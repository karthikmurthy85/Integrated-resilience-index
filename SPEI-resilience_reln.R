library(raster)
library(rgdal)
library(sf)
library(rgeos)
library(ggplot2)
library(reshape2)
library(stringr)

spei_na_fill <- function(speival){
  gf <- subset(speival, y %in% NA)
  if(nrow(gf) == 0){
    return(speival)
  }
  
  else{
    for(k in 1:nrow(gf)){
      gid <- as.numeric(rownames(gf)[k])
      speival$y[gid] <- mean(c(speival$y[gid-1],speival$y[gid+1]), na.rm=T)
    }
    return(speival)
  }
}

spei_xtrme_index <- function(spei_val, qq){
  df1 <- data.frame(x = 0, y = 0) 
  df3 <- data.frame(x = (length(spei_val)+1), y= 0)
  df2 <- data.frame(x = 1:length(spei_val), y = spei_val)
  df <- rbind(df1, df2, df3)
  curv_wet <- spei_na_fill(df)
  curv_wet$y[curv_wet$y < 0] <- 0
  #plot(curv_wet, ty='l')
  curv_wet1 <- rbind(curv_wet, curv_wet[1,])
  #plot(curv_wet1, ty='l')
  pw = Polygon(curv_wet1)
  pws = Polygons(list(pw),1)
  sps_wet = SpatialPolygons(list(pws))
  #plot(sps_wet)
  gwet <- gArea(sps_wet, byid = T)
  
  dry_df <- spei_val[spei_val < 0]
  dry_qq <- quantile(abs(dry_df), probs=qq, na.rm=T)
  dff1 <- data.frame(x = 0, y = -dry_qq) 
  dff3 <- data.frame(x = (length(spei_val)+1), y= -dry_qq)
  dff2 <- data.frame(x = 1:length(spei_val), y = spei_val)
  dff <- rbind(dff1, dff2, dff3)
  rownames(dff) <- 1:nrow(dff)
  curv_dry <- spei_na_fill(dff)
  curv_dry$y[curv_dry$y >= -dry_qq] <- -dry_qq
  #plot(curv_dry, ty='l')
  curv_dry1 <- rbind(curv_dry, curv_dry[1,])
  #plot(curv_dry1, ty='l')
  pw1 = Polygon(curv_dry1)
  pws1 = Polygons(list(pw1),1)
  sps_dry = SpatialPolygons(list(pws1))
  #plot(sps_dry)
  gdry <- gArea(sps_dry, byid = T)
  
  return(gdry/gwet)
  
}




spei1 <- brick('global_veg/rasters/CRU_SPEI.tif')
d1 <- read.csv('/media/karthik/ADATA HD720/global_veg/datasets/NDVIt_equi_ncycl_ts.csv',
               h=F)
colnames(d1) <- c('id','nodata','x','y', 'Ncycle','CC','aic', 'Radj',
                      'model_arb', 'model_CC_aic', 'model_Radj',
                      'mn_rg', 'mn_rr', 'mn_rrt', 'mn_rtp', 'rtn',
                      'sd_rg', 'sd_rr', 'sd_rrt', 'sd_rtp')

library(foreach)
library(doParallel)
registerDoParallel(4)
foreach(i = 1:nrow(d1))%dopar%{
  library(reshape2);library(rgeos)
  cord <- matrix(as.numeric(d1[i, 3:4]), ncol=2)
  speidt <- melt(extract(spei1, cord))
  dry_indx075 <- spei_xtrme_index(speidt$value, 0.75)
  dry_indx080 <- spei_xtrme_index(speidt$value, 0.8)
  dry_indx085 <- spei_xtrme_index(speidt$value, 0.85)
  dry_indx090 <- spei_xtrme_index(speidt$value, 0.9)
  dry_indx095 <- spei_xtrme_index(speidt$value, 0.95)
  
  dry_indx_dat <- matrix(c(i,cord,dry_indx075,dry_indx080,
                           dry_indx085,dry_indx090,dry_indx095), byrow=T, nrow = 1)
  path <- '/media/karthik/ADATA HD720/global_veg/datasets/water_stress_new_dat.csv'
  write.table(dry_indx_dat, path, sep=',', 
              append=TRUE, row.names=FALSE, col.names=FALSE)
}

spei_dat <- read.csv('/media/karthik/ADATA HD720/global_veg/datasets/water_stress_new_dat.csv', h=F)
head(spei_dat)
ids <- spei_dat$V1
ids1 <- as.numeric(rownames(d1))
ids2 <- ids1[!ids1 %in% ids]

for(i in 3:length(ids2)){
  library(reshape2);library(rgeos)
  cord <- matrix(as.numeric(d1[ids2[i], 3:4]), ncol=2)
  speidt <- melt(extract(spei1, cord))
  dry_indx075 <- spei_xtrme_index(speidt$value, 0.75)
  dry_indx080 <- spei_xtrme_index(speidt$value, 0.8)
  dry_indx085 <- spei_xtrme_index(speidt$value, 0.85)
  dry_indx090 <- spei_xtrme_index(speidt$value, 0.9)
  dry_indx095 <- spei_xtrme_index(speidt$value, 0.95)
  
  dry_indx_dat <- matrix(c(i,cord,dry_indx075,dry_indx080,
                           dry_indx085,dry_indx090,dry_indx095), byrow=T, nrow = 1)
  path <- '/media/karthik/ADATA HD720/global_veg/datasets/water_stress_new_dat.csv'
  write.table(dry_indx_dat, path, sep=',', 
              append=TRUE, row.names=FALSE, col.names=FALSE)
  print(i)
}

spei_xtrm <- read.csv( '/media/karthik/ADATA HD720/global_veg/datasets/water_stress_new_dat.csv', h=F)
nrow(spei_xtrm)
whb <- brick('global_veg/rasters/whittaker_worldclim_ecoregion.tif')
plot(whb)

head(d1)
d2 <- d1
d2$id <- d2$Ncycle <- d2$CC <- d2$aic <- d2$Radj <- NULL
d2$model_arb <- d2$model_CC_aic <- d2$model_Radj <- NULL

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

d2$mn_rg <- 1 / d1$mn_rg
d2$sd_rg <- 1 / d1$sd_rg
d2 <- subset(d2, nodata < 0.30)
d3 <- d2
d2$nodata <- d2$x <- d2$y <- NULL

rsl_idx_dt <- apply(d2, 2, function(x) scale(x))
rt_idx <- apply(rsl_idx_dt, 1, function(x) mean(x, na.rm = T))
rsl_idx <- (exp(-rt_idx))

##rownames of d1 is the V1 column of spei_xtrm
spei_xtrm1 <- spei_xtrm[!duplicated(spei_xtrm),]
spei_xtrm1$V1 <- paste(spei_xtrm1$V2, spei_xtrm1$V3, sep='_')
spei_xtrm1$resl_indx <- NA

df_spei_resl <- data.frame(matrix(NA, ncol = ncol(spei_xtrm1), nrow = nrow(spei_xtrm1)))
for(i in 1:nrow(d3)){
  xy <- paste(d3$x[i], d3$y[i], sep='_')
  df <- subset(spei_xtrm1, V1 %in% xy)
  df$resl_indx <- rsl_idx[i]
  df_spei_resl[i,] <- df
  print(i)
}
head(df_spei_resl)
df_spei_resl$biome <- extract(whb, df_spei_resl[,c(2,3)])
df_spei_resl$X1 <- NULL

resl_clim_xtrem_rel <- function(extr_qq_id, extr_dat){
  prb <- seq(0, 1, 0.025)
  bm_slp <- matrix(NA, ncol=20, nrow=10)
  for(i in 1:10){
    df_whb <- subset(extr_dat, biome %in% i)
    qq_spei <- melt(quantile(df_whb[,extr_qq_id], probs=prb))
    resl075 <- vector()
    for(j in 1:(nrow(qq_spei)-1)){
      spdf1 <- subset(df_whb, df_whb[,extr_qq_id] > qq_spei$value[j] &
                        df_whb[,extr_qq_id] <= qq_spei$value[j+1])
      resl075[j] <- mean(spdf1$X9)
    }
    mdf <- data.frame(mn_resl = resl075, t = 1:length(resl075))
    md <- confint(lm(mn_resl ~ t, data = mdf), level = 0.9)[2,]
    bm_slp[i,] <- runif(20, md[1], md[2])
  }
  return(bm_slp)
}

slp075_spei <- resl_clim_xtrem_rel(3, df_spei_resl)
slp080_spei <- resl_clim_xtrem_rel(4, df_spei_resl)
slp085_spei <- resl_clim_xtrem_rel(5, df_spei_resl)
slp090_spei <- resl_clim_xtrem_rel(6, df_spei_resl)
slp095_spei <- resl_clim_xtrem_rel(7, df_spei_resl)

slp_spei_xtr <- data.frame(cbind(slp090_spei, slp095_spei))
slp_spei_xtr$bm <- 1:10
slp_spei_xtr1 <- melt(slp_spei_xtr, id.vars = 'bm') 

df1_spei_resl_mn <- melt(with(slp_spei_xtr1,
                             tapply(value, factor(bm), function(x) mean(x,na.rm=T))))[,2]
df1_spei_resl_sd <- melt(with(slp_spei_xtr1,
                             tapply(value, factor(bm), function(x) sd(x,na.rm=T))))[,2]
df1_spei_resl <- data.frame(spei_xtr_mn = df1_spei_resl_mn,
                               spei_xtr_sd = df1_spei_resl_sd,
                              biome = factor(str_pad(1:10, width=2, pad='0')))
y_lab <- expression(Delta ~ "resilience /" ~ Delta ~ "water stress extreme")
xt2 <- ggplot(df1_spei_resl, aes(y = spei_xtr_mn, x =biome,
                            colour = biome))+
  geom_hline(yintercept=0, linetype=2)+
  geom_point(size=4)+
  geom_errorbar(aes(ymin=spei_xtr_mn-spei_xtr_sd,
                    ymax= spei_xtr_mn+spei_xtr_sd), width=.2,
                position=position_dodge(.9))+
  #coord_cartesian(ylim=c(0,4))+ 
  ylab(y_lab)+ xlab("")+
  ggtitle("b) Resilience ~ Water stress")+
  scale_x_discrete(breaks = str_pad(1:10, width=2, pad='0'), 
                   label = c("Tropical \n rainforest","Tropical \n deciduous", 
                             "Tropical \n grassland", "Hot \n desert",
                             "Temperate \n rainforest","Temperate \n deciduous",
                             "Cold \n desert","Temperate \n grassland", "Taiga", "Tundra"))+
  scale_colour_manual(values= c("green4", "olivedrab3", "tan4", "red4", "cyan",
                              "royalblue4", "grey70", "yellow", "sandybrown", "steelblue3"))+
  theme_bw()+
  theme(
    legend.position = "none",
    panel.grid = element_blank(), axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle= 90, size=14, colour='black',vjust=0.70),
    axis.text.y = element_text(size=10, colour='black'),
    axis.title.y = element_text(size=14, colour='black'),
    plot.title = element_text(size=14, colour = 'black')
  )###800, 400, water_stress_biomes
xt2

matl <- matrix(c(1,1,1,2,2,2,2))
gridExtra::grid.arrange(xt1, xt2, layout_matrix=matl)

spei_xtrm090_mn <- melt(with(df_spei_resl,
                     tapply(X7, biome, mean)))[,2]
spei_xtrm090_sd <- melt(with(df_spei_resl,
                             tapply(X7, biome, sd)))[,2]
spei_xtrm_df090 <- data.frame(mn = spei_xtrm090_mn,
                              sd = spei_xtrm090_sd,
                              biome = factor(str_pad(1:10, pad = '0', width = 2)))
wxt1 <- ggplot(spei_xtrm_df090, aes(x = biome, y = mn, colour=biome))+
  geom_hline(yintercept = 0, linetype=2)+
  geom_point(size=4)+
  geom_errorbar(aes(ymin=mn-sd, ymax=mn+sd),width=.2,
                position=position_dodge(.9))+
  ylab("water extreme index")+ xlab("")+
  ggtitle("c) water extremes (90%) ")+
  scale_x_discrete(breaks = str_pad(1:10, width=2, pad='0'), 
                   label = c("Tropical \n rainforest","Tropical \n deciduous", 
                             "Tropical \n grassland", "Hot \n desert",
                             "Temperate \n rainforest","Temperate \n deciduous",
                             "Cold \n desert","Temperate \n grassland", "Taiga", "Tundra"))+
  scale_colour_manual(values= c("green4", "olivedrab3", "tan4", "red4", "cyan",
                                "royalblue4", "grey70", "yellow", "sandybrown", "steelblue3"))+
  theme_bw()+
  theme(
    legend.position = "none",
    panel.grid = element_blank(), axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle= 90, size=14, colour='black',vjust=0.70),
    axis.text.y = element_text(size=10, colour='black'),
    axis.title.y = element_text(size=14, colour='black'),
    plot.title = element_text(size=14, colour = 'black')
  )
  

spei_xtrm090_mn <- melt(with(df_spei_resl,
                             tapply(X8, biome, mean)))[,2]
spei_xtrm090_sd <- melt(with(df_spei_resl,
                             tapply(X8, biome, sd)))[,2]
spei_xtrm_df095 <- data.frame(mn = spei_xtrm090_mn,
                              sd = spei_xtrm090_sd,
                              biome = factor(str_pad(1:10, pad = '0', width = 2)))
wxt2 <- ggplot(spei_xtrm_df095, aes(x = biome, y = mn, colour=biome))+
  geom_hline(yintercept = 0, linetype=2)+
  geom_point(size=4)+
  geom_errorbar(aes(ymin=mn-sd, ymax=mn+sd),width=.2,
                position=position_dodge(.9))+
  ylab("water extreme index")+ xlab("")+
  ggtitle("d) water extremes (95%) ")+
  scale_x_discrete(breaks = str_pad(1:10, width=2, pad='0'), 
                   label = c("Tropical \n rainforest","Tropical \n deciduous", 
                             "Tropical \n grassland", "Hot \n desert",
                             "Temperate \n rainforest","Temperate \n deciduous",
                             "Cold \n desert","Temperate \n grassland", "Taiga", "Tundra"))+
  scale_colour_manual(values= c("green4", "olivedrab3", "tan4", "red4", "cyan",
                                "royalblue4", "grey70", "yellow", "sandybrown", "steelblue3"))+
  theme_bw()+
  theme(
    legend.position = "none",
    panel.grid = element_blank(), axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle= 90, size=14, colour='black',vjust=0.70),
    axis.text.y = element_text(size=10, colour='black'),
    axis.title.y = element_blank(),
    plot.title = element_text(size=14, colour = 'black')
  )
wxt2
gridExtra::grid.arrange(wxt1, wxt2, ncol=1)
mtl1 <- matrix(c(1, 2, 1, 2,3, 4, 3,4,3,4), byrow=T, ncol=2)
gridExtra::grid.arrange(txtr1, txtr2, wxt1, wxt2, layout_matrix=mtl1)


spei_mod090 <- summary(lm(X7~factor(biome), data = df_spei_resl)) 
spei_mod095 <- summary(lm(X8~factor(biome), data = df_spei_resl)) 
