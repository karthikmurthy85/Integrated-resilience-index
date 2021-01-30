plot(speidt$value, ty='l')
v1 <- speidt$value
v1[v1 > 0] <- 0
plot(v1, ty='l')
v2 <- c(0, v1, 0)

head(v1)

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


spei_seas_xtrme_index <- function(spei_val, qq){
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
  
  
  ar1 <- 0
  v1 <- spei_val
  v1[v1 > 0] <- 0
  v2 <- c(0, v1, 0)
  for(i in 1:12){
    mn1 <- v1[seq(i, length(v1), 12)]
    options(digits = 10)
    mnx <- (-1*quantile(abs(mn1[mn1 < 0]),
                                   probs=c(qq), na.rm=T))
    mnx1 <- v2[seq(i, length(v1), 12)[which(mn1 < mnx)]+1]
    if(length(mnx1) > 0){
      for(j in 1:length(mnx1)){
        poly <- Polygon(data.frame(y = c(mnx, mnx1[j], mnx), x = c(1:3)))
        ar1 <- ar1 + poly@area
      }
    }
  
    print(i)
  }
  
  return(ar1/gwet)
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
 
  dry_indx090 <- spei_seas_xtrme_index(speidt$value, 0.9)
  dry_indx095 <- spei_seas_xtrme_index(speidt$value, 0.95)
  
  dry_seas_indx_dat <- matrix(c(i,cord,dry_indx090,dry_indx095), byrow=T, nrow = 1)
  path <- '/media/karthik/ADATA HD720/global_veg/datasets/water_stress_seasonal_new_dat.csv'
  write.table(dry_seas_indx_dat, path, sep=',', 
              append=TRUE, row.names=FALSE, col.names=FALSE)
}


library(ggplot2)
library(raster)
spei_xtrm <- read.csv( 'D:/global_veg/datasets/water_stress_seasonal_new_dat.csv', h=F)
nrow(spei_xtrm)
whb <- brick('D:/global_veg/rasters/whittaker_worldclim_ecoregion.tif')
plot(whb)

d1 <- read.csv('D:/global_veg/datasets/NDVIt_equi_ncycl_ts.csv',
               h=F)
colnames(d1) <- c('id','nodata','x','y', 'Ncycle','CC','aic', 'Radj',
                  'model_arb', 'model_CC_aic', 'model_Radj',
                  'mn_rg', 'mn_rr', 'mn_rrt', 'mn_rtp', 'rtn',
                  'sd_rg', 'sd_rr', 'sd_rrt', 'sd_rtp')
head(d1)
d2 <- d1
d2$id <- d2$Ncycle <- d2$CC <- d2$aic <- d2$Radj <- NULL
d2$model_arb <- d2$model_CC_aic <- d2$model_Radj <- NULL

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
  library(reshape2)
  prb <- seq(0, 1, 0.025)
  bm_slp <- matrix(NA, ncol=20, nrow=10)
  for(i in 1:10){
    df_whb <- subset(extr_dat, biome %in% i)
    qq_spei <- melt(quantile(df_whb[,extr_qq_id], probs=prb))
    resl075 <- vector()
    for(j in 1:(nrow(qq_spei)-1)){
      spdf1 <- subset(df_whb, df_whb[,extr_qq_id] > qq_spei$value[j] &
                        df_whb[,extr_qq_id] <= qq_spei$value[j+1])
      resl075[j] <- mean(spdf1$X6)
    }
    mdf <- data.frame(mn_resl = resl075, t = 1:length(resl075))
    md <- confint(lm(mn_resl ~ t, data = mdf), level = 0.9)[2,]
    bm_slp[i,] <- runif(20, md[1], md[2])
  }
  return(bm_slp)
}


slp090_spei <- resl_clim_xtrem_rel(3, df_spei_resl)
slp095_spei <- resl_clim_xtrem_rel(4, df_spei_resl)

slp_spei_xtr <- data.frame(cbind(slp090_spei, slp095_spei))
slp_spei_xtr$bm <- 1:10
slp_spei_xtr1 <- melt(slp_spei_xtr, id.vars = 'bm') 

library(stringr)
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

