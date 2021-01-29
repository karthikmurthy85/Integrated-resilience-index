library(raster)
library(rgdal)

tmp_xtrme_index <- function(tmp_val, qq){
  tmp_qq <- quantile(tmp_val, probs=qq, na.rm=T)
  tff1 <- data.frame(x = 0, y = tmp_qq) 
  tff3 <- data.frame(x = (length(tmp_val)+1), y= tmp_qq)
  tff2 <- data.frame(x = 1:length(tmp_val), y = tmp_val)
  tff <- rbind(tff1, tff2, tff3)
  rownames(tff) <- 1:nrow(tff)
  curv_tmp <- spei_na_fill(tff)
  curv_tmp$y[curv_tmp$y <= tmp_qq] <- tmp_qq
  #plot(curv_tmp, ty='l')
  curv_tmp1 <- rbind(curv_tmp, curv_tmp[1,])
  #plot(curv_tmp1, ty='l')
  pw1 = Polygon(curv_tmp1)
  pws1 = Polygons(list(pw1),1)
  tmpp = SpatialPolygons(list(pws1))
  #plot(tmpp)
  txtr <- gArea(tmpp, byid = T)
  return(txtr)
}

d1 <- read.csv('/media/karthik/ADATA HD720/global_veg/datasets/NDVIt_equi_ncycl_ts.csv',
               h=F)
colnames(d1) <- c('id','nodata','x','y', 'Ncycle','CC','aic', 'Radj',
                  'model_arb', 'model_CC_aic', 'model_Radj',
                  'mn_rg', 'mn_rr', 'mn_rrt', 'mn_rtp', 'rtn',
                  'sd_rg', 'sd_rr', 'sd_rrt', 'sd_rtp')

tmpmx <- brick('global_veg/rasters/CRU_Tmax.tif')

library(foreach)
library(doParallel)
registerDoParallel(4)
foreach(i = 1:nrow(d1))%dopar%{
  library(reshape2);library(rgeos)
  cord <- matrix(as.numeric(d1[i, 3:4]), ncol=2)
  tmpdt <- melt(extract(tmpmx, cord))

  tmp_indx090 <- tmp_xtrme_index(tmpdt$value, 0.9)
  tmp_indx095 <- tmp_xtrme_index(tmpdt$value, 0.95)
  
  tmp_indx_dat <- matrix(c(i,cord,tmp_indx090,tmp_indx095), byrow=T, nrow = 1)
  path <- '/media/karthik/ADATA HD720/global_veg/datasets/Tmax_extr_new_dat.csv'
  write.table(tmp_indx_dat, path, sep=',', 
              append=TRUE, row.names=FALSE, col.names=FALSE)
}

tmpxt_indx_df <- read.csv('/media/karthik/ADATA HD720/global_veg/datasets/Tmax_extr_new_dat.csv', h=F)
head(tmpxt_indx_df)
whb <- brick('global_veg/rasters/whittaker_worldclim_ecoregion.tif')
plot(whb)

d2 <- d1
d2$id <- d2$Ncycle <- d2$CC <- d2$aic <- d2$Radj <- NULL
d2$model_arb <- d2$model_CC_aic <- d2$model_Radj <- NULL

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

d2$mn_rg <- 1 / d1$mn_rg
d2$sd_rg <- 1 / d1$sd_rg
d2 <- subset(d2, nodata < 0.40)
d3 <- d2
d2$nodata <- d2$x <- d2$y <- NULL

rsl_idx_dt <- apply(d2, 2, function(x) scale(x))
rt_idx <- apply(rsl_idx_dt, 1, function(x) mean(x, na.rm = T))
rsl_idx <- (exp(-rt_idx))

##rownames of d1 is the V1 column of tmpxt_indx_df
tmpxt_indx_df1 <- tmpxt_indx_df[!duplicated(tmpxt_indx_df),]
tmpxt_indx_df1$V1 <- paste(tmpxt_indx_df1$V2, tmpxt_indx_df1$V3, sep='_')
tmpxt_indx_df1$resl_indx <- NA

df_tmp_resl <- data.frame(matrix(NA, ncol = ncol(tmpxt_indx_df1), nrow = nrow(tmpxt_indx_df1)))
for(i in 1:nrow(d2)){
  xy <- paste(d3$x[i], d3$y[i], sep='_')
  df <- subset(tmpxt_indx_df1, V1 %in% xy)
  df$resl_indx <- rsl_idx[i]
  df_tmp_resl[i,] <- df
  print(i)
}
head(df_tmp_resl)
df_tmp_resl$biome <- extract(whb, df_tmp_resl[,c(2,3)])
df_tmp_resl$X1 <- NULL

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
      resl075[j] <- mean(spdf1$X6) ##X6 is the resilience column
    }
    mdf <- data.frame(mn_resl = resl075, t = 1:length(resl075))
    md <- confint(lm(mn_resl ~ t, data = mdf), level = 0.9)[2,]
    bm_slp[i,] <- runif(20, md[1], md[2])
  }
  return(bm_slp)
}
library(reshape2); library(stringr)
slp090_tmx <- resl_clim_xtrem_rel(3, df_tmp_resl)
slp095_tmx <- resl_clim_xtrem_rel(4, df_tmp_resl)

slp_tmx_xtr <- data.frame(cbind(slp090_tmx, slp095_tmx))
slp_tmx_xtr$bm <- 1:10
slp_tmx_xtr1 <- melt(slp_tmx_xtr, id.vars = 'bm') 
head(slp_tmx_xtr1)

slp_tmx_xtr2 <- melt(unlist(with(slp_tmx_xtr1, 
                            tapply(value, bm,
                                          function(x) c(mean(x), sd(x))))))[,1]
slp_tmx_xtr_mn <- slp_tmx_xtr2[seq(1, 20, 2)]
slp_tmx_xtr_sd <- slp_tmx_xtr2[seq(2, 20, 2)]
slp_tmx_xtr_df <- data.frame(beta_mean = slp_tmx_xtr_mn, 
                             beta_sd = slp_tmx_xtr_sd,
                             biome = factor(str_pad(1:10, width=2, pad='0')))

y_lab <- expression(Delta ~ "resilience /" ~ Delta ~ "temperature extreme")
library(ggplot2)
xt1 <- ggplot(slp_tmx_xtr_df, aes(y = beta_mean, x = biome,
                                 fill = biome, colour= biome))+
  geom_hline(yintercept=0, linetype=2)+
  geom_point(size=4)+
  geom_errorbar(aes(ymin=beta_mean-beta_sd,
                    ymax= beta_mean+beta_sd), width=.2,
                position=position_dodge(.9))+
  #coord_cartesian(ylim=c(0,4))+ 
  ylab(y_lab)+ xlab("")+
  ggtitle("a) Resilience ~ temperature extreme ")+
  scale_x_discrete(label = c("Tropical \n rainforest","Tropical \n deciduous", 
                             "Tropical \n grassland", "Hot \n desert",
                             "Temperate \n rainforest","Temperate \n deciduous",
                             "Cold \n desert","Temperate \n grassland", "Taiga", "Tundra"))+
  scale_colour_manual(values= c("green4", "olivedrab3", "tan3", "red4", "cyan",
                              "royalblue4", "grey70", "yellow4", "sandybrown", "steelblue3"))+
  theme_bw()+
  theme(
    legend.position = "none",
    panel.grid = element_blank(), axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    #axis.text.x = element_text(angle= 90, size=14, colour='black',vjust=0.70),
    axis.text.y = element_text(size=10, colour='black'),
    axis.title.y = element_text(size=14, colour='black'),
    plot.title = element_text(size=14, colour = 'black')
  )###800, 400, water_stress_biomes
xt1
nrow(d2)


df1_tmp_resl_mn <- melt(with(df_tmp_resl,
                     tapply(X4, factor(biome), function(x) mean(x,na.rm=T))))[,2]
df1_tmp_resl_sd <- melt(with(df_tmp_resl,
                        tapply(X4, factor(biome), function(x) sd(x,na.rm=T))))[,2]
df1_tmp_resl090 <- data.frame(tmpx_xtr_mn = df1_tmp_resl_mn,
                           tmpx_xtr_sd = df1_tmp_resl_sd,
                           biome = factor(str_pad(1:10, width=2, pad='0')))

txtr1 <- ggplot(df1_tmp_resl090, aes(y = tmpx_xtr_mn, x = biome,
                           fill = biome, colour= biome))+
  geom_hline(yintercept=0, linetype=2)+
  geom_point(size=4)+
  geom_errorbar(aes(ymin=tmpx_xtr_mn-tmpx_xtr_sd,
                    ymax= tmpx_xtr_mn+tmpx_xtr_sd), width=.2,
                position=position_dodge(.9))+
  #coord_cartesian(ylim=c(0,4))+ 
  ylab("Temperature extreme index")+ xlab("")+
  ggtitle("a) Temperature extremes (90%) ")+
  scale_x_discrete(label = c("Tropical \n rainforest","Tropical \n deciduous", 
                             "Tropical \n grassland", "Hot \n desert",
                             "Temperate \n rainforest","Temperate \n deciduous",
                             "Cold \n desert","Temperate \n grassland", "Taiga", "Tundra"))+
  scale_colour_manual(values= c("green4", "olivedrab3", "tan3", "red4", "cyan",
                                "royalblue4", "grey70", "yellow4", "sandybrown", "steelblue3"))+
  theme_bw()+
  theme(
    legend.position = "none",
    panel.grid = element_blank(), axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    #axis.text.x = element_text(angle= 90, size=14, colour='black',vjust=0.70),
    axis.text.y = element_text(size=10, colour='black'),
    axis.title.y = element_text(size=14, colour='black'),
    plot.title = element_text(size=14, colour = 'black')
  )


df1_tmp_resl_mn <- melt(with(df_tmp_resl,
                             tapply(X5, factor(biome), function(x) mean(x,na.rm=T))))[,2]
df1_tmp_resl_sd <- melt(with(df_tmp_resl,
                             tapply(X5, factor(biome), function(x) sd(x,na.rm=T))))[,2]
df1_tmp_resl095 <- data.frame(tmpx_xtr_mn = df1_tmp_resl_mn,
                              tmpx_xtr_sd = df1_tmp_resl_sd,
                              biome = factor(str_pad(1:10, width=2, pad='0')))

txtr2 <- ggplot(df1_tmp_resl095, aes(y = tmpx_xtr_mn, x = biome,
                            fill = biome, colour= biome))+
  geom_hline(yintercept=0, linetype=2)+
  geom_point(size=4)+
  geom_errorbar(aes(ymin=tmpx_xtr_mn-tmpx_xtr_sd,
                    ymax= tmpx_xtr_mn+tmpx_xtr_sd), width=.2,
                position=position_dodge(.9))+
  #coord_cartesian(ylim=c(0,4))+ 
  ylab("Temperature extreme index")+ xlab("")+
  ggtitle("b) Temperature extremes (95%) ")+
  scale_x_discrete(label = c("Tropical \n rainforest","Tropical \n deciduous", 
                             "Tropical \n grassland", "Hot \n desert",
                             "Temperate \n rainforest","Temperate \n deciduous",
                             "Cold \n desert","Temperate \n grassland", "Taiga", "Tundra"))+
  scale_colour_manual(values= c("green4", "olivedrab3", "tan3", "red4", "cyan",
                                "royalblue4", "grey70", "yellow4", "sandybrown", "steelblue3"))+
  theme_bw()+
  theme(
    legend.position = "none",
    panel.grid = element_blank(), axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    #axis.text.x = element_text(angle= 90, size=14, colour='black',vjust=0.70),
    axis.text.y = element_text(size=10, colour='black'),
    axis.title.y = element_blank(),
    plot.title = element_text(size=14, colour = 'black')
  )

gridExtra::grid.arrange(txtr1, txtr2, nrow=1)

head(df_tmp_resl)
evgr <- subset(df_tmp_resl, biome %in% 1)
nrow(evgr)
evgr_tmxpr <- melt(quantile(evgr$X4, probs=c(seq(0, 1, 0.025))))[,1]
resl <- matrix(NA, ncol=3, nrow=40)
for(i in 1:(length(evgr_tmxpr)-1)){
  resl[i,] <- c(mean(subset(evgr, 
                           X4 >= evgr_tmxpr[i] &
                             X4 < evgr_tmxpr[i+1])[,5]), evgr_tmxpr[i], evgr_tmxpr[i+1])
}
resl<- data.frame(resl)
resl$id <- 1:40
resl$rng <- paste(round(resl[,2],2),  round(resl[,3],2), sep='-')
demo1 <- ggplot(resl, aes(x = id, y = X1))+
  geom_point(size=1) + 
  geom_smooth(method = 'lm')+
  xlab("Extreme temperature index") +
  ylab("mean resilience")+
  scale_x_continuous(breaks= c(1, 10, 20, 30, 40),
                     labels=resl$rng[ c(1, 10, 20, 30, 40)])+
  theme_bw()+
  ggtitle("a) Tropical Evergreen: Resilience ~ Temperature extreme(90%)")+
  theme(
    axis.text = element_text(size=12, colour='black'),
    axis.title = element_text(size=14, colour='black'),
    plot.title = element_text(size=14, colour='black'),
    panel.grid = element_blank()
  )

evgr_tmxpr <- melt(quantile(evgr$X5, probs=c(seq(0, 1, 0.025))))[,1]
resl <- matrix(NA, ncol=3, nrow=40)
for(i in 1:(length(evgr_tmxpr)-1)){
  resl[i,] <- c(mean(subset(evgr, 
                            X5 >= evgr_tmxpr[i] &
                              X5 < evgr_tmxpr[i+1])[,5]), evgr_tmxpr[i], evgr_tmxpr[i+1])
}
resl<- data.frame(resl)
resl$id <- 1:40
resl$rng <- paste(round(resl[,2],2),  round(resl[,3],2), sep='-')
demo2 <- ggplot(resl, aes(x = id, y = X1))+
  geom_point(size=1) + 
  geom_smooth(method = 'lm')+
  xlab("Extreme temperature index") +
  ylab("mean resilience")+
  scale_x_continuous(breaks= c(1, 10, 20, 30, 40),
                     labels=resl$rng[ c(1, 10, 20, 30, 40)])+
  theme_bw()+
  ggtitle("b) Tropical Evergreen: Resilience ~ Temperature extreme(95%)")+
  theme(
    axis.text = element_text(size=12, colour='black'),
    axis.title = element_text(size=14, colour='black'),
    plot.title = element_text(size=14, colour='black'),
    panel.grid = element_blank()
  )
gridExtra::grid.arrange(demo1, demo2, nrow=1)


##Convert to 0 to 1 range##
tmp_mod090 <- summary(lm(X4~factor(biome), data = df_tmp_resl)) 
tmp_mod095 <- summary(lm(X5~factor(biome), data = df_tmp_resl)) 
