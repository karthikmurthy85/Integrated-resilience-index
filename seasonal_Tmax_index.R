library(raster)
library(rgdal)

tmp_seas_xtrme_index <- function(tmp_val, qq){
  at1 <- 0
  v1 <- tmp_val
  v2 <- c(0, v1, 0)
  for(i in 1:12){
    mn1 <- v1[seq(i, length(v1), 12)]
    options(digits = 10)
    mnx <- quantile(mn1,probs=c(qq), na.rm=T)
    mnx1 <- v2[seq(i, length(v1), 12)[which(mn1 > mnx)]+1]
    if(length(mnx1) > 0){
      for(j in 1:length(mnx1)){
        poly <- Polygon(data.frame(y = c(mnx, mnx1[j], mnx), x = c(1:3)))
        at1 <- at1 + poly@area
      }
    }
    
    print(i)
  }
  return(at1)
}


d1 <- read.csv('D:/global_veg/datasets/NDVIt_equi_ncycl_ts.csv',
               h=F)
colnames(d1) <- c('id','nodata','x','y', 'Ncycle','CC','aic', 'Radj',
                  'model_arb', 'model_CC_aic', 'model_Radj',
                  'mn_rg', 'mn_rr', 'mn_rrt', 'mn_rtp', 'rtn',
                  'sd_rg', 'sd_rr', 'sd_rrt', 'sd_rtp')

tmpmx <- brick('D:/global_veg/rasters/CRU_Tmax.tif')

library(foreach)
library(doParallel)
registerDoParallel(8)
foreach(i = 1:nrow(d1))%dopar%{
  library(reshape2);library(rgeos); library(raster)
  cord <- matrix(as.numeric(d1[i, 3:4]), ncol=2)
  tmpdt <- melt(extract(tmpmx, cord))
  
  tmp_indx090 <- tmp_seas_xtrme_index(tmpdt$value, 0.9)
  tmp_indx095 <- tmp_seas_xtrme_index(tmpdt$value, 0.95)
  
  tmp_indx_dat <- matrix(c(i,cord,tmp_indx090,tmp_indx095), byrow=T, nrow = 1)
  path <- 'D:/global_veg/datasets/Tmax_extr_seasonal_new_dat.csv'
  write.table(tmp_indx_dat, path, sep=',', 
              append=TRUE, row.names=FALSE, col.names=FALSE)
}



tmpxt_indx_df <- read.csv('D:/global_veg/datasets/Tmax_extr_seasonal_new_dat.csv', h=F)
head(tmpxt_indx_df)
whb <- brick('D:/global_veg/rasters/whittaker_worldclim_ecoregion.tif')
plot(whb)

d1 <- read.csv('D:/global_veg/datasets/NDVIt_equi_ncycl_ts.csv',
               h=F)
colnames(d1) <- c('id','nodata','x','y', 'Ncycle','CC','aic', 'Radj',
                  'model_arb', 'model_CC_aic', 'model_Radj',
                  'mn_rg', 'mn_rr', 'mn_rrt', 'mn_rtp', 'rtn',
                  'sd_rg', 'sd_rr', 'sd_rrt', 'sd_rtp')

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

##rownames of d1 is the V1 column of tmpxt_indx_df
tmpxt_indx_df1 <- tmpxt_indx_df[!duplicated(tmpxt_indx_df),]
tmpxt_indx_df1$V1 <- paste(tmpxt_indx_df1$V2, tmpxt_indx_df1$V3, sep='_')
tmpxt_indx_df1$resl_indx <- NA

df_tmp_resl <- data.frame(matrix(NA, ncol = ncol(tmpxt_indx_df1), nrow = nrow(tmpxt_indx_df1)))
for(i in 1:nrow(d2)){
  xy <- paste(d3$x[i], d3$y[i], sep='_')
  df <- subset(tmpxt_indx_df1, V1 %in% xy)
  if(nrow(df)>0){
    df$resl_indx <- rsl_idx[i]
    df_tmp_resl[i,] <- df
  }
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
  scale_x_discrete(breaks = str_pad(1:10, width=2, pad='0'), 
    label = c("Tropical \n rainforest","Tropical \n deciduous", 
                             "Tropical \n grassland", "Hot \n desert",
                             "Temperate \n rainforest","Temperate \n deciduous",
                             "Cold \n desert","Temperate \n grassland", "Taiga", "Tundra"))+
  scale_colour_manual(values= c("green4", "olivedrab3", "tan3", "red4", "cyan",
                                "royalblue4", "grey70", "yellow3", "sandybrown", "steelblue3"))+
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

