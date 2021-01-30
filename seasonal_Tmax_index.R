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


