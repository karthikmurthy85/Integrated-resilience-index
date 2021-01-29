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

