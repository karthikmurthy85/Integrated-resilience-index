library(raster)
library(ggplot2)
library(reshape2)
library(stringr)
library(zoo)
library(ggpmisc)
library(reconPlots)

##Aggregate the GIMMS NDVI data to 55km
f1 <-  list.files('/media/karthik/ADATA HD720/global_veg/rasters/GIMMS_3gV1/', full.names = T)
f1 <- f1[-(1:2)]
#ndras <- aggregate(stack(f1), fact=6)
ndras <- brick('/media/karthik/ADATA HD720/global_veg/GIMMS_NDVI_55km.grd')

##GEt data on global vegetated area 
ndvi_area <- read.csv('/media/karthik/ADATA HD720/global_veg/datasets/veg_clim_rel_lm.csv')
head(ndvi_area)
ndvi_area[,1:9] <- NULL

##Annual mean NDVI of the long-term NDVI time series
mean_fun <- function(ts_ndvi){
  yr_id <- c(seq(1, 816, 24), 816)
  mn_yri <- NA
  for(j in 1:34){
    yr_ndvi <- ts_ndvi[yr_id[j]:(yr_id[j+1]-1)]
    mn_yri[j] <- mean(yr_ndvi)
  }
  return(mn_yri)
}

##Function to select the trend window of the stl function such that 
##the number of cycles in the trend NDVI and mean annual NDVI time series are equal

twindow_sel <- function(ndvi_ts, meanNDVI_ts){
  maxs <- match(meanNDVI_ts[ggpmisc:::find_peaks(meanNDVI_ts)], meanNDVI_ts) ##time -  Maxima occurred in time series
  mins <- match(meanNDVI_ts[ggpmisc:::find_peaks(-meanNDVI_ts)], meanNDVI_ts) ##time -  Minima occurred in time series
  pk <- c(mins, maxs) ## Time - all the extrema in the time series
  pk1 <- pk[order(pk)] 
  cyc_num1 <- floor(length(pk1)/2) ##Number of cycles in the mean annual NDVI time series data
  
  cyc_num <- NA
  tws <- seq(1,5,0.25) Different trend window values
  for(i in 1:length(tws)){
    ndt1 <- ts(ndvi_ts, frequency = 24) ##GIMMS data is bi-monthly : 2 NDVI values per month; hence 24 data points annually 
    ndt2 <- stl(ndt1, s.window = 24, t.window = 24*tws[i]) ##Conduct stl decomposition
    plot(ndt2, main = "NDVI time-series decomposition")
    ts5 <- ndt2$time.series[,2]
    maxs <- match(ts5[ggpmisc:::find_peaks(ts5)], ts5) ##time -  Maxima occurred in time series
    mins <- match(ts5[ggpmisc:::find_peaks(-ts5)], ts5) ##time -  Minima occurred in time series
    pks <- c(mins, maxs)  ## Time - all the extrema in the time series
    pks1 <- pks[order(pks)]
    cyc_num[i] <- floor(length(pks1)/2) ##Number of cycles in the trend NDVI time series data
    print(i)
  }
  cyc_df <- data.frame(cyc_num=cyc_num, tws = tws)
  cyc_df1 <- subset(cyc_df, cyc_num %in% cyc_num1)  ##Get the trend window value which gives equal number of cycles in mean annual NDVI
  
  if(nrow(cyc_df1)==0){
    cyc_df$diff <- abs(cyc_df$cyc_num - cyc_num1)
    cydf <- subset(cyc_df, diff %in% min(diff))
    return(mean(cydf$tws))
  }
  
  else{
    return(cyc_df1$tws[1])
  }
  
}

library(foreach)
library(doParallel)
registerDoParallel(4)
foreach(k = 1:nrow(ndvi_area))%dopar%{
#for(k in 1:2000){
  library(raster)
  library(ggplot2)
  library(reshape2)
  library(stringr)
  library(zoo)
  library(ggpmisc)
  library(reconPlots)
  library(nlme)
  library(stringr)
  library(segmented)
  
  source('global_veg/new_codes/RR_indices.R')
  
  lip <- k
  coord <- matrix(as.numeric(ndvi_area[k, 1:2]), ncol=2)
  ts1 <- extract(ndras, coord)/10000
  ts1[ts1 < 0] <- 0
  nodata <- length(ts1[ts1==0])/length(ts1)
  ts2 <- melt(ts1)
  
  mn <- mean_fun(ts2$value)
  mn[is.na(mn)] <- 0
  plot(mn, ty='o')
  
  
  
  mn_dlog <- mn_dlog_inv <- mn_log <- mn_log_inv <- NA
  mn_mod <- summary(lm(mn~c(1:34)))
  mn_mod
  
  ##GEt the details regarding the trajectory of mean annual NDVI#
  source('global_veg/codes/dlogistic.R')
  source('global_veg/codes/dlogistic_inv.R')
  source('global_veg/codes/logistic.R')
  source('global_veg/codes/logistic_inv.R')
  source('global_veg/new_codes/predict_index.R')
  source('global_veg/new_codes/model_sel_functions.R')
  
  dlog_mod <- mn_dlog <- mn_aic_dlog <- mn_radj_dlog <- NA
  try(dlog_mod <- dLogistic(mn))
  try(mn_dlog <- cc_score(mn, dlog_mod))
  try(mn_aic_dlog <- aic_calc(mn, mn - dlog_mod, 7))
  try(mn_radj_dlog <- R_adj_fun(mn, dlog_mod,7))
  
  dlog_inv_mod <- mn_dlog_inv <- mn_aic_dlog_inv <-mn_radj_dlog_inv <- NA
  try(dlog_inv_mod <- dLogistic_inv(mn))
  try(mn_dlog_inv <- cc_score(mn, dlog_inv_mod))
  try(mn_aic_dlog_inv <- aic_calc(mn, mn - dlog_inv_mod, 7))
  try(mn_radj_dlog_inv <- R_adj_fun(mn, dlog_inv_mod, 7))
  
  log_mod <- mn_log <- mn_aic_log <- mn_radj_log <- NA
  try(log_mod <- logistic(mn))
  try(mn_log <- cc_score(mn, log_mod))
  try(mn_aic_log <- aic_calc(mn, mn - log_mod, 4))
  try(mn_radj_log <- R_adj_fun(mn, log_mod, 4))
  
  log_inv_mod <- mn_log_inv <- mn_aic_log_inv <- mn_radj_log_inv <- NA
  try(log_inv_mod <- logistic_inv(mn))
  try(mn_log_inv <- cc_score(mn, log_inv_mod))
  try(mn_aic_log_inv <- aic_calc(mn, mn - log_inv_mod, 4))
  try(mn_radj_log_inv <- R_adj_fun(mn, log_inv_mod, 4))
  
  mn_lin <- cc_score(mn, mn_mod$coefficients[1] + mn_mod$coefficients[2]*(1:34))
  mn_aic_lin <- aic_calc(mn, mn_mod$residuals, 2)
  mn_radj_lin <- mn_mod$adj.r.squared
  
  
  mn_tr_dt <- data.frame(cc_score = c(mn_lin, mn_log, mn_log_inv, 
                                      mn_dlog, mn_dlog_inv),
                         aic_score = c(mn_aic_lin, mn_aic_log, mn_aic_log_inv, 
                                       mn_aic_dlog, mn_aic_dlog_inv),
                         radj_score = c(mn_radj_lin, mn_radj_log, mn_radj_log_inv, 
                                        mn_radj_dlog, mn_radj_dlog_inv),
                         models = c("linear", "logistic", "logistic_inverse",
                                    "dlogistic", "dlogistic_inverse"),
                         trend = c("gradual", "abrupt", "abrupt", 
                                   "reverse", "reverse"), 
                         rank =c(1, 2, 2, 3, 3))
  mn_tr_dt <- subset(mn_tr_dt, !cc_score %in% NA)
  
  ##Runing the model selection to identify the best fit##
  modsel <- model_sel_fun(mn_tr_dt)
  modsel
  mn_tr_dt$cc_score <- round(mn_tr_dt$cc_score,3)
  mn_tr_dt
  modsel1 <- model_sel_fun2(mn_tr_dt)
  modsel1
  modsel2 <- model_sel_fun3(mn_tr_dt)
  modsel2
  
  ##Ananlysis of trend NDVI
  tw_num <- twindow_sel(ts2$value,mn) ##Identify the trend window value
  ndt1 <- ts(ts2$value, frequency = 24)
  ndt2 <- stl(ndt1, s.window = 24, t.window = 24*tw_num)
  plot(ndt2, main = "NDVI time-series decomposition")
  ts5 <- ndt2$time.series[,2]
  ts6 <- data.frame(ts = ts5, id = 1:length(ts5))
  plot(ts6$id, ts6$ts, ty='l')
  plot(mn, ty='o')
  plot(ts6$id, ts6$ts, ty='l')
  
  ##Identify cycles
  fitt <- ts5
  maxs <- match(fitt[ggpmisc:::find_peaks(fitt)], fitt)
  mins <- match(fitt[ggpmisc:::find_peaks(-fitt)], fitt)
  pks <- c(mins, maxs)
  pks1 <- pks[order(pks)]
  abline(v=pks1, col='blue')
  abline(h=mean(mn), col ='red4')
  
  cyc_num <- floor(length(pks1)/2)
  rem <- length(pks1)%%2
  
  ##Calculate the return indices for each cycle
  ifelse(maxs[1]<mins[1], ky <- 1, ky <- -1)
  source('global_veg/new_codes/RR_indices.R')
  res <- matrix(data = NA, nrow = cyc_num, ncol = 4)
  lindt1 <-  data.frame(x = 1:34, y = mean(mn))
  
  if(ky == 1){
    num_cyc <- cyc_num 
    ifelse(rem == 0, pks2 <- c(pks1, 34), pks2 <- pks1)
    for(i in 1:num_cyc){
      cyc <- pks2[c((2*i-1), (2*i+1))]
      cyci <- data.frame(nd = c(fitt[cyc[1] : cyc[2]]),
                         t = c(cyc[1] : cyc[2]),
                         id = 1:length(c(cyc[1] : cyc[2])))
      cyci1 <- subset(cyci, !nd %in% NA)
      res[i,] <- rr_fun(cyci1,ky,lindt1)
      print(i)
    }
  }
  
  if(ky == -1){
    
    num_cyc <- cyc_num 
    ifelse(rem == 0, pks2 <- c(pks1, 34), pks2 <- pks1)
    
    for(i in 1:num_cyc){
      cyc <- pks2[c((2*i-1), (2*i+1))]
      cyci <- data.frame(nd = c(fitt[cyc[1] : cyc[2]]),
                         t = c(cyc[1] : cyc[2]),
                         id = 1:length(c(cyc[1] : cyc[2])))
      cyci1 <- subset(cyci, !nd %in% NA)
      res[i,] <- rr_fun(cyci1,ky,lindt1)
      print(i)
    }
  }
  
  ##Calculate multi-return index based on the mean and standard deviation of the return indices 
  ##that are calculated for all cycles in the trend NDVI time series
  ifelse(maxs[1]<mins[1], rtp <- pk_dif(maxs), rtp <- pk_dif(mins))
  ifelse(maxs[1]<mins[1], rr <- rrg_fun(fitt[maxs], mean(mn)), 
         rr <- rrg_fun(fitt[mins], mean(mn)))  
  
  tim <- c(res[,3], res[,4])
  t1 <- tim[order(tim)]
  t1 <- t1[!is.na(t1)]
  t2 <- NA
  for(j in 1:(length(t1)-1)){t2[j] <- t1[j+1]-t1[j]}
  rg1 <- res[,1]
  rg1 <- rg1[!rg1 %in% 0]
  mn_rg <- mean(rg1) ##RG : rate of return
  mn_rr <- mean(rr) ##RR : distance b/w steady state and equllibrium
  mn_rrt <- mean(t2) ## Time difference b/w the intersection points of time series with long-term mean
  mn_rtp <- mean(rtp) ##Time difference between the steady states
  rtn <- 1 - (length(t1)/length(tim))
  sd_rg <- sd(rg1)
  sd_rr <- sd(rr)
  sd_rrt <- sd(t2, na.rm = T)
  sd_rtp <- sd(rtp)
  
  mod_var <- NA
  try(mod_var <- (c(melt(modsel[1:3])[,2],as.character(melt(modsel[4])[1,]),modsel1,  modsel2)))
  ifelse(mod_var %in% NA, 
         mod_var <- c(modsel[1], NA, NA, modsel[2], modsel1,  modsel2),
         mod_var <- mod_var)
  mnr_var <- c(mn_rg, mn_rr, mn_rrt, mn_rtp, rtn)
  sd_var <- c(sd_rg, sd_rr, sd_rrt, sd_rtp)
  
  vardt <- matrix(data = c(lip, nodata,  coord, num_cyc, mod_var, mnr_var, sd_var),
                  nrow=1)
  
  ##Write the data into csv file
  path = "/media/karthik/ADATA HD720/global_veg/datasets/NDVIt_equi_ncycl_ts.csv"  
  write.table(vardt, path, sep=',',
              append=TRUE, row.names=FALSE, col.names=FALSE)
  #print(vardt)
}

