
res_mn <- read.csv('/media/karthik/ADATA HD720/global_veg/datasets/NDVIt_equi_ncycl_ts.csv',
                   h=F)
head(res_mn)
colnames(res_mn) <- c('id','nodata','x','y', 'Ncycle','CC','aic', 'Radj',
                      'model_arb', 'model_CC_aic', 'model_Radj',
                      'mn_rg', 'mn_rr', 'mn_rrt', 'mn_rtp', 'rtn',
                      'sd_rg', 'sd_rr', 'sd_rrt', 'sd_rtp')
null_df <- subset(res_mn, model_Radj %in% 'null' & CC < 0.2)
dlog_df <- subset(res_mn, model_Radj %in% 'dlogistic' & 
                    model_arb %in% 'dlogistic'&
                    model_CC_aic %in% 'dlogistic' & CC > 0.6)
dlog_inv_df <- subset(res_mn, model_Radj %in% 'dlogistic_inverse' & 
                        model_arb %in% 'dlogistic_inverse'&
                        model_CC_aic %in% 'dlogistic_inverse' & CC > 0.6)
log_df <- subset(res_mn, model_Radj %in% 'logistic' & 
                   model_arb %in% 'logistic'&
                   model_CC_aic %in% 'logistic' & CC > 0.65)
log_inv_df <- subset(res_mn, model_Radj %in% 'logistic_inverse' & 
                       model_arb %in% 'logistic_inverse'&
                       model_CC_aic %in% 'logistic_inverse' & CC > 0.65)
lin_df <- subset(res_mn, model_Radj %in% 'linear' & 
                   model_arb %in% 'linear'&
                   model_CC_aic %in% 'linear' & CC > 0.65)

revr_df <- rbind(dlog_df, dlog_inv_df)
abr_df <- rbind(log_df, log_inv_df)

sm_fun <- function(num_rows, nd){
  sm1 <- ceiling(runif(nd, 1, num_rows))
  sm2 <- ceiling(runif(nd, 1, num_rows))
  sm3 <- ceiling(runif(nd, 1, num_rows))
  
  sm4 <- unique(c(sm1, sm2, sm3))
  sm5 <- sample(sm4, nd, replace = F)
  return(sm5)
}

nsam <- 40
null1 <- null_df[sm_fun(nrow(null_df), nsam),]
abr1 <- abr_df[sm_fun(nrow(abr_df), nsam),]
revr1 <- revr_df[sm_fun(nrow(revr_df), nsam),]
grd1 <- lin_df[sm_fun(nrow(lin_df), nsam),]

com_dt <- rbind(null1, abr1, revr1, grd1)
com_dt$mov_mod <- rep(c('null','abrupt','reverse','gradual'), each= nsam)

# ggplot(com_dt, aes( y = (1/mn_rg), x = mov_mod, fill = mov_mod))+
#   geom_boxplot()+theme_bw() + ylim(0,1000)
# 
# ggplot(com_dt, aes( y = (1/sd_rg), x = mov_mod, fill = mov_mod))+
#   geom_boxplot()+theme_bw() + ylim(0,5000)
# 
# ggplot(com_dt, aes( y = mn_rtp, x = mov_mod, fill = mov_mod))+
#   geom_boxplot()+theme_bw()
# 
# ggplot(com_dt, aes( y = sd_rtp, x = mov_mod, fill = mov_mod))+
#   geom_boxplot()+theme_bw()
# 

# inds <- com_dt[,c(12:21)]
# inds$mn_rtp <- NULL
# inds$sd_rtp <- NULL
# inds$mn_rg <- 1 /inds$mn_rg
# inds$sd_rg <- 1 /inds$sd_rg
# for(i in 1:(ncol(inds)-1)){inds[,i] <- (exp(scale(inds[,i])))}
# inds <- inds[complete.cases(inds),]
# library(vegan)
# rtnmds<-metaMDS(inds[,-ncol(inds)],k=2,trymax=30, dist="euclid")
# rtnmds1<-metaMDS(inds[,-ncol(inds)],k=2,trymax=30, dist="bray")
# inds1 <- cbind(inds, rtnmds$points)
# inds2 <- cbind(inds, rtnmds1$points)
# 
# ellipse_eucl<-ggplot(inds1, aes(x=MDS1, y=MDS2, fill=mov_mod, colour=mov_mod)) +
#   geom_jitter(size=1, width=0.05, alpha=0.5)+ theme_bw()+ 
#   #ylim(-0.1,0.1)+xlim(-0.25,0.25)+
#   stat_ellipse(type = "t", level= 0.5, linetype=1, lwd=1.2, segments=6) + 
#   ylab("MDS2") + xlab("MDS1")+
#   ggtitle("a) Ordinance return time index")
# ellipse_eucl
# 
# ellipse_bry<-ggplot(inds2, aes(x=MDS1, y=MDS2, fill=mov_mod, colour=mov_mod)) +
#   geom_jitter(size=1, width=0.05, alpha=0.5)+ theme_bw()+ 
#   #ylim(-0.1,0.25)+xlim(-0.25,0.25)+
#   stat_ellipse(type = "t", level= 0.5, linetype=1, lwd=1.2, segments=6) + 
#   ylab("MDS2") + xlab("MDS1")+
#   ggtitle("a) Ordinance return time index")
# ellipse_bry

mind <- com_dt[,c(12:21)]
# mind$mn_rtp <- NULL
# mind$sd_rtp <- NULL
mind$mn_rrt <- NULL
mind$sd_rrt <- NULL
mind$mn_rg <- 1 / mind$mn_rg
mind$sd_rg <- 1 / mind$sd_rg
for(i in 1:(ncol(mind)-1)){mind[,i] <- ((scale(mind[,i])))}

minds <- data.frame(multi = apply(mind[,-ncol(mind)], 1, mean),
                    mov_mod = mind$mov_mod)
head(minds)

# ggplot(minds, aes(x = mov_mod, y = (multi),  fill = mov_mod))+
#   geom_boxplot(outlier.shape = NA)+ #ylim(0, 1.5) + 
#   theme_bw()
# 
# ggplot(minds, aes(x = mov_mod, y = 1/exp(multi),  fill = mov_mod))+
#   geom_boxplot(outlier.shape = NA)+ ylim(0, 3) + theme_bw()


com_dt$mov_mod <- as.factor(com_dt$mov_mod)
levels(com_dt$mov_mod)

minds$mov_mod <- relevel(minds$mov_mod, ref="null")
md <- lm((1/exp(multi)) ~ mov_mod, data = minds)
summary(lm((1/exp(multi)) ~ mov_mod, data = minds))
summary(lm((multi) ~ mov_mod, data = minds))
confint(md)
md
