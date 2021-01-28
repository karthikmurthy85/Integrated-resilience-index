cc_score <- function(dist, fitted){
  
  num<-sum((dist-fitted)^2)
  denom1<-sum((dist-mean(fitted))^2)
  denom2<-sum((fitted-mean(fitted))^2)
  denom3<-length(dist)*(mean(dist)-mean(fitted))^2
  CC_score<- 1-((num)/(denom1+denom2+denom3))
  return(CC_score)
}


##Input raw data, fitted data and number of parameters
aic_calc <- function(raw, residuals, parameters){
  aic <- length(raw)*(log(2*pi)+1+log((sum(residuals^2)/length(raw)))) +((length(parameters)+1)*2)
  return(aic)
}


##R-adj 
R_adj_fun <- function(raw, pred, k){
  lm1 <- summary(lm(raw~pred))
  Rsq <- lm1$r.squared
  n = length(raw)
  Radj <- 1 - ((n-1)/(n-k-1))*(1-Rsq)
  return(Radj)
}

