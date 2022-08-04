##Inverse double logistic function - parameter k is negative
dLogistic_inv <- function(ndvi){
  
  time <- 1:length(ndvi)
  
  Wx <- function(x){
    erg <- sum(((x[1] + (x[3]/(1+exp(-x[6]*(time-x[4])))) - 
                   ((x[3]+x[1]-x[2])/(1+exp(-x[5]*(time-x[7]))))) - (ndvi))^2)
    return(ifelse((is.infinite(erg)||is.nan(erg)),1e50,erg))
  }
  
  Wt <- function(vb,ve,k,p,d,c,q) {
    erg <- vb + (k/(1+exp(-c*(time-p)))) - ((k+vb-ve)/(1+exp(-d*(time-q))))
    return(erg)
  }
  
  pk <- match(min(ndvi), ndvi)
  
  vb <- ndvi[1]
  ve <- ndvi[length(ndvi)]
  c <- d <- 0.0001
  k <- min(ndvi)- ndvi[1]
  p <- pk-(pk/2)
  q <- pk +(pk/2)
  ifelse(q > length(ndvi), q <- pk+1,  q <- q )
  ifelse(p <0, p <- 2,  p<-p )
  
  lowc <- lowd <- c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1 )
  
  fitdat <- data.frame(matrix(NA, ncol = 2, nrow = length(lowc)))
  
  for(i in 1:length(lowc)){
    
    optimal <- optim(c(vb,ve,k,p,d,c,q),fn=Wx,gr=NULL,method="L-BFGS-B",
                     lower=c(0, (0), (min(ndvi)- ndvi[1]- 0.1), 2, lowc[i], lowd[i], (pk+1)),
                     upper=c((vb + 0.001), (ve + 0.002), 0, (pk-5), 0.3, 0.3, 
                             (length(ndvi)-2)),
                     control=list(maxit = 5000, pgtol = 1e-10, 
                                  ndeps = c(1e-10, 1e-10, 1e-10, 1, 1e-10, 1e-10,1), 
                                  lmm = 200))
    
    vb1 <- optimal$par[1]
    ve1 <- optimal$par[2]
    k1 <- optimal$par[3]
    p1 <- optimal$par[4]
    d1 <- optimal$par[5]
    c1 <- optimal$par[6]
    q1 <- optimal$par[7]
    
    fit <- Wt(vb1,ve1,k1,p1,d1,c1,q1)
    plot(ndvi, ty='l')
    points(fit, ty='l')
    
    fitdat[i, ] <- c(lowc[i], cor(fit, ndvi))
    
  }
  
  fitdat1 <- subset(fitdat, X2 %in% max(X2))
  
  optimal <- optim(c(vb,ve,k,p,d,c,q),fn=Wx,gr=NULL,method="L-BFGS-B",
                   lower=c(0, (0), (min(ndvi)- ndvi[1]- 0.1), 2, fitdat1$X1[1], fitdat1$X1[1], (pk+1)),
                   upper=c((vb + 0.001), (ve + 0.002), 0, (pk-5), 0.3, 0.3, 
                           (length(ndvi)-2)),
                   control=list(maxit = 5000, pgtol = 1e-10, 
                                ndeps = c(1e-10, 1e-10, 1e-10, 1, 1e-10, 1e-10,1), 
                                lmm = 200))
  
  vb1 <- optimal$par[1]
  ve1 <- optimal$par[2]
  k1 <- optimal$par[3]
  p1 <- optimal$par[4]
  d1 <- optimal$par[5]
  c1 <- optimal$par[6]
  q1 <- optimal$par[7]
  
  fit <- Wt(vb1,ve1,k1,p1,d1,c1,q1)
  plot(ndvi, ty='l')
  points(fit, ty='l')
  return(fit)
  
}
