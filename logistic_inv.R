logistic_inv <- function(ndvi){
  source('global_veg/movement models/cc_score.R')
  time <- 1:length(ndvi)
  
  Wx <- function(x){
    erg <- sum( ((x[1] + (x[2]/(1 + exp((x[3]-time)/x[4])))) - (ndvi))^2)
    return(ifelse((is.infinite(erg)||is.nan(erg)),1e50,erg))
  }
  
  Wt <- function(Init, Asym, xmid, scal) {
    erg <-  Init + (Asym/(1 + exp((xmid-time)/scal)))
    return(erg)
  }
  
  Init <- 0.1
  Asym <- 0.1
  xmid <- 6
  scal <- -0.9
  
  xmids <- seq(2, 33, by= 3)
  par_dat <- matrix(data=NA, nrow=length(xmids), ncol=6)
  for(z in 1:length(xmids)){
    optimal1 <- try(optim(c(Init, Asym, xmid, scal),fn=Wx,gr=NULL,method="L-BFGS-B",
                          lower=c(0.001,0, xmids[z], -7), upper=c(1, 1,34, -0.1),
                          control=list(maxit = 5000, pgtol = 1e-10, 
                                       ndeps = c(1e-10, 1e-10, 1e-10, 1e-10), 
                                       lmm = 200)))
    Init1 <- optimal1$par[1]
    Asym1 <- optimal1$par[2]
    xmid1 <- optimal1$par[3]
    scal1 <- optimal1$par[4]
    cc1 <- cc_score(ndvi, Wt(Init1, Asym1, xmid1, scal1))
    par_dat[z,] <- c(z, cc1, c(Init1, Asym1, xmid1, scal1))
  }
  
  id1 <- match(max(par_dat[,2]), par_dat[,2])
  par <- par_dat[id1, 3:6]
  
  fit <- Wt(par[1], par[2], par[3], par[4])
  return(fit)
}
