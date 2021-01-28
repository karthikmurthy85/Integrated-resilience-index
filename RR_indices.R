rg_fun <- function(df){
  sl <- NA
  for(i in 1:(nrow(df)-1)){
    sl[i] <- df$y[i+1] - df$y[i]
  }
  return(mean(sl))
}


rr_fun <- function(cycdt, key, lindt){
  ifelse(key == -1, k1 <- parse(text="max(nd)"), 
         k1 <- parse(text="min(nd)"))
  
  p1 <- subset(cycdt, nd %in% eval(k1))
  lin1_x <- cycdt$t[1:p1$id]
  lin1_y <- cycdt$nd[1:p1$id]
  curv1 <- data.frame(x = lin1_x, y = lin1_y)
  g1 <- NA
  try(g1 <- curve_intersect(curv1, lindt,  empirical = TRUE))
  ifelse(length(g1) == 2, ints1 <- g1$x, ints1 <- NA)
  
  #print("intersection 1 done")
  lin2_x <- cycdt$t[p1$id:nrow(cycdt)]
  lin2_y <- cycdt$nd[p1$id:nrow(cycdt)]
  curv2 <- data.frame(x = lin2_x,y = lin2_y)
  g2 <- NA
  try(g2 <- curve_intersect(curv2, lindt,  empirical = TRUE))
  ifelse(length(g2) == 2, ints2 <- g2$x, ints2 <- NA)
  
  #print("intersection 2 done")
  
  rg1 <- rg_fun(curv1); rg2 <- rg_fun(curv2)
 
  rg <- sum(c(abs(rg1), abs(rg2)), na.rm=T)
  rr <- sqrt((cycdt$nd[1] - lindt$y[1])^2)
  
  rt1 <- c(ints1, ints2)
  
  return(c(rg, rr, rt1))
}


rrg_fun <- function(xtr_y, amn){
  dt <- NA
  for(i in 1:length(xtr_y)){
    y1 <- xtr_y[i]
    y2 <- amn
    dt[i] <- sqrt((y1-y2)^2)
  }
  return(dt)
}

pk_dif <- function(xtrms){
  tdif <- NA
  for(i in 1:(length(xtrms)-1)){
    tdif[i] <- xtrms[i+1] - xtrms[i]
  }
  return(tdif)
}

