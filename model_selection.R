model_sel_fun <- function(model_dat){
  
  fita <- subset(model_dat, cc_score > (max(cc_score) - 0.075))
  if(nrow(fita)==1){
    fitm <- fita
  }
  
  else{
    fitb <- fita[order(fita$cc_score, decreasing = T),]
    fitb1 <- fitb[order(fitb$rank),]
    
    ifelse(fitb1$rank[1] == 1, fitm <- fitb1[1, ],
           ifelse(fitb1$rank[1] == fitb1$rank[2], fitm <- fitb1[1, ], 
                  ifelse(abs(fitb1$cc_score[1] - fitb1$cc_score[2]) < 0.03,
                         fitm <- c(max(fitb1$cc_score), 
                                   paste(fitb1$models[1], 
                                         fitb1$models[2], sep='_'), "unresolved"),
                         ifelse(fitb1$cc_score[1] - fitb1$cc_score[2] < -0.075,
                                fitm <- fitb1[2, ], fitm <- fitb1[1, ] ))))
  }
  
  return(fitm)
}



model_sel_fun2 <- function(x){
  x$models <- as.character(x$models)
  ##cc score and aic
  sel1 <- x[order(x$cc_score, decreasing = T),]
  sel2 <- sel1[1:2,]
  delcc <- sel2$cc_score[1] - sel2$cc_score[2]
  
  if(delcc > 0.02){selmod <- sel2$models[1]}
  
  else{
    if(delcc < 0.02){
      aicdat <- sel2[order(sel2$aic_score, decreasing = F),]
      delaic <- aicdat$aic_score[1] - aicdat$aic_score[2]
      selmod <- ifelse(delaic < -3, aicdat$models[1], 
                       paste(aicdat$models[1],aicdat$models[2], sep='-'))
    } 
  }
  return(selmod)
}


model_sel_fun3 <- function(x){
  x$models <- as.character(x$models)
  ##Rsq-adj
  sel1 <- x[order(x$radj_score, decreasing = T),]
  sel2 <- sel1[1:2,]
  
  Radj <- sel2$radj_score
  if(length(Radj[Radj < 0.05])==2){
    selmod <- 'null'
  }
  
  else{
    delRsq <- sel2$radj_score[1] - sel2$radj_score[2]
    
    selmod <- ifelse(delRsq > 0.01, sel2$models[1],
                     paste(sel2$models[1],sel2$models[2], sep='-'))
  }
  return(selmod)
}



