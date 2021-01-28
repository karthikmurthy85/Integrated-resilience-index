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



