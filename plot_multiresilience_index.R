beta <- matrix(NA, ncol=4,nrow=100)
for(j in 1:100){
  source('global_veg/new_codes/ultiresilience_calc_ndvit.R')
  b1 <- melt(md$coefficients)
  b1$beta <- b1$value + b1$value[1]
  beta[j,] <- c(b1$value[1],b1$beta[2:4])
  print(j)
}
bmn <- apply(beta, 2, mean)
bsd <- apply(beta, 2, sd)

df <- data.frame(model = c('stable','abrupt','gradual','reversible'),
                 mean = bmn, sd = bsd)

library(ggplot2)
ggplot(df, aes(x = model, y = bmn, colour = model))+
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mean-sd, ymax= mean+sd), width=.2,
                position=position_dodge(.9))+
  scale_x_discrete(limits = c("stable", "reversible", "abrupt", "gradual"))+
  scale_color_manual(values=c("tan4","red4","royalblue4","green4"))+
  theme_bw() + xlab("") + ylab("Resilience index") + 
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.ticks = element_blank(),
    axis.title = element_text(size=16, colour = 'black'),
    axis.text = element_text(size=14, colour = 'black')
  )


##MeanNDVI time series##
beta <- matrix(NA, ncol=4,nrow=100)
for(j in 1:100){
  source('global_veg/new_codes/multiresilience_calc_meanNDVI.R')
  b1 <- melt(md$coefficients)
  b1$beta <- b1$value + b1$value[1]
  beta[j,] <- c(b1$value[1],b1$beta[2:4])
  print(j)
}
bmn <- apply(beta, 2, mean)
bsd <- apply(beta, 2, sd)

df <- data.frame(model = c('stable','abrupt','gradual','reversible'),
                 mean = bmn, sd = bsd)

ggplot(df, aes(x = model, y = bmn, colour = model))+
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mean-sd, ymax= mean+sd), width=.2,
                position=position_dodge(.9))+
  scale_x_discrete(limits = c("stable", "reversible", "abrupt", "gradual"))+
  scale_color_manual(values=c("tan4","red4","grey40","green4"))+
  theme_bw() + xlab("") + ylab("Resilience index") + 
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.ticks = element_blank(),
    axis.title = element_text(size=16, colour = 'black'),
    axis.text = element_text(size=14, colour = 'black')
  )


