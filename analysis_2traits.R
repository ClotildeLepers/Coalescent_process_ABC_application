##########################################################################################################
###             Analysis of simulations with 2 traits: location and reproductive strategy              ###
##########################################################################################################

# install.packages("abc")
library(abc)

##########################################################################################################
###---- Load simulated data + check the data.
##########################################################################################################
setwd("C://Users//Clotilde//Dropbox//Lille 1//Analyse//Resultats_2traits")
simulsEqualb_0 <-as.vector(read.table("simulations_4_e.csv", header=T,check.names = TRUE,sep = ";",fill=TRUE))
simulsEqualb_0 = simulsEqualb_0[,1:69]
simulsEqualb_0$pxb10=0
simulsEqualb_0=simulsEqualb_0[,c(1,69,2:68)]
simulsDiffb_0 <-as.vector(read.table("simulations_4_d.csv", header=T,check.names = TRUE,sep = ";",fill=TRUE))
simulsDiffb_0 = simulsDiffb_0[,1:68]
simulsDiffb_0$pxb10=0
simulsDiffb_0=simulsDiffb_0[,c(1,69,2:68)]


simulsEqualb$diffb=0
simulsDiffb$diffb=1

simuls_2x=rbind(simulsDiffb[1:10000,],simulsEqualb[1:10000,])


##########################################################################################################
###---- posterior distribution for each abc
##########################################################################################################

###------ Choosing data from the simulations
###---- with choice of the number of sub-populations

simuls <- simuls_2x

summary(simuls$nbSubPop)
nbS = 10
diffb = 0
simulsNbS = c()
for ( row in 1:nrow(simuls)) {
  if (simuls$nbSubPop[row] >= nbS & is.na(simuls$nbSubPop[row]) == FALSE & (simuls$diffb[row]== diffb) & (is.na(simuls$varLocIntraX[row]) == FALSE & simuls$meanXPop[row] >0.1 & simuls$meanXPop [row]<0.9 )) {
    simulsNbS = c(simulsNbS,row)
  }
}
iObs = sample(simulsNbS,100,replace = FALSE)

nbS = 10
diffb = 1
simulsNbS = c()
for ( row in 1:nrow(simuls)) {
  if (simuls$nbSubPop[row] >= nbS & is.na(simuls$nbSubPop[row]) == FALSE & (simuls$diffb[row]== diffb) & (is.na(simuls$varLocIntraX[row]) == FALSE & simuls$meanXPop[row] >0.1 & simuls$meanXPop [row]<0.9 )) {
    simulsNbS = c(simulsNbS,row)
  }
}

iObs = c(iObs,sample(simulsNbS,100,replace = FALSE))

posteriorWeights = simuls_2x[,c(1:11,70)]
names(posteriorWeights)

###---- For each pseudo-data : abc + keep weight of simuls
for(obsNb in 1:length(iObs)){
  simuls <- simuls_2x
  ###---- delete the data row from simulations
  obs <-  simuls[iObs[obsNb],]
  simuls <- simuls[c(1:(iObs[obsNb]-1),(iObs[obsNb]+1):nrow(simuls)),]
  
  ###---- abcRes_3 : statistics obtained on data
  abcRes_Data3=abc(data.frame(c(obs[c(25,27,28,32,35:55,57,59,61,63)])),
                   data.frame(c(simuls[c(1:11)])),
                   data.frame(c(simuls[c(25,27,28,32,35:55,57,59,61,63)])),
                   tol=0.01,method="neuralnet")
  
  ###---- save posterior weights
  weights = abcRes_Data3$region
  weights[weights==TRUE]=abcRes_Data3$weights/sum(abcRes_Data3$weights)
  weights =c(weights[1:(iObs[obsNb]-1)],NA,weights[(iObs[obsNb]):length(weights)])
  
  posteriorWeights = cbind(posteriorWeights,weights)
  
}



simuls <- simuls_2x

HeyerObs<-as.vector(read.table("variables_Heyer.csv", header=T,check.names = TRUE,sep = ";",fill=TRUE))

abcRes_Data3=abc(data.frame(c(HeyerObs[c(25,27,28,32,35:55,57,59,61,63)])),
                 data.frame(c(simuls[c(1:11)])),
                 data.frame(c(simuls[c(25,27,28,32,35:55,57,59,61,63)])),
                 tol=0.05,method="neuralnet")

weightsHeyer = abcRes_Data3$region
weightsHeyer[weightsHeyer==TRUE]=abcRes_Data3$weights/sum(abcRes_Data3$weights)

posteriorWeights = cbind(posteriorWeights,weightsHeyer)


par(mfrow = c(2,6))
names(simuls)
names(simuls)[11]="q"

par(mfrow = c(2,6))
for(param in 1:11){
  ymax=max(
    density(simuls[,param],na.rm=TRUE)$y,
    density(abcRes_Data3$adj.values[,param])$y
  )
  xmin=max(0,min(
    quantile(simuls[,param],na.rm=TRUE,probs=0.05),
    quantile(abcRes_Data3$adj.values[,param],probs=0.05))
  )
  xmax=max(
    quantile(simuls[,param],na.rm=TRUE,probs=0.95),
    quantile(abcRes_Data3$adj.values[,param],probs=0.95)
  )
  plot(density(simuls[,param],na.rm=TRUE),col="black", lty=2, xlab="", main=names(simuls)[param], ylim=c(0,ymax+ymax/10),xlim=c(xmin-xmin/10,xmax+xmax/10))
  lines(density(abcRes_Data3$adj.values[,param]),col="red")
}

