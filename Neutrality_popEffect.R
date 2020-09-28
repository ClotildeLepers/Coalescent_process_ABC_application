##########################################################################################################
### 2017/07/25
### Analysis of the neutrality of simulations: step 1-3
###   1) abc on 10 pseudo-data per number of sub-populations
###   2) selection of 25 parameter sets for each abc run (in the posterior, depending on probability)
###   3) export the parameter sets and launch 500 simulations per set
###   4) analysis of neutrality of each pack of 500 simulations 
##########################################################################################################





# install.packages("abc")
library(abc)

##########################################################################################################
###---- Load simulated data: 
##########################################################################################################

###------ Loading of the simulated data
load("simuls.RData")

abcSaving=c()

##########################################################################################################
###---- Select post-ABC simulations: 1 sub-pop (meanFST = varFST = NA...)
##########################################################################################################

nbS = 1
#---- Select simulations with a given 2 sub-populations
simulsNbS = c()
for ( row in 2:nrow(simuls)) {
  if (simuls$nbSubPop[row] == nbS & is.na(simuls$nbSubPop[row]) == FALSE) {
    simulsNbS = c(simulsNbS,row)
  }
}
#---- Step 1: Choose 10 of those simulations and process ABC analysis using them as speudo-data
for(dataNb in 1:10){
  iObs = sample(simulsNbS,1,replace = FALSE)
  obs <-  simuls[iObs,]
  simulsABC <- simuls[c(1:(iObs-1),(iObs+1):nrow(simuls)),]
  
  #---- abc on the simulation => NOT TAKING INTO ACCOUNT: nbCherries, tMRCA, branch length
  abcRes=abc(data.frame(c(obs[24],obs[26:27],obs[30:34])),
             data.frame(c(simulsABC[1:7])),
             data.frame(c(simulsABC[24],simulsABC[26:27],simulsABC[30:34])),
             tol=0.01,method="neuralnet")
  
  #---- Step 2: keep 25 parameter sets from the posterior distribution
  nbSimulsKept = length(abcRes$weights)
  echParam = sample(nbSimulsKept,size=25,replace=TRUE, prob=abcRes$weights)
  
  #---- Step 3: Store the parameters of the selected simulations => used to run 500 simulations with those parameters
  params=cbind(abcRes$adj.values[echParam,],rep(nbS,length(echParam)))
  write.table(params, file="1trait_2017_07_25_automaticLauncher/parametersPostABC_pop1-5.txt",append = TRUE,row.names=F,col.names=F)
  abcSaving=append(abcSaving, abcRes)
  
}

##########################################################################################################
###---- Select post-ABC simulations: 2 sub-pop (varFST = NA...) => not taking into account pop neutral
##########################################################################################################

nbS = 2
#---- Select simulations with a given 2 sub-populations
simulsNbS = c()
for ( row in 2:nrow(simuls)) {
  if (simuls$nbSubPop[row] == nbS & is.na(simuls$nbSubPop[row]) == FALSE) {
    simulsNbS = c(simulsNbS,row)
  }
}

#---- Step 1: Choose 10 of those simulations and process ABC analysis using them as speudo-data
for(dataNb in 1:10){
  iObs = sample(simulsNbS,1,replace = FALSE)
  obs <-  simuls[iObs,]
  simulsABC <- simuls[c(1:(iObs-1),(iObs+1):nrow(simuls)),]
  
  #---- abc on the simulation => NOT TAKING INTO ACCOUNT: nbCherries, tMRCA, branch length
  # abcRes=abc(data.frame(c(obs[24],obs[26:27],obs[30:32],obs[c(33,35,37,39,41,43,45,47,49)])),
  #            data.frame(c(simulsABC[1:7])),
  #            data.frame(c(simulsABC[24],simulsABC[26:27],simulsABC[30:32],simulsABC[c(33,35,37,39,41,43,45,47,49)])),
  #            tol=0.01,method="neuralnet")
  abcRes=abc(data.frame(c(obs[24],obs[26:27],obs[30:34])),
             data.frame(c(simulsABC[1:7])),
             data.frame(c(simulsABC[24],simulsABC[26:27],simulsABC[30:34])),
             tol=0.01,method="neuralnet")
  
  #---- Step 2: keep 25 parameter sets from the posterior distribution
  nbSimulsKept = length(abcRes$weights)
  echParam = sample(nbSimulsKept,size=25,replace=TRUE, prob=abcRes$weights)
  
  #---- Step 3: Store the parameters of the selected simulations => used to run 500 simulations with those parameters
  params=cbind(abcRes$adj.values[echParam,],rep(nbS,length(echParam)))
  write.table(params, file="1trait_2017_07_25_automaticLauncher/parametersPostABC_pop1-5.txt",append = TRUE,row.names=F,col.names=F)
  abcSaving=append(abcSaving, abcRes)
  
}



##########################################################################################################
###---- Select post-ABC simulations: > 2 sub-pop (varFST != NA...)
##########################################################################################################

for(nbS in 3:5){
  #---- Select simulations with a given number of sub-populations (3 to 5)
  simulsNbS = c()
  for ( row in 2:nrow(simuls)) {
    if (simuls$nbSubPop[row] == nbS & is.na(simuls$nbSubPop[row]) == FALSE) {
      simulsNbS = c(simulsNbS,row)
    }
  }
  #---- Step 1: Choose 10 of those simulations and process ABC analysis using them as speudo-data
  for(dataNb in 1:10){
    iObs = sample(simulsNbS,1,replace = FALSE)
    obs <-  simuls[iObs,]
    simulsABC <- simuls[c(1:(iObs-1),(iObs+1):nrow(simuls)),]
    
    #---- abc on the simulation => NOT TAKING INTO ACCOUNT: nbCherries, tMRCA, branch length
    abcRes=abc(data.frame(c(obs[24],obs[26:27],obs[30:50])),
               data.frame(c(simulsABC[1:7])),
               data.frame(c(simulsABC[24],simulsABC[26:27],simulsABC[30:50])),
               tol=0.01,method="neuralnet")
    
    #---- Step 2: keep 25 parameter sets from the posterior distribution
    nbSimulsKept = length(abcRes$weights)
    echParam = sample(nbSimulsKept,size=25,replace=TRUE, prob=abcRes$weights)
    
    #---- Step 3: Store the parameters of the selected simulations => used to run 500 simulations with those parameters
    params=cbind(abcRes$adj.values[echParam,],rep(nbS,length(echParam)))
    write.table(params, file="1trait_2017_07_25_automaticLauncher/parametersPostABC_pop1-5.txt",append = TRUE,row.names=F,col.names=F)
    abcSaving=append(abcSaving, abcRes)
  }

}

#---- Problem of encoding when reading the file in calc => change the separator + add lineNb
param=read.table(file="1trait_2017_07_25_automaticLauncher/parametersPostABC_pop1-5.csv",sep=",",fileEncoding = "")
write.table(param, file="1trait_2017_07_25_automaticLauncher/parametersPostABC_pop1-5.txt",append = FALSE,row.names=F,col.names=F)



##########################################################################################################
###     Step 4: Analysis of packs of 500 simulations made with the same parameter set
###        Comparison of the distribution of the statistics of the 500 simulations
###              vs the distribution of the statistics of 5000 Kingman trees
##########################################################################################################


##########################################################################################################
###---- Load simulated data: 
##########################################################################################################

simulsPostABC <- as.vector(read.table("1trait_2017_07_25_automaticLauncher/selectStat_postABC_pop1-5.csv", header=T,sep = ",",fill=TRUE))
statKingCoeff <-as.vector(read.table("Kingman_coalescent/kingman_stat_coeff", header=T,sep = ",",fill=TRUE))



##########################################################################################################
###---- Distribution of statistics under Kingman's coalescents
##########################################################################################################

#---- number of cherries
sampleSize=1000
CnKing=(statKingCoeff$nbCherries - sampleSize /3)/sqrt(2*sampleSize/45)


#---- branch length
meanExtBranchKing = weighted.mean(statKingCoeff$extBranch)
varExtBranchKing = weighted.mean((statKingCoeff$extBranch)*(statKingCoeff$extBranch)) - (meanExtBranchKing*meanExtBranchKing)
LnKing=(statKingCoeff$extBranch-meanExtBranchKing)/(sqrt(0.5*varExtBranchKing))

#---- tMRCA is calculated for each parameter set (nx/bx vary)


##########################################################################################################
###---- Analysis of packs of 500 simulations made with the same parameter set: all simulations
###    Comparison of the distribution of the statistics of the 500 simulations
###          vs the distribution of the statistics of 5000 Kingman trees
##########################################################################################################

sampleSize = 1000


#---- initialization of tables
LnKS = c()
LnWilcox = c()
CnKS=c()
CnWilcox=c()
tmrcaWilcox=c()
tmrcaKS=c()

LnKSstat = c()
LnWilcoxStat = c()
CnKSstat=c()
CnWilcoxStat=c()
tmrcaWilcoxStat=c()
tmrcaKSstat=c()


for(paramSet in unique(simulsPop$lineNb)){
  #---- Selection of simulations with the right parameter set
  simulsParam = simulsPop[simulsPop$lineNb == paramSet,]
  
  
  #---- External Branch length
  try({
    meanExtBranch = weighted.mean(simulsParam$sumCoalExternalBranchLength)
    varExtBranch = weighted.mean((simulsParam$sumCoalExternalBranchLength)*(simulsParam$sumCoalExternalBranchLength)) - (meanExtBranch*meanExtBranch)
    
    LnSimuls=(simulsParam$sumCoalExternalBranchLength-meanExtBranch)/(sqrt(0.5*varExtBranch))
    
    LnKS = c(LnKS,ks.test(LnSimuls, LnKing)$p.value)
    LnWilcox = c(LnWilcox, wilcox.test(LnSimuls, LnKing,paired=FALSE,"two.sided")$p.value)
    
    LnKSstat = c(LnKSstat,ks.test(LnSimuls, LnKing)$statistic)
    LnWilcoxStat = c(LnWilcoxStat, wilcox.test(LnSimuls, LnKing,paired=FALSE,"two.sided")$statistic)
  })
  
  
  
  #---- Number of cherries
  try({
    CnSimuls=(simulsParam$nbCherries - sampleSize /3)/sqrt(2*sampleSize/45)
    
    CnKS=c(CnKS, ks.test(CnSimuls, CnKing)$p.value)
    CnWilcox=c(CnWilcox,wilcox.test(CnSimuls, CnKing,paired=FALSE,"two.sided")$p.value)
    
    CnKSstat=c(CnKSstat, ks.test(CnSimuls, CnKing)$statistic)
    CnWilcoxStat=c(CnWilcoxStat,wilcox.test(CnSimuls, CnKing,paired=FALSE,"two.sided")$statistic)
    
  })
  
  
  #---- tMRCA
  try({
    bx= exp(-(simulsParam$meanXind/simulsParam$sigB)^2)/2
    nx=(bx/(simulsParam$etaC * exp(-(0/simulsParam$sigC)^2)/2))
    nbSimuls = length(nx)
    tmrcaKing=rep(0,nbSimuls)
    for(i in 2:sampleSize)
    {
      tmrcaKing=tmrcaKing+rexp(nbSimuls,rate=(i*(i-1)/2))
    }
    tmrcaKing = tmrcaKing*nx/bx
    
    tmrcaKing[tmrcaKing>mean(simulsParam$tsim)] = mean(simulsParam$tsim)
    
    tmrcaKS=c(tmrcaKS, ks.test(tmrcaKing,simulsParam$tMRCA)$p.value)
    tmrcaWilcox=c(tmrcaWilcox,wilcox.test(tmrcaKing,simulsParam$tMRCA,paired=FALSE,"two.sided")$p.value)
    
    
    tmrcaKSstat=c(tmrcaKSstat, ks.test(tmrcaKing,simulsParam$tMRCA)$statistic)
    tmrcaWilcoxStat=c(tmrcaWilcoxStat,wilcox.test(tmrcaKing,simulsParam$tMRCA,paired=FALSE,"two.sided")$statistic)
  })

  
}


par(mfrow = c(3,2))

boxplot(LnKS,names=c(""),#
        main="normalized external branch length",
        xlab="number of sub-population in the reference simulation",
        ylab="p-value of Kolmogorov-Smirnov tests",
        ylim=c(0,1))
abline(h=0.05, col='red',lty=3)


boxplot(LnWilcox,names=c(""),#
        main="normalized external branch length",
        xlab="number of sub-population in the reference simulation",
        ylab="p-value of Mann-Whitney tests",
        ylim=c(0,1))
abline(h=0.05, col='red',lty=3)

boxplot(CnKS,names=c(""),#
        main="normalized number of cherries",
        xlab="number of sub-population in the reference simulation",
        ylab="p-value of Kolmogorov-Smirnov tests",
        ylim=c(0,1))
abline(h=0.05, col='red',lty=3)

boxplot(CnWilcox,names=c(""),#
        main="normalized number of cherries",
        xlab="number of sub-population in the reference simulation",
        ylab="p-value of the Mann-Whitney tests",
        ylim=c(0,1))
abline(h=0.05, col='red',lty=3)

boxplot(tmrcaKS,names=c(""),#
        main="time before the MRCA",
        xlab="number of sub-population in the reference simulation",
        ylab="p-value of Kolmogorov-Smirnov tests",
        ylim=c(0,1))
abline(h=0.05, col='red',lty=3)

boxplot(tmrcaWilcox,names=c(""),#
        main="time before the MRCA",
        xlab="all simulations",
        ylab="p-value of the Mann-Whitney tests",
        ylim=c(0,1))
abline(h=0.05, col='red',lty=3)






##########################################################################################################
###---- Analysis of packs of 500 simulations made with the same parameter set: as a function of sub-pop nb
###           1. Determine the mean/mode number of sub-populations
###                 of simulations made with the parameter set
##########################################################################################################


#---- Add columns with the mean and mode number of subpopulations of all simulations with the same parameters
simulsPostABCmodif<-simulsPostABC
simulsPostABCmodif<-cbind(simulsPostABCmodif,rep(0, nrow(simulsPostABC)))
simulsPostABCmodif<-cbind(simulsPostABCmodif,rep(0, nrow(simulsPostABC)))
simulsPostABCmodif<-cbind(simulsPostABCmodif,rep(0, nrow(simulsPostABC)))
simulsPostABCmodif<-cbind(simulsPostABCmodif,rep(0, nrow(simulsPostABC)))
names(simulsPostABCmodif)[22] = "meanPopSimuls"
names(simulsPostABCmodif)[23] = "modePopSimuls"
names(simulsPostABCmodif)[24] ="nbSimulsParam"
names(simulsPostABCmodif)[25] ="maxPopSimuls"

nbSimuls=c()
for(paramSet in unique(simulsPostABCmodif$lineNb)){
  #---- Selection of simulations with the right parameter set
  simulsParam = simulsPostABCmodif[simulsPostABCmodif$lineNb == paramSet,]
  meanPopNb = mean(simulsParam$nbSubPop)
  modePopNb = mean(as.numeric(names(which.max(table(simulsParam$nbSubPop)))))
  maxPopNb = max(simulsParam$nbSubPop)
  nbSimulsParam = nrow(simulsParam)
  nbSimuls = cbind(nbSimuls, nbSimulsParam)
  try({
    simulsPostABCmodif[simulsPostABCmodif$lineNb == paramSet,]$meanPopSimuls = meanPopNb
    simulsPostABCmodif[simulsPostABCmodif$lineNb == paramSet,]$modePopSimuls = modePopNb
    simulsPostABCmodif[simulsPostABCmodif$lineNb == paramSet,]$nbSimulsParam = nbSimulsParam
    simulsPostABCmodif[simulsPostABCmodif$lineNb == paramSet,]$maxPopSimuls = maxPopNb
  })
}
summary(simulsPostABCmodif)

table(nbSimuls)
ncol(nbSimuls)
table(param$V8)

plot(modeParamPopNb~nbSubPopInit)

mean(varParamPopNb,na.rm=TRUE)

par(mfrow = c(1,1))
plot(simulsPostABCmodif$modePopSimuls~simulsPostABCmodif$nbSubPopData)





##########################################################################################################
###---- Analysis of packs of 500 simulations made with the same parameter set: as a function of sub-pop nb
###       2. Compare the distributions of the statistics (Ln, Cn, tMRCA)
###            500 simuls with same parameters vs 5000 Kingman trees
##########################################################################################################

table(simulsPostABCmodif$modePopSimuls)
summary(simulsPostABCmodif$meanPopSimuls)
summary(simulsPostABCmodif$modePopSimuls)
length(unique(simulsPostABCmodif[(simulsPostABCmodif$modePopSimuls == 10),]$lineNb))

nbPopMin = 1
nbPopMax = 4

branchKSpop4 = c()
branchWilcoxPop4 =c()

cherriesKSpop4=c()
cherriesWilcoxPop4=c()

tmrcaKSpop4=c()
tmrcaWilcoxPop4=c()
tmrcaKSpopTronc4=c()
tmrcaWilcoxPopTronc4=c()

sampleSize = 1000

for(nbPop in nbPopMin:nbPopMax){
  
  #---- Selection of the simulations with the right number of populations
  simulsPop = simulsPostABCmodif[(simulsPostABCmodif$modePopSimuls == nbPop),]
  if(nbPop == 4){
    simulsPostABCmodif[(simulsPostABCmodif$modePopSimuls >= nbPop),]
  }
  
  #---- initialization of tables
  LnKS = c()
  LnWilcox = c()
  CnKS=c()
  CnWilcox=c()
  tmrcaWilcox=c()
  tmrcaKS=c()
  tmrcaWilcoxTronc=c()
  tmrcaKSTronc=c()
  nbSimulsNbPop = 0
  
  #---- select simulations with the same parameters
  for(paramSet in unique(simulsPop$lineNb)){
    #---- Selection of simulations with the right parameter set
    simulsParam = simulsPop[simulsPop$lineNb == paramSet,]
    
    # cat(summary(simulsParam))
    
    #---- Verify that the number of simulations with the same parameter set is high enough
    if(mean(simulsParam$nbSimulsParam) > 300){
      nbSimulsNbPop = nbSimulsNbPop +1
      #---- External Branch length
      try({
        meanExtBranch = weighted.mean(simulsParam$sumCoalExternalBranchLength)
        varExtBranch = weighted.mean((simulsParam$sumCoalExternalBranchLength)*(simulsParam$sumCoalExternalBranchLength)) - (meanExtBranch*meanExtBranch)
        
        LnSimuls=(simulsParam$sumCoalExternalBranchLength-meanExtBranch)/(sqrt(0.5*varExtBranch))
        
        # qqplot(LnSimuls,LnKing)
        
        LnKS = c(LnKS,ks.test(LnSimuls, LnKing)$p.value)
        LnWilcox = c(LnWilcox, wilcox.test(LnSimuls, LnKing,paired=FALSE,"two.sided")$p.value)
      })
      
      
      
      #---- Number of cherries
      try({
        CnSimuls=(simulsParam$nbCherries - sampleSize /3)/sqrt(2*sampleSize/45)
        
        # qqplot(CnSimuls,CnKing)
        
        CnKS=c(CnKS, ks.test(CnSimuls, CnKing)$p.value)
        CnWilcox=c(CnWilcox,wilcox.test(CnSimuls, CnKing,paired=FALSE,"two.sided")$p.value)
      })
      
      
      #---- tMRCA
      try({
        bx= exp(-(simulsParam$meanXind/simulsParam$sigB)^2)/2
        nx=(bx/(simulsParam$etaC * exp(-(0/simulsParam$sigC)^2)/2))
        nbSimuls = length(nx)
        tmrcaKing=rep(0,nbSimuls)
        for(i in 2:sampleSize)
        {
          tmrcaKing=tmrcaKing+rexp(nbSimuls,rate=(i*(i-1)/2))
        }
        tmrcaKing = tmrcaKing*nx/bx
        
        tmrcaKingTronc = tmrcaKing
        tmrcaKingTronc[tmrcaKingTronc>mean(simulsParam$tsim)] = mean(simulsParam$tsim)
        
        # qqplot(simulsParam$tMRCA,tmrcaKing)
        # qqplot(simulsParam$tMRCA,tmrcaKingTronc)
        
        tmrcaKS=c(tmrcaKS, ks.test(tmrcaKing,simulsParam$tMRCA)$p.value)
        tmrcaWilcox=c(tmrcaWilcox,wilcox.test(tmrcaKing,simulsParam$tMRCA,paired=FALSE,"two.sided")$p.value)
        
        
        tmrcaKSTronc=c(tmrcaKSTronc, ks.test(tmrcaKingTronc,simulsParam$tMRCA)$p.value)
        tmrcaWilcoxTronc=c(tmrcaWilcoxTronc,wilcox.test(tmrcaKingTronc,simulsParam$tMRCA,paired=FALSE,"two.sided")$p.value)
      })      
      
      cat("param ",paramSet, "done ")
      
      
    }
    
    
  }
  
  cat("\n pop ",nbPop, "done.", nbSimulsNbPop ,"paramater sets")
  branchKSpop4 = cbind(branchKSpop4, LnKS)
  branchWilcoxPop4 = cbind(branchWilcoxPop4,LnWilcox)
  
  cherriesKSpop4 = cbind(cherriesKSpop4, CnKS)
  cherriesWilcoxPop4 = cbind(cherriesWilcoxPop4,CnWilcox)
  
  tmrcaKSpop4=cbind(tmrcaKSpop4,tmrcaKS)
  tmrcaWilcoxPop4=cbind(tmrcaWilcoxPop4,tmrcaWilcox)
  
  tmrcaKSpopTronc4=cbind(tmrcaKSpopTronc4,tmrcaKSTronc)
  tmrcaWilcoxPopTronc4=cbind(tmrcaWilcoxPopTronc4,tmrcaWilcoxTronc)
  
}



par(mfrow = c(4,2))

boxplot(branchKSpop4,names=c("1","2","3","4 or +"),
        main="Ln",
        xlab="mode number of sub-populations",
        ylab="p-value K-S tests",
        ylim=c(0,1))
abline(h=0.05, col='red',lty=3)

boxplot(branchWilcoxPop4,names=c("1","2","3","4 or +"),
        main="Ln",
        xlab="mode number of sub-populations",
        ylab="p-value M-W tests",
        ylim=c(0,1))
abline(h=0.05, col='red',lty=3)


boxplot(cherriesKSpop4,names=c("1","2","3","4 or +"),
        main="Cn",
        xlab="mode number of sub-populations",
        ylab="p-value K-S tests",
        ylim=c(0,1))
abline(h=0.05, col='red',lty=3)

boxplot(cherriesWilcoxPop4,names=c("1","2","3","4 or +"),
        main="Cn",
        xlab="mode number of sub-populations",
        ylab="p-value M-W tests",
        ylim=c(0,1))
abline(h=0.05, col='red',lty=3)



boxplot(tmrcaKSpop4,names=c("1","2","3","4 or +"),
        main="tMRCA",
        xlab="mode number of sub-populations",
        ylab="p-value K-S",
        ylim=c(0,1))
abline(h=0.05, col='red',lty=3)

boxplot(tmrcaWilcoxPop4,names=c("1","2","3","4 or +"),
        main="tMRCA",
        xlab="mode number of sub-populations",
        ylab="p-value M-W tests",
        ylim=c(0,1))
abline(h=0.05, col='red',lty=3)



boxplot(tmrcaKSpopTronc4,names=c("1","2","3","4 or +"),
        main="trunctaed tMRCA",
        xlab="mode number of sub-populations",
        ylab="p-value K-S tests",
        ylim=c(0,1))
abline(h=0.05, col='red',lty=3)

boxplot(tmrcaWilcoxPopTronc4,names=c("1","2","3","4 or +"),
        main="truncated tMRCA",
        xlab="mode number of sub-populations",
        ylab="p-value M-W tests",
        ylim=c(0,1))
abline(h=0.05, col='red',lty=3)





##########################################################################################################
### qq plot
##########################################################################################################

par(mfrow = c(2,2))

simulsPop = simulsPostABCmodif[(simulsPostABCmodif$modePopSimuls == 5),]
iObs = sample(unique(simulsPop$lineNb),1,replace = FALSE)
simulsParam = simulsPop[simulsPop$lineNb == iObs,]

#---- External Branch length
meanExtBranch = weighted.mean(simulsParam$sumCoalExternalBranchLength)
varExtBranch = weighted.mean((simulsParam$sumCoalExternalBranchLength)*(simulsParam$sumCoalExternalBranchLength)) - (meanExtBranch*meanExtBranch)

LnSimuls=(simulsParam$sumCoalExternalBranchLength-meanExtBranch)/(sqrt(0.5*varExtBranch))

qqplot(LnSimuls,LnKing)

ks.test(LnSimuls, LnKing)
wilcox.test(LnSimuls, LnKing,paired=FALSE,"two.sided")

#---- Number of cherries
CnSimuls=(simulsParam$nbCherries - sampleSize /3)/sqrt(2*sampleSize/45)

qqplot(CnSimuls,CnKing)

ks.test(CnSimuls, CnKing)
wilcox.test(CnSimuls, CnKing,paired=FALSE,"two.sided")

#---- tMRCA
bx= exp(-(simulsParam$meanXind/simulsParam$sigB)^2)/2
nx=(bx/(simulsParam$etaC * exp(-(0/simulsParam$sigC)^2)/2))
nbSimuls = length(nx)
tmrcaKing=rep(0,nbSimuls)
for(i in 2:sampleSize)
{
  tmrcaKing=tmrcaKing+rexp(nbSimuls,rate=(i*(i-1)/2))
}
tmrcaKing = tmrcaKing*nx/bx

tmrcaKingTronc = tmrcaKing
tmrcaKingTronc[tmrcaKingTronc>mean(simulsParam$tsim)] = mean(simulsParam$tsim)

qqplot(simulsParam$tMRCA,tmrcaKing)
qqplot(simulsParam$tMRCA,tmrcaKingTronc)

ks.test(tmrcaKing,simulsParam$tMRCA)
wilcox.test(tmrcaKing,simulsParam$tMRCA,paired=FALSE,"two.sided")


ks.test(tmrcaKingTronc,simulsParam$tMRCA)
wilcox.test(tmrcaKingTronc,simulsParam$tMRCA,paired=FALSE,"two.sided")
