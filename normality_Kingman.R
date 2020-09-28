##########################################################################################################
### 2017/07/25
### Analysis of Kingman's coalescents (simulated with package Kingman).
### The normalized number of cherries Cn is expected to follow a N(0,1)
### The normalized external branch length of the coalescent Ln is expected to follow a N(0,1)
###  The time before the MRCA is expected to follow the same distribution than 
###        Sum(i in 2:nbInd){exp(nbSimuls,rate=(coefficient * i*(i-1))) }
###
### From this analysis, we can see that Cn and Ln of Kingman coalescent do not follow N(0,1).
### Comparison of the simulations will be with the distribution obtained with Kingman's coalescent, not with N(0,1)
##########################################################################################################





###------ Loading of the simulated data of Kingman's coalescents
statKing <-as.vector(read.table("kingman_stat", header=T,sep = ",",fill=TRUE))
statKingCoeff <-as.vector(read.table("kingman_stat_coeff", header=T,sep = ",",fill=TRUE))
summary(statKing)






##########################################################################################################
#################                    Tests of normality                     ##############################
#################                      nb cherries Cn                       ##############################
##########################################################################################################



####----   All simulations
sampleSize=1000
ech=(statKing$nbCherries - sampleSize /3)/sqrt(2*sampleSize/45)
par(mfrow = c(1,1))
val=0.3
hist(ech,prob=TRUE,
     main="",
     xlab = "normalized number of cherries",
     breaks = seq(min(ech)-val,max(ech)+val,val)
     # xlim=c(-4,4),
     # ylim=c(0,0.5)
)
abs=seq(-10,10,0.001)
lines(abs,dnorm(abs),col='red',lty=3,lwd=6)
abline(v=mean(ech),col='red',lwd=2)

ks.test(ech, "pnorm")
# shapiro.test(ech)
shapiro.test(sample(ech,size=5000, replace=FALSE))
wilcox.test(ech, rnorm(length(ech)),paired=FALSE,"two.sided")


####----   Sampling from all simulations
cherriesWiloxPvalKing1000 = c()
cherriesKSPvalKing1000 = c()
cherriesShapiroPvalKing1000 = c()
cherriesWiloxPvalKing500 = c()
cherriesKSPvalKing500 = c()
cherriesShapiroPvalKing500 = c()
cherriesWiloxPvalKing100 = c()
cherriesKSPvalKing100 = c()
cherriesShapiroPvalKing100 = c()

for(repet in 1:100){
  ech1000=sample((statKing$nbCherries - sampleSize /3)/sqrt(2*sampleSize/45),size=1000,replace=FALSE)
  ech500=sample((statKing$nbCherries - sampleSize /3)/sqrt(2*sampleSize/45),size=500,replace=FALSE)
  ech100=sample((statKing$nbCherries - sampleSize /3)/sqrt(2*sampleSize/45),size=100,replace=FALSE)
  
  #---- tests of normality
  cherriesWiloxPvalKing1000 = c(cherriesWiloxPvalKing1000, wilcox.test(ech1000, rnorm(length(ech1000)),paired=FALSE,"two.sided")$p.value )
  cherriesKSPvalKing1000 = c(cherriesKSPvalKing1000, ks.test(ech1000, "pnorm")$p.value)
  cherriesShapiroPvalKing1000 = c(cherriesShapiroPvalKing1000, shapiro.test(ech1000)$p.value)
  
  #---- tests of normality
  cherriesWiloxPvalKing500 = c(cherriesWiloxPvalKing500, wilcox.test(ech500, rnorm(length(ech500)),paired=FALSE,"two.sided")$p.value )
  cherriesKSPvalKing500 = c(cherriesKSPvalKing500, ks.test(ech500, "pnorm")$p.value)
  cherriesShapiroPvalKing500 = c(cherriesShapiroPvalKing500, shapiro.test(ech500)$p.value)
  
  #---- tests of normality
  cherriesWiloxPvalKing100 = c(cherriesWiloxPvalKing100, wilcox.test(ech100, rnorm(length(ech100)),paired=FALSE,"two.sided")$p.value )
  cherriesKSPvalKing100 = c(cherriesKSPvalKing100, ks.test(ech100, "pnorm")$p.value)
  cherriesShapiroPvalKing100 = c(cherriesShapiroPvalKing100, shapiro.test(ech100)$p.value)
}

par(mfrow=c(1,3))
boxplot(cherriesWiloxPvalKing1000,cherriesKSPvalKing1000,cherriesShapiroPvalKing1000,
        names=c("M-W","K-S","Sh"),
        main="Cn sample = 1000",
        xlab="",
        ylab="p-value of the tests",
        ylim=c(0,1))
abline(h=0.05, col='red',lty=3,lwd=3)

boxplot(cherriesWiloxPvalKing500,cherriesKSPvalKing500,cherriesShapiroPvalKing500,
        names=c("M-W","K-S","Sh"),
        main="Cn sample = 500",
        xlab="",
        ylab="p-value of the tests",
        ylim=c(0,1))
abline(h=0.05, col='red',lty=3,lwd=3)

boxplot(cherriesWiloxPvalKing100,cherriesKSPvalKing100,cherriesShapiroPvalKing100,
        names=c("M-W","K-S","Sh"),
        main="Cn sample = 100",
        xlab="",
        ylab="p-value of the tests",
        ylim=c(0,1))
abline(h=0.05, col='red',lty=3,lwd=3)





##########################################################################################################
#################                    Tests of normality                     ##############################
#################                  External Branch Length Ln                ##############################
##########################################################################################################

###---- All simulations
meanExtBranch = weighted.mean(statKing$extBranch)
varExtBranch = weighted.mean((statKing$extBranch)*(statKing$extBranch)) - (meanExtBranch*meanExtBranch)
distrib=(statKing$extBranch-meanExtBranch)/(sqrt(0.5*varExtBranch))


par(mfrow = c(1,1))
val=0.3
hist(distrib, prob=TRUE,
     breaks = seq(min(distrib)-val, max(distrib)+val,val),
     xlab="normalized branch length",
     xlim=c(-5,10),
     main="")
abs=seq(-10,10,0.001)
lines(abs,dnorm(abs),col='red',lty=3,lwd=2)
abline(v=mean(distrib),col='red',lwd=2)


ks.test(distrib, "pnorm")
shapiro.test(sample(distrib,size=5000,replace = FALSE))
wilcox.test(distrib, rnorm(length(distrib)),paired=FALSE,"two.sided")





###---- Sampling from simulations

branchWiloxPvalKing1000 = c()
branchKSPvalKing1000 = c()
branchShapiroPvalKing1000 = c()
branchWiloxPvalKing500 = c()
branchKSPvalKing500 = c()
branchShapiroPvalKing500 = c()
branchWiloxPvalKing100 = c()
branchKSPvalKing100 = c()
branchShapiroPvalKing100 = c()

for(repet in 1:100){
  meanExtBranch = weighted.mean(statKing$extBranch)
  varExtBranch = weighted.mean((statKing$extBranch)*(statKing$extBranch)) - (meanExtBranch*meanExtBranch)
  
  distrib1000=sample((statKing$extBranch-meanExtBranch)/(sqrt(0.5*varExtBranch)),size=1000,replace=FALSE)
  distrib500=sample((statKing$extBranch-meanExtBranch)/(sqrt(0.5*varExtBranch)),size=500,replace=FALSE)
  distrib100=sample((statKing$extBranch-meanExtBranch)/(sqrt(0.5*varExtBranch)),size=100,replace=FALSE)
  
  #---- tests of normality
  branchWiloxPvalKing1000 = c(branchWiloxPvalKing1000, wilcox.test(distrib1000, rnorm(length(distrib1000)),paired=FALSE,"two.sided")$p.value )
  branchKSPvalKing1000 = c(branchKSPvalKing1000, ks.test(distrib1000, "pnorm")$p.value)
  branchShapiroPvalKing1000 = c(branchShapiroPvalKing1000, shapiro.test(distrib1000)$p.value)
  
  #---- tests of normality
  branchWiloxPvalKing500 = c(branchWiloxPvalKing500, wilcox.test(distrib500, rnorm(length(distrib500)),paired=FALSE,"two.sided")$p.value )
  branchKSPvalKing500 = c(branchKSPvalKing500, ks.test(distrib500, "pnorm")$p.value)
  branchShapiroPvalKing500 = c(branchShapiroPvalKing500, shapiro.test(distrib500)$p.value)
  
  #---- tests of normality
  branchWiloxPvalKing100 = c(branchWiloxPvalKing100, wilcox.test(distrib100, rnorm(length(distrib100)),paired=FALSE,"two.sided")$p.value )
  branchKSPvalKing100 = c(branchKSPvalKing100, ks.test(distrib100, "pnorm")$p.value)
  branchShapiroPvalKing100 = c(branchShapiroPvalKing100, shapiro.test(distrib100)$p.value)
}


par(mfrow=c(1,3))

boxplot(branchWiloxPvalKing1000,branchKSPvalKing1000,branchShapiroPvalKing1000,
        names=c("M-W","K-S","Sh"),
        main="Ln sample = 1000",
        xlab="",
        ylab="p-value of the tests",
        ylim=c(0,1))
abline(h=0.05, col='red',lty=3,lwd=3)

boxplot(branchWiloxPvalKing500,branchKSPvalKing500,branchShapiroPvalKing500,
        names=c("M-W","K-S","Sh"),
        main="Ln sample = 500",
        xlab="",
        ylab="p-value of the tests",
        ylim=c(0,1))
abline(h=0.05, col='red',lty=3,lwd=3)

boxplot(branchWiloxPvalKing100,branchKSPvalKing100,branchShapiroPvalKing100,
        names=c("Mann-Whitney test","K-S test","Shapiro test"),
        main="Ln sample = 100",
        xlab="",
        ylab="p-value of the tests",
        ylim=c(0,1))
abline(h=0.05, col='red',lty=3,lwd=3)




##########################################################################################################
#################                    Tests of normality                     ##############################
#################                           tMRCA                           ##############################
##########################################################################################################
sampleSize = 1000
nbSimuls = nrow(statKing)
tmrcaNeutral=rep(0,nbSimuls)
for(i in 2:sampleSize)
{
  tmrcaNeutral=tmrcaNeutral+rexp(nbSimuls,rate=(i*(i-1)))
}

nbSimuls = nrow(statKingCoeff)
tmrcaNeutralCoeff=rep(0,nbSimuls)
for(i in 2:sampleSize)
{
  tmrcaNeutralCoeff=tmrcaNeutralCoeff+rexp(nbSimuls,rate=(0.00765*i*(i-1)))
}

par(mfrow=c(1,1))
hist(statKing$tMRCA,
     main="",
     xlab="time before the MRCA",
     xlim=c(0,max(statKing$tMRCA)),
     prob=TRUE)
lines(density(tmrcaNeutral,na.rm = TRUE),col="red",lwd=2,lty=3)

qqplot(x = statKing$tMRCA, y= tmrcaNeutral)

par(mfrow=c(1,1))
hist(statKingCoeff$tMRCA,
     main="",
     xlab="time before the MRCA",
     xlim=c(0,max(statKingCoeff$tMRCA)),
     prob=TRUE)
lines(density(tmrcaNeutralCoeff,na.rm = TRUE),col="red",lwd=2,lty=3)



wilcox.test(tmrcaNeutral,statKing$tMRCA,paired=FALSE,"two.sided")
ks.test(tmrcaNeutral,statKing$tMRCA)

tmrcaTronc = statKing$tMRCA[statKing$tMRCA<335.2]
wilcox.test(tmrcaNeutral,tmrcaTronc,paired=FALSE,"two.sided")
ks.test(tmrcaNeutral,tmrcaTronc)

tmrcaTronc = statKing$tMRCA
for(i in 1:length(tmrcaTronc)){
  if(tmrcaTronc[i] > 335.2){
    tmrcaTronc[i] = 335.2
  }
}
wilcox.test(tmrcaNeutral,tmrcaTronc,paired=FALSE,"two.sided")
ks.test(tmrcaNeutral,tmrcaTronc)

summary(simuls)
###---- Sampling from simulations
tMRCAWiloxPvalKing1000 = c()
tMRCAKSPvalKing1000 = c()
tMRCAWiloxPvalKing500 = c()
tMRCAKSPvalKing500 = c()
tMRCAWiloxPvalKing100 = c()
tMRCAKSPvalKing100 = c()


sampleSize = 1000
for(repet in 1:100){
  
  tmrca1000=rep(0,1000)
  for(i in 2:sampleSize)
  {
    tmrca1000=tmrca1000+rexp(1000,rate=(i*(i-1)))
  }
  
  tmrca500=rep(0,500)
  for(i in 2:sampleSize)
  {
    tmrca500=tmrca500+rexp(500,rate=(i*(i-1)))
  }
  
  tmrca100=rep(0,100)
  for(i in 2:sampleSize)
  {
    tmrca100=tmrca100+rexp(100,rate=(i*(i-1)))
  }
  
  
  #---- tests of normality
  tMRCAWiloxPvalKing1000 = c(tMRCAWiloxPvalKing1000, wilcox.test(tmrca1000,sample(statKing$tMRCA,size=1000,replace = FALSE),paired=FALSE,"two.sided")$p.value)
  tMRCAKSPvalKing1000 = c(tMRCAKSPvalKing1000, ks.test(tmrca1000,sample(statKing$tMRCA,size=1000,replace = FALSE))$p.value)
  tMRCAWiloxPvalKing500 = c(tMRCAWiloxPvalKing500, wilcox.test(tmrca500,sample(statKing$tMRCA,size=500,replace = FALSE),paired=FALSE,"two.sided")$p.value)
  tMRCAKSPvalKing500 = c(tMRCAKSPvalKing500, ks.test(tmrca500,sample(statKing$tMRCA,size=500,replace = FALSE))$p.value)
  tMRCAWiloxPvalKing100 = c(tMRCAWiloxPvalKing100, wilcox.test(tmrca100,sample(statKing$tMRCA,size=100,replace = FALSE),paired=FALSE,"two.sided")$p.value)
  tMRCAKSPvalKing100 = c(tMRCAKSPvalKing100, ks.test(tmrca100,sample(statKing$tMRCA,size=100,replace = FALSE))$p.value)
}


par(mfrow=c(1,3))

boxplot(tMRCAWiloxPvalKing1000, tMRCAKSPvalKing1000,
        names=c("M-W","K-S"),
        main="tMRCA sample = 1000",
        xlab="",
        ylab="p-value of the tests",
        ylim=c(0,1))
abline(h=0.05, col='red',lty=3,lwd=3)

boxplot(tMRCAWiloxPvalKing500, tMRCAKSPvalKing500,
        names=c("M-W","K-S"),
        main="tMRCA sample = 500",
        xlab="",
        ylab="p-value of the tests",
        ylim=c(0,1))
abline(h=0.05, col='red',lty=3,lwd=3)

boxplot(tMRCAWiloxPvalKing100, tMRCAKSPvalKing100,
        names=c("M-W","K-S"),
        main="tMRCA sample = 100",
        xlab="",
        ylab="p-value of the tests",
        ylim=c(0,1))
abline(h=0.05, col='red',lty=3,lwd=3)






###---- Sampling from simulations
tMRCACoeffWiloxPvalKing1000 = c()
tMRCACoeffKSPvalKing1000 = c()
tMRCACoeffWiloxPvalKing500 = c()
tMRCACoeffKSPvalKing500 = c()
tMRCACoeffWiloxPvalKing100 = c()
tMRCACoeffKSPvalKing100 = c()


sampleSize = 1000
for(repet in 1:100){
  
  tMRCACoeff1000=rep(0,1000)
  for(i in 2:sampleSize)
  {
    tMRCACoeff1000=tMRCACoeff1000+rexp(1000,rate=(0.00765*i*(i-1)))
  }
  
  tMRCACoeff500=rep(0,500)
  for(i in 2:sampleSize)
  {
    tMRCACoeff500=tMRCACoeff500+rexp(500,rate=(0.00765*i*(i-1)))
  }
  
  tMRCACoeff100=rep(0,100)
  for(i in 2:sampleSize)
  {
    tMRCACoeff100=tMRCACoeff100+rexp(100,rate=(0.00765*i*(i-1)))
  }
  
  
  #---- tests of normality
  tMRCACoeffWiloxPvalKing1000 = c(tMRCACoeffWiloxPvalKing1000, wilcox.test(tMRCACoeff1000,sample(statKingCoeff$tMRCA,size=1000,replace = FALSE),paired=FALSE,"two.sided")$p.value)
  tMRCACoeffKSPvalKing1000 = c(tMRCACoeffKSPvalKing1000, ks.test(tMRCACoeff1000,sample(statKingCoeff$tMRCA,size=1000,replace = FALSE))$p.value)
  tMRCACoeffWiloxPvalKing500 = c(tMRCACoeffWiloxPvalKing500, wilcox.test(tMRCACoeff500,sample(statKingCoeff$tMRCA,size=500,replace = FALSE),paired=FALSE,"two.sided")$p.value)
  tMRCACoeffKSPvalKing500 = c(tMRCACoeffKSPvalKing500, ks.test(tMRCACoeff500,sample(statKingCoeff$tMRCA,size=500,replace = FALSE))$p.value)
  tMRCACoeffWiloxPvalKing100 = c(tMRCACoeffWiloxPvalKing100, wilcox.test(tMRCACoeff100,sample(statKingCoeff$tMRCA,size=100,replace = FALSE),paired=FALSE,"two.sided")$p.value)
  tMRCACoeffKSPvalKing100 = c(tMRCACoeffKSPvalKing100, ks.test(tMRCACoeff100,sample(statKingCoeff$tMRCA,size=100,replace = FALSE))$p.value)
}


par(mfrow=c(1,3))

boxplot(tMRCACoeffWiloxPvalKing1000, tMRCACoeffKSPvalKing1000,
        names=c("M-W","K-S"),
        main="tMRCACoeff sample = 1000",
        xlab="",
        ylab="p-value of the tests",
        ylim=c(0,1))
abline(h=0.05, col='red',lty=3,lwd=3)

boxplot(tMRCACoeffWiloxPvalKing500, tMRCACoeffKSPvalKing500,
        names=c("M-W","K-S"),
        main="tMRCACoeff sample = 500",
        xlab="",
        ylab="p-value of the tests",
        ylim=c(0,1))
abline(h=0.05, col='red',lty=3,lwd=3)

boxplot(tMRCACoeffWiloxPvalKing100, tMRCACoeffKSPvalKing100,
        names=c("M-W","K-S"),
        main="tMRCACoeff sample = 100",
        xlab="",
        ylab="p-value of the tests",
        ylim=c(0,1))
abline(h=0.05, col='red',lty=3,lwd=3)
