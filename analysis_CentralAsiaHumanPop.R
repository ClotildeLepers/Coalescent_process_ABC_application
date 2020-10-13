#########################################################
# ABC exploitation of the data of Chaix and Austerlitz
#########################################################

rm(list=ls())
library(abc)

donnees <- read.csv("CentralAsia_posteriorWeights.csv",sep=",")

dim(donnees)
head(donnees)

# indicator of the simulations where b0<b1

table(donnees$diffb)
bool=(donnees$diffb==1)
sum(bool)

  # in each column "weight.xx", we have the weights computed from the ABC for one line corresponding to the simulation used as data for computing these weights.
  # For this simulation, the corresponding line contains an NA. 100 of the simulations used as data are chosen in the set of simulations where b0=b1, and 100 in the set b0<b1.

sum(is.na(donnees))


# computation of the posterior probabilities of the alternatives H0 : b0=b1 and Ha : b0<b1
posteriors=rep(0,201)

posteriors=colSums(donnees[bool,14:214],na.rm=TRUE)
abs=seq(1,200,1)
plot(abs[1:100],posteriors[1:100],ylim=c(0,1),xlim=c(0,200),ylab="Pr(Ha | Sobs)",xlab="simulated sample")
points(abs[101:200],posteriors[101:200],col='red')
lines(c(100,100),c(0,1),lty=3)

legend(x="topright", legend=c("b0=b1", "b0<b1"), col=c("black", "red"),lwd=1,cex=0.8)

mean(posteriors[1:100])
median(posteriors[1:100])
lines(c(0,100),c(median(posteriors[1:100]),median(posteriors[1:100])))
text(x=45,y=0.05,"median prob = 0.4335",cex=0.75)

mean(posteriors[101:200])
median(posteriors[101:200])
lines(c(100,200),c(median(posteriors[101:200]),median(posteriors[101:200])),col="red")
text(x=155,y=0.05,"median prob = 0.5843",cex=0.75)

# choice of the threshold alpha defining the critical region {P(H_a | Sobs)>alpha}.

alpha=seq(0,1,0.01)
length(alpha)
type1=rep(0,101)
type2=rep(0,101)

for(i in 1:101)
{
  type1[i]=sum(posteriors[1:100]>alpha[i])/100
  type2[i]=sum(posteriors[101:200]<=alpha[i])/100
}

# Graphs of the two errors means (the a priori distribution gives probabilities 1/2,1/2 to the hypotheses "b0=b1" and "b0<b1")
error=(type1+type2)/2
alpha[which.min(error)]
which.min(error)
type1[which.min(error)]
type2[which.min(error)]
(median(posteriors[1:100])+median(posteriors[101:200]))/2
median(posteriors[1:200])

plot(alpha,error,ylim=c(0,1),xlim=c(0,1),ylab="P(Ha | Sobs, H0)+P(H0 | Sobs,Ha)",type='l')
lines(c(alpha[which.min(error)],alpha[which.min(error)]),c(0,1),lty=2)
lines(c((median(posteriors[1:100])+median(posteriors[101:200]))/2,(median(posteriors[1:100])+median(posteriors[101:200]))/2),c(0,1),lty=3)

legend(x="topright",legend=c("Sum of errors of types 1 and 2 minimal", "Mean between the medians under H0 and Ha"),col=c("black","black"),lty=c(2,3),cex=0.6)

#

plot(alpha,type1, type='l', ylim=c(0,1),ylab="Missclassification rate")
lines(alpha,type2,lty=3,col="red")

lines(c(alpha[which.min(error)],alpha[which.min(error)]),c(0,1),lty=2)
lines(c(0.525,0.525),c(0,1),lty=3)
legend(x=0.45,y=0.95,legend=c("Prob. of type 1 error", "Prob. of type 2 error"),col=c("black","red"),lty=c(1,3),cex=0.65)
text(x=0.7,y=0.5,"intersection at alpha=0.53",cex=0.7)
text(x=0.7,y=0.4,"min of sum at alpha=0.43",cex=0.7)

type1[40:60]
type2[40:60]
type1[54]
type2[54]
alpha[54]


####### Test with the real data from Heyer et al. (Central Asia)

  # the column 201 corresponds to the weights computed for the summary statistics from Chaix et al. 2007 and Heyer et al. 2015 (see file data_central_Asia.csv)

posteriors[201]
  # 0.4517<0.5, so we do not reject the null hypothesis H0:b0=b1. 
  # Hence,the test concludes that there is no significant difference of fertility rates between the two social orga-nizations: patrilineal and cognatic.

sum(posteriors[1:100]>posteriors[201])/100

summary(donnees$b0)
summary(donnees$b1)
plot(density(donnees$b0,weights=donnees$weightsHeyer), type='l',main="",xlab="",ylab="")
lines(density(donnees$b1,weights=donnees$weightsHeyer),col='red',lty=2)
legend(x="topright",legend=c("b0","b1"),col=c("black", "red"),lty=c(1,2),cex=2)


plot(density(donnees$b1-donnees$b0,weights=donnees$weightsHeyer), type='l',xlim=c(-0.1,2), main="",xlab="",ylab="")
m01=sum((donnees$b1-donnees$b0)*donnees$weightsHeyer)
lines(c(m01,m01),c(0,6.5),col="red", lty=3)

sd01=sqrt(sum((donnees$b1-donnees$b0)^2*donnees$weightsHeyer)-m01^2)

o01=order(donnees$b1-donnees$b0)
fdr01=cumsum(donnees$weightsHeyer[o01])
min(which(fdr01>0.95))
fdr01[18794:18795]
srank=donnees$b1[o01]-donnees$b0[o01]
srank[18794:18795]

length(donnees$b0)
sum(donnees$b0==donnees$b1)
sum(donnees$b0<donnees$b1)
bool=(donnees$diffb==1)
sum(bool)




####### Additional graphs for the case b0<b1:

head(donnees)
dim(donnees)
summary(donnees$b0)
summary(donnees$b1)

summary(donnees$b1[bool]-donnees$b0[bool])
hist(donnees$b1[bool]-donnees$b0[bool],main="",xlab="b1-b0",freq=FALSE)
lines(density(ecartsb),col="red")
temp=density(ecartsb)
temp$x
temp$y

plot(temp$x,temp$y)

ecartsb=rep(0,100)

for(i in 1:100)
{
  booli=(is.na(donnees[,112+i]))
  ecartsb[i]=donnees$b1[booli]-donnees$b0[booli]
}
o=order(ecartsb)
posteriorsdiff=posteriors[101:200]
pourcentagetype2=cumsum(posteriorsdiff[o]<=0.29)/100
plot(ecartsb[o],pourcentagetype2, type="l",ylim=c(0,1),xlab="b1-b0",ylab="Cum. freq. of Type II err.")
lines(density(ecartsb),col="red")
