
library(abc)

simuls10neutrals <-as.vector(read.table("simulations_10neutral.csv", header=F,sep = ",",fill=TRUE))
simuls20neutrals <-as.vector(read.table("simulations_20neutral.csv", header=F,sep = ",",fill=TRUE))
simuls30neutrals <-as.vector(read.table("simulations_30neutral.csv", header=F,sep = ",",fill=TRUE))
simuls40neutrals <-as.vector(read.table("simulations_40neutral.csv", header=F,sep = ",",fill=TRUE))
simuls50neutrals <-as.vector(read.table("simulations_50neutral.csv", header=F,sep = ",",fill=TRUE))
simuls60neutrals <-as.vector(read.table("simulations_60neutral.csv", header=F,sep = ",",fill=TRUE))
simuls70neutrals <-as.vector(read.table("simulations_70neutral.csv", header=F,sep = ",",fill=TRUE))
simuls80neutrals <-as.vector(read.table("simulations_80neutral.csv", header=F,sep = ",",fill=TRUE))
simuls90neutrals <-as.vector(read.table("simulations_90neutral.csv", header=F,sep = ",",fill=TRUE))
simuls100neutrals <-as.vector(read.table("simulations_100neutral.csv", header=F,sep = ",",fill=TRUE))

nomsColonnes=c("p","sigB","sigC","sigM","etaC","tsim","theta",
               "nbSubPop","meanN","varN","meanXpop","varXpop","meanXind","varXind",
               "nbSubPopSampled","varNrelSample","meanXpopSample","varXpopSample","meanXindSample","varXindSample",
               "meanNbAllleles","meanUnbiasedGeneDiversity","meanGeneDiversity","meanVarAlleles","meanMIndex",
               "meanNeiDa","varNeiDa","meanNeiDs","varNeiDs","meanDeltaMu","varDeltaMu",
               "meanFst","varFst","meanWeightedDa","varWeightedDa","meanWeightedDs","varWeightedDs",
               "meanWeightedDeltaMu","varWeightedDeltaMu","meanWeightedFst","varWeightedFst",
               "nbEvents","sackinIS","sumCoalExternalBranchLength","nbCherries","tMRCA")
colnames(simuls10neutrals)<-nomsColonnes
colnames(simuls20neutrals)<-nomsColonnes
colnames(simuls30neutrals)<-nomsColonnes
colnames(simuls40neutrals)<-nomsColonnes
colnames(simuls50neutrals)<-nomsColonnes
colnames(simuls60neutrals)<-nomsColonnes
colnames(simuls70neutrals)<-nomsColonnes
colnames(simuls80neutrals)<-nomsColonnes
colnames(simuls90neutrals)<-nomsColonnes
colnames(simuls100neutrals)<-nomsColonnes

simuls10neutrals<-simuls10neutrals[1:100000,]
simuls20neutrals<-simuls20neutrals[1:100000,]
simuls30neutrals<-simuls30neutrals[1:100000,]
simuls40neutrals<-simuls40neutrals[1:100000,]
simuls50neutrals<-simuls50neutrals[1:100000,]
simuls60neutrals<-simuls60neutrals[1:100000,]
simuls70neutrals<-simuls70neutrals[1:100000,]
simuls80neutrals<-simuls80neutrals[1:100000,]
simuls90neutrals<-simuls90neutrals[1:100000,]
simuls100neutrals<-simuls100neutrals[1:100000,]

simulsData10 <-as.vector(read.table("simulations_10neutral_A-D.csv", header=F,sep = ",",fill=TRUE))
simulsData20 <-as.vector(read.table("simulations_20neutral_A-D.csv", header=F,sep = ",",fill=TRUE))
simulsData30 <-as.vector(read.table("simulations_30neutral_A-D.csv", header=F,sep = ",",fill=TRUE))
simulsData40 <-as.vector(read.table("simulations_40neutral_A-D.csv", header=F,sep = ",",fill=TRUE))
simulsData50 <-as.vector(read.table("simulations_50neutral_A-D.csv", header=F,sep = ",",fill=TRUE))
simulsData60 <-as.vector(read.table("simulations_60neutral_A-D.csv", header=F,sep = ",",fill=TRUE))
simulsData70 <-as.vector(read.table("simulations_70neutral_A-D.csv", header=F,sep = ",",fill=TRUE))
simulsData80 <-as.vector(read.table("simulations_80neutral_A-D.csv", header=F,sep = ",",fill=TRUE))
simulsData90 <-as.vector(read.table("simulations_90neutral_A-D.csv", header=F,sep = ",",fill=TRUE))
simulsData100 <-as.vector(read.table("simulations_100neutral_A-D.csv", header=F,sep = ",",fill=TRUE))

colnames(simulsData10)<-nomsColonnes
colnames(simulsData20)<-nomsColonnes
colnames(simulsData30)<-nomsColonnes
colnames(simulsData40)<-nomsColonnes
colnames(simulsData50)<-nomsColonnes
colnames(simulsData60)<-nomsColonnes
colnames(simulsData70)<-nomsColonnes
colnames(simulsData80)<-nomsColonnes
colnames(simulsData90)<-nomsColonnes
colnames(simulsData100)<-nomsColonnes

#Realisation of abc on pseudo-data A to D (i0bs in 1:4), depending on the number of microsatellites
for(iObs in 1:4){
  for(nbMicrosat in c(10,100)){
    switch(nbMicrosat/10,
           {simuls<-simuls10neutrals
            simulsData = simulsData10},
           {simuls<-simuls20neutrals
           simulsData = simulsData20},
           {simuls<-simuls30neutrals
           simulsData = simulsData30},
           {simuls<-simuls40neutrals
           simulsData = simulsData40},
           {simuls<-simuls50neutrals
           simulsData = simulsData50},
           {simuls<-simuls60neutrals
           simulsData = simulsData60},
           {simuls<-simuls70neutrals
           simulsData = simulsData70},
           {simuls<-simuls80neutrals
           simulsData = simulsData80},
           {simuls<-simuls90neutrals
           simulsData = simulsData90},
           {simuls<-simuls100neutrals
           simulsData = simulsData100}
    )
    obs <-  simulsData[iObs,]
    
    abcRes=abc(data.frame(c(obs[c(8:14,15:20,21:25,26:41,42:46)])),
               data.frame(c(simuls[1:7])),
               data.frame(c(simuls[c(8:14,15:20,21:25,26:41,42:46)])),
               tol=0.01,method="neuralnet")
    
    abcRes_sample=abc(data.frame(c(obs[c(15:20,21:25,26:41)])),
                      data.frame(c(simuls[1:7])),
                      data.frame(c(simuls[c(15:20,21:25,26:41)])),
                      tol=0.01,method="neuralnet")
    
    abcRes_Data=abc(data.frame(c(obs[c(15,21:25,26:41)])),
                    data.frame(c(simuls[1:7])),
                    data.frame(c(simuls[c(15,21:25,26:41)])),
                    tol=0.01,method="neuralnet")
    
    abcRes_Data3=abc(data.frame(c(obs[c(15,17,18,21:25,26:41)])),
                     data.frame(c(simuls[1:7])),
                     data.frame(c(simuls[c(15,17,18,21:25,26:41)])),
                     tol=0.01,method="neuralnet")
    
    nomExport=paste0("2020_03_17_prior_posterior_L0",iObs,"_nbMicrosat_",nbMicrosat,".pdf")
    pdf(nomExport,width=6.2,height=7.2)
    ###---- plot of prior/posterior distributions of parameters
    par(mfrow = c(4,2))
    
    ymax=max(
      density(simuls$p,na.rm=TRUE)$y,
      density(abcRes$adj.values[,1])$y,
      density(abcRes_sample$adj.values[,1])$y,
      density(abcRes_Data$adj.values[,1])$y,
      density(abcRes_Data3$adj.values[,1])$y
    )
    xmin=min(
      quantile(simuls$p,na.rm=TRUE,probs=0.05),
      quantile(abcRes$adj.values[,1],probs=0.05),
      quantile(abcRes_sample$adj.values[,1],probs=0.05),
      quantile(abcRes_Data$adj.values[,1],probs=0.05),
      quantile(abcRes_Data3$adj.values[,1],probs=0.05)
    )
    xmax=max(
      quantile(simuls$p,na.rm=TRUE,probs=0.95),
      quantile(abcRes$adj.values[,1],probs=0.95),
      quantile(abcRes_sample$adj.values[,1],probs=0.95),
      quantile(abcRes_Data$adj.values[,1],probs=0.95),
      quantile(abcRes_Data3$adj.values[,1],probs=0.95)
    )
    plot(density(simuls$p,na.rm=TRUE),col="black", lty=2, xlab="", main="p", ylim=c(0,ymax+ymax/10),xlim=c(0,xmax))
    abline(v=obs$p, col="red")
    lines(density(abcRes$adj.values[,1]),col="blue")
    lines(density(abcRes_sample$adj.values[,1]),col="magenta")
    lines(density(abcRes_Data$adj.values[,1]),col="orange")
    lines(density(abcRes_Data3$adj.values[,1]),col="red")
    
    
    ymax=max(
      density(simuls$sigM,na.rm=TRUE)$y,
      density(abcRes$adj.values[,4])$y,
      density(abcRes_sample$adj.values[,4])$y,
      density(abcRes_Data$adj.values[,4])$y,
      density(abcRes_Data3$adj.values[,4])$y
    )
    xmin=min(
      quantile(simuls$sigM,na.rm=TRUE,probs=0.05),
      quantile(abcRes$adj.values[,4],probs=0.05),
      quantile(abcRes_sample$adj.values[,4],probs=0.05),
      quantile(abcRes_Data$adj.values[,4],probs=0.05),
      quantile(abcRes_Data3$adj.values[,4],probs=0.05)
    )
    xmax=max(
      quantile(simuls$sigM,na.rm=TRUE,probs=0.95),
      quantile(abcRes$adj.values[,4],probs=0.95),
      quantile(abcRes_sample$adj.values[,4],probs=0.95),
      quantile(abcRes_Data$adj.values[,4],probs=0.95),
      quantile(abcRes_Data3$adj.values[,4],probs=0.95)
    )
    plot(density(simuls$sigM,na.rm=TRUE),col="black", lty=2, xlab="", main="sigM", ylim=c(0,ymax+ymax/10),xlim=c(0,xmax))
    abline(v=obs$sigM, col="red")
    lines(density(abcRes$adj.values[,4]),col="blue")
    lines(density(abcRes_sample$adj.values[,4]),col="magenta")
    lines(density(abcRes_Data$adj.values[,4]),col="orange")
    lines(density(abcRes_Data3$adj.values[,4]),col="red")
    
    
    ymax=max(
      density((simuls$p * simuls$sigM)^2 ,na.rm=TRUE)$y,
      density((abcRes$adj.values[,1] * abcRes$adj.values[,4])^2)$y,
      density((abcRes_sample$adj.values[,1]*abcRes_sample$adj.values[,4])^2)$y,
      density((abcRes_Data$adj.values[,1]*abcRes_Data$adj.values[,4])^2)$y,
      density((abcRes_Data3$adj.values[,1]*abcRes_Data3$adj.values[,4])^2)$y
    )
    xmin=min(
      quantile((abcRes$adj.values[,1] * abcRes$adj.values[,4])^2,probs=0.05),
      quantile((abcRes_sample$adj.values[,1]*abcRes_sample$adj.values[,4])^2,probs=0.05),
      quantile((abcRes_Data$adj.values[,1]*abcRes_Data$adj.values[,4])^2,probs=0.05),
      quantile((abcRes_Data3$adj.values[,1]*abcRes_Data3$adj.values[,4])^2,probs=0.05)
    )
    xmax=max(
      quantile((simuls$p * simuls$sigM)^2,na.rm=TRUE,probs=0.95),
      quantile((abcRes$adj.values[,1] * abcRes$adj.values[,4])^2,probs=0.95),
      quantile((abcRes_sample$adj.values[,1]*abcRes_sample$adj.values[,4])^2,probs=0.95),
      quantile((abcRes_Data$adj.values[,1]*abcRes_Data$adj.values[,4])^2,probs=0.95),
      quantile((abcRes_Data3$adj.values[,1]*abcRes_Data3$adj.values[,4])^2,probs=0.95)
    )
    plot(density((simuls$p * simuls$sigM)^2,na.rm=TRUE),col="black", lty=2, xlab="", main="p*sigM", ylim=c(0,ymax+ymax/10),xlim=c(0,xmax))
    abline(v=(obs$p*obs$sigM)^2, col="red")
    lines(density((abcRes$adj.values[,1] * abcRes$adj.values[,4])^2),col="blue")
    lines(density((abcRes_sample$adj.values[,1]*abcRes_sample$adj.values[,4])^2),col="magenta")
    lines(density((abcRes_Data$adj.values[,1]*abcRes_Data$adj.values[,4])^2),col="orange")
    lines(density((abcRes_Data3$adj.values[,1]*abcRes_Data3$adj.values[,4])^2),col="red")
    
    
    ymax=max(
      density(simuls$sigB,na.rm=TRUE)$y,
      density(abcRes$adj.values[,2])$y,
      density(abcRes_sample$adj.values[,2])$y,
      density(abcRes_Data$adj.values[,2])$y,
      density(abcRes_Data3$adj.values[,2])$y
    )
    xmin=min(
      quantile(simuls$sigB,na.rm=TRUE,probs=0.05),
      quantile(abcRes$adj.values[,2],probs=0.05),
      quantile(abcRes_sample$adj.values[,2],probs=0.05),
      quantile(abcRes_Data$adj.values[,2],probs=0.05),
      quantile(abcRes_Data3$adj.values[,2],probs=0.05)
    )
    xmax=max(
      quantile(simuls$sigB,na.rm=TRUE,probs=0.95),
      quantile(abcRes$adj.values[,2],probs=0.95),
      quantile(abcRes_sample$adj.values[,2],probs=0.95),
      quantile(abcRes_Data$adj.values[,2],probs=0.95),
      quantile(abcRes_Data3$adj.values[,2],probs=0.95)
    )
    plot(density(simuls$sigB,na.rm=TRUE),col="black", lty=2, xlab="", main="sigB", ylim=c(0,ymax+ymax/10),xlim=c(0,xmax))
    abline(v=obs$sigB, col="red")
    lines(density(abcRes$adj.values[,2]),col="blue")
    lines(density(abcRes_sample$adj.values[,2]),col="magenta")
    lines(density(abcRes_Data$adj.values[,2]),col="orange")
    lines(density(abcRes_Data3$adj.values[,2]),col="red")
    
    
    ymax=max(
      density(simuls$tsim,na.rm=TRUE)$y,
      density(abcRes$adj.values[,6])$y,
      density(abcRes_sample$adj.values[,6])$y,
      density(abcRes_Data$adj.values[,6])$y,
      density(abcRes_Data3$adj.values[,6])$y
    )
    xmin=min(
      quantile(simuls$tsim,na.rm=TRUE,probs=0.05),
      quantile(abcRes$adj.values[,6],probs=0.05),
      quantile(abcRes_sample$adj.values[,6],probs=0.05),
      quantile(abcRes_Data$adj.values[,6],probs=0.05),
      quantile(abcRes_Data3$adj.values[,6],probs=0.05)
    )
    xmax=max(
      quantile(simuls$tsim,na.rm=TRUE,probs=0.95),
      quantile(abcRes$adj.values[,6],probs=0.95),
      quantile(abcRes_sample$adj.values[,6],probs=0.95),
      quantile(abcRes_Data$adj.values[,6],probs=0.95),
      quantile(abcRes_Data3$adj.values[,6],probs=0.95)
    )
    plot(density(simuls$tsim,na.rm=TRUE),col="black", lty=2, xlab="", main="tsim", ylim=c(0,ymax+ymax/10),xlim=c(0,xmax))
    abline(v=obs$tsim, col="red")
    lines(density(abcRes$adj.values[,6]),col="blue")
    lines(density(abcRes_sample$adj.values[,6]),col="magenta")
    lines(density(abcRes_Data$adj.values[,6]),col="orange")
    lines(density(abcRes_Data3$adj.values[,6]),col="red")
    
    
    ymax=max(
      density(simuls$theta,na.rm=TRUE)$y,
      density(abcRes$adj.values[,7])$y,
      density(abcRes_sample$adj.values[,7])$y,
      density(abcRes_Data$adj.values[,7])$y,
      density(abcRes_Data3$adj.values[,7])$y
    )
    xmin=min(
      quantile(simuls$theta,na.rm=TRUE,probs=0.05),
      quantile(abcRes$adj.values[,7],probs=0.05),
      quantile(abcRes_sample$adj.values[,7],probs=0.05),
      quantile(abcRes_Data$adj.values[,7],probs=0.05),
      quantile(abcRes_Data3$adj.values[,7],probs=0.05)
    )
    xmax=max(
      quantile(simuls$theta,na.rm=TRUE,probs=0.95),
      quantile(abcRes$adj.values[,7],probs=0.95),
      quantile(abcRes_sample$adj.values[,7],probs=0.95),
      quantile(abcRes_Data$adj.values[,7],probs=0.95),
      quantile(abcRes_Data3$adj.values[,7],probs=0.95)
    )
    plot(density(simuls$theta,na.rm=TRUE),col="black", lty=2, xlab="", main="q", ylim=c(0,ymax+ymax/10),xlim=c(0,xmax))
    abline(v=obs$theta, col="red")
    lines(density(abcRes$adj.values[,7]),col="blue")
    lines(density(abcRes_sample$adj.values[,7]),col="magenta")
    lines(density(abcRes_Data$adj.values[,7]),col="orange")
    lines(density(abcRes_Data3$adj.values[,7]),col="red")
    
    
    ymax=max(
      density(simuls$sigC,na.rm=TRUE)$y,
      density(abcRes$adj.values[,3])$y,
      density(abcRes_sample$adj.values[,3])$y,
      density(abcRes_Data$adj.values[,3])$y,
      density(abcRes_Data3$adj.values[,3])$y
    )
    xmin=min(
      quantile(simuls$sigC,na.rm=TRUE,probs=0.05),
      quantile(abcRes$adj.values[,3],probs=0.05),
      quantile(abcRes_sample$adj.values[,3],probs=0.05),
      quantile(abcRes_Data$adj.values[,3],probs=0.05),
      quantile(abcRes_Data3$adj.values[,3],probs=0.05)
    )
    xmax=max(
      quantile(simuls$sigC,na.rm=TRUE,probs=0.95),
      quantile(abcRes$adj.values[,3],probs=0.95),
      quantile(abcRes_sample$adj.values[,3],probs=0.95),
      quantile(abcRes_Data$adj.values[,3],probs=0.95),
      quantile(abcRes_Data3$adj.values[,3],probs=0.95)
    )
    plot(density(simuls$sigC,na.rm=TRUE),col="black", lty=2, xlab="", main="sigC", ylim=c(0,ymax+ymax/10),xlim=c(0,xmax))
    abline(v=obs$sigC, col="red")
    lines(density(abcRes$adj.values[,3]),col="blue")
    lines(density(abcRes_sample$adj.values[,3]),col="magenta")
    lines(density(abcRes_Data$adj.values[,3]),col="orange")
    lines(density(abcRes_Data3$adj.values[,3]),col="red")
    
    
    ymax=max(
      density(simuls$etaC,na.rm=TRUE)$y,
      density(abcRes$adj.values[,5])$y,
      density(abcRes_sample$adj.values[,5])$y,
      density(abcRes_Data$adj.values[,5])$y,
      density(abcRes_Data3$adj.values[,5])$y
    )
    xmin=min(
      quantile(simuls$etaC,na.rm=TRUE,probs=0.05),
      quantile(abcRes$adj.values[,5],probs=0.05),
      quantile(abcRes_sample$adj.values[,5],probs=0.05),
      quantile(abcRes_Data$adj.values[,5],probs=0.05),
      quantile(abcRes_Data3$adj.values[,5],probs=0.05)
    )
    xmax=max(
      quantile(simuls$etaC,na.rm=TRUE,probs=0.95),
      quantile(abcRes$adj.values[,5],probs=0.95),
      quantile(abcRes_sample$adj.values[,5],probs=0.95),
      quantile(abcRes_Data$adj.values[,5],probs=0.95),
      quantile(abcRes_Data3$adj.values[,5],probs=0.95)
    )
    plot(density(simuls$etaC,na.rm=TRUE),col="black", lty=2, xlab="", main="etaC", ylim=c(0,ymax+ymax/10),xlim=c(0,xmax))
    abline(v=obs$etaC, col="red")
    lines(density(abcRes$adj.values[,5]),col="blue")
    lines(density(abcRes_sample$adj.values[,5]),col="magenta")
    lines(density(abcRes_Data$adj.values[,5]),col="orange")
    lines(density(abcRes_Data3$adj.values[,5]),col="red")
    
    dev.off()
    
  }
}


# For a given pseudo-data, realize abc depending on the number of microsattelites and evaluate the square deviation.
sommesRes=c()
sommesRes_sample=c()
sommesRes_Data=c()
sommesRes_Data3=c()
iObs=4
for(nbMicrosat in seq(10,100,10)){
  switch(nbMicrosat/10,
         {simuls<-simuls10neutrals
         simulsData = simulsData10},
         {simuls<-simuls20neutrals
         simulsData = simulsData20},
         {simuls<-simuls30neutrals
         simulsData = simulsData30},
         {simuls<-simuls40neutrals
         simulsData = simulsData40},
         {simuls<-simuls50neutrals
         simulsData = simulsData50},
         {simuls<-simuls60neutrals
         simulsData = simulsData60},
         {simuls<-simuls70neutrals
         simulsData = simulsData70},
         {simuls<-simuls80neutrals
         simulsData = simulsData80},
         {simuls<-simuls90neutrals
         simulsData = simulsData90},
         {simuls<-simuls100neutrals
         simulsData = simulsData100}
  )
  obs <-  simulsData[iObs,]
  abcRes=abc(data.frame(c(obs[c(8:14,15:20,21:25,26:41,42:46)])),
             data.frame(c(simuls[1:7])),
             data.frame(c(simuls[c(8:14,15:20,21:25,26:41,42:46)])),
             tol=0.01,method="neuralnet")
  ecartsRes=cbind((abcRes$adj.values[,1]-obs[1,1])^2*abcRes$weights,
                  (abcRes$adj.values[,4]-obs[1,4])^2*abcRes$weights,
                  ((abcRes$adj.values[,4]*abcRes$adj.values[,1])-(obs[1,4]*obs[1,1]))^2*abcRes$weights,
                  (abcRes$adj.values[,2]-obs[1,2])^2*abcRes$weights,
                  (abcRes$adj.values[,6]-obs[1,6])^2*abcRes$weights,
                  (abcRes$adj.values[,7]-obs[1,7])^2*abcRes$weights,
                  (abcRes$adj.values[,3]-obs[1,3])^2*abcRes$weights,
                  (abcRes$adj.values[,5]-obs[1,5])^2*abcRes$weights)
  sommesRes=rbind(sommesRes,colSums(ecartsRes))
  abcRes_sample=abc(data.frame(c(obs[c(15:20,21:25,26:41)])),
                    data.frame(c(simuls[1:7])),
                    data.frame(c(simuls[c(15:20,21:25,26:41)])),
                    tol=0.01,method="neuralnet")
  ecartsRes_sample=cbind((abcRes_sample$adj.values[,1]-obs[1,1])^2*abcRes_sample$weights,
                         (abcRes_sample$adj.values[,4]-obs[1,4])^2*abcRes_sample$weights,
                         ((abcRes_sample$adj.values[,4]*abcRes_sample$adj.values[,1])-(obs[1,4]*obs[1,1]))^2*abcRes_sample$weights,
                         (abcRes_sample$adj.values[,2]-obs[1,2])^2*abcRes_sample$weights,
                         (abcRes_sample$adj.values[,6]-obs[1,6])^2*abcRes_sample$weights,
                         (abcRes_sample$adj.values[,7]-obs[1,7])^2*abcRes_sample$weights,
                         (abcRes_sample$adj.values[,3]-obs[1,3])^2*abcRes_sample$weights,
                         (abcRes_sample$adj.values[,5]-obs[1,5])^2*abcRes_sample$weights)
  sommesRes_sample=rbind(sommesRes_sample,colSums(ecartsRes_sample))
  abcRes_Data=abc(data.frame(c(obs[c(15,21:25,26:41)])),
                  data.frame(c(simuls[1:7])),
                  data.frame(c(simuls[c(15,21:25,26:41)])),
                  tol=0.01,method="neuralnet")
  ecartsRes_Data=cbind((abcRes_Data$adj.values[,1]-obs[1,1])^2*abcRes_Data$weights,
                       (abcRes_Data$adj.values[,4]-obs[1,4])^2*abcRes_Data$weights,
                       ((abcRes_Data$adj.values[,4]*abcRes_Data$adj.values[,1])-(obs[1,4]*obs[1,1]))^2*abcRes_Data$weights,
                       (abcRes_Data$adj.values[,2]-obs[1,2])^2*abcRes_Data$weights,
                       (abcRes_Data$adj.values[,6]-obs[1,6])^2*abcRes_Data$weights,
                       (abcRes_Data$adj.values[,7]-obs[1,7])^2*abcRes_Data$weights,
                       (abcRes_Data$adj.values[,3]-obs[1,3])^2*abcRes_Data$weights,
                       (abcRes_Data$adj.values[,5]-obs[1,5])^2*abcRes_Data$weights)
  sommesRes_Data=rbind(sommesRes_Data,colSums(ecartsRes_Data))
  abcRes_Data3=abc(data.frame(c(obs[c(15,17,18,21:25,26:41)])),
                   data.frame(c(simuls[1:7])),
                   data.frame(c(simuls[c(15,17,18,21:25,26:41)])),
                   tol=0.01,method="neuralnet")
  ecartsRes_Data3=cbind((abcRes_Data3$adj.values[,1]-obs[1,1])^2*abcRes_Data3$weights,
                        (abcRes_Data3$adj.values[,4]-obs[1,4])^2*abcRes_Data3$weights,
                        ((abcRes_Data3$adj.values[,4]*abcRes_Data3$adj.values[,1])-(obs[1,4]*obs[1,1]))^2*abcRes_Data3$weights,
                        (abcRes_Data3$adj.values[,2]-obs[1,2])^2*abcRes_Data3$weights,
                        (abcRes_Data3$adj.values[,6]-obs[1,6])^2*abcRes_Data3$weights,
                        (abcRes_Data3$adj.values[,7]-obs[1,7])^2*abcRes_Data3$weights,
                        (abcRes_Data3$adj.values[,3]-obs[1,3])^2*abcRes_Data3$weights,
                        (abcRes_Data3$adj.values[,5]-obs[1,5])^2*abcRes_Data3$weights)
  sommesRes_Data3=rbind(sommesRes_Data3,colSums(ecartsRes_Data3))
  
}

nomExport=paste0("2020_03_17_square_deviation_L0",iObs,".pdf")
pdf(nomExport,width=6.2,height=7.2)
par(mfrow = c(4,2))
microsatNb=c(10,20,30,40,50,60,70,80,90,100)
plot(sommesRes[,1]~microsatNb,xlab="number of microsatellites",ylab="sum of square deviations",main="p",col="blue",
     ylim=c(0,max(sommesRes[,1],sommesRes_sample[,1],sommesRes_Data[,1],sommesRes_Data3[,1])))
lines(sommesRes_sample[,1]~microsatNb,col="magenta",type="p")
lines(sommesRes_Data[,1]~microsatNb,col="orange",type="p")
lines(sommesRes_Data3[,1]~microsatNb,col="red",type="p")
plot(sommesRes[,2]~microsatNb,xlab="number of microsatellites",ylab="sum of square deviations",main="sigM",col="blue",
     ylim=c(0,max(sommesRes[,2],sommesRes_sample[,2],sommesRes_Data[,2],sommesRes_Data3[,2])))
lines(sommesRes_sample[,2]~microsatNb,col="magenta",type="p")
lines(sommesRes_Data[,2]~microsatNb,col="orange",type="p")
lines(sommesRes_Data3[,2]~microsatNb,col="red",type="p")
plot(sommesRes[,3]~microsatNb,xlab="number of microsatellites",ylab="sum of square deviations",main="p*sigM",col="blue",
     ylim=c(0,max(sommesRes[,3],sommesRes_sample[,3],sommesRes_Data[,3],sommesRes_Data3[,3])))
lines(sommesRes_sample[,3]~microsatNb,col="magenta",type="p")
lines(sommesRes_Data[,3]~microsatNb,col="orange",type="p")
lines(sommesRes_Data3[,3]~microsatNb,col="red",type="p")
plot(sommesRes[,4]~microsatNb,xlab="number of microsatellites",ylab="sum of square deviations",main="sigB",col="blue",
     ylim=c(0,max(sommesRes[,4],sommesRes_sample[,4],sommesRes_Data[,4],sommesRes_Data3[,4])))
lines(sommesRes_sample[,4]~microsatNb,col="magenta",type="p")
lines(sommesRes_Data[,4]~microsatNb,col="orange",type="p")
lines(sommesRes_Data3[,4]~microsatNb,col="red",type="p")
plot(sommesRes[,5]~microsatNb,xlab="number of microsatellites",ylab="sum of square deviations",main="tsim",col="blue",
     ylim=c(0,max(sommesRes[,5],sommesRes_sample[,5],sommesRes_Data[,5],sommesRes_Data3[,5])))
lines(sommesRes_sample[,5]~microsatNb,col="magenta",type="p")
lines(sommesRes_Data[,5]~microsatNb,col="orange",type="p")
lines(sommesRes_Data3[,5]~microsatNb,col="red",type="p")
plot(sommesRes[,6]~microsatNb,xlab="number of microsatellites",ylab="sum of square deviations",main="q",col="blue",
     ylim=c(0,max(sommesRes[,6],sommesRes_sample[,6],sommesRes_Data[,6],sommesRes_Data3[,6])))
lines(sommesRes_sample[,6]~microsatNb,col="magenta",type="p")
lines(sommesRes_Data[,6]~microsatNb,col="orange",type="p")
lines(sommesRes_Data3[,6]~microsatNb,col="red",type="p")
plot(sommesRes[,7]~microsatNb,xlab="number of microsatellites",ylab="sum of square deviations",main="sigC",col="blue",
     ylim=c(0,max(sommesRes[,7],sommesRes_sample[,7],sommesRes_Data[,7],sommesRes_Data3[,7])))
lines(sommesRes_sample[,7]~microsatNb,col="magenta",type="p")
lines(sommesRes_Data[,7]~microsatNb,col="orange",type="p")
lines(sommesRes_Data3[,7]~microsatNb,col="red",type="p")
plot(sommesRes[,8]~microsatNb,xlab="number of microsatellites",ylab="sum of square deviations",main="etaC",col="blue",
     ylim=c(0,max(sommesRes[,8],sommesRes_sample[,8],sommesRes_Data[,8],sommesRes_Data3[,8])))
lines(sommesRes_sample[,8]~microsatNb,col="magenta",type="p")
lines(sommesRes_Data[,8]~microsatNb,col="orange",type="p")
lines(sommesRes_Data3[,8]~microsatNb,col="red",type="p")


dev.off()


nomExport=paste0("2020_03_17_square_deviation_all_parameters_L0",iObs,".pdf")
pdf(nomExport,width=3.1,height=3.1)
par(mfrow = c(1,1))

plot(rowSums(sommesRes[,c(1,2,4:8)])~microsatNb,xlab="number of microsatellites",ylab="sum of square deviations",main="all parameters",col="blue",
     ylim=c(0,max(sommesRes[,c(1,2,4:8)],sommesRes_sample[,c(1,2,4:8)],sommesRes_Data[,c(1,2,4:8)],sommesRes_Data3[,c(1,2,4:8)])))
lines(rowSums(sommesRes_sample[,c(1,2,4:8)])~microsatNb,col="magenta",type="p")
lines(rowSums(sommesRes_Data[,c(1,2,4:8)])~microsatNb,col="orange",type="p")
lines(rowSums(sommesRes_Data3[,c(1,2,4:8)])~microsatNb,col="red",type="p")

dev.off()


# for simulations with 100microsattelites, evaluate the difference between our model and normality.
simuls<-simuls100neutrals
simulsData = simulsData100
for(iObs in c(1:6)){
  obs <-  simulsData[iObs,]
  abcRes=abc(data.frame(c(obs[c(8:14,15:20,21:25,26:41,42:46)])),
             data.frame(c(simuls[1:7])),
             data.frame(c(simuls[c(8:14,15:20,21:25,26:41,42:46)])),
             tol=0.01,method="neuralnet")
  #---- hist number of cherries vs normality
  abcRes$names$statistics.names[38]
  sampleSize = 1000
  abs=seq(-10,10,0.001)
  
  nomExport=paste0("2020_03_17_hist_Cn_L0",iObs,".pdf")
  pdf(nomExport,width=4.5,height=5.5)
  par(mfrow = c(1,1))
  hist((abcRes$ss[,38] - sampleSize /3)/sqrt(2*sampleSize/45), weights=abcRes$weights ,probability = TRUE,xlim=c(-4,4),ylim=c(0,0.5),main="number of cherries",xlab="")
  lines(abs,dnorm(abs),col='black',lty=2)
  dev.off()

  #---- hist external branch length vs normality
  abcRes$names$statistics.names[37]
  sampleSize = 1000
  abs=seq(-10,10,0.001)
  meanExtBranch = weighted.mean(abcRes$ss[,37] , w=abcRes$weights/sum(abcRes$weights))
  varExtBranch = weighted.mean((abcRes$ss[,37])*(abcRes$ss[,37]) , w=abcRes$weights/sum(abcRes$weights)) - (meanExtBranch*meanExtBranch)
  
  distrib=(abcRes$ss[,37]-meanExtBranch)/(sqrt(0.5*varExtBranch))
  
  nomExport=paste0("2020_03_17_hist_Ln_L0",iObs,".pdf")
  pdf(nomExport,width=4.5,height=5.5)
  par(mfrow = c(1,1))
  hist((abcRes$ss[,37] -meanExtBranch)/(sqrt(0.5*varExtBranch)) ,probability = TRUE,breaks = seq(-4,8,0.5),main="coalescent external branch length",xlim=c(-4,8),xlab="",ylim=c(0.0,0.5))
  lines(abs,dnorm(abs),col='black',lty=2)
  dev.off()
  
  #----- hist tmrca vs normality
  nomExport=paste0("2020_03_17_hist_tmrca_L0",iObs,".pdf")
  pdf(nomExport,width=4.5,height=5.5)
  par(mfrow = c(1,1))
  abcRes$names$statistics.names[39]
  hist(abcRes$ss[,39],weights=abcRes$weights,probability = TRUE, na.rm = TRUE,breaks = seq(0,5500,100), main="time before the MRCA",xlab="",xlim=c(0,2500),ylim=c(0,0.01),prob=TRUE)
  #---- estimations de la distribution de tmrca sous un modèle neutre
  abcRes$names$parameter.names[2]
  val_sigB = abcRes$unadj.values[,2]
  abcRes$names$statistics.names[4]
  val_meanXpop = abcRes$ss[,4]
  abcRes$names$parameter.names[5]
  val_etaC = abcRes$unadj.values[,5]
  abcRes$names$parameter.names[3]
  val_sigC = abcRes$unadj.values[,3]
  
  bx= exp(-(val_meanXpop/val_sigB)^2)/2
  nx=(bx/(val_etaC * exp(-(val_meanXpop/val_sigC)^2)/2) )
  
  nbSimuls = nrow(abcRes$unadj.values)
  tmrcaNeutral=rep(0,nbSimuls)
  sampleSize = 1000
  for(i in 2:sampleSize)
  {
    tmrcaNeutral=tmrcaNeutral+rexp(nbSimuls,rate=(weighted.mean(bx/(2*nx), w=abcRes$weights/sum(abcRes$weights))*i*(i-1)))
  }
  lines(density(tmrcaNeutral,na.rm = TRUE),col="black",lty=2)
  dev.off()
  
}