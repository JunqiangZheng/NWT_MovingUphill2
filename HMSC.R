
library(HMSC)
library(vegan)
library(corrplot)
library(circlize)


####### Real data #####

comm.dataALL

write.csv(comm.dataALL,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/comm.dataALL.csv",row.names=F)

#get plant cover data to use as "light" and merge with biogeo
plantcov<-read.csv("/Users/farrer/Dropbox/Niwot Moving Uphill/Analysis/Niwot_MovingUpHill_plots2015.csv")
names(plantcov)[1]<-"Sample_name"
plantcov$plantcover<-(plantcov$MOSS+plantcov$VEG)/plantcov$TOTAL
plantcov$plantcover[which(plantcov$Sample_name==64)]<-0 #this I'm 99% sure has zero plants, not sure why bare and rock were recorded as 0s. 1% chance that we just forgot to sample it, but it doesn't have nematodes either so I think it was bare
head(plantcov)

head(biogeo2)
biogeo3<-merge(biogeo2,plantcov,"Sample_name") #with merge, they got put in the right order
cbind(biogeo3$Sample_name,row.names(hmscY))
#biogeo4<-biogeo3[order(biogeo3$Sample_name),]#now they are in the same order

write.csv(biogeo3,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/biogeo3.csv",row.names=F)


##### Start from here: read in files #####
biogeo3<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/biogeo3.csv")

comm.dataALL<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/comm.dataALL.csv")

hmscY<-comm.dataALL[,32:6704]
hmscY<-comm.dataALL[,32:1000]
hmscY<-comm.dataALL[,c(32:1000,5250:6000)]
hmscY<-comm.dataALL[,32:200] #for practice

rownames(hmscY)<-comm.dataALL$X.SampleID

hmscX<-data.frame(inter=rep(1,90),pH=biogeo3$pH.x,moisture=biogeo3$moisture.x,snowdepth=biogeo3$snowdepth,TN=biogeo3$TN,TC=biogeo3$TC,plantcov=biogeo3$plantcov,lomehi=biogeo3$lomehi) #there are NAs in the pH data and possibly more
rownames(hmscX)<-biogeo3$X.SampleID

#take out any rows with NAs
ind<-which(!is.na(rowSums(hmscX[,1:dim(hmscX)[2]-1])))
#ind<-which(!is.na(hmscX[,1]))
hmscXb<-hmscX[ind,]
hmscYb<-hmscY[ind,]

#readjust the lo me hi categories based on NAs and # samples in each category
XXXX

#select lo/hi
ind<-which(hmscXb$lomehi=="hi")
hmscXc<-hmscXb[ind,]
hmscYc<-hmscYb[ind,]

#select species with greater than 10 (11 or more) occurrences and remove lo me hi (since you can't have text in a matrix or the whole matrix becomes character)
ind<-which(colSums(hmscYc>0)>10)
hmscYd<-hmscYc[,ind]
hmscXd<-hmscXc[,1:dim(hmscXc)[2]-1]#
dim(hmscYd)

#testing a really strong correlation with pH
#hmscYd[,2]<-rnorm(29,mean=2*formdata$X[,2]+1,sd=.02)

#select only species 2 and 3
#hmscYd<-cbind(hmscYd[,3],hmscYd[,3])

#the y data are not normal (the only options I have are normal, binary, poisson, overdispersed poisson), so I could do a sqrt transformation on Y (log(0) is -Inf). log(x+1) doesn't work since the proportions are so low, could do log(x*100+1) but the sqrt actually makes it more normal
hmscYe<-sqrt(hmscYd*100)
#hmscYe<-log(hmscYd*100+1)
#hist(hmscYd[,6])
#hist(hmscYe[,12])

#check if the values are too low that some tolerance is messing up the CI estimates, yes important to scale y
hmscYf<-scale(d)

#make them matrices
hmscXe<-as.matrix(hmscXd)
hmscYg<-as.matrix(hmscYe)#hmscYf



#make the random (residual matrix)
pimat<-data.frame(plot=1:dim(hmscYe)[1])
pimat$plot<-as.factor(pimat$plot)
rownames(pimat)<-rownames(hmscYe)

formdata <- as.HMSCdata(Y = hmscYg, X = hmscXe, Random=pimat,interceptX = F, scaleX=T)
formprior <- as.HMSCprior(formdata,family="gaussian") #not necessary, this just generates flat priors
formparam <- as.HMSCparam(formdata, formprior) #not necessary, this just generates random staring parameters

model <- hmsc(formdata, family = "gaussian", niter = 10000, nburn = 1000, thin = 10)

mixing <- as.mcmc(model,parameters = "paramX")
temp<-as.vector(mixing[1:900,3])
#hist(temp)
plot(temp,type="l")

### Convert the mixing object to a matrix
#mixingDF <- as.data.frame(mixing)
#boxplot(mixingDF[,2], las = 2)

#CI
average <- apply(model$results$estimation$paramX, 1:2, mean)
### 95% confidence intervals
CI.025 <- apply(model$results$estimation$paramX, 1:2, quantile, probs = 0.025)
CI.975 <- apply(model$results$estimation$paramX, 1:2, quantile, probs = 0.975)
CI <- cbind(as.vector(CI.025), as.vector(CI.975))

plot(0, 0, xlim = c(1, nrow(CI)), ylim = range(CI), type = "n", xlab = "", ylab = "", main="paramX")
abline(h = 0,col = "grey")
arrows(x0 = 1:nrow(CI), x1 = 1:nrow(CI), y0 = CI[, 1], y1 = CI[, 2], code = 3, angle = 90, length = 0.05)
points(1:nrow(CI), average, pch = 15, cex = 1.5)

#checking that intercepts make sense, yes more or less
bdenovo90996 intercept is .72
mean(hmscYe$bdenovo90996) is .81
summary(lm(hmscYg[,27]~decostand(hmscXe[,2],method="standardize"))) #yes the coefficients are very close, bdenovo90996

### Summary table
paramXCITable <- cbind(unlist(as.data.frame(average)),
                       unlist(as.data.frame(CI.025)),
                       unlist(as.data.frame(CI.975)))
colnames(paramXCITable) <- c("paramX", "lowerCI", "upperCI")
rownames(paramXCITable) <- paste(rep(colnames(average),
                                     each = nrow(average)), "_",
                                 rep(rownames(average),
                                     ncol(average)), sep="")


###### Variance partitioning######
#this code is ssentially going through and calculating the parameters*Xdata for each group. it is stll variance though so I'm not sure how it relates to R2
variationPart <- variPart(model,c(rep("abiotic",5),rep("biotic",2)))#"climate","habitat"
barplot(t(variationPart), legend.text=colnames(variationPart),
        args.legend=list(y=1.1, x=nrow(variationPart)/2, xjust=0.5, horiz=T))
#effect of abiotic:
mean(variationPart[,1])
#biotic:
mean(variationPart[,2])
#plot
mean(variationPart[,3])

colnames(hmscYe)
#lo: bact:1:256  fun:267:289   euk: 259:261  nem:262-267
#hi bact:1:273  fun:293:311  euk:274-278  nem:279-292 

###### R2 ######
#I don't understand the difference b/t variance partitioning and R2, the numbers are not the same. This calcualtes R2 <- 1 - ssRes/ssY. Ok I think the difference is that the total height of the bar (which is 100% in the variance partitining) should be scaled to this R2 value; because above is the contribution of each parameter to the predicted value and here is the difference b/t the predicted value and the observed value
Ymean <- apply(model$data$Y,2,mean)
R2 <- Rsquared(model, averageSp=FALSE)
mean(R2)
plot(Ymean,R2,pch=19)

#checking that R2 makes sense. this doesn't work b/c I guess R2 is relative to a model with only an intercept. but it is good b/c bdenovo990996 has high correlations with other taxa
bdenovo990996 R2 is 0.49353501
summary(lm(hmscYe$bdenovo90996~decostand(hmscXe[,2],method="standardize")))$r.squared #r2 is .25 but that doesn't include a random plot effect, and I'm not sure if it is relevant to varipart b/c th random+fixed effects R2 adds up to 100, so it is relative??
str(summary(lm(hmscYe$bdenovo90996~1)))
test=data.frame(y=hmscYe$bdenovo90996,x=decostand(hmscXe[,2],method="standardize"),r=1:29)

require(MuMIn)
#The marginal R squared values are those associated with your fixed effects, the conditional ones are those of your fixed effects plus the random effects. 
m1<-lme(y~x,random=~1|r,data=test)
r.squaredGLMM(m1)
#marginal (fixed only) R2: .244, conditional (f+r): .907
#these r2 are not matching up with the ones from hmsc, i can't figure out why, unless it is because the hmsc r2 consider the effect of interactions with other species??? I could try fitting the hmsc model with more or fewer species to see if the R2 changes


###### Species co-occurrances #####

##### Covariances #####
#Confidence intervals for pairwise covariances. Covariances mcmc are much more normally distributed, so I feel more comfortable using them for calculating z statsitics and p values. They aren't perfectly normally distributed (so it would probably be more accurrate to use some kind of on sided z-statistic) but it is probably ok and would be easier than explaining all that.
corMat <- corRandomEff(model,cor=F) #cor=F gives the covariance. These covariances are not actualy fit by the model, they are calcuated post hoc based on the latent variables at each iteration, start 3:05, end 3:15ish

#ltri <- lower.tri(apply(corMat[, , , 1], 1:2, quantile, probs = 0.025),diag=F)#start 3:20, end waited until 4:20 but it hadn't stopped so I stopped it #used to be diat=T
#apply function takes a ridiculousy long time, need to vectorize things
ltri<-lower.tri(array(NA,dim=c(dim(corMat)[1],dim(corMat)[2])),diag=F) #does the same thing and is way way faster

### Average
#averagec <- as.vector(apply(corMat[, , , 1], 1:2, mean)[ltri])
averagec<- as.vector(rowSums(corMat[,,,1],dims=2)[ltri])/dim(corMat)[3] #faster
medianc <- as.vector(apply(corMat[, , , 1], 1:2, median)[ltri])
head(cbind(averagec,medianc))
hist(corMat[1,2,,1]) #dims 1 and 2 are the matrix species (35) by species, 3 is the iteration, 4 is the plot level (I only have 1 level in that dimension), bdenovo32932 and bdenovo195709
plot(density(corMat[1,2,,1]))
### 95% confidence intervals
corMat.025 <- as.vector(apply(corMat[, , , 1], 1:2, quantile,probs = 0.025)[ltri]) 
corMat.975 <- as.vector(apply(corMat[, , , 1], 1:2, quantile,probs=0.975)[ltri])
#try something like this to make it faster???
#a<-array(data=1:90,dim=c(3,3,10,1))
#colSums(matrix(a,nrow=10,ncol=9,byrow=T))/dim(a)[3]

CICov <- cbind(corMat.025, corMat.975)
head(CICov)
# plot(0, 0, xlim = c(1, nrow(CI)), ylim = range(CI), type = "n", xlab = "", main = "cov(paramLatent[[1, 1]])")
# abline(h = 0, col = "grey")
# arrows(x0 = 1:nrow(CI), x1 = 1:nrow(CI), y0 = CI[, 1], y1 = CI[, 2], code = 3, angle = 90, length = 0.05)
# points(1:nrow(CI), average, pch = 15,cex = 1.5)

#put labels on the pairwise correlations
rownames(corMat)
CorSp<-data.frame(sp1=rep(NA,length(averagec)),sp2=rep(NA,length(averagec)))
r=1
for(i in 1:(length(rownames(corMat))-1)){
  for(j in (i+1):length(rownames(corMat))){
    CorSp[r,1]<-rownames(corMat)[i]
    CorSp[r,2]<-rownames(corMat)[j]
    r=r+1
  }
}
head(CorSp)
tail(CorSp)
CovsCI<-cbind(CorSp,averagec,CICov)
head(CovsCI)

ind<-which(CovsCI$corMat.025<0&CovsCI$corMat.975<0|CovsCI$corMat.025>0&CovsCI$corMat.975>0)
CovsCI[ind,]#there are 12 that don't overlap 0!!! need to see if my code below is correct
length(ind)
unique(c(CovsCI$sp1[ind],CovsCI$sp2[ind]))


###### Correlations #####
#Confidence intervals for pairwise correlations, these distributions are very odd, they are bimodal with a lot at -1 and a lot at 1, so the mean is somewhere in the middle
corMat <- corRandomEff(model, cor = TRUE)
hist(corMat[15,14 , , 1])
ltri <- lower.tri(apply(corMat[, , , 1], 1:2, quantile, probs = 0.025),diag=F) #originally diag=T
### Average
averageCor <- as.vector(apply(corMat[, , , 1], 1:2, mean)[ltri])
medianCor <- as.vector(apply(corMat[, , , 1], 1:2, median)[ltri])
### 95% confidence intervals
corMat.025 <- as.vector(apply(corMat[, , , 1], 1:2, quantile,probs = 0.025)[ltri])
corMat.975 <- as.vector(apply(corMat[, , , 1], 1:2, quantile,probs=0.975)[ltri])
CICor <- cbind(corMat.025, corMat.975)
head(CICor)
#which(CI[,1]>0&CI[,2]>0)

#put labels on the pairwise correlations
rownames(corMat)
CorSp<-data.frame(sp1=rep(NA,length(averageCor)),sp2=rep(NA,length(averageCor)))
r=1
for(i in 1:(length(rownames(corMat))-1)){
  for(j in (i+1):length(rownames(corMat))){
    CorSp[r,1]<-rownames(corMat)[i]
    CorSp[r,2]<-rownames(corMat)[j]
    r=r+1
  }
}
head(CorSp)
tail(CorSp)
CorsCI<-cbind(CorSp,averageCor,CICor)
head(CorsCI)

#thse are the same 12 as above with the Cov matrix, so I guess that's ok, even though the distributions are kind of wonky
ind<-which(CorsCI$corMat.025<0&CorsCI$corMat.975<0|CorsCI$corMat.025>0&CorsCI$corMat.975>0)
CorsCI[ind,]#there are 12 that don't overlap 0!!! need to see if my code below is correct
hist(CorsCI$averageCor)

ind<-which(CorsCI$averageCor>.7)
CorsCI[ind,]

hist(corMat["bdenovo193772","bdenovo188874" , , 1])
hist(corMat["bdenovo117183","bdenovo85595" , , 1])



#looking at some of the correlations
CorsCI3<-subset(CorsCI2,qval<.01)
plot(hmscYe[,"bdenovo193772"],hmscYe[,"bdenovo145037"])
abline(lm(hmscYe[,"bdenovo145037"]~hmscYe[,"bdenovo193772"]))
summary(lm(hmscYe[,"bdenovo145037"]~hmscYe[,"bdenovo193772"]))


###### Plotting matrix #####
averageCor2 <- apply(corMat[, , , 1], 1:2, mean) #this is the full matrix, not just the lower triangle, for the mean it might not matter, but for CI it might? but no I don't think it does for the CI, because the CI is taking the confidence interval over the mcmc chain, it just does the calcualtion twice if you do it on the whole matrix
corrplot(averageCor2, method = "color", col = colorRampPalette(c("blue", "white", "red"))(200))


###### Circle network ######
corMat2 <- corRandomEff(model, cor = TRUE)
averageCor <- apply(corMat2[, , , 1], 1:2, mean)
colMat <- matrix(NA, nrow = nrow(averageCor), ncol = ncol(averageCor))
colMat[which(averageCor > 0.7, arr.ind = TRUE)] <- "red"
colMat[which(averageCor < -0.7, arr.ind = TRUE)] <- "blue"
chordDiagram(averageCor, symmetric = TRUE,
             annotationTrack = c("name", "grid"),
             grid.col = "grey",col=colMat)

ltri<-lower.tri(array(NA,dim=c(dim(corMat2)[1],dim(corMat2)[2])),diag=F)
averageCor <- as.vector(rowSums(corMat2[,,,1],dims=2)[ltri])/dim(corMat2)[3] #faster but this is a vector not a matrix, but I can use it to count
length(which(averageCor>.7|averageCor<(-.7)))
length(which(averageCor>.8|averageCor<(-.8)))
length(which(averageCor>.9|averageCor<(-.9)))
which(averageCor<(-.7))
which(averageCor>(.7))
hist(averageCor)

library(tidyr)
library(dplyr)

averageCor2<-data.frame(taxon1=rownames(averageCor),averageCor)
averageCor3<-gather(averageCor2,taxon2,cor,bdenovo195709:idenovo18189) #not labeling like I want
head(averageCor3)
subset(averageCor3,averageCor3$cor>.9&averageCor3$cor<1)

#checking, for species 2 and 3 intercept and slope should be significant
plot(formdata$X[,2],formdata$Y[,2])
abline(a=-2.571260,b=0.5421437)
abline(a=-7.1447,b=1.4906,col=2)
summary(lm(formdata$Y[,1]~formdata$X[,2]))
summary(lm(formdata$Y[,4]~formdata$X[,2]))


#trying just calculating intercept on simulated data, everything checks out, confidence intervals are correct and standard error in lm is the same as the calucated SE and you can calculate the CI from the standard error in lm()
hmscYd<-as.matrix(data.frame(y1=rnorm(1000,mean=1,sd=.2),y2=rnorm(1000,mean=1,sd=.2)))
hmscXd<-as.matrix(rep(1,1000))

formdata <- as.HMSCdata(Y = hmscYd, X = hmscXd, interceptX = F, scaleX=T)
model <- hmsc(formdata, family = "gaussian", niter = 10000, nburn = 1000, thin = 10)
0.99505+1.96*.2/sqrt(1000)
summary(lm(hmscYd[,1]~1))
.2/sqrt(1000)
sd(hmscYd)/sqrt(1000)
0.995050-0.006077*1.96
0.995050+0.006077*1.96

#trying to calculate likelihood
-sum(dnorm(formdata$Y[,2],mean=-100+1.4906*formdata$X[,2],sd=1,log=T)) 
-sum(dnorm(formdata$Y[,2],mean=-7.1447+1.4906*formdata$X[,2],sd=1,log=T)) #I don't know the sd, but I think it doesn't matter if I'm only looking at relative differences in loglik
-sum(dnorm(formdata$Y[,2],mean=-2.5712604+0.5421437*formdata$X[,2],sd=1,log=T))
#in terms of likelihood, the lm() model estimates are better




#upshot
#mcmc is very sensitive to the range in the y variables, if the numbers are too low, the CI will be huge b/c somehow it takes huge jumps in the mcmc - solution: scale the Y
#mcmc is also sensitive to the absolute value of the x variables (but I did this on my original data trial so this was not the problem then). even if it ranges from 6-8 the mcmc gives non-optimal results compared if the range is -1 to 1 








#####Old code#####

# Covariances
#translate to p values
#P from CI for a difference
#If the upper and lower limits of a 95% CI are u and l respectively:
#  1 calculate the standard error: SE = (u − l)/(2*1.96)
#  2 calculate the test statistic: z = Est/SE
#  3 calculate the P value2: P = exp(−0.717×z − 0.416×z2). #I changed to using pnorm(), more accurate
#Even with a model with only an intercept, the lowest p value is .001, which is a qval of .8, so nothing is significant
#Thoughts - since these are all estimated simulteneously, I don't know if I really need to account for multiple comparions. Also, I could mayb just use the cutoff "significance is when th 95% CI don't cross zero" this is the bayesian credible interval, right? yes, and a 95% CI should correspond to a 0.05 p value
CovsCI$SE<-(CovsCI$corMat.975-CovsCI$corMat.025)/(2*1.96)
CovsCI$z<-CovsCI$averagec/CovsCI$SE
CovsCI$absz<-abs(CovsCI$z)
#CovsCI$Pone<-pnorm(-CovsCI$absz) #one tailed test, don't use
CovsCI$P<-2*pnorm(-CovsCI$absz) #two tailed test. it should be a two tailed test b/c it could be higher or lower than 0 (and that's what lme does for example)
#CovsCI$Pold<-exp(-.717*CovsCI$absz-.416*CovsCI$absz^2)
which(CovsCI$Ptwo<.1)
CovsCI$qval<-p.adjust(CovsCI$P,method="fdr")#
sort(CovsCI$qval)
sort(CovsCI$P) 
CovsCI2<-subset(CovsCI,averagec>0)
which(CovsCI2$P<.1)
CovsCI2[37,]#8412,3112
CovsCI2[135,]#8412,3112

CovsCI[which(CovsCI$sp1=="bdenovo193772"&CovsCI$sp2=="bdenovo20595"),]
hist(corMat["bdenovo193772","bdenovo20595",,1])
hist(corMat["bdenovo117183","bdenovo85595",,1])
#so I think my p values are high b/c I'm dealing with slighly non normal data

#looking into one species pair with strong positive covariance
bdenovo88234 bdenovo137544
bdenovo184998 bdenovo117183
bdenovo198139 bdenovo51656
bdenovo193772 bdenovo20595
hmscXe
hmscYe[,"bdenovo193772"]
hmscYe[,"bdenovo20595"]
cov(hmscYe[,"bdenovo193772"],hmscYe[,"bdenovo20595"])
cor(hmscYe[,"bdenovo193772"],hmscYe[,"bdenovo20595"])
#library(Hmisc)
#rcorr(cbind(hmscYe[,"bdenovo193772"],hmscYe[,"bdenovo20595"]))$P #same P value as regression model
plot(hmscYe[,"bdenovo193772"],hmscYe[,"bdenovo20595"])
summary(lm(hmscYe[,"bdenovo193772"]~hmscYe[,"bdenovo20595"]))
abline(lm(hmscYe[,"bdenovo193772"]~hmscYe[,"bdenovo20595"]))
lm1<-lm(hmscYe[,"bdenovo193772"]~hmscXe[,"pH"])
lm2<-lm(hmscYe[,"bdenovo20595"]~hmscXe[,"pH"])
plot(resid(lm1),resid(lm2))
abline(lm(resid(lm1)~resid(lm2)))
summary(lm(resid(lm1)~resid(lm2)))





# correlations
#translate to p values
#P from CI for a difference
#If the upper and lower limits of a 95% CI are u and l respectively:
#  1 calculate the standard error: SE = (u − l)/(2*1.96)
#  2 calculate the test statistic: z = Est/SE
#  3 calculate the P value2: P = exp(−0.717×z − 0.416×z2).
#This only works when the CIs are symmetric around the estimate, as shown below, many of my CIs are not asymmetric b/c correlations can't be larger than 1 so the CI is squinched at that side
CorsCI$SE<-(CorsCI$corMat.975-CorsCI$corMat.025)/(2*1.96)
CorsCI$z<-CorsCI$averageCor/CorsCI$SE
CorsCI$absz<-abs(CorsCI$z)
CorsCI$P<-exp(-.717*CorsCI$absz-.416*CorsCI$absz^2)
which(CorsCI$P<.001)
CorsCI$qval<-p.adjust(CorsCI$P,method="fdr")
CorsCI2<-subset(CorsCI,averageCor>0)

#trying it by taking the "one-sided" p value
CorsCI$SE2<-ifelse(CorsCI$averageCor>0,(CorsCI$averageCor-CorsCI$corMat.025)/(1*1.96),(CorsCI$corMat.975-CorsCI$averageCor)/(1*1.96)) #this way is still doing a "two tailed test" b/c I'm using 1.96 sd away from the mean, if I were doing a one tailed test I would use 1.645 sd. In other words this way is like doing a one-tailed test with an alpha=0.025, and I'm using a one-tailed test simply b/c I don't have the proper data (i.e. asymmetrical CIs) to calcualted a true two tailed test. however, it seems that the z statistic calculation is based on the type of confidence interval you have, so if you have 95% CI, then you should use 1.96, it doesnt make sense to use a different number. I think the alpha cutoff is for translating the z value to a pvalue and giving you a cut off. Another thought is that I maybe should be using the median instead of the mean here 
CorsCI$z2<-CorsCI$averageCor/CorsCI$SE2
CorsCI$absz<-abs(CorsCI$z)
CorsCI$P<-exp(-.717*CorsCI$absz-.416*CorsCI$absz^2)
which(CorsCI$P<.001)
CorsCI$qval<-p.adjust(CorsCI$P,method="fdr")
CorsCI2<-subset(CorsCI,averageCor>0)


length(which(CorsCI2$P<.05))
length(which(CorsCI2$qval<.01))
hist(CorsCI2$averageCor[which(CorsCI2$qval<.01)])







# Mixing object, I don't think this is useful. I don't totally understand latent variables- they somehow help with fitting a model with so many parameters (the covoariance matrix has 630 parameters in this example case accounting for all the species pairwise correaltions.) when you run a model with no env factors, the latent variables are essentially ordination axes
mixing <- as.mcmc(model, parameters = "paramLatent")
### Draw trace and density plots for all combination of paramters
dim(mixing) #900 x 140, I don't undrsand what this is, what the 140 refers to, it looks like there are about 3-5 latent vaiables per species
plot(mixing[,140])
mean(mixing[,1])
hist(mixing[,1])
cbind(corMat[1,2,,1],mixing[,1],corMat[1,2,,1]/mixing[,1])
### Convert the mixing object to a matrix
mixingDF <- as.data.frame(mixing[[2]])
### Draw boxplot for each parameters
par(mar = c(7, 4, 4, 2))
boxplot(mixingDF, las = 2)
### Draw beanplots
library(beanplot)
par(mar = c(7, 4, 4, 2))
beanplot(mixingDF, las = 2)

