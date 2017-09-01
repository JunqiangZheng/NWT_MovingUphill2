#HMSC

# To install HMSC: 
#   install Xcode
#   make sure command line tools are installed by typing "make" in the terminal
#   install gfortran 6.3 for sierra
#   follow instructions on http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/ for doing the following in terminal:
#   mkdir ~/.R
#   cat << EOF >> ~/.R/Makevars
#   FLIBS=-L/usr/local/gfortran/lib/gcc/x86_64-apple-darwin16/6.3.0 -L/usr/local/gfortran/lib -lgfortran -lquadmath -lm
#   EOF
#   install dependencies and HMSC follwoing GitHub instructions: https://github.com/guiblanchet/HMSC

library(HMSC)

data("simulEx1")
data("simulParamEx1") #true parameters
model <- hmsc(simulEx1, family = "probit", niter = 10000, nburn = 1000, thin = 10)

#save results outside of the work env, but it just saves it as an R workspace
save(model, file = "/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/HMSC/SimulatedExample1/model.RData")
load(file = "/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/HMSC/SimulatedExample1/model.RData")

mixing <- as.mcmc(model, parameters = "paramX")
plot(mixing)

### Convert the mixing object to a matrix
mixingDF <- as.data.frame(mixing)
### Draw beanplots
library(beanplot)
par(mar = c(7, 4, 4, 2))
beanplot(mixingDF, las = 2)
### Draw boxplot for each parameters
par(mar = c(7, 4, 4, 2))
boxplot(mixingDF, las = 2)
### True values
truth <- as.vector(simulParamEx1$param$paramX)
### Average
average <- apply(model$results$estimation$paramX, 1:2, mean)
### 95% confidence intervals
CI.025 <- apply(model$results$estimation$paramX, 1:2, quantile,
                probs = 0.025)
CI.975 <- apply(model$results$estimation$paramX, 1:2, quantile,
                probs = 0.975)
CI <- cbind(as.vector(CI.025), as.vector(CI.975))
### Draw confidence interval plots
plot(0, 0, xlim = c(1, nrow(CI)), ylim = range(CI, truth), type = "n",
     xlab = "", ylab = "", main="paramX")
abline(h = 0,col = "grey")
arrows(x0 = 1:nrow(CI), x1 = 1:nrow(CI), y0 = CI[, 1], y1 = CI[, 2],
       code = 3, angle = 90, length = 0.05)
points(1:nrow(CI), average, pch = 15, cex = 1.5)
points(1:nrow(CI), truth, col = "red", pch = 19)

paramXCITable <- cbind(unlist(as.data.frame(average)),
                       unlist(as.data.frame(CI.025)),
                       unlist(as.data.frame(CI.975)))
colnames(paramXCITable) <- c("paramX", "lowerCI", "upperCI")
rownames(paramXCITable) <- paste(rep(colnames(average),
                                     each = nrow(average)), "_",
                                 rep(rownames(average),
                                     ncol(average)), sep="")
#variance partitioning
variationPart <- variPart(model,c(rep("climate",2),"habitat"))
barplot(t(variationPart), legend.text=colnames(variationPart),
        args.legend=list(y=1.1, x=nrow(variationPart)/2, xjust=0.5, horiz=T))

#association matrix
### Plot random effect estimation through correlation matrix
corMat <- corRandomEff(model,cor=FALSE)
#--------------------------------------------------------------------------
  ### Sampling units level
#----------------------------------------------------------------------------
### Isolate the values of interest
ltri <- lower.tri(apply(corMat[, , , 1], 1:2, quantile, probs = 0.025),diag=TRUE)
### True values
truth <- as.vector(tcrossprod(simulParamEx1$param$paramLatent[[1]])[ltri])
### Average
average <- as.vector(apply(corMat[, , , 1], 1:2, mean)[ltri])
### 95% confidence intervals
corMat.025 <- as.vector(apply(corMat[, , , 1], 1:2, quantile,
                              probs = 0.025)[ltri])
corMat.975 <- as.vector(apply(corMat[, , , 1], 1:2, quantile,
                              probs=0.975)[ltri])
CI <- cbind(corMat.025, corMat.975)


### Plot the results
plot(0, 0, xlim = c(1, nrow(CI)), ylim = range(CI, truth), type = "n",
     xlab = "", main = "cov(paramLatent[[1, 1]])")
abline(h = 0, col = "grey")
arrows(x0 = 1:nrow(CI), x1 = 1:nrow(CI), y0 = CI[, 1], y1 = CI[, 2],
       code = 3, angle = 90, length = 0.05)
points(1:nrow(CI), average, pch = 15,cex = 1.5)
points(1:nrow(CI), truth, col = "red", pch=19)
### Mixing object
mixing <- as.mcmc(model, parameters = "paramLatent")
### Draw trace and density plots for all combination of parameters
plot(mixing[[1]])
### Convert the mixing object to a matrix
mixingDF <- as.data.frame(mixing[[1]])
### Draw boxplot for each parameters
par(mar=c(7, 4, 4, 2))
boxplot(mixingDF, las = 2)
### Draw beanplots
library(beanplot)
par(mar = c(7, 4, 4, 2))
beanplot(mixingDF, las = 2)

### Draw estimated correlation matrix
library(corrplot)
corMat <- corRandomEff(model, cor = TRUE)
averageCor <- apply(corMat[, , , 1], 1:2, mean)
corrplot(averageCor, method = "color",col = colorRampPalette(c("blue", "white", "red"))(200))
### Draw chord diagram
library(circlize)
corMat <- corRandomEff(model, cor = TRUE)
averageCor <- apply(corMat[, , , 1], 1:2, mean)
colMat <- matrix(NA, nrow = nrow(averageCor), ncol = ncol(averageCor))
colMat[which(averageCor > 0.4, arr.ind = TRUE)] <- "red"
colMat[which(averageCor < -0.4, arr.ind = TRUE)] <- "blue"
chordDiagram(averageCor, symmetric = TRUE,
             annotationTrack = c("name", "grid"),
             grid.col = "grey",col=colMat)

#----------------------------------------------------------------------------
### Plot level
#----------------------------------------------------------------------------
### Isolate the values of interest
ltri <- lower.tri(apply(corMat[, , , 2], 1:2, quantile, probs=0.025),diag=TRUE)
### True values
truth <- as.vector(tcrossprod(simulParamEx1$param$paramLatent[[2]])[ltri])
### Average
average <- as.vector(apply(corMat[, , , 2], 1:2, mean)[ltri])
### 95% confidence intervals
corMat.025 <- as.vector(apply(corMat[, , , 2], 1:2, quantile,
                              probs = 0.025)[ltri])
corMat.975 <- as.vector(apply(corMat[, , , 2], 1:2, quantile,
                              probs = 0.975)[ltri])
CI <- cbind(corMat.025, corMat.975)
### Plot the results
plot(0, 0, xlim = c(1, nrow(CI)), ylim = range(CI, truth), type = "n",
     xlab = "", main = "cov(paramLatent[[1,2]])")
abline(h = 0, col = "grey")
arrows(x0 = 1:nrow(CI), x1 = 1:nrow(CI), y0 = CI[, 1], y1 = CI[, 2],
       code = 3, angle = 90, length = 0.05)
points(1:nrow(CI), average, pch = 15, cex = 1.5)
points(1:nrow(CI), truth, col = "red", pch = 19)
### Mixing object
mixing <- as.mcmc(model, parameters = "paramLatent")
### Draw trace and density plots for all combination of paramters
plot(mixing[[2]])

### Convert the mixing object to a matrix
mixingDF <- as.data.frame(mixing[[2]])
### Draw boxplot for each parameters
par(mar = c(7, 4, 4, 2))
boxplot(mixingDF, las = 2)
### Draw beanplots
library(beanplot)
par(mar = c(7, 4, 4, 2))
beanplot(mixingDF, las = 2)
### Draw estimated correlation matrix
library(corrplot)
corMat <- corRandomEff(model, cor = TRUE)
averageCor <- apply(corMat[, , , 2], 1:2, mean)
corrplot(averageCor, method = "color",
         col = colorRampPalette(c("blue", "white", "red"))(200))

### Draw chord diagram
library(circlize)
corMat <- corRandomEff(model, cor = TRUE)
averageCor <- apply(corMat[, , , 2], 1:2, mean)
colMat <- matrix(NA, nrow = nrow(averageCor), ncol = ncol(averageCor))
colMat[which(averageCor>0.4, arr.ind = TRUE)] <- "red"
colMat[which(averageCor< -0.4, arr.ind = TRUE)] <- "blue"
chordDiagram(averageCor, symmetric = TRUE,
             annotationTrack = c("name", "grid"), grid.col = "grey",
             col = colMat)

# R squared
Ymean <- apply(model$data$Y,2,mean)
R2 <- Rsquared(model, averageSp=FALSE)
plot(Ymean,R2,pch=19)


### For a particular set of parameter, this given directly by the hmsc function. For example for the beta (paramX)
model$results$estimation$paramX
### For the full joint probability distribution
fullPost <- jposterior(model)

#training data
predTrain <- predict(model)

#----------------------------------------------------------------------------
### Simulating "validation" data
#----------------------------------------------------------------------------
X <- matrix(nrow = 10, ncol = 3)
colnames(X) <- colnames(simulEx1$X)
X[, 1] <- 1
X[, 2] <- rnorm(10)
X[, 3] <- rnorm(10)
RandomSel <- sample(200, 10)
Random <- simulEx1$Random[RandomSel, ]
for(i in 1:ncol(Random)){
  Random[, i] <- as.factor(as.character(Random[, i]))
}
colnames(Random) <- colnames(simulEx1$Random)
dataVal <- as.HMSCdata(X = X, Random = Random)
#----------------------------------------------------------------------------
### Prediction for a new set of values
#----------------------------------------------------------------------------
predVal <- predict(model, dataVal)







#for real data

comm.dataALL

hmscY<-comm.dataALL[,32:6704]
hmscY<-comm.dataALL[,32:1000] #for practice
rownames(hmscY)<-comm.dataALL$X.SampleID

biogeo3<-biogeo2[order(biogeo2$Sample_name),]#now they are in the same order
hmscX<-data.frame(inter=rep(1,90),pH=biogeo3$pH.x,lomehi=biogeo3$lomehi) #there are NAs in the pH data and possibly more, moisture=biogeo3$moisture.x,TN=biogeo3$TN,TC=biogeo3$TC,
rownames(hmscX)<-biogeo3$X.SampleID

#take out any rows with NAs
ind<-which(!is.na(rowSums(hmscX[,1:2])))
hmscXb<-hmscX[ind,]
hmscYb<-hmscY[ind,]

#readjust the lo me hi categories based on NAs and # samples in each category
#get plant cover data to use as "light"

#select lo
ind<-which(hmscXb$lomehi=="lo")
hmscXc<-hmscXb[ind,]
hmscYc<-hmscYb[ind,]

#select species with greater than 10 (11 or more) and remove lo me hi (since you can't have text in a matrix or the whole matrix becomes character)
ind<-which(colSums(hmscYc>0)>10)
hmscYd<-hmscYc[,ind]
hmscXd<-hmscXc[,1:2]

#testing a really strong correlation with pH
#hmscYd[,2]<-rnorm(29,mean=2*formdata$X[,2]+1,sd=.02)

#select only species 2 and 3
#hmscYd<-cbind(hmscYd[,3],hmscYd[,3])

#check if the values are too low that some tolerance is messing up the CI estimates, yes important to scale y
hmscYe<-scale(hmscYd)

#make them matrices
hmscXe<-as.matrix(hmscXd)
hmscYe<-as.matrix(hmscYe)

#the y data are not normal (the only options I have are normal, binary, poisson, overdispersed poisson), so I could do a sqrt transformation on Y (log(0) is -Inf)
#hmscYf<-sqrt(hmscYe)

#make the random (residual matrix)
pimat<-data.frame(plot=1:dim(hmscYe)[1])
pimat$plot<-as.factor(pimat$plot)
rownames(pimat)<-rownames(hmscYe)

formdata <- as.HMSCdata(Y = hmscYe, X = hmscXe, Random=pimat,interceptX = F, scaleX=T)
#formdata <- as.HMSCdata(Y = hmscYe, X = hmscXe,interceptX = F, scaleX=T)
formprior <- as.HMSCprior(formdata,family="gaussian") #not necessary, this just generates flat priors
formparam <- as.HMSCparam(formdata, formprior) #not necessary, this just generates random staring parameters

model <- hmsc(formdata, family = "gaussian", niter = 10000, nburn = 1000, thin = 10)

mixing <- as.mcmc(model,parameters = "paramX")
temp<-as.vector(mixing[1:900,1])
hist(temp)
plot(temp,type="l")

### Convert the mixing object to a matrix
mixingDF <- as.data.frame(mixing)
boxplot(mixingDF[,2], las = 2)

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


### Summary table
paramXCITable <- cbind(unlist(as.data.frame(average)),
                       unlist(as.data.frame(CI.025)),
                       unlist(as.data.frame(CI.975)))
colnames(paramXCITable) <- c("paramX", "lowerCI", "upperCI")
rownames(paramXCITable) <- paste(rep(colnames(average),
                                     each = nrow(average)), "_",
                                 rep(rownames(average),
                                     ncol(average)), sep="")

#variance partitioning
variationPart <- variPart(model,c("climate","habitat"))
barplot(t(variationPart), legend.text=colnames(variationPart),
        args.legend=list(y=1.1, x=nrow(variationPart)/2, xjust=0.5, horiz=T))

#R2
Ymean <- apply(model$data$Y,2,mean)
R2 <- Rsquared(model, averageSp=FALSE)
plot(Ymean,R2,pch=19)

#species correlations
library(corrplot)

corMat <- corRandomEff(model,cor=FALSE) #cor=F gives the covariance
ltri <- lower.tri(apply(corMat[, , , 1], 1:2, quantile, probs = 0.025),diag=TRUE)
### Average
average <- as.vector(apply(corMat[, , , 1], 1:2, mean)[ltri])
### 95% confidence intervals
corMat.025 <- as.vector(apply(corMat[, , , 1], 1:2, quantile,probs = 0.025)[ltri])
corMat.975 <- as.vector(apply(corMat[, , , 1], 1:2, quantile,probs=0.975)[ltri])
CI <- cbind(corMat.025, corMat.975)
plot(0, 0, xlim = c(1, nrow(CI)), ylim = range(CI), type = "n", xlab = "", main = "cov(paramLatent[[1, 1]])")
abline(h = 0, col = "grey")
arrows(x0 = 1:nrow(CI), x1 = 1:nrow(CI), y0 = CI[, 1], y1 = CI[, 2], code = 3, angle = 90, length = 0.05)
points(1:nrow(CI), average, pch = 15,cex = 1.5)

corMat <- corRandomEff(model, cor = TRUE)
averageCor <- apply(corMat[, , , 1], 1:2, mean)
corrplot(averageCor, method = "color", col = colorRampPalette(c("blue", "white", "red"))(200))

(a) P from CI for a difference
If the upper and lower limits of a 95% CI are u and l respectively:
1 calculate the standard error: SE = (u − l)/(2*1.96)
2 calculate the test statistic: z = Est/SE
3 calculate the P value2: P = exp(−0.717×z − 0.416×z2).

library(circlize)
corMat <- corRandomEff(model, cor = TRUE)
averageCor <- apply(corMat[, , , 1], 1:2, mean)
colMat <- matrix(NA, nrow = nrow(averageCor), ncol = ncol(averageCor))
colMat[which(averageCor > 0.6, arr.ind = TRUE)] <- "red"
colMat[which(averageCor < -0.6, arr.ind = TRUE)] <- "blue"
chordDiagram(averageCor, symmetric = TRUE,
             annotationTrack = c("name", "grid"),
             grid.col = "grey",col=colMat)



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