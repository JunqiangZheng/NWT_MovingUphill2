#HMSCtutorial

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
plot(mixing[[1]]) #this doesn't work, it just extracts the first number
plot(mixing[,41])
### Convert the mixing object to a matrix
mixingDF <- as.data.frame(mixing[,13])
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


