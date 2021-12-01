library(ggplot2)
library(ggpubr)
library(ks)
library(rjags)
library(runjags)
library(dplyr)
library(tidyr)
library(Hmisc)
setwd("F:\\Subbu\\RMIT\\sem 2\\Applied Bayesian Statistics\\Final Project\\")
source("DBDA2E-utilities.R") 


#===============PRELIMINARY FUNCTIONS FOR POSTERIOR INFERENCES====================

smryMCMC_HD = function(  codaSamples , compVal = NULL,  saveName=NULL) {
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  paramName = colnames(mcmcMat)
  for ( pName in paramName ) {
    if (pName %in% colnames(compVal)){
      if (!is.na(compVal[pName])) {
        summaryInfo = rbind( summaryInfo , summarizePost( paramSampleVec = mcmcMat[,pName] , 
                                                          compVal = as.numeric(compVal[pName]) ))
      }
      else {
        summaryInfo = rbind( summaryInfo , summarizePost( paramSampleVec = mcmcMat[,pName] ) )
      }
    } else {
      summaryInfo = rbind( summaryInfo , summarizePost( paramSampleVec = mcmcMat[,pName] ) )
    }
  }
  rownames(summaryInfo) = paramName
  
  # summaryInfo = rbind( summaryInfo , 
  #                      "tau" = summarizePost( mcmcMat[,"tau"] ) )
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}

#===============================================================================

plotMCMC_HD = function( codaSamples , data , xName="x" , yName="y", preds = FALSE ,
                        showCurve=FALSE ,  pairsPlot=FALSE , compVal = NULL, 
                        saveName=NULL , saveType="jpg" ) {
   #-----------------------------------------------------------------------------
  y = data[,yName]
  x = as.matrix(data[,xName])
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  beta0 = mcmcMat[,"beta0"]
  beta  = mcmcMat[,grep("^beta$|^beta\\[",colnames(mcmcMat))]
  if ( ncol(x)==1 ) { beta = matrix( beta , ncol=1 ) }
  if (preds){
    pred = mcmcMat[,grep("^pred$|^pred\\[",colnames(mcmcMat))]
  }
  guess = mcmcMat[,"guess"]
  #-----------------------------------------------------------------------------
  # Compute R^2 for credible parameters:
  YcorX = cor( y , x ) # correlation of y with each x predictor
  Rsq = beta %*% matrix( YcorX , ncol=1 )
  #-----------------------------------------------------------------------------
  if ( pairsPlot ) {
    # Plot the parameters pairwise, to see correlations:
    openGraph()
    nPtToPlot = 1000
    plotIdx = floor(seq(1,chainLength,by=chainLength/nPtToPlot))
    panel.cor = function(x, y, digits=2, prefix="", cex.cor, ...) {
      usr = par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      r = (cor(x, y))
      txt = format(c(r, 0.123456789), digits=digits)[1]
      txt = paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
      text(0.5, 0.5, txt, cex=1.25 ) # was cex=cex.cor*r
    }
    pairs( cbind( beta0 , beta , tau )[plotIdx,] ,
           labels=c( "beta[0]" , 
                     paste0("beta[",1:ncol(beta),"]\n",xName) , 
                     expression(tau) ) , 
           lower.panel=panel.cor , col="skyblue" )
    if ( !is.null(saveName) ) {
      saveGraph( file=paste(saveName,"PostPairs",sep=""), type=saveType)
    }
  }
  #-----------------------------------------------------------------------------
  # Marginal histograms:
  
  decideOpenGraph = function( panelCount , saveName , finished=FALSE , 
                              nRow=5 , nCol=4 ) {
    # If finishing a set:
    if ( finished==TRUE ) {
      if ( !is.null(saveName) ) {
        saveGraph( file=paste0(saveName,ceiling((panelCount-1)/(nRow*nCol))), 
                   type=saveType)
      }
      panelCount = 1 # re-set panelCount
      return(panelCount)
    } else {
      # If this is first panel of a graph:
      if ( ( panelCount %% (nRow*nCol) ) == 1 ) {
        # If previous graph was open, save previous one:
        if ( panelCount>1 & !is.null(saveName) ) {
          saveGraph( file=paste0(saveName,(panelCount%/%(nRow*nCol))), 
                     type=saveType)
        }
        # Open new graph
        openGraph(width=nCol*7.0/3,height=nRow*2.0)
        layout( matrix( 1:(nRow*nCol) , nrow=nRow, byrow=TRUE ) )
        par( mar=c(4,4,2.5,0.5) , mgp=c(2.5,0.7,0) )
      }
      # Increment and return panel count:
      panelCount = panelCount+1
      return(panelCount)
    }
  }
  
  # Original scale:
  panelCount = 1
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( beta0 , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(beta[0]) , main="Intercept", compVal = as.numeric(compVal["beta0"] ))
  for ( bIdx in 1:ncol(beta) ) {
    panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
    if (!is.na(compVal[paste0("beta[",bIdx,"]")])){
      histInfo = plotPost( beta[,bIdx] , cex.lab = 1.75 , showCurve=showCurve ,
                           xlab=bquote(beta[.(bIdx)]) , main=xName[bIdx],
                           compVal = as.numeric(compVal[paste0("beta[",bIdx,"]")]))
    } else{
      histInfo = plotPost( beta[,bIdx] , cex.lab = 1.75 , showCurve=showCurve ,
                           xlab=bquote(beta[.(bIdx)]) , main=xName[bIdx])
    }
  }
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( Rsq , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(R^2) , main=paste("Prop Var Accntd") , finished=FALSE )
  
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( guess , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab="Guess parameter" , main=paste("Prop Var Accntd") , finished=TRUE )
  
  panelCount = 1
  if ( pred){
    
    for ( pIdx in 1:ncol(pred) ) {
      panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
      histInfo = plotPost( pred[,pIdx] , cex.lab = 1.75 , showCurve=showCurve ,
                           xlab=bquote(pred[.(pIdx)]) , main=paste0("Prediction ",pIdx) ) 
    }
  }
  panelCount = 1
 
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMargZ") )
  histInfo = plotPost( Rsq , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(R^2) , main=paste("Prop Var Accntd") )
  panelCount = decideOpenGraph( panelCount , finished=TRUE , saveName=paste0(saveName,"PostMargZ") )
  
  #-----------------------------------------------------------------------------
}

#===============PRELIMINARY FUNCTIONS FOR POSTERIOR INFERENCES====================


assign_ds <- read.csv("weatherAUS.csv")
assign_ds <- filter(assign_ds, Location == "Melbourne")
assign_ds <- assign_ds[,-23] #Removing RISK_MM variable

colSums(is.na(assign_ds))
str(assign_ds)
summary(assign_ds$RainTomorrow)
summary(assign_ds)
 
summary(assign_ds)

# Relace NAs with mean+error
assign_ds$MinTemp[is.na(assign_ds$MinTemp)] = round(mean(assign_ds$MinTemp, na.rm = TRUE),2)
assign_ds$MaxTemp[is.na(assign_ds$MaxTemp)] = round(mean(assign_ds$MaxTemp, na.rm = TRUE),2)

#rainfall - replacing NA values with normally distributed values with variance 1
mean_rainfall <- mean(assign_ds$Rainfall,na.rm = TRUE)

assign_ds[which(is.na(assign_ds[,'Rainfall'])),'Rainfall'] <- rnorm(length(which(is.na(assign_ds[,'Rainfall']))),mean_rainfall,1)

#sunshine - replacing NA values with normally distributed values with variance 1
assign_ds$Sunshine[is.na(assign_ds$Sunshine)] = round(mean(assign_ds$Sunshine, na.rm = TRUE),1)

#WindGustDir
assign_ds$WindGustDir <- impute(assign_ds$WindGustDir,fun = mean)

#WindGustSpeed - replacing NA values with normally distributed values with variance 1
mean_wgspeed <- mean(assign_ds$WindGustSpeed,na.rm = TRUE)

assign_ds[which(is.na(assign_ds[,'WindGustSpeed'])),'WindGustSpeed'] <- rnorm(length(which(is.na(assign_ds[,'WindGustSpeed']))),mean_wgspeed,5)

#WindDir9am
assign_ds$WindDir9am <- impute(assign_ds$WindDir9am, fun = median)

#WindDir3pm
assign_ds$WindDir3pm <- impute(assign_ds$WindDir3pm, fun = median)

#WindSpeed9am - replacing NA values with normally distributed values with variance 1
mean_wspd9 <- mean(assign_ds$WindSpeed9am,na.rm = TRUE)

assign_ds[which(is.na(assign_ds[,'WindSpeed9am'])),'WindSpeed9am'] <- rnorm(length(which(is.na(assign_ds[,'WindSpeed9am']))),mean_wspd9,2)

#WindSpeed3pm - replacing NA values with normally distributed values with variance 1
mean_wspd3 <- mean(assign_ds$WindSpeed3pm,na.rm = TRUE)

assign_ds[which(is.na(assign_ds[,'WindSpeed3pm'])),'WindSpeed3pm'] <- rnorm(length(which(is.na(assign_ds[,'WindSpeed3pm']))),mean_wspd3,2)

#Humidity9am - replacing NA values with normally distributed values with variance 1
mean_humid9 <- mean(assign_ds$Humidity9am,na.rm = TRUE)

assign_ds[which(is.na(assign_ds[,'Humidity9am'])),'Humidity9am'] <- rnorm(length(which(is.na(assign_ds[,'Humidity9am']))),mean_humid9,5)

#Humidity3pm - replacing NA values with normally distributed values with variance 1
mean_humid3 <- mean(assign_ds$Humidity3pm,na.rm = TRUE)

assign_ds[which(is.na(assign_ds[,'Humidity3pm'])),'Humidity3pm'] <- rnorm(length(which(is.na(assign_ds[,'Humidity3pm']))),mean_humid3,5)



#pressure9am - replacing NA values with normally distributed values with variance 1
mean_p9 <- mean(assign_ds$Pressure9am,na.rm = TRUE)
assign_ds[which(is.na(assign_ds[,'Pressure9am'])),'Pressure9am'] <- rnorm(length(which(is.na(assign_ds[,'Pressure9am']))),mean_p9,1)


#pressure3pm - replacing NA values with normally distributed values with variance 1
mean_p3 <- mean(assign_ds$Pressure3pm,na.rm = TRUE)

assign_ds[which(is.na(assign_ds[,'Pressure3pm'])),'Pressure3pm'] <- rnorm(length(which(is.na(assign_ds[,'Pressure3pm']))),mean_p3,1)


#cloud9am - replacing NA values with normally distributed values with variance 1
mean_c9 <- mean(assign_ds$Cloud9am,na.rm = TRUE)

assign_ds[which(is.na(assign_ds[,'Cloud9am'])),'Cloud9am'] <- rnorm(length(which(is.na(assign_ds[,'Cloud9am']))),mean_c9,1)


#cloud3pm - replacing NA values with normally distributed values with variance 1
mean_c3 <- mean(assign_ds$Cloud3pm,na.rm = TRUE)

assign_ds[which(is.na(assign_ds[,'Cloud3pm'])),'Cloud3pm'] <- rnorm(length(which(is.na(assign_ds[,'Cloud3pm']))),mean_c3,1)


#temp9am - replacing NA values with normally distributed values with variance 1
mean_t9 <- mean(assign_ds$Temp9am,na.rm = TRUE)

assign_ds[which(is.na(assign_ds[,'Temp9am'])),'Temp9am'] <- rnorm(length(which(is.na(assign_ds[,'Temp9am']))),mean_t9,1)


#temp3pm - replacing NA values with normally distributed values with variance 1
mean_t3 <- mean(assign_ds$Temp3pm,na.rm = TRUE)

assign_ds[which(is.na(assign_ds[,'Temp3pm'])),'Temp3pm'] <- rnorm(length(which(is.na(assign_ds[,'Temp3pm']))),mean_t3,1)

colSums(is.na(assign_ds))

assign_ds <- na.omit(assign_ds)


#=== 1. Descriptive look ===
# Scatter plots
p1 <- ggplot(assign_ds, aes(x=MinTemp, y = RainTomorrow)) +
  geom_point()
p1

p2 <- ggplot(assign_ds, aes(x=MaxTemp, y = RainTomorrow)) +
  geom_point()
p2
p3 <- ggplot(assign_ds, aes(x=Rainfall, y = RainTomorrow)) +
  geom_point()
p3
p4 <- ggplot(assign_ds, aes(x=Evaporation, y = RainTomorrow)) +
  geom_point()
p4
p5 <- ggplot(assign_ds, aes(x=Sunshine, y = RainTomorrow)) +
  geom_point()
p5
p6 <- ggplot(assign_ds, aes(x=WindGustDir, y = RainTomorrow)) +
  geom_point()
p6
p7 <- ggplot(assign_ds, aes(x=WindGustSpeed, y = RainTomorrow)) +
  geom_point()
p7
p8 <- ggplot(assign_ds, aes(x=WindDir9am, y = RainTomorrow)) +
  geom_point()
p8
p9 <- ggplot(assign_ds, aes(x=WindDir3pm, y = RainTomorrow)) +
  geom_point()
p9
p10 <- ggplot(assign_ds, aes(x=WindSpeed9am, y = RainTomorrow)) +
  geom_point()
p10

#Scatter plots
p11 <- ggplot(assign_ds, aes(x=WindSpeed3pm, y = RainTomorrow)) +
  geom_point()
p11

p12 <- ggplot(assign_ds, aes(x=Humidity9am, y = RainTomorrow)) +
  geom_point()
p12

p13 <- ggplot(assign_ds, aes(x=Humidity3pm, y = RainTomorrow)) +
  geom_point()
p13

p14 <- ggplot(assign_ds, aes(x=Pressure9am, y = RainTomorrow)) +
  geom_point()
p14

p15 <- ggplot(assign_ds, aes(x=Pressure3pm, y = RainTomorrow)) +
  geom_point()
p15

p16 <- ggplot(assign_ds, aes(x=Cloud9am, y = RainTomorrow)) +
  geom_point()
p16

p17 <- ggplot(assign_ds, aes(x=Cloud3pm, y = RainTomorrow)) +
  geom_point()
p17

p18 <- ggplot(assign_ds, aes(x=Temp9am, y = RainTomorrow)) +
  geom_point()
p18

p19 <- ggplot(assign_ds, aes(x=Temp3pm, y = RainTomorrow)) +
  geom_point()
p19

p20 <- ggplot(assign_ds, aes(x=RainToday, y = RainTomorrow)) +
  geom_point()
p20

figure <- ggarrange(p2, p3, p4,p6,p7,p9,p11,p13,p15,p17,p19,p20, nrow = 6, ncol = 2)
figure
# assign_bak <- assign_ds
assign_ds <- assign_bak


table1 <- table(assign_ds$WindGustDir,assign_ds$RainTomorrow)
#Check for which factor has higher influence
prop.table(table1,margin = 2)
assign_ds$WindGustDir <- as.numeric(as.factor(assign_ds$WindGustDir))
#converting as indicating variable
assign_ds$WindGustDir[assign_ds$WindGustDir != 4] <- 0
assign_ds$WindGustDir[assign_ds$WindGustDir == 4] <- 1


table2 <- table(assign_ds$WindDir9am,assign_ds$RainTomorrow)
#Check for which factor has higher influence
prop.table(table2,margin = 2)
assign_ds$WindDir9am <- as.numeric(as.factor(assign_ds$WindDir9am))
#converting as indicating variable
assign_ds$WindDir9am[assign_ds$WindDir9am != 4] <- 0
assign_ds$WindDir9am[assign_ds$WindDir9am == 4] <- 1

table3 <- table(assign_ds$WindDir3pm,assign_ds$RainTomorrow)
#Check for which factor has higher influence
prop.table(table3,margin = 2)

#All the factors contribute equally and hence assumed insignificant

str(assign_ds)

assign_ds$RainTomorrow <- as.numeric(as.factor(assign_ds$RainTomorrow)) - 1 #To obtain 0/1 instead of 1/2;
assign_ds$RainToday <- as.numeric(as.factor(assign_ds$RainToday)) - 1 # To obtain 0/1 instead of 1/2; 
assign_bak1 <- assign_ds

assign_ds <- assign_ds[,c(1,3:10,12:23)] #removing insignificant variables


# THE DATA.
y = assign_ds[,"RainTomorrow"]
x = as.matrix(assign_ds[,c(2:20)])
str(assign_ds)
x

show(round(cor(x[,1:19]),3))

sample_size = floor(0.75*nrow(assign_ds))
training_val <- sample(1:nrow(assign_ds),sample_size,replace = FALSE) # Find the indexes of observations going into the training set
test_val <- setdiff(1:nrow(assign_ds),training_val)               # Find the indexes of observations going into the test set
training_Data <- assign_ds[training_val,]
test_Data <- assign_ds[test_val,]
test_Data1 <- test_Data[,c(-1,-21)]

xPred = as.matrix(test_Data1)

Nx = ncol(x)

# Specify the data in a list, for later shipment to JAGS:
dataList <- list(
  x = x ,
  y = y ,
  xPred = xPred ,
  Ntotal = length(y),
  Nx = Nx, 
  Npred = nrow(xPred)
)

# First run without initials!
initsList <- list(
  beta0 = 0,
  beta1 = 0,
  beta2 = 0,
  beta3 = 0,
  beta4 = 0
)

modelString = "
data {
  for ( j in 1:Nx ) {
    xm[j]  <- mean(x[,j])
    xsd[j] <-   sd(x[,j])
    for ( i in 1:Ntotal ) {
      zx[i,j] <- ( x[i,j] - xm[j] ) / xsd[j]
    }
  }
}

model {
  for ( i in 1:Ntotal ) {
    # In JAGS, ilogit is logistic:
    y[i] ~ dbern( mu[i] )
      mu[i] <- ( guess*(1/2) + (1.0-guess)*ilogit(beta0+sum(beta[1:Nx]*x[i,1:Nx])) )
  }
  # Priors vague on standardized scale:
  zbeta0 ~ dnorm( 0 , 1/2^2 )
  # non-informative run
  for ( j in 1:Nx ) {
    zbeta[j] ~ dnorm( 0 , 1/2^2 )
  }
  guess ~ dbeta(1,9)
  # Transform to original scale:
  beta[1:Nx] <- zbeta[1:Nx] / xsd[1:Nx]
  beta0 <- zbeta0 - sum( zbeta[1:Nx] * xm[1:Nx] / xsd[1:Nx] )

  # Compute predictions at every step of the MCMC
  for ( k in 1:Npred){
    pred[k] <- ilogit(beta0 + sum(beta[1:Nx] * xPred[k,1:Nx]))
  }
}
" # close quote for modelString
# Write out modelString to a text file
writeLines( modelString , con="TEMPmodel.txt" )

# parameters = c( "zbeta0" , "beta0")
parameters = c( "beta0")
for ( i in 1:Nx){
  # parameters = c(parameters, paste0("zbeta[",i,"]"), paste0("beta[",i,"]")) 
  parameters = c(parameters, paste0("beta[",i,"]")) 
}
for ( i in 1:nrow(xPred)){
  parameters = c(parameters, paste0("pred[",i,"]")) 
}

parameters = c(parameters, "guess")
adaptSteps = 2000  # Number of steps to "tune" the samplers
burnInSteps = 7000
nChains = 4
thinSteps = 13
numSavedSteps = 6000
nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )

runJagsOut <- run.jags( method="parallel" ,
                        model="TEMPmodel.txt" ,
                        monitor=parameters  ,
                        data=dataList ,
                        n.chains=nChains ,
                        adapt=adaptSteps ,
                        burnin=burnInSteps ,
                        sample=numSavedSteps ,
                        thin=thinSteps , summarise=FALSE , plots=FALSE )
codaSamples = as.mcmc.list( runJagsOut )

diagMCMC( codaSamples , parName="beta0" )
for ( i in 1:Nx){
  diagMCMC( codaSamples , parName=paste0("beta[",i,"]") )
}
# diagMCMC( codaSamples , parName="zbeta0" )
# for ( i in 1:Nx){
#   diagMCMC( codaSamples , parName=paste0("zbeta[",i,"]") )
# }
# for ( i in 1:nrow(xPred)){
#   diagMCMC( codaSamples , parName=paste0("pred[",i,"]") )
# }

compVal <- data.frame("beta0" = 15, "beta[1]" = 0, "beta[2]" = 0, "beta[3]" = 0, "beta[4]" =  0,  "beta[5]" =  0, 
                      "beta[6]" =  0, "beta[7]" = 0,"beta[8]" = 0, "beta[9]" = 0,"beta[10]" = 0,
                      "beta[11]" = 0,"beta[12]" = 0, "beta[13]" = 0, "beta[14]" = 0, "beta[15]" = 0,
                      "beta[16]" = 0, "beta[17]" = 0, "beta[18]" = 0, "beta[19]" = 0, check.names=FALSE)

summaryInfo <- smryMCMC_HD( codaSamples = codaSamples , compVal = compVal )
print(summaryInfo)

plotMCMC_HD( codaSamples = codaSamples , data = assign_ds, xName=c("MinTemp","MaxTemp","Rainfall","Evaporation","Sunshine","WindGustDir","WindGustSpeed","WindDir9am","WindSpeed9am","WindSpeed3pm","Humidity9am","Humidity3pm","Pressure9am","Pressure3pm","Cloud9am","Cloud3pm","Temp9am","Temp3pm","RainToday") ,
             yName="RainTomorrow", compVal = compVal, preds = FALSE)

# Predictions for full records in training set
preds <- data.frame(date = test_Data[,1], PredProb = summaryInfo[22:596,3], actual = test_Data[,21] )

threshold <- 0.5# summaryInfo[427,3]
preds[which(preds[,2]<threshold),3] <- 0
preds[which(preds[,2]>threshold),3] <- 1

predsSorted <- preds[order(preds$date),]

table(preds$actual)

actualrain <- test_Data[which(test_Data$Date %in% predsSorted$date),21]


# ============ Predictive check ============

confusionMatrix <- function(resp, pred){
  classRes <- data.frame(response = resp , predicted = pred)
  conf = xtabs(~ predicted + response, data = classRes)
  
  accuracy = sum(diag(conf))/sum(conf)
  accuracy
  precision = conf[1,1]/(conf[1,1]+conf[1,2])
  precision
  recall = conf[1,1]/(conf[1,1]+conf[2,1])
  recall
  Fscore = 2*((precision*recall)/(precision+recall))
  Fscore
  return(list(Accuracy = accuracy, Precision = precision, Recall = recall, Fscore = Fscore, Confusion_Matrix = conf))
}

confusionMatrix(resp = actualrain, pred = predsSorted[,3])


# Saving Data
save.image(file = "firstrun.RData")
save.image(file = "secondrun.RData")
save.image(file = "thirdrun.RData")
save.image(file = "fourthrun.RData")


