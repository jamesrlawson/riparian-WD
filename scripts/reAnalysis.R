require(plyr)
require(ggplot2)
require(moments)
require(Hmisc)
require(ape)
require(ade4)
source("scripts/functions.R")


percentcover <- read.csv("data/percentcover.csv", header=TRUE)
data <- read.csv("data/data1.csv", header=TRUE)

# remove observations that average less than 1% cover
percentcover <- subset(percentcover, avgcover > 1)

# remove groundcover stratum and non-woody / freestanding species remaining in other strata
percentcover <- subset(percentcover, stratum != "groundcover")
percentcover <- subset(percentcover, species != "Trophis scandens") 
percentcover <- subset(percentcover, species != "Ripogonum album") 
percentcover <- subset(percentcover, species != "Rubus fruticosus agg.") 
percentcover <- subset(percentcover, species != "Sida rhombifolia") 
percentcover <- subset(percentcover, species != "Solanum americanum") 
percentcover <- subset(percentcover, species != "Senecio prenanthoides")
percentcover <- subset(percentcover, species != "Pteridium esculentum") 


# generate df with total % cover (across all strata) of species at each plot
plotsums <- ddply(percentcover, .(plotID, species), summarise, speciescover = sum(avgcover))
plotcover <- ddply(percentcover, .(plotID), summarise, plotCover = sum(avgcover))

#data <- read.csv("data/data.csv", header=TRUE)

data.fieldP <- subset(data, fieldPRESABS == "P")
data.allP <- subset(data, allPRESABS == "P")

plotcover.fieldP <- ddply(data.fieldP, .(plotID), summarise, plotCover = sum(speciescover))
plotcover.allP <- ddply(data.allP, .(plotID), summarise, plotCover = sum(speciescover))

data.fieldP  <- merge(plotcover.fieldP, data.fieldP)
data.allP  <- merge(plotcover.allP, data.allP)

data.fieldP$relcover <- data.fieldP$speciescover / data.fieldP$plotCover
data.allP$relcover <- data.allP$speciescover / data.allP$plotCover


dataDensity <- data.frame(cbind(plotcover.fieldP$plotID,
                                ddply(data.fieldP, .(plotID), summarise, length = length(allWD))$length,  
                                plotcover.fieldP$plotCover, 
                                plotcover.fieldP$plotCover / plotcover$plotCover,
                                ddply(data.allP, .(plotID), summarise, length = length(allWD))$length,                                  
                                plotcover.allP$plotCover, 
                                plotcover.allP$plotCover / plotcover$plotCover,
                                plotcover$plotCover))
                                
                                
colnames(dataDensity) <- c("plotID",
                           "fieldSampledNbsp",
                           "fieldSampledCover", 
                           "fieldSampledProportion", 
                           "allValuesNbsp",
                           "allValuesCover", 
                           "allValuesProportion", 
                           "totalCover")



CWM <- ddply(data.allP, .(plotID), function(X) data.frame(CWM=weighted.mean(X$allWD, X$speciescover),
                                                                       UWM=mean(X$allWD)))

med <-ddply(data.allP, .(plotID), summarise, med = median(allWD))

range <- ddply(data.allP, .(plotID), summarise, range = max(allWD) - min(allWD))

WDIQR <- ddply(data.allP, .(plotID), summarise, IQR = IQR(allWD))

CV <- ddply(data.allP, .(plotID), summarise, CV = CV(allWD))

SD <- ddply(data.allP, .(plotID), summarise, SD = sd(allWD))

kurt <- ddply(data.allP, .(plotID), summarise, kurt = kurtosis(allWD))

skew <- ddply(data.allP, .(plotID), summarise, skew = skewness(allWD))

var <- ddply(data.allP, .(plotID), summarise, Var = var(allWD))

hydro <- read.csv("data/hydro_fixed.csv", header=TRUE)
hydroCWM <- cbind(hydro, CWM)
getStats(hydroCWM, hydroCWM$CWM, CWM)
getStats(hydroCWM, hydroCWM$UWM, UWM)


plot.linear2(hydroCWM, hydroCWM$CWM, CWM)
plot.quad2(hydroCWM, hydroCWM$CWM, CWM)

hydroMedian <- cbind(hydro, med)
getStats(hydroMedian, hydroMedian$med, median)

hydroRange <- cbind(hydro, range)
getStats(hydroRange, hydroRange$range, range)

FDs <- WDIQR$IQR / range$range #Schleuter et al's FDS
hydroFDs <- cbind(FDs, hydro)
getStats(hydroFDs, hydroFDs$FDs, FDs)

hydroCV <- cbind(hydro, CV)
getStats(hydroCV, hydroCV$CV, CV)
plot.linear2(hydroCV, hydroCV$CV, CV)
plot.quad2(hydroCV, hydroCV$CV, CV)


hydroSD <- cbind(hydro, SD)
getStats(hydroSD, hydroSD$SD, SD)

hydroKurtosis<- cbind(hydro, kurt)
getStats(hydroKurtosis, hydroKurtosis$kurt, kurtosis)

hydroSkew<- cbind(hydro, skew)
getStats(hydroSkew, hydroSkew$skew, skew)
plot.linear2(hydroSkew, hydroSkew$skew, skew)
plot.quad2(hydroSkew, hydroSkew$skew, skew)

hist(subset(data.allP, plotID == 1)$allWD) # smallest skewness
hist(subset(data.allP, plotID == 2)$allWD) # 
hist(subset(data.allP, plotID == 3)$allWD) #
hist(subset(data.allP, plotID == 4)$allWD) # 
hist(subset(data.allP, plotID == 5)$allWD) # 
hist(subset(data.allP, plotID == 6)$allWD) # 
hist(subset(data.allP, plotID == 7)$allWD) # 
hist(subset(data.allP, plotID == 8)$allWD) # 
hist(subset(data.allP, plotID == 9)$allWD) # most negative skewness
hist(subset(data.allP, plotID == 10)$allWD) # 
hist(subset(data.allP, plotID == 11)$allWD) # 
hist(subset(data.allP, plotID == 12)$allWD) # 
hist(subset(data.allP, plotID == 13)$allWD) # 
hist(subset(data.allP, plotID == 14)$allWD) # 
hist(subset(data.allP, plotID == 15)$allWD) # most positive skewness

hydroVar <- cbind(hydro, var)
getStats(hydroVar, hydroVar$Var, Var)


########## PCA  ###########

# create dataframe of hydro metrics for which we find significant relationships with WD CWMs

hydro.signif <- cbind(hydroCWM["AS20YrARI"],
                      hydroCWM["LSPeak"],
                      hydroCWM["M_MDFM"],
                      hydroCWM["BFI"],
                      hydroCWM["C_MDFM"],
                      hydroCWM["CVAnnBFI"],
                      hydroCWM["CVAnnMRateFall"],
                      hydroCWM["CVAnnMRateRise"],
                      hydroCWM["MRateRise"],
                      hydroCWM["HSPeak"],
                      hydroCWM["C_MinM"])


catname <- data.frame("catnamesfull" = as.factor(c(
  "unpredictable intermittent",
  "unpredictable baseflow",
  "unpredictable baseflow",
  "unpredictable intermittent",
  "unpredictable intermittent",
  "unpredictable baseflow",
  "stable winter baseflow",
  "stable winter baseflow",
  "stable winter baseflow",
  "stable winter baseflow",
  "unpredictable baseflow",
  "unpredictable baseflow",
  "unpredictable intermittent",
  "stable winter baseflow",
  "unpredictable intermittent")))

hydro.signif.pca <- prcomp(hydro.signif, scale.=TRUE, centre=TRUE)
summary(hydro.signif.pca)

hydro.signif.pca$rotatio[,1:3] # loadings
hydro.signif.pca$x[,1] # PC1 site scores
PC1 <- hydro.signif.pca$x[,1]
PC2 <- hydro.signif.pca$x[,2]


PC_metrics <- data.frame(PC1,PC2,CWM,skew,CV)
plot(PC_metrics$PC1, PC_metrics$CWM)
plot(PC_metrics$PC2, PC_metrics$CWM)
plot(PC_metrics$PC1, PC_metrics$skew)
plot(PC_metrics$PC2, PC_metrics$skew)
plot(PC_metrics$PC1, range$range)
plot(PC_metrics$PC1, CV$CV)
plot(PC_metrics$PC1, SD$SD)
plot(PC_metrics$PC1, med$med)
plot(PC_metrics$PC1, kurt$kurt)

PC_metrics.na <- na.omit(PC_metrics)

cor(PC_metrics$PC1, PC_metrics$CWM)
cor(PC_metrics$PC2, PC_metrics$CWM)
cor(PC_metrics.na$PC1, PC_metrics.na$skew)
cor(PC_metrics.na$PC1, PC_metrics.na$CV)


PC1_CWM.lm <- lm(CWM$CWM ~ PC1)
summary(PC1_CWM.lm)
#plot(PC1_CWM.lm)

PC1_var.lm <- lm(var$Var ~ PC_metrics$PC1)
summary(PC1_var.lm)
#plot(PC1_var.lm)
plot(var$Var ~ PC_metrics$PC1)

PC1_range.lm <- lm(range$range ~ PC_metrics$PC1)
summary(PC1_range.lm)
plot(PC_metrics$PC1, range$range)

PC1_CV.lm <- lm(CV$CV ~ PC_metrics$PC1)
summary(PC1_CV.lm)
#plot(PC1_var.lm)
plot(CV$CV ~ PC_metrics$PC1)

PC1_skew.lm <- lm(skew$skew ~ PC1, data = PC_metrics)
summary(PC1_skew.lm)
#plot(PC1_skew.lm)
plot(PC_metrics.na$PC1, PC_metrics.na$skew)

PC1_kurtosis.lm <- lm(kurt$kurt ~ PC_metrics$PC1)
summary(PC1_kurtosis.lm)
plot(PC_metrics$PC1, kurt$kurt)


#### ABUNDANCE WEIGHTED METRICS ####


# abundance weighted variance #

getCWV <- function(df, number) {
  x <- subset(df, plotID == number)
  wtd.var(x$allWD, x$speciescover)
}

CWV <- data.frame()

for(i in 1:15) {
  
  x <- data.frame(getCWV(data.allP, i))
  CWV <- rbind(CWV,x) 
}

rm(x)

colnames(CWV) <- c("CWV")
CWV.na <- na.omit(CWV)

plot(PC1, CWV$CWV) 


hydroCWV <- data.frame(hydro, CWV)
getStats(hydroCWV, hydroCWV$CWV, CWV)
plot(CWV$CWV ~ hydro$C_MDFM)

CWV.C_MDFM.lm <- lm(CWV$CWV ~ hydro$C_MDFM)
summary(CWV.C_MDFM.lm)

CWV.lm <- lm(CWV$CWV ~ PC1)
summary(CWV.lm)

# abundance weighted CV #

getCWCV <- function(df, number) {
  x <- subset(df, plotID == number)
  sqrt(wtd.var(x$allWD, x$speciescover)) / wtd.mean(x$allWD, x$speciescover)
}

CWCV <- data.frame()

for(i in 1:15) {
  
  x <- data.frame(getCWCV(data.allP, i))
  CWCV <- rbind(CWCV,x) 
}

rm(x)

colnames(CWCV) <- c("CWCV")
CWCV.na <- na.omit(CWCV)

plot(PC_metrics.na$PC1, CWCV.na$CWCV) 
plot(PC_metrics.na$PC1, PC_metrics.na$CV) 

hydroCWCV <- data.frame(hydro, CWCV)
getStats(hydroCWCV, hydroCWCV$CWCV, CWCV)





############

## PCA over all hydro metrics ##

hydro.pca <- prcomp(hydro[2:25], scale.=TRUE, centre=TRUE)
summary(hydro.pca)
plot(hydro.pca)
biplot(hydro.pca)

hydro.pca1 <- hydro.pca$x[,1] # PC1 site scores
hydro.pca2 <- hydro.pca$x[,2] # PC2 site scores
hydro.pca3 <- hydro.pca$x[,3] # PC3 site scores
hydro.pca4 <- hydro.pca$x[,4] # PC4 site scores

hydro.pca$rotatio[,1:4] # loadings

cor(PC1, hydro.pca1) # PC1 of hydro.signif.pca is the same as PC1 of hydro.pca

cor(CWM, hydro.pca1)
cor(CWM, hydro.pca2)
cor(CWM, hydro.pca3)
cor(CWM, hydro.pca4)



plot(PC_metrics$PC2,range$range)



## SPATIAL AUTOCORRELATION TEST ##

require(ade4)
require(ape)
spatial <- read.csv("data/spatial.csv", header=TRUE)

# for hydro metrics #

spatial.dist <- dist(spatial[,2:3])
hydro.dist <- dist(hydro[,2:25])
hydro.signif.dist <- dist(hydro.signif)

mantel.rtest(spatial.dist, hydro.dist, nrepet=9999)
mantel.rtest(spatial.dist, hydro.signif.dist, nrepet=9999)

# for CWM #

spatial.dist.inv <- as.matrix(1/spatial.dist)
diag(spatial.dist.inv) <- 0

Moran.I(CWM$CWM, spatial.dist.inv)

# remove plotID 2 for CWV analysis #

spatial.na <- spatial
spatial.na <- spatial.na[-2,]
hydro.na <- hydro[-2,]
spatial.na.dist <- dist(spatial.na[,2:3])

spatial.na.dist.inv <- as.matrix(1/spatial.na.dist)
diag(spatial.na.dist.inv) <- 0

Moran.I(CWV.na$CWV, spatial.na.dist.inv)

# for range #

Moran.I(range$range, spatial.dist.inv)

