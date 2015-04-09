require(plyr)
require(ggplot2)
require(moments)
require(Hmisc)
require(ape)
require(ade4)
source("scripts/functions.R")


percentcover <- read.csv("data/percentcover.csv", header=TRUE)
data <- read.csv("data/data.csv", header=TRUE)
hydro <- read.csv("data/hydro_fixed.csv", header=TRUE)
spatial <- read.csv("data/spatial.csv", header=TRUE)


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

# find data density

data.fieldP <- subset(data, fieldPRESABS == "P")
data.allFieldP <- subset(data, allFieldPRESABS == "P")
data.allP <- subset(data, allPRESABS == "P")

plotcover.fieldP <- ddply(data.fieldP, .(plotID), summarise, plotCover = sum(speciescover))
plotcover.allFieldP <- ddply(data.allFieldP, .(plotID), summarise, plotCover = sum(speciescover))
plotcover.allP <- ddply(data.allP, .(plotID), summarise, plotCover = sum(speciescover))

data.fieldP  <- merge(plotcover.fieldP, data.fieldP)
data.allFieldP  <- merge(plotcover.allFieldP, data.allFieldP)
data.allP  <- merge(plotcover.allP, data.allP)

data.fieldP$relcover <- data.fieldP$speciescover / data.fieldP$plotCover
data.allFieldP$relcover <- data.allFieldP$speciescover / data.allFieldP$plotCover
data.allP$relcover <- data.allP$speciescover / data.allP$plotCover


dataDensity <- data.frame(cbind(plotcover.fieldP$plotID,
                                ddply(data.fieldP, .(plotID), summarise, length = length(allWD))$length,  
                                plotcover.fieldP$plotCover, 
                                plotcover.fieldP$plotCover / plotcover$plotCover,
                                ddply(data.allFieldP, .(plotID), summarise, length = length(allWD))$length,  
                                plotcover.allFieldP$plotCover, 
                                plotcover.allFieldP$plotCover / plotcover$plotCover,
                                ddply(data.allP, .(plotID), summarise, length = length(allWD))$length,                                  
                                plotcover.allP$plotCover, 
                                plotcover.allP$plotCover / plotcover$plotCover,
                                plotcover$plotCover))


colnames(dataDensity) <- c("plotID",
                           "fieldSampledNbsp",
                           "fieldSampledCover", 
                           "fieldSampledProportion", 
                           "allFieldSampledNbsp",
                           "allFieldSampledCover", 
                           "allFieldSampledProportion", 
                           "allValuesNbsp",
                           "allValuesCover", 
                           "allValuesProportion", 
                           "totalCover")

# calculate CWM for different datasets

CWM.field <- ddply(data.fieldP, .(plotID), function(X) data.frame(CWM=weighted.mean(X$fieldWD, X$speciescover),
                                                                  UWM=mean(X$fieldWD)))
CWM.allField <- ddply(data.allFieldP, .(plotID), function(X) data.frame(CWM=weighted.mean(X$allFieldWD, X$speciescover),
                                                                        UWM=mean(X$allFieldWD)))
CWM.all <- ddply(data.allP, .(plotID), function(X) data.frame(CWM=weighted.mean(X$allWD, X$speciescover),
                                                              UWM=mean(X$allWD)))

dataDensity <- cbind(dataDensity, CWM.field$CWM, CWM.allField$CWM, CWM.all$CWM)

cor(cbind(CWM.field$CWM, CWM.allField$CWM, CWM.all$CWM))
plot(data.frame(cbind(CWM.field$CWM, CWM.allField$CWM, CWM.all$CWM)))


# calculate CWM, CWV 

CWM <- ddply(data.allP, .(plotID), function(X) data.frame(CWM=weighted.mean(X$allWD, X$speciescover),
                                                          UWM=mean(X$allWD)))
CWV <- data.frame()

for(i in 1:15) {
  
  x <- data.frame(getCWV(data.allP, i))
  CWV <- rbind(CWV,x) 
}

rm(x)

colnames(CWV) <- c("CWV")

hydroCWM <- cbind(hydro, CWM$CWM)
colnames(hydroCWM)[26] <- "CWM"

hydroCWV <- cbind(hydro, CWV$CWV)
colnames(hydroCWV)[26] <- "CWV"


# get stats!

#getStats(hydroCWM, hydroCWM$CWM, CWM)
#getStats(hydroCWV, hydroCWV$CWV, CWV)

## get real p.adjusted values (using the correct vector with linear and quad vals )##
# this step isn't automated (vals contained in csv file) but is important. 
# Values given by BH p adjustment depend on other values in the vector. 
# So the vector needs to be composed of p values for the model you want. I have 4 quadratic models here.

pvec <- read.csv("output/CWM_sig2.csv", header=TRUE)
pvec$padj <- p.adjust(pvec$p, method="BH")
write.csv(pvec, "output/pvec.csv")

# plot graphs!

#plot.quad(hydroCWM, hydroCWM$CWM, CWM)
#plot.linear(hydroCWM, hydroCWM$CWM, CWM)


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
                      hydroCWM["MA.7daysMinMean"],
                      hydroCWM["HSPeak"],
                      hydroCWM["C_MinM"],
                      hydroCWM["MDFAnnLSNum"])

hydro.signif.pca <- prcomp(hydro.signif, scale.=TRUE, centre=TRUE)

summary(hydro.signif.pca)
hydro.signif.pca$rotatio[,1:3] # loadings
PC1 <- hydro.signif.pca$x[,1] # PC1 site scores


## SPATIAL AUTOCORRELATION TEST ##

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
CWV.na <- na.omit(CWV)

spatial.na.dist <- dist(spatial.na[,2:3])

spatial.na.dist.inv <- as.matrix(1/spatial.na.dist)
diag(spatial.na.dist.inv) <- 0

Moran.I(CWV.na$CWV, spatial.na.dist.inv)


## single species regressions

WDdata <- read.csv("data/WDdata.csv", header=TRUE)
WDraw_hydro <- merge(WDdata, hydro)

WDraw_hydro_cas <- WDraw_hydro[WDraw_hydro$species =="Casuarina cunninghamiana",]
WDraw_hydro_cas$species <- NULL
WDraw_hydro_cas$catname <- NULL
WDraw_hydro_cas$category <- NULL

############

## PCA over all hydro metrics ##

hydro.pca <- prcomp(hydro[2:25], scale.=TRUE, centre=TRUE)
summary(hydro.pca)
plot(hydro.pca)
biplot(hydro.pca)

hydro.pc1 <- hydro.pca$x[,1] # PC1 site scores
cor(hydro.pc1, PC1)

