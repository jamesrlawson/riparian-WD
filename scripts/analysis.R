require(plyr)
require(ggplot2)
require(moments)
require(Hmisc)
require(ape)
require(ade4)
source("scripts/functions.R")


percentcover <- read.csv("data/percentcover.csv", header=TRUE)
data <- read.csv("data/data.csv", header=TRUE)

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

cor(cbind(CWV.field$CWV, CWV.allField$CWV, CWV.all$CWV))


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


blah <- lm(CWV$CWV ~ PC1 + I(PC1^2))
summary(blah)

hydro <- read.csv("data/hydro_fixed.csv", header=TRUE)

hydroCWM <- cbind(hydro, CWM$CWM)
colnames(hydroCWM)[26] <- "CWM"


getStats(hydroCWM, hydroCWM$CWM, CWM)

#hydroCWV <- cbind(hydro, CWV)

# import hydrological data


## remove snowy creek as outlier (site 10)
#hydroCWM <- hydroCWM[-10,]

#hydroCWM$catname <- NULL # we won't use this just yet

# get pvals and rsquared vals for linear vs quadratic vs exponential models (outputs to file)

ldply(hydroCWM, fitmodels)

# remember to delete fitmodels.csv if you're going to run this multiple times  
# otherwise the function will keep appending
fitmodels.output <- read.csv("output/fitmodels.csv", header=FALSE)

hydronames <- colnames(hydroCWM)
modelcols <- c("linear.pval", "linear.r2", "quad.pval", "quad.r2", "exp.pval", "exp.r2")

rownames(fitmodels.output) <- hydronames
colnames(fitmodels.output) <- modelcols

rm(hydronames)
rm(modelcols)

#not using these rows
fitmodels.output <- fitmodels.output[-26,] # CWM
fitmodels.output <- fitmodels.output[-1,] # plotID

fitmodels.output <- format(fitmodels.output, digits=3, trim=TRUE) 

fitmodels.output$names <- rownames(fitmodels.output)
fitmodels.output <- fitmodels.output[order(fitmodels.output$names),]
fitmodels.output$names <- NULL

# read in comparison groups
comparisonGroups <- read.csv("data/comparisonGroups.csv", header=TRUE)
comparisonGroups <- comparisonGroups[order(comparisonGroups$names),]

# merge
fitmodels.output$comparisonGroups <- comparisonGroups[[2]]

floodmodels <- fitmodels.output[fitmodels.output$comparisonGroups == "1",]
unpredictablemodels <- fitmodels.output[fitmodels.output$comparisonGroups == "2",]

# adjust pvals for multiple comparison
floodmodels$linear.padj <- p.adjust(floodmodels$linear.pval, method="BH")
floodmodels$quad.padj <- p.adjust(floodmodels$quad.pval, method="BH")
floodmodels$exp.padj <- p.adjust(floodmodels$exp.pval, method="BH")
unpredictablemodels$linear.padj <- p.adjust(unpredictablemodels$linear.pval, method="BH")
unpredictablemodels$quad.padj <- p.adjust(unpredictablemodels$quad.pval, method="BH")
unpredictablemodels$exp.padj <- p.adjust(unpredictablemodels$exp.pval, method="BH")

write.csv(floodmodels, file="output/floodmodels.csv")
write.csv(unpredictablemodels, file="output/unpredictablemodels.csv")

       
# now, inspect output and manually determine best fitting model 
# values 0, 1, 2, 3 
# (corresponding to whether no significant models, linear, quadratic or exponential)

bestmodels.flood <- as.data.frame(rownames(floodmodels))
bestmodels.flood$bestmodel <- c(1,
                                0,
                                0,
                                2,
                                2,
                                2,
                                0,
                                0,
                                2
)
colnames(bestmodels.flood)[1] <- c("metric")
bestmodels.flood$metric <- as.character(bestmodels.flood$metric)

bestmodels.unpredictable <- as.data.frame(rownames(unpredictablemodels))
bestmodels.unpredictable$bestmodel <- c(2,
                                        2,
                                        2,
                                        2,
                                        0,
                                        0,
                                        0,
                                        0,
                                        1,
                                        1,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0
)
colnames(bestmodels.unpredictable)[1] <- c("metric")
bestmodels.unpredictable$metric <- as.character(bestmodels.unpredictable$metric) #because the metrics were reading in as factors

# now combine bestmodels with padj

padj <- rbind(floodmodels, unpredictablemodels)
padj$metric <- rownames(padj)
bestmodels <- rbind(bestmodels.flood, bestmodels.unpredictable)
padj <- merge(padj, bestmodels)

# now make dataframes that contain the right hydro parameters for each model
# hydro.quad
# hydro.linear
# hydro.exp
# hydro.nonsignif, uses plot.quad

hydro.quad <- as.data.frame(cbind(hydroCWM$CWM,
                                        hydroCWM$CVAnnMRateRise,	
                                        hydroCWM$CVAnnMRateFall,
                                        hydroCWM$MRateRise,
                                        hydroCWM$HSPeak,
                                        hydroCWM$BFI,  
                                        hydroCWM$CVAnnBFI,	
                                        hydroCWM$C_MDFM,
                                        hydroCWM$C_MinM
                                  ))
                                  
colnames(hydro.quad)       <- c("zCWM", # z so that CWM ends up last when sorted alphabetically
                                "CVAnnMRateRise", 
                                "CVAnnMRateFall",
                                "MRateRise",
                                "HSPeak",
                                "BFI",  
                                "CVAnnBFI",	
                                "C_MDFM",
                                "C_MinM"
                                )

 #sort columns alphabetically to match subsetting output later on

hydro.quad <- hydro.quad[,order(names(hydro.quad))]

                                            
hydro.linear <- as.data.frame(cbind(hydroCWM$CWM,
                                    hydroCWM$AS20YrARI,
                                    hydroCWM$M_MDFM,
                                    hydroCWM$LSPeak                                    
                                    ))
                                                
colnames(hydro.linear) <-       c("zCWM",
                                  "AS20YrARI", 
                                  "M_MDFM",
                                  "LSPeak"
                                  )

hydro.linear <- hydro.linear[,order(names(hydro.linear))]

# and all the remainders plotted as quads so I can look at outliers etc.

hydro.nonsignif      <- as.data.frame(cbind(hydroCWM$CWM,
                                            hydroCWM$MDFAnnHSNum,	
                                            hydroCWM$CVAnnHSNum,	
                                            hydroCWM$CVAnnHSPeak,	
                                            hydroCWM$MRateFall,	
                                            hydroCWM$MDFAnnZer,	
                                            hydroCWM$MDFAnnUnder0.1,	
                                            hydroCWM$LSMeanDur,	
                                            hydroCWM$MDFAnnLSNum,	
                                            hydroCWM$CVAnnLSNum,	
                                            hydroCWM$CVAnnLSPeak,	
                                            hydroCWM$CVAnnLSMeanDur,	
                                            hydroCWM$M_MinM,
                                            hydroCWM$MA.7daysMinMean                                               
                                            )) 

colnames(hydro.nonsignif    ) <- c("zCWM", 
                                    "MDFAnnHSNum",
                                    "CVAnnHSNum",
                                    "CVAnnHSPeak",
                                    "MRateFall",
                                    "MDFAnnZer",
                                    "MDFAnnUnder0.1",
                                    "LSMeanDur",
                                    "MDFAnnLSNum",
                                    "CVAnnLSNum",
                                    "CVAnnLSPeak",
                                    "CVAnnLSMeanDur",
                                    "M_MinM",
                                    "MA.7daysMinMean"
                                    )

hydro.nonsignif <- hydro.nonsignif[,order(names(hydro.nonsignif))]

## the next thing to do is subset padj according to $bestmodel and use the resulting
## df as the input df for padj in the function

hydro.linear.padj <- padj[padj$bestmodel == 1,]
hydro.quad.padj <- padj[padj$bestmodel == 2,]
hydro.nonsignif.padj <- padj[padj$bestmodel == 0,]

# add figure panel labels (a,b,c etc.) to padj, so they can be used in the graph (sloppy but neat...)

hydro.quad.padj$figlabel <- c("(a)","(c)","(e)","(f)","(g)","(f)","(a)","(d)")
hydro.linear.padj$figlabel <- c("(b)","(g)","(b)")
hydro.nonsignif.padj$figlabel <- c("NA","(c)","(NA)","NA","NA","NA","(d)","(h)","NA","NA","NA","NA", "(e)")

# get stats!

getStats(hydroCWM, hydroCWM$CWM, CWM)
getStats(hydroCWV, hydroCWV$CWV, CWV)
#getStats(hydroRange, hydroRange$range, range)

# get real p.adjusted values (using the correct vector with linear and quad vals)

pvec <- read.csv("output/CWM_sig2.csv", header=TRUE)
pvec$padj <- p.adjust(pvec$p, method="BH")
write.csv(pvec, "output/pvec.csv")

# plot graphs!

plot.quad(hydro.quad, hydro.quad.padj)
#plot.exp(hydro.exp, hydro.exp.padj)
plot.linear(hydro.linear, hydro.linear.padj)
plot.quad(hydro.nonsignif, hydro.nonsignif.padj)


hydroCWM$zCWM <- hydroCWM$CWM

plot.linear(hydroCWM, hydro.quad.padj)

### now lets get some graphs for means with error bars for standard error

# for raw data

#WDraw <- read.csv("data/WDdata.csv", header=TRUE)
#cats <- read.csv("data/categories.csv", header=TRUE)
#WDraw_cats <- merge(WDraw, cats)
#WDraw_cats$cats <- as.factor(WDraw_cats$cats)

#raw_cats.mean <- tapply(WDraw_cats$WD, WDraw_cats$cats, mean)
#raw_cats.stderr <- tapply(WDraw_cats$WD, WDraw_cats$cats, stderr)

#raw_cats.df <- as.data.frame(cbind("category"=c(1,2,3), mean = raw_cats.mean, stderr = raw_cats.stderr))
#raw_cats.df$category <- as.factor(raw_cats.df$category)
#raw_cats.df$labels <- c("Hydrological class", "(b)", "x")

#raw_cats.aov <- aov(WD ~ cats, data = WDraw_cats)
#summary(raw_cats.aov)
#TukeyHSD(raw_cats.aov)


# for CWMs

cats <- read.csv("data/categories.csv", header=TRUE)

WDCWM <- as.data.frame(cbind("plotID" = hydroCWM$plotID, "cats" = cats$cats, "CWM" = hydroCWM$CWM))
WDCWM$cats <- as.factor(WDCWM$cats)

CWM_cats.mean <- tapply(WDCWM$CWM, WDCWM$cats, mean)
CWM_cats.stderr <- tapply(WDCWM$CWM, WDCWM$cats, stderr)

CWM_cats.df <- as.data.frame(cbind("category"=c(1,2,3), mean = CWM_cats.mean, stderr = CWM_cats.stderr))
CWM_cats.df$category <- as.factor(CWM_cats.df$category)
CWM_cats.df$labels <- c("Hydrological class", "(a)", "x") # x is a placeholder


CWM_cats.aov <- aov(CWM ~ cats, data = WDCWM)
summary(CWM_cats.aov)
TukeyHSD(CWM_cats.aov)

plot.means2(CWM_cats.df)


# for CWV 

WDCWM$CWV <- CWV$CWV
#CWV[2,] <- c(0) # replacing NaN with zero for site 2
#WDCWM.na <- na.omit(WDCWM)
CWV_cats.aov <- aov(CWV ~ cats, data = WDCWM)
summary(CWV_cats.aov)
TukeyHSD(CWV_cats.aov)

CWV_cats.mean <- tapply(WDCWM.na$CWV, WDCWM.na$cats, mean)
CWV_cats.stderr <- tapply(WDCWM.na$CWV, WDCWM.na$cats, stderr)

CWV_cats.df <- as.data.frame(cbind("category"=c(1,2,3), mean = CWV_cats.mean, stderr = CWV_cats.stderr))
CWV_cats.df$category <- as.factor(CWV_cats.df$category)
CWV_cats.df$labels <- c("Hydrological class", "(a)", "x") # x is a placeholder

plot.means2(CWV_cats.df)


# for range

WDCWM$range <- range$range
range_cats.aov <- aov(range ~ cats, data = WDCWM)
summary(range_cats.aov)
TukeyHSD(range_cats.aov)

range_cats.mean <- tapply(WDCWM$range, WDCWM$cats, mean)
range_cats.stderr <- tapply(WDCWM$range, WDCWM$cats, stderr)

range_cats.df <- as.data.frame(cbind("category"=c(1,2,3), mean = range_cats.mean, stderr = range_cats.stderr))
range_cats.df$category <- as.factor(range_cats.df$category)
range_cats.df$labels <- c("Hydrological class", "(a)", "x") # x is a placeholder

plot.means2(range_cats.df)



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
PC1 <- hydro.signif.pca$x[,1] # PC1 site scores


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
CWV.na <- na.omit(CWV)

spatial.na.dist <- dist(spatial.na[,2:3])

spatial.na.dist.inv <- as.matrix(1/spatial.na.dist)
diag(spatial.na.dist.inv) <- 0

Moran.I(CWV.na$CWV, spatial.na.dist.inv)






## Permutational Multivariate Analysis of Variance Using Distance Matrices ##
# in vegan #

library(vegan)

hydrocats <- hydro
hydrocats$plotID <- NULL

cats <- read.csv("data/categories.csv", header=TRUE)

hydrocats$cats <- cats$cats

hydrocats12 <- subset(hydrocats, cats < 3)
hydrocats13 <- subset(hydrocats, cats != 2)
hydrocats23 <- subset(hydrocats, cats > 1)

hydrocats12.cats <- hydrocats12$cats
hydrocats13.cats <- hydrocats13$cats
hydrocats23.cats <- hydrocats23$cats
hydrocats12$cats <- NULL
hydrocats13$cats <- NULL
hydrocats23$cats <- NULL

hydro12.dist <- vegdist(hydrocats12, method = "bray")
hydro13.dist <- vegdist(hydrocats13, method = "bray")
hydro23.dist <- vegdist(hydrocats23, method = "bray")

hydrocats12$cats <- as.factor(hydrocats12.cats)
hydrocats13$cats <- as.factor(hydrocats13.cats)
hydrocats23$cats <- as.factor(hydrocats23.cats)

hydro12.adonis <- adonis(hydro12.dist ~ cats, data = hydrocats12)
hydro13.adonis <- adonis(hydro13.dist ~ cats, data = hydrocats13)
hydro23.adonis <- adonis(hydro23.dist ~ cats, data = hydrocats23)

print(hydro12.adonis)
print(hydro13.adonis)
print(hydro23.adonis)







## single species regressions

WDdata <- read.csv("data/WDdata.csv", header=TRUE)
WDraw_hydro <- merge(WDdata, hydro)

WDraw_hydro_cas <- WDraw_hydro[WDraw_hydro$species =="Casuarina cunninghamiana",]
WDraw_hydro_cas$species <- NULL
WDraw_hydro_cas$catname <- NULL
WDraw_hydro_cas$category <- NULL

#WDraw_hydro_tris <- WDraw_hydro[WDraw_hydro$species =="Tristaniopsis laurina",]
#WDraw_hydro_tris$species <- NULL
#WDraw_hydro_tris$catname <- NULL
#WDraw_hydro_tris$category <- NULL

#WDraw_hydro_lep <- WDraw_hydro[WDraw_hydro$species =="Leptospermum brevipes",]
#WDraw_hydro_lep$species <- NULL
#WDraw_hydro_lep$catname <- NULL
#WDraw_hydro_lep$category <- NULL

#WDraw_hydro_aca <- WDraw_hydro[WDraw_hydro$species =="Acacia dealbata",]
#WDraw_hydro_aca$species <- NULL
#WDraw_hydro_aca$catname <- NULL
#WDraw_hydro_aca$category <- NULL

#fitspecies(hydro, WDraw_hydro_tris, trilau)
