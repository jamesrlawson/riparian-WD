library(plyr)
library(ggplot2)


source("scripts/functions.R")

CWM <- read.csv("data/CWM.csv", header=TRUE)
hydro <- read.csv("data/hydro.csv", header=TRUE)

hydroCWM <- merge(hydro, CWM)
hydroCWM$catname <- NULL # we won't use this just yet

# get pvals and rsquared vals for linear vs quadratic vs exponential models 

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
fitmodels.output <- fitmodels.output[-27,] # CWM
fitmodels.output <- fitmodels.output[-26,] # category
fitmodels.output <- fitmodels.output[-1,] # plotID


fitmodels.output <- format(fitmodels.output, digits=2, trim=TRUE) 

# read in comparison groups
comparisonGroups <- read.csv("data/comparisonGroups.csv", header=FALSE)

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
bestmodels.flood$bestmodel <- c(0,0,0,2,2,2,0,0,3)
colnames(bestmodels.flood)[1] <- c("metric")
bestmodels.flood$metric <- as.character(bestmodels.flood$metric)

bestmodels.unpredictable <- as.data.frame(rownames(unpredictablemodels))
bestmodels.unpredictable$bestmodel <- c(0,0,0,0,0,0,0,2,2,2,1,0,2,1,0)
colnames(bestmodels.unpredictable)[1] <- c("metric")
bestmodels.unpredictable$metric <- as.character(bestmodels.unpredictable$metric) #because the metrics were reading in as factors

#put catname in because we need it now for our plots

#hydroCWM$catname <- hydro$catname

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
                                        hydroCWM$HSPeaknorm,
                                        hydroCWM$BFI,  
                                        hydroCWM$CVAnnBFI,	
                                        hydroCWM$C_MDFM,
                                        hydroCWM$M_MinM,
                                        hydro$catname))
                                  
colnames(hydro.quad)       <- c("zCWM", # z so that CWM ends up last when sorted alphabetically
                                "CVAnnMRateRise", 
                                "CVAnnMRateFall",
                                "HSPeaknorm",
                                "BFI",  
                                "CVAnnBFI",	
                                "C_MDFM",
                                "M_MinM",
                                "catname")
#hydro.quad$catname <- as.factor(hydro.quad$catname)
hydro.quad$catname <- NULL

 #sort columns alphabetically to match subsetting output later on

hydro.quad <- hydro.quad[,order(names(hydro.quad))]

hydro.exp <- as.data.frame(cbind(hydroCWM$CWM,
                                 hydroCWM$AS20YrARInorm,
                                 hydro$catname))
colnames(hydro.exp) <-       c("zCWM",
                               "AS20YrARInorm",
                               "catname")
#hydro.exp$catname <- as.factor(hydro.exp$catname)
hydro.exp$catname <- NULL
hydro.exp <- hydro.exp[,order(names(hydro.exp))]

                                            
hydro.linear <- as.data.frame(cbind(hydroCWM$CWM,
                                    hydroCWM$LSPeaknorm,
                                    hydroCWM$M_MDFM,
                                    hydro$catname))
                                                
colnames(hydro.linear) <-       c("zCWM",
                                  "LSPeaknorm", 
                                  "M_MDFM",
                                  "catname")
#hydro.linear$catname <- as.factor(hydro.linear$catname)
hydro.linear$catname <- NULL
hydro.linear <- hydro.linear[,order(names(hydro.linear))]

# and all the remainders plotted as quads so I can look at outliers etc.

hydro.nonsignif      <- as.data.frame(cbind(hydroCWM$CWM,
                                            hydroCWM$MDFAnnHSNum,	
                                            hydroCWM$CVAnnHSNum,	
                                            hydroCWM$CVAnnHSPeak,	
                                            hydroCWM$MRateRisenorm,	
                                            hydroCWM$MRateFallnorm,	
                                            hydroCWM$MDFAnnZer,	
                                            hydroCWM$MDFAnnUnder0.1,	
                                            hydroCWM$LSMeanDur,	
                                            hydroCWM$MDFAnnLSNum,	
                                            hydroCWM$CVAnnLSNum,	
                                            hydroCWM$CVAnnLSPeak,	
                                            hydroCWM$CVAnnLSMeanDur,	
                                            hydroCWM$C_MinM,	
                                            hydroCWM$MA.7daysMinMeannorm,
                                            hydro$catname)) 
colnames(hydro.nonsignif    ) <- c("zCWM", 
                                    "MDFAnnHSNum",
                                    "CVAnnHSNum",
                                    "CVAnnHSPeak",
                                    "MRateRisenorm",
                                    "MRateFallnorm",
                                    "MDFAnnZer",
                                    "MDFAnnUnder0.1",
                                    "LSMeanDur",
                                    "MDFAnnLSNum",
                                    "CVAnnLSNum",
                                    "CVAnnLSPeak",
                                    "CVAnnLSMeanDur",
                                    "C_MinM",
                                    "MA.7days.MinMeannorm",
                                    "catname")
#hydro.nonsignif$catname <- as.factor(hydro.nonsignif$catname)
hydro.nonsignif$catname <- NULL
hydro.nonsignif <- hydro.nonsignif[,order(names(hydro.nonsignif))]

## the next thing to do is subset padj according to $bestmodel and use the resulting
## df as the input df for padj in the function

hydro.linear.padj <- padj[padj$bestmodel == 1,]
hydro.quad.padj <- padj[padj$bestmodel == 2,]
hydro.exp.padj <- padj[padj$bestmodel == 3,]
hydro.nonsignif.padj <- padj[padj$bestmodel == 0,]

# plot graphs!

plot.quad(hydro.quad, hydro.quad.padj)
plot.exp(hydro.exp, hydro.exp.padj)
plot.linear(hydro.linear, hydro.linear.padj)
plot.quad(hydro.nonsignif, hydro.nonsignif.padj)

### now lets get some graphs for means with error bars for standard error

# for raw data

WDraw <- read.csv("data/WDraw.csv", header=TRUE)
cats <- as.data.frame(cbind(hydro$plotID, hydro$category))
names(cats) <- c("plotID","cats")
WDraw_cats <- merge(WDraw, cats)

raw_cats.mean <- tapply(WDraw_cats$heart.avg, WDraw_cats$cats, mean)
raw_cats.stderr <- tapply(WDraw_cats$heart.avg, WDraw_cats$cats, stderr)

raw_cats.df <- as.data.frame(cbind("category"=c(1,2,3), mean = raw_cats.mean, stderr = raw_cats.stderr))
raw_cats.df$category <- as.factor(raw_cats.df$category)

raw_cats.df$labels <- c("Difference in mean wood density between classes (raw values)", "hydrological class", "mean wood density (g/cm^3)")

plot.means(raw_cats.df)

# for CWMs

WDCWM <- as.data.frame(cbind("plotID" = hydroCWM$plotID, "cats" = hydroCWM$category, "CWM" = hydroCWM$CWM))

CWM_cats.mean <- tapply(WDCWM$CWM, WDCWM$cats, mean)
CWM_cats.stderr <- tapply(WDCWM$CWM, WDCWM$cats, stderr)

CWM_cats.df <- as.data.frame(cbind("category"=c(1,2,3), mean = CWM_cats.mean, stderr = CWM_cats.stderr))
CWM_cats.df$category <- as.factor(CWM_cats.df$category)
raw_cats.df$labels <- c("Difference in mean wood density between classes (abundance weighted)", "hydrological class", "mean wood density (g/cm^3)")

plot.means(CWM_cats.df)
