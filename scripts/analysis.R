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

# now, nspect output and manually determine best fitting model 
# values 0, 1, 2, 3 
# (corresponding to whether no significant models, linear, quadratic or exponential)

bestmodels.flood <- as.data.frame(rownames(floodmodels))
bestmodels.flood$bestmodel <- c(0,0,0,2,2,2,0,0,2)
colnames(bestmodels.flood)[1] <- c("metric")
bestmodels.flood$metric <- as.character(bestmodels.flood$metric)

bestmodels.unpredictable <- as.data.frame(rownames(unpredictablemodels))
bestmodels.unpredictable$bestmodel <- c(0,0,0,0,0,0,0,2,2,2,1,0,2,2,0)
colnames(bestmodels.unpredictable)[1] <- c("metric")
bestmodels.unpredictable$metric <- as.character(bestmodels.unpredictable$metric) #because the metrics were reading in as factors

#put catname in because we need it now for our plots

hydroCWM$catname <- hydro$catname

plot.quad(hydroCWM)
plot.linear(hydroCWM)
plot.exp(hydroCWM)



