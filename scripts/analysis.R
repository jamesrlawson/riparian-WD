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

comparisonGroups <- read.csv("data/comparisonGroups.csv", header=FALSE)

fitmodels.output$comparisonGroups <- comparisonGroups[[2]]


floodmodels <- fitmodels.output[fitmodels.output$comparisonGroups == "1",]
unpredictablemodels <- fitmodels.output[fitmodels.output$comparisonGroups == "2",]




