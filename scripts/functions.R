fitmodels <- function(x) {
    
 if(any(hydroCWM$x == 0)) {
   
   fit.linear  <- lm(CWM ~ x, data = hydroCWM)
   fit.quad    <- lm(CWM ~ x + I(x ^2), data = hydroCWM) 
   
   linear.pval <- anova(fit.linear)[1,"Pr(>F)"]
   quad.pval   <- anova(fit.quad)[1,"Pr(>F)"]
   exp.pval <- c(0)
   
   linear.r2   <- summary(fit.linear)$r.squared
   quad.r2     <- summary(fit.quad)$r.squared
   exp.r2      <- c(0)  
 }
 
 else {
   
  fit.linear  <- lm(CWM ~ x, data = hydroCWM)
  fit.quad    <- lm(CWM ~ x + I(x ^2), data = hydroCWM)  
  fit.exp     <- lm(CWM ~ x, data = hydroCWM)
  
  # p-values
  linear.pval <- anova(fit.linear)[1,"Pr(>F)"]
  quad.pval   <- anova(fit.quad)[1,"Pr(>F)"]
  exp.pval    <- anova(fit.exp)[1,"Pr(>F)"]
  
  # r-squared
  linear.r2   <- summary(fit.linear)$r.squared
  quad.r2     <- summary(fit.quad)$r.squared
  exp.r2      <- summary(fit.exp)$r.squared  
  
  }
  
  fitmodels.df <- cbind(linear.pval, linear.r2, quad.pval, quad.r2, exp.pval, exp.r2)
  
  write.table(fitmodels.df, file="output/fitmodels.csv", append=TRUE, sep=",", col.names=FALSE, row.names=FALSE)
 
}

#columns <- c("linear.pval", "linear.r2", "quad.pval", "quad.r2", "exp.pval", "exp.r2")
#rownames(fitmodels.df) <- colnames(hydroCWM)

