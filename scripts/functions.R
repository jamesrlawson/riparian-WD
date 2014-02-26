fitmodels <- function(x) {
  
  fit.linear  <- lm(CWM ~ x, data = hydroCWM)
  fit.quad    <- lm(CWM ~ x + I(x ^2), data = hydroCWM)
  
  # can't do lm's if there are any NaN values, so we'll omit them
  
  log10x <- log10(x)
  hydroCWM_naomit <- hydroCWM
  hydroCWM_naomit$log10x <- log10x
  hydroCWM_naomit <- na.omit(hydroCWM_naomit)
  
  # but let's get a readout of which hydro metrics were affected
  
 # hydroname <- names(hydroCWM[x])
  
 # if(any(hydroCWM == 0)) {
 #   NaNs <- rbind(hydroname)
 # }
          
  fit.exp     <- lm(CWM ~ log10x, data = hydroCWM_naomit)
  
  # p-values
  linear.pval <- anova(fit.linear)[1,"Pr(>F)"]
  quad.pval   <- anova(fit.quad)[1,"Pr(>F)"]
  exp.pval    <- anova(fit.exp)[1,"Pr(>F)"]
  
  # r-squared
  linear.r2   <- summary(fit.linear)$r.squared
  quad.r2     <- summary(fit.quad)$r.squared
  exp.r2      <- summary(fit.exp)$r.squared  
  
  fitmodels.df <- cbind(linear.pval, linear.r2, quad.pval, quad.r2, exp.pval, exp.r2)
  colnames(fitmodels.df) <- c("linear.pval", "linear.r2", "quad.pval", "quad.r2", "exp.pval", "exp.r2")
  
 # write.csv(NaNs, file="output/NaNs.csv")
  write.csv(fitmodels.df, file="output/NaNs.csv")
  
}


