fitmodels <- function(x) {
  
  fit.linear  <- lm(CWM ~ x, data = CWM)
  fit.quad    <- lm(CWM ~ x + I(x ^2), data = CWM)
  fit.exp     <- lm(CWM ~ log10(x))
  
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
  
}


