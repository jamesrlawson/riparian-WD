fitmodels <- function(x) {
     
  fit.linear  <- lm(CWM ~ x, data = hydroCWM)
  fit.quad    <- lm(CWM ~ x + I(x ^2), data = hydroCWM)  
  fit.exp     <- lm(CWM ~ log10(x), subset = x > 0, data = hydroCWM) # THIS IS THE LINE THAT HANDLES INF's!
  
  # p-values
  linear.pval <- anova(fit.linear)[1,"Pr(>F)"]
  quad.pval   <- anova(fit.quad)[1,"Pr(>F)"]
  exp.pval    <- anova(fit.exp)[1,"Pr(>F)"]
  
  # r-squared
  linear.r2   <- summary(fit.linear)$r.squared
  quad.r2     <- summary(fit.quad)$r.squared
  exp.r2      <- summary(fit.exp)$r.squared  
  
  
  fitmodels.df <- cbind(linear.pval, linear.r2, quad.pval, quad.r2, exp.pval, exp.r2)
  
  write.table(fitmodels.df, file="output/fitmodels.csv", append=TRUE, sep=",", col.names=FALSE, row.names=FALSE)
 
}


bestmodel <- function(data) { #group is the dataset 
  
  for (i in 1:length(data$bestmodel))  {
    
    hydro <- data$metric[i]
    
    if (data$bestmodel[i] == 1) {
      plot.linear(hydro)
    }
    else {
      if (data$bestmodel[i] == 2) {
        plot.quad(hydro)
      }
      else {
        if (data$bestmodel[i] ==3) {
          plot.exp(hydro)
        }
      }
    }
  }

}

plot.linear <- function(x) {
  p <- qplot(x, CWM, data = hydroCWM)
  p <- p + geom_point()
  p <- p + stat_smooth(method = "lm", formula = y ~ x, se=TRUE)   
  print(p)  
}

plot.quad <- function(x) {
  p <- qplot(x, CWM, data = hydroCWM)
  p <- p + geom_point()
  p <- p + stat_smooth(method = "lm", formula = y ~ x + I(x^2), se=TRUE)   
  print(p)  
}

plot.exp <- function(x) {
  p <- qplot(x, CWM, data = hydroCWM)
  p <- p + geom_point()
  p <- p + stat_smooth(method = "lm", formula = y ~ log10(x), se=TRUE)   
  print(p)  
}
