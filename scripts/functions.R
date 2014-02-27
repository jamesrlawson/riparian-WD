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

# using a loop allows us to access colnames, where plyr does not

plot.linear <- function(df) {
  
  for(i in 1:ncol(df)) {
    hydro <- df[[i]]  
    hydroname <- as.expression(names(df[i]))   # could also ask hydroname to refer to a vector of proper label names
    
    fit.linear <- lm(CWM ~ hydro, data = df)
    
    p <- qplot(hydro, CWM, data = df) 
    p = p + geom_point(aes(shape = catname))
    p <- p + scale_shape_discrete(name = "Hydrological \n class", labels = c("stable winter baseflow", "unpredictable baseflow", "unpredictable intermittent"))
    p <- p + stat_smooth(method = "lm", formula = y ~ x, se=TRUE, col="black") 
    p = p + xlab(hydroname)
    p = p + ylab("AWM wood density (g/cm^3)")
    p = p + ggtitle(paste("Adj R2 = ",signif(summary(fit.linear)$adj.r.squared, 5),
                          "; p =",signif(summary(fit.linear)$coef[2,4], 5)))
    
    p = p + theme(panel.grid.major = element_line(size = .5, color = "grey"),
                  #increase size of axis lines
                  axis.line = element_line(size=.7, color = "black"),
                  legend.position = "bottom",
                  panel.background = element_blank(),      
                  plot.title = element_text(size=10),
                  text = element_text(size=12))   
    
    print(p)
    
  }
}



plot.quad <- function(df, pvals) {
  
  for(i in 2:ncol(df)) {
    hydro <- df[[i]]  
    hydroname <- as.expression(names(df[i]))   # could also ask hydroname to refer to a vector of proper label names
    
    fit.quad <- lm(CWM ~ hydro + I(hydro^2), data = df)
    
    p <- qplot(hydro, CWM, data = df) 
    p = p + geom_point(aes(shape = catname))
    p <- p + scale_shape_discrete(name = "Hydrological \n class", labels = c("stable winter baseflow", "unpredictable baseflow", "unpredictable intermittent"))
    p <- p + stat_smooth(method = "lm", formula = y ~ x + I(x^2), se=TRUE, col="black") 
    p = p + xlab(hydroname)
    p = p + ylab("AWM wood density (g/cm^3)")
    p = p + ggtitle(paste("Adj. R2 = ",signif(summary(fit.quad)$adj.r.squared, 5),
 #                         "; p =",signif(summary(fit.quad)$coef[2,4], 5)
                          "; p =",pvals$quad.padj[i]  ))
    
    p = p + theme(panel.grid.major = element_line(size = .5, color = "grey"),
                  axis.line = element_line(size=.7, color = "black"),
                  legend.position = "bottom",
                  panel.background = element_blank(),      
                  plot.title = element_text(size=10),
                  text = element_text(size=12))   
    
    print(p)
    
  }
}

plot.exp <- function(df) {
  
  for(i in 1:ncol(df)) {
    hydro <- df[[i]]  
    hydroname <- as.expression(names(df[i]))   # could also ask hydroname to refer to a vector of proper label names
    
    fit.exp <- lm(CWM ~ hydro + log10(hydro), data = df)
    
    p <- qplot(hydro, CWM, data = df) 
    p = p + geom_point(aes(shape = catname))
    p <- p + scale_shape_discrete(name = "Hydrological \n class", labels = c("stable winter baseflow", "unpredictable baseflow", "unpredictable intermittent"))
    p <- p + stat_smooth(method = "lm", formula = y ~ log10(x), se=TRUE, col="black") 
    p = p + xlab(hydroname)
    p = p + ylab("AWM wood density (g/cm^3)")
    p = p + ggtitle(paste("Adj. R2 = ",signif(summary(fit.exp)$adj.r.squared, 5),
                          #                         "; p =",signif(summary(fit.quad)$coef[2,4], 5)
                          "; p =",padj$quad.padj[i]  ))
    
    p = p + theme(panel.grid.major = element_line(size = .5, color = "grey"),
                  axis.line = element_line(size=.7, color = "black"),
                  legend.position = "bottom",
                  panel.background = element_blank(),      
                  plot.title = element_text(size=10),
                  text = element_text(size=12))   
    
    print(p)
    
  }
}
