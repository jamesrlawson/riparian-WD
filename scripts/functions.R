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


plot.linear <- function(df, pvals) {
  
  for(i in 1:ncol(df)) {
    hydro <- df[[i]]  
    hydroname <- as.expression(names(df[i]))   # could also ask hydroname to refer to a vector of proper label names
    
    fit.linear <- lm(zCWM ~ hydro, data = df)
    catname <- as.factor(c(3,2,2,3,3,2,1,1,1,1,2,2,3,1,3))
    
    p <- qplot(hydro, zCWM, data = df) 
    p = p + geom_point(aes(shape = catname), size =3)
    p <- p + scale_shape_discrete(name = "Hydrological \n class", labels = c("stable winter baseflow", "unpredictable baseflow", "unpredictable intermittent"))
    p <- p + stat_smooth(method = "lm", formula = y ~ x, se=TRUE, col="black") 
    p = p + xlab(hydroname)
    p = p + ylab("AWM wood density (g/cm^3)")
    p = p + annotate("text",                    
                     x=max(hydro)/1.5, y=0.5,
                     label=paste("R^2 = ",signif(summary(fit.linear)$r.squared, 5),
                                 "\np.adj =",pvals$linear.padj[i]),
                     size = 4)
    
    p = p + theme(panel.grid.major = element_line(size = .5, color = "grey"),
                  axis.line = element_line(size=.7, color = "black"),
                  legend.position = "bottom",
                  panel.background = element_blank(),      
                  plot.title = element_text(size=12),
                  axis.text = element_text(size=12),
                  text = element_text(size=12))   
    
    print(p)
    
  }
}


plot.quad <- function(df, pvals) {
  
  for(i in 1:ncol(df)) {
    hydro <- df[[i]]  
    hydroname <- as.expression(names(df[i]))   # could also ask hydroname to refer to a vector of proper label names
    
    fit.quad <- lm(zCWM ~ hydro + I(hydro^2), data = df)
    catname <- as.factor(c(3,2,2,3,3,2,1,1,1,1,2,2,3,1,3))
    
    p <- qplot(hydro, zCWM, data = df) 
    p = p + geom_point(aes(shape = catname), size =3)
    p <- p + scale_shape_discrete(name = "Hydrological \n class", labels = c("stable winter baseflow", "unpredictable baseflow", "unpredictable intermittent"))
    p <- p + stat_smooth(method = "lm", formula = y ~ x + I(x^2), se=TRUE, col="black") 
    p = p + xlab(hydroname)
    p = p + ylab("AWM wood density (g/cm^3)")
    p = p + annotate("text",                    
                    x=max(hydro)/1.5, y=0.5,
                    label=paste("R^2 = ",signif(summary(fit.quad)$r.squared, 5),
                                "\np.adj =",pvals$quad.padj[i]),
                    size = 4)
    
    p = p + theme(panel.grid.major = element_line(size = .5, color = "grey"),
                  axis.line = element_line(size=.7, color = "black"),
                  legend.position = "bottom",
                  panel.background = element_blank(),      
                  plot.title = element_text(size=12),
                  axis.text = element_text(size=12),
                  text = element_text(size=12))   
    
    print(p)
    
  }
}

plot.exp <- function(df, pvals) {
  
  for(i in 1:ncol(df)) {
    hydro <- df[[i]]  
    hydroname <- as.expression(names(df[i]))   # could also ask hydroname to refer to a vector of proper label names
    
    fit.exp <- lm(zCWM ~ log10(hydro), data = df)
    catname <- as.factor(c(3,2,2,3,3,2,1,1,1,1,2,2,3,1,3))
    
    p <- qplot(hydro, zCWM, data = df) 
    p = p + geom_point(aes(shape = catname), size =3)
    p <- p + scale_shape_discrete(name = "Hydrological \n class", labels = c("stable winter baseflow", "unpredictable baseflow", "unpredictable intermittent"))
    p <- p + stat_smooth(method = "lm", formula = y ~ log10(x), se=TRUE, col="black") 
    p = p + xlab(hydroname)
    p = p + ylab("AWM wood density (g/cm^3)")
    p = p + annotate("text",                    
                     x=max(hydro)/1.5, y=0.5,
                     label=paste("R^2 = ",signif(summary(fit.exp)$r.squared, 5),
                                 "\np.adj =",pvals$exp.padj[i]),
                     size = 4)
    
    p = p + theme(panel.grid.major = element_line(size = .5, color = "grey"),
                  axis.line = element_line(size=.7, color = "black"),
                  legend.position = "bottom",
                  panel.background = element_blank(),      
                  plot.title = element_text(size=12),
                  axis.text = element_text(size=12),
                  text = element_text(size=12))   
    #p = p + theme_tufte()
    print(p)
    
  }
}
