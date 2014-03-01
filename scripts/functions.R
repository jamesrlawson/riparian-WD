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





#mainDir <- "C:/Users/JLawson/Desktop/stuff/data/analysis/R"
#subDir <- "TGA all traits/output_seedmass/tp"

#dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

#setwd("C:/Users/JLawson/Desktop/stuff/data/analysis/R")

# use sprintf within the jpeg graphics device to add pvals and R2's to our file names  

#jpeg(sprintf("TGA all traits/output_seedmass/tp/tp_%s_p-%s_r2-%s.jpg", k, pvaluerounded, Rsquared), quality = 100, width = 600, height = 500)






plot.linear <- function(df, pvals) {

  setwd("C:/Users/JLawson/Desktop/stuff/data/analysis/R/WDmeans")
  
  for(i in 1:ncol(df)) {
    hydro <- df[[i]]  
    hydroname <- as.expression(names(df[i]))   # could also ask hydroname to refer to a vector of proper label names
        
    fit.linear <- lm(zCWM ~ hydro, data = df)
    catname <- as.factor(c(3,2,2,3,3,2,1,1,1,1,2,2,3,1,3))
    
    padj <- pvals$linear.padj[i]
    r2 <- signif(summary(fit.linear)$r.squared, 5)

    png(sprintf("output/figures/%s_p-%s_r2-%s.png", hydroname, padj, r2), width = 600, height = 500)
    #on.exit(dev.off())
    
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
    dev.off()
  }
}


plot.quad <- function(df, pvals) {
  
  for(i in 1:ncol(df)) {
    hydro <- df[[i]]  
    hydroname <- as.expression(names(df[i]))   # could also ask hydroname to refer to a vector of proper label names
    
    fit.quad <- lm(zCWM ~ hydro + I(hydro^2), data = df)
    catname <- as.factor(c(3,2,2,3,3,2,1,1,1,1,2,2,3,1,3))
    
    padj <- pvals$quad.padj[i]
    r2 <- signif(summary(fit.quad)$r.squared, 5)
    
    png(sprintf("output/figures/%s_p-%s_r2-%s.png", hydroname, padj, r2), width = 600, height = 500)
    #on.exit(dev.off())
    
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
    dev.off()
  }
}

plot.exp <- function(df, pvals) {
  
  for(i in 1:ncol(df)) {
    hydro <- df[[i]]  
    hydroname <- as.expression(names(df[i]))   # could also ask hydroname to refer to a vector of proper label names
    
    fit.exp <- lm(zCWM ~ log10(hydro), data = df)
    catname <- as.factor(c(3,2,2,3,3,2,1,1,1,1,2,2,3,1,3))
    
    padj <- pvals$exp.padj[i]
    r2 <- signif(summary(fit.exp)$r.squared, 5)
    
    pdf(sprintf("output/figures/%s_p-%s_r2-%s.pdf", hydroname, padj, r2), width = 5.5, height = 3.2)
    
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
    
    p = p + theme(panel.grid.major = element_blank(),
                  axis.line = element_line(size=.2, color = "black"),
                  legend.position = "bottom",
                  legend.text = element_text(size=8),
                  legend.title = element_text(size=10),
                  panel.background = element_blank(),      
                  axis.text = element_text(size=8),
                  text = element_text(size=8))   
    #p = p + theme_tufte()
    print(p)
    dev.off()
   #ggsave(sprintf("output/figures/%s_p-%s_r2-%s.pdf", hydroname, padj, r2), width = 50, height = 40, units=c("mm"), p, scale=1)
    
  }
}


stderr <- function(x) sd(x)/sqrt(length(x))

plot.means <- function(df) {
  pd <- position_dodge(.001)
#  png(sprintf("output/figures/categories/%s.png", df), width = 480, height = 480, units="px")
  p <- ggplot(df, aes(x=category, y=mean)) + 
    geom_errorbar(aes(ymin=mean-stderr, ymax=mean+stderr), colour="black", width=.1, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
    xlab(df$labels[2]) +
    ylab(df$labels[3]) +
    ggtitle(df$labels[1]) +
    theme_minimal() +
    theme(legend.justification=c(1,0), legend.position=c(1,0)) # Position legend in bottom right
  print(p)
#  dev.off()
  
}
