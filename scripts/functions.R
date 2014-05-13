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

  setwd("C:/Users/JLawson/Desktop/stuff/data/analysis/R/WDmeans")
    
  for(i in 1:ncol(df)) {
    hydro <- df[[i]]  
    hydroname <- as.expression(names(df[i]))   # could also ask hydroname to refer to a vector of proper label names
        
    fit.linear <- lm(zCWM ~ hydro, data = df)
    catname <- as.factor(c(3,2,2,3,3,2,1,1,1,1,2,2,3,1,3))
    
    padj <- pvals$linear.padj[i]
    r2 <- signif(summary(fit.linear)$r.squared, 5)
    
    label <- substitute(deparse(pvals$figlabel[i]))

    pdf(sprintf("output/figures/%s_p-%s_r2-%s.pdf", hydroname, padj, r2), width = 2.795, height = 2.329)
    #on.exit(dev.off())
    
    df$hydro <- hydro
    df$catname <- catname
    
    p <- ggplot(df, aes(x = hydro, y = zCWM))
    p <- p + geom_point(aes(shape = catname), size = 2)
    p <- p + stat_smooth(size = 0.2, fullrange = TRUE, method = "lm", formula = y ~ x, se=TRUE, col="black", alpha = 0.2) 
    p <- p + xlab(hydroname)
    p <- p + ylab("Mean wood density (g/cm^3)")
    p <- p + ylim(0.45, 0.75)
    p <- p + annotate("text", x = min(hydro) * 1.1, y = 0.74, label = sprintf("%s", pvals$figlabel[i]), size = 3)
    p <- p + theme_bw() 
    p <- p + theme_set(theme_bw(base_size = 8))
    p <- p + theme(legend.position = "none",
                   axis.text = element_text(size = rel(0.8)),
                   axis.title.y = element_text(hjust=0.35),
                   axis.title.x = element_text(vjust=0.35),
                   panel.border = element_blank(),
                   panel.grid.minor = element_blank()
    )
    
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
    
    pdf(sprintf("output/figures/%s_p-%s_r2-%s.pdf", hydroname, padj, r2), width = 2.795, height = 2.329)
    #on.exit(dev.off())
    df$hydro <- hydro
    df$catname <- catname
    
    p <- ggplot(df, aes(x = hydro, y = zCWM))
    p <- p + geom_point(aes(shape = catname), size = 2)
    p <- p + stat_smooth(size = 0.2, fullrange = TRUE, method = "lm", formula = y ~ x + I(x^2), se=TRUE, col="black", alpha = 0.2) 
    p <- p + xlab(hydroname)
    p <- p + ylab("Mean wood density (g/cm^3)")
    p <- p + ylim(0.45, 0.75)
    p <- p + annotate("text", x = min(hydro) * 1.1, y = 0.74, label = sprintf("%s", pvals$figlabel[i]), size = 3)
    p <- p + theme_bw() 
    p <- p + theme_set(theme_bw(base_size = 8))
    p <- p + theme(legend.position = "none",
                   axis.text = element_text(size = rel(0.8)),
                   axis.title.y = element_text(hjust=0.35),
                   axis.title.x = element_text(vjust=0.35),
                   panel.border = element_blank(),
                   panel.grid.minor = element_blank()
    )
    
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
    
    pdf(sprintf("output/figures/%s_p-%s_r2-%s.pdf", hydroname, padj, r2), width = 2.795, height = 2.329)
    #on.exit(dev.off())
    
    df$hydro <- hydro
    df$catname <- catname
    
    p <- ggplot(df, aes(x = hydro, y = zCWM))
    p <- p + geom_point(aes(shape = catname), size = 2)
    p <- p + stat_smooth(size = 0.2, fullrange = TRUE, method = "lm", formula = y ~ log10(x), se=TRUE, col="black", alpha = 0.2) 
    p <- p + xlab(hydroname)
    p <- p + ylab("Mean wood density (g/cm^3)")
    p <- p + ylim(0.45, 0.75)
    p <- p + annotate("text", x = min(hydro) * 1.1, y = 0.74, label = sprintf("%s", pvals$figlabel[i]), size = 3)
    p <- p + theme_bw() 
    p <- p + theme_set(theme_bw(base_size = 8))
    p <- p + theme(legend.position = "none",
                   axis.text = element_text(size = rel(0.8)),
                   axis.title.y = element_text(hjust=0.35),
                   axis.title.x = element_text(vjust=0.35),
                   panel.border = element_blank(),
                   panel.grid.minor = element_blank()
              )
    

    print(p)
    dev.off()
   #ggsave(sprintf("output/figures/%s_p-%s_r2-%s.pdf", hydroname, padj, r2), width = 50, height = 40, units=c("mm"), p, scale=1)
    
  }
}


stderr <- function(x) sd(x)/sqrt(length(x))



plot.means <- function(df) {
   
  outDir = "output/figures/categories"
  dir.create(outDir, recursive=TRUE, showWarnings = FALSE)
  
  pdf(sprintf("%s/comparemeans.pdf", outDir), width = 4.29, height = 2.329)
    
  p <-ggplot(df, aes(x=category, y=mean)) + 
      geom_point(stat = "identity", size=2) + 
      geom_errorbar(aes(ymin=mean-stderr, ymax=mean+stderr), size = 0.1, width=.1) +
      facet_grid(. ~ datatype) +
      xlab(df$labels[2]) +
      ylab(df$labels[3]) +
      scale_y_continuous(limits = c(0.5, 0.7)) +
      ggtitle(df$labels[1]) +
      theme_bw() +
      theme_set(theme_bw(base_size = 8)) +
      theme(legend.position = "none",
                   axis.text = element_text(size = rel(0.8)),
                   axis.title.y = element_text(hjust=0.35),
                   axis.title.x = element_text(vjust=0.35),
                   panel.border = element_blank(),
                   panel.grid.minor = element_blank(),
                   strip.background = element_blank()
            )
  print(p)
  dev.off()
  
}

ggbiplotshape <- function (pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE,
                          obs.scale = 1 - scale, var.scale = scale, groups = NULL,
                          ellipse = FALSE, ellipse.prob = 0.68, labels = NULL, labels.size = 3,
                          alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69,
                          varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE,
                          ...) {

  library(ggplot2)
  library(plyr)
  library(scales)
  library(grid)
  stopifnot(length(choices) == 2)
  if (inherits(pcobj, "prcomp")) {
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$rotation
  }
  else if (inherits(pcobj, "princomp")) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$loadings
  }
  else if (inherits(pcobj, "PCA")) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1/(d * nobs.factor), FUN = "*")
    v <- sweep(pcobj$var$coord, 2, sqrt(pcobj$eig[1:ncol(pcobj$var$coord),
                                                  1]), FUN = "/")
  }
  else {
    stop("Expected a object of class prcomp, princomp or PCA")
  }
  df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale,
                              FUN = "*"))
  v <- sweep(v, 2, d^var.scale, FUN = "*")
  df.v <- as.data.frame(v[, choices])
  names(df.u) <- c("xvar", "yvar")
  names(df.v) <- names(df.u)
  if (pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  v.scale <- rowSums(v^2)
  df.v <- r * df.v/sqrt(max(v.scale))
  if (obs.scale == 0) {
    u.axis.labs <- paste("standardized PC", choices, sep = "")
  }
  else {
    u.axis.labs <- paste("PC", choices, sep = "")
  }
  u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained var.)",
                                            100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  if (!is.null(labels)) {
    df.u$labels <- labels
  }
  if (!is.null(groups)) {
    df.u$groups <- groups
  }
  if (varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  }
  else {
    df.v$varname <- rownames(v)
  }
  df.v$angle <- with(df.v, (180/pi) * atan(yvar/xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar))/2)
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + xlab(u.axis.labs[1]) +
    ylab(u.axis.labs[2]) + coord_equal()
  if (var.axes) {
    if (circle) {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi,
                                                length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r *
                             sin(theta))
      g <- g + geom_path(data = circle, color = muted("white"),
                         size = 1/2, alpha = 1/3)
    }
    g <- g + geom_segment(data = df.v, aes(x = 0, y = 0,
                                           xend = xvar, yend = yvar), arrow = arrow(length = unit(1/2,
                                                                                                  "picas")), color = muted("red"))
  }
  if (!is.null(df.u$labels)) {
    if (!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups),
                         size = labels.size)
    }
    else {
      g <- g + geom_text(aes(label = labels), size = labels.size)
    }
  }
  else {
    if (!is.null(df.u$groups)) {
      g <- g + geom_point(aes(shape = groups, size = 1), alpha = alpha)      
    }
    else {
      g <- g + geom_point(aes(size=1), alpha = alpha)
    }
  }
  if (!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    ell <- ddply(df.u, "groups", function(x) {
      if (nrow(x) < 2) {
        return(NULL)
      }
      else if (nrow(x) == 2) {
        sigma <- var(cbind(x$xvar, x$yvar))
      }
      else {
        sigma <- diag(c(var(x$xvar), var(x$yvar)))
      }
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2,
                       mu, FUN = "+"), groups = x$groups[1])
    })
    names(ell)[1:2] <- c("xvar", "yvar")
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }
#  if (var.axes) {
#    g <- g + geom_text(data = df.v, aes(label = varname,
#                                        x = xvar, y = yvar, angle = angle, hjust = hjust),
#                       color = "black", size = varname.size,
#                       position=position_jitter(h=0.5,w=0))
#  }
  return(g)
}


plot.species <- function(df, species, labels) {
  
  setwd("C:/Users/JLawson/Desktop/stuff/data/analysis/R/WDmeans")
  
  for(i in 1:ncol(df)) {
    hydro <- df[[i]]  
    hydroname <- as.expression(names(df[i]))   # could also ask hydroname to refer to a vector of proper label names
    
    fit.linear <- lm(heart.avg ~ hydro, data = df)
    catname <- as.factor(c(3,2,2,3,3,2,1,1,1,1,2,2,3,1,3))
    
    pval<- anova(fit.linear)[1,"Pr(>F)"]
    r2 <- signif(summary(fit.linear)$r.squared, 5)
    
    stats <- cbind(as.character(hydroname), pval, r2)
    
    sp = deparse(substitute(species))
        
    outDir = sprintf("output/figures/species/%s", sp)
    dir.create(outDir, recursive=TRUE, showWarnings = FALSE)
                
    pdf(sprintf("%s/%s.pdf", outDir, hydroname), width = 2.795, height = 2.329)
    #on.exit(dev.off())
    
    df$hydro <- hydro
    
    p <- ggplot(df, aes(x = hydro, y = heart.avg))
    p <- p + geom_point(ysize = 2)
    p <- p + stat_smooth(size = 0.2, method = "lm", formula = y ~ x, se=TRUE, col="black", alpha = 0.2) 
    p <- p + xlab(hydroname)
    p <- p + ylab("Mean wood density (g/cm^3)")
    p <- p + ylim(0.45, 0.75)
    p <- p + annotate("text", x = min(hydro) * 1.1, y = 0.74, label = sprintf("%s", labels[i]), size = 3)
    p <- p + theme_bw() 
    p <- p + theme_set(theme_bw(base_size = 8))
    p <- p + theme(legend.position = "none",
                   axis.text = element_text(size = rel(0.8)),
                   axis.title.y = element_text(hjust=0.35),
                   axis.title.x = element_text(vjust=0.35),
                   panel.border = element_blank(),
                   panel.grid.minor = element_blank()
    )
    
    print(p) 
    dev.off()
  }
  
}

fitspecies <- function(hydro, df, species) {
  
  fit.linear  <- lm(heart.avg ~ hydro, data = df)
  fit.quad    <- lm(heart.avg ~ hydro + I(hydro ^2), data = df)  
#  fit.exp     <- lm(heart.avg ~ log10(hydro), subset = hydro > 0, data = df) # THIS IS THE LINE THAT HANDLES INF's!
  
  # p-values
  linear.pval <- anova(fit.linear)[1,"Pr(>F)"]
  quad.pval   <- anova(fit.quad)[1,"Pr(>F)"]
#  exp.pval    <- anova(fit.exp)[1,"Pr(>F)"]
  
  # r-squared
  linear.r2   <- summary(fit.linear)$r.squared
  quad.r2     <- summary(fit.quad)$r.squared
#  exp.r2      <- summary(fit.exp)$r.squared  
  
  
  fitspecies.df <- cbind(linear.pval, linear.r2, quad.pval, quad.r2)
  
  sp = deparse(substitute(species))
  
  write.table(fitspecies.df, file=sprintf("output/figures/species/%s/stats.csv", sp), append=TRUE, sep=",", col.names=FALSE, row.names=FALSE)
  
}

