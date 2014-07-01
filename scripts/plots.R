library(ggplot2)
library(ggbiplot)
library(plyr)

source("scripts/functions.R")

# plot CWM/hydro relationships 
# my plot functions are spitting out the wrong R2 values in file names. 
# this happened suddenly and unexpectedly and I can't work out what is happening or why

plot.quad(hydro.quad, hydro.quad.padj)
plot.exp(hydro.exp, hydro.exp.padj)
plot.linear(hydro.linear, hydro.linear.padj)
plot.quad(hydro.nonsignif, hydro.nonsignif.padj) 

# pca (not really worth it's own function...

pdf("output/figures/categories/PCAbiplot.pdf", width = 4.29, height = 2.329)

g <- ggbiplotshape(hydro.signif.pca, obs.scale = 1, var.scale = 1, 
              groups = catname$catnamesfull, ellipse = FALSE)
g <- g + theme_bw()
#p <- p + theme_set(theme_bw(base_size = 8))
g <- g + theme( 
               legend.position = 'none',
               panel.grid.major = element_blank(), # switch off major gridlines
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               axis.line = element_line(size=.2, color = "black")
               )
g <- g + scale_x_continuous(limits = c(-4, 5))
g <- g + scale_y_continuous(limits = c(-2, 3))

print(g)
rm(g)

dev.off()

# generate plots for individual species (this will just dump plots into the working directory, because I'm lazy)

acadeb <- data.frame()

plot.species(WDraw_hydro_aca, acadea, labels = c(
                                                  "NA",
                                                  "NA",
                                                  "NA",
                                                  "NA",
                                                  "NA",
                                                  "NA",
                                                  "NA",
                                                  "NA",
                                                  "NA",
                                                  "NA",
                                                  "NA",
                                                  "NA",
                                                  "NA",
                                                  "NA",
                                                  "NA",
                                                  "NA",
                                                  "NA",
                                                  "NA",
                                                  "NA",
                                                  "NA",
                                                  "NA",
                                                  "NA",
                                                  "(a)",
                                                  "NA",
                                                  "NA",
                                                  "NA")
   )                                               

plot.species(WDraw_hydro_tris, trilau, labels = c(
                                                  "NA",
                                                  "NA",
                                                  "(a)",
                                                 "(b)",
                                                 "NA",
                                                 "NA",
                                                 "NA",
                                                 "NA",
                                                 "NA",
                                                 "NA",
                                                 "NA",
                                                 "NA",
                                                 "NA",
                                                 "NA",
                                                 "NA",
                                                 "NA",
                                                 "NA",
                                                 "NA",
                                                 "NA",
                                                 "NA",
                                                 "NA",
                                                 "NA",
                                                 "NA",
                                                 "NA",
                                                 "NA",
                                                 "NA")                                                 
  )
  
plot.species(WDraw_hydro_cas, cascun, labels = c(
                                                "NA",
                                                "NA",
                                                "NA",
                                                 "NA",
                                                 "NA",
                                                 "NA",
                                                 "NA",
                                                 "NA",
                                                 "NA",
                                                 "NA",
                                                 "NA",
                                                 "NA",
                                                 "NA",
                                                 "NA",
                                                 "NA",
                                                 "(a)",
                                                 "NA",
                                                 "(b)",
                                                 "NA",
                                                 "NA",
                                                 "NA",
                                                 "NA",
                                                 "(c)",
                                                 "NA",
                                                 "NA",
                                                 "NA")
             )

ldply(WDraw_hydro_tris, fitspecies, WDraw_hydro_tris, trilau)
ldply(WDraw_hydro_cas, fitspecies, WDraw_hydro_cas, cascun)
ldply(WDraw_hydro_lep, fitspecies, WDraw_hydro_lep, lepbre)
ldply(WDraw_hydro_aca, fitspecies, WDraw_hydro_aca, acadea)

test <- lm(heart.avg ~ MRateRisenorm, data = WDraw_hydro_aca)
summary(test)


