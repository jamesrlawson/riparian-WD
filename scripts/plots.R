library(ggplot2)
library(ggbiplot)

source("scripts/functions.R")

# plot CWM/hydro relationships 
# my plot functions are spitting out the wrong R2 values in file names. 
# this happened suddenly and unexpectedly and I can't work out what is happening or why

plot.quad(hydro.quad, hydro.quad.padj)
plot.exp(hydro.exp, hydro.exp.padj)
plot.linear(hydro.linear, hydro.linear.padj)
plot.quad(hydro.nonsignif, hydro.nonsignif.padj)


# plot categorical comparisons

plot.means(compareMeans)


# pca (not really worth it's own function...

pdf("output/figures/categories/PCAbiplot.pdf", width = 4.29, height = 2.329)

g <- ggbiplotshape(hydro.signif.pca, obs.scale = 1, var.scale = 1, 
              groups = catname$catnamesfull, ellipse = FALSE)
g <- g + theme_bw()
p <- p + theme_set(theme_bw(base_size = 8))
g <- g + theme( 
               legend.position = 'none',
               panel.grid.major = element_blank(), # switch off major gridlines
               panel.grid.minor = element_blank())
g <- g + scale_x_continuous(limits = c(-4, 5))
g <- g + scale_y_continuous(limits = c(-2, 3))

print(g)
rm(g)

dev.off()

# generate plots for individual species (this will just dump plots into the working directory, because I'm lazy)

plot.species(WDraw_hydro_tris, trilau)
plot.species(WDraw_hydro_cas, cascun)
plot.species(WDraw_hydro_lep, lepbre)
plot.species(WDraw_hydro_aca, acadeb)







