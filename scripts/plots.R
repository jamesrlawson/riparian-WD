library(ggplot2)
library(ggbiplot)

source("scripts/functions.R")

# plot CWM/hydro relationships

plot.quad(hydro.quad, hydro.quad.padj)
plot.exp(hydro.exp, hydro.exp.padj)
plot.linear(hydro.linear, hydro.linear.padj)
plot.quad(hydro.nonsignif, hydro.nonsignif.padj)

# plot categorical comparisons

plot.means(raw_cats.df)
plot.means(CWM_cats.df)

# pca (not really worth it's own function...

g <- ggbiplotshape(hydro.signif.pca, obs.scale = 1, var.scale = 1, 
              groups = catname$catnamesfull, ellipse = TRUE)
g <- g + theme_minimal()
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'bottom')
g <- g + scale_x_continuous(limits = c(-4, 5))

print(g)
rm(g)

# generate plots for individual species (this will just dump plots into the working directory, because I'm lazy)

plot.species(WDraw_hydro_tris)
plot.species(WDraw_hydro_cas)
plot.species(WDraw_hydro_lep)
plot.species(WDraw_hydro_aca)