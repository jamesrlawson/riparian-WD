library(ggplot2)
library(ggbiplot)
library(plyr)

source("scripts/functions.R")

# plot CWM/hydro relationships 

plot.quad(hydro.quad, hydro.quad.padj)
plot.linear(hydro.linear, hydro.linear.padj)
plot.quad(hydro.nonsignif, hydro.nonsignif.padj) 

# pca (not really worth it's own function...)

svg("output/figures/PCAbiplot2.svg", width = 4.29, height = 2.329, pointsize=8)

#g <- ggbiplotshape(hydro.signif.pca, obs.scale = 1, var.scale = 1, groups = catname$catnamesfull, ellipse = FALSE)
g <- ggbiplotshape(hydro.signif.pca, obs.scale = 1, var.scale = 1, ellipse = FALSE)
g <- g + theme_bw()
g <- g + theme_set(theme_bw(base_size = 8))
g <- g + theme( 
               legend.position = 'none',
               axis.text = element_text(size = rel(0.8)),
               axis.title = element_text(size = rel(0.8)),
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


# plot for CWV vs PC1 'unimodal distribution'

CWV.lm <- lm(CWV$CWV ~ PC1)
summary(CWV.lm)
CWV.lm.quad<- lm(CWV$CWV ~ PC1 + I(PC1^2))
summary(CWV.lm.quad)

CWV_PC1 <- data.frame(CWV$CWV, PC1)
colnames(CWV_PC1) <- c("CWV", "PC1")

svg("output/figures/CWV.svg", width = 2.795, height = 2.329, pointsize=8)

p <- ggplot(CWV_PC1, aes(x = PC1, y = CWV))
p <- p + geom_point(size = 2)
p <- p + stat_smooth(size = 0.2, fullrange = TRUE, method = "lm", formula = y ~ x + I(x^2), se=TRUE, col="black", alpha = 0.2) 
p <- p + xlab("PC1")
p <- p + ylab(expression(paste("CWV wood density ", "(g / ", cm^3,")")))
#p <- p + ylim(0.45, 0.75)
p <- p + theme_bw() 
p <- p + theme_set(theme_bw(base_size = 8))
p <- p + theme(legend.position = "none",
               axis.text = element_text(size = rel(0.8)),
               axis.title.y = element_text(hjust=0.35),
               axis.title.x = element_text(vjust=0.35),
               panel.border = element_blank(),
               panel.grid.minor = element_blank(),
               panel.grid.major = element_blank(),
               axis.line = element_line(size=.2, color = "black")
)

print(p)
dev.off()




svg("output/figures/HSPeak_linear.svg", width = 2.795, height = 2.329, pointsize=8)

p <- ggplot(hydroCWM, aes(x = HSPeak, y = CWM))
#p <- p + geom_point(aes(shape = catname$catnamesfull), size = 2)
p <- p + geom_point(size = 2)
p <- p + stat_smooth(size = 0.2, fullrange = TRUE, method = "lm", formula = y ~ x, se=TRUE, col="black", alpha = 0.2) 
p <- p + xlab("HSPeak")
p <- p + ylab(expression(paste("Mean wood density ", "(g / ", cm^3,")")))
#p <- p + ylim(0.45, 0.75)
#    p <- p + annotate("text", x = min(hydro) * 1.1, y = 0.74, label = sprintf("%s", pvals$figlabel[i]), size = 3)
p <- p + theme_bw() 
p <- p + theme_set(theme_bw(base_size = 8))
p <- p + theme(legend.position = "none",
               axis.text = element_text(size = rel(0.8)),
               axis.title.y = element_text(hjust=0.35),
               axis.title.x = element_text(vjust=0.35),
               panel.border = element_blank(),
               panel.grid.minor = element_blank(),
               panel.grid.major = element_blank(),
               axis.line = element_line(size=.2, color = "black")
)

print(p)
dev.off()




# generate plots for individual species (this will just dump plots into the working directory, because I'm lazy)

WDdata <- read.csv("data/WDdata.csv", header=TRUE)
WDraw_hydro <- merge(WDdata, hydro)

WDraw_hydro_cas <- WDraw_hydro[WDraw_hydro$species =="Casuarina cunninghamiana",]
WDraw_hydro_cas$species <- NULL
WDraw_hydro_cas$catname <- NULL
WDraw_hydro_cas$category <- NULL

plot.species(WDraw_hydro_cas, cascun)

#acadeb <- data.frame()

#plot.species(WDraw_hydro_aca, acadea, labels = c(
#                                                  "NA",
#                                                  "NA",
#                                                  "NA",
#                                                  "NA",
#                                                  "NA",
#                                                  "NA",
#                                                 "NA",
#                                                  "NA",
#                                                  "NA",
#                                                  "NA",
#                                                  "NA",
#                                                  "NA",
#                                                 "NA",
#                                                  "NA",
 #                                                 "NA",
  #                                                "(a)",
   #                                               "NA",
    #                                              "NA",
     #                                             "NA",
      #                                            "NA",
       #                                           "NA",
        #                                          "NA",
         #                                         "(a)",
          #                                        "NA",
           #                                       "NA",
            #                                      "NA")
   #)                                               

#plot.species(WDraw_hydro_tris, trilau, labels = c(
 #                                                 "NA",
  #                                                "NA",
   #                                               "NA",
    #                                             "NA",
     #                                            "(a)",
      #                                           "NA",
       #                                          "NA",
        #                                         "NA",
         #                                        "NA",
          #                                       "NA",
           #                                      "NA",
            #                                     "NA",
             #                                    "NA",
              #                                   "NA",
               #                                  "NA",
                #                                 "NA",
                 #                                "NA",
                  #                               "NA",
                   #                              "NA",
                    #                             "NA",
                     #                            "NA",
                      #                           "NA",
                       #                          "NA",
                        #                         "NA",
                         #                        "NA",
                          #                       "NA")                                                 
  #)
  
WDdata <- read.csv("data/WDdata.csv", header=TRUE)
WDraw_hydro <- merge(WDdata, hydro)

WDraw_hydro_cas <- WDraw_hydro[WDraw_hydro$species =="Casuarina cunninghamiana",]
WDraw_hydro_cas$species <- NULL
WDraw_hydro_cas$catname <- NULL
WDraw_hydro_cas$category <- NULL

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
                                                 "(d)",
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
                                                 "(c)",
                                                 "NA")
             )

#ldply(WDraw_hydro_tris, fitspecies, WDraw_hydro_tris, trilau)
ldply(WDraw_hydro_cas, fitspecies, WDraw_hydro_cas, cascun)
#ldply(WDraw_hydro_lep, fitspecies, WDraw_hydro_lep, lepbre)
#ldply(WDraw_hydro_aca, fitspecies, WDraw_hydro_aca, acadea)

#test <- lm(WD ~ MRateRise, data = WDraw_hydro_aca)
#summary(test)


