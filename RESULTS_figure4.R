library(data.table)
library(ggplot2)

## get function createFigure/createFigureBIS/groupFigures,
##     objects greek.label.MM, greek.label.factor, greek.label.lvm, label.statistic, name.statistic, n.statistic
source("FCT.R") 

export <- TRUE

## * Path
path.results <- "./Results"

## * Load data
dtLS.sim.MMtype1 <- readRDS(file.path(path.results,"type1error-simulation-mixedModel.rds"))
dtLS.sim.factortype1 <- readRDS(file.path(path.results,"type1error-simulation-factorModel.rds"))
dtLS.sim.lvmtype1 <- readRDS(file.path(path.results,"type1error-simulation-lvmModel.rds"))

## * Name parameters
dtLS.sim.MMtype1[, link.txt := factor(link,
                                      levels = names(greek.label.MM),
                                      labels = as.character(greek.label.MM))]


dtLS.sim.factortype1[link %in% names(greek.label.factor), link.txt := factor(link,
                                                                             levels = names(greek.label.factor),
                                                                             labels = as.character(greek.label.factor))]



dtLS.sim.lvmtype1[link %in% names(greek.label.lvm), link.txt := factor(link,
                                                                       levels = names(greek.label.lvm),
                                                                       labels = as.character(greek.label.lvm))]



## * Figure 3
## unique(dtLS.sim.factortype1$link)
gg.mm <- createFigure(dtLS.sim.MMtype1,
                      robust = FALSE, link = c("Y2","eta~Gene1Y"),
                      vec.name = name.statistic,
                      vec.label = label.statistic)
gg.factor <- createFigure(dtLS.sim.factortype1,
                          robust = FALSE, link = c("Y2","Y4~eta","eta~Gene1Y","Y1~Gene2Y"),
                          vec.name = name.statistic,
                          vec.label = label.statistic)
gg.lvm <- createFigure(dtLS.sim.lvmtype1,
                       robust = FALSE, link = c("Y2","Y1~Gene2Y","Y4~eta1","eta1~Gene1Y","eta1~eta2","Y1~~Y2"),
                       vec.name = name.statistic,
                       vec.label = label.statistic)
## groupFigures(gg.mm,gg.factor,gg.lvm)

## postscript(file.path(path.figures,"type1error-Wald.eps"), height = 9.5)
pdf(file.path("Figures","type1error-Wald.pdf"), height = 9.25)
groupFigures(gg.mm,gg.factor,gg.lvm)
dev.off()

## * Figure 3 bis (for review only)
ggCoef.mm <- createFigureBIS(data = dtLS.sim.MMtype1,
                             n = c(20,50),
                             robust = FALSE, link = NULL,
                             vec.name = name.statistic,
                             vec.label = label.statistic)
ggCoef.factor <- createFigureBIS(data = dtLS.sim.factortype1,
                                 n = c(20,50),
                                 robust = FALSE, link = NULL,
                                 vec.name = name.statistic,
                                 vec.label = label.statistic)
ggCoef.lvm <- createFigureBIS(data = dtLS.sim.lvmtype1,
                              n = c(20,50),
                              robust = FALSE, link = NULL,
                              vec.name = name.statistic,
                              vec.label = label.statistic)

pdf(file.path("Figures","type1error-Wald2.pdf"), height = 9.25)
groupFigures(ggCoef.mm,ggCoef.factor,ggCoef.lvm, reduce.x = FALSE)
dev.off()

