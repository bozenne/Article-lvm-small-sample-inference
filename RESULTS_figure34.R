## path <- "P:/Cluster/LVMproject/article-smallSampleInference"
## setwd(path) 

library(data.table)
library(ggplot2)
source("FCT.R") ## get function createFigure/createFigureBIS/groupFigures
export <- TRUE

## TODO decrease thickness of lines figures BIS

## * path
path.results <- "./Results"
path.figures <- "./output"

## * load data
dtLS.sim.MMtype1 <- readRDS(file.path(path.results,"type1error-simulation-mixedModel.rds"))
## unique(dtLS.sim.MMtype1$link)
## dtLS.sim.MMtype1[link == "eta~GenderF", link := "eta~Gene1"]
## dtLS.sim.MMtype1[link == "eta~G", linkn := "eta~Age"]

dtLS.sim.factortype1 <- readRDS(file.path(path.results,"type1error-simulation-factorModel.rds"))
## dtLS.sim.factortype1[link == "Y1~Gene1Y", link := "Y1~Gene2Y"]
## dtLS.sim.factortype1[link == "eta~Gene2Y", link := "eta~Gene1Y"]
## dtLS.sim.factortype1[link == "eta~G", link := "eta~Age"]

dtLS.sim.lvmtype1 <- readRDS(file.path(path.results,"type1error-simulation-lvmModel.rds"))
## dtLS.sim.lvmtype1[link == "Y1~Age", link := "Y1~Gene2Y"]
## dtLS.sim.lvmtype1[link == "eta1~G1", link := "eta1~Gene1Y"]
## dtLS.sim.lvmtype1[link == "eta2~G2", link := "eta1~Gender"]

## * name parameters
greek.label.MM <- c(eta = expression(paste(eta, "   (within)")),
                    Y2 = expression(paste(nu[2], "   (within)")),
                    Y3 = expression(paste(nu[3], "   (within)")),
                    "eta~Gene1Y" = expression(paste(gamma[1], "   (between)")),
                    "eta~Age" = expression(paste(gamma[2], "   (between)")))
dtLS.sim.MMtype1[, link.txt := factor(link,
                                      levels = names(greek.label.MM),
                                      labels = as.character(greek.label.MM))]

greek.label.factor <- c("eta" = expression(paste(eta, "   (intercept)")),
                        "Y2" = expression(paste(nu[2], "   (intercept)")),
                        "Y3" = expression(paste(nu[3], "   (intercept)")),
                        "Y4" = expression(paste(nu[4], "   (intercept)")),
                        "Y1~Gene2Y" = expression(paste(k[1], "   (reg. OV~OV)")),
                        "eta~Gene1Y" = expression(paste(gamma[1], "   (reg. LV~OV)")),
                        "eta~Age" = expression(paste(gamma[2], "   (reg. LV~OV)")),
                        "Y2~eta" = expression(paste(lambda[2], "   (loading)")),
                        "Y3~eta" = expression(paste(lambda[3], "   (loading)")),
                        "Y4~eta" = expression(paste(lambda[4], "   (loading)"))
                        )

dtLS.sim.factortype1[link %in% names(greek.label.factor), link.txt := factor(link,
                                                                             levels = names(greek.label.factor),
                                                                             labels = as.character(greek.label.factor))]


greek.label.lvm <- c(Y2 = expression(paste(nu[2], "   (intercept)")),
                     Y3 = expression(paste(nu[3], "   (intercept)")),
                     Y4 = expression(paste(nu[4], "   (intercept)")),
                     Y5 = expression(paste(nu[5], "   (intercept)")),
                     eta1 = expression(paste(eta[1], "   (intercept)")),
                     Z2 = expression(paste(Z[2], "   (intercept)")),
                     Z3 = expression(paste(Z[3], "   (intercept)")),
                     Z4 = expression(paste(Z[4], "   (intercept)")),
                     Z5 = expression(paste(Z[5], "   (intercept)")),
                     eta = expression(paste(eta[2], "   (intercept)")),
                     "Y1~Gene2Y" = expression(paste(k[1], "   (reg. OV~OV)")),
                     "Y2~eta1" = expression(paste(lambda[2], "   (loading)")),
                     "Y3~eta1" = expression(paste(lambda[3], "   (loading)")),
                     "Y4~eta1" = expression(paste(lambda[4], "   (loading)")),
                     "Y5~eta1" = expression(paste(lambda[5], "   (loading)")),
                     "eta1~eta2" = expression(paste(b[1], "   (reg. on LV~LV)")),
                     "eta1~Age" = expression(paste(gamma[2], "   (reg. OV~LV)")),
                     "eta1~Gene1Y" = expression(paste(gamma[1], "   (reg. OV~LV)")),
                     "Z2~eta1" = expression(paste(lambda[6], "   (loading)")),
                     "Z3~eta1" = expression(paste(lambda[7], "   (loading)")),
                     "Z4~eta1" = expression(paste(lambda[8], "   (loading)")),
                     "Z5~eta1" = expression(paste(lambda[9], "   (loading)")),
                     "eta2~Gender" = expression(paste(gamma[3], "   (reg. OV~LV)")),
                     "Y1~~Y2" = expression(paste(sigma[12], "   (covariance)"))
                     )

dtLS.sim.lvmtype1[link %in% names(greek.label.lvm), link.txt := factor(link,
                                                                       levels = names(greek.label.lvm),
                                                                       labels = as.character(greek.label.lvm))]

label.statistic <- c( 
    Ztest = "Gaussian approx.",
    Satt = "Satterthwaite approx.",
    SSC = "bias correction",
    KR = "Satterthwaite approx. with bias correction"
)
name.statistic <- names(label.statistic)
n.statistic <- length(name.statistic)


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

if(export){
    ## postscript(file.path(path.figures,"type1error-Wald.eps"), height = 9.5)
    pdf(file.path(path.figures,"type1error-Wald.pdf"), height = 9.25)
    groupFigures(gg.mm,gg.factor,gg.lvm)
    dev.off()
}

## * Figure 3 bis
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

if(export){
    pdf(file.path(path.figures,"type1error-Wald2.pdf"), height = 9.25)
    groupFigures(ggCoef.mm,ggCoef.factor,ggCoef.lvm, reduce.x = FALSE)
    dev.off()
}

## * Figure 4
ggR.mm <- createFigure(dtLS.sim.MMtype1,
                       robust = TRUE, link = c("Y2","eta~Gene1Y"),
                       vec.name = name.statistic,
                       vec.label = label.statistic)
ggR.factor <- createFigure(dtLS.sim.factortype1,
                           robust = TRUE, link = c("Y2","Y4~eta","eta~Gene1Y","Y1~Gene2Y"),
                           vec.name = name.statistic,
                           vec.label = label.statistic) + scale_y_continuous(breaks = scales::pretty_breaks(n = 5))
ggR.lvm <- createFigure(dtLS.sim.lvmtype1,
                        robust = TRUE, link = c("Y2","Y1~Gene2Y","Y4~eta1","eta1~Gene1Y","eta1~eta2","Y1~~Y2"),
                        vec.name = name.statistic,
                        vec.label = label.statistic) + scale_y_continuous(breaks = scales::pretty_breaks(n = 5))

if(export){
    ## postscript(file.path(path.figures,"type1error-Wald-robust.eps"))
    pdf(file.path(path.figures,"type1error-Wald-robust.pdf"), height = 9.25)
    groupFigures(ggR.mm,ggR.factor,ggR.lvm)
    dev.off()
}

##  groupFigures(ggR.mm,ggR.factor,ggR.lvm)

## * Figure 4 bis
ggCoefR.mm <- createFigureBIS(data = dtLS.sim.MMtype1,
                              n = c(20,50),
                              robust = TRUE, link = NULL,
                              vec.name = name.statistic,
                              vec.label = label.statistic)
ggCoefR.factor <- createFigureBIS(data = dtLS.sim.factortype1,
                                  n = c(20,50),
                                  robust = TRUE, link = NULL,
                                  vec.name = name.statistic,
                                  vec.label = label.statistic)
ggCoefR.lvm <- createFigureBIS(data = dtLS.sim.lvmtype1,
                               n = c(20,50),
                               robust = TRUE, link = NULL,
                               vec.name = name.statistic,
                               vec.label = label.statistic)

if(export){
    pdf(file.path(path.figures,"type1error-Wald2-robust.pdf"), height = 9.25)
    groupFigures(ggCoefR.mm,ggCoefR.factor,ggCoefR.lvm, reduce.x = FALSE)
    dev.off()
}

## * Tables reviewers
greek.label.factor2 <- c("eta" = "$\\alpha$",
                         "Y2" = "$\\nu_2$",
                         "Y3" = "$\\nu_3$",
                         "Y4" = "$\\nu_4$",
                         "Y2~Gene1Y" = "$k_1$",
                         "eta~Gene2Y" = "$\\gamma_1$",
                         "eta~Age" = "$\\gamma_2$",
                         "Y2~eta" = "$\\lambda_2$",
                         "Y3~eta" = "$\\lambda_3$",
                         "Y4~eta" = "$\\lambda_4$")

## ** comparisons


dtLS.sim.Comptype1 <- readRDS(file.path(path.results,"type1error-simulation-comparison.rds"))

dtLS.sim.Comptype1[link %in% names(greek.label.factor2), link.txt := factor(link,
                                                                            levels = names(greek.label.factor2),
                                                                            labels = as.character(greek.label.factor2))]
dtLS.sim.Comptype1[estimator == "IV", estimator := "IV (lava)"]
dtLS.sim.Comptype1[estimator == "IVlavaan", estimator := "IV (lavaan)"]
dtLS.sim.Comptype1[estimator == "robustML", estimator := "robust ML"]


addtorow <- list()
addtorow$pos <- list(10)
addtorow$command <- "[4mm]"
table.R1 <- dcast(dtLS.sim.Comptype1,
                  formula = link.txt+n ~ estimator, value.var = "type1")[n%in%c(20,50)]
setnames(table.R1, old = c("link.txt"), new = c("parameter"))
table.R1 <- table.R1[order(table.R1$n),]
table.R1$n[duplicated(table.R1$n)] <- ""

order.col <- c("parameter","n","ML","robust ML","GLS","WLS","IV (lava)","IV (lavaan)")
print(xtable::xtable(table.R1[,.SD,.SDcols = order.col], label = "tab:comparison", caption = "Comparison of the type 1 error of Wald tests for various estimation methods in scenario (b) under a correctly specified model. No small sample correction is used. Robust ML corresponds to the use of robust Wald tests. NA indicates that the estimator never converged in the simulations and so the Wald statistic could not be computed. The R packages lava and MIIVsem were used for IV estimation. GLS and WLS were carried out using the R package lavaan.", digits = 3),
      add.to.row =  addtorow, NA.string="NA",
      include.rownames = FALSE, booktabs = TRUE, sanitize.text.function = function(x){x})


## ** misspecification correction small samples (student)
dtLS.sim.IV1type1 <- readRDS(file = file.path("Results","type1error-simulation-IV1.rds"))

dtLS.sim.IV1type1[link %in% names(greek.label.factor2), link.txt := factor(link,
                                                                            levels = names(greek.label.factor2),
                                                                           labels = as.character(greek.label.factor2))]
dtLS.sim.IV1type1[estimator == "IV", estimator := "IV (lava)"]
dtLS.sim.IV1type1[estimator == "IVlavaan", estimator := "IV (lavaan)"]
dtLS.sim.IV1type1[estimator == "robustML", estimator := "robust ML"]
dtLS.sim.IV1type1[estimator == "robustMLC", estimator := "corrected robust ML"]

table.R2 <- dcast(dtLS.sim.IV1type1, formula = link.txt+n ~ estimator, value.var = "type1")[n==20]

names(table.R2)[1] <- "parameter"
order.col <- c("parameter","robust ML","corrected robust ML","IV (lava)","IV (lavaan)")
table.R2[,n:=NULL]

print(xtable::xtable(table.R2[,.SD,.SDcols = order.col], label = "tab:student", caption = "Type 1 error of Wald tests in a misspecified latent variable model (residuals following a Student's t-distribution) for a sample size of 20.", digits = 3), NA.string="NA",
      include.rownames = FALSE, booktabs = TRUE, sanitize.text.function = function(x){x})

## ** misspecification correction small samples (misspecified covariance structure)
dtLS.sim.IV3type1 <- readRDS(file = file.path("Results","type1error-simulation-IV3.rds"))

dtLS.sim.IV3type1[link %in% names(greek.label.factor2), link.txt := factor(link,
                                                                           levels = names(greek.label.factor2),
                                                                           labels = as.character(greek.label.factor2))]
dtLS.sim.IV3type1[estimator == "IV", estimator := "IV (lava)"]
dtLS.sim.IV3type1[estimator == "IVlavaan", estimator := "IV (lavaan)"]
dtLS.sim.IV3type1[estimator == "robustML", estimator := "robust ML (uncorrected)"]
dtLS.sim.IV3type1[estimator == "robustMLC", estimator := "robust ML (df=non-robust)"]
dtLS.sim.IV3type1[estimator == "robustMLC2", estimator := "robust ML (df=Pan)"]

table.R3 <- dcast(dtLS.sim.IV3type1, formula = link.txt+n ~ estimator, value.var = "type1")[n %in% c(20,100)]
setkeyv(table.R3, "n")
names(table.R3)[1] <- "parameter"
order.col <- c("parameter","n","robust ML (uncorrected)","robust ML (df=non-robust)","robust ML (df=Pan)","IV (lava)","IV (lavaan)")
table.R3[,n := as.character(n)]
table.R3[duplicated(n),n := ""]

addtorow <- list()
addtorow$pos <- list(0,0,10)
addtorow$command <- c("\\multirow{2}{*}{parameter} & \\multirow{2}{*}{n} & robust ML & robust ML & robust ML & IV & IV \\\\",
                      " & & (uncorrected) & (df=non-robust) & (df=Pan) &  (lava) &  (lavaan) \\\\",
                      "[4mm]")

print(xtable::xtable(table.R3[,.SD,.SDcols = order.col], label = "tab:miscov", caption = "Type 1 error of Wald tests in a misspecified latent variable model (incorrect covariance structure)", digits = 3), NA.string="NA",
      add.to.row = addtorow,
      include.colnames = FALSE,
      include.rownames = FALSE,
      booktabs = TRUE, sanitize.text.function = function(x){x})


## ** timings
library(microbenchmark)
bench.MM <- readRDS(file.path(path.results,"speed-Algo2-MM.rds"))
bench.Factor <- readRDS(file.path(path.results,"speed-Algo2-Factor.rds"))
bench.Lvm <- readRDS(file.path(path.results,"speed-Algo2-LVM.rds"))
class(bench.MM)

## args(microbenchmark:::summary.microbenchmark)
speed.Algo2 <- rbind(cbind(scenario = "a", summary(bench.MM, unit = "s")),
                     cbind(scenario = "b", summary(bench.Factor, unit = "s")),
                     cbind(scenario = "c", summary(bench.Lvm, unit = "s")))
names(speed.Algo2)[names(speed.Algo2)=="expr"] <- "n"
keep.col <- c("scenario","n","min","mean","median","max")
speed.Algo2 <- speed.Algo2[,keep.col]
speed.Algo2$scenario <- as.character(speed.Algo2$scenario)
speed.Algo2$scenario[duplicated(speed.Algo2$scenario)] <- ""


addtorow <- list()
addtorow$pos <- list(9,18)
addtorow$command <- c("[4mm]","[4mm]")
print(xtable::xtable(speed.Algo2, label = "tab:speed", caption = "Computation time for algorithm 2 over 50 repetitions (1 CPU: AMD Opteron(tm) Processor 6380, 2.44 GHz).", digits = 3),
      add.to.row =  addtorow, 
      include.rownames = FALSE, booktabs = TRUE, sanitize.text.function = function(x){x})

