## path <- "P:/Cluster/LVMproject/article-smallSampleInference"
## setwd(path)

library(data.table)
source("FCT.R") ## get function sinkDirectory


## * path
path.results <- "./Results"
path.simulation.mixedModel <- "./Results/simulation-mixedModel"
path.simulation.factorModel <- "./Results/simulation-factorModel"
path.simulation.lvm <- "./Results/simulation-lvm"
path.simulation.comparison <- "./Results/comparison-ML-IV-GLS"
path.simulation.IV1 <- "./Results/IV-non-normal"
path.simulation.IV2 <- "./Results/IV-non-normal2"
path.simulation.IV3 <- "./Results/IV-non-normal3"
path.illustration.mixedModel <- "./Results/illustration-mixedModel"
path.illustration.factorModel <- "./Results/illustration-factorModel"
path.illustration.lvm <- "./Results/illustration-lvm"

## * Simulation
cat("Simulation: ")

## ** Mixed model - type 1 error
cat(" - mixed model (type 1 error) \n")

dt.sim.MMtype1 <- sinkDirectory(path.simulation.mixedModel,
                                string.keep = "type1error", string.exclude = "(tempo)")
## head(dt.sim.MMtype1)
## dt.sim.MMtype1[link == link[1], .(nb=.N,pc=100*.N/20000,warning=sum(warning),niterMax=max(niter)), by = "n"]
##      n    nb      pc warning niterMax
## 1:  20 19929  99.645       0        8
## 2:  30 19993  99.965       0        6
## 3:  50 20000 100.000       0        5
## 4:  75 20000 100.000       0        5
## 5: 100 20000 100.000       0        4
## 6: 150 20000 100.000       0        4
## 7: 200 20000 100.000       0        4
## 8: 300 20000 100.000       0        3
## 9: 500 20000 100.000       0        3

dtL.sim.MMtype1 <- melt(dt.sim.MMtype1,
                        id.vars = c("n","rep","iFile","link"),
                        measure.vars = grep("^p.",names(dt.sim.MMtype1),value = TRUE),
                        value.name = "p.value", variable.name = "method")

dtLS.sim.MMtype1 <- dtL.sim.MMtype1[,.(n.rep=.N,
                                       type1=mean(p.value<=0.05,na.rm=TRUE)
                                       ), by = c("n","method","link")]

dtLS.sim.MMtype1[, robust := grepl("robust",method)]
dtLS.sim.MMtype1[, correction := gsub("p.|robust","",method)]

## ** Mixed model - bias
cat(" - mixed model (bias) \n")

dt.sim.MMestimate <- sinkDirectory(path.simulation.mixedModel,
                                   string.keep = "estimate", string.exclude = "(tempo)")
dt.sim.MMestimate[, bias.ML := estimate.truth - estimate.ML]
dt.sim.MMestimate[, bias.MLcorrected := estimate.truth - estimate.MLcorrected]

dtL.sim.MMestimate <- melt(dt.sim.MMestimate,
                       id.vars = c("n","rep","seed","iFile","type","name"),
                       measure.vars = c("bias.ML","bias.MLcorrected"),
                       value.name = "bias", variable.name = "corrected")
dtL.sim.MMestimate[,corrected := factor(corrected, levels = c("bias.ML","bias.MLcorrected"), labels = c("FALSE","TRUE"))]

dtLS.sim.MMestimate <- dtL.sim.MMestimate[,.(rep = .N,
                                     inf = quantile(bias,probs = 0.05),
                                     lower = quantile(bias,probs = 0.25),
                                     middle = quantile(bias,probs = 0.5),
                                     upper = quantile(bias,probs = 0.75),
                                     sup = quantile(bias,probs = 0.95),
                                     mean = mean(bias),
                                     sd = sd(bias)),
                                  by = c("n","corrected","type")]

## ** factor model - type 1 error
cat(" - factor model (type 1 error) \n")

dt.sim.factortype1 <- sinkDirectory(path.simulation.factorModel,
                                    string.keep = "type1error", string.exclude = "(tempo)")

## dt.sim.factortype1[link == link[1], .(nb=.N,pc=100*.N/20000,warning=sum(warning),niterMax=max(niter)), by = "n"]
##      n    nb      pc warning niterMax
## 1:  20 17852  89.260       0        9
## 2:  30 19572  97.860       0        7
## 3:  50 19969  99.845       0        6
## 4:  75 20000 100.000       0        5
## 5: 100 20000 100.000       0        5
## 6: 150 20000 100.000       0        4
## 7: 200 20000 100.000       0        4
## 8: 300 20000 100.000       0        4
## 9: 500 20000 100.000       0        4

dtL.sim.factortype1 <- melt(dt.sim.factortype1,
                            id.vars = c("n","rep","iFile","link"),
                            measure.vars = grep("^p.",names(dt.sim.factortype1),value = TRUE),
                            value.name = "p.value", variable.name = "method")

dtLS.sim.factortype1 <- dtL.sim.factortype1[,.(n.rep=.N,
                                       type1=mean(p.value<=0.05,na.rm=TRUE)
                                       ), by = c("n","method","link")]

dtLS.sim.factortype1[, robust := grepl("robust",method)]
dtLS.sim.factortype1[, correction := gsub("p.|robust","",method)]

## ** factor model - bias
cat(" - factor model (bias) \n")

dt.sim.factorestimate <- sinkDirectory(path.simulation.factorModel,
                                       string.keep = "estimate", string.exclude = "(tempo)")
dt.sim.factorestimate[, bias.ML := estimate.truth - estimate.ML]
dt.sim.factorestimate[, bias.MLcorrected := estimate.truth - estimate.MLcorrected]

dtL.sim.factorestimate <- melt(dt.sim.factorestimate,
                       id.vars = c("n","rep","seed","iFile","type","name"),
                       measure.vars = c("bias.ML","bias.MLcorrected"),
                       value.name = "bias", variable.name = "corrected")
dtL.sim.factorestimate[,corrected := factor(corrected, levels = c("bias.ML","bias.MLcorrected"), labels = c("FALSE","TRUE"))]

dtLS.sim.factorestimate <- dtL.sim.factorestimate[,.(rep = .N,
                                             inf = quantile(bias,probs = 0.05),
                                             lower = quantile(bias,probs = 0.25),
                                             middle = quantile(bias,probs = 0.5),
                                             upper = quantile(bias,probs = 0.75),
                                             sup = quantile(bias,probs = 0.95),
                                             mean = mean(bias),
                                             sd = sd(bias)),
                                          by = c("n","corrected","type")]

## ** lvm - type 1 error
cat(" - lvm (type 1 error) \n")

dt.sim.lvmtype1 <- sinkDirectory(path.simulation.lvm,
                                 string.keep = "type1error", string.exclude = "(tempo)")
## dt.sim.lvmtype1[link == link[1], .(nb=.N,pc=100*.N/20000,warning=sum(warning),niterMax=sum(niter==20),niter=mean(niter)), by = "n"]
##      n    nb      pc warning niterMax    niter
## 1:  20 16796  83.980       7       11 7.332579
## 2:  30 19482  97.410       0        0 6.370085
## 3:  50 19979  99.895       0        0 5.745833
## 4:  75 20000 100.000       0        0 5.000000
## 5: 100 20000 100.000       0        0 5.000000
## 6: 150 20000 100.000       0        0 4.022950
## 7: 200 20000 100.000       0        0 4.000000
## 8: 300 20000 100.000       0        0 4.000000
## 9: 500 20000 100.000       0        0 4.000000

dtL.sim.lvmtype1 <- melt(dt.sim.lvmtype1,
                         id.vars = c("n","rep","iFile","link"),
                         measure.vars = grep("^p.",names(dt.sim.lvmtype1),value = TRUE),
                         value.name = "p.value", variable.name = "method")

dtLS.sim.lvmtype1 <- dtL.sim.lvmtype1[,.(n.rep=.N,
                                       type1=mean(p.value<=0.05,na.rm=TRUE)
                                       ), by = c("n","method","link")]

dtLS.sim.lvmtype1[, robust := grepl("robust",method)]
dtLS.sim.lvmtype1[, correction := gsub("p.|robust","",method)]

## ** lvm - bias
cat(" - lvm (bias) \n")

dt.sim.lvmestimate <- sinkDirectory(path.simulation.lvm,
                               string.keep = "estimate", string.exclude = "(tempo)")
dt.sim.lvmestimate[, bias.ML := estimate.truth - estimate.ML]
dt.sim.lvmestimate[, bias.MLcorrected := estimate.truth - estimate.MLcorrected]

dtL.sim.lvmestimate <- melt(dt.sim.lvmestimate,
                       id.vars = c("n","rep","seed","iFile","type","name"),
                       measure.vars = c("bias.ML","bias.MLcorrected"),
                       value.name = "bias", variable.name = "corrected")
dtL.sim.lvmestimate[,corrected := factor(corrected, levels = c("bias.ML","bias.MLcorrected"), labels = c("FALSE","TRUE"))]

dtLS.sim.lvmestimate <- dtL.sim.lvmestimate[,.(rep = .N,
                                     inf = quantile(bias,probs = 0.05),
                                     lower = quantile(bias,probs = 0.25),
                                     middle = quantile(bias,probs = 0.5),
                                     upper = quantile(bias,probs = 0.75),
                                     sup = quantile(bias,probs = 0.95),
                                     mean = mean(bias),
                                     sd = sd(bias)),
                                  by = c("n","corrected","type")]


## ** export
saveRDS(dtLS.sim.MMtype1, file = file.path(path.results,"type1error-simulation-mixedModel.rds"))
saveRDS(dtLS.sim.factortype1, file = file.path(path.results,"type1error-simulation-factorModel.rds"))
saveRDS(dtLS.sim.lvmtype1, file = file.path(path.results,"type1error-simulation-lvmModel.rds"))

saveRDS(dtLS.sim.MMestimate, file = file.path(path.results,"bias-simulation-mixedModel.rds"))
saveRDS(dtLS.sim.factorestimate, file = file.path(path.results,"bias-simulation-factorModel.rds"))
saveRDS(dtLS.sim.lvmestimate, file = file.path(path.results,"bias-simulation-lvmModel.rds"))

saveRDS(unique(dt.sim.MMestimate$name), file = file.path(path.results,"param-simulation-mixedModel.rds"))
saveRDS(unique(dt.sim.factorestimate$name), file = file.path(path.results,"param-simulation-factorModel.rds"))
saveRDS(unique(dt.sim.lvmestimate$name), file = file.path(path.results,"param-simulation-lvmModel.rds"))

keep.col <- c("n","name","type","estimate.MLcorrected","se.MLcorrected","se.robustMLcorrected","df.MLcorrected","df.robustMLcorrected")
saveRDS(list(MM = dt.sim.MMestimate[,.SD, .SDcols = keep.col],
             factor = dt.sim.factorestimate[,.SD, .SDcols = keep.col],
             lvm = dt.sim.lvmestimate[,.SD, .SDcols = keep.col]),
        file = file.path(path.results,"dist-simulation.rds"))





## * Illustration
cat("Illustration: \n")
## ** Mixed model - type 1 error
cat(" - mixed model (type 1 error) \n")

dt.ill.MMtype1 <- sinkDirectory(path.illustration.mixedModel,
                                string.keep = "type1error", string.exclude = "(tempo)")

dtL.ill.MMtype1 <- melt(dt.ill.MMtype1,
                        id.vars = c("n","rep","iFile","link"),
                        measure.vars = grep("^p.",names(dt.ill.MMtype1),value = TRUE),
                        value.name = "p.value", variable.name = "method")
dtLS.ill.MMtype1 <- dtL.ill.MMtype1[,.(n.rep=.N,
                                       type1=mean(p.value<=0.05,na.rm=TRUE)
                                       ), by = c("n","method","link")]

dtLS.ill.MMtype1[, robust := grepl("robust",method)]
dtLS.ill.MMtype1[, correction := gsub("p.|robust","",method)]

## dtLS.ill.MMtype1[link == "global" & robust == FALSE]
## Ztest: 0.103
## SSC  : 0.057
## KR   : 0.045

## ** factor model - type 1 error
cat(" - factor model (type 1 error) \n")

dt.ill.factortype1 <- sinkDirectory(path.illustration.factorModel,
                                    string.keep = "type1error", string.exclude = "(tempo)")
## dt.ill.factortype1[, length(unique(rep)), by = iFile]

dtL.ill.factortype1 <- melt(dt.ill.factortype1,
                        id.vars = c("n","rep","iFile","link"),
                        measure.vars = grep("^p.",names(dt.ill.factortype1),value = TRUE),
                        value.name = "p.value", variable.name = "method")
dtLS.ill.factortype1 <- dtL.ill.factortype1[,.(n.rep=.N,
                                       type1=mean(p.value<=0.05,na.rm=TRUE)
                                       ), by = c("n","method","link")]

dtLS.ill.factortype1[, robust := grepl("robust",method)]
dtLS.ill.factortype1[, correction := gsub("p.|robust","",method)]

## dtLS.ill.factortype1[link %in% c("u~bdnf2","neo~httlpr2") & robust == TRUE & correction %in% c("Ztest","KR")]
## Ztest: 0.0695, 0.0700  ---> 0.0741, 0.0613
## KR   : 0.0555, 0.0602  ---> 0.0559, 0.0478

## ** lvm model - type 1 error
cat(" - lvm (type 1 error) \n")

dt.ill.lvmtype1 <- sinkDirectory(path.illustration.lvm,
                                 string.keep = "type1error", string.exclude = "(tempo)")

dt.ill.lvmtype1[,length(unique(rep)),by = "iFile"]

dtL.ill.lvmtype1 <- melt(dt.ill.lvmtype1,
                        id.vars = c("n","rep","iFile","link"),
                        measure.vars = grep("^p.",names(dt.ill.lvmtype1),value = TRUE),
                        value.name = "p.value", variable.name = "method")
dtLS.ill.lvmtype1 <- dtL.ill.lvmtype1[,.(n.rep=.N,
                                       type1=mean(p.value<=0.05,na.rm=TRUE)
                                       ), by = c("n","method","link")]

dtLS.ill.lvmtype1[, robust := grepl("robust",method)]
dtLS.ill.lvmtype1[, correction := gsub("p.|robust","",method)]

## dtLS.ill.lvmtype1[robust == FALSE & correction %in% c("Ztest","KR")]
## Ztest: 0.104, 0.112, 0.076 --> 0.0633, 0.084, 0.085
## KR   : 0.048, 0.058, 0.034 --> 0.02, 0.034, 0.037

## ** export
saveRDS(dtLS.ill.MMtype1, file = file.path(path.results,"type1error-illustration-mixedModel.rds"))
saveRDS(dtLS.ill.factortype1, file = file.path(path.results,"type1error-illustration-factorModel.rds"))
saveRDS(dtLS.ill.lvmtype1, file = file.path(path.results,"type1error-illustration-lvmModel.rds"))

## * IV
cat("IV: \n")
## ** type 1 error in small sample
cat(" - small sample \n")

dtL.sim.Comptype1 <- sinkDirectory(path.simulation.comparison,
                                   string.keep = "type1error", string.exclude = "(tempo)")
## dtL.sim.IVtype1[, link := name]
## dt.sim.IVtype1[warning==TRUE]
dtLS.sim.Comptype1 <- dtL.sim.Comptype1[,.(n.rep=.N,
                                           type1=mean(p.value<=0.05,na.rm=TRUE)
                                           ), by = c("n","estimator","link")]
dcast(dtLS.sim.Comptype1, formula = link+n ~ estimator, value.var = "type1")[n==20]
##           link  n        GLS      IV IVlavaan         ML WLS robustML
##  1:         Y2 20 0.06734401 0.10130  0.08075 0.08235824 NaN  0.09335
##  2:  Y2~Gene1Y 20 0.09734560 0.08825  0.08080 0.11495000 NaN  0.10805
##  3:     Y2~eta 20 0.11186765 0.16125  0.12825 0.08895890 NaN  0.09395
##  4:         Y3 20 0.05282196 0.07500  0.06955 0.06240624 NaN  0.06775
##  5:     Y3~eta 20 0.11356987 0.17995  0.14320 0.08560856 NaN  0.09275
##  6:         Y4 20 0.05729028 0.07720  0.07340 0.07260726 NaN  0.07385
##  7:     Y4~eta 20 0.08766424 0.10105  0.06655 0.07480748 NaN  0.08645
##  8:        eta 20 0.07479121 0.10805  0.08680 0.08825000 NaN  0.09390
##  9:    eta~Age 20 0.13543274 0.12555  0.08810 0.10681068 NaN  0.11440
## 10: eta~Gene2Y 20 0.08133411 0.09375  0.08495 0.09635964 NaN  0.08950

## dcast(dtLS.sim.Comptype1, formula = link+n ~ estimator, value.var = "type1")[n==30]
## dcast(dtLS.sim.Comptype1, formula = link+n ~ estimator, value.var = "type1")[n==50]
## dcast(dtLS.sim.Comptype1, formula = link+n ~ estimator, value.var = "type1")[n==100]
## ** type 1 error with student distributed residuals 
cat(" - student \n")

dtL.sim.IV1type1 <- sinkDirectory(path.simulation.IV1,
                                  string.keep = "type1error", string.exclude = "(tempo)")
100*mean(dtL.sim.IV1type1$shapiroMax<=1e-3)
## 99.61833
dtLS.sim.IV1type1 <- dtL.sim.IV1type1[,.(n.rep=.N,
                                         type1=mean(p.value<=0.05,na.rm=TRUE)
                                         ), by = c("n","estimator","link")]

dcast(dtLS.sim.IV1type1, formula = link+n ~ estimator, value.var = "type1")[n==20]
dcast(dtLS.sim.IV1type1, formula = link+n ~ estimator, value.var = "type1")[n==100]


## ** type 1 error with chi2 distributed residuals 
cat(" - chi2 \n")

dtL.sim.IV2type1 <- sinkDirectory(path.simulation.IV2,
                                  string.keep = "type1error", string.exclude = "(tempo)")
100*mean(dtL.sim.IV2type1$shapiroMax<=1e-3)
## 100
dtLS.sim.IV2type1 <- dtL.sim.IV2type1[,.(n.rep=.N,
                                         type1=mean(p.value<=0.05,na.rm=TRUE)
                                         ), by = c("n","estimator","link")]

dcast(dtLS.sim.IV2type1, formula = link+n ~ estimator, value.var = "type1")[n==20]


## ** type 1 error with correlated residuals 
cat(" - correlation \n")

dtL.sim.IV3type1 <- sinkDirectory(path.simulation.IV3,
                                  string.keep = "type1error", string.exclude = "(tempo)")
dtLS.sim.IV3type1 <- dtL.sim.IV3type1[,.(n.rep=.N,
                                         type1=mean(p.value<=0.05,na.rm=TRUE)
                                         ), by = c("n","estimator","link")]

dcast(dtLS.sim.IV3type1, formula = link+n ~ estimator, value.var = "type1")[n==1000]
##           link    n      IV IVlavaan      ML robustML robustMLC
##  1:         Y2 1000 0.05040  0.05035 0.05305  0.04980   0.04980
##  2:  Y2~Gene1Y 1000 0.05040  0.05155 0.06070  0.05015   0.05030
##  3:     Y2~eta 1000 0.04890  0.31255 0.04870  0.04880   0.04880
##  4:         Y3 1000 0.05105  0.04875 0.05105  0.05095   0.05095
##  5:     Y3~eta 1000 0.05095  0.30005 0.05135  0.05150   0.05130
##  6:         Y4 1000 0.04960  0.05105 0.04950  0.05015   0.05000
##  7:     Y4~eta 1000 0.04850  0.31305 0.04945  0.05105   0.05095
##  8:        eta 1000 0.31345  0.05045 0.32175  0.05160   0.05145
##  9:    eta~Age 1000 0.30115  0.04955 0.31785  0.04970   0.04960
## 10: eta~Gene2Y 1000 0.31335  0.04770 0.33235  0.05160   0.05145

dcast(dtLS.sim.IV3type1, formula = link+n ~ estimator, value.var = "type1")[n==20]

## ** export
saveRDS(dtLS.sim.Comptype1, file = file.path(path.results,"type1error-simulation-comparison.rds"))
saveRDS(dtLS.sim.IV1type1, file = file.path(path.results,"type1error-simulation-IV1.rds"))
saveRDS(dtLS.sim.IV3type1, file = file.path(path.results,"type1error-simulation-IV3.rds"))
