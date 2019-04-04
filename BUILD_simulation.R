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
cat("Simulation: \n")

## ** Mixed model - type 1 error
cat(" - mixed model (type 1 error) \n")

dt.sim.MMtype1 <- sinkDirectory(path.simulation.mixedModel,
                                string.keep = "type1error", string.exclude = "(tempo)")

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
dt.sim.MMestimate[, bias.ML := estimate.ML - estimate.truth]
dt.sim.MMestimate[, bias.MLcorrected := estimate.MLcorrected - estimate.truth]

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
dt.sim.factorestimate[, bias.ML := estimate.ML - estimate.truth]
dt.sim.factorestimate[, bias.MLcorrected := estimate.MLcorrected - estimate.truth]

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
dt.sim.lvmestimate[, bias.ML := estimate.ML - estimate.truth]
dt.sim.lvmestimate[, bias.MLcorrected := estimate.MLcorrected - estimate.truth]

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

## ** factor model - type 1 error
cat(" - factor model (type 1 error) \n")

dt.ill.factortype1 <- sinkDirectory(path.illustration.factorModel,
                                    string.keep = "type1error", string.exclude = "(tempo)")

dtL.ill.factortype1 <- melt(dt.ill.factortype1,
                        id.vars = c("n","rep","iFile","link"),
                        measure.vars = grep("^p.",names(dt.ill.factortype1),value = TRUE),
                        value.name = "p.value", variable.name = "method")
dtLS.ill.factortype1 <- dtL.ill.factortype1[,.(n.rep=.N,
                                       type1=mean(p.value<=0.05,na.rm=TRUE)
                                       ), by = c("n","method","link")]

dtLS.ill.factortype1[, robust := grepl("robust",method)]
dtLS.ill.factortype1[, correction := gsub("p.|robust","",method)]

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

dtLS.sim.Comptype1 <- dtL.sim.Comptype1[,.(n.rep=.N,
                                           type1=mean(p.value<=0.05,na.rm=TRUE)
                                           ), by = c("n","estimator","link")]
dcast(dtLS.sim.Comptype1, formula = link+n ~ estimator, value.var = "type1")[n==20]

## ** type 1 error with student distributed residuals 
cat(" - student \n")

dtL.sim.IV1type1 <- sinkDirectory(path.simulation.IV1,
                                  string.keep = "type1error", string.exclude = "(tempo)")
100*mean(dtL.sim.IV1type1$shapiroMax<=1e-3)

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

dcast(dtLS.sim.IV3type1, formula = link+n ~ estimator, value.var = "type1")[n==20]
dcast(dtLS.sim.IV3type1, formula = link+n ~ estimator, value.var = "type1")[n==20]

## ** export
saveRDS(dtLS.sim.Comptype1, file = file.path(path.results,"type1error-simulation-comparison.rds"))
saveRDS(dtLS.sim.IV1type1, file = file.path(path.results,"type1error-simulation-IV1.rds"))
saveRDS(dtLS.sim.IV2type1, file = file.path(path.results,"type1error-simulation-IV2.rds"))
saveRDS(dtLS.sim.IV3type1, file = file.path(path.results,"type1error-simulation-IV3.rds"))
