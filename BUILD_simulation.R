## path <- "P:/Cluster/LVMproject/article-smallSampleInference"
## setwd(path)

library(data.table)
source("FCT.R") ## get function sinkDirectory

## * path
path.results <- "./Results"
path.simulation.mixedModel <- "./Results/simulation-mixedModel"
path.simulation.factorModel <- "./Results/simulation-factorModel"
path.simulation.lvm <- "./Results/simulation-lvm"
path.simulation.IV <- "./Results/IV"
path.illustration.mixedModel <- "./Results/illustration-mixedModel"
path.illustration.factorModel <- "./Results/illustration-factorModel"
path.illustration.lvm <- "./Results/illustration-lvm"


digit.table <- 4
vec.greek <- c("alpha" = expression(alpha),
               "Gamma" = expression(Gamma),
               "Lambda" = expression(lambda),
               "K" = "K",
               "nu" = expression(nu),
               "Psi_var" = expression(Sigma[zeta]),
               "Psi_cov" = expression(Sigma[list(zeta,tilde(zeta))]),
               "Sigma_var" = expression(Sigma[epsilon]),
               "Sigma_cov" = expression(Sigma[list(epsilon,tilde(epsilon))])
               )

## * Simulation
## ** Mixed model - type 1 error
dt.sim.MMtype1 <- sinkDirectory(path.simulation.mixedModel,
                                string.keep = "type1error", string.exclude = "(tempo)")
## dt.sim.MMtype1[warning==TRUE]
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
dt.sim.factortype1 <- sinkDirectory(path.simulation.factorModel,
                                    string.keep = "type1error", string.exclude = "(tempo)")
## dt.sim.factortype1[warning==TRUE]
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
dt.sim.lvmtype1 <- sinkDirectory(path.simulation.lvm,
                                 string.keep = "type1error", string.exclude = "(tempo)")
## dt.sim.lvmtype1[warning==TRUE]
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



## ** fit distribution
if(FALSE){
    ## iDT <- dtL.sim.factorestimate[name == "eta1~eta2" & n == 20,]
    ## iDT <- dt.sim.factorestimate[name == "Y4~eta" & n == 20,]
    iDT <- dt.sim.factorestimate[name == "Y1~Gene2Y" & n == 20,]
    gg <- ggplot(iDT,aes(x = estimate.ML))
    gg <- gg + geom_histogram(aes(y = ..density..), binwidth = 0.05)
    gg <- gg + stat_function(fun = dnorm, n = 100, args = list(mean = mean(iDT$estimate.ML), sd = sd(iDT$estimate.ML)), color = "red")
    gg <- gg + facet_wrap(~n) + coord_cartesian(xlim = c(-1,1))
    fitStudent <- function(x,
                           mu = mean(x), se, df,
                           hist = TRUE, breaks = 100, xlim = c(-7,7), n.points = 100){ ##x <- iDT$bias
        var.x <- var(x)
        mean.x <- mean(x)
        out <- MASS::fitdistr(x = x, densfun = function(x, m, s, df) dt((x-m)/s, df)/s,
                              start = list(m = mean.x, s = var.x, df = 5), lower = c(-Inf, 0, 1))
        if(hist){
            x2 <- (x-out$estimate["m"])/out$estimate["s"]
            hist(x2, breaks = breaks, freq = FALSE, xlim = xlim)
            seqX <- seq(min(x2),max(x2),length.out=n.points) 
            points(seqX, y = dt(x = seqX, df = 3), type = "l", col = "red")
        }
        return(out)
    }
    fitStudent(iDT[,estimate.ML], se = mean(iDT$se.MLcorrected), df = mean(iDT$df.MLcorrected))

    iDT[,.(se = mean(se.ML), seSSC = mean(se.MLcorrected), df = mean(df.ML), df.MLcorrected = mean(df.MLcorrected))]
    
}

## * Illustration
## ** Mixed model - type1
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

## ** factor model - type1
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

## ** lvm model - type1
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
dt.sim.IVtype1 <- sinkDirectory(ath.simulation.IV,
                                string.keep = "type1error", string.exclude = "(tempo)")
## dt.sim.IVtype1[warning==TRUE]
dtL.sim.IVtype1 <- melt(dt.sim.IVtype1,
                        id.vars = c("n","rep","iFile","link"),
                        measure.vars = grep("^p.",names(dt.sim.IVtype1),value = TRUE),
                        value.name = "p.value", variable.name = "method")

dtLS.sim.IVtype1 <- dtL.sim.IVtype1[,.(n.rep=.N,
                                       type1=mean(p.value<=0.05,na.rm=TRUE)
                                       ), by = c("n","method","link")]

dtLS.sim.IVtype1[, robust := grepl("robust",method)]
dtLS.sim.IVtype1[, correction := gsub("p.|robust","",method)]

saveRDS(dtLS.sim.IVtype1, file = file.path(path.results,"type1error-simulation-IV.rds"))
