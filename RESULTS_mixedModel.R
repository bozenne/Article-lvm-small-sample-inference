## path <- "P:/Cluster/LVMproject/article-smallSampleInference"
## setwd(path)
path.results <- "./Results"
path.data <- "./data/"

## * packages
devtools::load_all("lavaSearch2")
library(lme4)
library(lmerTest)
library(data.table)
library(ggplot2)

## * Import

## ** simulation results
dttype1.sim.mm <- readRDS(file.path(path.results, "type1error-simulation-mixedModel.rds"))
dtbias.sim.mm <- readRDS(file.path(path.results, "bias-simulation-mixedModel.rds"))
param.sim.mm <- readRDS(file.path(path.results, "param-simulation-mixedModel.rds"))
dttype1.ill.mm <- readRDS(file.path(path.results, "type1error-illustration-mixedModel.rds"))

## ** real data
dtL.vitamin <- fread(file.path(path.data, "vitamin.txt"), header = TRUE)

## * Analysis of the dataset

## ** data processing
dtL.vitamin[, animal := as.factor(animal)]
dtL.vitamin[, grp := factor(grp, levels = 1:2, labels = c("C","T"))]
dtL.vitamin[, treat := factor(grp, levels = c("C","T"), labels = c("No","Yes"))]
dtL.vitamin[week<5, treat := "No"] ## no treatment before week five
dtL.vitamin[, week.factor := paste0("w",as.factor(week))]
dtL.vitamin[, weight0 := scale(weight)]
dtL.vitamin[,interaction := paste0(treat,week.factor)]
dtL.vitamin[treat=="No" | week <5, interaction := "base"]

dtW.vitamin <- dcast(dtL.vitamin, value.var = "weight0",
                     formula = grp + animal ~ week.factor)

## ** random intercept model
vitamin.REML <- lmer(weight0 ~ -1 + week.factor + interaction + (1|animal),
                     data = dtL.vitamin, REML = TRUE)
## summary(vitamin.REML)

coef.vitamin.REML <- c(fixef(vitamin.REML), 
                       sigma2 = sigma(vitamin.REML)^2, 
                       tau = as.double(summary(vitamin.REML)$varcor))
## coef.vitamin.REML

## ** latent variable model 
m.vitamin <- lvm(w1[mu1:sigma] ~ 1*eta,
                 w3[mu3:sigma] ~ 1*eta,
                 w4[mu4:sigma] ~ 1*eta,
                 w5[mu5:sigma] ~ grp + 1*eta,
                 w6[mu6:sigma] ~ grp + 1*eta,
                 w7[mu7:sigma] ~ grp + 1*eta,
                 eta ~ 0)
latent(m.vitamin) <- ~eta

vitamin.ML <- estimate(m.vitamin, data = dtW.vitamin)
## summary(vitamin.ML)

coef.vitamin.ML <- coef(vitamin.ML)
## coef.vitamin.ML

## * Simulation study (section 7.1)

## convergence
dttype1.sim.mm[,100*min(n.rep)/20000, by = "n"]

## number of parameters
length(param.sim.mm)

## inflation of the type 1 error without correction
dttype1.sim.mm[n==20 & method %in% c("p.Ztest"), .(link,inflation = type1-0.05)]
## nu2 is Y2
## gamma1 is eta Gene1Y

## type 1 error after correction
dttype1.sim.mm[n==20 & method %in% c("p.KR")]

## * Illustration (section 8.1)

## REML estimates
Ftest.vitamin.REML <- anova(vitamin.REML)
cat("\n")
c(statistic = round(Ftest.vitamin.REML[["F value"]][2],2), 
  pvalue = round(Ftest.vitamin.REML[["Pr(>F)"]][2],4))

## REML vs. ML
rdiff <- 100*(coef.vitamin.ML - coef.vitamin.REML)/abs(coef.vitamin.REML)
cbind("rdiff (%)" = round(rdiff,2))

## ML estimates
Ftest.vitamin.ML <- lava::compare(vitamin.ML, par = c("w5~grpT","w6~grpT","w7~grpT")) ## 15.072/3
cat("\n")
c(statistic = round(Ftest.vitamin.ML[["statistic"]]/3, 2),
  pvalue = round(Ftest.vitamin.ML[["p.value"]], 5))

## type 1 error found in the simulation
dttype1.ill.mm[link == "global" & method == "p.Ztest",type1]


## corrected ML estimates
sCorrect(vitamin.ML) <- TRUE
coef.vitamin.MLc <- vitamin.ML$sCorrect$param
Mcompare <- cbind("ML-correct" = coef.vitamin.MLc,
                  "ML" = coef.vitamin.ML, 
                  "rdiff ML (%)" = 100*(coef.vitamin.MLc-coef.vitamin.ML)/abs(coef.vitamin.ML),
                  "REML" = coef.vitamin.REML,
                  "rdiff REML (%)" = 100*(coef.vitamin.MLc-coef.vitamin.REML)/abs(coef.vitamin.REML)
                  )
Mcompare[c("w1~~w1","eta~~eta"),]

## REML vs. corrected ML
Ftest.vitamin.MLc <- compare2(vitamin.ML, par = c("w5~grpT","w6~grpT","w7~grpT")) ## 15.072/3
cat("\n")
c(statistic = round(Ftest.vitamin.MLc[["statistic"]], 2),
  pvalue = round(Ftest.vitamin.MLc[["p.value"]], 5))

## type 1 error found in the simulation
dttype1.ill.mm[link == "global" & method %in% c("p.SSC","p.KR"),.(method,link,type1)]


## ** Table 2

row.ML <- c("residual variance" = as.double(coef.vitamin.ML["w1~~w1"]),
            "variance random intercept" = as.double(coef.vitamin.ML["eta~~eta"]),
            "statistic" = as.double(Ftest.vitamin.ML[["statistic"]])/3,
            "degrees of freedom" = Inf,
            "p-value" = as.double(Ftest.vitamin.ML[["p.value"]]))
row.correctedML <- c("residual variance" = as.double(coef.vitamin.MLc["w1~~w1"]),
                     "variance random intercept" = as.double(coef.vitamin.MLc["eta~~eta"]),
                     "statistic" = as.double(Ftest.vitamin.MLc[["statistic"]]),
                     "degrees of freedom" = as.double(Ftest.vitamin.REML[["DenDF"]][2]),
              "p-value" = as.double(Ftest.vitamin.REML[["Pr(>F)"]][2]))
M <- cbind("ML" = row.ML, 
           "ML with correction" = row.correctedML,
           "REML" = row.REML)
M.txt <- rbind(format(M[1:3,], digits = 3, nsmall = 3),
               format(M[4,,drop=FALSE], digits = 2, nsmall = 2),
               format(M[5,,drop=FALSE], digits = 2, nsmall = 2))
M.txt[is.infinite(M)] <- "$\\infty$"

print(xtable::xtable(M.txt,
                     label = "tab:mixed"),
      include.colnames = TRUE,
      include.rownames = TRUE,
      sanitize.text.function = function(x){x},
      booktabs = TRUE)

## ** figure (Appendix A.1)

dtL.vitamin[, group := factor(grp, levels = c("C","T"), labels = c("Control","Treatment"))]
gg.spaguetti <- ggplot(dtL.vitamin, aes(x = as.factor(week), y = weight0, group = animal,
                                        color = animal))
gg.spaguetti <- gg.spaguetti + geom_line(size = 2) + geom_point(size = 3)
gg.spaguetti <- gg.spaguetti + facet_grid(~group, labeller = label_both)
gg.spaguetti <- gg.spaguetti + xlab("week") + ylab("weight")
gg.spaguetti <- gg.spaguetti + theme(text = element_text(size=20))
