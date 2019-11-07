path.results <- "./Results"

## * Packages
library(data.table)

## * Import

## ** Simulation results
dttype1.ill.mm <- readRDS(file.path(path.results, "type1error-illustration-mixedModel.rds"))
dttype1.ill.factor <- readRDS(file.path(path.results, "type1error-illustration-factorModel.rds"))
dttype1.ill.lvm <- readRDS(file.path(path.results, "type1error-illustration-factorModel.rds"))

## ** Real data
source("ANALYSIS.R")

## * Application A
## ** Number of parameter of the LVM
length(coef.vitamin.ML)
## [1] 7

## ** Wald test mixed model
c(statistic = round(Ftest.vitamin.REML[["F value"]][2],2), 
  pvalue = round(Ftest.vitamin.REML[["Pr(>F)"]][2],4))
## statistic    pvalue 
##    4.2500    0.0102 

## ** Wald test LVM
c(statistic = round(Ftest.vitamin.ML[["statistic"]]/3, 2),
  pvalue = round(Ftest.vitamin.ML[["p.value"]], 5))
## statistic.chisq       pvalue.df 
##         5.02000         0.00176 

## ** Type 1 error by simulation
dttype1.ill.mm[link == "global" & method == "p.Ztest",type1]
## [1] 0.1022204

## * Application B
## ** Number of patients/observations
c(nrow = NROW(dt.bdnf), n.id = length(unique(dt.bdnf$cimbi.id)))
## nrow n.id 
##   73   68 

## ** Number of parameter of the LVM
length(coef.bdnf.ML)
## [1] 29

## ** Estimates
summary(e.bdnf2)$coef[bdnf.null,]
##                  Estimate  Robust SE   Naive SE      P-value
## u~bdnf2mx      0.07399285 0.02642470 0.02590882 5.108053e-03
## neo~httlpr2sx -0.07294591 0.01625799 0.01761417 7.231043e-06

## ** Type 1 error
dttype1.ill.factor[link %in% c("u~bdnf2","neo~httlpr2") & method %in% c("p.robustZtest"), .(link, type1 = type1)]
##           link      type1
## 1:     u~bdnf2 0.07411511
## 2: neo~httlpr2 0.06137384




## * Application C
## ** Number of patients/observations
c(nrow = NROW(dtG.memory), id = length(unique(dtG.memory$cimbi.id)))
## nrow   id 
##   24   24 

## ** Number of parameter of the LVM
length(coef(e.memory))
## [1] 48

## ** Estimates
summary(e.memory)$coef[memory.null,]

## ** Type 1 error
dttype1.ill.lvm[method == "p.Ztest",.(link,type1)][match(link,memory.null)]
