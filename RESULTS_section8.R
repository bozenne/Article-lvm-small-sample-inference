path.results <- "./Results"

## * Packages
library(data.table)

## * Import

## ** Simulation results
dttype1.ill.mm <- readRDS(file.path(path.results, "type1error-illustration-mixedModel.rds"))
dttype1.ill.factor <- readRDS(file.path(path.results, "type1error-illustration-factorModel.rds"))
dttype1.ill.lvm <- readRDS(file.path(path.results, "type1error-illustration-lvmModel.rds"))

## ** Real data
source("ANALYSIS.R")

## * Application A
## ** Corrected ML vs REML and corrected ML vs. ML
cbind("ML-correct" = coef.vitamin.MLc,
      "ML" = coef.vitamin.ML, 
      "rdiff ML (%)" = 100*(coef.vitamin.MLc-coef.vitamin.ML)/abs(coef.vitamin.ML),
      "REML" = coef.vitamin.REML,
      "rdiff REML (%)" = 100*(coef.vitamin.MLc-coef.vitamin.REML)/abs(coef.vitamin.REML)
      )[c("w1~~w1","eta~~eta"),]
##          ML-correct        ML rdiff ML (%)      REML rdiff REML (%)
## w1~~w1    0.1768755 0.1482961     19.27189 0.1762872     0.33369954
## eta~~eta  0.3888395 0.3488731     11.45584 0.3885398     0.07711967

## ** Corrected stat and p-value
c(statistic = round(Ftest.vitamin.MLc[["statistic"]], 2),
  pvalue = round(Ftest.vitamin.MLc[["p.value"]], 5))

## ** Type 1 error found in the simulation
dttype1.ill.mm[link == "global" & method %in% c("p.KR"),.(method,link,type1)]
##    method   link      type1
## 1:   p.KR global 0.04490898

## * Application B
## ** Corrected ML vs. ML
cbind("ML-correct" = coef.bdnf.MLc,
      "ML" = coef.bdnf.ML, 
      "rdiff ML (%)" = 100*(coef.bdnf.MLc-coef.bdnf.ML)/abs(coef.bdnf.ML)
      )[grep("~~",names(coef.bdnf.MLc)),]
##            ML-correct          ML rdiff ML (%)
## neo~~neo 0.0004038557 0.000325987    23.887049
## u~~u     0.0117522161 0.010594862    10.923726
## cau~~cau 0.0117670745 0.011259093     4.511744
## put~~put 0.0049439930 0.004755134     3.971691
## hip~~hip 0.0120897817 0.011644897     3.820427
## amy~~amy 0.0249588302 0.024045687     3.797536
## cau~~put 0.0041899347 0.004030385     3.958668
## hip~~amy 0.0112839263 0.010866509     3.841319

## ** Corrected stat and p-values
summary2(e.bdnf2, robust = TRUE)$coef[bdnf.null,]
##                  Estimate  robust SE   t-value      P-value       df
## u~bdnf2mx      0.07399285 0.02749638  2.691004 9.034888e-03 65.37332
## neo~httlpr2sx -0.07294591 0.01678952 -4.344729 4.802431e-05 67.50082

## ** Compare p-values
cbind(pvalue = summary(e.bdnf2, robust = TRUE)$coef[bdnf.null,"P-value"],
      pvalueC = summary2(e.bdnf2, robust = TRUE)$coef[bdnf.null,"P-value"],
      pc.increase =  100*(summary2(e.bdnf2, robust = TRUE)$coef[bdnf.null,"P-value"] / summary(e.bdnf2, robust = TRUE)$coef[bdnf.null,"P-value"]-1)
      )
##                     pvalue      pvalueC pc.increase
## u~bdnf2mx     5.108053e-03 9.034888e-03    76.87536
## neo~httlpr2sx 7.231043e-06 4.802431e-05   564.14085

## ** Type 1 error found in the simulation
dttype1.ill.factor[link %in% c("u~bdnf2","neo~httlpr2") & method %in% c("p.robustKR"), .(link,type1)]
##           link      type1
## 1:     u~bdnf2 0.05586332
## 2: neo~httlpr2 0.04775789
## * Application C
## ** Correct ML vs. ML
coef.type <- as.data.table(coefType(e.memory, as.lava = FALSE))
vec.type <- coef.type[!is.na(lava),setNames(detail,name)]

Mcompare <- data.table("link" = names(coef.memory.MLc),
                       "type" = vec.type[names(coef.memory.MLc)],
                       "ML-correct" = coef.memory.MLc,
                       "ML" = coef.memory.ML, 
                       "rdiff ML (%)" = 100*(coef.memory.MLc-coef.memory.ML)/abs(coef.memory.ML)
                       )
Mcompare[,.(min = min(`rdiff ML (%)`), max = max(`rdiff ML (%)`)),by=type]
##         type      min       max
## 1:     alpha 0.000000  0.000000
## 2:        nu 0.000000  0.000000
## 3:         K 0.000000  0.000000
## 4:     Gamma 0.000000  0.000000
## 5:    Lambda 0.000000  0.000000
## 6:         B 0.000000  0.000000
## 7: Sigma_var 4.779645 21.199375
## 8:   Psi_var 9.090919 10.150635
## 9: Sigma_cov 4.832805  4.832805
## ** Corrected stat and p-values
summary2(e.memory)$coef[memory.null,"P-value",drop=FALSE]
##             P-value
## m.pos~u 0.004528687
## m.neu~u 0.019622746
## m.neg~u 0.124485463
## ** Compare p-values
cbind(pvalue = summary(e.memory)$coef[memory.null,"P-value"],
      pvalueC = summary2(e.memory)$coef[memory.null,"P-value"],
      pc.increase =  100*(summary2(e.memory)$coef[memory.null,"P-value"] / summary(e.memory)$coef[memory.null,"P-value"]-1)
      )
##               pvalue     pvalueC pc.increase
## m.pos~u 0.0005447762 0.004528687   731.29297
## m.neu~u 0.0042110757 0.019622746   365.97941
## m.neg~u 0.0720243779 0.124485463    72.83796
## ** Type 1 error found in the simulation
dttype1.ill.lvm[method == "p.KR",.(link,type1MLc = type1)][match(link,memory.null)]
##       link   type1MLc
## 1: m.pos~u 0.02023424
## 2: m.neu~u 0.03709611
## 3: m.neg~u 0.03422504
