path.results <- "./Results"

## * Packages
library(data.table)

## * Import

## ** Simulation results
dttype1.sim.mm <- readRDS(file.path(path.results, "type1error-simulation-mixedModel.rds"))
dtbias.sim.mm <- readRDS(file.path(path.results, "bias-simulation-mixedModel.rds"))
param.sim.mm <- readRDS(file.path(path.results, "param-simulation-mixedModel.rds"))

dttype1.sim.factor <- readRDS(file.path(path.results, "type1error-simulation-factorModel.rds"))
dtbias.sim.factor <- readRDS(file.path(path.results, "bias-simulation-factorModel.rds"))
param.sim.factor <- readRDS(file.path(path.results, "param-simulation-factorModel.rds"))

dttype1.sim.lvm <- readRDS(file.path(path.results, "type1error-simulation-lvmModel.rds"))
dtbias.sim.lvm <- readRDS(file.path(path.results, "bias-simulation-lvmModel.rds"))
param.sim.lvm <- readRDS(file.path(path.results, "param-simulation-lvmModel.rds"))


## * Application A
## nu2 is Y2
## gamma1 is eta Gene1Y

## ** Convergence
dttype1.sim.mm[,100*min(n.rep)/20000, by = "n"]
##      n      V1
## 1:  20  99.660
## 2:  30  99.975
## 3:  50 100.000
## 4:  75 100.000
## 5: 100 100.000
## 6: 150 100.000
## 7: 200 100.000
## 8: 300 100.000
## 9: 500 100.000

## ** Number of parameters
length(param.sim.mm)
## [1] 7

## ** Inflation of the type 1 error without correction
dttype1.sim.mm[n==20 & method %in% c("p.Ztest"), .(link,inflation = type1-0.05)]
##          link  inflation
## 1:        eta 0.02736303
## 2:         Y2 0.01456954
## 3:         Y3 0.01306442
## 4:    eta~Age 0.03965483
## 5: eta~Gene1Y 0.03674493

## type 1 error after correction
dttype1.sim.mm[n==20 & method %in% c("p.KR")]
##     n method       link n.rep      type1 robust correction
## 1: 20   p.KR        eta 19932 0.05057194  FALSE         KR
## 2: 20   p.KR         Y2 19932 0.05037126  FALSE         KR
## 3: 20   p.KR         Y3 19932 0.04866546  FALSE         KR
## 4: 20   p.KR    eta~Age 19932 0.05338150  FALSE         KR
## 5: 20   p.KR eta~Gene1Y 19932 0.05017058  FALSE         KR

## * Application B
## ** Convergence
dttype1.sim.factor[,100*min(n.rep)/20000, by = "n"]
##      n      V1
## 1:  20  89.260
## 2:  30  97.860
## 3:  50  99.845
## 4:  75 100.000
## 5: 100 100.000
## 6: 150 100.000
## 7: 200 100.000
## 8: 300 100.000
## 9: 500 100.000

## ** Number of parameters
length(param.sim.factor)
## [1] 15

## ** Bias of ML
dtbias.sim.factor[n==20 & type == "Sigma_var" & corrected == FALSE,mean]
dtbias.sim.factor[n==20 & type == "Psi_var" & corrected == FALSE,mean]

## ** Bias of corrected ML
dtbias.sim.factor[n==20 & type == "Sigma_var" & corrected == TRUE,mean]
dtbias.sim.factor[n==20 & type == "Psi_var" & corrected == TRUE,mean]

## ** Inflation of the type 1 error without correction
dttype1.sim.factor[n == 20 & link %in% c("Y2","Y4~eta","eta~Gene1Y","Y1~Gene2Y") & method %in% c("p.robustZtest"), .(link,type1 = type1, inflation = type1-0.05)]
##          link      type1  inflation
## 1:         Y2 0.08811338 0.03811338
## 2:  Y1~Gene2Y 0.11656957 0.06656957
## 3: eta~Gene1Y 0.10889536 0.05889536
## 4:     Y4~eta 0.11259243 0.06259243

## ** Type 1 error after correction
dttype1.sim.factor[n == 20 & link %in% c("Y2","Y4~eta","eta~Gene1Y","Y1~Gene2Y") & method %in% c("p.robustKR"), .(link,type1)]
##          link      type1
## 1:         Y2 0.04105983
## 2:  Y1~Gene2Y 0.07410934
## 3: eta~Gene1Y 0.05640825
## 4:     Y4~eta 0.03601837

## * Application C
## ** Number of parameters in the model
length(param.sim.lvm)
## [1] 36

## ** Type 1 error 
keep.link <- c("Y2", "Y4~eta1","Y1~Gene2Y", "eta1~Gene1Y","eta1~eta2","Y1~~Y2")
cbind(dttype1.sim.lvm[n==20 & link %in% keep.link & method %in% c("p.Ztest"), .(link,type1ML = type1)],
      dttype1.sim.lvm[n==20 & link %in% keep.link & method %in% c("p.SSC"),.(type1MLssc = type1)],
      dttype1.sim.lvm[n==20 & link %in% keep.link & method %in% c("p.KR"),.(type1MLc = type1)])
##           link    type1ML type1MLssc   type1MLc
## 1:          Y2 0.07662539 0.06114915 0.03170063
## 2:   Y1~Gene2Y 0.12258871 0.09913075 0.07752769
## 3:     Y4~eta1 0.07489879 0.06823459 0.02269343
## 4:   eta1~eta2 0.10014289 0.07325789 0.01293052
## 5: eta1~Gene1Y 0.12080257 0.09074670 0.06123786
## 6:      Y1~~Y2 0.09371279 0.06085713 0.02390320
