## path <- "P:/Cluster/LVMproject/article-smallSampleInference"
## setwd(path)

library(data.table)
source("FCT.R") ## get function createTable
export <- TRUE

## * path
path.results <- "./Results"

## * load data
dtLS.sim.MMbias <- readRDS(file.path(path.results,"bias-simulation-mixedModel.rds"))

dtLS.sim.factorbias <- readRDS(file.path(path.results,"bias-simulation-factorModel.rds"))

dtLS.sim.lvmbias <- readRDS(file.path(path.results,"bias-simulation-lvmModel.rds"))

## * table 1
table1 <- createTable(dt = dtLS.sim.MMbias, seqN = c(20,30,50,100), seqType = c("Sigma_var","Psi_var"), digit = 3, convert2latex = TRUE)
dtLS.sim.MMbias[corrected == FALSE & n==20,.(type,mean,middle)]
table2 <- createTable(dt = dtLS.sim.factorbias, seqN = c(20,30,50,100), seqType = c("Sigma_var","Psi_var"), digit = 3, convert2latex = TRUE)
dtLS.sim.factorbias[corrected == FALSE & n==20,.(type,mean,middle)]
table3 <- createTable(dt = dtLS.sim.lvmbias, seqN = c(20,30,50,100), seqType = c("Sigma_var","Psi_var","Sigma_cov"), digit = 3, convert2latex = TRUE)
dtLS.sim.lvmbias[corrected == FALSE & n==20,.(type,mean,middle)]
