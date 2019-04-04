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
cat(table1)

table2 <- createTable(dt = dtLS.sim.factorbias, seqN = c(20,30,50,100), seqType = c("Sigma_var","Psi_var"), digit = 3, convert2latex = TRUE)
cat(table2)

table3 <- createTable(dt = dtLS.sim.lvmbias, seqN = c(20,30,50,100), seqType = c("Sigma_var","Psi_var","Sigma_cov"), digit = 3, convert2latex = TRUE)
cat(table3)
