## path <- "P:/Cluster/LVMproject/article-smallSampleInference"
## setwd(path)
## source("BATCH_simulation-mixedModel.R")

rm(list = ls())
gc()

## * seed
iter_sim <- as.numeric(Sys.getenv("SGE_TASK_ID"))
n.iter_sim <- as.numeric(Sys.getenv("SGE_TASK_LAST"))
if(is.na(iter_sim)){iter_sim <- 89}
if(is.na(n.iter_sim)){n.iter_sim <- 200}
cat("iteration ",iter_sim," over ",n.iter_sim,"\n",sep="")

set.seed(1)
seqSeed <- sample(1:max(1e5,n.iter_sim),size=n.iter_sim,replace=FALSE)
iSeed <- seqSeed[iter_sim]
set.seed(iSeed)

cat("seed: ",iSeed,"\n")

## * path
path <- "."
path.res <- file.path(path,"Results","simulation-mixedModel")
if(dir.exists(path.res)==FALSE){
    dir.create(path.res)
}
path.output <- file.path(path,"output","simulation-mixedModel")
if(dir.exists(path.output)==FALSE){
    dir.create(path.output)
}

## * libraries
library(lava)
library(data.table)
devtools::load_all("lavaSearch2")

## * settings
seqN <- c(20,30,50,75,100,150,200,300,500)
n.boot <- 0# 1e3
n.rep <- 5e2

## * model

## ** generative model
m.generative <- lvm(c(Y1[mu1:sigma]~1*eta,
                      Y2[0:sigma]~1*eta,
                      Y3[mu3:sigma]~1*eta,
                      eta~beta1*Age))
latent(m.generative) <- ~eta
categorical(m.generative, labels = c("N","Y")) <- ~Gene1
 
## plot(m)

## ** investigator model
m.fit <- lvm(c(Y1[mu1:sigma]~1*eta,
               Y2[mu2:sigma]~1*eta,
               Y3[mu3:sigma]~1*eta,
               eta~beta1*Age + beta2*Gene1))
latent(m.fit) <- ~eta
categorical(m.fit, labels = c("N","Y")) <- ~Gene1
## lava::estimate(m.fit, lava::sim(m.generative, n = 1e5))

## ** true value of the coefficients
if(FALSE){ ## create true.coef
    true.coef <- round(coef(lava::estimate(m.fit, lava::sim(m.generative, n = 1e3))))
    Mtrue.coef <- cbind(names(true.coef)," = ",as.double(true.coef),",")
    Mtrue.coef[1,1] <- paste0("c(",Mtrue.coef[1,1])
    Mtrue.coef[NROW(Mtrue.coef),4] <- ")"
    Mtrue.coef <- cbind(apply(Mtrue.coef, 1, paste, collapse = ""))
    print(as.data.table(Mtrue.coef), quote = FALSE, row.names = FALSE)
}

true.coef <- c(eta = 0,
               Y2 = 0,
               Y3 = 0,
               "eta~Age" = 1,
               "eta~Gene1Y" = 0,
               "Y1~~Y1" = 1,
               "eta~~eta" = 1)

## * simulation
out <- calibrateType1(m.fit, true.coef = true.coef,
                      param = names(true.coef[1:5]),
                      null = true.coef[1:5],
                      n = seqN, n.rep = n.rep,
                      generative.object = m.generative,
                      dir.save = path.res, label.file = iter_sim,
                      bootstrap = FALSE, seed = NULL, trace = 2)

## * display
print(sessionInfo())
