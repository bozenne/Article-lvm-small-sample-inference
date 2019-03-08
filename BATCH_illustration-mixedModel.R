## path <- "P:/Cluster/LVMproject/article-smallSampleInference"
## setwd(path)
## source("BATCH_illustration-mixedModel.R")

rm(list = ls())
gc()

## * seed
iter_sim <- as.numeric(Sys.getenv("SGE_TASK_ID"))
n.iter_sim <- as.numeric(Sys.getenv("SGE_TASK_LAST"))
if(is.na(iter_sim)){iter_sim <- 2}
if(is.na(n.iter_sim)){n.iter_sim <- 180}
cat("iteration ",iter_sim," over ",n.iter_sim,"\n",sep="")

set.seed(2)
seqSeed <- sample(1:max(1e5,n.iter_sim),size=n.iter_sim,replace=FALSE)
iSeed <- seqSeed[iter_sim]
set.seed(iSeed)

cat("seed: ",iSeed,"\n")

## * path
path <- "."
path.res <- file.path(path,"Results","illustration-mixedModel")
if(dir.exists(path.res)==FALSE){
    dir.create(path.res)
}
path.output <- file.path(path,"output","illustration-mixedModel")
if(dir.exists(path.output)==FALSE){
    dir.create(path.output)
}

## * libraries
library(lava)
library(data.table)
## devtools::load_all("lavaSearch2")
library(lavaSearch2)

## * settings
seqN <- c(10)
n.boot <- 0# 1e3
n.rep <- 5e2


## * model

## ** generative model
m.generative <- lvm(w1[mu1:sigma] ~ 1*eta,
                    w3[mu3:sigma] ~ 1*eta,
                    w4[mu4:sigma] ~ 1*eta,
                    w5[mu5:sigma] ~ grp + 1*eta,
                    w6[mu6:sigma] ~ grp + 1*eta,
                    w7[mu7:sigma] ~ grp + 1*eta,
                    eta ~ 0, u ~ 1)
latent(m.generative) <- ~eta+u
transform(m.generative,grp~u) <- function(x){
  n <- NROW(x)
  return(c(rep(1,n/2),rep(0,n/2)))
}
 
## plot(m)

generative.coef <- c("mu1" = -1.20756,
                     "mu3" = -0.328751,
                     "mu4" = 0.253379,
                     "mu5" = 0.246786,
                     "mu6" = 0.006236,
                     "mu7" = 0.413567,
                     "w5~grp" = 0,
                     "w6~grp" = 0,
                     "w7~grp" = 0,
                     "sigma" = 0.148296,
                     "eta~~eta" = 0.348873)

## ** fit model
m.fit <- lvm(w1[mu1:sigma] ~ 1*eta,
             w3[mu3:sigma] ~ 1*eta,
             w4[mu4:sigma] ~ 1*eta,
             w5[mu5:sigma] ~ grp + 1*eta,
             w6[mu6:sigma] ~ grp + 1*eta,
             w7[mu7:sigma] ~ grp + 1*eta,
             eta ~ 0)
latent(m.fit) <- ~eta      
## lava::estimate(m.fit, lava::sim(m.generative, n = 1e5, p = generative.coef))
## vcov(lava::estimate(m.fit, lava::sim(m.generative, n = 10, p = generative.coef)))

## ** true value of the coefficients
if(FALSE){ ## create true.coef

    ## perform analysis under the null
    dtL.vitamin <- fread(file.path(path.data, "vitamin.txt"), header = TRUE)

    dtL.vitamin[, animal := as.factor(animal)]
    dtL.vitamin[, grp := factor(grp, levels = 1:2, labels = c("C","T"))]
    dtL.vitamin[, week.factor := paste0("w",as.factor(week))]
    dtL.vitamin[, weight0 := scale(weight)]

    dtW.vitamin <- dcast(dtL.vitamin, value.var = "weight0",
                         formula = grp + animal ~ week.factor)
    
    m0.vitamin <- lvm(w1[mu1:sigma] ~ 1*eta,
                      w3[mu3:sigma] ~ 1*eta,
                      w4[mu4:sigma] ~ 1*eta,
                      w5[mu5:sigma] ~ grp + 1*eta,
                      w6[mu6:sigma] ~ grp + 1*eta,
                      w7[mu7:sigma] ~ grp + 1*eta,
                      eta ~ 0)
    latent(m0.vitamin) <- ~eta

    vitamin0.ML <- estimate(m0.vitamin, data = dtW0.vitamin)

    ## extract coef and display
    true.coef <- round(coef(vitamin0.ML), 6)
    true.coef[grep("grpT",names(true.coef))] <- 0 ## under the null
    Mtrue.coef <- cbind(names(true.coef)," = ",as.double(true.coef),",")
    Mtrue.coef[1,1] <- paste0("c(",Mtrue.coef[1,1])
    Mtrue.coef[NROW(Mtrue.coef),4] <- ")"
    Mtrue.coef <- cbind(apply(Mtrue.coef, 1, paste, collapse = ""))
    print(as.data.table(Mtrue.coef), quote = FALSE, row.names = FALSE)

}
true.coef <- c("w1" = -1.20756,
               "w3" = -0.328751,
               "w4" = 0.253379,
               "w5" = 0.246786,
               "w6" = 0.006236,
               "w7" = 0.413567,
               "w5~grp" = 0,
               "w6~grp" = 0,
               "w7~grp" = 0,
               "w1~~w1" = 0.148296,
               "eta~~eta" = 0.348873)

## * simulation
## n.rep <- 10
out <- calibrateType1(m.fit, true.coef = true.coef,
                      param = c("w5~grp","w6~grp","w7~grp"), F.test = TRUE,
                      n = seqN, n.rep = n.rep,
                      generative.object = m.generative, generative.coef = generative.coef,
                      dir.save = path.res, label.file = iter_sim,
                      bootstrap = FALSE, seed = NULL, trace = 2)

## out$p.value

## * display
print(sessionInfo())
summary(out)

