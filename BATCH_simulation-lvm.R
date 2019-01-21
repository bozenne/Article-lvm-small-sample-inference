## path <- "P:/Cluster/LVMproject/article-smallSampleInference"
## setwd(path)
## source("BATCH_simulation-lvm.R")

rm(list = ls())
gc()

## * seed
iter_sim <- as.numeric(Sys.getenv("SGE_TASK_ID"))
n.iter_sim <- as.numeric(Sys.getenv("SGE_TASK_LAST"))
if(is.na(iter_sim)){iter_sim <- 15}
if(is.na(n.iter_sim)){n.iter_sim <- 40}
cat("iteration ",iter_sim," over ",n.iter_sim,"\n",sep="")

set.seed(1)
seqSeed <- sample(1:max(1e5,n.iter_sim),size=n.iter_sim,replace=FALSE)
iSeed <- seqSeed[iter_sim]
set.seed(iSeed)

cat("seed: ",iSeed,"\n")

## * path
path <- "."
path.res <- file.path(path,"Results","simulation-lvm")
if(dir.exists(path.res)==FALSE){
    dir.create(path.res)
}
path.output <- file.path(path,"output","simulation-lvm")
if(dir.exists(path.output)==FALSE){
    dir.create(path.output)
}

## * libraries
library(lava)
library(data.table)
devtools::load_all("lavaSearch2")

## * settings
seqN <- c(20,30,50,75,100,150,200,300,500)
n.boot <- 0
n.rep <- 5e2

## * model

## ** generative model
m.generative <- lvm(c(Y1,Y2,Y3,Y5) ~ eta1,
                    Y4 ~ 1,
                    c(Z1,Z2,Z3,Z4,Z5) ~ eta2,
                    eta1 ~ Age,
                    eta2 ~ Gender)
latent(m.generative) <- ~eta1+eta2
categorical(m.generative, labels = c("N","Y")) <- ~Gene1
categorical(m.generative, labels = c("N","Y")) <- ~Gene2
categorical(m.generative, labels = c("F","M")) <- ~Gender

## plot(m.generative)

## ** fit model
m.fit <- lvm(c(Y1,Y2,Y3,Y4,Y5) ~ eta1,
             c(Z1,Z2,Z3,Z4,Z5) ~ eta2,
             Y1 ~ Gene2,
             eta1 ~ Gene1 + Age,
             eta2 ~ Gender,
             eta1 ~ eta2)
covariance(m.fit) <- Y1 ~ Y2 
latent(m.fit) <- ~eta1+eta2
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

true.coef <- c(Y2 = 0,
               Y3 = 0,
               Y4 = 0,
               Y5 = 0,
               eta1 = 0,
               Z2 = 0,
               Z3 = 0,
               Z4 = 0,
               Z5 = 0,
               eta2 = 0,
               "Y1~Gene2Y" = 0,
               "Y2~eta1" = 1,
               "Y3~eta1" = 1,
               "Y4~eta1" = 0,
               "Y5~eta1" = 1,
               "eta1~eta2" = 0,
               "eta1~Age" = 1,
               "eta1~Gene1Y" = 0,
               "Z2~eta2" = 1,
               "Z3~eta2" = 1,
               "Z4~eta2" = 1,
               "Z5~eta2" = 1,
               "eta2~GenderM" = 1,
               "Y1~~Y1" = 1,
               "Y2~~Y2" = 1,
               "Y3~~Y3" = 1,
               "Y4~~Y4" = 1,
               "Y5~~Y5" = 1,
               "eta1~~eta1" = 1,
               "Z1~~Z1" = 1,
               "Z2~~Z2" = 1,
               "Z3~~Z3" = 1,
               "Z4~~Z4" = 1,
               "Z5~~Z5" = 1,
               "eta2~~eta2" = 1,
               "Y1~~Y2" = 0)

## * simulation
out <- calibrateType1(m.fit, true.coef = true.coef,
                      param = names(true.coef),
                      null = true.coef,
                      n = seqN, n.rep = n.rep, 
                      generative.object = m.generative, 
                      dir.save = path.res, label.file = iter_sim,
                      bootstrap = FALSE, seed = NULL, trace = 2)

## * display
print(sessionInfo())
