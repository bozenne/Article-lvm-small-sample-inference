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
library(lavaSearch2)

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
 
## ** investigator model
m.fit <- lvm(c(Y1[mu1:sigma]~1*eta,
               Y2[mu2:sigma]~1*eta,
               Y3[mu3:sigma]~1*eta,
               eta~beta1*Age + beta2*Gene1))
latent(m.fit) <- ~eta
categorical(m.fit, labels = c("N","Y")) <- ~Gene1

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

## * running time
if(FALSE){

    ls.data <- lapply(seqN, function(iN){lava::sim(m.generative, n = iN, latent = FALSE)})
    names(ls.data) <- seqN
    ls.lvmfit <- lapply(ls.data, function(iData){lava::estimate(m.fit, data = iData)})
    names(ls.lvmfit) <- seqN

    ## algo 2
    if(packageVersion("lavaSearch2")>="2.0.0"){
        speedMM <- microbenchmark::microbenchmark("20" = estimate2(ls.lvmfit[["20"]]),
                                                  "30" = estimate2(ls.lvmfit[["30"]]),
                                                  "50" = estimate2(ls.lvmfit[["50"]]),
                                                  "75" = estimate2(ls.lvmfit[["75"]]),
                                                  "100" = estimate2(ls.lvmfit[["100"]]),
                                                  "150" = estimate2(ls.lvmfit[["150"]]),
                                                  "200" = estimate2(ls.lvmfit[["200"]]),
                                                  "300" = estimate2(ls.lvmfit[["300"]]),
                                                  "500" = estimate2(ls.lvmfit[["500"]]),
                                                  times = 50
                                                  )
    }else{
        speedMM <- microbenchmark::microbenchmark("20" = sCorrect(ls.lvmfit[["20"]]),
                                                  "30" = sCorrect(ls.lvmfit[["30"]]),
                                                  "50" = sCorrect(ls.lvmfit[["50"]]),
                                                  "75" = sCorrect(ls.lvmfit[["75"]]),
                                                  "100" = sCorrect(ls.lvmfit[["100"]]),
                                                  "150" = sCorrect(ls.lvmfit[["150"]]),
                                                  "200" = sCorrect(ls.lvmfit[["200"]]),
                                                  "300" = sCorrect(ls.lvmfit[["300"]]),
                                                  "500" = sCorrect(ls.lvmfit[["500"]]),
                                                  times = 50
                                                  )
    }

    saveRDS(speedMM, file = file.path("Results","speed-Algo2-MM.rds"))
    ##  Unit: milliseconds
    ## expr       min        lq      mean    median        uq       max neval     cld
    ##   20  69.39570  69.81790  77.56782  70.57196  83.42209 106.99486    50 a c
    ##   30  65.93456  66.25498  74.50335  75.90799  79.75036 109.47599    50 a
    ##   50  67.34745  67.79331  75.19469  74.24777  81.23589 110.43351    50 ab
    ##   75  73.80876  74.39447  80.65022  74.81652  87.78234  97.88324    50 c
    ##  100  73.41202  74.04703  80.18060  74.54868  87.18580 100.43682    50 bc
    ##  150  84.22918  84.84177  91.80059  88.15417  98.20942 110.56206    50 d
    ##  200  95.16414  96.03317 103.78675 102.60769 109.54751 124.72909    50 e
    ##  300 100.85332 101.92437 110.57415 114.17656 115.49912 127.33436    50 f
    ##  500 135.84538 136.51482 146.19712 147.00547 150.89283 168.85675    50 g
}
