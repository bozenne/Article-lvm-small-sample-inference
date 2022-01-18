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
library(lavaSearch2)

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

## ** fit model
m.fit <- lvm(c(Y1,Y2,Y3,Y4,Y5) ~ eta1,
             c(Z1,Z2,Z3,Z4,Z5) ~ eta2,
             Y1 ~ Gene2,
             eta1 ~ Gene1 + Age,
             eta2 ~ Gender,
             eta1 ~ eta2)
covariance(m.fit) <- Y1 ~ Y2 
latent(m.fit) <- ~eta1+eta2

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

## * running time
if(FALSE){

    ls.data <- lapply(seqN, function(iN){lava::sim(m.generative, n = iN, latent = FALSE)})
    names(ls.data) <- seqN
    ls.lvmfit <- lapply(ls.data, function(iData){lava::estimate(m.fit, data = iData)})
    names(ls.lvmfit) <- seqN

    ## algo 2
    if(packageVersion("lavaSearch2")>="2.0.0"){
        speedLvm <- microbenchmark::microbenchmark("20" = estimate2(ls.lvmfit[["20"]]),
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
        speedLvm <- microbenchmark::microbenchmark("20" = sCorrect(ls.lvmfit[["20"]]),
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
    
    saveRDS(speedLvm, file = file.path("Results","speed-Algo2-Lvm.rds"))
    ## Unit: seconds
    ## expr      min       lq     mean   median       uq      max neval     cld
    ##   20 1.334072 1.349604 1.445554 1.359906 1.387825 2.460707    50 a
    ##   30 1.400533 1.419628 1.505777 1.436489 1.453078 2.192798    50 ab
    ##   50 1.457068 1.484657 1.555015 1.494763 1.506217 2.418348    50 ab
    ##   75 1.540719 1.559081 1.663383 1.583262 1.695684 2.307196    50  bc
    ##  100 1.669573 1.704736 1.780416 1.715103 1.799912 2.331430    50   c
    ##  150 1.848124 1.877930 2.029438 1.907948 2.082418 2.650552    50    d
    ##  200 2.121476 2.161624 2.354556 2.183057 2.346099 3.517504    50     e
    ##  300 2.641857 2.689073 2.964438 2.733649 3.096009 4.339638    50      f
    ##  500 3.724511 3.814207 4.264608 3.952911 4.711809 6.242196    50       g
}
