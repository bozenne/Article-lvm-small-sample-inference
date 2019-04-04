## path <- "P:/Cluster/LVMproject/article-smallSampleInference"
## setwd(path)
## source("BATCH_simulation-factorModel.R")

rm(list = ls())
gc()

## * seed
iter_sim <- as.numeric(Sys.getenv("SGE_TASK_ID"))
n.iter_sim <- as.numeric(Sys.getenv("SGE_TASK_LAST"))
if(is.na(iter_sim)){iter_sim <- 2}
if(is.na(n.iter_sim)){n.iter_sim <- 10}
cat("iteration ",iter_sim," over ",n.iter_sim,"\n",sep="")

set.seed(1)
seqSeed <- sample(1:max(1e5,n.iter_sim),size=n.iter_sim,replace=FALSE)
iSeed <- seqSeed[iter_sim]
set.seed(iSeed)

cat("seed: ",iSeed,"\n")

## * path
path <- "."
path.res <- file.path(path,"Results","simulation-factorModel")
if(dir.exists(path.res)==FALSE){
    dir.create(path.res)
}
path.output <- file.path(path,"output","simulation-factorModel")
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
m.generative <- lvm(c(Y1~eta,
                      Y2~eta,
                      Y3~eta,
                      Y4~1,
                      eta~Age))
latent(m.generative) <- ~eta
categorical(m.generative, labels = c("N","Y")) <- ~Gene1
categorical(m.generative, labels = c("N","Y")) <- ~Gene2

## ** fit model
m.fit <- lvm(c(Y1~eta+Gene2,
               Y2~eta,
               Y3~eta,
               Y4~eta,
               eta~Age+Gene1))
latent(m.fit) <- ~eta

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
               Y4 = 0,
               "Y1~Gene2Y" = 0,
               "eta~Age" = 1,
               "eta~Gene1Y" = 0,
               "Y2~eta" = 1,
               "Y3~eta" = 1,
               "Y4~eta" = 0,
               "Y1~~Y1" = 1,
               "eta~~eta" = 1,
               "Y2~~Y2" = 1,
               "Y3~~Y3" = 1,
               "Y4~~Y4" = 1)

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
    speedFactor <- microbenchmark::microbenchmark("20" = sCorrect(ls.lvmfit[["20"]]),
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
    saveRDS(speedFactor, file = file.path("Results","speed-Algo2-Factor.rds"))
 ##    Unit: milliseconds
 ## expr      min       lq     mean   median       uq      max neval cld
 ##   20 153.4920 165.4882 167.9470 168.0757 173.8985 181.4775    50 b
 ##   30 139.3625 143.3120 153.0009 154.8625 156.3171 176.4130    50 a
 ##   50 140.7812 143.2224 155.0999 155.5626 157.3227 193.3208    50 a
 ##   75 154.9952 170.4748 174.3150 174.5028 180.2789 210.7192    50 b
 ##  100 170.1527 183.2981 192.2278 186.0499 195.3349 389.3620    50 c
 ##  150 175.7185 189.7478 195.4327 194.0275 202.5155 223.1499    50 c
 ##  200 202.3027 215.0806 220.8639 220.7766 225.8146 252.8478    50 d
 ##  300 262.3011 266.8939 272.6470 271.2441 275.5871 292.8013    50 e
 ##  500 360.1440 371.8377 384.6029 380.2178 390.5027 556.9517    50 f


    ## algo 1
    microbenchmark::microbenchmark("20" = sCorrect(ls.lvmfit[["20"]], adjust.n = FALSE),
                                   "30" = sCorrect(ls.lvmfit[["30"]], adjust.n = FALSE),
                                   "50" = sCorrect(ls.lvmfit[["50"]], adjust.n = FALSE),
                                   "75" = sCorrect(ls.lvmfit[["75"]], adjust.n = FALSE),
                                   "100" = sCorrect(ls.lvmfit[["100"]], adjust.n = FALSE),
                                   "150" = sCorrect(ls.lvmfit[["150"]], adjust.n = FALSE),
                                   "200" = sCorrect(ls.lvmfit[["200"]], adjust.n = FALSE),
                                   "300" = sCorrect(ls.lvmfit[["300"]], adjust.n = FALSE),
                                   "500" = sCorrect(ls.lvmfit[["500"]], adjust.n = FALSE),
                                   times = 50
                                   )

 ##    Unit: milliseconds
 ## expr      min       lq     mean   median       uq      max neval cld
 ##   20 137.4206 139.4769 151.5129 151.5326 159.3462 179.8258    50 cd
 ##   30 124.6302 125.9618 136.9098 138.2271 144.7811 164.9541    50 ab
 ##   50 120.9232 122.0929 131.6263 132.6343 135.1917 166.9360    50 a
 ##   75 126.1521 127.8582 139.4041 140.5283 146.9436 166.2210    50 ab
 ##  100 131.2951 141.6610 145.3657 145.1715 153.7632 157.3701    50 bd
 ##  150 131.5896 132.9748 142.7944 145.3925 146.9813 159.4823    50 bc
 ##  200 141.1600 151.5561 154.4614 154.8508 157.2333 171.8212    50 d
 ##  300 159.7682 170.4607 179.4216 175.0345 183.6101 356.9537    50 e
 ##  500 207.7388 221.7226 229.4150 223.6882 229.9727 410.6021    50 f


}
