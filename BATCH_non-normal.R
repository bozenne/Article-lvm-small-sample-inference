## path <- "P:/Cluster/LVMproject/article-smallSampleInference"
## setwd(path)
## source("BATCH_non-normal.R")

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
path.res <- file.path(path,"Results","non-normal")
if(dir.exists(path.res)==FALSE){
    dir.create(path.res)
}
path.output <- file.path(path,"output","non-normal")
if(dir.exists(path.output)==FALSE){
    dir.create(path.output)
}

## * libraries
library(lava)
library(data.table)
devtools::load_all("lavaSearch2")

## * settings
seqN <- c(500,1000,2000)
n.boot <- 0# 1e3
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

# plot(m.generative)

## ** fit model
m.fit <- lvm(c(Y1~eta,
               Y2~eta+Gene1,
               Y3~eta,
               Y4~eta,
               eta~Age+Gene2))
latent(m.fit) <- ~eta
distribution(m.generative,~Y1) <- lava:::student.lvm(df = 2)
distribution(m.generative,~Y2) <- lava:::student.lvm(df = 3)
distribution(m.generative,~Y3) <- lava:::student.lvm(df = 4)
distribution(m.generative,~Y4) <- lava:::student.lvm(df = 5)

## lava::estimate(m.fit, lava::sim(m.generative, n = 1e5))

## * simulation
coefnull <- c("Y2","Y2~Gene1Y","eta~Gene2Y","Y4~eta")
dt.res <- NULL

for(iN in 1:length(seqN)){ ## iN <- 1
    cat("sample size: ",seqN[iN],"\n")
    for(iSim in 1:n.rep){ ## iSim <- 1
        cat("*")
        d <- lava::sim(m.generative, n = seqN[iN])
        e.ML <- estimate(m.fit, data = d, estimator = "Gaussian")    
        e.IV <- estimate(m.fit, data = d, estimator = "IV")

        eS.ML <- summary(e.ML)$coef
        eS.RML <- estimate(e.ML)$coefmat
        eS.IV <- summary(e.IV)$coef

        iRes <- data.table(seed = iSeed,
                           n = seqN[iN],
                           rep = iSim,
                           coef = coefnull,                           
                           p.ML = eS.ML[coefnull, "P-value"],
                           p.RML = eS.RML[coefnull, "P-value"],
                           p.IV = eS.IV[coefnull, "P-value"])
        dt.res <- rbind(dt.res,iRes)
    }
    cat("\n")
    filename <- paste0("type1error-S",iter_sim,"(tempo).rds")
    saveRDS(dt.res, file = file.path(path.res,filename))

}

filename <- paste0("type1error-S",iter_sim,".rds")
saveRDS(dt.res, file = file.path(path.res,filename))

## * display
print(sessionInfo())


## * extra (lm)
if(FALSE){
    ## ** generative model
    m.generative <- lvm(Y~1,X1+X2+X3~1)
    ## distribution(m.generative,~Y) <- Gamma.lvm(rate=1,shape=0.5)
    distribution(m.generative,~Y) <- lava:::student.lvm(df = 2)

    ## dd <- sim(m.generative, 1e3)
    ## qqtest::qqtest(residuals(lm(Y~X1+X2+X3, data = dd)))

    ## ** investigator model
    m.fit <- lvm(Y~X1+X2+X3)

    ## ** simulation
    n.sim <- 4000
    n.obs <- 1000
    cpus <- 4

    cl <- snow::makeSOCKcluster(cpus)
    doSNOW::registerDoSNOW(cl)

    pb <- txtProgressBar(max = n.sim, style=3)
    opts <- list(progress = function(n) setTxtProgressBar(pb, n))


    ls.res <- foreach::`%dopar%`(
                           foreach::foreach(i=1:n.sim, .options.snow=opts, .packages = "lava"), {
                               d <- lava::sim(m.generative, n = n.obs)
                               e <- estimate(m.fit, data = d)
                               c("model" = summary(e)$coef[c("Y~X1","Y~X2","Y~X3"),"P-value"],
                                 "robust" = estimate(e)$coefmat[c("Y~X1","Y~X2","Y~X3"),"P-value"])
                           })

    vec.type1 <- apply(do.call(rbind,ls.res),2, function(vecP){mean(vecP<=0.05)})
    print(vec.type1)
    ##  model.Y~X1  model.Y~X2  model.Y~X3 robust.Y~X1 robust.Y~X2 robust.Y~X3 
    ## 0.05400     0.05625     0.05575     0.04875     0.04900     0.04550 

}
