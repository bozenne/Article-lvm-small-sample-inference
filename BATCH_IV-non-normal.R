## path <- "P:/Cluster/LVMproject/article-smallSampleInference"
## setwd(path)
## source("BATCH_IV-non-normal.R")

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
path.res <- file.path(path,"Results","IV-non-normal")
if(dir.exists(path.res)==FALSE){
    dir.create(path.res)
}
path.output <- file.path(path,"output","IV-non-normal")
if(dir.exists(path.output)==FALSE){
    dir.create(path.output)
}

## * libraries
library(lava)
library(data.table)
library(multcomp)
library(MIIVsem)
## devtools::load_all("lavaSearch2")
library(lavaSearch2)

## * settings
seqN <- c(20,30,50,100,300,500,1000)
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
transform(m.generative, Gene1num~Gene1) <- function(x){as.numeric(x[,1])-1}
transform(m.generative, Gene2num~Gene2) <- function(x){as.numeric(x[,1])-1}
distribution(m.generative, ~Y1) <- student.lvm(df = 2)
distribution(m.generative, ~Y2) <- student.lvm(df = 2)
distribution(m.generative, ~Y3) <- student.lvm(df = 2)
distribution(m.generative, ~Y4) <- student.lvm(df = 2)

# plot(m.generative)

## ** fit model
m.fit <- lvm(c(Y1~eta,
               Y2~eta+Gene1,
               Y3~eta,
               Y4~eta,
               eta~Age+Gene2))
latent(m.fit) <- ~eta
## ee <- lava::estimate(m.fit, lava::sim(m.generative, n = 1e3))
## butils::qqplot2(ee)

m.fit2 <- '
    eta =~ Y1 + Y2 + Y3 + Y4  
    eta ~ Age + Gene2num
    Y2 ~ Gene1num
  '

## * simulation
dt.res <- NULL
## paste(names(coef(e.LM)),collapse=",")
keep.coef <- c("eta","Y2","Y3","Y4","Y2~Gene1Y","eta~Age","eta~Gene2Y","Y2~eta","Y3~eta","Y4~eta")
keep.type <- c("alpha","nu","nu","nu","K","Gamma","Gamma","Lambda","Lambda","Lambda")
null <- c("Y2" = 0,
          "Y2~Gene1Y" = 0,
          "Y2~eta" = 1,
          "Y3" = 0,
          "Y3~eta" = 1,
          "Y4" = 0,
          "Y4~eta" = 0,
          "eta" = 0,
          "eta~Age" = 1,
          "eta~Gene2Y" = 0,
          "Y1~~Y1" = 1,  
          "eta~~eta" = 1,
          "Y2~~Y2" = 1,
          "Y3~~Y3" = 1,
          "Y4~~Y4" = 1   
          )

for(iN in 1:length(seqN)){ ## iN <- 1
    
    cat("sample size: ",seqN[iN],"\n",sep = "")
    
    for(iRep in 1:n.rep){
        cat("*")
        ## simulation
        iD <- sim(m.generative, n = seqN[iN], latent = FALSE)

        ## model fit        
        e.ML <- estimate(m.fit, data = iD, estimator = "Gaussian")
        if(e.ML$opt$convergence>0){
            e.ML <- estimate(m.fit, data = iD, estimator = "Gaussian", control = list(constrain = TRUE))
        }
        if(e.ML$opt$convergence>0){next}        
        if(min(eigen(information(e.ML))$value)<1e-5){next}

        e.IV <- estimate(m.fit, data = iD, estimator = "IV")
        e.IVlavaan <- miive(model = m.fit2, data = iD)
        ## coef(e.IV)[c("Y2","eta~Gene2Y","Y2~Gene1Y","Y4~eta")]
        ## coef(e.IVlavaan)[c("Y2~1","eta~Gene2num","Y2~Gene1num","Y4~eta")]
        ## diag(e.IVlavaan$coefCov)[c("Y2~1","eta~Gene2num","Y2~Gene1num","Y4~eta")]
        ## diag(vcov(e.IV))[c("Y2","eta~Gene2Y","Y2~Gene1Y","Y4~eta")]

        ## diagnostics
        p.testNorm <- max(apply(residuals(e.ML),2,function(x){shapiro.test(x)$p.value}))

        ## ML estimates
        eS.ML <- summary(e.ML)$coef
        eS.ML[names(null),"Z-value"] <- (eS.ML[names(null),"Estimate"]-null)/eS.ML[names(null),"Std. Error"]
        eS.ML[names(null),"P-value"] <- 2*(1-pnorm(abs(eS.ML[names(null),"Z-value"])))
        
        ## ML robust estimates
        eS.robustML <- cbind(estimate(e.ML)$coefmat, "Z-value" = NA)
        eS.robustML[names(null),"Z-value"] <- (eS.robustML[names(null),"Estimate"]-null)/eS.robustML[names(null),"Std.Err"]
        eS.robustML[names(null),"P-value"] <- 2*(1-pnorm(abs(eS.robustML[names(null),"Z-value"])))

        ## corrected ML robust estimates
        if(seqN[iN]<200){
            eS.robustMLC <- summary2(e.ML, robust = TRUE)$coef
            eS.robustMLC[names(null),"t-value"] <- (eS.robustMLC[names(null),"Estimate"]-null)/eS.robustMLC[names(null),"robust SE"]
            eS.robustMLC[names(null),"P-value"] <- 2*(1-pt(abs(eS.robustMLC[names(null),"t-value"]), df = eS.robustMLC[names(null),"df"]))
        }else{
            eS.robustMLC <- eS.robustML
            eS.robustMLC[] <- NA
        }
        
        ## IV estimates
        eS.IV <- summary(e.IV)$coef
        eS.IV[names(null),"Z-value"] <- (eS.IV[names(null),"Estimate"]-null)/eS.IV[names(null),"Std. Error"]
        eS.IV[names(null),"P-value"] <- 2*(1-pnorm(abs(eS.IV[names(null),"Z-value"])))

        eS.IVlavaan <- cbind("Estimate" = c(coef(e.IVlavaan), rep(NA, 5)),
                             "Std. Error" = c(sqrt(diag(e.IVlavaan$coefCov)), rep(NA,5)),
                             "Z-value" = NA,
                             "P-value" = NA)
        eS.IVlavaan[,"Z-value"] <- (eS.IVlavaan[,"Estimate"]-null)/eS.IVlavaan[,"Std. Error"]
        eS.IVlavaan[,"P-value"] <- 2*(1-pnorm(abs(eS.IVlavaan[,"Z-value"])))
        

        
        ## if(p.testNorm>0.05){browser()}
        ## qqtest::qqtest(residuals(e.ML)[,4])
        iDT.ML <- data.table(seed = iSeed,
                             n = seqN[iN],
                             rep = iRep,
                             link = keep.coef,
                             type = keep.type,
                             estimate = eS.ML[keep.coef,"Estimate"],
                             p.value = eS.ML[keep.coef,"P-value"],
                             estimator = "ML")
        iDT.robustML <- data.table(seed = iSeed,
                                   n = seqN[iN],
                                   rep = iRep,
                                   link = keep.coef,
                                   type = keep.type,
                                   estimate = eS.robustML[keep.coef,"Estimate"],
                                   p.value = eS.robustML[keep.coef,"P-value"],
                                   estimator = "robustML")
        iDT.robustMLC <- data.table(seed = iSeed,
                                    n = seqN[iN],
                                    rep = iRep,
                                    link = keep.coef,
                                    type = keep.type,
                                    estimate = eS.robustMLC[keep.coef,"Estimate"],
                                    p.value = eS.robustMLC[keep.coef,"P-value"],
                                    estimator = "robustMLC")

        iDT.IV <- data.table(seed = iSeed,
                             n = seqN[iN],
                             rep = iRep,
                             link = keep.coef,
                             type = keep.type,
                             estimate = eS.IV[keep.coef,"Estimate"],
                             p.value = eS.IV[keep.coef,"P-value"],
                             estimator = "IV")
        iDT.IVlavaan <- data.table(seed = iSeed,
                                   n = seqN[iN],
                                   rep = iRep,
                                   link = keep.coef,
                                   type = keep.type,
                                   estimate = eS.IVlavaan[1:length(keep.coef),"Estimate"],
                                   p.value = eS.IVlavaan[1:length(keep.coef),"P-value"],
                                   estimator = "IVlavaan")
        iDT <- rbind(iDT.ML,
                     iDT.robustML,
                     iDT.robustMLC,
                     iDT.IV,
                     iDT.IVlavaan)

        dt.res <- rbind(dt.res,
                        cbind(iDT,
                              shapiroMax = p.testNorm,
                              cvML = e.ML$opt$convergence)
                        )
    }

    filename <- paste0("type1error-S",iter_sim,"(tempo).rds")
    saveRDS(dt.res, file = file.path(path.res,filename))

}
filename <- paste0("type1error-S",iter_sim,".rds")
saveRDS(dt.res, file = file.path(path.res,filename))

## * display
print(sessionInfo())


## * extra
if(FALSE){
    seqN <- c(20,30,50,100)

    ## ** model
    library(lava)
    m1 <- lvm(c(Y1,Y2,Y3,Y4) ~ eta)
    regression(m1) <- eta ~ X1 + X2
    regression(m1) <- Y1 ~ X1
    latent(m1) <- ~eta

    m2 <- lvm(c(Y1,Y2,Y3,Y4) ~ eta)
    regression(m2) <- eta ~ X1 + X2
    regression(m2) <- Y1 ~ X3
    regression(m2) <- Y2 ~ X4
    regression(m2) <- Y3 ~ X5
    regression(m2) <- Y4 ~ X6
    latent(m2) <- ~eta

    ## ** warper
    runML1 <- function(..., n){
        n <- 1000
        d <- lava::sim(m1, n = n, latent = FALSE)
        e <- estimate(m1, data = d, estimator = "gaussian")
        e2 <- estimate(m1, data = d, estimator = "IV")
        summary(e2)
        return(c("ML"=coef(e)[c("Y4~X1","eta~X1","eta~X2")],
                 "IV"=coef(e2)[c("Y4~X1","eta~X1","eta~X2")]))
    }
    runML2 <- function(..., n){
        d <- lava::sim(m2, n = n, latent = FALSE)
        e <- estimate(m2, data = d, estimator = "gaussian")
        e2 <- estimate(m2, data = d, estimator = "IV")
        return(c("ML"=coef(e)[c("Y1~X3","Y","eta~X1","eta~X2")],
                 "IV"=coef(e2)[c("Y1~X3","eta~X1","eta~X2")]))
    }

    ## ** run
    set.seed(17)
    ls.sim <- lapply(seqN, function(iN){
        sim(runML, m = m1, n = iN, mc.cores = 4)
    })
    names(ls.sim) <- seqN
    lapply(ls.sim,summary)

set.seed(17)
simres <- lava::sim(runML2, n = 100, mc.cores = 4)


library(MIIVsem)
model1 <- '
    eta =~ y1 + y2  + y3  + y4  
    eta ~ X1 + X2
    y1 ~ X1
  '

model2 <- '
    eta =~ y1 + y2 + y3 + y4  
    eta ~ X1 + X2
    y1 ~ X3
    y2 ~ X4
    y3 ~ X5
    y4 ~ X6
  '

dd1 <- sim(m1, n = 100)
e.lava <- estimate(m1, data = dd1, estimator = "IV")
e.lavaan <- miive(model = model1, data = dd1)


dd2 <- sim(m2, n = 100)
e.lavaan <- miive(model = model2, data = dd2)


## from the help of MIIVsem:::miive
## Certain model specifications are not currently supported. For example, the scaling indicator of
## a latent variable is not permitted to cross-load on another latent variable. 

## In addition, MIIVsem does not currently support relations where the scaling indicator of a
## latent variable is also the dependent variable in a regression equation. The model below would
## not be valid under the current algorithm
}
