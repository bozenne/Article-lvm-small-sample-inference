## path <- "P:/Cluster/LVMproject/article-smallSampleInference"
## setwd(path)
## source("BATCH_IV-non-normal3.R")

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
path.res <- file.path(path,"Results","IV-non-normal3")
if(dir.exists(path.res)==FALSE){
    dir.create(path.res)
}
path.output <- file.path(path,"output","IV-non-normal3")
if(dir.exists(path.output)==FALSE){
    dir.create(path.output)
}

## * libraries
library(lava)
library(data.table)
library(multcomp)
library(MIIVsem)
library(lavaSearch2) ##

## * settings
seqN <- c(20,30,50,100,300,500,1000)
n.boot <- 0
n.rep <- 5e2

## * model

## ** generative model
m.generative <- lvm(c(Y1.aa~eta.aa,
                      Y2.aa~eta.aa,
                      Y3.aa~eta.aa,
                      Y4.aa~1,
                      eta.aa~Age+2*eta.id,
                      Y1.bb~eta.bb,
                      Y2.bb~eta.bb,
                      Y3.bb~eta.bb,
                      Y4.bb~1,
                      eta.bb~Age+2*eta.id,
                      Y1.cc~eta.cc,
                      Y2.cc~eta.cc,
                      Y3.cc~eta.cc,
                      Y4.cc~1,
                      eta.cc~Age+2*eta.id,
                      Y1.dd~eta.dd,
                      Y2.dd~eta.dd,
                      Y3.dd~eta.dd,
                      Y4.dd~1,
                      eta.dd~Age+2*eta.id,
                      Y1.ee~eta.ee,
                      Y2.ee~eta.ee,
                      Y3.ee~eta.ee,
                      Y4.ee~1,
                      eta.ee~Age+2*eta.id))
latent(m.generative) <- ~eta.aa + eta.bb + eta.cc + eta.dd + eta.ee + eta.id
categorical(m.generative, labels = c("N","Y")) <- ~Gene1
categorical(m.generative, labels = c("N","Y")) <- ~Gene2
transform(m.generative, Gene1num~Gene1) <- function(x){as.numeric(x[,1])-1}
transform(m.generative, Gene2num~Gene2) <- function(x){as.numeric(x[,1])-1}
transform(m.generative, id~eta.id) <- function(x){paste0("Id",1:NROW(x))}

sim2 <- function(n){
    dW <- as.data.table(sim(m.generative, n = n, latent = FALSE))
    dL <- melt(dW,
               measure=patterns("^Y1", "^Y2","^Y3","^Y4"),
               id.vars = c("Age","Gene1","Gene2","Gene1num","Gene2num","id"),
               value.name = c("Y1","Y2","Y3","Y4")
               )
    return(dL)
}
## ggplot(sim2(100), aes(x=id, y = Y1)) + geom_boxplot()

## ** fit model
m.fit <- lvm(c(Y1~eta,
               Y2~eta+Gene1,
               Y3~eta,
               Y4~eta,
               eta~Age+Gene2))
latent(m.fit) <- ~eta
## ee <- lava::estimate(m.fit, sim2(n = 1e3))
## butils::qqplot2(ee)
## hist(lava::sim(m.generative, n = 1e3)$Y1)

m.fit2 <- '
    eta =~ Y1 + Y2 + Y3 + Y4  
    eta ~ Age + Gene2num
    Y2 ~ Gene1num
  '


## * simulation
dt.res <- NULL
df.null <- rbind(data.frame(lava.name = "eta", lavaan.name = "eta~1", type = "alpha", null = 0, stringsAsFactors = FALSE),
                 data.frame(lava.name = "Y2", lavaan.name = "Y2~1", type = "nu", null = 0, stringsAsFactors = FALSE),
                 data.frame(lava.name = "Y3", lavaan.name = "Y3~1", type = "nu", null = 0, stringsAsFactors = FALSE),
                 data.frame(lava.name = "Y4", lavaan.name = "Y4~1", type = "nu", null = 0, stringsAsFactors = FALSE),
                 data.frame(lava.name = "eta~Age", lavaan.name = "eta~Age", type = "Gamma", null = 1, stringsAsFactors = FALSE),
                 data.frame(lava.name = "eta~Gene2Y", lavaan.name = "eta~Gene2num", type = "Gamma", null = 0, stringsAsFactors = FALSE),
                 data.frame(lava.name = "Y2~eta", lavaan.name = "Y2~eta", type = "Lambda", null = 1, stringsAsFactors = FALSE),
                 data.frame(lava.name = "Y2~Gene1Y", lavaan.name = "Y2~Gene1num", type = "K", null = 0, stringsAsFactors = FALSE),
                 data.frame(lava.name = "Y3~eta", lavaan.name = "Y3~eta", type = "Lambda", null = 1, stringsAsFactors = FALSE),
                 data.frame(lava.name = "Y4~eta", lavaan.name = "Y4~eta", type = "Lambda", null = 0, stringsAsFactors = FALSE)
                 )

for(iN in 1:length(seqN)){ ## iN <- 1
    
    cat("sample size: ",seqN[iN],"\n",sep = "")
    
    for(iRep in 1:n.rep){
        cat("*")
        ## simulation
        iD <- sim2(n = seqN[iN])

        ## model fit
        e.ML <- estimate(m.fit, data = iD, estimator = "Gaussian")
        if(e.ML$opt$convergence>0){
            e.ML <- estimate(m.fit, data = iD, estimator = "Gaussian", control = list(constrain = TRUE))
        }
        e.IV <- estimate(m.fit, data = iD, estimator = "IV")
        e.IVlavaan <- miive(model = m.fit2, data = iD)
        ## coef(e.IV)[c("Y2","eta~Gene2Y","Y2~Gene1Y","Y4~eta")]
        ## coef(e.IVlavaan)[c("Y2~1","eta~Gene2num","Y2~Gene1num","Y4~eta")]
        ## diag(e.IVlavaan$coefCov)[c("Y2~1","eta~Gene2num","Y2~Gene1num","Y4~eta")]
        ## diag(vcov(e.IV))[c("Y2","eta~Gene2Y","Y2~Gene1Y","Y4~eta")]

        ## diagnostics
        p.testNorm <- max(apply(residuals(e.ML),2,function(x){shapiro.test(x)$p.value}))

        ## ML estimates
        eS.ML <- summary(e.ML)$coef[df.null$lava.name,]
        eS.ML[,"Z-value"] <- (eS.ML[,"Estimate"]-df.null$null)/eS.ML[,"Std. Error"]
        eS.ML[,"P-value"] <- 2*(1-pnorm(abs(eS.ML[,"Z-value"])))
        
        ## ML robust estimates
        eS.robustML <- cbind(estimate(e.ML, id = iD$id)$coefmat, "Z-value" = NA)[df.null$lava.name,]
        eS.robustML[,"Z-value"] <- (eS.robustML[,"Estimate"]-df.null$null)/eS.robustML[,"Std.Err"]
        eS.robustML[,"P-value"] <- 2*(1-pnorm(abs(eS.robustML[,"Z-value"])))

        ## corrected ML robust estimates
        if(seqN[iN]<200){
            eS.robustMLC <- summary2(e.ML, cluster = iD$id, robust = TRUE)$coef[df.null$lava.name,]
            eS.robustMLC[,"t-value"] <- (eS.robustMLC[,"Estimate"]-df.null$null)/eS.robustMLC[,"robust SE"]
            eS.robustMLC[,"P-value"] <- 2*(1-pt(abs(eS.robustMLC[,"t-value"]), df = eS.robustMLC[,"df"]))
        }else{
            eS.robustMLC <- eS.robustML
            eS.robustMLC[] <- NA
        }
        ## IV estimates
        eS.IV <- summary(e.IV)$coef[df.null$lava.name,]
        eS.IV[,"Z-value"] <- (eS.IV[,"Estimate"]-df.null$null)/eS.IV[,"Std. Error"]
        eS.IV[,"P-value"] <- 2*(1-pnorm(abs(eS.IV[,"Z-value"])))
        
        eS.IVlavaan <- cbind("Estimate" = coef(e.IVlavaan)[df.null$lavaan.name],
                             "Std. Error" = sqrt(diag(e.IVlavaan$coefCov))[df.null$lavaan.name],
                             "Z-value" = NA,
                             "P-value" = NA)
        eS.IVlavaan[,"Z-value"] <- (eS.IVlavaan[,"Estimate"]-df.null$null)/eS.IVlavaan[,"Std. Error"]
        eS.IVlavaan[,"P-value"] <- 2*(1-pnorm(abs(eS.IVlavaan[,"Z-value"])))
        

        
        ## if(p.testNorm>0.05){browser()}
        ## qqtest::qqtest(residuals(e.ML)[,4])
        iDT.ML <- data.table(seed = iSeed,
                             n = seqN[iN],
                             rep = iRep,
                             link = df.null$lava.name,
                             type = df.null$type,
                             estimate = eS.ML[,"Estimate"],
                             p.value = eS.ML[,"P-value"],
                             estimator = "ML")
        iDT.robustML <- data.table(seed = iSeed,
                                   n = seqN[iN],
                                   rep = iRep,
                                   link = df.null$lava.name,
                                   type = df.null$type,
                                   estimate = eS.robustML[,"Estimate"],
                                   p.value = eS.robustML[,"P-value"],
                                   estimator = "robustML")
        iDT.robustMLC <- data.table(seed = iSeed,
                                    n = seqN[iN],
                                    rep = iRep,
                                    link = df.null$lava.name,
                                    type = df.null$type,
                                    estimate = eS.robustMLC[,"Estimate"],
                                    p.value = eS.robustMLC[,"P-value"],
                                    estimator = "robustMLC")
        iDT.IV <- data.table(seed = iSeed,
                             n = seqN[iN],
                             rep = iRep,
                             link = df.null$lava.name,
                             type = df.null$type,
                             estimate = eS.IV[,"Estimate"],
                             p.value = eS.IV[,"P-value"],
                             estimator = "IV")
        iDT.IVlavaan <- data.table(seed = iSeed,
                                   n = seqN[iN],
                                   rep = iRep,
                                   link = df.null$lava.name,
                                   type = df.null$type,
                                   estimate = eS.IVlavaan[,"Estimate"],
                                   p.value = eS.IVlavaan[,"P-value"],
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


