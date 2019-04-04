## path <- "P:/Cluster/LVMproject/article-smallSampleInference"
## setwd(path)
## source("BATCH_comparison-ML-IV-GLS.R")

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
path.res <- file.path(path,"Results","comparison-ML-IV-GLS")
if(dir.exists(path.res)==FALSE){
    dir.create(path.res)
}
path.output <- file.path(path,"output","comparison-ML-IV-GLS")
if(dir.exists(path.output)==FALSE){
    dir.create(path.output)
}

## * libraries
library(lava)
library(lavaan)
library(data.table)
library(MIIVsem)

## * settings
seqN <- c(20,30,50,75,100)
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

# plot(m.generative)

## ** fit model
m.fit <- lvm(c(Y1~eta,
               Y2~eta+Gene1,
               Y3~eta,
               Y4~eta,
               eta~Age+Gene2))
latent(m.fit) <- ~eta

m.fit2 <- '
    eta =~ Y1 + Y2 + Y3 + Y4  
    eta ~ Age + Gene2num
    Y2 ~ Gene1num
  '

m.fit3 <- '
    eta =~ Y1 + Y2 + Y3 + Y4  
    eta ~ 1 + Age + Gene2num
    Y1 ~ 0
    Y2 ~ 1 + Gene1num
    Y3 ~ 1
    Y4 ~ 1
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
        iD <- sim(m.generative, n = seqN[iN], latent = FALSE)

        ## model fit
        e.ML <- estimate(m.fit, data = iD, estimator = "Gaussian")
        ## logLik(cfa(m.fit3, data = iD, estimator = "ML"))-logLik(e.ML)
        if(e.ML$opt$convergence>0){
            e.ML <- estimate(m.fit, data = iD, estimator = "Gaussian", control = list(constrain = TRUE))
        }
        e.IV <- estimate(m.fit, data = iD, estimator = "IV")
        e.IVlavaan <- miive(model = m.fit2, data = iD)
        e.GLS <- cfa(m.fit3, data = iD, estimator = "GLS")
        e.WLS <- suppressWarnings(try(cfa(m.fit3, data = iD, estimator = "WLS"), silent = TRUE))
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
        eS.robustML <- cbind(estimate(e.ML)$coefmat[df.null$lava.name,], "Z-value" = NA)
        eS.robustML[,"Z-value"] <- (eS.robustML[,"Estimate"]-df.null$null)/eS.robustML[,"Std.Err"]
        eS.robustML[,"P-value"] <- 2*(1-pnorm(abs(eS.robustML[,"Z-value"])))

        ## GLS estimates
        eS0.GLS <- parameterEstimates(e.GLS)
        eS0.GLS[grep("=~",eS0.GLS$op),c("lhs","rhs")] <- cbind(eS0.GLS[grep("=~",eS0.GLS$op),c("rhs","lhs")])
        rownames(eS0.GLS) <- paste0(eS0.GLS$lhs,gsub("=~","~",eS0.GLS$op), eS0.GLS$rhs)
        
        eS.GLS <- cbind("Estimate" = eS0.GLS[df.null$lavaan.name,"est"],
                        "Std. Error" = eS0.GLS[df.null$lavaan.name,"se"],
                        "Z-value" = NA,
                        "P-value" = NA)
        rownames(eS.GLS) <- df.null$lava.name
        eS.GLS[,"Z-value"] <- (eS.GLS[,"Estimate"]-df.null$null)/eS.GLS[,"Std. Error"]
        eS.GLS[,"P-value"] <- 2*(1-pnorm(abs(eS.GLS[,"Z-value"])))

        ## WLS estimates
        if(inherits(e.WLS,"try-error")){
            eS.WLS <- cbind("Estimate" = rep(NA, NROW(df.null)),
                            "P-value" = rep(NA, NROW(df.null))
                            )
            rownames(eS.WLS) <- df.null$lava.name            
        }else{
            eS0.WLS <- parameterEstimates(e.WLS)
            eS0.WLS[grep("=~",eS0.WLS$op),c("lhs","rhs")] <- cbind(eS0.WLS[grep("=~",eS0.WLS$op),c("rhs","lhs")])
            rownames(eS0.WLS) <- paste0(eS0.WLS$lhs,gsub("=~","~",eS0.WLS$op), eS0.WLS$rhs)
        
            eS.WLS <- cbind("Estimate" = eS0.WLS[df.null$lavaan.name,"est"],
                            "Std. Error" = eS0.WLS[df.null$lavaan.name,"se"],
                            "Z-value" = NA,
                            "P-value" = NA)
            rownames(eS.WLS) <- df.null$lava.name
            eS.WLS[,"Z-value"] <- (eS.WLS[,"Estimate"]-df.null$null)/eS.WLS[,"Std. Error"]
            eS.WLS[,"P-value"] <- 2*(1-pnorm(abs(eS.WLS[,"Z-value"])))
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
        iDT.GLS <- data.table(seed = iSeed,
                              n = seqN[iN],
                              rep = iRep,
                              link = df.null$lava.name,
                              type = df.null$type,
                              estimate = eS.GLS[,"Estimate"],
                              p.value = eS.GLS[,"P-value"],
                              estimator = "GLS")
        iDT.WLS <- data.table(seed = iSeed,
                              n = seqN[iN],
                              rep = iRep,
                              link = df.null$lava.name,
                              type = df.null$type,
                              estimate = eS.WLS[,"Estimate"],
                              p.value = eS.WLS[,"P-value"],
                              estimator = "WLS")
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
                     iDT.GLS,
                     iDT.WLS,
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


