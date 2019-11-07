path.data <- "./data/"
path.data2 <- "./data2/"

## * Packages
library(lavaSearch2)
library(lme4)
library(lmerTest)
library(data.table)
library(foreign)

## * Application A

## ** Import
dtL.vitamin <- fread(file.path(path.data, "vitamin.txt"), header = TRUE)

## ** Data processing
dtL.vitamin[, animal := as.factor(animal)]
dtL.vitamin[, grp := factor(grp, levels = 1:2, labels = c("C","T"))]
dtL.vitamin[, treat := factor(grp, levels = c("C","T"), labels = c("No","Yes"))]
dtL.vitamin[week<5, treat := "No"] ## no treatment before week five
dtL.vitamin[, week.factor := paste0("w",as.factor(week))]
dtL.vitamin[, weight0 := scale(weight)]
dtL.vitamin[,interaction := paste0(treat,week.factor)]
dtL.vitamin[treat=="No" | week <5, interaction := "base"]

dtW.vitamin <- dcast(dtL.vitamin, value.var = "weight0",
                     formula = grp + animal ~ week.factor)

## ** Random intercept model
vitamin.REML <- lmer(weight0 ~ -1 + week.factor + interaction + (1|animal),
                     data = dtL.vitamin, REML = TRUE)
## summary(vitamin.REML)
coef.vitamin.REML <- c(fixef(vitamin.REML), 
                       sigma2 = sigma(vitamin.REML)^2, 
                       tau = as.double(summary(vitamin.REML)$varcor))
## coef.vitamin.REML
Ftest.vitamin.REML <- anova(vitamin.REML)

## ** Latent variable model 
m.vitamin <- lvm(w1[mu1:sigma] ~ 1*eta,
                 w3[mu3:sigma] ~ 1*eta,
                 w4[mu4:sigma] ~ 1*eta,
                 w5[mu5:sigma] ~ grp + 1*eta,
                 w6[mu6:sigma] ~ grp + 1*eta,
                 w7[mu7:sigma] ~ grp + 1*eta,
                 eta ~ 0)
latent(m.vitamin) <- ~eta

vitamin.ML <- estimate(m.vitamin, data = dtW.vitamin)
## summary(vitamin.ML)

coef.vitamin.ML <- coef(vitamin.ML)
## coef.vitamin.ML
Ftest.vitamin.ML <- lava::compare(vitamin.ML, par = c("w5~grpT","w6~grpT","w7~grpT")) ## 15.072/3

## ** Small sample correction
sCorrect(vitamin.ML) <- TRUE

coef.vitamin.MLc <- vitamin.ML$sCorrect$param

Ftest.vitamin.MLc <- compare2(vitamin.ML, par = c("w5~grpT","w6~grpT","w7~grpT")) ## 15.072/3

## * Application B

## ** Import
if(dir.exists(path.data2)){
    dt0.bdnf <- fread(file.path(path.data2,"SB_5HTTLPR_Patrick_Klaus.csv"))
}else{
    stop("Data relative to example 2 is under the Danish rules on data protection which do not allow to freely share person-sensitive health
care related data - contact the corresponding author of the paper for more information \n")
}

## ** Data processing
keep.variables <- c("cimbi.id","bmi","age","gender","scanner", ## clinical variables
                    "sb.injmass","sb.wa.injmass", ## injected mass
                    "httlpr","bdnf", ## genotype
                    "neo.bpnd","cau.bpnd","amy.bpnd","put.bpnd","hip.bpnd"
                    )

## Exclude indiv. w. high injected mass
## Exclude one "large old lady". She was excluded from Neuroimage manuscript analysis.
dt.bdnf <- dt0.bdnf[sb.injmass<6 & bmi<32,.SD,.SDcols = keep.variables]

## Mean-centering of continuous data
dt.bdnf[,age0 := age - mean(age,na.rm=TRUE)]
dt.bdnf[,inj0 := sb.wa.injmass - mean(sb.wa.injmass,na.rm=TRUE)]
dt.bdnf[,bmi0 := bmi - mean(bmi,na.rm=TRUE)]

## Create two-group factors for genotypes
dt.bdnf[,httlpr2 := factor(httlpr=="ll", labels=c("sx","ll"))]
dt.bdnf[,httlpr2 := relevel(httlpr2, ref = "ll")]
dt.bdnf[,bdnf2 := factor(bdnf=="val/val",labels=c("mx","vv"))]
dt.bdnf[,bdnf2 := relevel(bdnf2, ref = "vv")]

## Log transform bpnd
dt.bdnf[,neo := log(neo.bpnd)]
dt.bdnf[,cau := log(cau.bpnd)]
dt.bdnf[,amy := log(amy.bpnd)]
dt.bdnf[,put := log(put.bpnd)]
dt.bdnf[,hip := log(hip.bpnd)]

## ** Latent variable model
## definition (as in the original paper)
m.bdnf0 <- lvm(neo ~ u + httlpr2 + scanner,
               cau ~ u + scanner,
               put ~ u + scanner,
               hip ~ u + scanner,
               amy ~ u + scanner,
               u ~ age0 + gender + bdnf2 + bmi0 + inj0)
latent(m.bdnf0) <- ~u

m.bdnf1 <- m.bdnf0
covariance(m.bdnf1) <- cau ~ put
covariance(m.bdnf1) <- amy ~ hip

m.bdnf2 <- m.bdnf1
regression(m.bdnf2) <- cau ~ age0

## estimation
e.bdnf0 <- estimate(m.bdnf0, data = dt.bdnf, cluster = "cimbi.id",
                    control = list(constrain=TRUE))
e.bdnf1 <- estimate(m.bdnf1, data = dt.bdnf, cluster = "cimbi.id",
                    control = list(constrain=FALSE, start = coef(e.bdnf0)))
e.bdnf2 <- estimate(m.bdnf2, data = dt.bdnf, cluster = "cimbi.id",
                    control = list(constrain=FALSE, start = coef(e.bdnf1)))
## summary(e.bdnf2)

## coefficients
bdnf.null <- c('$\\gamma_2$' = 'u~bdnf2mx', '$k_1$' = 'neo~httlpr2sx')
coef.bdnf.ML <- coef(e.bdnf2)


## ** Small sample correction
sCorrect(e.bdnf2) <- TRUE
coef.bdnf.MLc <- e.bdnf2$sCorrect$param

## ** Permutation test [not in the article]
if(FALSE){
    library(pbapply)
    data0 <- dt.bdnf[,.SD,.SDcols = c("cimbi.id",manifest(m.bdnf2))]

    warper_resample <- function(var){
        var <- match.arg(var, c("httlpr2","bdnf2"))
    
        data <- data.table::copy(data0)
        data$repetition <- duplicated(data$cimbi.id)+1
        data[repetition==1, c(var) := .SD[sample.int(.N, replace = FALSE)], .SDcols = var]
        rep.id <- data[repetition==2,cimbi.id]
        data[cimbi.id %in% rep.id, c(var) := .SD[1], by = "cimbi.id", .SDcols = var]
        return(data)
    }

    warper_fit <- function(data){
        suppressWarnings(e0 <- estimate(m.bdnf0, data = data, cluster = "cimbi.id",
                                        control = list(constrain=TRUE)))
        suppressWarnings(e1 <- estimate(m.bdnf1, data = data, cluster = "cimbi.id",
                                        control = list(constrain=FALSE, start = coef(e0))))
        e2 <- estimate(m.bdnf2, data = data, cluster = "cimbi.id",
                       control = list(constrain=FALSE, start = coef(e1)))
        return(e2)
    }

## names(dt.bdnf)
## summary(warper_fit(data0))$coef[c("neo~httlpr2sx","u~bdnf2mx"),]

## dd <- warper_resample("bdnf2")
## summary(warper_fit(dd))$coef[c("neo~httlpr2sx","u~bdnf2mx"),]

## dd <- warper_resample("httlpr2")
## summary(warper_fit(dd))$coef[c("neo~httlpr2sx","u~bdnf2mx"),]

    n.perm <- 1e4

    httlpr.tps <- system.time({
        httlpr.perm <- pblapply(1:n.perm, function(x){
            iOut <- suppressWarnings(warper_fit(warper_resample(var = "httlpr")))
            return(iOut)
        })
    })
    attr(httlpr.perm,"tps") <- httlpr.tps
    saveRDS(httlpr.perm, file.path("Results","permTest-httlpr-factor.rds"))

    bdnf.tps <- system.time({
        bdnf.perm <- pblapply(1:n.perm, function(x){
            iOut <- suppressWarnings(warper_fit(warper_resample(var = "bdnf")))
            return(iOut)
        })
    })
    attr(bdnf.perm,"tps") <- bdnf.tps
    saveRDS(bdnf.perm, file.path("Results","permTest-bdnf-factor.rds"))


    ##summary(httlpr.perm[[1]])$coef[c("neo~httlpr2sx","u~bdnf2mx"),]
    ##summary(bdnf.perm[[1]])$coef[c("neo~httlpr2sx","u~bdnf2mx"),]

    httlpr.cv <- sapply(httlpr.perm,function(x){x$opt$convergence})
    bdnf.cv <- sapply(bdnf.perm,function(x){x$opt$convergence})
    ## mean(httlpr.cv)
    ## mean(bdnf.cv)

    mean( abs(do.call(rbind,lapply(httlpr.perm[httlpr.cv==0],coef))[,"neo~httlpr2sx"]) > abs(coef(e.bdnf2)["neo~httlpr2sx"]) )
    mean( abs(do.call(rbind,lapply(bdnf.perm[bdnf.cv==0],coef))[,"u~bdnf2mx"]) > abs(coef(e.bdnf2)["u~bdnf2mx"]) )
}


## * Application C
## ** Import
if(dir.exists(path.data2)){
    dt0.memory <- as.data.table(read.spss(file.path(path.data2,"SPSSmedscan_v2.sav"), to.data.frame = TRUE))
    dt0.gene <- fread(file.path(path.data2,"SB_5HTTLPR_Patrick_Klaus.csv"))
}else{ 
    stop("Data relative to example 3 is under the Danish rules on data protection which do not allow to freely share person-sensitive health
care related data - contact the corresponding author of the paper for more information \n")
}

## ** data processing
keep.col <- c(cimbi.id = "CIMBIID",
              amy = "Total_amygdala_SB_BPnd_NonPV_GM",
              acc = "Total_antcin_SB_BPnd_NonPV_GM",
              hip = "Total_hippocampus_SB_BPnd_NonPV_GM", 
              frc = "FrontalCortex_SB_BPnd_NonPV_GM", 
              pos.imm = "CAMTImmediaterecallofpositivewordssumA1A5", 
              pos.shr = "CAMTShorttermmemoryofpositivewordsA6", 
              pos.del = "CAMTDelayedrecallofpositivewordsA7", 
              neg.imm = "CAMTImmediaterecallofnegativewordssumA1A5", 
              neg.shr = "CAMTShorttermmemoryofnegativewordsA6", 
              neg.del = "CAMTDelayedrecallofnegativewordsA7", 
              neu.imm = "CAMTImmediaterecallofneutralwordssumA1A5", 
              neu.shr = "CAMTShorttermmemoryofneutralwordsA6", 
              neu.del = "CAMTDelayedrecallofneutralwordsA7", 
              age = "Age", 
              sex = "Gender", 
              bmi = 'BMI', 
              wa.inj = 'SBinjectedmassperkgmicrogramkg', 
              inj = 'SBinjectedmassmicrogram')

dt.memory <- dt0.memory[CIMBIID != 51906,.SD,.SDcols = as.character(keep.col)]
setnames(dt.memory, old = as.character(keep.col), new = names(keep.col))

## Align scaling of immediate memory performance with short-term and delayed.
dt.memory[,pos.imm5 := pos.imm/5]
dt.memory[,neg.imm5 := neg.imm/5]
dt.memory[,neu.imm5 := neu.imm/5]

## Mean-centered and log-transformed data.
dt.memory[,age0 := age - mean(age)]
dt.memory[,wa.inj0 := wa.inj - mean(wa.inj)]
dt.memory[,inj := inj - mean(inj)]

setkeyv(dt.memory, cols = "cimbi.id")

dt0.gene[, bdnf2 := factor(bdnf == "val/val", labels = c("mx","vv"))]
dt0.gene[, httlpr2 := factor(httlpr == "ll", labels = c("sx","ll"))]
dt.gene <- dt0.gene[cimbi.id %in% dt.memory$cimbi.id,.(cimbi.id,bdnf2,httlpr2)]
setkeyv(dt.gene, cols = "cimbi.id")

dtG.memory <- merge(dt.gene,dt.memory, by = "cimbi.id")

## ** latent variable model
## definition (as in the original paper)
m.memory0 <- lvm(frc ~ httlpr2 + u , amy ~ u, acc ~ u, hip ~ u)
regression(m.memory0) <- c(pos.imm5, pos.shr, pos.del) ~ m.pos
regression(m.memory0) <- c(neg.imm5, neg.shr, neg.del) ~ m.neg
regression(m.memory0) <- c(neu.imm5, neu.shr, neu.del) ~ m.neu
regression(m.memory0) <- u~age0
regression(m.memory0) <- m.pos~age0
regression(m.memory0) <- m.neg~age0
regression(m.memory0) <- m.neu~age0
latent(m.memory0) <- ~u + m.pos + m.neg + m.neu

m.memory <- m.memory0
covariance(m.memory) <- amy~hip
regression(m.memory) <- c(m.pos ~ u, m.neg ~ u, m.neu ~ u)

## estimation0
e.memory0 <- estimate(m.memory0, data = dtG.memory)

e.memory <- estimate(m.memory, data = dtG.memory,
                     control = list(constrain = TRUE, start = coef(e.memory0)))
## summary(e.memory)

## coefficients
memory.null <- c('$b_1$' = 'm.pos~u', '$b_2$' = 'm.neu~u', '$b_3$' = 'm.neg~u')
coef.memory.ML <- coef(e.memory)

## ** small sample correction
sCorrect(e.memory) <- TRUE
coef.memory.MLc <- e.memory$sCorrect$param
