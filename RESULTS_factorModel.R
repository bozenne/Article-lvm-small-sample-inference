path.results <- "./Results"
path.data <- "./data2/"

## * packages
library(lavaSearch2)
library(data.table)

## * Import

## ** simulation results
dttype1.sim.factor <- readRDS(file.path(path.results, "type1error-simulation-factorModel.rds"))
dtbias.sim.factor <- readRDS(file.path(path.results, "bias-simulation-factorModel.rds"))
param.sim.factor <- readRDS(file.path(path.results, "param-simulation-factorModel.rds"))
dttype1.ill.factor <- readRDS(file.path(path.results, "type1error-illustration-factorModel.rds"))

## ** real data
if(dir.exists(path.data)){
    dt0.bdnf <- fread(file.path(path.data,"SB_5HTTLPR_Patrick_Klaus.csv"))
}else{
    stop("Data relative to example 2 is under the Danish rules on data protection which do not allow to freely share person-sensitive health
care related data - contact the corresponding author of the paper for more information \n")
}
## * Analysis of the dataset

## ** data processing
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

## ** latent variable model
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
bdnf.null <- c('u~bdnf2mx', 'neo~httlpr2sx')
coef.bdnf.ML <- coef(e.bdnf2)


## * Simulation study (section 8.2) 

## convergence
dttype1.sim.factor[,100*min(n.rep)/20000, by = "n"]

## number of parameters
length(param.sim.factor)

## inflation of the type 1 error without correction
dttype1.sim.factor[n == 20 & link %in% c("Y2","Y4~eta","eta~Gene1Y","Y1~Gene2Y") & method %in% c("p.robustZtest"), .(link,type1 = type1, inflation = type1-0.05)]

## type 1 error after correction
dttype1.sim.factor[n == 20 & link %in% c("Y2","Y4~eta","eta~Gene1Y","Y1~Gene2Y") & method %in% c("p.robustKR"), .(link,type1)]

## * Illustration (section 9.2)

## number of patients/observations
c(nrow = NROW(dt.bdnf), n.id = length(unique(dt.bdnf$cimbi.id)))

## number of parameters
length(coef(e.bdnf2))

## genetic effect
summary(e.bdnf2)$coef[bdnf.null,]

## type 1 error from the simulation
dttype1.ill.factor[method == "p.robustZtest",.(link,type1)]

## corrected ML estimates
sCorrect(e.bdnf2) <- TRUE
coef.bdnf.MLc <- e.bdnf2$sCorrect$param
Mcompare <- cbind("ML-correct" = coef.bdnf.MLc,
                  "ML" = coef.bdnf.ML, 
                  "rdiff ML (%)" = 100*(coef.bdnf.MLc-coef.bdnf.ML)/abs(coef.bdnf.ML)
                  )
Mcompare

## genetic effect
summary2(e.bdnf2, robust = TRUE)$coef[bdnf.null,]

cbind(pvalue = summary(e.bdnf2, robust = TRUE)$coef[bdnf.null,"P-value"],
      pvalueC = summary2(e.bdnf2, robust = TRUE)$coef[bdnf.null,"P-value"],
      pc.increase =  100*(summary2(e.bdnf2, robust = TRUE)$coef[bdnf.null,"P-value"] / summary(e.bdnf2, robust = TRUE)$coef[bdnf.null,"P-value"]-1)
      )

## type 1 error found in the simulation
cbind(dttype1.ill.factor[method == "p.robustKR",.(link,type1MLc = type1)],
      dttype1.ill.factor[method == "p.robustZtest",.(typeML=type1)])



