## path <- "P:/Cluster/LVMproject/article-smallSampleInference"
## setwd(path)
## source("BATCH_illustration-factorModel.R")

rm(list = ls())
gc()

## * seed
iter_sim <- as.numeric(Sys.getenv("SGE_TASK_ID"))
n.iter_sim <- as.numeric(Sys.getenv("SGE_TASK_LAST"))
if(is.na(iter_sim)){iter_sim <- 14}
if(is.na(n.iter_sim)){n.iter_sim <- 40}
cat("iteration ",iter_sim," over ",n.iter_sim,"\n",sep="")

set.seed(3)
seqSeed <- sample(1:max(1e5,n.iter_sim),size=n.iter_sim,replace=FALSE)
iSeed <- seqSeed[iter_sim]
set.seed(iSeed)

cat("seed: ",iSeed,"\n",sep="")

## * path
path <- "."
path.res <- file.path(path,"Results","illustration-factorModel")
if(dir.exists(path.res)==FALSE){
    dir.create(path.res)
}
path.output <- file.path(path,"output","illustration-factorModel")
if(dir.exists(path.output)==FALSE){
    dir.create(path.output)
}

## * libraries
library(lava)
library(data.table)
devtools::load_all("lavaSearch2")

## * settings
seqN <- c(73)
n.boot <- 0
n.rep <- 5e2 * 3

## * model

## ** generative model
m.generative <- lvm(neo ~ u + httlpr2 + scanner,
                    cau ~ u + age0 + scanner,
                    put ~ u + scanner,
                    hip ~ u + scanner,
                    amy ~ u + scanner,
                    u ~ age0 + gender + bdnf2 + bmi0 + inj0)
latent(m.generative) <- ~u
covariance(m.generative) <- cau ~ put
covariance(m.generative) <- amy ~ hip
distribution(m.generative, ~scanner) <- lava::binomial.lvm(p=0.6575)
distribution(m.generative, ~gender) <- lava::binomial.lvm(p=0.7671)
distribution(m.generative, ~bdnf2) <- lava::binomial.lvm(p=0.5753)
distribution(m.generative, ~httlpr2) <- lava::binomial.lvm(p=0.5890)
distribution(m.generative, ~age0) <- lava::gaussian.lvm(sd=16.09)
distribution(m.generative, ~inj0) <- lava::gaussian.lvm(sd=2.479)
distribution(m.generative, ~bmi0) <- lava::gaussian.lvm(sd=0.02276)
transform(m.generative, cimbi.id~u) <- function(x){
    vec.id <- c(10753, 11074, 11131, 30006, 30009, 30012, 30018, 30027, 30030, 30036, 30045, 30054, 30060, 30060, 30066, 30075, 30075, 30081, 30081, 30087, 30087, 30114, 30114, 30120, 30126, 30141, 30156, 30159, 50054, 50062, 50066, 50068, 50151, 50228, 50235, 50246, 50249, 50251, 50285, 50299, 50313, 50315, 50334, 50335, 50336, 50337, 50469, 50524, 50528, 50532, 50542, 50547, 50548, 50562, 50677, 50678, 50689, 50701, 50763, 50804, 50819, 51081, 51325, 51742, 51745, 51756, 51777, 51906, 51907, 51911, 51915, 52325, 52334)
    return(vec.id[1:NROW(x)])
}

## m.generative <- lvm(c(neo, cau, amy, put, hip)~u)
## m.generative <- regression(m.generative, to = endogenous(m.generative), from = "scanner")
## m.generative <- regression(m.generative, u ~ age0 + gender + bdnf2 + bmi0 + inj0)
## latent(m.generative) <- ~u
## covariance(m.generative) <- cau ~ put
## covariance(m.generative) <- amy ~ hip
## regression(m.generative) <- neo ~ httlpr2
## distribution(m.generative, ~scanner+gender+bdnf2+httlpr2) <- lava::binomial.lvm()
## butils::object2script(dt.bdnf$cimbi.id)
generative.coef <- c("u" = -1.015588,
                     "cau" = 1.17927,
                     "put" = 1.466333,
                     "hip" = 0.357585,
                     "amy" = 0.446344,
                     "neo~httlpr2" = 0,
                     "neo~scanner" = 0.49859,
                     "u~age0" = 0,
                     "u~bmi0" = 0,
                     "u~inj0" = -2.116476,
                     "u~gender" = 0,
                     "u~bdnf2" = 0,
                     "cau~u" = 0.392194,
                     "cau~age0" = -0.004079,
                     "cau~scanner" = 0.401596,
                     "put~u" = 0.635584,
                     "put~scanner" = 0.347252,
                     "hip~u" = 0.599209,
                     "hip~scanner" = 0.194701,
                     "amy~u" = 0.717358,
                     "amy~scanner" = 0.0362,
                     "neo~~neo" = 0.000326,
                     "u~~u" = 0.010595,
                     "cau~~cau" = 0.011259,
                     "put~~put" = 0.004755,
                     "hip~~hip" = 0.011645,
                     "amy~~amy" = 0.024046,
                     "cau~~put" = 0.00403,
                     "hip~~amy" = 0.010866)

# plot(m.generative)

## ** fit model
## m.fit0 <- lvm(neo ~ u + httlpr2 + scanner,
##              cau ~ u + age0 + scanner,
##              put ~ u + scanner,
##              hip ~ u + scanner,
##              amy ~ u + scanner,
##              u ~ age0 + gender + bdnf2 + bmi0 + inj0)
## latent(m.fit0) <- ~u

## m.fit1 <- m.fit0
## covariance(m.fit1) <- cau ~ put
## covariance(m.fit1) <- amy ~ hip

m.fit0 <- lvm(neo ~ u + httlpr2 + scanner,
              cau ~ u + scanner,
              put ~ u + scanner,
              hip ~ u + scanner,
              amy ~ u + scanner,
              u ~ age0 + gender + bdnf2 + bmi0 + inj0)
latent(m.fit0) <- ~u

m.fit1 <- m.fit0
covariance(m.fit1) <- cau ~ put
covariance(m.fit1) <- amy ~ hip

m.fit2 <- m.fit1
regression(m.fit2) <- cau ~ u + age0

## lava::estimate(m.fit, lava::sim(m.generative, n = 1e4, p = generative.coef))
## aa1bbbbabbbbbbbaaabba1

## * simulation
out <- calibrateType1(m.fit2, warmup = list(m.fit0, m.fit1) ,cluster = "cimbi.id", true.coef = generative.coef,
                      param =  c('u~bdnf2', 'neo~httlpr2', 'u~age0', 'u~bmi0', 'u~gender'),
                      n = seqN, n.rep = n.rep, 
                      generative.object = m.generative, generative.coef = generative.coef,
                      dir.save = path.res, label.file = iter_sim,
                      bootstrap = FALSE, seed = NULL, trace = 2, control = list(constrain = TRUE))

## * display
print(sessionInfo())
 

if(FALSE){
    ## ** generative model
m.generative <- lvm(neo ~ u + httlpr2 + scanner,
                    cau ~ u + age0 + scanner,
                    put ~ u + scanner,
                    hip ~ u + scanner,
                    amy ~ u + scanner,
                    u ~ age0 + gender + bdnf2 + bmi0 + inj0)
latent(m.generative) <- ~u
covariance(m.generative) <- cau ~ put
covariance(m.generative) <- amy ~ hip
distribution(m.generative, ~scanner) <- lava::binomial.lvm(p=0.6575)
distribution(m.generative, ~gender) <- lava::binomial.lvm(p=0.7671)
distribution(m.generative, ~bdnf2) <- lava::binomial.lvm(p=0.5753)
distribution(m.generative, ~httlpr2) <- lava::binomial.lvm(p=0.5890)
distribution(m.generative, ~age0) <- lava::gaussian.lvm(sd=16.09)
distribution(m.generative, ~inj0) <- lava::gaussian.lvm(sd=2.479)
distribution(m.generative, ~bmi0) <- lava::gaussian.lvm(sd=0.02276)
transform(m.generative, cimbi.id~u) <- function(x){
    vec.id <- c(10753, 11074, 11131, 30006, 30009, 30012, 30018, 30027, 30030, 30036, 30045, 30054, 30060, 30060, 30066, 30075, 30075, 30081, 30081, 30087, 30087, 30114, 30114, 30120, 30126, 30141, 30156, 30159, 50054, 50062, 50066, 50068, 50151, 50228, 50235, 50246, 50249, 50251, 50285, 50299, 50313, 50315, 50334, 50335, 50336, 50337, 50469, 50524, 50528, 50532, 50542, 50547, 50548, 50562, 50677, 50678, 50689, 50701, 50763, 50804, 50819, 51081, 51325, 51742, 51745, 51756, 51777, 51906, 51907, 51911, 51915, 52325, 52334)
    return(vec.id[1:NROW(x)])
}

## m.generative <- lvm(c(neo, cau, amy, put, hip)~u)
## m.generative <- regression(m.generative, to = endogenous(m.generative), from = "scanner")
## m.generative <- regression(m.generative, u ~ age0 + gender + bdnf2 + bmi0 + inj0)
## latent(m.generative) <- ~u
## covariance(m.generative) <- cau ~ put
## covariance(m.generative) <- amy ~ hip
## regression(m.generative) <- neo ~ httlpr2
## distribution(m.generative, ~scanner+gender+bdnf2+httlpr2) <- lava::binomial.lvm()
## butils::object2script(dt.bdnf$cimbi.id)
generative.coef <- c("u" = -1.015588,
                     "cau" = 1.17927,
                     "put" = 1.466333,
                     "hip" = 0.357585,
                     "amy" = 0.446344,
                     "neo~httlpr2" = 0,
                     "neo~scanner" = 0.49859,
                     "u~age0" = 0,
                     "u~bmi0" = 0,
                     "u~inj0" = -2.116476,
                     "u~gender" = 0,
                     "u~bdnf2" = 0,
                     "cau~u" = 0.392194,
                     "cau~age0" = -0.004079,
                     "cau~scanner" = 0.401596,
                     "put~u" = 0.635584,
                     "put~scanner" = 0.347252,
                     "hip~u" = 0.599209,
                     "hip~scanner" = 0.194701,
                     "amy~u" = 0.717358,
                     "amy~scanner" = 0.0362,
                     "neo~~neo" = 0.000326,
                     "u~~u" = 0.010595,
                     "cau~~cau" = 0.011259,
                     "put~~put" = 0.004755,
                     "hip~~hip" = 0.011645,
                     "amy~~amy" = 0.024046,
                     "cau~~put" = 0.00403,
                     "hip~~amy" = 0.010866)

# plot(m.generative)

    ## ** fit model
    ## m.fit0 <- lvm(neo ~ u + httlpr2 + scanner,
    ##              cau ~ u + age0 + scanner,
    ##              put ~ u + scanner,
    ##              hip ~ u + scanner,
    ##              amy ~ u + scanner,
    ##              u ~ age0 + gender + bdnf2 + bmi0 + inj0)
    ## latent(m.fit0) <- ~u

    ## m.fit1 <- m.fit0
    ## covariance(m.fit1) <- cau ~ put
    ## covariance(m.fit1) <- amy ~ hip

    m.fit0 <- lvm(neo ~ u + httlpr2 + scanner,
                  cau ~ u + scanner,
                  put ~ u + scanner,
                  hip ~ u + scanner,
                  amy ~ u + scanner,
                  u ~ age0 + gender + bdnf2 + bmi0 + inj0)
    latent(m.fit0) <- ~u

    m.fit1 <- m.fit0
    covariance(m.fit1) <- cau ~ put
    covariance(m.fit1) <- amy ~ hip

    m.fit2 <- m.fit1
    regression(m.fit2) <- cau ~ u + age0

    ## lava::estimate(m.fit, lava::sim(m.generative, n = 1e4, p = generative.coef))
    ## aa1bbbbabbbbbbbaaabba1

    }
