## path <- "P:/Cluster/LVMproject/article-smallSampleInference"
## setwd(path)
## source("BATCH_illustration-lvm.R")

rm(list = ls())
gc()

## * seed
iter_sim <- as.numeric(Sys.getenv("SGE_TASK_ID"))
n.iter_sim <- as.numeric(Sys.getenv("SGE_TASK_LAST"))
if(is.na(iter_sim)){iter_sim <- 15}
if(is.na(n.iter_sim)){n.iter_sim <- 40}
cat("iteration ",iter_sim," over ",n.iter_sim,"\n",sep="")

set.seed(3)
seqSeed <- sample(1:max(1e5,n.iter_sim),size=n.iter_sim,replace=FALSE)
iSeed <- seqSeed[iter_sim]
set.seed(iSeed)

cat("seed: ",iSeed,"\n",sep="")

## * path
path <- "."
path.res <- file.path(path,"Results","illustration-lvm")
if(dir.exists(path.res)==FALSE){
    dir.create(path.res)
}
path.output <- file.path(path,"output","illustration-lvm")
if(dir.exists(path.output)==FALSE){
    dir.create(path.output)
}

## * libraries
library(lava)
library(data.table)
library(lavaSearch2)

## * settings
seqN <- c(24)
n.boot <- 0
n.rep <- 5e2 * 10

## * model

## ** generative model
m.generative <- lvm(frc ~ httlpr2 + u , amy ~ u, acc ~ u, hip ~ u)
regression(m.generative) <- c(pos.imm5, pos.shr, pos.del) ~ m.pos
regression(m.generative) <- c(neg.imm5, neg.shr, neg.del) ~ m.neg
regression(m.generative) <- c(neu.imm5, neu.shr, neu.del) ~ m.neu
regression(m.generative) <- u ~ age0
regression(m.generative) <- m.pos ~ u + age0
regression(m.generative) <- m.neg ~ u + age0
regression(m.generative) <- m.neu ~ u  + age0
covariance(m.generative) <- amy ~ hip
latent(m.generative) <- ~u + m.pos + m.neg + m.neu

distribution(m.generative, ~httlpr2) <- lava::binomial.lvm(p=0.4166)
distribution(m.generative, ~age0) <- lava::gaussian.lvm(sd=6.926)

generative.coef <- c("u" = 0.665862,
                     "amy" = -0.124351,
                     "acc" = -0.082397,
                     "hip" = 0.327872,
                     "pos.shr" = -5.860545,
                     "pos.del" = -6.363676,
                     "m.pos" = 9.640682,
                     "neg.shr" = -3.138381,
                     "neg.del" = -4.851905,
                     "m.neg" = 7.120127,
                     "neu.shr" = -1.187194,
                     "neu.del" = -2.133435,
                     "m.neu" = 9.846,
                     "frc~httlpr2" = 0.030707,
                     "u~age0" = -0.002664,
                     "amy~u" = 1.520811,
                     "acc~u" = 1.357993,
                     "hip~u" = 1.065282,
                     "pos.shr~m.pos" = 2.110027,
                     "pos.del~m.pos" = 2.284593,
                     "m.pos~u" = 0,
                     "m.pos~age0" = -0.06504,
                     "neg.shr~m.neg" = 1.731918,
                     "neg.del~m.neg" = 2.107359,
                     "m.neg~u" = 0,
                     "m.neg~age0" = -0.075826,
                     "neu.shr~m.neu" = 1.333535,
                     "neu.del~m.neu" = 1.494576,
                     "m.neu~u" = 0,
                     "m.neu~age0" = -0.057311,
                     "frc~~frc" = 0.000346,
                     "u~~u" = 0.004845,
                     "amy~~amy" = 0.021234,
                     "acc~~acc" = 0.001896,
                     "hip~~hip" = 0.009439,
                     "pos.imm5~~pos.imm5" = 0.517555,
                     "pos.shr~~pos.shr" = 0.913716,
                     "pos.del~~pos.del" = 0.32748,
                     "m.pos~~m.pos" = 0.158005,
                     "neg.imm5~~neg.imm5" = 0.31036,
                     "neg.shr~~neg.shr" = 0.81319,
                     "neg.del~~neg.del" = 0.156672,
                     "m.neg~~m.neg" = 0.393181,
                     "neu.imm5~~neu.imm5" = 0.304089,
                     "neu.shr~~neu.shr" = 0.725098,
                     "neu.del~~neu.del" = 0.562479,
                     "m.neu~~m.neu" = 0.401224,
                     "amy~~hip" = 0.010661)

## ** fit model
m.fit0 <- lvm(frc ~ u , amy ~ u, acc ~ u, hip ~ u)
regression(m.fit0) <- c(pos.imm5, pos.shr, pos.del) ~ m.pos
regression(m.fit0) <- c(neg.imm5, neg.shr, neg.del) ~ m.neg
regression(m.fit0) <- c(neu.imm5, neu.shr, neu.del) ~ m.neu
latent(m.fit0) <- ~u + m.pos + m.neg + m.neu

m.fit1 <- m.fit0
regression(m.fit1) <- u ~ age0
regression(m.fit1) <- frc ~ httlpr2
regression(m.fit1) <- m.pos ~ age0
regression(m.fit1) <- m.neg ~ age0
regression(m.fit1) <- m.neu ~ age0

m.fit2 <- m.fit1
covariance(m.fit2) <- amy ~ hip
regression(m.fit2) <- m.pos ~ u
regression(m.fit2) <- m.neg ~ u
regression(m.fit2) <- m.neu ~ u

## * simulation
out <- calibrateType1(m.fit2, warmup = list(m.fit1), true.coef = generative.coef,
                      param =  c("m.pos~u","m.neg~u","m.neu~u"),
                      n = seqN, n.rep = n.rep, 
                      generative.object = m.generative, generative.coef = generative.coef,
                      dir.save = path.res, label.file = iter_sim,
                      bootstrap = FALSE, seed = NULL, trace = 2, constrain = FALSE)

## * display
print(sessionInfo())
