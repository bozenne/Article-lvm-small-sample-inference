## path <- "P:/Cluster/LVMproject/article-smallSampleInference"
## setwd(path)

path.results <- "./Results"
path.data <- "./data2/"

## * packages
devtools::load_all("lavaSearch2")
library(data.table)
library(ggplot2)
library(foreign)

## * Import
## ** simulation results
dttype1.sim.lvm <- readRDS(file.path(path.results, "type1error-simulation-lvmModel.rds"))
dtbias.sim.lvm <- readRDS(file.path(path.results, "bias-simulation-lvmModel.rds"))
dttype1.ill.lvm <- readRDS(file.path(path.results, "type1error-illustration-lvmModel.rds"))

## m.generative.lvm <- lvm(c(Y1,Y2,Y3,Y5) ~ eta1,
                        ## Y4 ~ 1,
                        ## Age ~ 1,
                        ## c(Z1,Z2,Z3,Z4,Z5) ~ eta2,
                        ## eta1 ~ G1,
                        ## eta2 ~ G2)
## latent(m.generative.lvm) <- ~eta1+eta2

## m.fit.lvm <- lvm(c(Y1,Y2,Y3,Y4,Y5) ~ eta1,
                 ## c(Z1,Z2,Z3,Z4,Z5) ~ eta2,
                 ## Y1 ~ Age,
                 ## eta1 ~ G1 + Age,
                 ## eta2 ~ G2,
                 ## eta1 ~ eta2)
## covariance(m.fit.lvm) <- Y1 ~ Y2 
## latent(m.fit.lvm) <- ~eta1+eta2

## ** real data

dt0.memory <- as.data.table(read.spss(file.path(path.data,"SPSSmedscan_v2.sav"), to.data.frame = TRUE))
dt0.gene <- fread(file.path(path.data,"SB_5HTTLPR_Patrick_Klaus.csv"))

## * Analysis of the dataset

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

## estimation
e.memory0 <- estimate(m.memory0, data = dtG.memory)

e.memory <- estimate(m.memory, data = dtG.memory,
                     control = list(constrain = TRUE, start = coef(e.memory0)))
## summary(e.memory)

## coefficients
memory.null <- c( c('m.pos~u', 'm.neg~u', 'm.neu~u'))
coef.memory.ML <- coef(e.memory)

## * ML estimates (section 2.1.3)

## number of patients/observations
c(nrow = NROW(dtG.memory), id = length(unique(dtG.memory$cimbi.id)))

## loadings
summary(e.memory)$coef[memory.null,]

## type 1 error from the simulation
dttype1.ill.lvm[method == "p.Ztest",.(link,type1)]

## * corrected type 1 error in the simulation study (section 7.3)

## chunk 57
e.true.lvm <- lava::estimate(m.fit.lvm, lava::sim(m.generative.lvm, n = 1e5))
length(coef(e.true.lvm))

## chunk 58
cbind(dtSS.coverage.lvm[n==20 & method %in% c("p.Ztest"), .(link,type1ML = type1)],
      dtSS.coverage.lvm[n==20 & method %in% c("p.SSC"),.(type1MLssc = type1)],
      dtSS.coverage.lvm[n==20 & method %in% c("p.KR"),.(type1MLc = type1)])

## * corrected ML estimates (section 8.3)

## corrected ML estimates
sCorrect(e.memory) <- TRUE
coef.memory.MLc <- e.memory$sCorrect$param
coef.type <- as.data.table(coefType(e.memory, as.lava = FALSE))
vec.type <- coef.type[!is.na(lava),setNames(detail,name)]

Mcompare <- data.table("link" = names(coef.memory.MLc),
                       "type" = vec.type[names(coef.memory.MLc)],
                       "ML-correct" = coef.memory.MLc,
                       "ML" = coef.memory.ML, 
                       "rdiff ML (%)" = 100*(coef.memory.MLc-coef.memory.ML)/abs(coef.memory.ML)
                       )
Mcompare[,.(min = min(`rdiff ML (%)`), max = max(`rdiff ML (%)`)),by=type]

## loadings
summary2(e.memory)$coef[memory.null,"P-value"]

cbind(pvalue = summary(e.memory)$coef[memory.null,"P-value"],
      pvalueC = summary2(e.memory)$coef[memory.null,"P-value"],
      pc.increase =  100*(summary2(e.memory)$coef[memory.null,"P-value"] / summary(e.memory)$coef[memory.null,"P-value"]-1)
      )

## type 1 error found in the simulation
cbind(dttype1.ill.lvm[method == "p.KR",.(link,type1MLc = type1)],
      dttype1.ill.lvm[method == "p.Ztest",.(typeML=type1)])

## * Appendix E - table 3

## chunk 64
table3 <- createTable(dtSS.bias.lvm, 
                      seqN = c(20,30,50,100), 
                      digit = 3, 
                      seqType = c("Sigma_var","Sigma_cov","Psi_var"),
                      convert2latex = FALSE)
table3
