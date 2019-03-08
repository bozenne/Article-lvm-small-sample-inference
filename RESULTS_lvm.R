## path <- "P:/Cluster/LVMproject/article-smallSampleInference"
## setwd(path)

path.results <- "./Results"
path.data <- "./data2/"

## * packages
devtools::load_all("lavaSearch2")
library(data.table)
library(foreign)

## * Import
## ** simulation results
dttype1.sim.lvm <- readRDS(file.path(path.results, "type1error-simulation-lvmModel.rds"))
dtbias.sim.lvm <- readRDS(file.path(path.results, "bias-simulation-lvmModel.rds"))
dttype1.ill.lvm <- readRDS(file.path(path.results, "type1error-illustration-lvmModel.rds"))
param.sim.lvm <- readRDS(file.path(path.results, "param-simulation-lvmModel.rds"))

## ** real data
if(dir.exists(path.data)){
    dt0.memory <- as.data.table(read.spss(file.path(path.data,"SPSSmedscan_v2.sav"), to.data.frame = TRUE))
    dt0.gene <- fread(file.path(path.data,"SB_5HTTLPR_Patrick_Klaus.csv"))
}else{ 
    stop("Data relative to example 3 is under the Danish rules on data protection which do not allow to freely share person-sensitive health
care related data - contact the corresponding author of the paper for more information \n")
}


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
memory.null <- c( c('m.pos~u', 'm.neu~u', 'm.neg~u'))
coef.memory.ML <- coef(e.memory)

## * Simulation study (section 7.3)

## number of parameters in the model
length(param.sim.lvm)

## type 1 error 
keep.link <- c("Y2", "Y4~eta1","Y1~Gene2Y", "eta1~Gene1Y","eta1~eta2","Y1~~Y2")
cbind(dttype1.sim.lvm[n==20 & link %in% keep.link & method %in% c("p.Ztest"), .(link,type1ML = type1)],
      dttype1.sim.lvm[n==20 & link %in% keep.link & method %in% c("p.SSC"),.(type1MLssc = type1)],
      dttype1.sim.lvm[n==20 & link %in% keep.link & method %in% c("p.KR"),.(type1MLc = type1)])

## * Illustration (section 8.3)

## number of patients/observations
c(nrow = NROW(dtG.memory), id = length(unique(dtG.memory$cimbi.id)))

## loadings
summary(e.memory)$coef[memory.null,]

## type 1 error from the simulation
dttype1.ill.lvm[method == "p.Ztest",.(link,type1)][match(link,memory.null)]

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
      dttype1.ill.lvm[method == "p.Ztest",.(typeML=type1)])[match(link,memory.null)]

