path.results <- "./Results"

## * packages
library(data.table)

## * Import

## ** Simulation results
dttype1.ill.mm <- readRDS(file.path(path.results, "type1error-illustration-mixedModel.rds"))
dttype1.ill.factor <- readRDS(file.path(path.results, "type1error-illustration-factorModel.rds"))
dttype1.ill.lvm <- readRDS(file.path(path.results, "type1error-illustration-lvmModel.rds"))

## ** Real data
source("ANALYSIS.R")

## * Application A
treatment.null <- c("$k_1,k_2,k_3$")
M.tab1a <- matrix(NA, nrow = 1, ncol = 8,
                  dimnames = list(treatment.null,as.vector(sapply(c("ML","MLc"),paste,c("sigma","df","p-value","type1"),sep =".")))
                  )
df.tab1a <- as.data.frame(M.tab1a)

df.tab1a[treatment.null,"ML.df"] <- "$\\infty$"
df.tab1a[treatment.null,"ML.p-value"] <- Ftest.vitamin.ML$p.value
df.tab1a[treatment.null,"ML.type1"] <- dttype1.ill.mm[link == "global" & method == "p.Ztest",type1]

df.tab1a[treatment.null,"MLc.df"] <- Ftest.vitamin.MLc$parameter
df.tab1a[treatment.null,"MLc.p-value"] <- Ftest.vitamin.MLc$p.value
df.tab1a[treatment.null,"MLc.type1"] <- dttype1.ill.mm[link == "global" & method == "p.KR",type1]

df.tab1a.print <- df.tab1a
df.tab1a.print[,c("ML.sigma","MLc.sigma")] <- ""
df.tab1a.print[,c("MLc.df")] <- round(df.tab1a[,c("MLc.df")],1)
df.tab1a.print[,c("ML.p-value")] <- paste0("\\num{",formatC(df.tab1a[,c("ML.p-value")], format = "e", digits = 2),"}")
df.tab1a.print[,c("MLc.p-value")] <- paste0("\\num{",formatC(df.tab1a[,c("MLc.p-value")], format = "e", digits = 2),"}")
df.tab1a.print[,c("ML.type1","MLc.type1")] <- sapply(round(100*df.tab1a[,c("ML.type1","MLc.type1")],2),paste0,"\\%")

df.tab1a.print <- cbind(cbind(treatment.null,df.tab1a.print))
names(df.tab1a.print) <- c("parameter",
                           "\\(\\sigma\\)","\\(\\df\\)","p-value","type 1 error",
                           "\\(\\sigma\\)","\\(\\df\\)","p-value","type 1 error")

## * Application B
M.tab1b <- matrix(NA, nrow = 2, ncol = 8,
                dimnames = list(bdnf.null,as.vector(sapply(c("ML","MLc"),paste,c("sigma","df","p-value","type1"),sep =".")))
                )
df.tab1b <- as.data.frame(M.tab1b)
df.tab1b[,c("ML.sigma","ML.p-value")] <- summary(e.bdnf2, robust = TRUE)$coef[bdnf.null,c("Robust SE","P-value")]
df.tab1b[,c("ML.df")] <- "$\\infty$"
df.tab1b["u~bdnf2mx","ML.type1"] <- dttype1.ill.factor[link == "u~bdnf2" & method == "p.robustZtest",.(typeML=type1)]  
df.tab1b["neo~httlpr2sx","ML.type1"] <- dttype1.ill.factor[link == "neo~httlpr2" & method == "p.robustZtest",.(typeML=type1)]  

df.tab1b[,c("MLc.sigma","MLc.df","MLc.p-value")] <- summary2(e.bdnf2, robust = TRUE)$coef[bdnf.null,c("robust SE","df","P-value")]
df.tab1b["u~bdnf2mx","MLc.type1"] <- dttype1.ill.factor[link == "u~bdnf2" & method == "p.robustKR",.(typeML=type1)]  
df.tab1b["neo~httlpr2sx","MLc.type1"] <- dttype1.ill.factor[link == "neo~httlpr2" & method == "p.robustKR",.(typeML=type1)]  

df.tab1b.print <- df.tab1b
df.tab1b.print[,c("ML.sigma","MLc.sigma")] <- round(df.tab1b[,c("ML.sigma","MLc.sigma")],3)
df.tab1b.print[,c("MLc.df")] <- round(df.tab1b[,c("MLc.df")],1)
df.tab1b.print[,c("ML.p-value")] <- paste0("\\num{",formatC(df.tab1b[,c("ML.p-value")], format = "e", digits = 2),"}")
df.tab1b.print[,c("MLc.p-value")] <- paste0("\\num{",formatC(df.tab1b[,c("MLc.p-value")], format = "e", digits = 2),"}")
df.tab1b.print[,c("ML.type1","MLc.type1")] <- sapply(round(100*df.tab1b[,c("ML.type1","MLc.type1")],2),paste0,"\\%")

df.tab1b.print <- cbind(cbind(names(bdnf.null),df.tab1b.print))
names(df.tab1b.print) <- c("parameter",
                           "\\(\\sigma\\)","\\(\\df\\)","p-value","type 1 error",
                           "\\(\\sigma\\)","\\(\\df\\)","p-value","type 1 error")

## * Application C
M.tab1c <- matrix(NA, nrow = 3, ncol = 8,
                  dimnames = list(memory.null,as.vector(sapply(c("ML","MLc"),paste,c("sigma","df","p-value","type1"),sep =".")))
                  )
df.tab1c <- as.data.frame(M.tab1c)
df.tab1c[,c("ML.sigma","ML.p-value")] <- summary(e.memory)$coef[memory.null,c("Std. Error","P-value")]
df.tab1c[,c("ML.df")] <- "$\\infty$"
df.tab1c[memory.null[1],"ML.type1"] <- dttype1.ill.lvm[link == memory.null[1] & method == "p.Ztest",.(typeML=type1)]  
df.tab1c[memory.null[2],"ML.type1"] <- dttype1.ill.lvm[link == memory.null[2] & method == "p.Ztest",.(typeML=type1)]  
df.tab1c[memory.null[3],"ML.type1"] <- dttype1.ill.lvm[link == memory.null[3] & method == "p.Ztest",.(typeML=type1)]  

df.tab1c[,c("MLc.sigma","MLc.df","MLc.p-value")] <- summary2(e.memory)$coef[memory.null,c("Std. Error","df","P-value")]
df.tab1c[memory.null[1],"MLc.type1"] <- dttype1.ill.lvm[link == memory.null[1] & method == "p.KR",.(typeML=type1)]  
df.tab1c[memory.null[2],"MLc.type1"] <- dttype1.ill.lvm[link == memory.null[2] & method == "p.KR",.(typeML=type1)]  
df.tab1c[memory.null[3],"MLc.type1"] <- dttype1.ill.lvm[link == memory.null[3] & method == "p.KR",.(typeML=type1)]  

df.tab1c.print <- df.tab1c
df.tab1c.print[,c("ML.sigma","MLc.sigma")] <- round(df.tab1c[,c("ML.sigma","MLc.sigma")],3)
df.tab1c.print[,c("MLc.df")] <- round(df.tab1c[,c("MLc.df")],1)
df.tab1c.print[,c("ML.p-value")] <- paste0("\\num{",formatC(df.tab1c[,c("ML.p-value")], format = "e", digits = 2),"}")
df.tab1c.print[,c("MLc.p-value")] <- paste0("\\num{",formatC(df.tab1c[,c("MLc.p-value")], format = "e", digits = 2),"}")
df.tab1c.print[,c("ML.type1","MLc.type1")] <- sapply(round(100*df.tab1c[,c("ML.type1","MLc.type1")],2),paste0,"\\%")

df.tab1c.print <- cbind(cbind(names(memory.null),df.tab1c.print))
names(df.tab1c.print) <- c("parameter",
                           "\\(\\sigma\\)","\\(\\df\\)","p-value","type 1 error",
                           "\\(\\sigma\\)","\\(\\df\\)","p-value","type 1 error")

## * Table 1
df.tab1.print <- rbind(cbind(Application = "A",df.tab1a.print,"x" = "\\\\"),
                       cbind(Application = "B",df.tab1b.print,"x" = "\\\\"),
                       cbind(Application = "C",df.tab1c.print,"x" = "\\\\"))
df.tab1.print$x <- as.character(df.tab1.print$x)
df.tab1.print$Application <- as.character(df.tab1.print$Application)
df.tab1.print$Application[duplicated(df.tab1.print$Application)] <- " "
df.tab1.print[c(1,3),"x"] <- paste0(df.tab1.print[c(1,3),"x"]," [0.5cm]")

cat("\n",apply(df.tab1.print,1,function(x){paste(paste(x, collapse =  " & "),"\n")}))
 ## A & $k_1,k_2,k_3$ &  & $\infty$ & \num{1.76e-03} & 10.22\% &  & 49.1 & \num{1.23e-02} & 4.49\% & \\ [0.5cm] 
 ## B & $\gamma_2$ & 0.026 & $\infty$ & \num{5.11e-03} & 7.41\% & 0.027 & 65.4 & \num{9.03e-03} & 5.59\% & \\ 
 ##   & $k_1$ & 0.016 & $\infty$ & \num{7.23e-06} & 6.14\% & 0.017 & 67.5 & \num{4.80e-05} & 4.78\% & \\ [0.5cm] 
 ## C & $b_1$ & 2.106 & $\infty$ & \num{5.45e-04} & 6.32\% & 2.143 & 13.5 & \num{4.53e-03} & 2.02\% & \\ 
 ##   & $b_2$ & 2.355 & $\infty$ & \num{4.21e-03} & 8.5\% & 2.429 & 10.0 & \num{1.96e-02} & 3.71\% & \\ 
 ##   & $b_3$ & 2.041 & $\infty$ & \num{7.20e-02} & 8.41\% & 2.13 &  7.7 & \num{1.24e-01} & 3.42\% & \\ 
