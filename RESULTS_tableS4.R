path.data2 <- "./data2/"

## * package
library(data.table)
source("FCT.R") ## get function fitStudent

if(FALSE){ ## sanity check 
    set.seed(10)
    dt.test <- rbind(data.table(name = "rt", estimate.MLcorrected = rt(1e3, df = 3), se.MLcorrected = 1, df.MLcorrected = 3),
                     data.table(name = "rnorm", estimate.MLcorrected = rnorm(1e3), se.MLcorrected = 1, df.MLcorrected = Inf),
                     data.table(name = "runif", estimate.MLcorrected = runif(1e3), se.MLcorrected = 1, df.MLcorrected = NA))
    fitStudent(dt.test, robust = FALSE)
    ##       se.empirical df.empirical mean.se mean.df     Etype1  zdist.ksD    zdist.ksP  tdist.ksD    tdist.ksP empirical.cv
    ## rt       1.0243040 3.061981e+00       1       3 0.05157778 0.07855195 8.739392e-06 0.03930024 9.108476e-02            1
    ## rnorm    0.9802015 4.408587e+01       1     Inf 0.05173461 0.02132348 7.534688e-01 0.01818344 8.955425e-01            1
    ## runif    0.2918500 2.455879e+07       1     NaN        NaN 0.07207078 6.157547e-05 0.07218601 5.956193e-05            1

}

## * Import
## ** simulation results
if(dir.exists(path.data2)){
    dist.simulation <- readRDS(file.path("data2","dist-simulation.rds"))
    ## saveRDS(lapply(dist.simulation, function(x){x[n==20]}), file = file.path("data2","dist-simulation.rds"))
}else{ 
    stop("Data is too large to be stored on Github (30Mb). Contact the corresponding author to obtain the dataset. \n")
}

## * fit distributions
checkStud.MM <- fitStudent(dist.simulation$MM[n==20],
                           robust = FALSE,
                           seqLink = c("Y2","eta~Gene1Y"))

checkStud.robustMM <- fitStudent(dist.simulation$MM[n==20],
                                 robust = TRUE,
                                 seqLink = c("Y2","eta~Gene1Y"))

checkStud.factor <- fitStudent(dist.simulation$factor[n==20],
                               robust = FALSE,
                               seqLink = c("Y2","Y1~Gene2Y","eta~Gene1Y","Y4~eta"))

checkStud.robustfactor <- fitStudent(dist.simulation$factor[n==20],
                                     robust = TRUE, 
                                     seqLink = c("Y2","Y1~Gene2Y","eta~Gene1Y","Y4~eta"))
 
checkStud.lvm <- fitStudent(dist.simulation$lvm[n==20],
                            robust = FALSE, value.max = 10,
                            seqLink = c("Y2","Y1~Gene2Y","eta1~Gene1Y","Y4~eta1","eta1~eta2","Y1~~Y2"))

checkStud.robustlvm <- fitStudent(dist.simulation$lvm[n==20],
                                  robust = TRUE, value.max = 10,
                                  seqLink = c("Y2","Y1~Gene2Y","eta1~Gene1Y","Y4~eta1","eta1~eta2","Y1~~Y2"))

dt_gghist <- dist.simulation$lvm[n==20 & name %in% c("eta1~eta2","Y1~~Y2"), .(wald = na.omit(estimate.MLcorrected/se.MLcorrected)), by = "name"]

## gg_hist <- ggplot(dt_gghist, aes(x = wald))
## gg_hist <- gg_hist + geom_histogram(bins = 100, aes(y=..density..), color = "red") + facet_wrap(~name)
## gg_hist <- gg_hist + stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = 1.1), size = 2)


rm.col <- c("empirical.cv","zdist.ksD","zdist.ksP")
keep.col <- setdiff(colnames(checkStud.MM),rm.col)

greek <- list(a = c("$\\nu_2$","$\\gamma_1$"),
              b = c("$\\nu_2$","$\\lambda_4$","$\\gamma_1$","$k_1$"),
              c = c("$\\nu_2$","$k_1$","$\\lambda_4$","$\\gamma_1$","$b_1$","$\\sigma_{12}$"))
row <- list(a = c("Y2","eta~Gene1Y"),
            b = c("Y2","Y4~eta","eta~Gene1Y","Y1~Gene2Y"),
            c = c("Y2","Y1~Gene2Y","Y4~eta1","eta1~Gene1Y","eta1~eta2","Y1~~Y2"))

## * Table S4
df.table1 <- rbind(cbind(scenario = "A", parameter = greek$a, as.data.frame(checkStud.MM[row$a,keep.col])),
                   cbind(scenario = "B", parameter = greek$b, as.data.frame(checkStud.factor[row$b,keep.col])),
                   cbind(scenario = "C", parameter = greek$c, as.data.frame(checkStud.lvm[row$c,keep.col]))
                   )

df.table1$scenario <- as.character(df.table1$scenario)
rownames(df.table1) <- NULL
df.table1$se.empirical <- formatC(round(df.table1$se.empirical, digits = 3), format = "f", digits = 3)
df.table1$df.empirical <- formatC(round(df.table1$df.empirical, digits = 1), format = "f", digits = 1)
df.table1$mean.se <- formatC(round(df.table1$mean.se, digits = 3), format = "f", digits = 0)
df.table1$mean.df <- formatC(round(df.table1$mean.df, digits = 1), format = "f", digits = 1)
df.table1$Etype1 <- formatC(round(df.table1$Etype1, digits = 3), format = "f", digits = 3)
df.table1$tdist.ksD <- formatC(round(df.table1$tdist.ksD, digits = 3), format = "f", digits = 3)
df.table1$tdist.ksP <- format.pval(df.table1$tdist.ksP, digits = 3, eps = 1e-3)
df.table1[duplicated(df.table1$scenario),"scenario"] <- ""

addtorow <- list()
addtorow$pos <- list(0,0,2,6)
addtorow$command <- c("&&\\multicolumn{2}{c}{Empirical Student} & \\multicolumn{2}{c}{Modeled Student} & Expected
& \\multicolumn{2}{c}{KS-test} \\\\ \\cmidrule(lr){3-4} \\cmidrule(lr){5-6} \\cmidrule(lr){8-9} \n",
"Scenario & parameter & dispersion & df & dispersion & df & type 1 error & statistic & p-value  \\\\\n",
"[4mm] ","[4mm] ")
print(xtable::xtable(df.table1,
                     label = "tab:validation",
                     caption = paste("Comparison between the empirical distribution of the Wald statistic vs. the modeled distribution (after correction) for n=20. ",
                                     "Empirical Student: standard error and degrees of freedom of a non-standardized Student's t-distribution fitted to the empirical distribution.",
                                     "Modeled Student: average estimated degrees of freedom over the simulations.",
                                     "Expected type 1 error: rejection rate under the empirical Student when using the 2.5 and 97.5\\% quantile of the modeled Student.",
                                     "Kolmogorov Smirnov test: comparison between the empirical cumluative distribution function (cdf) and the cdf of the empirical Student."),
                     ),
      add.to.row =  addtorow,
      include.colnames = FALSE,
      include.rownames=FALSE,
      sanitize.text.function = function(x){x},
      booktabs = TRUE)

## \begin{table}[ht]
## \centering
## \begin{tabular}{lllllllll}
##   \toprule
##   &&\multicolumn{2}{c}{Empirical Student} & \multicolumn{2}{c}{Modeled Student} & Expected
## & \multicolumn{2}{c}{KS-test} \\ \cmidrule(lr){3-4} \cmidrule(lr){5-6} \cmidrule(lr){8-9} 
##  Scenario & parameter & dispersion & df & dispersion & df & type 1 error & statistic & p-value  \\
##  \midrule
## A & $\nu_2$ & 1.006 & 40.8 & 1 & 36.7 & 0.051 & 0.005 & 0.659 \\ 
##    & $\gamma_1$ & 1.004 & 17.9 & 1 & 18.3 & 0.051 & 0.004 & 0.961 \\ 
##    [4mm] B & $\nu_2$ & 0.998 & 39.2 & 1 & 9.1 & 0.029 & 0.006 & 0.513 \\ 
##    & $\lambda_4$ & 1.067 & 117.2 & 1 & 5.9 & 0.023 & 0.004 & 0.885 \\ 
##    & $\gamma_1$ & 1.088 & 311.8 & 1 & 10.5 & 0.043 & 0.006 & 0.455 \\ 
##    & $k_1$ & 1.063 & 14.8 & 1 & 17.1 & 0.066 & 0.005 & 0.792 \\ 
##    [4mm] C & $\nu_2$ & 1.026 & 45.5 & 1 & 10.1 & 0.035 & 0.006 & 0.536 \\ 
##    & $k_1$ & 1.106 & 12.6 & 1 & 17.1 & 0.080 & 0.004 & 0.936 \\ 
##    & $\lambda_4$ & 1.080 & 145.9 & 1 & 6.0 & 0.025 & 0.008 & 0.237 \\ 
##    & $\gamma_1$ & 1.152 & 49.6 & 1 & 10.1 & 0.059 & 0.009 & 0.135 \\ 
##    & $b_1$ & 1.127 & 783296.7 & 1 & 3.6 & 0.010 & 0.017 & <0.001 \\ 
##    & $\sigma_{12}$ & 1.074 & 934003.5 & 1 & 7.8 & 0.031 & 0.038 & <0.001 \\ 
##    \bottomrule
## \end{tabular}
## \caption{Comparison between the empirical distribution of the Wald statistic vs. the modeled distribution (after correction) for n=20.  Empirical Student: standard error and degrees of freedom of a non-standardized Student's t-distribution fitted to the empirical distribution. Modeled Student: average estimated degrees of freedom over the simulations. Expected type 1 error: rejection rate under the empirical Student when using the 2.5 and 97.5\% quantile of the modeled Student. Kolmogorov Smirnov test: comparison between the empirical cumluative distribution function (cdf) and the cdf of the empirical Student.} 
## \label{tab:validation}
## \end{table}
## * Additional table: robust (not in the article)
if(FALSE){
    df.table2 <- rbind(cbind(scenario = "A", parameter = greek$a, as.data.frame(checkStud.robustMM[row$a,keep.col])),
                       cbind(scenario = "B", parameter = greek$b, as.data.frame(checkStud.robustfactor[row$b,keep.col])),
                       cbind(scenario = "C", parameter = greek$c, as.data.frame(checkStud.robustlvm[row$c,keep.col]))
                       )

    df.table2$scenario <- as.character(df.table2$scenario)
    rownames(df.table2) <- NULL
    df.table2$se.empirical <- as.character(round(df.table2$se.empirical, digits = 3))
    df.table2$df.empirical <- as.character(round(df.table2$df.empirical, digits = 1))
    df.table2$mean.se <- as.character(round(df.table2$mean.se, digits = 3))
    df.table2$mean.df <- as.character(round(df.table2$mean.df, digits = 1))
    df.table2$Etype1 <- as.character(round(df.table2$Etype1, digits = 3))
    df.table2$tdist.ksD <- as.character(round(df.table2$tdist.ksD, digits = 3))
    df.table2$tdist.ksP <- format.pval(df.table2$tdist.ksP, digits = 3, eps = 1e-3)
    df.table2[duplicated(df.table2$scenario),"scenario"] <- ""

    addtorow <- list()
    addtorow$pos <- list(0,0,2,6)
    addtorow$command <- c("&&\\multicolumn{2}{c}{Empirical Student} & \\multicolumn{2}{c}{Modeled Student} & Expected
& \\multicolumn{2}{c}{KS-test} \\\\ \\cmidrule(lr){3-4} \\cmidrule(lr){5-6} \\cmidrule(lr){8-9} \n",
"Scenario & parameter & se & df & se & df & type 1 error & statistic & p-value  \\\\\n",
"[4mm] ","[4mm] ")
    print(xtable::xtable(df.table2,
                         label = "tab:robustvalidation",
                         caption = paste("Comparison between the empirical distribution of the robust Wald statistic vs. the modeled distribution (after correction) for n=20. ",
                                         "Empirical Student: standard error and degrees of freedom of a Student's t-distribution fitted on the empirical distribution.",
                                         "Modeled Student: average estimated degrees of freedom over the simulations.",
                                         "Expected type 1 error: rejection rate under the empirical Student when using the 2.5 and 97.5\\% quantile of the modeled Student.",
                                         "Kolmogorov Smirnov test: comparison between the empirical cumluative distribution function (cdf) and the cdf of the empirical Student."),
                         ),
          add.to.row =  addtorow,
          include.colnames = FALSE,
          include.rownames=FALSE,
          sanitize.text.function = function(x){x},
          booktabs = TRUE)
}
