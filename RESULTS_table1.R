## path <- "P:/Cluster/LVMproject/article-smallSampleInference"
## setwd(path)

path.results <- "./Results"
path.data <- "./data/"

library(data.table)
source("FCT.R") ## get function fitStudent

path.simulation.mixedModel <- "./Results/simulation-mixedModel"
path.simulation.factorModel <- "./Results/simulation-factorModel"
path.simulation.lvm <- "./Results/simulation-lvm"


## * load data
dist.simulation <- readRDS(file.path(path.results,"dist-simulation.rds"))

## * fit distributions
if(FALSE){ ## check
    set.seed(10)
    dt.test <- rbind(data.table(name = "rt", estimate.MLcorrected = rt(1e3, df = 3), se.MLcorrected = 1, df.MLcorrected = 3),
                     data.table(name = "rnorm", estimate.MLcorrected = rnorm(1e3), se.MLcorrected = 1, df.MLcorrected = Inf),
                     data.table(name = "runif", estimate.MLcorrected = runif(1e3), se.MLcorrected = NA, df.MLcorrected = NA))
    fitStudent(dt.test)
    ##       se.empirical df.empirical mean.se mean.df     Etype1  zdist.ksD    zdist.ksP  tdist.ksD   tdist.ksP empirical.cv
    ## rt       0.9691165 2.533170e+00       1       3 0.05892354 0.10524949 4.778306e-10 0.01482152 0.980536804            1
    ## rnorm    1.0146781 2.196487e+02       1     Inf 0.05469411 0.01906850 8.602872e-01 0.01870678 0.875256720            1
    ## runif    0.2856405 4.226251e+07     NaN     NaN        NaN 0.06152528 1.030642e-03 0.06164236 0.001001343            1

}

checkStud.MM <- fitStudent(dist.simulation$MM[n==20],
                           seqLink = c("Y2","eta~Gene1Y"))
##            se.empirical df.empirical   mean.se  mean.df     Etype1   zdist.ksD zdist.ksP   tdist.ksD tdist.ksP empirical.cv
## Y2            0.3148938     14497.46 0.3137851 36.66667 0.04343510 0.003385657 0.9763259 0.003371633 0.9772740            1
## eta~Gene1Y    0.5494933    183493.45 0.5394871 18.33333 0.03940187 0.005588690 0.5622115 0.005599168 0.5597803            1

checkStud.factor <- fitStudent(dist.simulation$factor[n==20],
                               seqLink = c("Y2","Y1~Gene2Y","eta~Gene1Y","Y4~eta"))
##            se.empirical df.empirical   mean.se   mean.df     Etype1   zdist.ksD    zdist.ksP   tdist.ksD tdist.ksP empirical.cv
## Y2            0.4109210     8.072971 0.4561871  9.109677 0.03630368 0.020565057 5.535884e-07 0.005092165 0.7436142            1
## Y1~Gene2Y     0.5559323    70.766115 0.5181065 17.079044 0.05327330 0.005865745 5.708261e-01 0.004985804 0.7665729            1
## eta~Gene1Y    0.5503309    22.004143 0.5268317 10.547688 0.04569939 0.009336589 8.898157e-02 0.004793287 0.8066604            1
## Y4~eta        0.1649963     5.241632 0.1860292  5.908417 0.03744640 0.043599304 0.000000e+00 0.004420286 0.8764785            1

checkStud.lvm <- fitStudent(dist.simulation$lvm[n==20],
                            seqLink = c("Y2","Y1~Gene2Y","eta1~Gene1Y","Y4~eta1","eta1~eta2","Y1~~Y2"))
##             se.empirical df.empirical   mean.se   mean.df     Etype1   zdist.ksD    zdist.ksP   tdist.ksD    tdist.ksP empirical.cv
## Y2             0.4036533     7.114034 0.4489476 10.134260 0.04206787 0.049700838 0.000000e+00 0.009800696 0.0793768135            1
## Y1~Gene2Y      0.5490824    49.027617 0.4902965 17.112820 0.06564178 0.005256191 7.422624e-01 0.003790330 0.9692818285            1
## eta1~Gene1Y    0.5514802    20.813113 0.4926512 10.114808 0.06019669 0.011302646 2.737222e-02 0.005955252 0.5906315106            1
## Y4~eta1        0.1669083     5.446397 0.1863481  6.015801 0.03774302 0.068214988 0.000000e+00 0.007035541 0.3766476984            1
## eta1~eta2      0.2451067     3.652047 0.2772663  3.554216 0.03419134 0.093078627 0.000000e+00 0.005478091 0.6946017516            1
## Y1~~Y2         0.3632980    14.206881 0.3440811  7.804898 0.04538257 0.021471003 3.762927e-07 0.015343410 0.0007354053


rm.col <- c("empirical.cv","zdist.ksD","zdist.ksP")
keep.col <- setdiff(colnames(checkStud.MM),rm.col)

greek <- list(a = c("$\\nu_2$","$\\gamma_1$"),
              b = c("$\\nu_2$","$\\lambda_4$","$\\gamma_1$","$k_1$"),
              c = c("$\\nu_2$","$k_1$","$\\lambda_4$","$\\gamma_1$","$b_1$","$\\sigma_{12}$"))
row <- list(a = c("Y2","eta~Gene1Y"),
            b = c("Y2","Y4~eta","eta~Gene1Y","Y1~Gene2Y"),
            c = c("Y2","Y1~Gene2Y","Y4~eta1","eta1~Gene1Y","eta1~eta2","Y1~~Y2"))

df.table1 <- rbind(cbind(scenario = "a", parameter = greek$a, as.data.frame(checkStud.MM[row$a,keep.col])),
                   cbind(scenario = "b", parameter = greek$b, as.data.frame(checkStud.factor[row$b,keep.col])),
                   cbind(scenario = "c", parameter = greek$c, as.data.frame(checkStud.lvm[row$c,keep.col]))
                   )
rownames(df.table1) <- NULL
df.table1$se.empirical <- as.character(round(df.table1$se.empirical, digits = 3))
df.table1$df.empirical <- as.character(round(df.table1$df.empirical, digits = 1))
df.table1$mean.se <- as.character(round(df.table1$mean.se, digits = 3))
df.table1$mean.df <- as.character(round(df.table1$mean.df, digits = 1))
df.table1$Etype1 <- as.character(round(df.table1$Etype1, digits = 3))
df.table1$tdist.ksD <- as.character(round(df.table1$tdist.ksD, digits = 3))
df.table1$tdist.ksP <- format.pval(df.table1$tdist.ksP, digits = 3, eps = 1e-3)
df.table1[duplicated(df.table1$scenario),"scenario"] <- ""

addtorow <- list()
addtorow$pos <- list(0,0,2,6)
addtorow$command <- c("&&\\multicolumn{2}{c}{empirical Student} & \\multicolumn{2}{c}{modeled Student} & expected
& \\multicolumn{2}{c}{KS-test} \\\\ \\cmidrule(lr){3-4} \\cmidrule(lr){5-6} \\cmidrule(lr){8-9} \n",
"scenario & parameter & se & df & se & df & type 1 error & statistic & p-value  \\\\\n",
"[4mm] ","[4mm] ")
print(xtable::xtable(df.table1,
                     label = "tab:validation",
                     caption = paste("Comparison between the empirical distribution of the LVM parameters vs. the modeled distribution (after correction) for n=20. ",
                                     "Empirical Student: standard error and degree of freedom of a Student's t-distribution fitted on the empirical distribution.",
                                     "Modeled Student: average standard error and degree of freedom for the parameter distribution estimated in each simulation.",
                                     "Expected type 1 error: rejection rate under the empirical Student when using the 2.5 and 97.5\\% quantile of the modeled Student.",
                                     "Kolmogorov Smirnov test: comparison between the empirical cumluative distribution function (cdf) and the cdf of the empirical Student."),
                     ),
      add.to.row =  addtorow,
      include.colnames = FALSE,
      include.rownames=FALSE,
      sanitize.text.function = function(x){x},
      booktabs = TRUE)
