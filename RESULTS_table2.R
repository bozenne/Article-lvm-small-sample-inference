source("ANALYSIS.R")

row.ML <- c("residual variance" = as.double(coef.vitamin.ML["w1~~w1"]),
            "variance random intercept" = as.double(coef.vitamin.ML["eta~~eta"]),
            "statistic" = as.double(Ftest.vitamin.ML[["statistic"]])/3,
            "degrees of freedom" = Inf,
            "p-value" = as.double(Ftest.vitamin.ML[["p.value"]]))
row.correctedML <- c("residual variance" = as.double(coef.vitamin.MLc["w1~~w1"]),
                     "variance random intercept" = as.double(coef.vitamin.MLc["eta~~eta"]),
                     "statistic" = as.double(Ftest.vitamin.MLc[["statistic"]]),
                     "degrees of freedom" = as.double(Ftest.vitamin.MLc[["parameter"]]),
                     "p-value" = as.double(Ftest.vitamin.MLc[["p.value"]]))
row.REML <- c("residual variance" = as.double(coef.vitamin.REML["sigma2"]),
              "variance random intercept" = as.double(coef.vitamin.REML["tau"]),
              "statistic" = as.double(Ftest.vitamin.REML[["F value"]][2]),
              "degrees of freedom" = as.double(Ftest.vitamin.REML[["DenDF"]][2]),
              "p-value" = as.double(Ftest.vitamin.REML[["Pr(>F)"]][2]))
M <- cbind("ML" = row.ML, 
           "ML with correction" = row.correctedML,
           "REML" = row.REML)
M.txt <- rbind(format(M[1:3,], digits = 3, nsmall = 3),
               format(M[4,,drop=FALSE], digits = 2, nsmall = 2),
               format(M[5,,drop=FALSE], digits = 2, nsmall = 2))
M.txt[is.infinite(M)] <- "$\\infty$"

print(xtable::`align<-`(xtable::xtable(M.txt, label = "tab:mixed"), "rccc"),
      include.colnames = TRUE,
      include.rownames = TRUE,
      sanitize.text.function = function(x){x},
      booktabs = TRUE)

## \begin{table}[ht]
## \centering
## \begin{tabular}{rccc}
##   \toprule
##  & ML & ML with correction & REML \\ 
##   \midrule
## residual variance & 0.148 & 0.177 & 0.176 \\ 
##   variance random intercept & 0.349 & 0.389 & 0.389 \\ 
##   statistic & 5.024 & 4.019 & 4.247 \\ 
##   degrees of freedom & $\infty$ & 49.11 & 43.59 \\ 
##   p-value & 0.0018 & 0.0123 & 0.0102 \\ 
##    \bottomrule
## \end{tabular}
## \label{tab:mixed}
## \end{table}
