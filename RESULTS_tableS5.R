path.results <- "./Results"

## * packages
library(data.table)

## * Import

## ** Simulation results
dtLS.sim.Comptype1 <- readRDS(file.path(path.results,"type1error-simulation-comparison.rds"))

## * Data processing
greek.label.factor2 <- c("eta" = "$\\alpha$",
                         "Y2" = "$\\nu_2$",
                         "Y3" = "$\\nu_3$",
                         "Y4" = "$\\nu_4$",
                         "Y2~Gene1Y" = "$k_1$",
                         "eta~Gene2Y" = "$\\gamma_1$",
                         "eta~Age" = "$\\gamma_2$",
                         "Y2~eta" = "$\\lambda_2$",
                         "Y3~eta" = "$\\lambda_3$",
                         "Y4~eta" = "$\\lambda_4$")

dtLS.sim.Comptype1[link %in% names(greek.label.factor2), link.txt := factor(link,
                                                                            levels = names(greek.label.factor2),
                                                                            labels = as.character(greek.label.factor2))]
dtLS.sim.Comptype1[estimator == "IV", estimator := "IV (lava)"]
dtLS.sim.Comptype1[estimator == "IVlavaan", estimator := "IV (lavaan)"]
dtLS.sim.Comptype1[estimator == "robustML", estimator := "robust ML"]

## * Table 5
addtorow <- list()
addtorow$pos <- list(10)
addtorow$command <- "[4mm]"
table.R1 <- dcast(dtLS.sim.Comptype1,
                  formula = link.txt+n ~ estimator, value.var = "type1")[n%in%c(20,50)]
setnames(table.R1, old = c("link.txt"), new = c("parameter"))
table.R1 <- table.R1[order(table.R1$n),]
table.R1$n[duplicated(table.R1$n)] <- ""

order.col <- c("parameter","n","ML","robust ML","GLS","WLS","IV (lava)","IV (lavaan)")
print(xtable::xtable(table.R1[,.SD,.SDcols = order.col], label = "tab:comparison", caption = "Comparison of the type 1 error of Wald tests for various estimation methods in Scenario B under a correctly specified model. No small sample correction is used. Robust ML corresponds to the use of robust Wald tests. NA indicates that the estimator never converged in the simulations and so the Wald statistic could not be computed. The R packages lava and MIIVsem were used for IV estimation. GLS and WLS were carried out using the R package lavaan.", digits = 3),
      add.to.row =  addtorow, NA.string="NA",
      include.rownames = FALSE, booktabs = TRUE, sanitize.text.function = function(x){x})

## \begin{tabular}{llrrrrrr}
##   \toprule
## parameter & n & ML & robust ML & GLS & WLS & IV (lava) & IV (lavaan) \\ 
##   \midrule
## $\alpha$ & 20 & 0.088 & 0.094 & 0.075 & NA & 0.108 & 0.087 \\ 
##   $\nu_2$ &  & 0.082 & 0.093 & 0.067 & NA & 0.101 & 0.081 \\ 
##   $\nu_3$ &  & 0.062 & 0.068 & 0.053 & NA & 0.075 & 0.070 \\ 
##   $\nu_4$ &  & 0.073 & 0.074 & 0.057 & NA & 0.077 & 0.073 \\ 
##   $k_1$ &  & 0.115 & 0.108 & 0.097 & NA & 0.088 & 0.081 \\ 
##   $\gamma_1$ &  & 0.096 & 0.089 & 0.081 & NA & 0.094 & 0.085 \\ 
##   $\gamma_2$ &  & 0.107 & 0.114 & 0.135 & NA & 0.126 & 0.088 \\ 
##   $\lambda_2$ &  & 0.089 & 0.094 & 0.112 & NA & 0.161 & 0.128 \\ 
##   $\lambda_3$ &  & 0.086 & 0.093 & 0.114 & NA & 0.180 & 0.143 \\ 
##   $\lambda_4$ &  & 0.075 & 0.086 & 0.088 & NA & 0.101 & 0.067 \\ 
##    [4mm]$\alpha$ & 50 & 0.062 & 0.064 & 0.061 & 0.046 & 0.068 & 0.062 \\ 
##   $\nu_2$ &  & 0.058 & 0.062 & 0.055 & 0.092 & 0.063 & 0.057 \\ 
##   $\nu_3$ &  & 0.054 & 0.055 & 0.049 & 0.108 & 0.056 & 0.054 \\ 
##   $\nu_4$ &  & 0.061 & 0.062 & 0.057 & 0.077 & 0.061 & 0.061 \\ 
##   $k_1$ &  & 0.065 & 0.063 & 0.061 & 0.062 & 0.059 & 0.056 \\ 
##   $\gamma_1$ &  & 0.061 & 0.057 & 0.069 & 0.015 & 0.064 & 0.061 \\ 
##   $\gamma_2$ &  & 0.070 & 0.080 & 0.075 & 0.108 & 0.077 & 0.062 \\ 
##   $\lambda_2$ &  & 0.063 & 0.067 & 0.062 & 0.062 & 0.092 & 0.079 \\ 
##   $\lambda_3$ &  & 0.059 & 0.065 & 0.063 & 0.246 & 0.100 & 0.086 \\ 
##   $\lambda_4$ &  & 0.060 & 0.066 & 0.068 & 0.185 & 0.070 & 0.054 \\ 
##    \bottomrule
## \end{tabular}
