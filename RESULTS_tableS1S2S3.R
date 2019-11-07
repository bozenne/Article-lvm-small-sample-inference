library(data.table)
source("FCT.R") ## get function createTable
export <- TRUE

## * path
path.results <- "./Results"

## * load data
dtLS.sim.MMbias <- readRDS(file.path(path.results,"bias-simulation-mixedModel.rds"))

dtLS.sim.factorbias <- readRDS(file.path(path.results,"bias-simulation-factorModel.rds"))

dtLS.sim.lvmbias <- readRDS(file.path(path.results,"bias-simulation-lvmModel.rds"))

## * table 1
table1 <- createTable(dt = dtLS.sim.MMbias, seqN = c(20,30,50,100), seqType = c("Sigma_var","Psi_var"), digit = 3, convert2latex = TRUE)
cat(table1)
## \begin{tabular}{c|c|cccc} 
##  \multirow{2}{*}{type} & \multirow{2}{*}{sample size}  & \multicolumn{2}{c}{ML}  & \multicolumn{2}{c}{ML-corrected}  \\ 
##  & &  mean (sd) & median [Q1;Q3] & mean (sd) & median [Q1;Q3] \\ \hline 
##  $\Sigma_{\varepsilon}$ & 20 & -0.050 (0.218) & -0.067 [-0.204;0.086] & <0.001 (0.229) & -0.018 [-0.162;0.143] \\ 
##  $\Sigma_{\zeta}$ &  & -0.182 (0.394) & -0.227 [-0.469;0.058] & 0.002 (0.463) & -0.052 [-0.335;0.282] \\ 
##  $\Sigma_{\varepsilon}$ & 30 & -0.034 (0.180) & -0.045 [-0.161;0.080] & -0.001 (0.187) & -0.012 [-0.132;0.117] \\ 
##  $\Sigma_{\zeta}$ &  & -0.128 (0.330) & -0.157 [-0.362;0.075] & -0.006 (0.366) & -0.039 [-0.266;0.219] \\ 
##  $\Sigma_{\varepsilon}$ & 50 & -0.020 (0.140) & -0.027 [-0.118;0.071] & <0.001 (0.143) & -0.007 [-0.100;0.092] \\ 
##  $\Sigma_{\zeta}$ &  & -0.071 (0.264) & -0.089 [-0.260;0.098] & 0.002 (0.281) & -0.017 [-0.198;0.182] \\ 
##  $\Sigma_{\varepsilon}$ & 100 & -0.011 (0.099) & -0.014 [-0.080;0.054] & -0.001 (0.100) & -0.004 [-0.071;0.065] \\ 
##  $\Sigma_{\zeta}$ &  & -0.037 (0.189) & -0.046 [-0.169;0.087] & <0.001 (0.194) & -0.010 [-0.136;0.127] \\ 
##  \end{tabular} 
table2 <- createTable(dt = dtLS.sim.factorbias, seqN = c(20,30,50,100), seqType = c("Sigma_var","Psi_var"), digit = 3, convert2latex = TRUE)
cat(table2)
## \begin{tabular}{c|c|cccc} 
##  \multirow{2}{*}{type} & \multirow{2}{*}{sample size}  & \multicolumn{2}{c}{ML}  & \multicolumn{2}{c}{ML-corrected}  \\ 
##  & &  mean (sd) & median [Q1;Q3] & mean (sd) & median [Q1;Q3] \\ \hline 
##  $\Sigma_{\varepsilon}$ & 20 & -0.125 (0.411) & -0.161 [-0.412;0.118] & -0.029 (0.454) & -0.071 [-0.346;0.238] \\ 
##  $\Sigma_{\zeta}$ &  & -0.150 (0.494) & -0.230 [-0.514;0.124] & 0.018 (0.585) & -0.075 [-0.411;0.343] \\ 
##  $\Sigma_{\varepsilon}$ & 30 & -0.085 (0.349) & -0.107 [-0.321;0.127] & -0.021 (0.372) & -0.046 [-0.273;0.204] \\ 
##  $\Sigma_{\zeta}$ &  & -0.096 (0.415) & -0.148 [-0.399;0.146] & 0.017 (0.464) & -0.041 [-0.320;0.287] \\ 
##  $\Sigma_{\varepsilon}$ & 50 & -0.052 (0.270) & -0.065 [-0.231;0.115] & -0.013 (0.280) & -0.027 [-0.200;0.161] \\ 
##  $\Sigma_{\zeta}$ &  & -0.055 (0.326) & -0.087 [-0.288;0.146] & 0.014 (0.348) & -0.020 [-0.235;0.228] \\ 
##  $\Sigma_{\varepsilon}$ & 100 & -0.026 (0.190) & -0.032 [-0.152;0.095] & -0.006 (0.194) & -0.013 [-0.135;0.117] \\ 
##  $\Sigma_{\zeta}$ &  & -0.026 (0.234) & -0.046 [-0.190;0.119] & 0.008 (0.242) & -0.012 [-0.161;0.158] \\ 
##  \end{tabular} 
table3 <- createTable(dt = dtLS.sim.lvmbias, seqN = c(20,30,50,100), seqType = c("Sigma_var","Psi_var","Sigma_cov"), digit = 3, convert2latex = TRUE)
cat(table3)
## \begin{tabular}{c|c|cccc} 
##  \multirow{2}{*}{type} & \multirow{2}{*}{sample size}  & \multicolumn{2}{c}{ML}  & \multicolumn{2}{c}{ML-corrected}  \\ 
##  & &  mean (sd) & median [Q1;Q3] & mean (sd) & median [Q1;Q3] \\ \hline 
##  $\Sigma_{\varepsilon}$ & 20 & -0.110 (0.402) & -0.148 [-0.394;0.131] & -0.034 (0.437) & -0.076 [-0.342;0.227] \\ 
##  $\Sigma_{\zeta}$ &  & -0.142 (0.513) & -0.231 [-0.522;0.145] & -0.004 (0.584) & -0.104 [-0.438;0.324] \\ 
##  $\Sigma_{\varepsilon,\varepsilon'}$ &  & 0.015 (0.359) & 0.004 [-0.221;0.237] & 0.017 (0.392) & 0.005 [-0.241;0.259] \\ 
##  $\Sigma_{\varepsilon}$ & 30 & -0.075 (0.335) & -0.099 [-0.305;0.132] & -0.024 (0.353) & -0.050 [-0.267;0.193] \\ 
##  $\Sigma_{\zeta}$ &  & -0.094 (0.428) & -0.151 [-0.405;0.152] & <0.001 (0.467) & -0.061 [-0.340;0.272] \\ 
##  $\Sigma_{\varepsilon,\varepsilon'}$ &  & 0.001 (0.296) & -0.006 [-0.198;0.190] & 0.002 (0.313) & -0.005 [-0.209;0.202] \\ 
##  $\Sigma_{\varepsilon}$ & 50 & -0.046 (0.259) & -0.060 [-0.224;0.117] & -0.016 (0.267) & -0.030 [-0.198;0.153] \\ 
##  $\Sigma_{\zeta}$ &  & -0.052 (0.336) & -0.088 [-0.292;0.152] & 0.005 (0.353) & -0.031 [-0.248;0.220] \\ 
##  $\Sigma_{\varepsilon,\varepsilon'}$ &  & -0.003 (0.223) & -0.007 [-0.152;0.143] & -0.002 (0.231) & -0.007 [-0.156;0.149] \\ 
##  $\Sigma_{\varepsilon}$ & 100 & -0.023 (0.181) & -0.031 [-0.147;0.093] & -0.008 (0.184) & -0.016 [-0.134;0.110] \\ 
##  $\Sigma_{\zeta}$ &  & -0.026 (0.237) & -0.043 [-0.195;0.122] & 0.003 (0.243) & -0.015 [-0.170;0.154] \\ 
##  $\Sigma_{\varepsilon,\varepsilon'}$ &  & -0.001 (0.155) & -0.001 [-0.106;0.102] & -0.001 (0.158) & -0.001 [-0.107;0.104] \\ 
##  \end{tabular} 

