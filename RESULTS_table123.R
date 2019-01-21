## path <- "P:/Cluster/LVMproject/article-smallSampleInference"
## setwd(path)

library(data.table)
export <- TRUE

## * path
path.results <- "./Results"

## * load data
dtLS.sim.MMbias <- readRDS(file.path(path.results,"bias-simulation-mixedModel.rds"))

dtLS.sim.factorbias <- readRDS(file.path(path.results,"bias-simulation-factorModel.rds"))

dtLS.sim.lvmbias <- readRDS(file.path(path.results,"bias-simulation-lvmModel.rds"))

## * function
createTable <- function(dt, seqN, seqType, digit, convert2latex = TRUE){

    tab1 <- dt[(corrected == FALSE) & (type %in% seqType),
               .("mean (sd)"=paste0(format(mean, digits = digit)," (",format(sd, digits = digit),")"),
                 "median [Q1;Q3]"=paste0(format(middle, digits = digit)," [",format(lower, digits = digit),";",format(upper, digits = digit),"]")),
               by = c("n","type") ]
    ## convert n to character without warning
    tab1[, "sample size" := as.character(n)]
    ## put sample size firt column
    setcolorder(tab1, c("sample size",setdiff(names(tab1),"sample size")))

    tab2 <- dt[(corrected == TRUE) & (type %in% seqType),
               .("mean (sd)"=paste0(format(mean, digits = digit)," (",format(sd, digits = digit),")"),
                 "median [Q1;Q3]"=paste0(format(middle, digits = digit)," [",format(lower, digits = digit),";",format(upper, digits = digit),"]")),
               by = c("n","type") ]

    n.n <- length(seqN)
    tab.tex <- NULL
    for(iN in 1:n.n){ ## iN <- 1
        iTab <- cbind(tab1[n==seqN[iN],],
                      tab2[n==seqN[iN],.SD,.SDcols = setdiff(names(tab2),c("n","type"))])

        ## only keep first n
        iTab[, "n" := NULL]
        iTab[-1, "sample size" := ""]

        if(convert2latex){
            ## scientific notation for the coefficients
            iTab$type[iTab$type=="Sigma_var"] <- "Sigma_{\\varepsilon}"
            iTab$type[iTab$type=="Sigma_cov"] <- "Sigma_{\\varepsilon,\\varepsilon'}"
            iTab$type[iTab$type=="Psi_var"] <- "Sigma_{\\zeta}"
            iTab$type[iTab$type=="Psi_var"] <- "Sigma_{\\zeta,\\zeta'}"
            iTab$type <- paste0("$\\",iTab$type,"$")

            ## create table with coefficients
            iTab.tex <- paste(apply(iTab,1,paste,collapse=" & "),"\\\\ \n")
        
            ## add column title
            if(iN==1){
                iTab.tex <- c(paste("& & ",paste(names(iTab)[-(1:2)],collapse=" & "),"\\\\ \\hline \n"),
                              iTab.tex)
                iTab.tex <- c("\\multirow{2}{*}{sample size} & \\multirow{2}{*}{type} & \\multicolumn{2}{c}{ML} & \\multicolumn{2}{c}{ML-corrected} \\\\ \n",
                              iTab.tex)
            }
            tab.tex <- c(tab.tex, iTab.tex)
        }else{
            tab.tex <- rbind(tab.tex, iTab)
        }
    }
    ## wrap in tabular envir
    if(convert2latex){
        tab.tex <- c("\\begin{tabular}{c|c|cccc} \n",
                     tab.tex,
                     "\\end{tabular} \n \n")
    }
    return(tab.tex)
}

## * table 1
table1 <- createTable(dt = dtLS.sim.MMbias, seqN = c(20,30,50,100), seqType = c("Sigma_var","Psi_var"), digit = 3, convert2latex = TRUE)
table2 <- createTable(dt = dtLS.sim.factorbias, seqN = c(20,30,50,100), seqType = c("Sigma_var","Psi_var"), digit = 3, convert2latex = TRUE)
table3 <- createTable(dt = dtLS.sim.lvmbias, seqN = c(20,30,50,100), seqType = c("Sigma_var","Psi_var","Sigma_cov"), digit = 3, convert2latex = TRUE)
