### FCT.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 17 2019 (14:14) 
## Version: 
## Last-Updated: mar 18 2019 (11:15) 
##           By: Brice Ozenne
##     Update #: 110
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * sinkDirectory
sinkDirectory <- function (path, string.keep = NULL, string.exclude = NULL, addMissingCol = FALSE, 
                           fixed = FALSE, trace = TRUE) 
{
    if (!dir.exists(path)) {
        possible.dirs <- setdiff(list.dirs(file.path(path, ".."), 
                                           full.names = FALSE), "")
        stop("Directory ", path, " does not exists \n", "Existing dir. in parent dir.: \"", 
             paste0(possible.dirs, collapse = "\" \""), "\"\n")
    }
    allFiles <- list.files(path)
    if (length(allFiles) == 0) {
        warning("Empty directory \n")
        return(NULL)
    }
    index.file <- 1:length(allFiles)
    if (!is.null(string.keep)) {
        index.file <- intersect(index.file, grep(string.keep, 
                                                 allFiles, fixed = fixed))
    }
    if (!is.null(string.exclude)) {
        index.file <- setdiff(index.file, grep(string.exclude, 
                                               allFiles, fixed = fixed))
    }
    n.files <- length(index.file)
    if (n.files == 0) {
        warning("No file matching the query \n")
        return(NULL)
    }
    dt.merge <- NULL
    if (trace) {
        cat("read ", n.files, " files \n", sep = "")
        pb <- txtProgressBar(max = n.files)
    }
    for (iFile in 1:n.files) {
        iFilename <- allFiles[index.file[iFile]]
        iRes <- data.table::as.data.table(readRDS(file = file.path(path, 
                                                                   iFilename)))
        if (NROW(iRes) == 0) {
            next
        }
        iRes[, `:=`(iFile, iFilename)]
        if (!is.null(dt.merge) && addMissingCol == TRUE && NCOL(dt.merge) != 
            NCOL(iRes)) {
            missing <- setdiff(names(dt.merge), names(iRes))
            if (!is.null(missing)) {
                for (iMiss in missing) {
                    vec.tempo <- dt.merge[1, .SD, .SDcols = iMiss][[1]]
                    vec.tempo[1] <- NA
                    iRes[, `:=`(c(iMiss), vec.tempo[1])]
                }
            }
            missing <- setdiff(names(iRes), names(dt.merge))
            if (!is.null(missing)) {
                for (iMiss in missing) {
                    vec.tempo <- iRes[1, .SD, .SDcols = iMiss][[1]]
                    vec.tempo[1] <- NA
                    dt.merge[, `:=`(c(iMiss), vec.tempo[1])]
                }
            }
        }
        dt.merge <- rbind(dt.merge, iRes)
        if (trace) {
            setTxtProgressBar(pb, iFile)
        }
    }
    if (trace) {
        close(pb)
    }
    return(dt.merge)
}

## * createFigure
createFigure <- function(data, robust, link, vec.name, vec.label){
    data <- data.table::copy(data)
    
    index.keep <- which(data$robust==robust)
    data <- data[index.keep]
    
    if(any(link %in% data$link == FALSE)){
        txt <- link[link %in% data$link == FALSE]
        stop("Incorrect link: \"",paste(txt, collapse = "\" \""),"\" \n")
    }
    index.keep <- which(data$link %in% link)
    data <- data[index.keep]

    link2link.txt <- data[!duplicated(link),setNames(link.txt,link)]
    data$link.txt <- factor(data$link.txt, levels = link2link.txt[link])
    data$correction <- factor(data$correction, levels = vec.name)
    setkeyv(data, cols = c("link.txt","correction"))
    if(any(data[,.N,by = c("correction","n","link.txt")][["N"]]!=1)){
        stop("The combination \"correction\" \"n\" \"link.txt\" must lead exactly to one point \n")
    }
    n.name <- length(vec.name)
    seq.shape <- seq(from = 15, by = 1, length.out = n.name)

    gg <- ggplot(data, aes(x=as.factor(n),y=type1,group=correction,color=correction,shape=correction))
    gg <- gg + geom_abline(intercept = 0.05, slope = 0, size = 2)
    gg <- gg + geom_point(size = 4) + geom_line(size = 2)
    gg <- gg + facet_grid(~link.txt, labeller = label_parsed)
    gg <- gg + xlab("sample size")
    gg <- gg + theme(text = element_text(size = 25)) + theme(legend.key.height = unit(0.05, "npc"), legend.key.width = unit(0.08, "npc"))
    gg <- gg + scale_color_manual("",
                                  breaks = vec.name,
                                  label = vec.label,
                                  values = ggthemes::colorblind_pal()(n.name+2)[c(6,2:4)])
    gg <- gg + scale_shape_manual("",
                                  breaks = vec.name,
                                  label = vec.label,
                                  values = seq.shape)
    gg <- gg + ylab("type 1 error rate") 
    gg <- gg + theme(legend.position = "bottom") + guides(color=guide_legend(nrow=2,byrow=TRUE))
    gg <- gg + scale_y_continuous(breaks = scales::pretty_breaks(n = 5))
    return(gg)
}

## * createFigureBIS
createFigureBIS <- function(data, n, robust, link = NULL, vec.name, vec.label){
    data <- data.table::copy(data)[!is.na(link.txt)]    
    
    index.keep <- which(data$robust==robust)
    data <- data[index.keep]

    index.keep <- which(data$n%in%n)
    data <- data[index.keep]
    data[,n := paste0("sample size = ",factor(data$n))]
    
    if(!is.null(link)){
        if(any(link %in% data$link == FALSE)){
            txt <- link[link %in% data$link == FALSE]
            stop("Incorrect link: \"",paste(txt, collapse = "\" \""),"\" \n")
        }
        index.keep <- which(data$link %in% link)
        data <- data[index.keep]
    }else{
        link <- unique(data$link)
    }
    
    link2link.txt <- data[!duplicated(link),setNames(link.txt,link)]
    data$link.txt <- factor(data$link.txt, levels = link2link.txt[link])
    data$correction <- factor(data$correction, levels = vec.name)
    setkeyv(data, cols = c("link.txt","correction"))

    if(any(data[,.N,by = c("correction","n","link.txt")][["N"]]!=1)){
        stop("The combination \"correction\" \"n\" \"link.txt\" must lead exactly to one point \n")
    }

    vecX.label <- sapply(unique(data$link.txt), function(iLink){
        text = strsplit(gsub(pattern = "paste(", replacement = "", x = as.character(iLink), fixed = TRUE),
                        split = ",")[[1]][1]
    })
    vecX.expr <- eval(parse(text = paste0("c(",paste(paste("expression(",vecX.label,")",sep=""),collapse=","),")")))
    
    n.name <- length(vec.name)
    seq.shape <- seq(from = 15, by = 1, length.out = n.name)
    gg <- ggplot(data, aes(x=link.txt,y=type1,group=correction,color=correction,shape=correction))
    gg <- gg + geom_abline(intercept = 0.05, slope = 0, size = 2)
    gg <- gg + geom_point(size = 2.5) + geom_line(size = 2)
    gg <- gg + facet_grid(~n)
    gg <- gg + theme(text = element_text(size = 25)) + theme(legend.key.height = unit(0.05, "npc"), legend.key.width = unit(0.08, "npc"))
    gg <- gg + xlab("") + scale_x_discrete(labels = vecX.expr)

    
    gg <- gg + scale_color_manual("",
                                  breaks = vec.name,
                                  label = vec.label,
                                  values = ggthemes::colorblind_pal()(n.name+2)[c(6,2:4)])
    gg <- gg + scale_shape_manual("",
                                  breaks = vec.name,
                                  label = vec.label,
                                  values = seq.shape)
    gg <- gg + ylab("type 1 error rate") 
    gg <- gg + theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1)) + guides(color=guide_legend(nrow=2,byrow=TRUE))
    gg <- gg + scale_y_continuous(breaks = scales::pretty_breaks(n = 5))
    return(gg)
}

## * groupFigures
groupFigures <- function(figure1, figure2, figure3, reduce.x = TRUE){
    x.factor <- c("20" = "black", "30" = "black", "50" = "black", "75" = "black", "100" = "black",
                  "150" = "transparent", "200" = "black", "300" = "transparent", "500" = "black")
    x.lvm <- c("20" = "black", "30" = "transparent", "50" = "black", "75" = "transparent", "100" = "black",
               "150" = "transparent", "200" = "black", "300" = "transparent", "500" = "black")

    figure1 <- figure1 + ggtitle("scenario (A): mixed model") + xlab("") + theme(text = element_text(size = 10), legend.position="none")
    figure2 <- figure2 + ggtitle("scenario (B): single factor model") + xlab("") + theme(text = element_text(size = 10), legend.position="none")
    if(reduce.x == TRUE){
        figure2 <- figure2 + theme(axis.text.x=element_text(color=x.factor))
    }
    figure3 <- figure3 + ggtitle("scenario (C): two latent variables model") + theme(text = element_text(size = 10))
    figure3 <- figure3 + guides(color=guide_legend(nrow=2,byrow=TRUE)) + theme(legend.margin=margin(t = -0.25, unit='cm'))
    if(reduce.x == TRUE){
        figure3 <- figure3 + theme(axis.text.x=element_text(color=x.lvm))
    }
    ls.figure <- list(figure1 + theme(plot.margin=unit(c(0,0.1,0,0),"cm")),
                      figure2 + theme(plot.margin=unit(c(0,0.1,0,0),"cm")),
                      figure3 + theme(plot.margin=unit(c(0,0.1,0,0),"cm")))
    arrange.figure <- gridExtra::grid.arrange(grobs=ls.figure[1:3], ncol=1, widths = 10, 
                                              heights=unit(c(7,7,9), c("cm", "cm")))

    return(arrange.figure)

}


## * createTable
createTable <- function(dt, seqN, seqType, digit, convert2latex = TRUE){

    reformat <- function(value, digit){ ## value <- 0 ; digit = 2
        value.round <- round(value, digits = digit)
        vec <- format(value.round, nsmall = digit)
        if(any(value.round==0)){
            vec[value.round==0] <- paste0("<0.",paste(rep(0,digit-1), collapse=""),"1",sep="")
        }
        return(vec)
    }
    ## reformat(0, digit = 2)


    
    tab1 <- dt[(corrected == FALSE) & (type %in% seqType),
               .("mean (sd)"=paste0(reformat(mean, digit = digit)," (",reformat(sd, digit = digit),")"),
                 "median [Q1;Q3]"=paste0(reformat(middle, digit = digit)," [",reformat(lower, digit = digit),";",reformat(upper, digit = digit),"]")),
               by = c("n","type") ]
    ## convert n to character without warning
    tab1[, "sample size" := as.character(n)]
    ## put sample size firt column
    setcolorder(tab1, c("sample size",setdiff(names(tab1),"sample size")))
    tab2 <- dt[(corrected == TRUE) & (type %in% seqType),
               .("mean (sd)"=paste0(reformat(mean, digit = digit)," (",reformat(sd, digit = digit),")"),
                 "median [Q1;Q3]"=paste0(reformat(middle, digit = digit)," [",reformat(lower, digit = digit),";",reformat(upper, digit = digit),"]")),
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

## * fitStudent
fitStudent <- function(data, robust, seqLink = NULL, trace = TRUE,
                       quantile.min = 0, quantile.max = 1, value.max = Inf){

    if(is.null(seqLink)){
        seqLink <- unique(data$name)
    }else if(any(seqLink %in% data$name == FALSE)){
        txt <- seqLink[seqLink %in% data$name == FALSE]
        stop("Unknown link(s): \"",paste(txt, collapse = "\" \""),"\"\n")
    }

    if(robust){
        data$wald <- data[["estimate.MLcorrected"]]/data[["se.robustMLcorrected"]]
        data$wald.df <- data[["df.robustMLcorrected"]]
    }else{
        data$wald <- data[["estimate.MLcorrected"]]/data[["se.MLcorrected"]]
        data$wald.df <- data[["df.MLcorrected"]]
    }
    data <- data[!is.na(wald) & !is.na(wald.df)]
    
    n.link <- length(seqLink)
    M.dist <- matrix(nrow = n.link, ncol = 10,
                     dimnames = list(seqLink, c("se.empirical","df.empirical","mean.se","mean.df","Etype1","zdist.ksD","zdist.ksP","tdist.ksD","tdist.ksP","empirical.cv")))
    
    if(trace>0){pb <- txtProgressBar(max = n.link)}
    for(iIndex in 1:n.link){ ## iLink <- 1
        iLink <- seqLink[iIndex]
        iData <- data[name == iLink]

        iMin <- quantile(iData$wald, probs = quantile.min)
        iMax <- quantile(iData$wald, probs = quantile.max)

        if(trace>0){setTxtProgressBar(pb, iIndex)}
        M.dist[iLink, "mean.se"] <- 1
        M.dist[iLink, "mean.df"] <- mean(iData$wald.df)

        iData2 <- iData[wald >= iMin & wald <= iMax & abs(wald)<value.max]
        iFitted <- try(QRM_fit.st(iData2$wald), silent = FALSE)

        if(!inherits(iFitted, "try-error") && iFitted$converged ){
            M.dist[iLink,"empirical.cv"] <- iFitted$converged
            M.dist[iLink,"se.empirical"] <- iFitted$par.est["sigma"]
            M.dist[iLink,"df.empirical"] <- iFitted$par.est["df"]

            iKS.z <- ks.test(iData2$wald, y = pnorm, mean = mean(iData2$wald), sd = sd(iData2$wald))
            M.dist[iLink,"zdist.ksD"] <- as.double(iKS.z$statistic)
            M.dist[iLink,"zdist.ksP"] <- as.double(iKS.z$p.value)
            
            iX.norm <- (iData2$wald - iFitted$par.est["mean"])/iFitted$par.est["sigma"]
            iKS.t <- ks.test(iX.norm, y = pt, df = iFitted$par.est["df"])
            M.dist[iLink,"tdist.ksD"] <- as.double(iKS.t$statistic)
            M.dist[iLink,"tdist.ksP"] <- as.double(iKS.t$p.value)
        }
        
        q_alpha <- qt(c(0.025,0.975), df = M.dist[iLink,"mean.df"], ncp = 0) * M.dist[iLink,"mean.se"]
        M.dist[iLink,"Etype1"] <- 1-diff(pt(q_alpha/M.dist[iLink,"se.empirical"], df = M.dist[iLink,"df.empirical"], ncp = 0))
        
    }
    if(trace>0){close(pb)}

    return(M.dist)
}

## QRM_fit.st(rt(1e4, df = 5))
## QRM_fit.st(rt(1e4, df = 5)*2+1)
## QRM_fit.st(rnorm(1e4))
QRM_fit.st <- function (data, ...) 
{
    ## modification of the fit.st function from the QRM pakcage
    
    ## extract statistics
    mu <- mean(data)
    m2 <- mean((data - mu)^2)
    m4 <- mean((data - mu)^4)
    nu <- max(3,4 + (6 * m2^2)/(m4 - 3 * m2^2)) ## avoid wrong initialization
    sigma <- sqrt((nu - 2) * m2/nu)

    ## initialize parameters
    theta <- c(sigma, nu)

    ## define objective function
    negloglik <- function(theta, y) {
        -sum(log(dt((y - mu)/abs(theta[1]), df = abs(theta[2]))) - 
             log(abs(theta[1])))
    }

    ## optimization
    optim.1 <- try(optim(theta, fn = negloglik, y = data, method = "L-BFGS-B", lower = c(1e-6, 1), upper = c(Inf, 200),
                         control = list(trace = 0)), silent = TRUE)

    optim.2 <- try(optim(theta, fn = negloglik, y = data, method = "Nelder-Mead", control = list(trace = 0)),
                   silent = TRUE)

    if(inherits(optim.1,"try-error") && inherits(optim.2,"try-error")){
        optimfit <- list(par = c(NA,NA), value = Inf, convergence = 1, method = "none")
    }else{
    
        if(inherits(optim.1,"try-error")){
            optim.1 <- list(par = c(NA,NA), value = Inf, convergence = 1, method = "none")
        }else{
            optim.1$method <- "L-BFGS-B"
        }
        if(inherits(optim.2,"try-error")){
            optim.2 <- list(par = c(NA,NA), value = Inf, convergence = 1, method = "none")
        }else{
            optim.2$method <- "Nelder-Mead"
        }

        if(optim.2$value>optim.1$value || (optim.2$par[2]<=1)){
            optimfit <- optim.1
        }else{
            optimfit <- optim.2
        }
    }

    ## epxort
    par.ests <- optimfit$par
    loglh.max <- -negloglik(par.ests, y = data)
    out <- list(converged = (optimfit$convergence==0),
                par.ests = setNames(c(mu, par.ests), c("mean","sigma","df")),
                ll.max = loglh.max)
    return(out)
}

######################################################################
### FCT.R ends here
