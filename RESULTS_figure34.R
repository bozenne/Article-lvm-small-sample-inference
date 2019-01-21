## path <- "P:/Cluster/LVMproject/article-smallSampleInference"
## setwd(path)

library(data.table)
library(ggplot2)
export <- TRUE

## * path
path.results <- "./Results"
path.figures <- "./output"

## * load data
dtLS.sim.MMtype1 <- readRDS(file.path(path.results,"type1error-simulation-mixedModel.rds"))
## unique(dtLS.sim.MMtype1$link)
## dtLS.sim.MMtype1[link == "eta~GenderF", link := "eta~Gene1"]
## dtLS.sim.MMtype1[link == "eta~G", linkn := "eta~Age"]

dtLS.sim.factortype1 <- readRDS(file.path(path.results,"type1error-simulation-factorModel.rds"))
dtLS.sim.factortype1[link == "Y1~Gene1Y", link := "Y1~Gene2Y"]
dtLS.sim.factortype1[link == "eta~Gene2Y", link := "eta~Gene1Y"]
## dtLS.sim.factortype1[link == "eta~G", link := "eta~Age"]

dtLS.sim.lvmtype1 <- readRDS(file.path(path.results,"type1error-simulation-lvmModel.rds"))
## dtLS.sim.lvmtype1[link == "Y1~Age", link := "Y1~Gene2Y"]
## dtLS.sim.lvmtype1[link == "eta1~G1", link := "eta1~Gene1Y"]
## dtLS.sim.lvmtype1[link == "eta2~G2", link := "eta1~Gender"]

## * name parameters
greek.label.MM <- c(eta = expression(paste(eta, "   (within)")),
                    Y2 = expression(paste(nu[2], "   (within)")),
                    Y3 = expression(paste(nu[3], "   (within)")),
                    "eta~Gene1Y" = expression(paste(gamma[1], "   (between)")),
                    "eta~Age" = expression(paste(gamma[2], "   (between)")))
dtLS.sim.MMtype1[, link.txt := factor(link,
                                      levels = names(greek.label.MM),
                                      labels = as.character(greek.label.MM))]

greek.label.factor <- c("eta" = expression(paste(eta, "   (intercept)")),
                        "Y2" = expression(paste(nu[2], "   (intercept)")),
                        "Y3" = expression(paste(nu[3], "   (intercept)")),
                        "Y4" = expression(paste(nu[4], "   (intercept)")),
                        "Y1~Gene2Y" = expression(paste(k[1], "   (reg. OV~OV)")),
                        "eta~Gene1Y" = expression(paste(gamma[1], "   (reg. LV~OV)")),
                        "eta~Age" = expression(paste(gamma[2], "   (reg. LV~OV)")),
                        "Y2~eta" = expression(paste(lambda[2], "   (loading)")),
                        "Y3~eta" = expression(paste(lambda[3], "   (loading)")),
                        "Y4~eta" = expression(paste(lambda[4], "   (loading)"))
                        )

dtLS.sim.factortype1[link %in% names(greek.label.factor), link.txt := factor(link,
                                                                             levels = names(greek.label.factor),
                                                                             labels = as.character(greek.label.factor))]


greek.label.lvm <- c(Y2 = expression(paste(nu[2], "   (intercept)")),
                     Y3 = expression(paste(nu[3], "   (intercept)")),
                     Y4 = expression(paste(nu[4], "   (intercept)")),
                     Y5 = expression(paste(nu[5], "   (intercept)")),
                     eta1 = expression(paste(eta[1], "   (intercept)")),
                     Z2 = expression(paste(Z[2], "   (intercept)")),
                     Z3 = expression(paste(Z[3], "   (intercept)")),
                     Z4 = expression(paste(Z[4], "   (intercept)")),
                     Z5 = expression(paste(Z[5], "   (intercept)")),
                     eta = expression(paste(eta[2], "   (intercept)")),
                     "Y1~Gene2Y" = expression(paste(k[1], "   (reg. OV~OV)")),
                     "Y2~eta1" = expression(paste(lambda[2], "   (loading)")),
                     "Y3~eta1" = expression(paste(lambda[3], "   (loading)")),
                     "Y4~eta1" = expression(paste(lambda[4], "   (loading)")),
                     "Y5~eta1" = expression(paste(lambda[5], "   (loading)")),
                     "eta1~eta2" = expression(paste(b[1], "   (reg. on LV~LV)")),
                     "eta1~Age" = expression(paste(gamma[2], "   (reg. OV~LV)")),
                     "eta1~Gene1Y" = expression(paste(gamma[1], "   (reg. OV~LV)")),
                     "Z2~eta1" = expression(paste(lambda[6], "   (loading)")),
                     "Z3~eta1" = expression(paste(lambda[7], "   (loading)")),
                     "Z4~eta1" = expression(paste(lambda[8], "   (loading)")),
                     "Z5~eta1" = expression(paste(lambda[9], "   (loading)")),
                     "eta2~Gender" = expression(paste(gamma[3], "   (reg. OV~LV)")),
                     "Y1~~Y2" = expression(paste(sigma[12], "   (covariance)"))
                     )

dtLS.sim.lvmtype1[link %in% names(greek.label.lvm), link.txt := factor(link,
                                                                       levels = names(greek.label.lvm),
                                                                       labels = as.character(greek.label.lvm))]

label.statistic <- c( 
    Ztest = "Gaussian approx.",
    Satt = "Satterthwaite approx.",
    SSC = "small sample correction",
    KR = "Satterthwaite approx. with small sample correction"
)
name.statistic <- names(label.statistic)
n.statistic <- length(name.statistic)

## * FUNCTION
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
    gg <- gg + geom_point(size = 4) + geom_line(size = 2)
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
groupFigures <- function(figure1, figure2, figure3, reduce.x = TRUE){
    x.factor <- c("20" = "black", "30" = "black", "50" = "black", "75" = "black", "100" = "black",
                  "150" = "transparent", "200" = "black", "300" = "transparent", "500" = "black")
    x.lvm <- c("20" = "black", "30" = "transparent", "50" = "black", "75" = "transparent", "100" = "black",
               "150" = "transparent", "200" = "black", "300" = "transparent", "500" = "black")

    figure1 <- figure1 + ggtitle("scenario (a): mixed model") + xlab("") + theme(text = element_text(size = 10), legend.position="none")
    figure2 <- figure2 + ggtitle("scenario (b): single factor model") + xlab("") + theme(text = element_text(size = 10), legend.position="none")
    if(reduce.x == TRUE){
        figure2 <- figure2 + theme(axis.text.x=element_text(color=x.factor))
    }
    figure3 <- figure3 + ggtitle("scenario (c): two latent variables model") + theme(text = element_text(size = 10))
    figure3 <- figure3 + guides(color=guide_legend(nrow=2,byrow=TRUE)) + theme(legend.margin=margin(t = -0.5, unit='cm'))
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

## * Figure 3
## unique(dtLS.sim.factortype1$link)
gg.mm <- createFigure(dtLS.sim.MMtype1,
                      robust = FALSE, link = c("Y2","eta~Gene1Y"),
                      vec.name = name.statistic,
                      vec.label = label.statistic)
gg.factor <- createFigure(dtLS.sim.factortype1,
                          robust = FALSE, link = c("Y2","Y4~eta","eta~Gene1Y","Y1~Gene2Y"),
                          vec.name = name.statistic,
                          vec.label = label.statistic)
gg.lvm <- createFigure(dtLS.sim.lvmtype1,
                       robust = FALSE, link = c("Y2","Y1~Gene2Y","Y4~eta1","eta1~Gene1Y","eta1~eta2","Y1~~Y2"),
                       vec.name = name.statistic,
                       vec.label = label.statistic)

## groupFigures(gg.mm,gg.factor,gg.lvm)

if(export){
    ## postscript(file.path(path.figures,"type1error-Wald.eps"), height = 9.5)
    pdf(file.path(path.figures,"type1error-Wald.pdf"), height = 9.25)
    groupFigures(gg.mm,gg.factor,gg.lvm)
    dev.off()
}

## * Figure 3 bis
ggCoef.mm <- createFigureBIS(data = dtLS.sim.MMtype1,
                             n = c(20,50),
                             robust = FALSE, link = NULL,
                             vec.name = name.statistic,
                             vec.label = label.statistic)
ggCoef.factor <- createFigureBIS(data = dtLS.sim.factortype1,
                                 n = c(20,50),
                                 robust = FALSE, link = NULL,
                                 vec.name = name.statistic,
                                 vec.label = label.statistic)
ggCoef.lvm <- createFigureBIS(data = dtLS.sim.lvmtype1,
                              n = c(20,50),
                              robust = FALSE, link = NULL,
                              vec.name = name.statistic,
                              vec.label = label.statistic)

if(export){
    pdf(file.path(path.figures,"type1error-Wald2.pdf"), height = 9.25)
    groupFigures(ggCoef.mm,ggCoef.factor,ggCoef.lvm, reduce.x = FALSE)
    dev.off()
}

## * Figure 4
ggR.mm <- createFigure(dtLS.sim.MMtype1,
                       robust = TRUE, link = c("Y2","eta~Gene1Y"),
                       vec.name = name.statistic,
                       vec.label = label.statistic)
ggR.factor <- createFigure(dtLS.sim.factortype1,
                           robust = TRUE, link = c("Y2","Y4~eta","eta~Gene1Y","Y1~Gene2Y"),
                           vec.name = name.statistic,
                           vec.label = label.statistic) + scale_y_continuous(breaks = scales::pretty_breaks(n = 5))
ggR.lvm <- createFigure(dtLS.sim.lvmtype1,
                        robust = TRUE, link = c("Y2","Y1~Gene2Y","Y4~eta1","eta1~Gene1Y","eta1~eta2","Y1~~Y2"),
                        vec.name = name.statistic,
                        vec.label = label.statistic) + scale_y_continuous(breaks = scales::pretty_breaks(n = 5))

if(export){
    ## postscript(file.path(path.figures,"type1error-Wald-robust.eps"))
    pdf(file.path(path.figures,"type1error-Wald-robust.pdf"), height = 9.25)
    groupFigures(ggR.mm,ggR.factor,ggR.lvm)
    dev.off()
}

##  groupFigures(ggR.mm,ggR.factor,ggR.lvm)

## * Figure 4 bis
ggCoefR.mm <- createFigureBIS(data = dtLS.sim.MMtype1,
                              n = c(20,50),
                              robust = TRUE, link = NULL,
                              vec.name = name.statistic,
                              vec.label = label.statistic)
ggCoefR.factor <- createFigureBIS(data = dtLS.sim.factortype1,
                                  n = c(20,50),
                                  robust = TRUE, link = NULL,
                                  vec.name = name.statistic,
                                  vec.label = label.statistic)
ggCoefR.lvm <- createFigureBIS(data = dtLS.sim.lvmtype1,
                               n = c(20,50),
                               robust = TRUE, link = NULL,
                               vec.name = name.statistic,
                               vec.label = label.statistic)

if(export){
    pdf(file.path(path.figures,"type1error-Wald2-robust.pdf"), height = 9.25)
    groupFigures(ggCoefR.mm,ggCoefR.factor,ggCoefR.lvm, reduce.x = FALSE)
    dev.off()
}

## * IV
dtLS.sim.IV <- readRDS(file.path(path.results,"type1error-simulation-IV.rds"))
##unique(dtLS.sim.IV$link)
dtLS.sim.IV[link == "Y2~X1", link := "Y1~Gene2Y"]
dtLS.sim.IV[link == "eta~GenderF", link := "eta~Gene1Y"]
dtLS.sim.IV[link == "eta~G", link := "eta~Age"]

dtLS.sim.IV[link %in% names(greek.label.factor), link.txt := factor(link,
                                                                    levels = names(greek.label.factor),
                                                                    labels = as.character(greek.label.factor))]
gg.IV <- createFigure(dtLS.sim.IV[correction == "Ztest"],
                      robust = FALSE, link = c("Y2","Y4~eta","eta~Gene1Y","Y1~Gene2Y"),
                      vec.name = name.statistic,
                      vec.label = label.statistic)

