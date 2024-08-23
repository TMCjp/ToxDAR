library(NormExpression)
library(ggplot2)

nonzeroRatio2AUCVC <- function (data, dataType = c("bk", "sc"), HG7 = NULL,
                                ERCC = NULL, TN = NULL, TC = NULL, CR = NULL, NR = NULL,
                                DESeq = NULL, UQ = NULL, TMM = NULL, TU = 0, GAPDH = NULL,
                                nonzeroRatio = NULL, cvNorm = TRUE, cvResolution = 0.005)
{
  nozeroIndex <- filteredZero(data, nonzeroRatio = nonzeroRatio)
  methodsList <- list(HG7 = HG7, ERCC = ERCC, TN = TN, TC = TC,
                      CR = CR, NR = NR, DESeq = DESeq, UQ = UQ, TMM = TMM,
                      TU = TU, GAPDH = GAPDH)
  specifiedMethods <- methodsList[!unlist(lapply(methodsList,
                                                 is.null))]
  if (length(TU) == 1 && TU == 0) {
    specifiedMethods$TU <- NULL
  }
  if (length(TU) == 1 && TU == 1) {
    if (dataType == "bk") {
      optimalPara <- optTU(data, nonzeroRatio = nonzeroRatio,
                           pre_ratio_range = c(1, 1), prResolution = 0.1,
                           lower_range = c(0.05, 0.4), upper_range = c(0.6,
                                                                       0.95), qResolution = 0.05, min_ubq = 1000,
                           cvNorm = cvNorm, cvResolution = cvResolution)
    }
    else {
      optimalPara <- optTU(data, nonzeroRatio = nonzeroRatio,
                           pre_ratio_range = c(0.2, 0.6), prResolution = 0.1,
                           lower_range = c(0.05, 0.4), upper_range = c(0.6,
                                                                       0.95), qResolution = 0.05, min_ubq = 100, cvNorm = cvNorm,
                           cvResolution = cvResolution)
    }
    optimalPara <- as.matrix(optimalPara)
    lower_trim <- optimalPara["lower", 1]
    upper_trim <- optimalPara["upper", 1]
    pre_ratio <- optimalPara["ratio", 1]
    para <- c(nonzeroRatio, pre_ratio, lower_trim, upper_trim)
    names(para)[1] <- "nonzeroRatio"
    paraMatrix <- t(as.matrix(para))
    colnames(paraMatrix) <- c("nonzeroRatio", "pre_ratio",
                              "lower_trim", "upper_trim")
    # write.table(paraMatrix, file = "bestPara.txt",
    #             sep = "\t", row.names = FALSE, col.names = FALSE,
    #             append = TRUE)
    TU.factors <- getFactors(data, method = "TU", lower_trim = lower_trim,
                             upper_trim = upper_trim, pre_ratio = pre_ratio, min_ubq = 100)
    norm.matrix <- getNormMatrix(data, TU.factors)
    dataUse2CV <- norm.matrix[nozeroIndex, ]
    cv.result <- getCV(dataUse2CV, cvNorm = cvNorm)
    TU.AUCVC <- CV2AUCVC(cv.result, cvResolution = cvResolution)
    specifiedMethods$TU <- NULL
  }
  numMethod <- length(specifiedMethods)
  if (numMethod >= 1) {
    method_range <- seq(1, numMethod, 1)
    for (i in method_range) {
      norm.matrix <- getNormMatrix(data, specifiedMethods[[i]])
      dataUse2CV <- norm.matrix[nozeroIndex, ]
      cv.result <- getCV(dataUse2CV, cvNorm = cvNorm)
      assign(names(specifiedMethods)[i], CV2AUCVC(cv.result,
                                                  cvResolution = cvResolution))
    }
    AUCVC.result <- NULL
    for (i in method_range) {
      AUCVC.result <- cbind(AUCVC.result, get(names(specifiedMethods)[i]))
    }
    colnames(AUCVC.result) <- names(specifiedMethods)
    if (length(TU) == 1 && TU == 1) {
      AUCVC.result <- cbind(AUCVC.result, TU.AUCVC)
      colnames(AUCVC.result) <- c(names(specifiedMethods),
                                  "TU")
    }
  }
  if (numMethod == 0 && TU == 0)
    stop("Please specify at least one method!")
  if (numMethod == 0 && TU == 1) {
    AUCVC.result <- as.matrix(TU.AUCVC)
    colnames(AUCVC.result) <- "TU"
  }
  return(rels <- list(auc = AUCVC.result, bp = paraMatrix))
}



gridAUCVC <- function (data, dataType = c("bk", "sc"), HG7 = NULL,
                       ERCC = NULL, TN = NULL, TC = NULL, CR = NULL, NR = NULL,
                       DESeq = NULL, UQ = NULL, TMM = NULL, TU = 0, GAPDH = NULL,
                       nonzeroRatios = c(0.7, 0.8, 0.9, 1), cvNorm = TRUE, cvResolution = 0.005)
{
  grid_result <- NULL
  # if (length(TU) == 1 && TU == 1) {
  #   colnames_paraMatrix <- c("nonzeroRatio", "pre_ratio",
  #                            "lower_trim", "upper_trim")
  #   # write.table(t(as.matrix(colnames_paraMatrix)), file = "bestPara.txt",
  #   #             sep = "\t", row.names = FALSE, col.names = FALSE)
  # }
  for (i in nonzeroRatios) {
    if (dataType == "sc") {
      if ((ncol(data) * i) <= 100) {
        cat("nonzeroRatio:", i, " is too small!\n")
        stop("We suggest that the minimal counts of nonzero samples should be greater than 100!")
      }
    }
    result <- nonzeroRatio2AUCVC(data = data, dataType = dataType,
                                 HG7 = HG7, ERCC = ERCC, TN = TN, TC = TC, CR = CR,
                                 NR = NR, DESeq = DESeq, UQ = UQ, TMM = TMM, TU = TU,
                                 GAPDH = GAPDH, nonzeroRatio = i, cvNorm = cvNorm,
                                 cvResolution = cvResolution)
    aucresult <- result$auc
    paraMatrix <- result$bp

    nonzeroM <- matrix(i, 1, 1, TRUE)
    colnames(nonzeroM) <- "NonzeroRatio"
    grid_record <- cbind(nonzeroM, aucresult)
    grid_result <- rbind(grid_result, grid_record)
  }
  return(rels <- list(grid=grid_result, bestPara=paraMatrix))
}


normfactors <- function(genexp, dataType="bk", method = c("sizefactor", "DESeq", "RLE", "UQ", "TMM", "TU")){

  if(method == "TU"){
    # Grid of non-zero ratios to produce AUCVCs for TU and it is time-consuming
    # nonzeroRatios can be set to 1 for bulk data to reduce computing time
    AUCVCs <- gridAUCVC(data= genexp, dataType=dataType, TU=1, nonzeroRatios= 1)
    # Find the parameters used by TU when it achieved the maximum AUCVC in bestPara.txt
    bestPara <- as.data.frame(AUCVCs$bestPara)

    # Using the parameters to produce the TU normalization factor
    facs <- getFactors(data = genexp, method = "TU",
                       pre_ratio = bestPara$pre_ratio,
                       lower_trim = bestPara$lower_trim, upper_trim = bestPara$upper_trim)
  } else {
    facs <- getFactors(data = genexp, method = method)
  }

  return(facs)
}





# ###
# # library(data.table)
# # genexp <- as.data.frame(fread('caseDATA/GSE121134_series_matrix.txt', skip = 60))
# # rownames(genexp) <- genexp$ID_REF
# # genexp <- as.matrix(genexp[, -1])
# # suppressWarnings(storage.mode(genexp) <- "numeric")
# load(file = 'RPackage/genexp.rdata')
#
# NR.norfac <- normfactors(genexp, dataType="bk", method = "sizefactor")
# TMM.norfac <- normfactors(genexp, dataType="bk", method = "TMM")
# DESeq.norfac <- normfactors(genexp, dataType="bk", method = "DESeq")
# TU.norfac <- normfactors(genexp, dataType="bk", method = "TU")
#
# # Produce the TU-normalized gene expression matrix
# norm.matrix <- getNormMatrix(data = genexp, norm.factors = NR.norfac)
#
#
# ### plot auc
# cv_uniform <- gatherCVs(data = genexp, nonzeroRatio = 1,
#                         NR = NR.norfac,
#                         DESeq = DESeq.norfac,
#                         TMM= TMM.norfac,
#                         TU= TU.norfac)
#
# cv_uniform <- as.data.frame(cv_uniform)
# cv_uniform$Cutoff <- as.numeric(cv_uniform$Cutoff)
# cv_uniform$Counts <- as.numeric(cv_uniform$Counts)
#
# # plot
# library(ggplot2)
# png(file = "cv.png", res=300, width=(1200*4.17), height=(960*4.17))
# plotCVs(cv_uniform, methods=c("NR", "DESeq", "TMM", "TU"),
#         legend.position=c(.85, .48))
# dev.off()


plotCVs <- function (data, methods = c("None", "HG7", "ERCC",
                            "TN", "TC", "CR", "NR", "DESeq",
                            "UQ", "TMM", "TU"), legend.position = c(1.12, 0.48))
{
  if (!is.data.frame(data))
    data <- data.frame(data)
  if (is.factor(data$Cutoff))
    data$Cutoff <- as.numeric(as.character(data$Cutoff))
  if (is.factor(data$Counts))
    data$Counts <- as.numeric(as.character(data$Counts))
  data$Methods <- factor(data$Methods, levels = methods, labels = methods)
  change_colours(ggplot(data = data, aes(x = Cutoff, y = Counts)) +
                   geom_line(aes(group = Methods, color = Methods), size = 3) +
                   xlab("Normalized CV cutoff") + ylab("Number of uniform genes") +
                   theme_bw() + theme(panel.grid.minor = element_blank(),
                                      axis.title.x = element_text(size = 48), axis.title.y = element_text(size = 48),
                                      axis.text.x = element_text(size = 38), axis.text.y = element_text(size = 38),
                                      legend.text = element_text(size = 39), legend.title = element_text(size = 43),
                                      legend.position = legend.position, legend.background = element_blank(),
                                      legend.key = element_blank(), legend.key.height = unit(1.8,
                                                                                             "cm"), plot.margin = unit(c(0.5, 7, 0.5,
                                                                                                                         0.5), "cm")) + scale_x_continuous(breaks = seq(0, 1, 0.2), expand = c(0, 0)) + scale_y_continuous(limits = c(NA, 1.05*max(data$Counts)), expand = c(0, 0)) + guides(color = guide_legend(title = NULL)),
                 c("olivedrab", "blue", "red", "violet",
                   "orange", "yellow", "magenta",
                   "peru", "black", "maroon", "lightblue",
                   "darkslateblue", "seashell4", "tan2",
                   "darkgreen", "springgreen"))
}



plotCors <- function (data, methods = c("None", "HG7", "ERCC",
                                        "TN", "TC", "CR", "NR", "DESeq",
                                        "UQ", "TMM", "TU"), legend.position = c(1.12,
                                                                                0.56))
{
  if (!is.data.frame(data))
    data <- data.frame(data)
  if (is.factor(data$Value))
    data$Value <- as.numeric(as.character(data$Value))
  data$Methods <- factor(data$Methods, levels = methods, labels = methods)
  change_colours(ggplot(data = data, aes(x = Value, y = (..count../sum(..count..))*100)) +
                   geom_freqpoly(aes(group = Methods, color = Methods),
                                 size = 3, bins = 50) + xlab("Spearman correlation coefficient") +
                   ylab("Fraction of gene pairs (%)") + theme_bw() + theme(panel.grid.minor = element_blank(),
                                                                       axis.title.x = element_text(size = 48), axis.title.y = element_text(size = 48),
                                                                       axis.text.x = element_text(size = 38), axis.text.y = element_text(size = 38),
                                                                       legend.text = element_text(size = 39), legend.title = element_text(size = 43),
                                                                       legend.position = legend.position, legend.background = element_blank(),
                                                                       legend.key = element_blank(), legend.key.height = unit(1.8,
                                                                                                                              "cm"), plot.margin = unit(c(0.5, 7, 0.5, 0.5),
                                                                                                                                                        "cm")) + scale_x_continuous(expand = c(0.01,
                                                                                                                                                                                               0.01), breaks = round(seq(-1, 1, 0.25), 2)) + scale_y_continuous(expand = c(0.01,
                                                                                                                                                                                                                                                                           0)) + guides(color = guide_legend(title = NULL)), c("olivedrab",
                                                                                                                                                                                                                                                                                                                               "blue", "red", "violet", "orange",
                                                                                                                                                                                                                                                                                                                               "yellow", "magenta", "peru", "black",
                                                                                                                                                                                                                                                                                                                               "maroon", "lightblue", "darkslateblue",
                                                                                                                                                                                                                                                                                                                               "seashell4", "tan2", "darkgreen", "springgreen"))
}
