library(ggord)
library(ggpubr)
library(ggplot2)


# 生成随机颜色的函数
# generate_random_color <- function(n = 1) {
#   r <- sample(0:255, n) / 255
#   g <- sample(0:255, n) / 255
#   b <- sample(0:255, n) / 255
#   random_color <- sapply(1:n, function(x) rgb(r[x], g[x], b[x]))
#   return(random_color)
# }

toxPCA <- function(genexp, stype = NULL, outpng = TRUE, legendtitle = 'Group'){
  if(is.null(stype)){stype <- colnames(genexp)}

  exp_data <- cbind(type=stype, t(genexp))

  pca_data <- exp_data[, -1]
  suppressWarnings(storage.mode(pca_data) <- "numeric")

  ord <- prcomp(pca_data, center = T, scale. = F)


  if(outpng == TRUE){
    tys <- as.data.frame(exp_data)$type
    #20240609 tys转为factor，且与stype顺序一致，保障图例的排序
    levels <- unique(stype)
    tys <- factor(tys, levels = levels)
    # rcol <- generate_random_color(n <- length(unique(tys)))

    if(length(which(table(tys) < 2)) == 0){el = TRUE} else {el = FALSE}

    png(file = "toxPCA.png", res=300, width=(1200*3), height=(960*3))
    #使用默认参数直接出图
    p <- suppressWarnings(
      ggord(ord, tys, #cols=rcol,
            size = 6, obslab = F,  poly = FALSE, polylntyp = "dashed", ellipse = el,
            ellipse_pro = 0.95, grp_title = legendtitle, arrow=NULL, txt = NULL)
    )
    #设定不同的分组形状
    # p <- p + theme_classic() +
    p <- p + theme_bw() +
      theme(axis.text = element_text(size=30), axis.title = element_text(size=30),
            legend.title = element_text(size=28),
            legend.text = element_text(size=18),
            legend.key.size = unit(28, "pt"),
            # axis.line = element_line(size = 2),
            axis.ticks = element_blank(),
            panel.border = element_rect(size = 2))
    print(p)
    dev.off()
  }

  print(summary(ord))
  return(ord)
}


# load(file = 'RPackage/genexp.rdata')
# pcaresult <- toxPCA(genexp, stype = c('c', 'c', 'c', 't', 't', 't'), outpng = T)
# # pcaresult <- toxPCA(genexp, outpng = T)


plotVio <- function(genexp, stype = NULL, xlabs = 'Group'){
  if(is.null(stype)){stype <- colnames(genexp)}

  exp_data <- as.data.frame(cbind(type=stype, t(genexp)))
  viomat <- reshape2::melt(exp_data, id.vars='type')

  levels <- unique(stype)
  viomat$type <- factor(viomat$type, levels = levels)
  viomat$value <- as.numeric(viomat$value)


  p <- ggviolin(viomat, x = "type", y = "value", fill = "type", xlab = xlabs, alpha = 0.8, size = 1,
                add = "boxplot", error.plot = "crossbar", trim = FALSE, ggtheme = theme_bw()) +
    theme(axis.text = element_text(size=20), axis.title = element_text(size=30),
          legend.position = "none",
          # axis.line = element_line(size = 2),
          axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5),
          panel.border = element_rect(size = 2))
  ### legend position
  # p <- p + theme(legend.position=c(0.95,0.95), legend.justification=c(0,1))

  return(p)
}
