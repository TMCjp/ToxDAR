library(ggord)
library(ggplot2)


# 生成随机颜色的函数
# generate_random_color <- function(n = 1) {
#   r <- sample(0:255, n) / 255
#   g <- sample(0:255, n) / 255
#   b <- sample(0:255, n) / 255
#   random_color <- sapply(1:n, function(x) rgb(r[x], g[x], b[x]))
#   return(random_color)
# }

toxPCA <- function(genexp, stype = NULL, outpng = TRUE){
  if(is.null(stype)){stype <- colnames(genexp)}

  exp_data <- cbind(type=stype, t(genexp))

  pca_data <- exp_data[, -1]
  suppressWarnings(storage.mode(pca_data) <- "numeric")

  ord <- prcomp(pca_data, center = T, scale. = F)


  if(outpng == TRUE){
    tys <- as.data.frame(exp_data)$type
    # rcol <- generate_random_color(n <- length(unique(tys)))

    if(length(which(table(tys) < 2)) == 0){el = TRUE} else {el = FALSE}

    png(file = "toxPCA.png", res=300, width=(1200*3), height=(960*3))
    #使用默认参数直接出图
    p <- suppressWarnings(
      ggord(ord, tys, #cols=rcol,
            size = 6, obslab = T,  poly = FALSE, polylntyp = "dashed", ellipse = el,
            ellipse_pro = 0.95, grp_title = "", arrow=NULL, txt = NULL)
    )
    #设定不同的分组形状
    p <- p + theme_classic() +
      theme(axis.text = element_text(size=30), axis.title = element_text(size=30),
            axis.line = element_line(size = 2), axis.ticks = element_blank())
    print(p)
    dev.off()
  }

  print(summary(ord))
  return(ord)
}


# load(file = 'RPackage/genexp.rdata')
# pcaresult <- toxPCA(genexp, stype = c('c', 'c', 'c', 't', 't', 't'), outpng = T)
# # pcaresult <- toxPCA(genexp, outpng = T)

