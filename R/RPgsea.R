### R包 毒理机制 敏感生物学通路分析及热图呈现
library(data.table)

### gsea function
gseafc <- function(geneList, term = "KEGG", minGSSize = 5, maxGSSize = 500){
  library(msigdbr)
  library(fgsea)

  m = msigdbr(species = "Homo sapiens", category = 'C2', subcategory = term)

  df.gene2go <- unique(data.frame(gene = m$entrez_gene, go = m$gs_id))
  term2gene <- with(df.gene2go, split(as.character(gene), as.character(go)))
  # m <- m[, c('gs_id', 'gs_cat', 'gs_subcat')]

  gseamat <- suppressWarnings(
    fgsea(pathways=term2gene,
          stats=geneList,
          minSize=minGSSize,
          maxSize=maxGSSize,
          gseaParam=1,
          eps=0,
          nproc = 0)
  )
  pwname <- unique(m[, c('gs_id', 'gs_exact_source', 'gs_description')])
  gseamat <- merge(pwname, gseamat, by.x='gs_id', by.y='pathway', all.y = T)

  gseamat <- gseamat[order(gseamat$ES, decreasing = T), -9]
  return(gseamat)
}



### plot ES result
plotgsea <- function(gsid, geneList, term = "KEGG", ticksSize = 0.2, lineSize = 2, axistext.size = 20, axistitle.size = 20, title.size = 30){
  m = msigdbr(species = "Homo sapiens", category = 'C2', subcategory = term)
  pwgs <- m$entrez_gene[which(m$gs_id == gsid)]
  pwnm <- m$gs_name[which(m$gs_id == gsid)]

  library(ggplot2)
  plotEnrichment(pwgs, geneList, ticksSize = ticksSize, lineSize = lineSize) + labs(title=pwnm, subtitle = "Gene set enrichment in gene expression data") +
    theme(axis.text = element_text(size=axistext.size), axis.title = element_text(size=axistitle.size), plot.title = element_text(size = title.size))
}



### pathway class
pathclass <- function(pwclass = c("Metabolism","Genetic Information Processing","Environmental Information Processing","Cellular Processes","Organismal Systems","Human Diseases","Drug Development")){
  library(rjson)
  library(magrittr)
  ko.in <- fromJSON(file = "https://www.kegg.jp/kegg-bin/download_htext?htext=br08901.keg&format=json&filedir=")

  ko <- data.frame()
  ko["K", "term"] <- "all"
  ko["K", "id"] <- "K"
  ko["K", "hie"] <- "K"
  ko["K", "level"] <- "0"

  #level 1
  for (i in 1:length(ko.in[["children"]])) {
    name <- ko.in[["children"]][[i]][["name"]]
    hie.1 <- ifelse(i<10,paste0("K", ".0", as.character(i)),paste0("K", ".", as.character(i)))
    ko[hie.1, "term"] <- name
    ko[hie.1, "hie"] <-hie.1
    ko[hie.1, "level"] <- 1
    ko[hie.1, "id"] <- hie.1
    for (j in 1:length(ko.in[["children"]][[i]][["children"]])) {
      name<-
        ko.in[["children"]][[i]][["children"]][[j]][["name"]]
      hie.2<-ifelse(j<10,paste0(hie.1,".0",j),paste0(hie.1,".",j))
      ko[hie.2, "hie"] <- hie.2
      ko[hie.2, "term"] <- name
      ko[hie.2, "level"] <- 2
      ko[hie.2, "id"] <- ko[hie.2, "hie"]
      for (k in 1:length(ko.in[["children"]][[i]][["children"]][[j]][["children"]])) {
        name <-
          ko.in[["children"]][[i]][["children"]][[j]][["children"]][[k]][["name"]]
        hie.3<-ifelse(k<10,paste0(hie.2,".0",k),paste0(hie.2,".",k))
        id.tmp <- unlist(strsplit(name, " "))
        ko[hie.3, "hie"] <- hie.3
        ko[hie.3, "id"] <- id.tmp[1]
        ko[hie.3, "term"] <-
          paste(id.tmp[-c(1, 2)], collapse = " ")
        ko[hie.3, "level"] <- 3
      }
    }
  }


  ### human kegg pathway
  hsapath <- read.table('https://rest.kegg.jp/list/pathway/hsa', sep = '\t', header = F)
  hsapath$V1 <- gsub('hsa', '', hsapath$V1)


  ### merge 获取具有层级信息的 hsa pathway
  hsako <- merge(ko, hsapath, by.x = 'id', by.y = 'V1', all = FALSE)

  kc <- ko$id[match(pwclass, ko$term)]
  nn <- unlist(lapply(kc, function(x) grep(x, hsako$hie)))
  hsako <- hsako[nn, -c(4,5)]
  hsako$id <- paste0('hsa', hsako$id)

  return(hsako)
}



### 归一化函数
interval <- function(vec, min=-2, max=2){
  k = (max-min)/(max(vec)-min(vec))
  vecmin = min(vec)

  for(i in 1:length(vec)){
    vec[i] <- min + k*(vec[i]-vecmin)
  }
  return(vec)
}



### ES heatmap
ESheatmap <- function(esmicro, Normalization=TRUE){
  if(Normalization){
    for(i in 1:nrow(esmicro)){
      nn <- which(is.na(esmicro[i, ]) == "TRUE")
      if(length(nn) > 0){esmicro[i, -nn] <- interval(esmicro[i, -nn])} else {esmicro[i, ] <- interval(esmicro[i, ])}
    }
    nan <- which(is.na(esmicro) == 'TRUE')
    esmicro[nan] <- 0
  }


  hsako <- fread(paste0(system.file(package = "ToxDAR"),"/extdata/pathwayclass.csv"))

  nn <- na.omit(match(hsako$id, rownames(esmicro)))
  esmicro <- as.matrix(esmicro[nn, ])

  hieinfo <- hsako[na.omit(match(rownames(esmicro), hsako$id)), ]$hie
  hieanno <- c(rep('Metabolism', length(grep('K.01', hieinfo))),
               rep('Genetic Information Processing', length(grep('K.02', hieinfo))),
               rep('Environmental Information Processing', length(grep('K.03', hieinfo))),
               rep('Cellular Processes', length(grep('K.04', hieinfo))),
               rep('Organismal System', length(grep('K.05', hieinfo))),
               rep('Human Diseases', length(grep('K.06', hieinfo)))
  )


  ### plot heatmap
  library(ComplexHeatmap)

  rownames(esmicro) <- hsako$term[match(rownames(esmicro), hsako$id)]
  # ### 抽提样本特征
  # Feature <- read.table('C:/Users/tmcjp/Desktop/DW1118/toxgroupinfo.csv', sep = ',', header = T)
  # colnames(Feature) <- c('samples','Stype')
  #
  # Feature <- Feature[order(Feature$Stype), ]
  # esmicro <- esmicro[, Feature$samples]
  #
  #
  #
  # ### 构建热图的注释条等
  # ha_mix_top = HeatmapAnnotation(anno1=Feature$Stype,
  #                                annotation_legend_param = list(anno1=list(title = 'SampleClass', nrow = 1)),
  #                                # col = list(anno1 = c("BLM" = "#6B7CBC", "BXQ" = "#FBAD1E", "JZQ" = "#F15F6C", "LJS" = "#94CCB1", "PQ" = "#9FC85E")),
  #                                # show_legend = T,
  #                                show_annotation_name = FALSE)
  row_ha = rowAnnotation(anno1 = hieanno, annotation_legend_param = list(anno1=list(title = 'PathwayClass', nrow = 2, direction = "vertical")),
                         show_legend = T, show_annotation_name = FALSE)


  color <- circlize::colorRamp2(c(min(esmicro), 0, max(esmicro)), c('#3C78B9', '#F2F2F2', "#DC143C"))

  ht <- Heatmap(esmicro,     #需要绘制热图的矩阵
                name = 'Enrichment scores',
                col = color,
                cluster_rows = F,
                cluster_columns = T,
                show_column_names = T,
                row_names_gp = gpar(fontsize = 7),
                column_names_gp = gpar(fontsize = 10, fontface = 'bold'),
                heatmap_legend_param = list(title_position = "leftcenter-rot", legend_height = unit(4, "cm"), gap = unit(5, "cm")),

                #heatmap的大小
                #width = unit(0.95, "npc"), #heatmap的宽度
                row_title = ' ',
                row_title_side = "right",

                left_annotation = row_ha
                # top_annotation = ha_mix_top,
                #oncoprint形式
                # rect_gp = gpar(col = "#FFFFFF", lwd = 5)
  )

  draw(ht, heatmap_legend_side = "left", annotation_legend_side = "right")
}





plotEnrichment <- function (pathway, stats, gseaParam = 1, ticksSize = 0.2, lineSize = 2)
{
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway,
                          returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  diff <- (max(tops) - min(bottoms))/8
  x = y = NULL
  g <- ggplot(toPlot, aes(x = x, y = y)) + geom_point(color = "green",
                                                      size = lineSize) + geom_hline(yintercept = max(tops), colour = "red",
                                                                               linetype = "dashed") + geom_hline(yintercept = min(bottoms),
                                                                                                                 colour = "red", linetype = "dashed") + geom_hline(yintercept = 0,
                                                                                                                                                                   colour = "black") + geom_line(color = "green", size = lineSize) +
    theme_bw() + geom_segment(data = data.frame(x = pathway),
                              mapping = aes(x = x, y = -diff/2, xend = x, yend = diff/2),
                              size = ticksSize) + theme(panel.border = element_blank(),
                                                        panel.grid.minor = element_blank()) + labs(x = "Rank",
                                                                                                   y = "Enrichment score")
  g
}




# ### GSEA analysis workflow
# # gsea analysis ES
# data(geneList, package="DOSE")
# gsearesult <- gseafc(geneList, term = "KEGG", minGSSize = 5, maxGSSize = 500)
#
#
# # plot gsea result
# gsid <- gsearesult$gs_id[1]
# plotgsea(gsid, geneList, term = "KEGG")
#
#
# # pc <- pathclass(pwclass = c("Metabolism","Genetic Information Processing","Environmental Information Processing","Cellular Processes"))
# # write.csv(pc, 'RPackage/pathwayclass.csv', row.names = F, quote = T)
#
#
# # plot ES heatmap
# esmicro <- as.matrix(gsearesult[, 'ES'])
# rownames(esmicro) <- gsearesult$gs_exact_source
#
# ESheatmap(esmicro, Normalization = F)
