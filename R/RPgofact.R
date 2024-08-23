library(GOplot)

orthconvert <- function(genes, org = "mmusculus"){
  orthc <- gorth(query = genes, source_organism = org, target_organism = "hsapiens", filter_na = TRUE)
  return(orthc$ortholog_name)
}


gofact <- function(genes, slim = 'go'){

  ### loading
  if(slim == "GO"|slim == 'go'){
    load(paste0(system.file(package = "ToxDAR"), '/extdata/gofact/go.Rdata'))
    load(paste0(system.file(package = "ToxDAR"), '/extdata/gofact/Hs.goa.Rdata'))
    load(paste0(system.file(package = "ToxDAR"), '/extdata/gofact/GO-slim.Rdata'))

    #ontology root
    o2g = list(BP = "GO:0008150", CC = "GO:0005575", MF = "GO:0003674")
    u.o2g <- unlist(o2g)
  }
  if(slim == "KEGG"|slim == "kegg"){
    load(paste0(system.file(package = "ToxDAR"), '/extdata/gofact/ko.Rdata'))
    load(paste0(system.file(package = "ToxDAR"), '/extdata/gofact/Hs.koa.Rdata'))
    slim.all <- ko

    #ontology root
    o2g = list(KEGG = "K", kegg = "K")
    u.o2g <- unlist(o2g)
  }
  if(slim == "DO"|slim == "do"){
    load(paste0(system.file(package = "ToxDAR"), '/extdata/gofact/do.Rdata'))
    load(paste0(system.file(package = "ToxDAR"), '/extdata/gofact/Hs.doa.Rdata'))
    slim.all <- do

    #ontology root
    o2g = list(DO = "DOID:4", DO = "DOID:4")
    u.o2g <- unlist(o2g)
  }



  ### function ###
  ######
  result <- slim.all[order(slim.all$hie), ]
  rownames(result) <- result$id

  # result$all_genes <- u.t2g[result$id]
  result$nall <- u.t2Ng[result$id]
  result <- result[which(result$nall > 0),]


  ### 超几何分布计算
  Qterm2gene <-
    lapply(term2gene, function(x)
      x[x %in% genes])

  u.N.Qterm2gene <-
    unlist(lapply(Qterm2gene, length))

  Qterm2gene = Qterm2gene[which(u.N.Qterm2gene != 0)]

  u.N.Qterm2gene <- u.N.Qterm2gene[u.N.Qterm2gene != 0]

  Qterm2gene.id <-
    lapply(Qterm2gene, function(x)
      paste(x, collapse = ";"))

  u.Qterm2gene.id <- unlist(Qterm2gene.id)


  # 修改为dataframe格式 并重置行名
  data.table::setDF(df.geneid2genesym)
  rownames(df.geneid2genesym) <- df.geneid2genesym$id
  Qterm2gene.sym <-
    lapply(Qterm2gene, function(x)
      paste(df.geneid2genesym[x, 'sym'], collapse = ";"))

  u.Qterm2gene.sym <- unlist(Qterm2gene.sym)
  u.N.Qterm2gene <- u.N.Qterm2gene[u.N.Qterm2gene > 0]


  result <-
    result[result$id %in% names(u.N.Qterm2gene), ]

  result$npro <- u.N.Qterm2gene[result$id]
  result$genes.id <- u.Qterm2gene.id[result$id]
  result$genes.sym <- u.Qterm2gene.sym[result$id]

  result$Npro <-  result[u.o2g[result$ontology], 'npro']
  result$Nall <-  result[u.o2g[result$ontology], 'nall']


  result$ratio <- round(result$npro / result$Npro, digits = 4)
  result$bg.ratio <-
    round(result$nall / result$Nall, digits = 4)
  result$ER <- round(result$ratio / result$bg.ratio, digits = 4)
  result$exp <-
    round(result$nall / result$Nall * result$Npro, digits = 4)


  for (goid in result$id) {
    npro <- result[goid, 'npro']
    nall <- result[goid, 'nall']
    Nall <- result[goid, 'Nall']
    Npro <- result[goid, 'Npro']
    Exp <- result[goid, 'exp']
    ER <- result[goid, 'ER']
    p.h <-
      phyper(npro - 1, nall, Nall - nall, Npro, lower.tail = F)
    p.h <- round(p.h, digits = 4)
    if (ER >= 1) {
      result[goid, 'p.value'] <- p.h
    } else{
      result[goid, 'p.value'] <- 1 - p.h
    }
    if ((ER > 1) &
        (result[goid, 'p.value'] < 0.05)) {
      result[goid, 'P.direction'] = "++"
    } else if ((ER > 1) &
               (result[goid, 'p.value'] >= 0.05)) {
      result[goid, 'P.direction'] = "+"
    } else if ((ER < 1) &
               (result[goid, 'p.value'] < 0.05)) {
      result[goid, 'P.direction'] = "--"
    } else if ((ER < 1) &
               (result[goid, 'p.value'] >= 0.05)) {
      result[goid, 'P.direction'] = "-"
    } else{
      result[goid, 'P.direction'] = "0"
    }
  }

  p.adj <- p.adjust(result$p.value, method = 'BH')  #BH
  p.adj <- round(p.adj, digits = 4)
  result$p.value.adj <- p.adj

  result <- result[, c("hie","id","term","npro","ratio","ER","exp","p.value","P.direction","p.value.adj","nall","bg.ratio","Npro","Nall","genes.id","genes.sym")]
  ######


  return(result)
}


# 定义一个函数来实现首字母大写
capitalize_first_letter <- function(s) {
  # 获取第一个字母并转换为大写
  first_letter <- toupper(substr(s, 1, 1))
  # 获取除第一个字母外的其余部分
  rest <- substr(s, 2, nchar(s))
  # 返回大写的首字母加上其余部分
  return(paste(first_letter, rest, sep = ""))
}


plotfact <- function(enrichdf, ntop = 5, fillby = "logP", tn = 15, collimit = c(log(0.05), log(0)), xtextsize = 30, stripsize = 20){
  require(dplyr)

  tys <- unlist(lapply(enrichdf$hie, function(x) strsplit(x, '.', fixed = T)[[1]][1]))
  enrichdf$types <- ''
  nn <- which(tys == 'BP')
  if(length(nn) > 0){enrichdf$types[nn] <- 'Biological Process'}
  nn <- which(tys == 'CC')
  if(length(nn) > 0){enrichdf$types[nn] <- 'Cellular Component'}
  nn <- which(tys == 'MF')
  if(length(nn) > 0){enrichdf$types[nn] <- 'Molecular Function'}

  topER <- enrichdf %>%
    group_by(types) %>%       # 按照分类变量进行分组
    arrange(p.value) %>%      # 按照数值变量升序排序
    filter(row_number() <= ntop)    # 过滤掉排名在前5之后的行

  nn <- unlist(lapply(topER$term, function(x)nchar(x)))
  topER$term <- unlist(lapply(topER$term, function(x) capitalize_first_letter(x)))
  topER$term[which(nn > tn)] <- paste0(substr(topER$term[which(nn > tn)], 1, tn), '...')

  if(fillby == "logP"){
    p <- ggplot(topER, aes(term, npro, fill=log(p.value)))
  }
  if(fillby == "Pval"){
    p <- ggplot(topER, aes(term, npro, fill=p.value))
  }
  p <- p + geom_bar(position='dodge',stat="identity") + scale_fill_gradient(low = "#8B0000", high = "#FFC0CB", limits = collimit) +
    theme_bw() + ylab('Gene number') + xlab('')
  # p <- ggplot(topER, aes(id, npro, fill=p.value)) + geom_bar(position='dodge',stat="identity") + scale_fill_gradient(low = "#8B0000", high = "#FFC0CB") +
  #   theme_bw() + ylab('enrichment number') + xlab('')
  p <- p + facet_wrap(~types, scales='free_x') + theme(axis.text.x=element_text(angle=60,hjust=1, size = xtextsize), axis.title.x=NULL, strip.text = element_text(size = stripsize)) +
     theme(panel.grid =element_blank(), plot.margin=unit(c(0.5,0.5,0.5,5),'lines'))

  return(p)
}






plotexpfunc <- function(enrichdf, ntop = 5, diffresult, plottype = c('chord','cluster'), chord.space = 0.05, gene.space = 0.5, gene.size = 5, process.title = 27, process.label = 24, termcol = 3){
  require(dplyr)

  tys <- unlist(lapply(enrichdf$hie, function(x) strsplit(x, '.', fixed = T)[[1]][1]))
  enrichdf$hie <- tys

  enrichdf <- enrichdf %>%
    group_by(hie) %>%       # 按照分类变量进行分组
    arrange(p.value) %>%      # 按照数值变量升序排序
    filter(row_number() <= ntop)    # 过滤掉排名在前5之后的行

  # select dataframe
  egomat <- as.data.frame(enrichdf[, c("hie","id","term","genes.sym","p.value.adj")])
  egomat$genes.sym <- unlist(lapply(egomat$genes.sym, function(x) toString(strsplit(x, split = ';', fixed = T)[[1]]) ))
  colnames(egomat) <- c('category', 'ID', 'Term', 'Genes', 'adj_pval')


  # diff function genes
  funcgenes <- unique(strsplit(toString(egomat$Genes), split = ', ', fixed = T)[[1]])
  # limma result
  result_limma <- diffresult[, -c(2,9,10)]

  limmadat <- result_limma[match(funcgenes, result_limma[, 1]), ]
  colnames(limmadat)[1] <- "ID"


  # goplot
  circ <- circle_dat(egomat, limmadat)
  chordrs <- chord_dat(circ, limmadat, egomat$Term)

  # pdf('GOBPChord.pdf', width = 15, height = 15)
  if(plottype == 'chord'){
    p <- GOChord(chordrs, space = chord.space, gene.order = 'logFC',
                 gene.space = gene.space, gene.size = gene.size,
                 process.label = process.label, termcol = termcol)
    # 主图设置参考：https://ggplot2.tidyverse.org/reference/theme.html
    p <- p + theme(legend.title = element_text(size = process.title), legend.box="vertical") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#D3D3D3", size = 2) +
      annotate("text", x = -1.2, y = 1.2, label = "Gene expression", color = "blue", size = 15, fontface = "bold") +
      annotate("text", x = 1, y = 1.2, label = "Gene Function", color = "blue", size = 15, fontface = "bold")
  }
  if(plottype == 'cluster'){
    if(nrow(egomat) > 12){terms <- egomat$Term[1:12]} else {terms <- egomat$Term}

    p <- GOCluster(circ, terms, chord = chordrs, clust.by = 'logFC', term.width = 2)
  }


  return(p)
}







GOChord <- function (data, title, space, gene.order, gene.size, gene.space, termcol,
          nlfc = 1, lfc.col, lfc.min, lfc.max, ribbon.col, border.size,
          process.label, limit)
{
  y <- id <- xpro <- ypro <- xgen <- ygen <- lx <- ly <- ID <- logFC <- NULL
  Ncol <- dim(data)[2]
  if (missing(title))
    title <- ""
  if (missing(space))
    space = 0
  if (missing(gene.order))
    gene.order <- "none"
  if (missing(gene.size))
    gene.size <- 3
  if (missing(gene.space))
    gene.space <- 0.2
  if (missing(lfc.col))
    lfc.col <- c("brown1", "azure", "cornflowerblue")
  if (missing(lfc.min))
    lfc.min <- -3
  if (missing(lfc.max))
    lfc.max <- 3
  if (missing(border.size))
    border.size <- 0.5
  if (missing(process.label))
    process.label <- 11
  if (missing(limit))
    limit <- c(0, 0)
  if (gene.order == "logFC")
    data <- data[order(data[, Ncol], decreasing = T), ]
  if (gene.order == "alphabetical")
    data <- data[order(rownames(data)), ]
  if (sum(!is.na(match(colnames(data), "logFC"))) > 0) {
    if (nlfc == 1) {
      cdata <- GOplot:::check_chord(data[, 1:(Ncol - 1)], limit)
      lfc <- sapply(rownames(cdata), function(x) data[match(x,
                                                            rownames(data)), Ncol])
    }
    else {
      cdata <- GOplot:::check_chord(data[, 1:(Ncol - nlfc)], limit)
      lfc <- sapply(rownames(cdata), function(x) data[,
                                                      (Ncol - nlfc + 1)])
    }
  }
  else {
    cdata <- GOplot:::check_chord(data, limit)
    lfc <- 0
  }
  if (missing(ribbon.col))
    colRib <- grDevices::rainbow(dim(cdata)[2])
  else colRib <- ribbon.col
  nrib <- colSums(cdata)
  ngen <- rowSums(cdata)
  Ncol <- dim(cdata)[2]
  Nrow <- dim(cdata)[1]
  colRibb <- c()
  for (b in 1:length(nrib)) colRibb <- c(colRibb, rep(colRib[b],
                                                      202 * nrib[b]))
  r1 <- 1
  r2 <- r1 + 0.1
  xmax <- c()
  x <- 0
  for (r in 1:length(nrib)) {
    perc <- nrib[r]/sum(nrib)
    xmax <- c(xmax, (pi * perc) - space)
    if (length(x) <= Ncol - 1)
      x <- c(x, x[r] + pi * perc)
  }
  xp <- c()
  yp <- c()
  l <- 50
  for (s in 1:Ncol) {
    xh <- seq(x[s], x[s] + xmax[s], length = l)
    xp <- c(xp, r1 * sin(x[s]), r1 * sin(xh), r1 * sin(x[s] +
                                                         xmax[s]), r2 * sin(x[s] + xmax[s]), r2 * sin(rev(xh)),
            r2 * sin(x[s]))
    yp <- c(yp, r1 * cos(x[s]), r1 * cos(xh), r1 * cos(x[s] +
                                                         xmax[s]), r2 * cos(x[s] + xmax[s]), r2 * cos(rev(xh)),
            r2 * cos(x[s]))
  }
  df_process <- data.frame(x = xp, y = yp, id = rep(c(1:Ncol),
                                                    each = 4 + 2 * l))
  xp <- c()
  yp <- c()
  logs <- NULL
  x2 <- seq(0 - space, -pi - (-pi/Nrow) - space, length = Nrow)
  xmax2 <- rep(-pi/Nrow + space, length = Nrow)
  for (s in 1:Nrow) {
    xh <- seq(x2[s], x2[s] + xmax2[s], length = l)
    if (nlfc <= 1) {
      xp <- c(xp, (r1 + 0.05) * sin(x2[s]), (r1 + 0.05) *
                sin(xh), (r1 + 0.05) * sin(x2[s] + xmax2[s]),
              r2 * sin(x2[s] + xmax2[s]), r2 * sin(rev(xh)),
              r2 * sin(x2[s]))
      yp <- c(yp, (r1 + 0.05) * cos(x2[s]), (r1 + 0.05) *
                cos(xh), (r1 + 0.05) * cos(x2[s] + xmax2[s]),
              r2 * cos(x2[s] + xmax2[s]), r2 * cos(rev(xh)),
              r2 * cos(x2[s]))
    }
    else {
      tmp <- seq(r1, r2, length = nlfc + 1)
      for (t in 1:nlfc) {
        logs <- c(logs, data[s, (dim(data)[2] + 1 - t)])
        xp <- c(xp, (tmp[t]) * sin(x2[s]), (tmp[t]) *
                  sin(xh), (tmp[t]) * sin(x2[s] + xmax2[s]),
                tmp[t + 1] * sin(x2[s] + xmax2[s]), tmp[t +
                                                          1] * sin(rev(xh)), tmp[t + 1] * sin(x2[s]))
        yp <- c(yp, (tmp[t]) * cos(x2[s]), (tmp[t]) *
                  cos(xh), (tmp[t]) * cos(x2[s] + xmax2[s]),
                tmp[t + 1] * cos(x2[s] + xmax2[s]), tmp[t +
                                                          1] * cos(rev(xh)), tmp[t + 1] * cos(x2[s]))
      }
    }
  }
  if (lfc[1] != 0) {
    if (nlfc == 1) {
      df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:Nrow),
                                                      each = 4 + 2 * l), logFC = rep(lfc, each = 4 +
                                                                                       2 * l))
    }
    else {
      df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:(nlfc *
                                                             Nrow)), each = 4 + 2 * l), logFC = rep(logs,
                                                                                                    each = 4 + 2 * l))
    }
  }
  else {
    df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:Nrow),
                                                    each = 4 + 2 * l))
  }
  aseq <- seq(0, 180, length = length(x2))
  angle <- c()
  for (o in aseq) if ((o + 270) <= 360)
    angle <- c(angle, o + 270)
  else angle <- c(angle, o - 90)
  df_texg <- data.frame(xgen = (r1 + gene.space) * sin(x2 +
                                                         xmax2/2), ygen = (r1 + gene.space) * cos(x2 + xmax2/2),
                        labels = rownames(cdata), angle = angle)
  df_texp <- data.frame(xpro = (r1 + 0.15) * sin(x + xmax/2),
                        ypro = (r1 + 0.15) * cos(x + xmax/2), labels = colnames(cdata),
                        stringsAsFactors = FALSE)
  cols <- rep(colRib, each = 4 + 2 * l)
  x.end <- c()
  y.end <- c()
  processID <- c()
  for (gs in 1:length(x2)) {
    val <- seq(x2[gs], x2[gs] + xmax2[gs], length = ngen[gs] +
                 1)
    pros <- which((cdata[gs, ] != 0) == T)
    for (v in 1:(length(val) - 1)) {
      x.end <- c(x.end, sin(val[v]), sin(val[v + 1]))
      y.end <- c(y.end, cos(val[v]), cos(val[v + 1]))
      processID <- c(processID, rep(pros[v], 2))
    }
  }
  df_bezier <- data.frame(x.end = x.end, y.end = y.end, processID = processID)
  df_bezier <- df_bezier[order(df_bezier$processID, -df_bezier$y.end),
  ]
  x.start <- c()
  y.start <- c()
  for (rs in 1:length(x)) {
    val <- seq(x[rs], x[rs] + xmax[rs], length = nrib[rs] +
                 1)
    for (v in 1:(length(val) - 1)) {
      x.start <- c(x.start, sin(val[v]), sin(val[v + 1]))
      y.start <- c(y.start, cos(val[v]), cos(val[v + 1]))
    }
  }
  df_bezier$x.start <- x.start
  df_bezier$y.start <- y.start
  df_path <- bezier(df_bezier, colRib)
  if (length(df_genes$logFC) != 0) {
    tmp <- sapply(df_genes$logFC, function(x) ifelse(x >
                                                       lfc.max, lfc.max, x))
    logFC <- sapply(tmp, function(x) ifelse(x < lfc.min,
                                            lfc.min, x))
    df_genes$logFC <- logFC
  }
  g <- ggplot() + geom_polygon(data = df_process, aes(x, y,
                                                      group = id), fill = "gray70", inherit.aes = F,
                               color = "black") + geom_polygon(data = df_process,
                                                               aes(x, y, group = id), fill = cols, inherit.aes = F,
                                                               alpha = 0.6, color = "black") + geom_point(aes(x = xpro,
                                                                                                              y = ypro, size = factor(labels, levels = labels), shape = NA),
                                                                                                          data = df_texp) + guides(size = guide_legend("Terms",
                                                                                                                                                       ncol = termcol, byrow = T, override.aes = list(shape = 22,
                                                                                                                                                                                                fill = unique(cols), size = 8))) + theme(legend.text = element_text(size = process.label)) +
    geom_text(aes(xgen, ygen, label = labels, angle = angle),
              data = df_texg, size = gene.size) + geom_polygon(aes(x = lx,
                                                                   y = ly, group = ID), data = df_path, fill = colRibb,
                                                               color = "black", size = border.size, inherit.aes = F) +
    labs(title = title) + GOplot:::theme_blank
  if (nlfc >= 1) {
    g + geom_polygon(data = df_genes, aes(x, y, group = id,
                                          fill = logFC), inherit.aes = F, color = "black") +
      scale_fill_gradient2("logFC", space = "Lab",
                           low = lfc.col[3], mid = lfc.col[2], high = lfc.col[1],
                           guide = guide_colorbar(title.position = "top",
                                                  title.hjust = 0.5), breaks = c(min(df_genes$logFC),
                                                                                 max(df_genes$logFC)), labels = c(round(min(df_genes$logFC)),
                                                                                                                  round(max(df_genes$logFC)))) + theme(legend.position = "bottom",
                                                                                                                                                       legend.background = element_rect(fill = "transparent"),
                                                                                                                                                       legend.box = "horizontal", legend.direction = "horizontal")
  }
  else {
    g + geom_polygon(data = df_genes, aes(x, y, group = id),
                     fill = "gray50", inherit.aes = F, color = "black") +
      theme(legend.position = "bottom", legend.background = element_rect(fill = "transparent"),
            legend.box = "horizontal", legend.direction = "horizontal")
  }
}




GOCluster <- function (data, process, metric, chord, clust, clust.by, nlfc, lfc.col,
                       lfc.min, lfc.max, lfc.space, lfc.width, term.col, term.space,
                       term.width)
{
  x <- y <- xend <- yend <- width <- space <- logFC <- NULL
  if (missing(metric))
    metric <- "euclidean"
  if (missing(clust))
    clust <- "average"
  if (missing(clust.by))
    clust.by <- "term"
  if (missing(nlfc))
    nlfc <- 0
  if (missing(lfc.col))
    lfc.col <- c("firebrick1", "white", "dodgerblue")
  if (missing(lfc.min))
    lfc.min <- -3
  if (missing(lfc.max))
    lfc.max <- 3
  if (missing(lfc.space))
    lfc.space <- (-0.5)
  else lfc.space <- lfc.space * (-1)
  if (missing(lfc.width))
    lfc.width <- (-1.6)
  else lfc.width <- lfc.space - lfc.width - 0.1
  if (missing(term.col))
    term.col <- brewer.pal(length(process), "Set3")
  if (missing(term.space))
    term.space <- lfc.space + lfc.width
  else term.space <- term.space * (-1) + lfc.width
  if (missing(term.width))
    term.width <- 2 * lfc.width + term.space
  else term.width <- term.width * (-1) + term.space
  if (clust.by == "logFC")
    distance <- stats::dist(chord[, dim(chord)[2]], method = metric)
  if (clust.by == "term")
    distance <- stats::dist(chord, method = metric)
  cluster <- stats::hclust(distance, method = clust)
  dendr <- dendro_data(cluster)
  y_range <- range(dendr$segments$y)
  x_pos <- data.frame(x = dendr$label$x, label = as.character(dendr$label$label))
  chord <- as.data.frame(chord)
  chord$label <- as.character(rownames(chord))
  all <- merge(x_pos, chord, by = "label")
  all$label <- as.character(all$label)
  if (nlfc) {
    lfc_rect <- all[, c(2, dim(all)[2])]
    for (l in 4:dim(data)[2]) lfc_rect <- cbind(lfc_rect,
                                                sapply(all$label, function(x) data[match(x, data$genes),
                                                                                   l]))
    num <- dim(data)[2] - 1
    tmp <- seq(lfc.space, lfc.width, length = num)
    lfc <- data.frame(x = numeric(), width = numeric(), space = numeric(),
                      logFC = numeric())
    for (l in 1:(length(tmp) - 1)) {
      tmp_df <- data.frame(x = lfc_rect[, 1], width = tmp[l +
                                                            1], space = tmp[l], logFC = lfc_rect[, l + 1])
      lfc <- rbind(lfc, tmp_df)
    }
  }
  else {
    lfc <- all[, c(2, dim(all)[2])]
    lfc$space <- lfc.space
    lfc$width <- lfc.width
  }
  term <- all[, c(2:(length(process) + 2))]
  color <- NULL
  termx <- NULL
  tspace <- NULL
  twidth <- NULL
  for (row in 1:dim(term)[1]) {
    idx <- which(term[row, -1] != 0)
    if (length(idx) != 0) {
      termx <- c(termx, rep(term[row, 1], length(idx)))
      color <- c(color, term.col[idx])
      tmp <- seq(term.space, term.width, length = length(idx) +
                   1)
      tspace <- c(tspace, tmp[1:(length(tmp) - 1)])
      twidth <- c(twidth, tmp[2:length(tmp)])
    }
  }
  tmp <- sapply(lfc$logFC, function(x) ifelse(x > lfc.max,
                                              lfc.max, x))
  logFC <- sapply(tmp, function(x) ifelse(x < lfc.min, lfc.min,
                                          x))
  lfc$logFC <- logFC
  term_rect <- data.frame(x = termx, width = twidth, space = tspace,
                          col = color)
  legend <- data.frame(x = 1:length(process), label = process)
  ggplot() + geom_segment(data = segment(dendr), aes(x = x,
                                                     y = y, xend = xend, yend = yend)) + geom_rect(data = lfc,
                                                                                                   aes(xmin = x - 0.5, xmax = x + 0.5, ymin = width, ymax = space,
                                                                                                       fill = logFC)) + scale_fill_gradient2("logFC",
                                                                                                                                             space = "Lab", low = lfc.col[3], mid = lfc.col[2],
                                                                                                                                             high = lfc.col[1], guide = guide_colorbar(title.position = "top",
                                                                                                                                                                                       title.hjust = 0.5), breaks = c(min(lfc$logFC), max(lfc$logFC)),
                                                                                                                                             labels = c(round(min(lfc$logFC)), round(max(lfc$logFC)))) +
    geom_rect(data = term_rect, aes(xmin = x - 0.5, xmax = x +
                                      0.5, ymin = width, ymax = space), fill = term_rect$col) +
    geom_point(data = legend, aes(x = x, y = 0.1, size = factor(label,
                                                                levels = label), shape = NA)) + guides(size = guide_legend("GO Terms",
                                                                                                                           ncol = 4, byrow = T, override.aes = list(shape = 22,
                                                                                                                                                                    fill = term.col, size = 8))) + coord_polar() + scale_y_reverse() +
    theme(legend.position = "bottom", legend.background = element_rect(fill = "transparent"),
          legend.box = "horizontal", legend.direction = "horizontal") +
    GOplot:::theme_blank
}





bezier <- function(data, process.col)
{
  x <- c()
  y <- c()
  Id <- c()
  sequ <- seq(0, 1, by = 0.01)
  N <- dim(data)[1]
  sN <- seq(1, N, by = 2)
  if (process.col[1] == "")
    col_rain <- grDevices::rainbow(N)
  else col_rain <- process.col
  for (n in sN) {
    xval <- c()
    xval2 <- c()
    yval <- c()
    yval2 <- c()
    for (t in sequ) {
      xva <- (1 - t) * (1 - t) * data$x.start[n] + t *
        t * data$x.end[n]
      xval <- c(xval, xva)
      xva2 <- (1 - t) * (1 - t) * data$x.start[n + 1] +
        t * t * data$x.end[n + 1]
      xval2 <- c(xval2, xva2)
      yva <- (1 - t) * (1 - t) * data$y.start[n] + t *
        t * data$y.end[n]
      yval <- c(yval, yva)
      yva2 <- (1 - t) * (1 - t) * data$y.start[n + 1] +
        t * t * data$y.end[n + 1]
      yval2 <- c(yval2, yva2)
    }
    x <- c(x, xval, rev(xval2))
    y <- c(y, yval, rev(yval2))
    Id <- c(Id, rep(n, 2 * length(sequ)))
  }
  df <- data.frame(lx = x, ly = y, ID = Id)
  return(df)
}




# ###
# diffgs <- na.omit(diffanno$Symbol[1:200])
# diffgs <- orthconvert(diffgs, org = "mmusculus")
#
#
# library(org.Hs.eg.db)
# genes <- mapIds(org.Hs.eg.db, keys=diffgs, keytype="SYMBOL", column="ENTREZID", multiVals="first")
#
# gomat <- gofact(genes, slim = 'go')
