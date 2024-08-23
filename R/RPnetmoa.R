### R包 毒理机制网络分析及可视化
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(igraph)


netmoa <- function(toxid = 'D002995', ngid, toppathway = 5, ppiscore = 900){
  ### 1.毒物与基因关系
  load(paste0(system.file(package = "ToxDAR"), '/extdata/CTD_toxge.Rdata'))
  dag <- toxge
  # dag <- fread('RPackage/CTD_chem_gene_ixns.csv', sep = ',')
  # dag <- dag[, c(2, 4, 5, 10, 11)]
  dag <- dag[, c(1, 2, 3, 4, 8)]

  # Clofibric Acid  D002995
  nn <- which(dag$ChemicalID == toxid)
  mathdag <- dag[nn, ]

  nn <- na.omit(match(ngid, mathdag$GeneID))
  mathdag <- mathdag[nn, ]



  ### 2.基因与通路关系
  ekegg <- enrichKEGG(gene = ngid,
                      organism = "hsa",
                      keyType = "kegg",
                      minGSSize = 1,
                      # universe      = names(geneList),
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 1)
  ekeggmat <- ekegg@result
  if(nrow(ekeggmat) > toppathway){ekeggmat <- ekeggmat[1:toppathway, ]}

  eklist <- lapply(ekeggmat$geneID, function(x) strsplit(x, split = '/', fixed = T)[[1]])
  symbols <- mapIds(org.Hs.eg.db, keys=unlist(eklist), keytype="ENTREZID", column="SYMBOL", multiVals="first")
  mathgap <- as.data.frame(symbols)

  path <- unlist( lapply(1:length(eklist), function(x) rep(ekeggmat$ID[x], length(eklist[[x]]))) )
  mathgap <- cbind(path, mathgap)
  colnames(mathgap) <- c('path', 'GeneSymbol')



  ### 3.PPI网络
  load(paste0(system.file(package = "ToxDAR"), '/extdata/STRING_ppi.Rdata'))
  # ppi <- fread('RPackage/PPIdata.csv', sep = ',')
  nn <- which(ppi$combined_score > ppiscore)
  ppi <- ppi[nn, ]

  nf <- unlist(lapply(ngid, function(x) which(ppi$protein1 == x)))
  nt <- unlist(lapply(ngid, function(x) which(ppi$protein2 == x)))
  nn <- intersect(nf, nt)

  mathppi <- ppi[nn, ]
  # 0702 add ev annotation
  mathppi$combined_score <- paste0('STRING_combined_score:', mathppi$combined_score)



  ### 4.网络构建
  keygenes <- union(mathppi$protein1_symbol, mathppi$protein2_symbol)
  mathdag <- mathdag[mathdag$GeneSymbol %in% keygenes]
  # 0702 add ev annotation
  mathdag$PubMedIDs <- paste0('PubMedID:', mathdag$PubMedIDs)


  rows_to_keep <- which(mathgap$GeneSymbol %in% keygenes)
  mathgap <- mathgap[rows_to_keep, ]


  mathdag <- mathdag[, c(1,3,5)]
  colnames(mathdag) <- c('From', 'To', 'EV')
  mathgap <- cbind(mathgap, EV='')
  colnames(mathgap) <- c('From', 'To', 'EV')
  mathppi <- mathppi[, c(3,4,5)]
  colnames(mathppi) <- c('From', 'To', 'EV')
  net <- rbind( rbind(mathdag, mathgap), mathppi )


  return(net)
}





netfunc <- function(ngid, ppiscore = 300){
  ### 1.基因与通路关系
  ekegg <- enrichKEGG(gene = ngid,
                      organism = "hsa",
                      keyType = "kegg",
                      minGSSize = 1,
                      # universe      = names(geneList),
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 1)
  ekeggmat <- ekegg@result
  if(nrow(ekeggmat) > 5){ekeggmat <- ekeggmat[1:5, ]}

  eklist <- lapply(ekeggmat$geneID, function(x) strsplit(x, split = '/', fixed = T)[[1]])
  symbols <- mapIds(org.Hs.eg.db, keys=unlist(eklist), keytype="ENTREZID", column="SYMBOL", multiVals="first")
  mathgap <- as.data.frame(symbols)

  path <- unlist( lapply(1:length(eklist), function(x) rep(ekeggmat$ID[x], length(eklist[[x]]))) )
  mathgap <- cbind(path, mathgap)
  colnames(mathgap) <- c('path', 'GeneSymbol')


  ### 2.PPI网络
  load(paste0(system.file(package = "ToxDAR"), '/extdata/STRING_ppi.Rdata'))
  # ppi <- fread('RPackage/PPIdata.csv', sep = ',')
  nn <- which(ppi$combined_score > ppiscore)
  ppi <- ppi[nn, ]

  nf <- unlist(lapply(ngid, function(x) which(ppi$protein1 == x)))
  nt <- unlist(lapply(ngid, function(x) which(ppi$protein2 == x)))
  nn <- intersect(nf, nt)

  mathppi <- ppi[nn, ]


  ### 3.网络构建
  colnames(mathgap) <- c('From', 'To')
  mathppi <- mathppi[, c(3,4)]
  colnames(mathppi) <- c('From', 'To')
  net <- rbind(mathgap, mathppi)


  return(net)
}






netplot <- function(netmoa){
  ig <- graph_from_edgelist(as.matrix(netmoa[, c(1,2)]), directed = F)

  ### 5.网络分析&可视化
  # sub_g <- cluster_walktrap(ig)
  # submat <- as.data.frame(cbind(names=sub_g$names, cluster=sub_g$membership))
  # submat$cluster <- as.numeric(submat$cluster)
  # submat <- submat[order(submat$cluster, decreasing = F), ]
  # plot(sub_g, ig)

  # lay <- layout_as_tree(ig)
  ## 布局、形状、颜色
  lay <- matrix(NA, nrow = length(V(ig)), ncol = 2)
  tys <- character(length(V(ig)))
  col <- character(length(V(ig)))

  lay[1, ] <- c(-100, 5)
  lay[-1, 2] <- runif(length(V(ig))-1, 0, 10)
  tys[1] <- 'raster'

  nn <- grep('^hsa', names(V(ig)))
  lay[nn, 1] <- 100
  tys[nn] <- 'vrectangle'
  col[nn] <- '#707DDB'

  lay[-c(1,nn), 1] <- runif(length(V(ig))-1-length(nn), -50, 50)



  tys[which(tys == '')] <- 'circle'
  col[which(col == '')] <- '#7EB897'



  ### ID转通路名称
  nn <- grep('^hsa', netmoa$From)

  library(KEGGREST)
  kegg_ids <- netmoa$From[nn]
  kegg_ids <- unique(kegg_ids)


  ### 2024.08.22
  if(length(kegg_ids) > 10){

    # 使用lapply将字符串向量按每10个一组进行处理
    processed <- lapply(seq(1, length(kegg_ids), by = 10), function(i) {
      end <- min(i + 9, length(kegg_ids))
      keggGet(kegg_ids[i:end])
    })

    for(nn in 1:(length(processed)-1)){
      yy <- nn + 1
      processed[[yy]] <- c(processed[[nn]], processed[[yy]])
    }
    kegg_info <- processed[[length(processed)]]

  } else {
    kegg_info <- keggGet(kegg_ids)
  }

  kegg_info <- unlist(lapply(kegg_info, function(x) x$NAME[1]))
  kegg_info <- gsub(' - Homo sapiens (human)', '', kegg_info, fixed = T)

  library(plyr)
  netmoa$From <- mapvalues(netmoa$From, kegg_ids, kegg_info, warn_missing = F)

  ## 根据映射完的通路名称重新生成网络
  ig <- graph_from_edgelist(as.matrix(netmoa[, c(1,2)]), directed = F)

  node_shapes <- tys
  names(node_shapes) <- V(ig)$name


  # netplot <- plot(ig, layout = lay, vertex.shape = node_shapes, vertex.color = col, vertex.frame.color = NA,
  #                 vertex.label.degree = pi/4, vertex.label.dist = 2)
  # return(print(netplot))

  plot(ig, layout = lay, vertex.shape = node_shapes, vertex.color = col, vertex.frame.color = NA,
       vertex.label.degree = pi/4, vertex.label.dist = 2)
}




# ###
# gis <- c("2","5243","8647","18669","5244","1244","10057","368","215","5825","10061","9619","64240","84836","10449","34","130013","48","641371","10965")
# netmat <- netmoa(toxid = 'D002995', ngid = gis, ppiscore = 900)
#
# netplot(netmoa = netmat)





### 0304 add custom layout
customlay <- function(graph, coord) {
  ### 更新坐标
  graph$x$nodes$x <- coord[, 1]
  graph$x$nodes$y <- coord[, 2]

  # https://search.r-project.org/CRAN/refmans/visNetwork/html/visIgraphLayout.html
  graph$x$igraphlayout <- list(type = "square")
  # graph %>% visNodes(physics = FALSE) %>% visEdges(smooth = FALSE) %>%
  #   visPhysics(stabilization = FALSE)
  graph
}



netplothtml <- function(netmoa, outtype = 'pdf'){
  ### 构建网络
  ig <- igraph::graph_from_edgelist(as.matrix(netmoa[, 1:2]), directed = F)


  # 计算布局坐标
  lay <- matrix(NA, nrow = length(V(ig)), ncol = 2)
  lay[1, ] <- c(-1, 0.5)
  lay[-1, 2] <- runif(length(V(ig))-1, 0, 1)
  nn <- grep('^hsa', names(V(ig)))
  lay[nn, 1] <- 1
  lay[-c(1,nn), 1] <- runif(length(V(ig))-1-length(nn), -0.5, 0.5)

  # 计算节点形状
  tys <- character(length(V(ig)))
  tys[1] <- 'diamond'
  nn <- grep('^hsa', names(V(ig)))
  tys[nn] <- 'square'

  # 计算节点颜色
  col <- character(length(V(ig)))
  col[1] <- '#F4B183'
  col[nn] <- '#707DDB'

  tys[which(tys == '')] <- 'dot'
  col[which(col == '')] <- '#7EB897'



  ### ID转通路名称
  nn <- grep('^hsa', netmoa$From)

  library(KEGGREST)
  kegg_ids <- netmoa$From[nn]
  kegg_ids <- unique(kegg_ids)

  ### 2024.08.22
  if(length(kegg_ids) > 10){

    # 使用lapply将字符串向量按每10个一组进行处理
    processed <- lapply(seq(1, length(kegg_ids), by = 10), function(i) {
      end <- min(i + 9, length(kegg_ids))
      keggGet(kegg_ids[i:end])
    })

    for(nn in 1:(length(processed)-1)){
      yy <- nn + 1
      processed[[yy]] <- c(processed[[nn]], processed[[yy]])
    }
    kegg_info <- processed[[length(processed)]]

  } else {
    kegg_info <- keggGet(kegg_ids)
  }

  kegg_info <- unlist(lapply(kegg_info, function(x) x$NAME[1]))
  kegg_info <- gsub(' - Homo sapiens (human)', '', kegg_info, fixed = T)

  library(plyr)
  netmoa$From <- mapvalues(netmoa$From, kegg_ids, kegg_info, warn_missing = F)

  ## 根据映射完的通路名称重新生成网络
  ig <- graph_from_edgelist(as.matrix(netmoa[, 1:2]), directed = F)



  ### network html
  library(visNetwork)
  visGraph <- toVisNetworkData(ig)
  visGraph$nodes$color <- col
  visGraph$nodes$shape <- tys

  vg <- visNetwork(visGraph$nodes, visGraph$edges, height = 800, width = 1200) %>%   #, main = mt, submain = st) %>%
    visNodes( scaling = list(label = list(enabled = T, min = 16, max = 18)),
              color = list(background = visGraph$nodes$color, border = NULL),
              shape = visGraph$nodes$shape) %>%
    visEdges(color = '#D3D3D3', width = 3)

  vg <- vg %>% customlay(coord = lay)

  vg <- vg %>%
    # visNodes( scaling = list(label = list(enabled = T, min = 8, max = 15)) ) %>%
    # visEdges(smooth = list(type = "curvedCW")) %>%
    visOptions(manipulation = TRUE, highlightNearest = FALSE, nodesIdSelection = list(enabled = FALSE)) %>%
    # visLegend(addNodes = lnodes, addEdges = ledges, position = "right", useGroups = FALSE, ncol = 1) %>%
    visInteraction(dragNodes = TRUE, dragView = TRUE, zoomView = TRUE, navigationButtons = FALSE) %>%
    visPhysics(maxVelocity = 0.1)

  if(outtype == 'png'){
    vg %>%
      visExport(type = "png", label = 'Download PNG', style= "color: #ffffff;background: #3C90BC;border:none;")
  } else {
    vg %>%
      visExport(type = "pdf", label = 'Download PDF', style= "color: #ffffff;background: #3C90BC;border:none;")
  }

}



###
savehtml <- function(network, htmlfile){
  visSave(network, file = htmlfile)
}

