### R包 毒理机制网络分析及可视化
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(igraph)


netmoa <- function(toxid = 'D002995', ngid, ppiscore = 900){
  ### 1.毒物与基因关系
  load(paste0(system.file(package = "ToxDAR"), '/extdata/CTD_toxge.Rdata'))
  dag <- toxge
  # dag <- fread('RPackage/CTD_chem_gene_ixns.csv', sep = ',')
  dag <- dag[, c(2, 4, 5, 10, 11)]

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
  if(nrow(ekeggmat) > 5){ekeggmat <- ekeggmat[1:5, ]}

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



  ### 4.网络构建
  mathdag <- mathdag[, c(1,2,5)]
  colnames(mathdag) <- c('From', 'To', 'EV')
  mathgap <- cbind(mathgap, EV='')
  colnames(mathgap) <- c('From', 'To', 'EV')
  mathppi <- mathppi[, c(3,4,5)]
  colnames(mathppi) <- c('From', 'To', 'EV')
  net <- rbind( rbind(mathdag, mathgap), mathppi )


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
  lay <- matrix(NA, nrow = length(V(ig)), ncol = 2)

  lay[1, ] <- c(-100, 5)
  lay[-1, 2] <- runif(length(V(ig))-1, 0, 10)

  nn <- grep('^hsa', names(V(ig)))
  lay[nn, 1] <- 100

  lay[-c(1,nn), 1] <- runif(length(V(ig))-1-length(nn), -50, 50)

  netplot <- plot(ig, layout = lay)

  return(print(netplot))
}




# ###
# gis <- c("2","5243","8647","18669","5244","1244","10057","368","215","5825","10061","9619","64240","84836","10449","34","130013","48","641371","10965")
# netmat <- netmoa(toxid = 'D002995', ngid = gis, ppiscore = 900)
#
# netplot(netmoa = netmat)
