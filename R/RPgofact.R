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





# ###
# diffgs <- na.omit(diffanno$Symbol[1:200])
# diffgs <- orthconvert(diffgs, org = "mmusculus")
#
#
# library(org.Hs.eg.db)
# genes <- mapIds(org.Hs.eg.db, keys=diffgs, keytype="SYMBOL", column="ENTREZID", multiVals="first")
#
# gomat <- gofact(genes, slim = 'go')
