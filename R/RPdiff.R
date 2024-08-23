library(limma)
library(gprofiler2)
library(ggplot2)
library(ggrepel)
library(DT)
library(htmlwidgets)
library(dplyr)


diffexp <- function(genexp, groupinfo, ntop = 1000){
  ### diff analysis
  # # 0831 设置差异方向，默认后者比前者
  gls <- unique(groupinfo)
  groupinfo <- plyr::mapvalues(groupinfo, from = gls, to = c('g1', 'g2'), warn_missing = F)

  group_list <- factor(groupinfo)
  design <- model.matrix(~group_list)
  colnames(design) <- levels(group_list)

  fit <- lmFit(genexp, design)
  fit <- eBayes(fit, trend=TRUE)
  result_limma <- topTable(fit, coef=2, n=Inf, sort.by = "P")
  result_limma <- cbind(ID=rownames(result_limma), result_limma)
  result_limma <- result_limma[1:ntop, ]

  return(result_limma)
}



diffanno <- function(diffresult, org = "hsapiens", toxid = 'C009277'){
  ### 标识符统一到symbol
  # https://biit.cs.ut.ee/gprofiler/page/organism-list
  # https://biit.cs.ut.ee/gprofiler/page/namespaces-list
  sym <- gorth(query = diffresult$ID,
               source_organism = org, target_organism = "hsapiens",
               numeric_ns = "ENTREZGENE_ACC", mthreshold = Inf, filter_na = TRUE)
  # sym <- gconvert(query = result_limma$ID, organism = org, numeric_ns = "ENTREZGENE_ACC",
  #                 target="ENTREZGENE", mthreshold = Inf, filter_na = TRUE)
  nn <- match(unique(sym$input), sym$input)
  sym <- sym[nn, c('input', 'ortholog_name')]
  colnames(sym) <- c('ID', 'Symbol')

  diffresult <- merge(sym, diffresult, by = 'ID', all.y = T)


  ### knowledge annotation
  # toxge <- fread('./RPackage/CTD_chem_gene_ixns.csv')
  load(paste0(system.file(package = "ToxDAR"), '/extdata/CTD_toxge.Rdata'))
  nn <- which(toxge$ChemicalID == toxid)
  toxge <- toxge[nn, c('GeneSymbol', 'PubMedIDs')]

  unig <- unique(toxge$GeneSymbol)
  pmid <- sapply(unig, function(x) gsub(', ', '|', toString(unique(unlist(strsplit(toxge$PubMedIDs[which(toxge$GeneSymbol == x)], '|', fixed = T))), fixed = T)))
  toxanno <- cbind(Symbol=unig, is.Evidence='Y', PMID=pmid)

  diffresult <- merge(diffresult, toxanno, by = 'Symbol', all.x = T)
  diffresult <- diffresult[order(diffresult$P.Value), ]


  ###
  return(diffresult)
}



volanno <- function(diffres, nodesize = 2.5, annosize = 5, themsize = 21, legend.position = c(0.8, 0.15), legend.size = 12){
  # Gather Log-fold change and FDR-corrected pvalues from DESeq2 results
  ## - Change pvalues to -log10 (1.3 = 0.05)
  data <- data.frame(gene = diffres$Symbol,
                     pval = -log10(diffres$P.Value),
                     lfc = diffres$logFC,
                     evi = diffres$is.Evidence)
  data$label <- ifelse(data$evi == 'Y', data$gene, "")

  # Remove any rows that have NA as an entry
  data <- na.omit(data)

  # Color the points which are up or down
  ## If fold-change > 0 and pvalue > 1.3 (Increased significant)
  ## If fold-change < 0 and pvalue > 1.3 (Decreased significant)
  data <- mutate(data, color = case_when(data$lfc > 0 & data$pval > 1.3 ~ "Significantly up-regulated",
                                         data$lfc < 0 & data$pval > 1.3 ~ "Significantly down-regulated",
                                         data$pval < 1.3 ~ "No significance"))
  data$color <- factor(data$color, levels = c("Significantly up-regulated", "Significantly down-regulated", "No significance"))

  # Make a basic ggplot2 object with x-y values
  vol <- ggplot(data, aes(x = lfc, y = pval, color = color))

  # Add ggplot2 layers
  p <- vol +
    # ggtitle(label = "Volcano Plot", subtitle = "Colored by fold-change direction") +
    geom_point(size = nodesize, alpha = 0.8, na.rm = T) +
    scale_color_manual(name = "",   #"Legend",
                       values = c('Significantly up-regulated' = "#CD4F39", 'Significantly down-regulated' = "#008B00", 'No significance' = "darkgray")) +
    theme_bw(base_size = themsize) + # change overall theme
    theme(legend.position = legend.position,
          legend.background = element_blank(),
          legend.text = element_text(size = legend.size)) + # change the legend
    xlab(expression(log[2]("Case" / "Control"))) + # Change X-Axis label
    ylab(expression(-log[10]("P-value"))) + # Change Y-Axis label   # adjusted
    geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add p-adj value cutoff line
    scale_y_continuous(trans = "log1p") # Scale yaxis due to large p-values

  p <- p + geom_text_repel(data = data, aes(x = lfc, y = pval, label = label),
                           size = annosize, box.padding = unit(0.2, "lines"),
                           point.padding = unit(0.1, "lines"),
                           segment.color = "black",
                           show.legend = FALSE)

  warning('The annotated information shown on the volcano map is limited, please refer to the results table for details')
  return(p)
}


IntTab <- function(diffres, outhtml = TRUE){
  nurl <- lapply(diffres$PMID, function(x) unlist(strsplit(x, '|', fixed = T)))
  for(i in 1:length(nurl)){
    if(!is.na(nurl[i])){
      nurl[[i]] <- sapply(nurl[[i]], function(x)paste0('<a href="https://pubmed.ncbi.nlm.nih.gov/', x, '" target="_blank">', x, '</a>'))
      nurl[[i]] <- gsub(', ', '|', toString(nurl[[i]]), fixed = T)
    }
  }
  diffres$PMID <- unlist(nurl)
  diffres$logFC <- round(diffres$logFC, 3)
  diffres$AveExpr <- round(diffres$AveExpr, 3)
  diffres$t <- round(diffres$t, 3)
  diffres$P.Value <- format(diffres$P.Value, scientific = TRUE, digits = 3)
  diffres$adj.P.Val <- format(diffres$adj.P.Val, scientific = TRUE, digits = 3)
  diffres$B <- round(diffres$B, 3)

  dt <- datatable(diffres, escape = FALSE, selection = 'none', rownames = FALSE,
                  options = list(pageLength = 20, dom = 'frtip', scrollX = T))

  ###
  if(outhtml){saveWidget(dt, "diffresult.html", selfcontained = TRUE)}
  return(dt)
}





# ### 差异分析
# load(file = 'RPackage/genexp.rdata')
# gi <- c(rep("control",3), rep("drug",3))
#
# diffanno <- diffexp(genexp, groupinfo = gi, org = "hsapiens", toxid = 'C009277')
# IntTab(diffres = diffanno, outhtml = T)
#
# volanno(diffres = diffanno)


