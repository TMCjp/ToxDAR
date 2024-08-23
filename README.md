ToxDAR
================

A Workflow Software for Analyzing Toxicologically Relevant Proteomic and Transcriptomic Data, From Data Preprocessing to Toxicological Mechanism Elucidation.

Overview
--------

Exploration of toxicological mechanisms is imperative for the assessment of potential adverse reactions to chemicals and pharmaceutical agents, the engineering of safer compounds, and the preservation of public health.

High-throughput proteomics and transcriptomics can accurately capture the body's response to toxins and have become key tools for revealing complex toxicological mechanisms. Recently, a vast amount of omics data related to toxicological mechanisms has been accumulated. However, analyzing and utilizing this data remains a major challenge for researchers.

To address this, the current study has developed ToxDAR, a workflow-oriented R package for preprocessing and analyzing toxicological multi-omics data. ToxDAR integrates packages like NormExpression, DESeq2, and igraph, and utilizes R functions such as prcomp and phyper. It supports data preprocessing, quality control, differential expression analysis, functional analysis, and network analysis. ToxDAR’s architecture also includes five major categories of mechanism-related biological entities and details fifteen types of interactions among them, providing comprehensive knowledge annotation for omics data analysis results. 


Features
--------

**ToxDAR, designed to enhance the capabilities of data analysis and annotation, thereby improving the efficiency and quality of toxicological mechanism elucidation.**

-   ToxDAR is consisted of four main modules:
    - Data Preparation
    - Quality Control
    - Data Analysis
    - Data Interpretation

<img src = "vignettes/framework.jpg" width = "80%" align = "center" />
<br></br>

**The ToxDAR software package extracts associations and supporting evidence between toxicants and biological entities such as genes, pathways, and phenotypes from multiple knowledge bases.**

<table style="width:80%;">
<colgroup>
<col width="16%" />
<col width="20%" />
<col width="19%" />
</colgroup>
<thead>
<tr class="header">
<th>Source</th>
<th>Version</th>
<th>URL</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>CTD</td>
<td>v2021-10</td>
<td>http://ctdbase.org/</td>
</tr>
<tr class="even">
<td>AOP-Wiki</td>
<td>v2022-12</td>
<td>https://aopwiki.org/</td>
</tr>
<tr class="odd">
<td>KEGG</td>
<td>v0.7.2</td>
<td>https://www.kegg.jp/</td>
</tr>
<tr class="even">
<td>DrugBank</td>
<td>v2020-12-15</td>
<td>https://go.drugbank.com/</td>
</tr>
<tr class="odd">
<td>DisGeNet</td>
<td>v7.0</td>
<td>https://www.disgenet.org/</td>
</tr>
<tr class="even">
<td>Disease Ontology</td>
<td>v2021-10-11</td>
<td>https://disease-ontology.org/</td>
</tr>
<tr class="odd">
<td>Human Phenotype Ontology</td>
<td>v2021-10-10</td>
<td>https://hpo.jax.org/app/</td>
</tr>
<tr class="even">
<td>PhosphoSitePlus</td>
<td>v6.6</td>
<td>https://www.phosphosite.org/</td>
</tr>
<tr class="odd">
<td>UbiBrowser</td>
<td>v2.0</td>
<td>http://ubibrowser.ncpsb.org.cn</td>
</tr>
<tr class="even">
<td>ENCODE</td>
<td>v120</td>
<td>https://www.encodeproject.org/</td>
</tr>
<tr class="odd">
<td>STRINGdb</td>
<td>v11.0</td>
<td>https://string-db.org/</td>
</tr>
</tbody>
</table>


Installation
------------

**Hardware requirements**
- PC with 8G RAM or above is recommended
- Passed test on Windows 10, Ubuntu 18.04.6, macOS 14.6

**Software requirements**
- Depends: R (>= 3.6.1)
- Imports: limma, gprofiler2, ComplexHeatmap, igraph, NormExpression


**Run in R session:**

``` r
install.packages("devtools")
install.packages("remotes")
devtools::install_github("TMCjp/ToxDAR")
```

If there's any problem during installation, please refer to <a href="#FAQ">FAQ</a>.

How to use
----------


#### Try the simplest run with default data
> To assist users in reproducing the results, we have integrated the sample data into the R package. This includes public transcriptome data from the L02 cell line following exposure to Triphenyl Phosphate (TPP). Users can load the data with a single command.

``` r
library(ToxDAR)
load(file = paste0(system.file(package = "ToxDAR"), '/extdata/genexp.rdata'))
```

#### Data Preparation
> We utilize the ToxDAR software package to normalize the transcriptome data post-TPP exposure using ten different methods, subsequently generating a report on the effects of normalization. These normalization methods are evaluated based on two metrics: AUCVC and mSCC, to select the appropriate normalization method.

``` r
### 多种预处理方法
HG7.norfac <- gatherFactors(genexp, methods = "HG7")[, 1]
TC.norfac <- gatherFactors(genexp, method = "TC")[, 1]
ERCC.norfac <- gatherFactors(genexp, method = "ERCC")[, 1]
TN.norfac <- gatherFactors(genexp, method = "TN")[, 1]
CR.norfac <- gatherFactors(genexp, method = "CR")[, 1]
UQ.norfac <- gatherFactors(genexp, method = "UQ")[, 1]
NR.norfac <- gatherFactors(genexp, method = "NR")[, 1]
TMM.norfac <- gatherFactors(genexp, method = "TMM")[, 1]
DESeq.norfac <- gatherFactors(genexp, method = "DESeq")[, 1]
TU.norfac <- gatherFactors(genexp, method = "TU")[, 1]

### 根据标准化计算的factor，将表达谱进行预处理
norm.matrix <- getNormMatrix(data = genexp, norm.factors = NR.norfac)

### plot auc  几种分析方法的auc比较
cv_uniform1 <- gatherCVs(data = genexp, nonzeroRatio = 0.5,
                        NR = NR.norfac,
                        DESeq = DESeq.norfac,
                        TMM= TMM.norfac,
                        HG7 = HG7.norfac,
                        TC = TC.norfac,
                        ERCC = ERCC.norfac,
                        TN = TN.norfac,
                        CR = CR.norfac,
                        UQ = UQ.norfac,
                        TU= TU.norfac)

cv_uniform2 <- gatherCVs4Matrices(None= genexp, raw_matrix =genexp, nonzeroRatio=1);
cv_uniform <- rbind(cv_uniform1, cv_uniform2);

cv_uniform <- as.data.frame(cv_uniform)
cv_uniform$Cutoff <- as.numeric(cv_uniform$Cutoff)
cv_uniform$Counts <- as.numeric(cv_uniform$Counts)

### plot
library(ggplot2)
png(file = "cv3.png", res=300, width=(1200*4.17), height=(960*4.17))
# plotCVs(cv_uniform, methods=c("NR", "DESeq", "TMM", "TU"),
#         legend.position=c(.85, .48))
plotCVs(cv_uniform, methods=c("None", "NR", "DESeq", "TMM", "HG7", "TC", "ERCC", "TN", "CR", "UQ", "TU"),
        legend.position=c(1.12, .48))
dev.off()
```
<img src = "vignettes/cv4.png" width = "80%" align = "center" />
<br></br>


``` r
### spearman correlation
spearman_corr1 <- gatherCors(data= genexp, cor_method="spearman",
                            NR = NR.norfac,
                            DESeq = DESeq.norfac,
                            TMM= TMM.norfac,
                            HG7 = HG7.norfac,
                            TC = TC.norfac,
                            ERCC = ERCC.norfac,
                            TN = TN.norfac,
                            CR = CR.norfac,
                            UQ = UQ.norfac,
                            TU = TU.norfac,
                            pre_ratio=1, lower_trim=0.2, upper_trim=0.6, rounds=10000)

spearman_corr2 <- gatherCors4Matrices(None= genexp, raw_matrix=genexp, cor_method="spearman",
                                      pre_ratio=1, lower_trim=0.2, upper_trim=0.6, rounds=10000)
spearman_corr <- rbind(spearman_corr1, spearman_corr2);

### string to numeric
spearman_corr <- as.data.frame(spearman_corr)
spearman_corr$Value <- as.numeric(spearman_corr$Value)

### mSCC
cor.medians <- getCorMedians(spearman_corr)

### plot mSCC
png(file = "cor3.png", res=300, width=(1200*4.17), height=(960*4.17))
plotCors(spearman_corr, methods=c("None", "NR", "DESeq", "TMM", "HG7", "TC", "ERCC", "TN", "CR", "UQ", "TU"),
         legend.position=c(1.12, .48));
dev.off()
```
<img src = "vignettes/cor4.png" width = "80%" align = "center" />
<br></br>



#### Quality Control
> We generated PCA plots and violin plots of gene expression distribution for these datasets, allowing us to visually observe the variability both between different experimental groups and within the same group across samples.

``` r
### PCA分析，参数 stype是样本分类信息，outpng是逻辑值表示是否输出PCA的结果图
gi <- c(rep("control",5), rep("55mg/kg",3), rep("110mg/kg",3), 
        rep("220mg/kg",3), rep("441mg/kg",3), rep("881mg/kg",3))
pcaresult <- toxPCA(genexp, stype = gi, outpng = T)
```
<img src = "vignettes/toxPCA.png" width = "80%" align = "center" />
<br></br>

``` r
### 箱线图
png(file = "toxVio.png", res=300, width=(1200*3), height=(960*3))
plotVio(genexp, stype = gi)
dev.off()
```
<img src = "vignettes/toxVio.png" width = "80%" align = "center" />
<br></br>



#### Data Analysis
> We utilized the differential analysis module of the ToxDAR software package to identify differentially expressed genes associated with TPP exposure at various concentrations and time points. Additionally, the integrated knowledge base within ToxDAR provided clear associations between differentially expressed genes and TPP exposure.

``` r
genexp <- expdata[, c(2:6, 19,20,21)]

### 差异分析，参数org表示物种信息，toxid表示毒物ID，火山图可以直接标注已知毒物相关基因。差异方向默认是后者比前者
gi <- c(rep("control",5), rep("881mg/kg",3))
diffmat <- diffexp(genexp, groupinfo = gi, ntop = 10000)

annomat <- diffanno(diffmat, org = "rnorvegicus", toxid = 'C005445')
# 交互结果表格
IntTab(diffres = annomat, outhtml = T)

# 火山图，标注了已知与所分析毒物相关的基因
library(dplyr)
png(file = "volanno.png", res=300, width=(960*3), height=(1000*3))
volanno(diffres = annomat, nodesize = 1.5, annosize = 5, themsize = 24, legend.position = c(0.82, 0.15), legend.size = 12)
dev.off()
```
<img src = "vignettes/volanno.png" width = "80%" align = "center" />
<br></br>

> To gain a comprehensive understanding of the functions of these genes, we proceeded to perform enrichment analysis on the differentially expressed genes following TPP exposure using multiple annotation datasets including GO, DO, KEGG, and others collected within the software package.
``` r
library(clusterProfiler)
# 加载数据
gsy <- annomat$Symbol[1:30]
gis <- gconvert(query = gsy, organism = "hsapiens",
                target="ENTREZGENE_ACC", mthreshold = Inf, filter_na = TRUE)
gis <- gis$target

# 功能分析
# gofact分析，slim可选 go/kegg/do
gomat <- gofact(genes = gis, slim = 'go')

# function bar
png(file = "funcbar.png", res=300, width=(1200*3), height=(960*3))
p <- plotfact(gomat, ntop = 10, fillby = "Pval", tn = 30, collimit = c(0, 0.05), xtextsize = 20, stripsize = 18)
p + theme(axis.text = element_text(size=20), axis.title = element_text(size=30), 
          legend.text = element_text(size=20), legend.title = element_text(size=30), legend.key.size = unit(20, "pt"),
          axis.ticks = element_blank())
dev.off()
```
<img src = "vignettes/funcbar.png" width = "80%" align = "center" />
<br></br>


``` r
# function circos
# select 功能条目某一大类
tys <- unlist(lapply(gomat$hie, function(x) strsplit(x, '.', fixed = T)[[1]][1]))
bpmat <- gomat[which(tys == 'BP'), ]

png(file = "funcchord.png", res=300, width=(2200*3), height=(2000*3))
plotexpfunc(bpmat, ntop = 5, diffresult = annomat, plottype = 'chord', chord.space = 0.05, gene.space = 0.25, gene.size = 9, process.title = 30, process.label = 28, termcol = 3)
dev.off()

png(file = "funccluster.png", res=300, width=(1200*3), height=(1200*3))
plotexpfunc(bpmat, ntop = 5, diffresult = annomat, plottype = 'cluster')
dev.off()
```
<img src = "vignettes/funcchord.png" width = "80%" align = "center" />
<br></br>
<img src = "vignettes/funccluster.png" width = "80%" align = "center" />
<br></br>


> The ssGSEA method quantifies changes in gene expression profiles as a result of toxic exposure into the activation or inhibition states of biological pathways, revealing the mechanisms by which the organism responds to exposure at a molecular level. The changes are presented visually through enrichment plots and heatmaps for intuitive representation.
``` r
### gsea analysis ES
nn <- which(is.na(annomat$Symbol))
gseamat <- annomat[-nn, ]

geneList <- gseamat$logFC
genename <- gseamat$Symbol

library(org.Hs.eg.db)
gid <- mapIds(org.Hs.eg.db, keys=genename, column="ENTREZID", keytype="SYMBOL", multiVals = "first")
names(geneList) <- gid

geneList <- geneList[order(geneList, decreasing = T)]

# gsea分析
gsearesult <- gseafc(geneList, term = "KEGG", minGSSize = 5, maxGSSize = 500)


# plot gsea result
gsid <- gsearesult$gs_id[which(gsearesult$gs_description == 'Apoptosis')]
png(file = "gsea_up2.png", res=300, width=(1200*3), height=(960*3))
plotgsea(gsid, geneList, term = "KEGG")
dev.off()
```
<img src = "vignettes/gsea_up2.png" width = "80%" align = "center" />
<br></br>


``` r
# plot ES heatmap
esmicro <- as.matrix(gsearesult[, 'ES'])
rownames(esmicro) <- gsearesult$gs_exact_source

ESheatmap(esmicro, Normalization = F)
```
<img src = "vignettes/KeggSig.png" width = "80%" align = "center" />
<br></br>



#### Data Interpretation
> We employed the ToxDAR network analysis and annotation module to map the associations between the toxin and key molecules, biological pathways, and phenotypes. Within this network, TPP is linked with biological pathways such as apoptosis and inflammatory responses, indicating that triphenyl phosphate (TPP) may act as a potential endocrine disruptor. It exerts its effects by activating the thyroid hormone receptor β (THRB) molecule, thereby influencing signaling pathways related to apoptosis and inflammatory responses, ultimately adversely affecting thyroid function. 

``` r
netmat <- netmoa(toxid = 'C005445', ngid = gis, toppathway = 30, ppiscore = 500)
# 绘制AOP结果
# png(file = "netmoa.png", res=600, width=(1600*5), height=(1200*5))
# netplot(netmoa = netmat)
# dev.off()

### html network
nethtml <- netplothtml(netmoa = netmat, outtype = 'pdf')
savehtml(network = nethtml, htmlfile = 'netmoa.html')
```
<img src = "vignettes/netmoa.png" width = "80%" align = "center" />
<br></br>


FAQ
---

If `devtools::install_github()` raise `Installation failed: Problem with the SSL CA cert (path? access rights?)` error, try:

``` r
install.packages(c("curl", "httr"))
```

During installation there may be some configuration error (lack of libraries):

``` pre
------------------------- ANTICONF ERROR ---------------------------
Configuration failed because libcurl was not found. Try installing:
 * deb: libcurl4-openssl-dev (Debian, Ubuntu, etc)
 * rpm: libcurl-devel (Fedora, CentOS, RHEL)
 * csw: libcurl_dev (Solaris)
If libcurl is already installed, check that 'pkg-config' is in your
PATH and PKG_CONFIG_PATH contains a libcurl.pc file. If pkg-config
is unavailable you can set INCLUDE_DIR and LIB_DIR manually via:
R CMD INSTALL --configure-vars='INCLUDE_DIR=... LIB_DIR=...'
--------------------------------------------------------------------
```

Just follow the instruction to satisfy the dependencies. For instance, you can run `sudo apt-get install libcurl4-openssl-dev` in *Ubuntu* to fix the problem above.

> ToxDAR is developed using the R programming language and integrates numerous commonly used packages provided by Bioconductor, so if you get `'BiocInstaller' must be installed to install Bioconductor packages.`, please choose `1 (Yes)`. Since then you may see `Installation failed: cannot open the connection to 'https://bioconductor.org/biocLite.R'`, run `source('http://bioconductor.org/biocLite.R')`, finally try the installation commands above again.

> Please wait for minutes then **try again** if solving some dependencies from *GitHub* fails with `Connection timed out after 100001 milliseconds`.

License
-------

This package is free and open source software, licensed under [GPL v3.0](LICENSE).
