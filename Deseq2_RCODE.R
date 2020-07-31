---
title: "R Notebook"
output: html_notebook
---

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you execute code within the notebook, the results appear beneath the code. 

Try executing this chunk by clicking the *Run* button within the chunk or by placing your cursor inside it and pressing *Cmd+Shift+Enter*. 
```{r}
#-- provide directory with htseq-count output files
dir = ("/Input_HTCOUNTS")
getwd()
setwd("/Input_HTCOUNTS")
#-- define the pattern of files to be analysed, the file shold end as _htcount.txt
pattern="*_htcount.txt"
currentjob <- "GSE116250"
#-- provide metadata file containg file names and their the replicate information as below
#metadata_file <- "metadata.txt"

outfile <- "_count_matrix.tab"

#' Title: deseq_from_htseqcount
#' deseq_from_htseqcount(dir, pattern)
#' @param dir 
#' @param pattern 
#' @param metadata_file 
#'treatment_set1	treated
# treatment_set2	treated
# control_set1	untreated
# control_set2	untreated
#' @param outfile 
#' @author ankushs  
#' University of Oslo
#' @return
#' @export
#'
#' @examples
deseq_from_htseqcount <- function(dir=dir, pattern, metadata_file,outfile)
  #install.packages("GGally")
  #install.packages("openxlsx")
  #if (!requireNamespace('BiocManager', quietly = TRUE))
    #install.packages('BiocManager')
#BiocManager::install('EnhancedVolcano')

  #--- Load package
  library(tidyverse)
  library(dplyr)
  library(data.table)
  library(DESeq2)
  library(purrr)
  library(GGally)
  library(openxlsx)
library(EnhancedVolcano)
library("AnnotationDbi")
library("org.Hs.eg.db")
library("dplyr")
#--- Load DESeq2 files
dd <- list.files(path = dir,
                   pattern=pattern, 
                   full.names = F)
  print(dd)
#--- get count matrix          
df_count_matrix <- data_frame(file_name = dd) %>% 
    mutate(file_cont = map(file_name,fread,data.table = F))  %>%
    unnest() %>% 
    mutate(file_name = gsub(pattern="_htcount.txt",replacement="",file_name))  %>% 
    spread(key = file_name , value = V3) %>%
    dplyr::slice(0:nrow(.)) #removing unwanted alignment summary present at top #change 0 to number of lines to be removed

 
#--- write count matrix in a file
write_delim(df_count_matrix, paste(outfile,"_count_matrix.tab", sep=""),col_names = TRUE,delim="\t")
setwd("ISCHEMIA_RNA_SEQ_116520/FOR_PAPER")
write.table(df_count_matrix, file="ISCHEMIA_RNA_SEQ_116520/FOR_PAPER/COUNT_MATRIX.TAB", sep="\t")

``` 


```{r}
Readcount <- read.table("COUNT_MATRIX.TAB", header=TRUE, sep = "\t", row.names = 1)
metaData <- read.csv('../samplename1.txt', header = TRUE, sep = "\t")
metaData
  
names<- make.unique(as.character(Readcount$ENSEMBL))
  #CHecking duplicated Rownames
  #unique(Readcount[,1])
  #DUPID = duplicated(Readcount[,1])
  #write.table(DUPID,"xew.txt")
  Readcount
  
  row.names(Readcount)
  colnames(Readcount)
  countdata1 <- Readcount
  countdata1
  #plot(Readcount)
  #png("readcount.png",1000, 1500, pointsize=14)
  #dev.off()
  
  #Detecting "NA" in data 
  #c= sum(is.na(countdata1))
  #summarise(Readcount)
  
  #FINDING DUPLICATE
  #x <- countdata1[ ,(1)]
  #duplicated(x)
  #y = x[duplicated(x)] 
  #write.table(y, file="duplicatelist.txt")
  #y

```


```{r Selecting  <columns, include=FALSE}
#for selecting required columns
#countdata1
#countdata <- countdata1[ ,3:ncol(countdata)]

countdata <- countdata1[ ,c(2:28)]
  # Convert to matrix
  countdata <- as.matrix(countdata)
  library("DESeq2")

countdata_symbol <- countdata
rownames(countdata_symbol) <- Readcount$GeneSymbol
```

```{r}
# Assign condition (first four are controls, second four contain the expansion)
x <-(condition <- factor(c(rep("NF", 14), rep("Ischemia", 13))))
(condition <- factor(c(rep("NF", 14), rep("Ischemia", 13))))
```
# Analysis with DESeq2 ----------------------------------------------------
```{r}
library(DESeq2)

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata_symbol), condition))
View(coldata)

dds <- DESeqDataSetFromMatrix(countData=countdata_symbol, colData=coldata, design=~condition)
dds

# Run the DESeq pipeline
dds <- DESeq(dds)

# Plot dispersions
png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))

# Colors for plots below
## Ugly:
## (mycols <- 1:length(unique(condition)))
## Use RColorBrewer, better
library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

# Principal components analysis
## Could do with built-in DESeq2 function:
## DESeq2::plotPCA(rld, intgroup="condition")
## I like mine better:
rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
  #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
}
png("qc-pca.png", 1000, 1000, pointsize=20)
rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-75, 35))
dev.off()


# Get differential expression results
res <- results(dds)
table(res$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
## Write results
write.csv(resdata, file="diffexpr-results.csv")

## Examine plot of p-values
hist(res$pvalue, breaks=50, col="grey")

## Examine independent filtering
#attr(res, "filterThreshold")
#plot(attr(res,"filterNumRej"), type="b", xlab="quantiles of baseMean", ylab="number of rejections")

## MA plot
## Could do with built-in DESeq2 function:
## DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)
## I like mine better:
maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
png("diffexpr-maplot.png", 1500, 1000, pointsize=20)
maplot(resdata, main="MA Plot")
dev.off()







res
  #columns(org.Hs.eg.db)
  #ens.str <- substr(rownames(res), 1, 15)
  
  #res$symbol <- mapIds(org.Hs.eg.db,
   #                  keys=ens.str,
    #                 column="SYMBOL",
     #                keytype="ENSEMBL",
      #               multiVals="first")
 #res$symbol
  #resOrdered <- res[order(res$pvalue),]
  #head(res)
  
  #rownames(resOrdered) <- symbols
  ############
#resOrderedDF <- as.data.frame(resOrdered)[1:58395, ]
  ########
#res
```

```{r}
## Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="wheat", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="wheat", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("blue","wheat","red"))
}
png("diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
dev.off()
  ########################################
```


```{r}
 
######################################
png("ENHANCED-volcanoplot.png", 1200, 1000, pointsize=20)
EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'pvalue',
    xlim = c(-6,6),
    xlab = bquote(~Log[2]~ 'fold change'),
    pCutoff = 10e-14,
    FCcutoff = 2.0,
    pointSize = 4.0,
    labSize = 4.0,
    colAlpha = 1,
    legendPosition = 'right',
    legendLabSize = 12,
    legendIconSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 0.2,
    colConnectors = 'grey30')
dev.off()
```



```{r}




 outputPrefix <- "NFVSIschemia"
library("pheatmap")
  mat = assay(rld)[ head(order(resOrdered$padj),60), ] # select the top 30 genes with the lowest padj
  mat
  mat = mat - rowMeans(mat) # Subtract the row means from each value
  # Optional, but to make the plot nicer:
  df = as.data.frame(colData(rld)[,c("condition")]) # Create a dataframe with a column of the conditions
  colnames(df) = "condition" # Rename the column header
  rownames(df) = colnames(mat) # add rownames
  rownames(df)
  # and plot the actual heatmap
  
png(paste("HEATMAP_TOPGENES.png", sep=""),width = 1200, height = 1200, units = "px", pointsize = 12)
  pheatmap(mat, annotation_col=df)
  dev.off()
  dev.set(dev.next())
  #correlation matrix 
  
  
  resdata1 <- merge(as.data.frame(resdata), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  
  names(resdata1)[1] <- "Gene"
  #resdata = data.frame(resdata,significance) 
  print(head(resdata1))
  
  expression_matrix <- resdata1[,c(1,8:ncol(resdata))]
  print(head(expression_matrix))
  
  #--- get correlation matrix
  gg <- ggpairs(expression_matrix, columns = 2:ncol(expression_matrix), upper = list(continuous = wrap("cor", size = 5, color="black"))) +
    theme_bw() +
    theme(legend.text = element_text(size = 25),
          legend.title = element_text(size = 20),         
          axis.title.x = element_text(size = 15),        
          axis.title.y = element_text(size = 15),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10)) 


 
  #--- write deseq output
  up_deg <- subset(resdata, resdata$log2FoldChange> 0.6 & pvalue <= 0.05)
  down_deg <- subset(resdata, resdata$log2FoldChange< -0.6 & pvalue <= 0.05)
  write.table(resdata, file="Ischemia Vs Nonfailing.tab")
  openxlsx::write.xlsx(resdata, file=paste("ALL_DE_Output.tab", sep="\t"),row.names = F, append=F)
  openxlsx::write.xlsx(up_deg,file=paste("up_deg.xlsx", sep=""),row.names = F, append=T)
  openxlsx::write.xlsx(down_deg,file=paste("down_deg.xlsx", sep=""),row.names = F, append=T)
  
  
```


```{r}
  




columns(org.Hs.eg.db)
  ens.str <- substr(rownames(res), 1, 15)
  
  res$symbol <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
  res$symbol
  resOrdered <- res[order(res$pvalue),]
  head(resOrdered)
  
resOrderedDF <- as.data.frame(resOrdered)[1:20000, ]
 write.csv(resOrderedDF, file = "results.csv")
```
Add a new chunk by clicking the *Insert Chunk* button on the toolbar or by pressing *Cmd+Option+I*.


When you save the notebook, an HTML file containing the code and output will be saved alongside it (click the *Preview* button or press *Cmd+Shift+K* to preview the HTML file). 

The preview shows you a rendered HTML copy of the contents of the editor. Consequently, unlike *Knit*, *Preview* does not run any R code chunks. Instead, the output of the chunk when it was last run in the editor is displayed.

