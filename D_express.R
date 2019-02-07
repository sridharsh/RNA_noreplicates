library(DESeq2)
D_express <- function(in_file, condition){
  countdata <- as.matrix(in_file)
  (coldata <- data.frame(row.names=colnames(countdata), condition))
  dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
  dds <- DESeq(dds)
  res<-results(dds)
  res_filt <- res[ which(res$padj < 0.0000000001), ] # filtering based on p-val cut off of 1e-10 and saving it into a variable
  resOrdered <- res_filt[order(-res_filt$log2FoldChange),] # ordering based on fold change values
  write.csv(as.data.frame(resOrdered),file="Sensitive_vs_Resistant.csv")
  dds
}

heat_mapping <- function(dds){
  select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:20]
  nt <- normTransform(dds) # defaults to log2(x+1)
  log2.norm.counts <- assay(nt)[select,]
  df <- as.data.frame(colData(dds)[,c("condition", "sizeFactor")]) 
  df$sizeFactor <- NULL # removing sizeFactor column out of colData(dds)
  library("pheatmap")
  quartz()
  pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=TRUE,
           cluster_cols=FALSE, annotation_col=df, filename = "heatmap_deseq.pdf", fontsize = 6)
}


#df <- read.table(args[1], header=TRUE)
setwd("/out_dir")
df <- read.table("./featureCounts/cleaned_counts.csv", 
                 sep = ",", header = TRUE, row.names="GeneID")
condition <- factor(c("Sensitive", "Resistant", "Resistant", "Resistant"))
heat_mapping(D_express(df, condition))
