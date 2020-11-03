setwd("/home/zsuzsi/Dropbox/egyeb/R/RNAseq/")

library(DESeq2)
library("pheatmap")
library("RColorBrewer")
library(tidyverse)
library("AnnotationDbi")
library(org.Mm.eg.db)

annot<-read.csv(file= "genes_annotation.csv",sep = ",",header = T)
raw_counts<-read.csv(file= "raw_counts_emptySE.csv",row.names=1,sep = ",",header = T)
raw_counts_empty_SE_stim<-raw_counts[ c(1,3,5,7,9,11) ]
metadata<-read.csv(file= "sampleinfo_emptySE.csv",row.names=1, header = T)

#metadata<-metadata[!(metadata$Status=="unstim"),]
rownames(metadata)<-metadata$FileName
metadata<-dplyr::select(metadata,CellType,Status)

#all(colnames(raw_counts) == rownames(metadata))
#match(colnames(raw_counts), rownames(metadata))

dds<-DESeqDataSetFromMatrix(countData =  raw_counts,
                       colData =  metadata,
                       design = ~ CellType + Status)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

#unsupervised clustering analyses
vsd<- vst(dds, blind=TRUE)
vsd_mat <- assay(vsd)
vsd_cor <- cor(vsd_mat) 
pheatmap(vsd_cor, annotation = dplyr::select(metadata, Status))

#PCA
pcaData <- plotPCA(vsd, intgroup=c ("CellType", "Status"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=CellType, shape=Status)) +  theme_linedraw() +
  geom_point(size=4) +
  scale_color_manual(values = c("black","blue")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  expand_limits(y = c(-40, 40)) +
  ggtitle("Principal component analysis") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed()


#differencial expression
dds<-DESeq(dds)
plotDispEsts(dds)

#contrasts
stim_emptySE_res <- results(dds,
                            contrast = c("CellType","SE","empty"),
                              alpha = 0.05)


#LFC shrinkage
library("apeglm")
stim_emptySE_res_apeglm <- lfcShrink(dds, coef="CellType_SE_vs_empty", type="apeglm")

#Vulcano plots
stim_emptySE_res_df<-data.frame(stim_emptySE_res_apeglm)
ggplot(stim_emptySE_res_df) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj))) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))



stim_emptySE_res_df$ensgene<-gsub("\\..*","",rownames(stim_emptySE_res_df))
stim_emptySE_res_df$symbol <- mapIds(org.Mm.eg.db,
                                    keys=stim_emptySE_res_df$ensgene,
                                    column="SYMBOL",
                                    keytype="ENSEMBL",
                                    multiVals="first")
stim_emptySE_res_df$entrez <- mapIds(org.Mm.eg.db,
                     keys=stim_emptySE_res_df$ensgene,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

genes_emptySE<- stim_emptySE_res_df %>% dplyr::filter(!is.na(stim_emptySE_res_df$symbol)) %>%
  arrange(desc(abs(log2FoldChange)))
stim_emptySE_res_sig<-subset(genes_emptySE,padj<0.05)

#visualization
heat_colors <- brewer.pal(11,"RdYlBu")
normcounts<-as.data.frame(normalized_counts)
normcounts$ensgene<-gsub("\\..*","",rownames(normcounts))
my.heatmap.data<-stim_emptySE_res_sig %>% left_join(normcounts)


rownames(my.heatmap.data) <- my.heatmap.data$symbol
my.heatmap.data <- my.heatmap.data[,9:20]
my.heatmap.data <- as.matrix(my.heatmap.data)
mode(my.heatmap.data)<-"numeric"

my_colour = list(
  CellType = c(empty = "black", SE = "blue"))
pheatmap(my.heatmap.data,
         color = heat_colors,
         cluster_rows = T,
         show_rownames = T,
         annotation = dplyr::select(metadata, CellType),
         annotation_colors = my_colour,
         scale = "row")

