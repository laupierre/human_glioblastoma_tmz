library (DESeq2)
library (ggplot2)
library(ggrepel)
library (dplyr)


## normal STAR results from the human RNA-Seq IIT pipeline
a <- read.delim ("star_gene_raw_counts.txt")
annot <- a[ ,grep ("TMZ", colnames (a), invert=TRUE)]
head (annot)

counts <- a[ ,grep ("Gene|TMZ", colnames (a))]
row.names (counts) <- counts$Geneid
counts <- counts[ ,-1]

samples <- data.frame (matrix (nrow=dim (counts)[2], ncol=3))
colnames (samples) <- c("sample", "cell", "condition")
samples$cell <- gsub ("_.*", "", gsub ("TMZ_", "", colnames (counts)) )
samples$condition <- gsub ("[0-9]", "", gsub (".*_", "", colnames (counts)))
samples$sample <- colnames (counts)
row.names (samples) <- colnames (counts)


# U87 cells
cell_type <- "U87"
# cell_type <- "GBM6"

samples.s <- samples[samples$cell == cell_type, ]
samples.s

counts.s <- counts[ ,colnames (counts) %in% samples.s$sample]
stopifnot (colnames (counts.s) == samples.s$sample)


## DESeq2 
dds <- DESeqDataSetFromMatrix(countData = round (counts.s), colData = samples.s, design = ~ condition)
                                 
keep <- rowSums(counts(dds) >= 10) >= dim (counts.s)[2]/2
dds <- dds[keep,]
dds

dds$condition <- relevel(dds$condition, "C")

dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, contrast=list("condition_T_vs_C"))

res <- merge (data.frame (res), round (counts (dds, normalized=TRUE)), by="row.names")
res <- merge (res, annot, by.x="Row.names", by.y="Geneid", all.x= TRUE)
colnames (res)[1] <- "Geneid"
res <- res[order (res$padj), ]

boxplot (res$log2FoldChange)
abline (h=0)

write.table (res, "TMZ_U87_differential_expression.txt", row.names=F, quote=F, sep="\t")
# write.table (res, "TMZ_GBM6_differential_expression.txt", row.names=F, quote=F, sep="\t")


# PCA plot
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

p1 <- ggplot(pcaData, aes(PC1, PC2, color=condition, shape=condition)) +
  		geom_point(size=3) +
  		xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  		ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
		coord_fixed () + geom_label_repel (aes(label = name))
p1
ggsave ("PCA plot.pdf", p1)




## PCA of all samples

dds <- DESeqDataSetFromMatrix(countData = round (counts), colData = samples, design = ~ condition)
                                 
keep <- rowSums(counts(dds) >= 10) >= dim (counts)[2]/2
dds <- dds[keep,]
dds

vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
pcaData$cell <- gsub ("_.*", "", gsub ("TMZ_", "", pcaData$name))

percentVar <- round(100 * attr(pcaData, "percentVar"))

p1 <- ggplot(pcaData, aes(PC1, PC2, color=cell, shape=condition)) +
  		geom_point(size=3) +
  		xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  		ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
		coord_fixed () # + geom_label_repel (aes(label = cell), max.overlaps = Inf)
p1
ggsave ("PCA all cells plot.pdf", p1)




