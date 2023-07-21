library (DESeq2)
library (ggplot2)
library(ggrepel)
library (dplyr)


## normal STAR results from the RNA-Seq IIT pipeline
a <- read.delim ("star_gene_raw_counts.txt")
annot <- a[ ,grep ("TMZ", colnames (a), invert=TRUE)]

counts <- a[ ,grep ("Gene|TMZ", colnames (a))]
row.names (counts) <- counts$Geneid
counts <- counts[ ,-1]
head (counts)

samples <- data.frame (matrix (nrow=dim (counts)[2], ncol=3))
colnames (samples) <- c("sample", "cell", "condition")
samples$cell <- gsub ("_.*", "", gsub ("TMZ_", "", colnames (counts)) )
samples$condition <- gsub ("0-9", "", gsub (".*_", "", colnames (counts)))
samples$sample <- row.names (samples)
row.names (samples) <- colnames (counts)


## See Github RNA-Seq_mouse/gene_annotation.R
#system ("cp /projects/ncrrbt_share_la/dev_pipe/gencode.vM32.annotation.txt .")

anno <- read.delim ("gencode.vM32.annotation.txt")

anno <- anno[ ,grep ("transcript_id", colnames (anno), invert=TRUE)]
anno <- unique (anno)



## normal STAR results from the RNA-Seq IIT pipeline
a <- read.delim ("subread.counts.txt", skip=1)
a <- a[ ,grep ("Gene|bam", colnames (a))]

a <- merge (a, anno, by.x="Geneid", by.y="gene_id", all.x=TRUE) 

a <- a[grep ("miRNA|Mt_tRNA|Mt_rRNA|rRNA|snRNA|snoRNA|scRNA|sRNA|misc_RNA|scaRNA|ribozyme|IG_|TR_", a$gene_type, invert=TRUE), ]
colnames (a) <- gsub ("Aligned.out.bam", "", colnames (a))
colnames (a) <- gsub ("star.", "", colnames (a))


counts <- annot <- a

annot <- annot[ ,c("Geneid", "gene_name", "gene_type", "mgi_id", "external_gene_name", "description")]

row.names (counts) <- counts$Geneid
counts <- counts[ ,grep ("IIT", colnames (counts))]

colnames (counts) <- gsub ("IIT_SHA_LINE1_KO_", "", colnames (counts))
colnames (counts) <- gsub ("_S.*", "", colnames (counts))

samples <- data.frame (matrix (nrow=dim (counts)[2], ncol=3))
colnames (samples) <- c("sample", "condition", "area")
samples$area <- "OB"
samples$condition <- factor (c(rep ("KO",4), rep ("WT", 4)))
row.names (samples) <- colnames (counts)
samples$sample <- row.names (samples)
