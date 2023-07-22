library (readxl)

raw <- data.frame (read_excel("TMZ_rnaseq_raw_counts.xlsx"))
raw <- raw[ ,grep ("TMZ", colnames (raw), invert=TRUE)]

u87 <- data.frame (read_excel("TMZ_U87_differential_expression.xlsx"))
gbm1 <- data.frame (read_excel("TMZ_GBM1_differential_expression.xlsx"))
gbm2 <- data.frame (read_excel("TMZ_GBM2_differential_expression.xlsx"))
gbm5 <- data.frame (read_excel("TMZ_GBM5_differential_expression.xlsx"))
gbm6 <- data.frame (read_excel("TMZ_GBM6_differential_expression.xlsx"))

u87 <- u87[ ,c("Geneid", "baseMean", "log2FoldChange", "padj")]
gbm1 <- gbm1[ ,c("Geneid", "baseMean", "log2FoldChange", "padj")]
gbm2 <- gbm2[ ,c("Geneid", "baseMean", "log2FoldChange", "padj")]
gbm5 <- gbm5[ ,c("Geneid", "baseMean", "log2FoldChange", "padj")]
gbm6 <- gbm6[ ,c("Geneid", "baseMean", "log2FoldChange", "padj")]

colnames (u87)[2:4] <- paste (colnames (u87)[2:4], "_U87", sep="")
colnames (gbm1)[2:4] <- paste (colnames (gbm1)[2:4], "_GBM1", sep="")
colnames (gbm2)[2:4] <- paste (colnames (gbm2)[2:4], "_GBM2", sep="")
colnames (gbm5)[2:4] <- paste (colnames (gbm5)[2:4], "_GBM5", sep="")
colnames (gbm6)[2:4] <- paste (colnames (gbm6)[2:4], "_GBM6", sep="")

mylist <- list (u87, gbm1, gbm2, gbm5, gbm6)

res <- Reduce (function (x,y) merge (x,y, by="Geneid", all.x=TRUE, all.y=TRUE), mylist)
colnames (res) <- gsub ("log2FoldChange", "log2FC", colnames (res))

res <- merge (res, raw, by="Geneid")

library (writexl)

write_xlsx (res, "TMZ_all_cells_differential_expression.xlsx")
               
