rna_counts <- readRDS("rna_counts_1095_patients.rds")
clinical <- readRDS("TCGA_BRCA_clinical.rds")

clinica.rna <- clinical[ clinical$submitter_id %in% colnames(rna_counts), ]

rownames(clinica.rna) <- clinica.rna$submitter_id

rna_counts_o <- rna_counts[, rownames(clinica.rna)]

table(rownames(clinica.rna) == colnames(rna_counts_o))

library(DESeq2)
library(edgeR)

dds = DESeqDataSetFromMatrix(countData = rna_counts_o, colData = clinica.rna, design = ~ 1)

keep = rowSums(cpm(dds) > 5) >= 0.05*ncol(rna_counts_o) ## at least 5% samples with count of 5 or higher

dds = dds[keep,]

nrow(dds) ## 47609


vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)

colData(vsd)

vst.mat <- assay(vsd)

gene.var <- apply( vst.mat, 1, var)

gene.var

summary(gene.var)

quantile(
  gene.var,
  probs = c(0.25, 0.5, 0.75, 0.90, 0.95, 0.99)
)

genes.var.75.quantile <- gene.var[gene.var > 0.6923132] ## using 75th quantile

genes.keep <- names(genes.var.75.quantile)

vsd.75.quantile <- vst.mat[ genes.keep , ]
  

write.table(vsd.75.quantile,file="TCGA_BRCA_VSD_transformed_expression_matrix_above_75th_quantile_variance.txt",col.names = T,row.names = T, sep="\t", quote = FALSE)



expr_mut_count <-  read.table(file="TCGA_BRCA_Expresison_Mutation_Count_Matrix.txt", header = T,sep="\t", stringsAsFactors = FALSE)

expr_mut_status <-  read.table(file="TCGA_BRCA_Expresison_Mutation_Status_Matrix.txt", header = T,sep="\t", stringsAsFactors = FALSE)



mut_count <- expr_mut_count[, grep("count", colnames(expr_mut_count))]

mut_status <- expr_mut_status[, grep("status", colnames(expr_mut_status))]

mut_count <- t(mut_count)

mut_status <- t(mut_status)

overlapping_patients <- intersect(colnames(mut_count), colnames(vsd.75.quantile))

mut_count_overlapping <- mut_count[, overlapping_patients]

mut_status_overlapping <- mut_status[, overlapping_patients]

vsd.75.quantile.overlapping <- vsd.75.quantile[, overlapping_patients]

clinical.overlapping <- clinical[ clinical$submitter_id %in% overlapping_patients, ]

clinical.overlapping.export <- clinical.overlapping

list.cols <- sapply(clinical.overlapping.export, is.list)

clinical.overlapping.export[list.cols] <- lapply(
  clinical.overlapping.export[list.cols],
  function(x) sapply(x, paste, collapse = ";")
)

write.csv(vsd.75.quantile,file="TCGA_BRCA_VSD_transformed_expression_matrix_above_75th_quantile_variance_overlapping_patients.csv",col.names = T,row.names = T, quote = FALSE)

write.csv(mut_count_overlapping,file="TCGA_BRCA_mutation_count_overlapping_patients.csv",col.names = T,row.names = T, quote = FALSE)

write.csv(mut_status_overlapping,file="TCGA_BRCA_mutation_status_overlapping_patients.csv",col.names = T,row.names = T, quote = FALSE)

write.table(clinical.overlapping.export,file="TCGA_BRCA_clinical_information_overlapping_patients.csv",col.names = T,row.names = T, quote = FALSE)


  
