setwd("~/Desktop/Related_to_Haider_Lab/TCGA_TNBC")

wes_maf <- readRDS("TCGA_BRCA_WES_protein_altering_MAF.rds")

rna_counts <- readRDS("rna_counts_1095_patients.rds")

clinical <- readRDS("TCGA_BRCA_clinical.rds")

## check type of mutations
table(wes_maf$Variant_Classification)

all_patient_ids <- wes_maf[,c("Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode")]

head(all_patient_ids)


wes_maf$Tumor_Patient_ID <- sapply(strsplit(wes_maf$Tumor_Sample_Barcode, "-"), `[`, 4)

wes_maf$Normal_Patient_ID <- sapply(strsplit(wes_maf$Matched_Norm_Sample_Barcode, "-"), `[`, 4)

wes_maf$Participant_ID <- sapply(strsplit(wes_maf$Tumor_Sample_Barcode, "-"), `[`, 1:3)

wes_maf$Patient_ID <- sapply(
  strsplit(wes_maf$Tumor_Sample_Barcode, "-"),
  function(x) paste(x[1:3], collapse = "-")
)


head(all_patient_ids)

## Primary Solid Tumor
wes_maf_use <- wes_maf[
  sapply(strsplit(wes_maf$Tumor_Sample_Barcode, "-"), `[`, 4) |>
    substr(1, 2) == "01",
]

#BiocManager::install(c(
#  "MutationalPatterns"))


unique_patients <- unique(wes_maf$Patient_ID)

a_patient <- unique_patients[1]


wes_maf_a_patient <- wes_maf_use[ wes_maf_use$patient_id %in% a_patient, ]
  
mutation_counts <- as.data.frame(
  table(
    wes_maf_a_patient$patient_id,
    wes_maf_a_patient$Hugo_Symbol
  )
)

colnames(mutation_counts) <- c("patient_id", "Hugo_Symbol", "n_mutations")

all_mutated_genes <- unique(wes_maf$Hugo_Symbol)


all_mutated_genes <- as.data.frame(all_mutated_genes)


rownames(all_mutated_genes) <- all_mutated_genes$all_mutated_genes

for(patient in unique_patients){
  
  wes_maf_a_patient <- wes_maf_use[ wes_maf_use$patient_id %in% patient, ]
  
  mutation_counts <- as.data.frame(
    table(
      wes_maf_a_patient$patient_id,
      wes_maf_a_patient$Hugo_Symbol
    )
  )
  
  colnames(mutation_counts) <- c("patient_id", "Hugo_Symbol", patient)
  
  rownames(mutation_counts) <- mutation_counts$Hugo_Symbol
  
  all_mutated_genes <- cbind(all_mutated_genes, mutation_counts[rownames(all_mutated_genes), 3])
  
  colnames(all_mutated_genes)[ncol(all_mutated_genes)] <- patient
  
  
}

#all_mutated_genes$all_mutated_genes <- NULL

all_mutated_genes[is.na(all_mutated_genes)] <- 0

all_mutated_genes$all_mutated_genes <- NULL

apply(all_mutated_genes, 2, max)[1:10]
apply(all_mutated_genes, 1, max)

all_mutated_genes <- t(all_mutated_genes)

all_mutated_genes[,c("TP53", "PIK3CA", "GATA3", "BRCA1", "BRCA2")]

mutation_status_matrix <- (all_mutated_genes > 0) * 1


intersect(rownames(mutation_matrix), colnames(rna_counts))


colData_rna <- colnames(rna_counts)

colData_rna <- as.data.frame(colData_rna)

rownames(colData_rna) <- colData_rna$colData_rna

colData_rna$RNA <- 1

library(DESeq2)
library(edgeR)
dds = DESeqDataSetFromMatrix(countData = rna_counts,
                             colData = colData_rna,
                             design = ~ 1)

keep = rowSums(cpm(dds) > 5) >= 0.05*ncol(rna_counts) ## at least 5% samples with count of 1 or higher

dds = dds[keep,]

vsd <- vst(dds, blind = TRUE)

head(assay(vsd), 3)

vst.t.mat <- assay(vsd) ## Expression matrix to use

expression_genes <- rownames(vst.t.mat)

expression_genes <- sub("\\..*", "", expression_genes)  

rownames(vst.t.mat) <-  expression_genes


library(org.Hs.eg.db)  


expression_genes_symbol <- select(
  org.Hs.eg.db,
  keys = expression_genes,
  columns = "SYMBOL",
  keytype = "ENSEMBL"
)

symbol.no.na <- expression_genes_symbol[!is.na(expression_genes_symbol$SYMBOL),]

vst.t.mat.keep <- vst.t.mat[ symbol.no.na$ENSEMBL, ]

stopifnot(rownames(vst.t.mat.keep) == symbol.no.na$ENSEMBL)

rownames(vst.t.mat.keep) <- symbol.no.na$SYMBOL


vst.t.mat.keep <- t(vst.t.mat.keep)

length(intersect(rownames(vst.t.mat.keep), rownames(mutation_status_matrix)))

common.patients <- intersect(
  rownames(vst.t.mat.keep),
  rownames(mutation_status_matrix)
)

vst.common <- vst.t.mat.keep[common.patients, , drop = FALSE]

mutation.status.common <- mutation_status_matrix[ common.patients, , drop = FALSE ]

colnames(mutation.status.common) <- paste(colnames(mutation.status.common), "status",sep="_")


all_mutated_genes_common <- all_mutated_genes[ common.patients, , drop = FALSE]


colnames(all_mutated_genes_common) <- paste(colnames(all_mutated_genes_common), "count",sep="_")


genes_mutation_status <- apply(all_mutated_genes_common, 2, function(x) sum(x > 0) / length(x) * 100)

thresholds <- c(1, 2, 5, 10)

sapply(
  thresholds,
  function(p) sum(genes_mutation_status >= p)
)

genes_mutation_status_min_2_pct <- genes_mutation_status[ genes_mutation_status > 2 ]

## keep genes with at least 25% not zero
#genes_mutated_2_pct <- names(genes_mutation_status)[ genes_mutation_status < 75 ]

genes_with_mutation <-sub("_status$", "", names(genes_mutation_status_min_2_pct))



mutation.status.common.min.2pct <- mutation.status.common[, paste(genes_with_mutation, "_status",sep="")]

all_mutated_genes_common.min.2pct <- all_mutated_genes_common[, paste(genes_with_mutation, "_count", sep="")]


expression.mutation.matrix.status <- cbind(vst.common, mutation.status.common.min.2pct[rownames(vst.common), ])

expression.mutation.matrix.count <- cbind(vst.common, all_mutated_genes_common.min.2pct[rownames(vst.common), ])

stopifnot(dim(expression.mutation.matrix.status) == dim(expression.mutation.matrix.count))

write.table(expression.mutation.matrix.status,file="TCGA_BRCA_Expresison_Mutation_Status_Matrix.txt",col.names = T,row.names = T, sep="\t", quote = FALSE)

write.table(expression.mutation.matrix.count,file="TCGA_BRCA_Expresison_Mutation_Count_Matrix.txt",col.names = T,row.names = T, sep="\t", quote = FALSE)





