library(TCGAbiolinks)
library(dplyr)

# 1) download and prepare expression (example: HTSeq counts)
query.rna <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  experimental.strategy = "RNA-Seq",
  sample.type = "Primary Tumor"
)

rna.results <- getResults(query.rna)

dim(rna.results)
table(rna.results$sample_type)
table(rna.results$experimental_strategy)


GDCdownload(
  query = query.rna,
  directory = "TCGA_BRCA"
)

rna <- GDCprepare(
  query = query.rna,
  directory = "TCGA_BRCA"
)

list.files(
  "TCGA_BRCA",
  recursive = TRUE,
  full.names = TRUE
)[1:20]


length(
  list.files(
    "TCGA_BRCA",
    recursive = TRUE,
    full.names = TRUE
  )
)

object.size(query.rna)


rna.files <- list.files(
  "TCGA_BRCA/TCGA-BRCA",
  recursive = TRUE,
  full.names = TRUE
)

length(rna.files)

rna.files[grep("RNA|Gene_Expression|STAR", rna.files)][1:20]


rna.files <- rna.files[
  grepl(
    "Transcriptome_Profiling|Gene_Expression_Quantification",
    rna.files
  )
]

length(rna.files)
head(rna.files)

list.dirs(
  "TCGA_BRCA/TCGA-BRCA",
  recursive = TRUE
)[1:10]

length(list.dirs(
  "TCGA_BRCA/TCGA-BRCA",
  recursive = TRUE
))

rna.test <- read.delim(
  rna.files[1],
  header = TRUE,
  sep = "\t",
  comment.char = "#",
  check.names = FALSE
)

dim(rna.test)
colnames(rna.test)
head(rna.test)


rna.results <- getResults(query.rna)

dim(rna.results)
colnames(rna.results)

grep(
  "file|case|sample|submitter|barcode",
  colnames(rna.results),
  ignore.case = TRUE,
  value = TRUE
)

rna.map <- rna.results[, c(
  "file_id",
  "file_name",
  "cases.submitter_id",
  "sample.submitter_id",
  "sample_type"
)]

head(rna.map)


length(unique(rna.map$cases.submitter_id))
length(unique(rna.map$sample.submitter_id))

table(duplicated(rna.map$cases.submitter_id))

table(table(rna.map$cases.submitter_id))

head(basename(rna.files))
head(rna.map$file_name)


sum(basename(rna.files) %in% rna.map$file_name)


wes.patients <- unique(maf$patient_id)
length(wes.patients)



rna.patients <- unique(rna.map$cases.submitter_id)

length(rna.patients)
length(intersect(wes.patients, rna.patients))

dup.patients <- names(
  table(rna.map$cases.submitter_id)[
    table(rna.map$cases.submitter_id) > 1
  ]
)

rna.map[
  rna.map$cases.submitter_id %in% dup.patients,
  c(
    "file_id",
    "file_name",
    "cases.submitter_id",
    "sample.submitter_id",
    "sample_type"
  )
]

table(rna.map$sample_type)


id="53184"
sample.counts <- table(rna.map$sample.submitter_id)

sample.counts[sample.counts > 1]


dup.samples <- names(sample.counts[sample.counts > 1])

rna.map[
  rna.map$sample.submitter_id %in% dup.samples,
  c(
    "file_id",
    "file_name",
    "cases.submitter_id",
    "sample.submitter_id"
  )
]

rna.map[
  rna.map$sample.submitter_id %in% dup.samples,
  c(
    "file_id",
    "file_name",
    "cases.submitter_id",
    "sample.submitter_id"
  )
]

## determine whether the duplicate files are identical
for (s in dup.samples) {
  
  rows <- rna.map$sample.submitter_id == s
  
  cat("\n", s, "\n")
  print(
    rna.map[
      rows,
      c(
        "file_id",
        "file_name",
        "cases.submitter_id",
        "sample.submitter_id"
      )
    ]
  )
}

for (s in dup.samples) {
  
  rows <- which(rna.map$sample.submitter_id == s)
  
  f1 <- rna.files[
    basename(rna.files) == rna.map$file_name[rows[1]]
  ]
  
  f2 <- rna.files[
    basename(rna.files) == rna.map$file_name[rows[2]]
  ]
  
  x1 <- read.delim(
    f1,
    header = TRUE,
    sep = "\t",
    comment.char = "#",
    check.names = FALSE
  )
  
  x2 <- read.delim(
    f2,
    header = TRUE,
    sep = "\t",
    comment.char = "#",
    check.names = FALSE
  )
  
  # Remove STAR summary rows
  x1 <- x1[!grepl("^N_", x1$gene_id), ]
  x2 <- x2[!grepl("^N_", x2$gene_id), ]
  
  same <- identical(
    x1$unstranded,
    x2$unstranded
  )
  
  cat(s, ":", same, "\n")
}

## compare their metadata and sequencing depth for duplicate files
for (s in dup.samples) {
  
  rows <- which(rna.map$sample.submitter_id == s)
  
  cat("\n====================\n")
  cat(s, "\n")
  
  for (i in rows) {
    
    f <- rna.files[
      basename(rna.files) == rna.map$file_name[i]
    ]
    
    x <- read.delim(
      f,
      header = TRUE,
      sep = "\t",
      comment.char = "#",
      check.names = FALSE
    )
    
    # STAR summary rows
    star <- x[grepl("^N_", x$gene_id), ]
    
    cat(
      "\nFile:", rna.map$file_name[i],
      "\nTotal unstranded counts:",
      sum(x$unstranded[!grepl("^N_", x$gene_id)]),
      "\n"
    )
    
    print(
      star[, c("gene_id", "unstranded")]
    )
  }
}

## For duplicate sample.submitter_id
qc <- data.frame()

for (s in dup.samples) {
  
  rows <- which(rna.map$sample.submitter_id == s)
  
  for (i in rows) {
    
    f <- rna.files[
      basename(rna.files) == rna.map$file_name[i]
    ]
    
    x <- read.delim(
      f,
      header = TRUE,
      sep = "\t",
      comment.char = "#",
      check.names = FALSE
    )
    
    star <- x[grepl("^N_", x$gene_id), ]
    gene <- x[!grepl("^N_", x$gene_id), ]
    
    qc <- rbind(
      qc,
      data.frame(
        sample = s,
        file_name = basename(f),
        assigned = sum(gene$unstranded),
        multimapping = star$unstranded[
          star$gene_id == "N_multimapping"
        ],
        noFeature = star$unstranded[
          star$gene_id == "N_noFeature"
        ],
        unmapped = star$unstranded[
          star$gene_id == "N_unmapped"
        ],
        ambiguous = star$unstranded[
          star$gene_id == "N_ambiguous"
        ]
      )
    )
  }
}

qc

qc$assigned_fraction <- qc$assigned / (
  qc$assigned +
    qc$multimapping +
    qc$noFeature +
    qc$unmapped +
    qc$ambiguous
)

best_dup <- qc[
  ave(
    qc$assigned_fraction,
    qc$sample,
    FUN = function(x) x == max(x)
  ),
]

bad_files <- qc$file_name[
  qc$assigned_fraction <
    ave(
      qc$assigned_fraction,
      qc$sample,
      FUN = max
    )
]

rna.map.clean <- rna.map[
  !rna.map$file_name %in% bad_files,
]

nrow(rna.map.clean)
length(unique(rna.map.clean$sample.submitter_id))

patient.counts <- table(rna.map.clean$cases.submitter_id)

multi.patients <- names(
  patient.counts[patient.counts > 1]
)

length(multi.patients)

rna.multi <- rna.map.clean[
  rna.map.clean$cases.submitter_id %in% multi.patients,
  c(
    "cases.submitter_id",
    "sample.submitter_id",
    "file_name"
  )
]

rna.multi[
  order(
    rna.multi$cases.submitter_id,
    rna.multi$sample.submitter_id
  ),
]

# 2) get clinical table (GDC clinical or legacy file)
clin <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")  # inspect columns

grep(
  "days|death|vital|survival",
  colnames(clin),
  ignore.case = TRUE,
  value = TRUE
)

clin[, grep(
  "days|death|vital|survival",
  colnames(clin),
  ignore.case = TRUE
)]

saveRDS(
  clin,
  file = "TCGA_BRCA_clinical.rds",
  compress = FALSE
)


query.wes <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)

getResults(query.wes)

getResults(query.wes)[, c(
  "file_name",
  "data_type",
  "experimental_strategy",
  "workflow_type"
)]

GDCdownload(
  query.wes,
  method = "api",
  files.per.chunk = 5,
  directory = "TCGA_BRCA"
)

maf <- GDCprepare(
  query = query.wes,
  directory = "TCGA_BRCA"
)


query.wes <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open",
  experimental.strategy = "WXS"
)

nrow(getResults(query.wes))


wes.results <- getResults(query.wes)

dim(wes.results)

table(wes.results$data_type)
table(wes.results$experimental_strategy)
table(wes.results$data_format)

length(unique(wes.results$project))


maf$patient_id <- substr(
  maf$Tumor_Sample_Barcode,
  1,
  12
)

dim(maf)

length(unique(maf$patient_id))

table(maf$Variant_Classification)

## keep protein coding changes
keep <- maf$Variant_Classification %in% c(
  "Frame_Shift_Del",
  "Frame_Shift_Ins",
  "Missense_Mutation",
  "Nonsense_Mutation",
  "Splice_Site",
  "In_Frame_Del",
  "In_Frame_Ins",
  "Translation_Start_Site",
  "Nonstop_Mutation"
)

maf2 <- maf[keep, ]

dim(maf2)

length(unique(maf2$patient_id))

table(maf2$Variant_Classification)


gene_counts <- sort(
  table(maf2$Hugo_Symbol),
  decreasing = TRUE
)

length(gene_counts)

head(gene_counts, 30)



## keep 1095 patients with uniqu WES and RNA
# Keep 01A when a patient has multiple primary-tumor samples
rna.final <- do.call(
  rbind,
  lapply(
    split(rna.map.clean, rna.map.clean$cases.submitter_id),
    function(x) {
      if (nrow(x) == 1) {
        return(x)
      }
      
      x[x$sample.submitter_id == 
          paste0(unique(x$cases.submitter_id), "-01A"), ]
    }
  )
)

rownames(rna.final) <- NULL

dim(rna.final)
length(unique(rna.final$cases.submitter_id))

table(table(rna.final$cases.submitter_id))

## overlap RNA and WES samples
rna.patients <- unique(rna.final$cases.submitter_id)
wes.patients <- unique(maf$patient_id)

matched.patients <- intersect(
  rna.patients,
  wes.patients
)

length(rna.patients)
length(wes.patients)
length(matched.patients)

##create the final 1,095-patient mapping
# Keep one sample per patient.
# For patients with multiple samples, retain 01A.

rna.final <- do.call(
  rbind,
  lapply(
    split(rna.map.clean, rna.map.clean$cases.submitter_id),
    function(x) {
      
      if (nrow(x) == 1) {
        return(x)
      }
      
      x[x$sample.submitter_id == 
          paste0(unique(x$cases.submitter_id), "-01A"), ]
    }
  )
)

rownames(rna.final) <- NULL

dim(rna.final)
length(unique(rna.final$cases.submitter_id))


## Build the expression matrix
# Get local file paths corresponding to the final 1095 samples
rna.final$file_path <- rna.files[
  match(
    rna.final$file_name,
    basename(rna.files)
  )
]

# Read first file to establish gene information
x <- read.delim(
  rna.final$file_path[1],
  header = TRUE,
  sep = "\t",
  comment.char = "#",
  check.names = FALSE
)

# Remove STAR summary rows
x <- x[!grepl("^N_", x$gene_id), ]

# Keep gene ID, gene name and raw unstranded counts
gene_id <- x$gene_id
gene_name <- x$gene_name

# Initialize matrix
rna.counts <- matrix(
  NA_integer_,
  nrow = length(gene_id),
  ncol = nrow(rna.final)
)

rownames(rna.counts) <- gene_id
colnames(rna.counts) <- rna.final$cases.submitter_id

# First sample
rna.counts[, 1] <- x$unstranded

rm(x)


for (i in 2:nrow(rna.final)) {
  
  if (i %% 100 == 0) {
    cat("Processing", i, "of", nrow(rna.final), "\n")
  }
  
  x <- read.delim(
    rna.final$file_path[i],
    header = TRUE,
    sep = "\t",
    comment.char = "#",
    check.names = FALSE
  )
  
  x <- x[!grepl("^N_", x$gene_id), ]
  
  # Make sure gene order is identical
  if (!identical(x$gene_id, gene_id)) {
    x <- x[match(gene_id, x$gene_id), ]
  }
  
  rna.counts[, i] <- x$unstranded
  
  rm(x)
}

dim(rna.counts)

sum(is.na(rna.counts))

length(unique(rownames(rna.counts)))

saveRDS(
  rna.counts,
  file = "rna_counts_1095_patients.rds",
  compress = FALSE
)

gene.annotation <- data.frame(
  gene_id = gene_id,
  gene_name = gene_name,
  stringsAsFactors = FALSE
)

saveRDS(
  gene.annotation,
  file = "na_gene_annotation.rds",
  compress = FALSE
)


saveRDS(
  maf,
  file = "TCGA_BRCA_WES_MAF.rds",
  compress = FALSE
)

# Save the protein-altering subset
saveRDS(
  maf2,
  file = "TCGA_BRCA_WES_protein_altering_MAF.rds",
  compress = FALSE
)




## download TCGA clinical supplement data
library(TCGAbiolinks)

query.brca.clin <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR Biotab"
)

GDCdownload(query.brca.clin)

clinical.BCR <- GDCprepare(query.brca.clin)

names(clinical.BCR)


brca.patient <- clinical.BCR$clinical_patient_brca

dim(brca.patient)

colnames(brca.patient)

brca.patient.use <- brca.patient[ -c(1:2),]

brca.patient.final.cols <- brca.patient.use[, c(
  "bcr_patient_barcode",
  "age_at_diagnosis",
  "er_status_by_ihc",
  "pr_status_by_ihc",
  "her2_status_by_ihc",
  "ajcc_pathologic_tumor_stage",
  "ajcc_tumor_pathologic_pt",
  "ajcc_nodes_pathologic_pn",
  "ajcc_metastasis_pathologic_pm",
  "vital_status",
  "death_days_to",
  "last_contact_days_to",
  "days_to_patient_progression_free",
  "days_to_tumor_progression"
)]

write.csv(brca.patient.final.cols,file="TCGA_BRCA_subset_of_clinical_supplement_information_all_patients.csv",col.names = T,row.names = F, quote = FALSE)

write.csv(brca.patient.use,file="TCGA_BRCA_all_clinical_supplement_information_patients.csv",col.names = T,row.names = F, quote = FALSE)


