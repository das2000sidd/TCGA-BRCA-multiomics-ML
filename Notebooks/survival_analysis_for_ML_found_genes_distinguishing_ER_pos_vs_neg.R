setwd("~/Desktop/Related_to_Haider_Lab/TCGA_TNBC")

library(dplyr)

clinical <- read.table(file="TCGA_BRCA_Clinica_information_with_vital_status_and_days_to_death.txt",header = T, stringsAsFactors = F, sep="\t")

#(base) siddharthadas@Mac TCGA_TNBC % grep -Ril "days_to_last_follow_up" .
#./TCGA_BRCA/TCGA_BRCA_clinical.rds
#./Fit_DL_models_TCGA_data.ipynb
#./TCGA_BRCA_Clinica_information_with_vital_status_and_days_to_death.csv
#./Clinical_data_for_1098_TCGA_BRCA_patients.txt
#./TCGA_BRCA_clinical.rds
# ./TCGA_BRCA_Clinica_information_with_vital_status_and_days_to_death.txt

relevant_genes <- read.csv(file="Potential_Candidate_genes_distinguishing_ER_positive_versus_ER_negative.csv",header = T, stringsAsFactors = F)


exp_data <- read.csv(file="TCGA_BRCA_VSD_transformed_expression_matrix_above_75th_quantile_variance_overlapping_patients.csv",header = T, stringsAsFactors = F)

rownames(exp_data) <- exp_data$X;

exp_data <- exp_data[, -c(1)];

head(exp_data)

rownames(exp_data)[1:5]

rownames(exp_data) <- sub("\\..*$", "", rownames(exp_data))

exp_relevant_genes <- exp_data[ relevant_genes$gene, ]


test_samples <- read.csv(file="Test_samples_IDs.csv",header = F)

test_samples <- test_samples$V1

test_samples[1:5]

test_samples <- gsub("-", ".", test_samples)


exp_data_test <- exp_data[, test_samples]

  
exp_data_test_relevant_genes <- exp_data_test[ relevant_genes$gene, ]


survival_data <- clinical[, c( "submitter_id","vital_status", "days_to_death.x", "days_to_last_follow_up")]


survival_data$OS_event <- ifelse(
  survival_data$vital_status == "Dead", 1, 0
)

survival_data$OS_time <- ifelse(
  survival_data$vital_status == "Dead",
  survival_data$days_to_death.x,
  survival_data$days_to_last_follow_up
)


head(survival_data)

survival_data$submitter_id_dot <- gsub( "-",  ".", survival_data$submitter_id)

survival_data_test <- survival_data[ survival_data$submitter_id_dot %in% colnames(exp_data_test_relevant_genes), ]

rownames(survival_data_test) <- survival_data_test$submitter_id_dot

survival_data_test_o <- survival_data_test[ colnames(exp_data_test_relevant_genes), ]

stopifnot(survival_data_test_o$submitter_id_dot == colnames(exp_data_test_relevant_genes))

exp_data_test_relevant_genes <- t(exp_data_test_relevant_genes)

survival_data_test_o_test <- cbind(survival_data_test_o, exp_data_test_relevant_genes[rownames(survival_data_test_o), "ENSG00000233622"])

colnames(survival_data_test_o_test)[8] <- "ENSG00000233622"

library(survival)

#fit <- coxph(
#  Surv(OS_time, OS_event) ~ ENSG00000233622,
#  data = survival_data_test_o_test
#)

#summary(fit)

coxph_genes_of_interest <- matrix(NA, nrow = 0, ncol = 5)
colnames(coxph_genes_of_interest) <- c(
  "Gene", "HR", "HR_lower", "HR_upper", "p_value"
)

for(gene in relevant_genes$gene){
  
  survival_data_test_o_test <- cbind(
    survival_data_test_o,
    expression = exp_data_test_relevant_genes[
      rownames(survival_data_test_o), gene
    ]
  )
  
  fit <- coxph(
    Surv(OS_time, OS_event) ~ expression,
    data = survival_data_test_o_test
  )
  
  s <- summary(fit)
  
  result <- data.frame(
    Gene = gene,
    HR = s$coefficients[, "exp(coef)"],
    HR_lower = s$conf.int[, "lower .95"],
    HR_upper = s$conf.int[, "upper .95"],
    p_value = s$coefficients[, "Pr(>|z|)"]
  )
  
  coxph_genes_of_interest <- rbind(
    coxph_genes_of_interest,
    result
  )
}

coxph_genes_of_interest <- as.data.frame(coxph_genes_of_interest)

coxph_genes_of_interest$p_value <- as.numeric(
  coxph_genes_of_interest$p_value
)

coxph_genes_of_interest$FDR <- p.adjust(
  coxph_genes_of_interest$p_value,
  method = "BH"
)


gene_symbol <- relevant_genes[,c("gene", "SYMBOL")]

gene_symbol <- gene_symbol[ gene_symbol$SYMBOL!= "", ]


coxph_genes_of_interest_symbol <- left_join(coxph_genes_of_interest,gene_symbol, by = c("Gene" = "gene"))


library(survival)

gene <- "ENSG00000185008"

survival_km <- survival_data_test_o

survival_km$expression <- exp_data_test_relevant_genes[
  rownames(survival_km),
  gene
]

survival_km$expression_group <- ifelse(
  survival_km$expression >= median(
    survival_km$expression,
    na.rm = TRUE
  ),
  "High",
  "Low"
)

survival_km$expression_group <- factor(
  survival_km$expression_group,
  levels = c("Low", "High")
)


km_fit <- survfit(
  Surv(OS_time, OS_event) ~ expression_group,
  data = survival_km
)

summary(km_fit)

library(survminer)

pdf(
  "ENSG00000185008_survival_curve.pdf",
  width = 7,
  height = 7
)

print(
  ggsurvplot(
    km_fit,
    data = survival_km,
    pval = TRUE,
    risk.table = TRUE,
    conf.int = TRUE,
    xlab = "Overall survival (days)",
    ylab = "Survival probability",
    legend.title = gene,
    legend.labs = c("Low", "High")
  )
)

dev.off()



write.table(coxph_genes_of_interest,"COX_PRPOPORTIONAL_HAZARD_MODEL_RESULST_FOR_GENES_EXPRESSION_DISTINGUSHING_ER_POS_VERSUS_NEG.txt", col.names = T,row.names = F,quote = F, sep="\t")




