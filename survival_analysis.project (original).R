# ---------------------- #
# Survival analysis: KM analtsis, logrank test, cox regression
# ---------------------- #

# Install packages and load librarys
BiocManager::install("ggsurvfit")
BiocManager::install("survminer")
BiocManager::install("survival", force = TRUE)
BiocManager::install("TCGAbiolinks")
install.packages("tidyverse")

library(tidyverse)
library(survival)
library(ggsurvfit)
library(survminer)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)


# Setting up environment ===================================================
# Clean environment
rm(list = ls(all.names = TRUE)) # clear all objects including hidden objects
gc() # free up memory and report the memory usagex

# avoid truncated output in R console and scientific notation
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) 


#Download clinical data in TCGA (Thyroid carcinoma: TCGA-THCA cohort)===========================================
getProjectSummary("TCGA-THCA")

data_clinic <- GDCquery_clinic("TCGA-THCA") # download the clinical data
any(colnames(data_clinic) %in% c("vital_status", "days_to_last_follow_up", "days_to_death", "age_at_index", "gender",
                                 "ajcc_pathologic_t", "ajcc_pathologic_m", "ajcc_pathologic_n","ajcc_pathologic_stage"))
which(colnames(data_clinic) %in% c("vital_status", "days_to_last_follow_up", "days_to_death", "age_at_index", "gender",
                                   "ajcc_pathologic_t", "ajcc_pathologic_m", "ajcc_pathologic_n","ajcc_pathologic_stage"))
data_clinic[,c(4,9,20,23,24,42,44,45,52)]


# Check the statics and missing values
table(data_clinic$vital_status, useNA = "ifany") #event size is too small, hard to apply the Cox regression(specially for multivariate)
table(data_clinic$days_to_last_follow_up, useNA = "ifany")  
table(data_clinic$days_to_death, useNA = "ifany")
table(data_clinic$gender, useNA = "ifany")
table(data_clinic$ajcc_pathologic_stage, useNA = "ifany")
summary(data_clinic$age_at_index)
summary(data_clinic$days_to_last_follow_up) # consider time of years
summary(data_clinic$days_to_death)


# Variable transformation
data_clinic$vital_binary <- ifelse(data_clinic$vital_status == "Alive", FALSE, TRUE)
data_clinic$Age <- data_clinic$age_at_index
data_clinic$Gender <- ifelse(data_clinic$gender == "male", 0, 1)
data_clinic$Stage <- as.numeric(factor(data_clinic$ajcc_pathologic_stage, 
                                       levels = c("Stage I", "Stage II", "Stage III", "Stage IV", "Stage IVA", "Stage IVC"),
                                       labels = c(1, 2, 3, 4, 5, 6)))
data_clinic$Drug_Treatment <- ifelse(data_clinic$treatments_pharmaceutical_treatment_or_therapy == "yes", 0,
                                     ifelse(data_clinic$treatments_pharmaceutical_treatment_or_therapy == "no", 1, NA))
data_clinic$overall_survival <- ifelse(data_clinic$vital_status == "Alive",
                                       data_clinic$days_to_last_follow_up,
                                       data_clinic$days_to_death)

table(data_clinic$vital_binary, useNA = "ifany")
table(data_clinic$Gender, useNA = "ifany")
table(data_clinic$Stage, useNA = "ifany")
table(data_clinic$Drug_Treatment, useNA = "ifany")
summary(data_clinic$overall_survival)


#Download Transcriptome Profiling data in TCGA (Thyroid carcinoma: TCGA-THCA cohort)===========================================
query_gene <- GDCquery(project = "TCGA-THCA", #Building a Query
                       data.category = "Transcriptome Profiling",
                       workflow.type = "STAR - Counts",
                       data.type = "Gene Expression Quantification",
                       sample.type = "Primary Tumor")

output_THCA <- getResults(query_gene)
summary(output_THCA)

downloading_directory <- file.path(getwd(), "GDCdata")
GDCdownload(query_gene, directory = "GDCdata", files.per.chunk = 50) # get gene expression data

data_dir <- "C://Users//user//OneDrive//바탕 화면//새 폴더//GDCdata"
data_gene <- GDCprepare(query_gene, directory = data_dir) # merge individual raw data to one object

#Extract gene 
THCAexpr_matrix <- assay(data_gene, "unstranded")
THCAexpr_matrix[1:10, 1:10]
gene_metadata <- as.data.frame(rowData(data_gene))
coldata <- as.data.frame(colData(data_gene))


# Transform counts using vst===========================================================================================
# Set up countData
dds <- DESeqDataSetFromMatrix(countData = THCAexpr_matrix,
                              colData = coldata,
                              design = ~1)

# Removing low expressed genes with sum total of reads (>=10)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# vst
vsd <- vst(dds, blind=FALSE)
THCAexpr_matrix_vst <- assay(vsd)
THCAexpr_matrix_vst[1:5, 1:5]


# Get data interested genes and gene metadata
THCA_AKT1 <- THCAexpr_matrix_vst %>%
  as.data.frame() %>%
  rownames_to_column(var = 'gene_id') %>%
  gather(key = 'case_id', value = 'counts', -gene_id) %>%
  left_join(., gene_metadata, by = "gene_id") %>%
  filter(gene_name == "AKT1")

THCA_TNF <- THCAexpr_matrix_vst %>%
  as.data.frame() %>%
  rownames_to_column(var = 'gene_id') %>%
  gather(key = 'case_id', value = 'counts', -gene_id) %>%
  left_join(., gene_metadata, by = "gene_id") %>%
  filter(gene_name == "TNF")

THCA_HSP90AA1 <- THCAexpr_matrix_vst %>%
  as.data.frame() %>%
  rownames_to_column(var = 'gene_id') %>%
  gather(key = 'case_id', value = 'counts', -gene_id) %>%
  left_join(., gene_metadata, by = "gene_id") %>%
  filter(gene_name == "HSP90AA1")

THCA_IL6 <- THCAexpr_matrix_vst %>%
  as.data.frame() %>%
  rownames_to_column(var = 'gene_id') %>%
  gather(key = 'case_id', value = 'counts', -gene_id) %>%
  left_join(., gene_metadata, by = "gene_id") %>%
  filter(gene_name == "IL6")

THCA_EGFR <- THCAexpr_matrix_vst %>%
  as.data.frame() %>%
  rownames_to_column(var = 'gene_id') %>%
  gather(key = 'case_id', value = 'counts', -gene_id) %>%
  left_join(., gene_metadata, by = "gene_id") %>%
  filter(gene_name == "EGFR")


# Get median value
median_AKT1 <- median(THCA_AKT1$counts)
median_TNF <- median(THCA_TNF$counts)
median_HSP90AA1 <- median(THCA_HSP90AA1$counts)
median_IL6 <- median(THCA_IL6$counts)
median_EGFR <- median(THCA_EGFR$counts)


# Classify the high and low expression groups using median value
THCA_AKT1$Strata <- as.factor(ifelse(THCA_AKT1$counts >= median_AKT1, "AKT1_HIGH", "AKT1_LOW"))
THCA_TNF$Strata <- as.factor(ifelse(THCA_TNF$counts >= median_TNF, "TNF_HIGH", "TNF_LOW"))
THCA_HSP90AA1$Strata <- as.factor(ifelse(THCA_HSP90AA1$counts >= median_HSP90AA1, "HSP90AA1_HIGH", "HSP90AA1_LOW"))
THCA_IL6$Strata <- as.factor(ifelse(THCA_IL6$counts >= median_IL6, "IL6_HIGH", "IL6_LOW"))
THCA_EGFR$Strata <- as.factor(ifelse(THCA_EGFR$counts >= median_EGFR, "EGFR_HIGH", "EGFR_LOW"))

# Chage the ref group (if you need)
THCA_AKT1$Strata <- relevel(THCA_AKT1$Strata, ref = "AKT1_LOW")
THCA_TNF$Strata <- relevel(THCA_TNF$Strata, ref = "TNF_LOW")
THCA_HSP90AA1$Strata <- relevel(THCA_HSP90AA1$Strata, ref = "HSP90AA1_LOW")
THCA_IL6$Strata <- relevel(THCA_IL6$Strata, ref = "IL6_LOW")
THCA_EGFR$Strata <- relevel(THCA_EGFR$Strata, ref = "EGFR_LOW")


# Add clinical data
THCA_AKT1$case_id <- gsub('-01.*', '', THCA_AKT1$case_id)
THCA_TNF$case_id <- gsub('-01.*', '', THCA_TNF$case_id)
THCA_HSP90AA1$case_id <- gsub('-01.*', '', THCA_HSP90AA1$case_id)
THCA_IL6$case_id <- gsub('-01.*', '', THCA_IL6$case_id)
THCA_EGFR$case_id <- gsub('-01.*', '', THCA_EGFR$case_id)

THCA_AKT1 <- merge(THCA_AKT1, data_clinic, by.x = 'case_id', by.y = 'submitter_id')
THCA_TNF <- merge(THCA_TNF, data_clinic, by.x = 'case_id', by.y = 'submitter_id')
THCA_HSP90AA1 <- merge(THCA_HSP90AA1, data_clinic, by.x = 'case_id', by.y = 'submitter_id')
THCA_IL6 <- merge(THCA_IL6, data_clinic, by.x = 'case_id', by.y = 'submitter_id')
THCA_EGFR <- merge(THCA_EGFR, data_clinic, by.x = 'case_id', by.y = 'submitter_id')

write.csv(THCA_AKT1, "THCA_AKT1.csv")
write.csv(THCA_TNF, "THCA_TNF.csv")
write.csv(THCA_HSP90AA1, "THCA_HSP90AA1.csv")
write.csv(THCA_IL6, "THCA_IL6.csv")
write.csv(THCA_EGFR, "THCA_EGFR.csv")

# Fitting the Curve=================================================================================
fit_AKT1 <- survfit(Surv(overall_survival, vital_binary) ~ Strata, data = THCA_AKT1)
fit_TNF <- survfit(Surv(overall_survival, vital_binary) ~ Strata, data = THCA_TNF)
fit_HSP90AA1 <- survfit(Surv(overall_survival, vital_binary) ~ Strata, data = THCA_HSP90AA1)
fit_IL6 <- survfit(Surv(overall_survival, vital_binary) ~ Strata, data = THCA_IL6)
fit_EGFR <- survfit(Surv(overall_survival, vital_binary) ~ Strata, data = THCA_EGFR)


# Log-rank test and extract p-value for annotation to KM plot
logrank_AKT1 <- survdiff(formula = Surv(overall_survival, vital_binary) ~ Strata, data = THCA_AKT1)
logrank_TNF <- survdiff(formula = Surv(overall_survival, vital_binary) ~ Strata, data = THCA_TNF)
logrank_HSP90AA1 <- survdiff(formula = Surv(overall_survival, vital_binary) ~ Strata, data = THCA_HSP90AA1)
logrank_IL6 <- survdiff(formula = Surv(overall_survival, vital_binary) ~ Strata, data = THCA_IL6)
logrank_EGFR <- survdiff(formula = Surv(overall_survival, vital_binary) ~ Strata, data = THCA_EGFR)

pval_AKT1 <- paste0("logrank P = ", signif(1 - pchisq(logrank_AKT1$chisq, length(logrank_AKT1$n) - 1), 3))
pval_TNF <- paste0("logrank P = ", signif(1 - pchisq(logrank_TNF$chisq, length(logrank_TNF$n) - 1), 3))
pval_HSP90AA1 <- paste0("logrank P = ", signif(1 - pchisq(logrank_HSP90AA1$chisq, length(logrank_HSP90AA1$n) - 1), 3))
pval_IL6 <- paste0("logrank P = ", signif(1 - pchisq(logrank_IL6$chisq, length(logrank_IL6$n) - 1), 3))
pval_EGFR <- paste0("logrank P = ", signif(1 - pchisq(logrank_EGFR$chisq, length(logrank_EGFR$n) - 1), 3))

# Cox regression: Univariate analysis =================================================================
cox_AKT1_uni <- coxph(Surv(overall_survival, vital_binary) ~ Strata, data = THCA_AKT1)
cox_TNF_uni <- coxph(Surv(overall_survival, vital_binary) ~ Strata, data = THCA_TNF)
cox_HSP90AA1_uni <- coxph(Surv(overall_survival, vital_binary) ~ Strata, data = THCA_HSP90AA1)
cox_IL6_uni <- coxph(Surv(overall_survival, vital_binary) ~ Strata, data = THCA_IL6)
cox_EGFR_uni <- coxph(Surv(overall_survival, vital_binary) ~ Strata, data = THCA_EGFR)

# Extract HR ratio for annotation to KM plot
hr_AKT1 <- paste0("HR = ", round(summary(cox_AKT1_uni)$coefficients[,"exp(coef)"], 2), " (", round(summary(cox_AKT1_uni)$conf.int[,"lower .95"], 2),
                  " - ", round(summary(cox_AKT1_uni)$conf.int[,"upper .95"], 2), ")")
hr_TNF <- paste0("HR = ", round(summary(cox_TNF_uni)$coefficients[,"exp(coef)"], 2), " (", round(summary(cox_TNF_uni)$conf.int[,"lower .95"], 2),
                 " - ", round(summary(cox_TNF_uni)$conf.int[,"upper .95"], 2), ")")
hr_HSP90AA1 <- paste0("HR = ", round(summary(cox_HSP90AA1_uni)$coefficients[,"exp(coef)"], 2), " (", round(summary(cox_HSP90AA1_uni)$conf.int[,"lower .95"], 2),
                      " - ", round(summary(cox_HSP90AA1_uni)$conf.int[,"upper .95"], 2), ")")
hr_IL6 <- paste0("HR = ", round(summary(cox_IL6_uni)$coefficients[,"exp(coef)"], 2), " (", round(summary(cox_IL6_uni)$conf.int[,"lower .95"], 2),
                 " - ", round(summary(cox_IL6_uni)$conf.int[,"upper .95"], 2), ")")
hr_EGFR <- paste0("HR = ", round(summary(cox_EGFR_uni)$coefficients[,"exp(coef)"], 2), " (", round(summary(cox_EGFR_uni)$conf.int[,"lower .95"], 2),
                  " - ", round(summary(cox_EGFR_uni)$conf.int[,"upper .95"], 2), ")")


# KM plot and saving the pdf files and annotation
plot_AKT1 <- ggsurvplot(fit_akt1,
                        data = THCA_AKT1,
                        xlab = "Time (days)",
                        pval = F,
                        risk.table = T,
                        legend.title = "",
                        legend.labs = c("AKT1_LOW", "AKT1_HIGH"),
                        palette = c("#00BFC4", "#F8766D"))
plot_AKT1$plot <- plot_AKT1$plot +
  ggplot2::annotate("text",
                    x = max(THCA_AKT1$overall_survival, na.rm = TRUE) * 0.75,
                    y = 0.15,
                    label = pval_AKT1,
                    size = 4,
                    hjust = 0) +
  ggplot2::annotate("text",
                    x = max(THCA_AKT1$overall_survival, na.rm = TRUE) * 0.75,
                    y = 0.10,
                    label = hr_AKT1,
                    size = 4,
                    hjust = 0)
pdf("KM plot_AKT1.pdf", width = 9, height = 8)
print(plot_AKT1)
dev.off()


plot_TNF <- ggsurvplot(fit_TNF,
                       data = THCA_TNF,
                       xlab = "Time (days)",
                       pval = F,
                       risk.table = T,
                       legend.title = "",
                       legend.labs = c("TNF_LOW", "TNF_HIGH"),
                       palette = c("#00BFC4", "#F8766D"))
plot_TNF$plot <- plot_TNF$plot +
  ggplot2::annotate("text",
                    x = max(THCA_TNF$overall_survival, na.rm = TRUE) * 0.75,
                    y = 0.15,
                    label = pval_TNF,
                    size = 4,
                    hjust = 0) +
  ggplot2::annotate("text",
                    x = max(THCA_TNF$overall_survival, na.rm = TRUE) * 0.75,
                    y = 0.10,
                    label = hr_TNF,
                    size = 4,
                    hjust = 0)
pdf("KM plot_TNF.pdf", width = 9, height = 8)
print(plot_TNF)
dev.off()


plot_HSP90AA1 <- ggsurvplot(fit_HSP90AA1,
                            data = THCA_HSP90AA1,
                            xlab = "Time (days)",
                            pval = F,
                            risk.table = T,
                            legend.title = "",
                            legend.labs = c("HSP90AA1_LOW", "HSP90AA1_HIGH"),
                            palette = c("#00BFC4", "#F8766D"))
plot_HSP90AA1$plot <- plot_HSP90AA1$plot +
  ggplot2::annotate("text",
                    x = max(THCA_HSP90AA1$overall_survival, na.rm = TRUE) * 0.75,
                    y = 0.15,
                    label = pval_HSP90AA1,
                    size = 4,
                    hjust = 0) +
  ggplot2::annotate("text",
                    x = max(THCA_HSP90AA1$overall_survival, na.rm = TRUE) * 0.75,
                    y = 0.10,
                    label = hr_HSP90AA1,
                    size = 4,
                    hjust = 0)
pdf("KM plot_HSP90AA1.pdf", width = 9, height = 8)
print(plot_HSP90AA1)
dev.off()


plot_IL6 <- ggsurvplot(fit_IL6,
                       data = THCA_IL6,
                       xlab = "Time (days)",
                       pval = F,
                       risk.table = T,
                       legend.title = "",
                       legend.labs = c("IL6_LOW", "IL6_HIGH"),
                       palette = c("#00BFC4", "#F8766D"))
plot_IL6$plot <- plot_IL6$plot +
  ggplot2::annotate("text",
                    x = max(THCA_IL6$overall_survival, na.rm = TRUE) * 0.75,
                    y = 0.15,
                    label = pval_IL6,
                    size = 4,
                    hjust = 0) +
  ggplot2::annotate("text",
                    x = max(THCA_IL6$overall_survival, na.rm = TRUE) * 0.75,
                    y = 0.10,
                    label = hr_IL6,
                    size = 4,
                    hjust = 0)
pdf("KM plot_IL6.pdf", width = 9, height = 8)
print(plot_IL6)
dev.off()


plot_EGFR <- ggsurvplot(fit_EGFR,
                        data = THCA_EGFR,
                        xlab = "Time (days)",
                        pval = F,
                        risk.table = T,
                        legend.title = "",
                        legend.labs = c("EGFR_LOW", "EGFR_HIGH"),
                        palette = c("#00BFC4", "#F8766D"))
plot_EGFR$plot <- plot_EGFR$plot +
  ggplot2::annotate("text",
                    x = max(THCA_EGFR$overall_survival, na.rm = TRUE) * 0.75,
                    y = 0.15,
                    label = pval_EGFR,
                    size = 4,
                    hjust = 0) +
  ggplot2::annotate("text",
                    x = max(THCA_EGFR$overall_survival, na.rm = TRUE) * 0.75,
                    y = 0.10,
                    label = hr_EGFR,
                    size = 4,
                    hjust = 0)
pdf("KM plot_EGFR.pdf", width = 9, height = 8)
print(plot_EGFR)
dev.off()


# Cox regression: Multivariate analysis =================================================================
cox_AKT1 <- coxph(Surv(overall_survival, vital_binary) ~ Gender + Age + Stage + Strata, data = THCA_AKT1)
cox_TNF <- coxph(Surv(overall_survival, vital_binary) ~ Gender + Age + Stage + Strata, data = THCA_TNF)
cox_HSP90AA1 <- coxph(Surv(overall_survival, vital_binary) ~ Gender + Age + Stage + Strata, data = THCA_HSP90AA1)
cox_IL6 <- coxph(Surv(overall_survival, vital_binary) ~ Gender + Age + Stage + Strata, data = THCA_IL6)
cox_EGFR <- coxph(Surv(overall_survival, vital_binary) ~ Gender + Age + Stage + Strata, data = THCA_EGFR)

summary(cox_AKT1)
summary(cox_TNF)
summary(cox_HSP90AA1)
summary(cox_IL6)
summary(cox_EGFR)


# Cox regression: proportional hazards assumption and tree_plot
pdf("PHA_forest_AKT1.pdf", width = 8, height = 8)
cox_test_AKT1 <- survival::cox.zph(cox_AKT1)
PHA_AKT1 <-survminer::ggcoxzph(cox_test_AKT1, font.axis=8,
                               font.tickslab = 6,
                               font.x = 8,
                               font.y = 8,
                               font.main = 8)
forest_AKT1 <- ggforest(cox_AKT1, data = THCA_AKT1)
print(forest_AKT1)
print(PHA_AKT1)
dev.off()
pdf("PHA_forest_TNF.pdf", width = 8, height = 8)
cox_test_TNF <- survival::cox.zph(cox_TNF)
PHA_TNF <-survminer::ggcoxzph(cox_test_TNF, font.axis=8,
                              font.tickslab = 6,
                              font.x = 8,
                              font.y = 8,
                              font.main = 8)
forest_TNF <- ggforest(cox_TNF, data = THCA_TNF)
print(forest_TNF)
print(PHA_TNF)
dev.off()
pdf("PHA_forest_HSP90AA1.pdf", width = 8, height = 8)
cox_test_HSP90AA1 <- survival::cox.zph(cox_HSP90AA1)
PHA_HSP90AA1 <-survminer::ggcoxzph(cox_test_HSP90AA1, font.axis=8,
                                   font.tickslab = 6,
                                   font.x = 8,
                                   font.y = 8,
                                   font.main = 8)
forest_HSP90AA1 <- ggforest(cox_HSP90AA1, data = THCA_HSP90AA1)
print(forest_HSP90AA1)
print(PHA_HSP90AA1)
dev.off()
pdf("PHA_forest_IL6.pdf", width = 8, height = 8)
cox_test_IL6 <- survival::cox.zph(cox_IL6)
PHA_IL6 <-survminer::ggcoxzph(cox_test_IL6, font.axis=8,
                              font.tickslab = 6,
                              font.x = 8,
                              font.y = 8,
                              font.main = 8)
forest_IL6 <- ggforest(cox_IL6, data = THCA_IL6)
print(forest_IL6)
print(PHA_IL6)
dev.off()
pdf("PHA_forest_EGFR.pdf", width = 8, height = 8)
cox_test_EGFR <- survival::cox.zph(cox_EGFR)
PHA_EGFR <-survminer::ggcoxzph(cox_test_EGFR, font.axis=8,
                               font.tickslab = 6,
                               font.x = 8,
                               font.y = 8,
                               font.main = 8)
forest_EGFR <- ggforest(cox_EGFR, data = THCA_EGFR)
print(forest_EGFR)
print(PHA_EGFR)
dev.off()


