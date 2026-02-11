# ---------------------- #
# Survival analysis: Kaplan_Meier analysis, logrank test, cox regression, ROC
# ---------------------- #
# Install packages and load librarys
BiocManager::install("ggsurvfit", dependency = TRUE, force = TRUE)
BiocManager::install("survminer", dependency = TRUE)
BiocManager::install("survival", dependency = TRUE, force = TRUE)
BiocManager::install("TCGAbiolinks", dependency = TRUE)
install.packages("tidyverse", dependency = TRUE)
install.packages("timeROC", dependency = TRUE)
install.packages("limma", dependency = TRUE)
install.packages("pheatmap", dependency = TRUE)

library(timeROC)
library(tidyverse)
library(survival)
library(ggsurvfit)
library(survminer)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(limma)
library(pheatmap)

# Download clinical data in TCGA (Lung carcinoma: TCGA-LUAD cohort)
getProjectSummary("TCGA-LUAD")
clinic <- GDCquery_clinic("TCGA-LUAD") # download the clinical data

# Handling Variable level
clinic$Vital_status <- ifelse(clinic$vital_status == "Alive", FALSE, TRUE)

clinic$Age <- as.numeric(clinic$age_at_index)

clinic$Gender <- factor(ifelse(clinic$gender=="male",0,
                        ifelse(clinic$gender=="female",1,NA_integer_)), levels=c(0,1), labels=c("Male","Female"))

clinic$Smoking <- factor(ifelse(clinic$tobacco_smoking_status=="Lifelong Non-Smoker","Never",
                              ifelse(clinic$tobacco_smoking_status=="Current Reformed Smoker for > 15 yrs","Former>15",
                              ifelse(clinic$tobacco_smoking_status=="Current Reformed Smoker for < or = 15 yrs","Former<=15",
                              ifelse(clinic$tobacco_smoking_status=="Current Smoker","Current",
                              ifelse(clinic$tobacco_smoking_status %in% c("Unknown","Not Reported", "Current Reformed Smoker, Duration Not Specified"),"Unknown",NA_character_))))), 
                              levels=c("Never","Former>15","Former<=15","Current","Unknown"))

clinic$Stage <- factor(ifelse(clinic$ajcc_pathologic_stage %in% c("Stage I","Stage IA","Stage IB"), "I",
                       ifelse(clinic$ajcc_pathologic_stage %in% c("Stage II","Stage IIA","Stage IIB"), "II",
                       ifelse(clinic$ajcc_pathologic_stage %in% c("Stage III","Stage IIIA","Stage IIIB","Stage IIIC"), "III",
                       ifelse(clinic$ajcc_pathologic_stage %in% c("Stage IV","Stage IVA","Stage IVB","Stage IVC"), "IV", NA_character_)))),
                       levels = c("I","II","III","IV"))

clinic$M_stage <- factor(ifelse(clinic$ajcc_pathologic_m=="M0","M0",
                         ifelse(clinic$ajcc_pathologic_m %in% c("M1","M1a","M1b"),"M1",
                         ifelse(clinic$ajcc_pathologic_m=="MX","MX",NA_character_))),
                         levels=c("M0","M1","MX"))


clinic$N_stage <- factor(ifelse(clinic$ajcc_pathologic_n=="N0","N0",
                         ifelse(clinic$ajcc_pathologic_n=="N1","N1",
                         ifelse(clinic$ajcc_pathologic_n %in% c("N2","N3"), "N2/N3",
                         ifelse(clinic$ajcc_pathologic_n=="NX","NX",NA_character_)))),
                         levels=c("N0","N1","N2/N3","NX"))

# T_stage
clinic$T_stage <- factor(ifelse(clinic$ajcc_pathologic_t %in% c("T1","T1a","T1b"),"T1",
                         ifelse(clinic$ajcc_pathologic_t %in% c("T2","T2a","T2b"),"T2",
                         ifelse(clinic$ajcc_pathologic_t %in% c("T3","T4"), "T3/T4",
                         ifelse(clinic$ajcc_pathologic_t=="TX","TX",NA_character_)))),
                         levels=c("T1","T2","T3/T4","TX"))

clinic$Overall_survival <- ifelse(clinic$vital_status == "Alive", clinic$days_to_last_follow_up, clinic$days_to_death)
clinic$Overall_survival <- round(clinic$Overall_survival / 365, 2)

vars <- c("Vital_status", "Gender", "Smoking", "Stage", "M_stage", "N_stage", "T_stage")

lapply(vars, function(v) {print(table(clinic[[v]], useNA = "ifany"))})
summary(clinic$Age)
summary(clinic$Overall_survival)

#Download Transcriptome Profiling data in TCGA
query_gene <- GDCquery(project = "TCGA-LUAD", #Building a Query
                       data.category = "Transcriptome Profiling",
                       workflow.type = "STAR - Counts",
                       data.type = "Gene Expression Quantification",
                       sample.type = "Primary Tumor")
LUAD <- getResults(query_gene)
summary(LUAD)

path <- file.path(getwd(), "GDCdata")
GDCdownload(query_gene, directory = "GDCdata", files.per.chunk = 50) # get gene expression data
data_gene <- GDCprepare(query_gene, directory = path) # merge individual raw data to one object

#Extract gene 
exprs <- assay(data_gene, "unstranded")
exprs[1:10, 1:10]
metadata <- as.data.frame(rowData(data_gene))
coldata <- as.data.frame(colData(data_gene))

# Transform counts using vst
# Set up countData
dds <- DESeqDataSetFromMatrix(countData = exprs, colData = coldata, design = ~1)

# Removing low expressed genes with sum total of reads (>=10)
keep <- rowSums(counts(dds) >= 10) >= (ncol(dds) / 2)
dds <- dds[keep,]

# vst
vsd <- vst(dds, blind=FALSE)
exprs_vst <- assay(vsd)
hist(exprs_vst)

# Get data interested genes and metadata
target_genes <- c("CCNB1", "CDK1", "CCNA2", "AURKA", "CHEK1")
target_meta <- exprs_vst %>%
  as.data.frame() %>%
  rownames_to_column(var = 'gene_id') %>%
  pivot_longer(-gene_id, names_to = "case_id", values_to = "counts") %>%
  left_join(., metadata, by = "gene_id") %>%
  filter(gene_name %in% target_genes) %>%
  mutate(case_id = sub("-01.*", "", case_id)) 

## Collapse duplicates by averaging to a single row per
# Mapping information of targets for collapsing the rows
map_target <- target_meta %>% distinct(gene_id, gene_name)
map_target <- setNames(map_target$gene_name, map_target$gene_id)

# Change rownames to gene names
expr_target <- exprs_vst[rownames(exprs_vst) %in% names(map_target), , drop = FALSE]
rownames(expr_target) <- map_target[rownames(expr_target)]
ncol(expr_target)

# Averaging using limma
pid <- sub("-01.*", "", colnames(expr_target))
expr_target <- t(limma::avereps(t(expr_target), ID = pid, FUN = mean))
colnames(expr_target) <- unique(pid)
ncol(expr_target)

# Change data frame to long format
target_meta <- as.data.frame(expr_target) %>%
  tibble::rownames_to_column("gene_name") %>%
  pivot_longer(-gene_name, names_to = "case_id", values_to = "counts")

# Add clinical information
target_meta <- merge(target_meta, clinic, by.x = "case_id", by.y = "submitter_id")
stopifnot(sum(duplicated(target_meta[, c("case_id","gene_name")])) == 0)

# Classify the high and low expression groups using median value
target_meta <- target_meta %>%
  group_by(gene_name) %>%
  mutate(
    median_exp = median(counts, na.rm=TRUE),
    exp_group = factor(ifelse(counts >= median_exp, "HIGH", "LOW"), levels = c("LOW", "HIGH"))
  ) %>%
  ungroup()


# ---------------------- #
# 1) Kaplan-meier curve with log-rank test
# ---------------------- #
# Path for saving plots
result_dir <- file.path(getwd(), "result_dir")
if (!dir.exists(result_dir)) {
  dir.create(result_dir, recursive = TRUE)
}

for (g in target_genes) {
  
  dat_g <- subset(target_meta, gene_name == g)
  
  # Fitting the Curve
  fit <- survfit(Surv(Overall_survival, Vital_status) ~ exp_group, data = dat_g)
  
  # Log-rank test and extract p-value for annotation to KM plot
  logrank <- survdiff(formula = Surv(Overall_survival, Vital_status) ~ exp_group, data = dat_g)
  pval <- paste0("logrank P = ",
                 signif(1 - pchisq(logrank$chisq, length(logrank$n) - 1), 3))
  
  # Cox regression: Univariate analysis
  cox_uni <- coxph(Surv(Overall_survival, Vital_status) ~ exp_group, data = dat_g)
  
  # Save results in Global Env
  assign(paste0("logrank_", g), logrank, envir = .GlobalEnv)
  assign(paste0("cox_", g), cox_uni, envir = .GlobalEnv)
  
  # Extract HR ratio for annotation to KM plot
  HR <- paste0("HR = ",
               round(summary(cox_uni)$coefficients[,"exp(coef)"], 2), " (",
               round(summary(cox_uni)$conf.int[,"lower .95"], 2), " - ",
               round(summary(cox_uni)$conf.int[,"upper .95"], 2), ")")
  
  # KM plot and saving the pdf files and annotation
  plot_surv <- ggsurvplot(fit,
                          data = dat_g,
                          xlab = "Time (years)",
                          pval = F,
                          risk.table = T,
                          legend.title = "",
                          legend.labs = c(paste0(g, "_LOW"), paste0(g, "_HIGH")),
                          palette = c("#00BFC4", "#F8766D"),
                          ggtheme = theme_classic(base_size = 15),        
                          font.x = c(15),                         
                          font.y = c(15),                        
                          font.tickslab = 15,                            
                          risk.table.fontsize = 5.4)
  
  pdf(file.path(result_dir, paste0("KM_plot_", g, ".pdf")), width = 9, height = 8)
  plot_surv$plot <- plot_surv$plot +
    ggplot2::annotate("text",
                      x = max(dat_g$Overall_survival, na.rm = TRUE) * 0.75,
                      y = 0.95,
                      label = pval,
                      size = 5,
                      hjust = 0) +
    ggplot2::annotate("text",
                      x = max(dat_g$Overall_survival, na.rm = TRUE) * 0.75,
                      y = 0.88,
                      label = HR,
                      size = 5,
                      hjust = 0)
  print(plot_surv)
  dev.off()
}


# ---------------------- #
# 2) Multivariate Cox with target genes for deriving β and compute Risk Score
# ---------------------- #
# Clinical data in wide format
clin_cols <- c("Overall_survival","Vital_status","Gender","Age","Stage", "Smoking","M_stage","N_stage","T_stage")
clin_wide <- target_meta %>%
  select(all_of(c("case_id", clin_cols))) %>%
  distinct(case_id, .keep_all = TRUE)

# Expression data in wide format
expr_wide <- target_meta %>%
  select(case_id, gene_name, counts) %>%
  pivot_wider(names_from = gene_name, values_from = counts)

# Merge Clinical data + Expression data by case_id
cox_df <- clin_wide %>% inner_join(expr_wide, by = "case_id")
rownames(cox_df) <- cox_df$case_id

# Build formula for multivariate Cox model (continuous expression values)
form_genes <- as.formula(paste0("Surv(Overall_survival, Vital_status) ~ ", paste(target_genes, collapse = " + ")))

# Fit the multivariate Cox model
cox_multi <- coxph(form_genes, data = cox_df, x = TRUE)
summary(cox_multi)

# Checking the proportional hazards (PH) assumption
ph <- cox.zph(cox_multi)
pdf(file.path(result_dir, "PH_sig.pdf"), width = 9, height = 8)
PH_sig <-survminer::ggcoxzph(ph, font.axis=8,
                                 font.tickslab = 6,
                                 font.x = 8,
                                 font.y = 8,
                                 font.main = 8)
print(PH_sig)
dev.off()

# Calculate Risk Score from β coefficients (recommended: standardize as z-score)
betas <- coef(cox_multi)[target_genes]
Risk_score <- as.numeric(as.matrix(cox_df[, target_genes, drop=FALSE]) %*% betas)
cox_df$Risk_score <- as.numeric(scale(Risk_score))  # z-score


# ---------------------- #
# 3) Cox model with Risk Score + clinical covariates (independence test)
# ---------------------- #
form_risk <- as.formula(
  "Surv(Overall_survival, Vital_status) ~ Risk_score + Gender + Age + Stage + Smoking + M_stage + N_stage + T_stage")

# Fit Cox model including risk score and routine clinical variables
cox_clinic <- coxph(form_risk, data=cox_df, x=TRUE)
summary(cox_clinic)

pdf(file.path(result_dir, "Cox_clinic_Schoenfeld_residual.pdf"), width = 8, height = 6)
plot(ph)
dev.off()
print(ph)

pdf(file.path(result_dir, "FOREST.pdf"), width = 11, height = 13)
p <- ggforest(
  cox_clinic,
  data      = model.frame(cox_clinic),
  fontsize  = 1.3,                  # ← 글자 20% 키우기
  main      = "Hazard ratio",
  refLabel  = "reference",
  noDigits  = 2
)
print(p)
dev.off()

pdf(file.path(result_dir, "FOREST.pdf"), width = 9, height = 8)
print( ggforest(cox_clinic, data = model.frame(cox_clinic)) )
dev.off()

# Checking the proportional hazards (PH) assumption
ph <- cox.zph(cox_clinic)
pdf(file.path(result_dir, "PH_clinic.pdf"), width = 9, height = 8)
PH_clinic <-survminer::ggcoxzph(ph, font.axis=8,
                             font.tickslab = 6,
                             font.x = 8,
                             font.y = 8,
                             font.main = 8)
print(PH_clinic)
dev.off()

# Export tidy summary (HR, 95% CI, p-values)
write.csv(broom::tidy(cox_clinic, exponentiate=TRUE, conf.int=TRUE),
          file.path(result_dir, "Cox_risk_vs_clinical.csv"), row.names=FALSE)


# ---------------------- #
# 4) Time-dependent ROC (Overall_survival is in years 1, 3, 5)
# ---------------------- #
times_vec <- c(1,3,5)

# AUC for risk score at 1/3/5 years
roc_obj <- timeROC(T=cox_df$Overall_survival,
                   delta=as.integer(cox_df$Vital_status), # 1=event, 0=censored
                   marker=cox_df$Risk_score,
                   cause=1,
                   times=times_vec,
                   iid=FALSE)

# Prepare ROC curves for ggplot (stack FP/TP per time point)
curve_df <- do.call(rbind, lapply(seq_along(times_vec), function(i){
  data.frame(time_label=paste0(times_vec[i]," year"),
             FP=roc_obj$FP[,i], TP=roc_obj$TP[,i])
}))
auc_labs <- sprintf("%d year AUC = %.3f", times_vec, roc_obj$AUC)

# Plot overlaid ROC curves (1/3/5y) with AUC in legend
p_roc <- ggplot(curve_df, aes(x=FP, y=TP, color=time_label)) +
  geom_line(size=1.0) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey60") +
  scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,1)) +
  scale_color_manual(values=c("1 year"="#F8766D","3 year"="#00B0F6","5 year"="#00BA38"),
                     labels=auc_labs) +
  labs(title="ROC curves", x="False positive rate", y="True positive rate", color=NULL) +
  theme_minimal(base_size=12) +
  theme(legend.position=c(0.8,0.2),
        legend.text=element_text(size=13),
        legend.background=element_blank(),
        legend.key=element_blank(),
        panel.border=element_rect(color="black", fill=NA, linewidth=0.8))

pdf(file.path(result_dir, "ROC curve.pdf"), width=8, height=6)
print(p_roc)
dev.off()


# ---------------------- #
# 5) Risk score distribution & survival status (patient-ordered)
# ---------------------- #\
# Cox risk score with scale of hazard ratio
df <- cox_df[complete.cases(cox_df), ]
df$Risk_score <- predict(cox_multi, newdata = df, type = "risk") 

# Binary event indicator (1=event/death, 0=censored/alive)
df$event <- as.integer(df$Vital_status)

# Risk group split using median value (High vs Low)
cutoff <- median(df$Risk_score, na.rm=TRUE)
df$risk_group <- ifelse(df$Risk_score >= cutoff, "High","Low")

# Order patients by increasing risk score
df <- df[order(df$Risk_score), ]
df$idx <- seq_len(nrow(df))
df$Status <- ifelse(df$event==1, "Dead", "Alive")

# Risk score distribution (with dashed line at the cut point)
p_risk <- ggplot(df, aes(x=idx, y=Risk_score, color=risk_group)) +
  geom_point(size=1.6) +
  geom_vline(xintercept = sum(df$Risk_score < cutoff) + 0.5, linetype="dashed") +
  scale_color_manual(values=c("Low"="#00BFC4","High"="#F8766D")) +
  labs(x="Patients (Increasing risk score)", y="Risk score",
       title="Risk score distribution") +
  theme_minimal(base_size=12) +
  theme(legend.title=element_blank(),
        panel.border=element_rect(color="black", fill=NA, linewidth=0.8))

pdf(file.path(result_dir, "Distribution of risk score.pdf"), width=9, height=6)
print(p_risk)
dev.off()

# Survival time/status distribution in the same patient order
p_status <- ggplot(df, aes(x=idx, y=Overall_survival, color=Status)) +
  geom_point(size=1.6) +
  geom_vline(xintercept = sum(df$Risk_score < cutoff) + 0.5, linetype="dashed") +
  scale_color_manual(values=c("Dead"="#F8766D","Alive"="#00BFC4")) +
  labs(x="Patients (Increasing risk score)", y="Survival time (years)",
       title="Survival status by patient order") +
  theme_minimal(base_size=12) +
  theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.8))

pdf(file.path(result_dir, "Distribution of survival status.pdf"), width=9, height=6)
print(p_status)
dev.off()


# ---------------------- #
# 6) Heatmap of expression of signature ordered by risk score
# ---------------------- #
# Column order = patients sorted by risk score
rownames(df) <- df$case_id
ids <- df$case_id[order(df$Risk_score)]

# Column annotation: risk group and event status
annot <- data.frame(RiskGroup = df[ids, "risk_group"], row.names = ids)
annot$risk_group <- factor(annot$risk_group, levels = c("Low","High"))
annot <- annot[ids, , drop = FALSE]
ann_cols <- list(RiskGroup = c(Low = "#00BFC4", High = "#F8766D"))

# Gene-wise z-score for visualization
expr_z <- as.matrix(expr_wide[ , -1, drop = FALSE])
rownames(expr_z) <- expr_wide$case_id 
expr_z <- t(scale(expr_z))
idx <- match(ids, colnames(expr_z))
expr_z <- expr_z[, idx[!is.na(idx)], drop = FALSE]

# Draw heatmap
pdf(file.path(result_dir, "Heatmap_by_RiskScore.pdf"), width=7, height=5)
pheatmap(expr_z1,
         annotation_col= annot,
         annotation_colors = ann_cols,
         show_colnames=FALSE,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         fontsize_row=10)
dev.off()