library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(ggplot2)
library(Signac)
library(readxl)
library(dplyr)
library(GGally)
library(SummarizedExperiment)
library(DESeq2)
library(ggpubr)
library(tibble)
library(ggthemes)
library(data.table)
library(stringr)
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(UpSetR)
library(gplots)
library(RColorBrewer)
library(edgeR)
library(reshape2)
library(rWikiPathways)
library(msigdbr)
library(tidyverse)
library(fgsea)
library(grid)  
library(GSVA)
library(ggrepel)
library(impute)
library(WGCNA) 
library(tinyarray)
library(rbioapi)
library(org.Mm.eg.db)
library(ggvenn)

use_pinboard("onedrive")

mouse_list <- "4443|4436|4433|4419|4329|4535|4506"



# ---- Helper functions --------------------------------------------------------

#' Filter lowly expressed genes/miRNAs by CPM threshold
#'
#' @param count_mat  Raw count matrix (genes x samples)
#' @param cpm_cutoff Minimum CPM to count as expressed (default 0.5)
#' @return Logical vector; TRUE = keep
filter_by_cpm <- function(count_mat, cpm_cutoff = 0.5) {
  cpm_vals <- count_mat / colSums(count_mat) * 1e6
  rowSums(cpm_vals >= cpm_cutoff) > 0
}


#' Mean-center a log-CPM matrix by row means
#'
#' @param log_cpm  Log-CPM matrix (genes x samples)
#' @param row_means Named vector of row means; if NULL, computed from log_cpm
#' @return Row-mean-centered matrix
mean_center <- function(log_cpm, row_means = NULL) {
  if (is.null(row_means)) row_means <- rowMeans(log_cpm)
  sweep(log_cpm, 1, row_means, FUN = "-")
}


#' Build SVD state space from a centered expression matrix
#'
#' @param centered_mat Row-mean-centered matrix (genes x samples)
#' @return List with bV (right singular vectors, rownames = genes),
#'         bU (left singular vectors), bD (singular values)
build_svd_space <- function(centered_mat) {
  decomp <- svd(t(centered_mat))
  bV <- decomp$v
  rownames(bV) <- rownames(centered_mat)
  list(bV = bV, bU = decomp$u, bD = decomp$d)
}


#' Project samples onto a pre-built SVD state space
#'
#' @param centered_mat  Row-mean-centered matrix (genes x samples),
#'                      rows must be a subset of rownames(svd_space$bV)
#' @param svd_space     Output of build_svd_space()
#' @param n_pcs         Number of PCs to retain (default 9)
#' @return Data frame of projected coordinates (samples x PCs) + library_id
project_onto_svd <- function(centered_mat, svd_space, n_pcs = 9) {
  shared_genes <- intersect(rownames(centered_mat), rownames(svd_space$bV))
  proj <- t(centered_mat[shared_genes, ]) %*% svd_space$bV[shared_genes, seq_len(n_pcs)]
  proj_df <- as.data.frame(proj)
  proj_df$library_id <- rownames(proj_df)
  proj_df
}

#' Standard ggplot2 theme used throughout the analysis
theme_aml <- function(base_size = 25) {
  theme(
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.background   = element_rect(fill = "white", color = "white"),
    panel.border       = element_rect(fill = "transparent"),
    axis.title.x       = element_text(size = base_size),
    axis.title.y       = element_text(size = base_size),
    axis.text.x        = element_text(size = base_size),
    axis.text.y        = element_text(size = base_size),
    legend.text        = element_text(size = base_size - 5),
    legend.title       = element_text(size = base_size - 5)
  )
}


#' Convert timepoint strings (e.g. "T8p5") to numeric weeks
#'
#' @param tp Character vector of timepoint labels
#' @return Numeric vector
parse_timepoint <- function(tp) {
  as.numeric(gsub("p", ".", sub("^T", "", tp)))
}


# =============================================================================
# SECTION 1 — 2018 Cohort (Baseline State Space)
# =============================================================================

# ----  mRNA 2018 -----------------------------------------------------------

AML18       <- readRDS("~/Downloads/AML.mRNA.2018.rds")
meta_18     <- data.frame(AML18@colData@listData)
counts_18   <- AML18@assays@data@listData[["counts"]]

# Filter lowly expressed genes and compute log-CPM
keep_genes_18  <- filter_by_cpm(counts_18)
logcpm_18_mRNA <- cpm(counts_18, log = TRUE)
rowmean_18_mRNA <- rowMeans(logcpm_18_mRNA)

# Build SVD state space on expressed genes
mRNA_18_centered <- mean_center(logcpm_18_mRNA[keep_genes_18, ],
                                row_means = rowmean_18_mRNA[keep_genes_18])
svd_mRNA_18 <- build_svd_space(mRNA_18_centered)

# Project all 2018 samples and attach metadata
proj_mRNA_18 <- project_onto_svd(mRNA_18_centered, svd_mRNA_18)
proj_mRNA_18 <- merge(proj_mRNA_18, meta_18, by = "library_id")

CM18_m_proj <- proj_mRNA_18[, c("library_id", "V2", "mouse_id",
                                "tissue", "timepoint", "percent_ckit", "treatment")]
colnames(CM18_m_proj) <- c("library_id", "mRNA_PC2", "mouse_id",
                           "tissue", "sample_weeks", "percent_ckit", "treatment")
CM18_m_proj$sample_weeks <- parse_timepoint(CM18_m_proj$sample_weeks)
CM18_m_proj$treat_group  <- "CM18"

# Trajectory plot — CM treatment arm only
ggplot(CM18_m_proj[CM18_m_proj$treatment == "CM", ],
       aes(x = sample_weeks, y = mRNA_PC2, color = mouse_id)) +
  geom_point(size = 4) +
  geom_line(size = 1) +
  labs(x = "Time (weeks)", y = "mRNA PC2") +
  theme_aml()


# ----  miRNA 2018 ----------------------------------------------------------

miRNA_exp_18 <- read.table("~/Downloads/mirna.AML.miRNA.2018.tsv",
                           header = TRUE, row.names = 1)
colnames(miRNA_exp_18) <- gsub("_seqcluster", "", colnames(miRNA_exp_18))

all_meta   <- get_pin("metadata_mmu.csv")
mi_meta_18 <- all_meta[all_meta$cohort == "AML.miRNA.2018", ]
mi_18_exp  <- miRNA_exp_18[, mi_meta_18$library_id]

# Filter and log-CPM
keep_mi_18     <- filter_by_cpm(mi_18_exp)
logcpm_18_miRNA <- cpm(mi_18_exp[keep_mi_18, ], log = TRUE)
rowmean_18_miRNA <- rowMeans(logcpm_18_miRNA)

# Build SVD state space
mi_18_centered <- mean_center(logcpm_18_miRNA)
svd_miRNA_18   <- build_svd_space(mi_18_centered)

# Project and attach metadata
proj_miRNA_18 <- project_onto_svd(mi_18_centered, svd_miRNA_18)
proj_miRNA_18 <- merge(proj_miRNA_18, mi_meta_18, by = "library_id")

CM18_mi_proj <- proj_miRNA_18[, c("library_id", "V1", "mouse_id",
                                  "tissue", "timepoint", "percent_ckit", "treatment")]
colnames(CM18_mi_proj) <- c("library_id", "miRNA_PC1", "mouse_id",
                            "tissue", "sample_weeks", "percent_ckit", "treatment")
CM18_mi_proj$sample_weeks <- parse_timepoint(CM18_mi_proj$sample_weeks)
CM18_mi_proj$treat_group  <- "CM18"

# Trajectory plot — CM treatment arm only
ggplot(CM18_mi_proj[CM18_mi_proj$treatment == "CM", ],
       aes(x = sample_weeks, y = -miRNA_PC1, color = mouse_id)) +
  geom_point(size = 4) +
  geom_line(size = 1) +
  labs(x = "Time (weeks)", y = "miRNA PC1 (inverted)") +
  theme_aml()


# ---- 2D mRNA–miRNA scatter (2018 cohort) ---------------------------------

CM18_m_proj$combined_column  <- paste(CM18_m_proj$mouse_id,  CM18_m_proj$sample_weeks,  sep = "_")
CM18_mi_proj$combined_column <- paste(CM18_mi_proj$mouse_id, CM18_mi_proj$sample_weeks, sep = "_")
CM18_m_mi <- merge(CM18_m_proj, CM18_mi_proj, by = "combined_column")

ggplot(CM18_m_mi, aes(x = mRNA_PC2, y = miRNA_PC1, color = treatment.x)) +
  geom_point(size = 3) +
  labs(title = "2018 cohort — mRNA PC2 vs miRNA PC1",
       x = "mRNA PC2", y = "miRNA PC1") +
  theme_minimal() +
  theme_aml()



# =============================================================================
# SECTION 2 — Rx Cohort (Projected onto 2018 State Space)
# =============================================================================


# Shared metadata for both Rx mRNA and miRNA
df_mRNA_miRNA <- read_excel("/Users/ziachen/Documents/df_mRNA_miRNA.xlsx")

# ----  miRNA Rx ------------------------------------------------------------

# Subset to mice of interest
mi_data_ess <- df_mRNA_miRNA[grepl(mouse_list, df_mRNA_miRNA$mouse_id_y), ]
mi_Rx_meta  <- mi_data_ess[, c("mouse_id_y", "library_id_y", "tissue_y",
                               "time_rx", "cohort_x", "treatment_y", "percent_ckit_x")]
colnames(mi_Rx_meta) <- c("mouse_id", "library_id", "tissue",
                          "sample_weeks", "project", "treatment", "percent_ckit")

# Load Rx count matrix
mi_count_Rx <- read.table("~/Downloads/AML.miRNA.2022.Rx.mirna 1.tsv", header = TRUE)
rownames(mi_count_Rx) <- mi_count_Rx$miRNA
colnames(mi_count_Rx)  <- gsub("_seqcluster$", "", colnames(mi_count_Rx))
mi_count_Rx <- na.omit(mi_count_Rx[, mi_Rx_meta$library_id])

# Log-CPM, mean-center using 2018 row means, project onto 2018 SVD space
logcpm_Rx_miRNA  <- cpm(mi_count_Rx, log = TRUE)
shared_mi        <- intersect(rownames(logcpm_Rx_miRNA), rownames(svd_miRNA_18$bV))
Rx_mi_centered   <- sweep(logcpm_Rx_miRNA[shared_mi, ], 1,
                          rowmean_18_miRNA[shared_mi], FUN = "-")
proj_Rx_miRNA    <- project_onto_svd(Rx_mi_centered, svd_miRNA_18)
proj_Rx_miRNA    <- merge(proj_Rx_miRNA, mi_Rx_meta, by = "library_id")

Rx_mi_proj <- proj_Rx_miRNA[, c("library_id", "V1", "mouse_id",
                                "tissue", "sample_weeks", "percent_ckit", "treatment")]
colnames(Rx_mi_proj) <- c("library_id", "miRNA_PC1", "mouse_id",
                          "tissue", "sample_weeks", "percent_ckit", "treatment")
Rx_mi_proj$treat_group <- "Rx"

# Trajectory plot
ggplot(Rx_mi_proj, aes(x = sample_weeks, y = -miRNA_PC1, color = mouse_id)) +
  geom_point(size = 4) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(-5, 20)) +
  labs(title = "miRNA trajectories — Rx cohort",
       x = "Time (weeks)", y = "miRNA PC1") +
  theme_aml()


# ----  mRNA Rx -------------------------------------------------------------

CM_data_ess <- df_mRNA_miRNA[grepl(mouse_list, df_mRNA_miRNA$mouse_id_x), ]
mRNA_Rx_meta <- CM_data_ess[, c("mouse_id_x", "timepoint_x", "time_rx",
                                "timepoint_y", "treatment_y", "library_id_x",
                                "mRNA_PC2", "tissue_x", "percent_ckit_x")]

# Load Rx mRNA counts
AML_Rx       <- readRDS("~/Downloads/AML.mRNA.2022.Rx.se.rds")
all_counts_Rx <- AML_Rx@assays@data@listData[["counts"]]
counts_Rx     <- all_counts_Rx[, mRNA_Rx_meta$library_id_x]
rownames(counts_Rx) <- sub("\\..*", "", rownames(counts_Rx))   # strip version suffix
counts_Rx <- na.omit(counts_Rx)

# Log-CPM, mean-center using 2018 row means, project onto 2018 SVD space
logcpm_Rx_mRNA <- cpm(counts_Rx, log = TRUE)
shared_mRNA    <- intersect(rownames(logcpm_Rx_mRNA), rownames(svd_mRNA_18$bV))
Rx_mRNA_centered <- sweep(logcpm_Rx_mRNA[shared_mRNA, ], 1,
                          rowmean_18_mRNA[shared_mRNA], FUN = "-")
proj_Rx_mRNA   <- project_onto_svd(Rx_mRNA_centered, svd_mRNA_18)
proj_Rx_mRNA   <- merge(proj_Rx_mRNA, mRNA_Rx_meta,
                        by.x = "library_id", by.y = "library_id_x")

Rx_m_proj <- proj_Rx_mRNA[, c("library_id", "V2", "mouse_id_x",
                              "tissue_x", "time_rx", "percent_ckit_x", "treatment_y")]
colnames(Rx_m_proj) <- c("library_id", "mRNA_PC2", "mouse_id",
                         "tissue", "sample_weeks", "percent_ckit", "treatment")
Rx_m_proj$sample_weeks <- round(Rx_m_proj$sample_weeks)
Rx_m_proj$treat_group  <- "Rx"

# Trajectory plot with grand mean overlay
ggplot(Rx_m_proj, aes(x = sample_weeks, y = mRNA_PC2, color = mouse_id)) +
  geom_point(size = 4) +
  geom_line(size = 1) +
  stat_summary(aes(group = 1), fun = mean, geom = "line",
               color = "black", size = 2) +
  scale_x_continuous(limits = c(-3, 10)) +
  labs(title = "mRNA trajectories — Rx cohort",
       x = "Time (weeks)", y = "mRNA PC2") +
  theme_aml()


# ----  2D mRNA–miRNA scatter (Rx cohort) -----------------------------------

Rx_m_proj$combined_column  <- paste(Rx_m_proj$mouse_id,  Rx_m_proj$sample_weeks,  sep = "_")
Rx_mi_proj$combined_column <- paste(Rx_mi_proj$mouse_id, Rx_mi_proj$sample_weeks, sep = "_")
Rx_m_mi <- merge(Rx_m_proj, Rx_mi_proj, by = "combined_column")

ggplot(Rx_m_mi, aes(x = mRNA_PC2, y = miRNA_PC1, color = treatment.x)) +
  geom_point(size = 3) +
  labs(title = "Rx cohort — mRNA PC2 vs miRNA PC1",
       x = "mRNA PC2", y = "miRNA PC1") +
  theme_minimal() +
  theme_aml()


# =============================================================================
# SECTION 3 — combine 2D state-space
# =============================================================================



CM18_m_mi_result$treat_group.x <- ifelse(CM18_m_mi_result$treatment.x == "CM", "18_CM", "18_Ctrl")
m_miRNA_combined_all <- rbind(CM18_m_mi_result[CM18_m_mi_result$treat_group.x == "18_CM",],Rx_m_mi_result)

#write.csv(m_miRNA_combined_all, file = "m_miRNA_combined_all.csv", row.names = FALSE)

m_miRNA_combined_all <- m_miRNA_combined_all %>%
  filter(!(sample_weeks.y > 10))
m_miRNA_combined_all <- m_miRNA_combined_all %>%
  filter(!(sample_weeks.y < 0))

custom_colors <- c(
  "Rx" = "#e7298a", 
  "18_CM" = "#7570b3",              
  "18_Ctrl" = "#d95f02"        
)

#[m_miRNA_combined_all$treatment.x !='CM-ST2KO',]
ggplot(m_miRNA_combined_all, aes(x = mRNA_PC2, y = -miRNA_PC1, color = treat_group.x)) +
  geom_point(size = 3) +
  scale_color_manual(values = custom_colors) +  # <-- assign your colors here
  labs(title = "2D Dot Plot", x = "mRNA_PC2", y = "miRNA_PC1") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = "white"),
    panel.border = element_rect(fill = "transparent"),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25), 
    axis.text.x = element_text(size = 25),
    axis.text.y = element_text(size = 25),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  ) +
  geom_smooth(
    data = subset(m_miRNA_combined_all, treat_group.x == "18_CM"),
    aes(x = mRNA_PC2, y = -miRNA_PC1),
    method = "lm",
    se = FALSE,  # Set to TRUE if you want to include the confidence interval
    color = "#7570b3",  # Choose the color for the regression line
    size = 1.5       # Adjust the size of the line
  ) +
  geom_smooth(
    data = subset(m_miRNA_combined_all, treat_group.x == "st2"),
    aes(x = mRNA_PC2, y = -miRNA_PC1),
    method = "lm",
    se = FALSE,  # Set to TRUE if you want to include the confidence interval
    color = "#a6761d",  # Choose the color for the regression line
    size = 1.5       # Adjust the size of the line
  )





# =============================================================================
# SECTION 4 — Pathway Figures (Figure 3)
# =============================================================================


# Rx meta
CM_data_ess <- df_mRNA_miRNA[grepl(mouse_list, df_mRNA_miRNA$mouse_id_x),]
Rx_meta <- CM_data_ess[, c("mouse_id_x", "timepoint_x",'time_rx','timepoint_y','treatment_y','library_id_x',"mRNA_PC2",
                            "tissue_x","percent_ckit_x")]
Rx_meta <- Rx_meta[,c("library_id_x","mouse_id_x","time_rx","mRNA_PC2")]
colnames(Rx_meta) <- c("library_id","mouse_id","sample_weeks","mRNA_PC2")


# Rx expression
Rx_exp <- counts_Rx
Rx_exp_sel <- Rx_exp[,Rx_meta$library_id]

# path_list_msigDB
{
  curpath = msigdbr(species = "mouse", category = "H")
  curpath$gs_name <- gsub("^HALLMARK_", "", curpath$gs_name)
  path_list_msigDB = split(x = curpath$ensembl_gene, f = curpath$gs_name)
}

re <- gsva(gsvaParam(Rx_exp_sel, path_list_msigDB)) # select database from above
re_t <- as.data.frame(t(re))
re_t$COHP <- rownames(re_t)
Rx_meta$mouse_time <- paste(Rx_meta$mouse_id, Rx_meta$sample_weeks, sep = "_")
re_t <- merge(re_t,Rx_meta[,c("library_id","mouse_time")], by.x ="COHP", by.y = "library_id")

all_pathways <- unique(curpath$gs_name)
final_result_sel <- re_t[,all_pathways]


r2_mat <- cor(final_result_sel, use = "pairwise.complete.obs")

dev.off()
heatmap(
  r2_mat,
  symm = TRUE,
  col = colorRampPalette(c("dodgerblue3", "#FFFFFF", "firebrick1"))(100)
)


# ----  select pathway set from heatmap -----------------------------------

# select pathway set
{
  path_sets_new1 <- list(
    "KRAS_SIGNALING_DN",                 
    "MYOGENESIS",                       
    "EPITHELIAL_MESENCHYMAL_TRANSITION", 
    "APICAL_SURFACE",                   
    "ANGIOGENESIS",                      
    "COAGULATION",                      
    "COMPLEMENT",                        
    "ESTROGEN_RESPONSE_LATE",           
    "ESTROGEN_RESPONSE_EARLY",           
    "IL2_STAT5_SIGNALING",              
    "INTERFERON_ALPHA_RESPONSE",         
    "INTERFERON_GAMMA_RESPONSE",        
    "HEDGEHOG_SIGNALING",                
    "ALLOGRAFT_REJECTION",              
    "KRAS_SIGNALING_UP",                 
    "INFLAMMATORY_RESPONSE",            
    "TNFA_SIGNALING_VIA_NFKB",           
    "APICAL_JUNCTION",                  
    "IL6_JAK_STAT3_SIGNALING"
  )
  
  
  path_sets_new2 <- list(
    "HEME_METABOLISM",
    "PANCREAS_BETA_CELLS",
    "BILE_ACID_METABOLISM",
    "XENOBIOTIC_METABOLISM",
    "UV_RESPONSE_DN",
    "NOTCH_SIGNALING",
    "HYPOXIA",
    "APOPTOSIS",
    "TGF_BETA_SIGNALING",
    "ANDROGEN_RESPONSE",
    "REACTIVE_OXYGEN_SPECIES_PATHWAY",
    "P53_PATHWAY",
    "CHOLESTEROL_HOMEOSTASIS",
    "UV_RESPONSE_UP")
  
  path_sets_new3 <- list(
    "WNT_BETA_CATENIN_SIGNALING",
    "MITOTIC_SPINDLE",
    "PROTEIN_SECRETION",
    "PI3K_AKT_MTOR_SIGNALING",
    "FATTY_ACID_METABOLISM",
    "PEROXISOME",
    "ADIPOGENESIS",
    "GLYCOLYSIS",
    "MYC_TARGETS_V2",
    "E2F_TARGETS",
    "G2M_CHECKPOINT",
    "OXIDATIVE_PHOSPHORYLATION",
    "MTORC1_SIGNALING",
    "UNFOLDED_PROTEIN_RESPONSE",
    "DNA_REPAIR",
    "MYC_TARGETS_V1","SPERMATOGENESIS")
}

# plot group dynamic(Figure 3B)
cols_to_plot <- path_sets_new1
cols_to_plot <- unlist(cols_to_plot)

pathway_table <- merge(re_t,Rx_meta[,c("library_id", "sample_weeks", "mouse_id")], by.x ="COHP", by.y = "library_id")
pathway_table$sample_weeks <- round(pathway_table$sample_weeks)
pathway_table_result <- pathway_table %>%
  group_by(sample_weeks) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
long_df <- pathway_table_result %>%
  pivot_longer(
    cols = all_of(cols_to_plot),
    names_to = "pathway",
    values_to = "value"
  )

mean_df <- long_df %>%
  group_by(sample_weeks) %>%
  summarise(mean_value = mean(value, na.rm = TRUE))

ggplot() +
  # grey individual trajectories
  geom_line(
    data = long_df,
    aes(x = sample_weeks, y = value, group = pathway),
    color = "grey70",
    linewidth = 1,
    alpha = 0.7
  ) +
  
  # red average line
  geom_line(
    data = mean_df,
    aes(x = sample_weeks, y = mean_value),
    color = "red",
    linewidth = 2
  ) +
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  scale_x_continuous(limits = c(-3, 10)) +
  
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = "white"),
    panel.border = element_rect(fill="transparent"), 
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25), 
    axis.text.x = element_text(size = 25),
    axis.text.y = element_text(size = 25),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  ) +
  labs(
    title = "",
    x = "Weeks",
    y = " "
  )


# =============================================================================
# SECTION 5 — mRNA Figures (Figure 4A, C, E, G, 5A)
# =============================================================================

CM_data_ess <- df_mRNA_miRNA[grepl(mouse_list, df_mRNA_miRNA$mouse_id_x),]
data_ess <- CM_data_ess[, c("mouse_id_x", "timepoint_x",'time_rx','timepoint_y','treatment_y','library_id_x',"mRNA_PC2",
                            "tissue_x","percent_ckit_x")]

mRNA_ess <- data_ess[data_ess$treatment_y == 'CM',]
mRNA_ess$new_group <- ifelse(mRNA_ess$timepoint_y == "END",
                             "END",
                             paste0("TIME", round(mRNA_ess$time_rx)))
meta_Rx_m <- mRNA_ess[,c('library_id_x','new_group')]
meta_Rx_m$new_group <- gsub("TIME-(\\d+)", "TIME_N\\1", meta_Rx_m$new_group)
meta_Rx_m$new_group[meta_Rx_m$new_group == 'END'] <- 'TLast'


inner_mRNA <- intersect(rownames(all_counts_Rx),rownames(m_bV))
mRNA_count <- ceiling(all_counts_Rx[inner_mRNA,meta_Rx_m$library_id_x])

dds <- DESeqDataSetFromMatrix(mRNA_count,
                              meta_Rx_m,
                              design = ~new_group)

dds <- DESeq(dds)

vsd <- varianceStabilizingTransformation(dds)
library(genefilter)
wpn_vsd <- getVarianceStabilizedData(dds)
rv_wpn <- rowVars(wpn_vsd)
q75_wpn <- quantile( rowVars(wpn_vsd), .75)  # <= original
expr_normalized <- wpn_vsd[ rv_wpn > q75_wpn, ]



expr_normalized_df <- data.frame(expr_normalized) %>%
  mutate(
    Gene_id = row.names(expr_normalized)
  ) %>%
  pivot_longer(-Gene_id)

input_mat = t(expr_normalized)
allowWGCNAThreads() 
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

sft = pickSoftThreshold(
  input_mat,             # <= Input data
  blockSize = 30,
  powerVector = powers,
  verbose = 5
)

#plot the result
{
  par(mfrow = c(1,2));
  cex1 = 0.9;
  
  plot(sft$fitIndices[, 1],
       -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       xlab = "Soft Threshold (power)",
       ylab = "Scale Free Topology Model Fit, signed R^2",
       main = paste("Scale independence")
  )
  text(sft$fitIndices[, 1],
       -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       labels = powers, cex = cex1, col = "red"
  )
  abline(h = 0.90, col = "red")
  plot(sft$fitIndices[, 1],
       sft$fitIndices[, 5],
       xlab = "Soft Threshold (power)",
       ylab = "Mean Connectivity",
       type = "n",
       main = paste("Mean connectivity")
  )
  text(sft$fitIndices[, 1],
       sft$fitIndices[, 5],
       labels = powers,
       cex = cex1, col = "red")
}

picked_power = 14
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)

netwk <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)
cor <- temp_cor     # Return cor function to original namespace

# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

mRNA_module_df <- module_df
# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$library_id_x = row.names(MEs0)
MEs0 <- merge(MEs0,meta_Rx_m, by = 'library_id_x')
MEs0_f <- subset(MEs0, select = -library_id_x)
mME = MEs0_f %>%
  pivot_longer(-new_group) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME$new_group <- factor(mME$new_group, 
                        levels = c("TIME_N5", "TIME_N4", "TIME_N3",
                                   'TIME_N2','TIME_N1',"TIME0",
                                   "TIME2","TIME3","TIME4","TIME5",
                                   "TIME6","TIME7","TIME8","TIME9",
                                   "TIME10","TIME12","TIME14",'TLast'))
mME %>% ggplot(., aes(x=new_group, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-time Relationships", y = "Modules", fill="corr")



# ----  Figure 4A  -----------------------------------

inter_all_mir <- intersect(mRNA_module_df$gene_id,rownames(m_bV))

sub_m_bV <- as.data.frame(m_bV[inter_all_mir,])
sub_m_bV$gene_id <- rownames(sub_m_bV)
group_m_bV <- merge(sub_m_bV,mRNA_module_df,by = "gene_id")


group_m_bV$colors <- factor(
  group_m_bV$colors,
  levels = c("pink","yellow","blue","red","black","brown","grey","turquoise","green")
)

color_map_B <- c(
  "pink"      = "yellow",  # bright pink
  "yellow"    = "#F781BF",  # vivid yellow
  "blue"      = "#1f7a47",  # medium blue (different from deep/navy)
  "red"       = "#fdcdac",  # true bright red (not purple-leaning)
  "black"     = "#4D4D4D",  # neutral dark gray (prints better than brownish)
  "brown"     = "#00008B",  # clear brown (more saturated)
  "grey"      = "#999999",  # mid gray (distinct from "black")
  "turquoise" = "#BF00FF",  # shifted to green-teal (clearly different)
  "green"     = "#A65628"   # bright cyan-blue (avoids overlap with green above)
)


ggplot(group_m_bV, aes(x = colors, y = V2, fill = colors)) +
  geom_boxplot() + 
  labs(title = "", y = "Aerage Loading Value", x = "Co-expressed Group") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = "white"),
    panel.border = element_rect(fill="transparent"),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 12)) +
  scale_fill_manual(values = color_map_B) 


# ----  Figure 4C  -----------------------------------

colors_to_process <- unique(mRNA_module_df$colors)
proj_results_list <- list()
for (color in colors_to_process){
  colored_group <- mRNA_module_df[mRNA_module_df$colors == color, ]
  inter_miR <- intersect(rownames(Rx_cpm_rowmean),colored_group$gene_id)
  mi.treat <- t(Rx_cpm_rowmean[inter_miR,])  %*% m_bV[inter_miR,c(1:9)]
  
  mi.treat <- as.data.frame(mi.treat)
  mi.treat$library_id_x <- rownames(mi.treat)
  proj_result <- merge(mi.treat, mRNA_ess, by = 'library_id_x')
  
  proj_result$mouse_id_x <- factor(proj_result$mouse_id_x , levels = c("4443", "4436", "4433", "4419"
                                                                       , "4329", 
                                                                       "4535", "4506"))
  plot_data <- proj_result[!is.na(proj_result$mouse_id), ]
  plot_data$color <- color
  proj_results_list[[color]] <- plot_data
}
mR_final_combined_table <- do.call(rbind, proj_results_list)
mR_final_combined_table$mouse_color <- paste(mR_final_combined_table$mouse_id, mR_final_combined_table$color, sep = "_")

input_data <- mR_final_combined_table
input_data$time_rx <- round(input_data$time_rx)
input_data <- input_data[input_data$time_rx <= 10,]
result <- input_data %>%
  group_by(color, time_rx) %>%
  summarise(mean_V2 = mean(V2), .groups = 'drop')

result$color <- factor(
  result$color,
  levels = c("pink","yellow","blue","red","black","brown","grey","turquoise","green")
)

ggplot(result, aes(x = time_rx, y = mean_V2, color = color)) +
  geom_point(size = 4)+
  geom_line(size = 1) +
  scale_x_continuous(limits = c(-3, 10),breaks = seq(-3, 10, by = 1)) +
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = "white"),
    panel.border = element_rect(fill="transparent"), 
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25), 
    axis.text.x = element_text(size = 25),
    axis.text.y = element_text(size = 25),
    
    # Increase legend text size (if you have a legend)
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  )  +
  labs(title = "mRNA mouse trajectories", x = "Weeks", y = "pc2") +          
  scale_color_manual(values = color_map_B) 


# ----  Figure 4E  -----------------------------------

mR_input <- mR_final_combined_table
mR_input$time_rx <- round(mR_input$time_rx)
mR_input <- mR_input[mR_input$time_rx <= 10,]
mR_input <- mR_input[mR_input$time_rx >= -3,]
mR_result <- mR_input %>%
  group_by(color, time_rx) %>%
  summarise(mean_V2 = mean(V2), .groups = 'drop')

mR_result_percet <- mR_result %>%
  group_by(time_rx) %>%
  mutate(Percent = mean_V2 / sum(mean_V2) * 100)

mR_result_percet$color <- factor(
  mR_result_percet$color,
  levels = rev(c("pink","yellow","blue","red","black","brown","grey","turquoise","green"))
)


# Plot 100% stacked bar chart
ggplot(mR_result_percet, aes(x = factor(time_rx), y = Percent, fill = color)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_map_B) +
  labs(title = "100% Stacked Bar Chart", y = "Percentage (%)", x = "Category") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = "white"),
    panel.border = element_rect(fill="transparent"), 
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25), 
    axis.text.x = element_text(size = 25),
    axis.text.y = element_text(size = 25),
    
    # Increase legend text size (if you have a legend)
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  )


# ----  Figure 4G  -----------------------------------

input_data <- mR_final_combined_table
input_data <- input_data[input_data$time_rx <= 10,]
input_data <- input_data[input_data$time_rx >= -3,]
input_data$time_rx <- round(input_data$time_rx)
input_data$mouse_time <- paste(input_data$mouse_id, input_data$time_rx, sep = "_")

Rx_cpm <- cpm(counts_Rx, log = T)
inner_mRNA <- intersect(rownames(Rx_cpm),rownames(m_bV))

# minus 18 mean
Rx_cpm_rowmean <- sweep(Rx_cpm[inner_mRNA,], 1, rowMean_18_mRNA[inner_mRNA], FUN="-")

# proj 18 bv
bU.treat <- t(Rx_cpm_rowmean)  %*% m_bV[inner_mRNA,c(1:9)]

# create result table
bU.treat <- as.data.frame(bU.treat)
bU.treat$library_id <- rownames(bU.treat)
proj_result <- merge(bU.treat, data_ess, by.x = 'library_id',by.y = "library_id_x")

# rename for standard output
Rx_m_proj <- proj_result[,c("library_id","V2","mouse_id_x","tissue_x","time_rx","percent_ckit_x","treatment_y")]
colnames(Rx_m_proj) <- c("library_id","mRNA_PC2","mouse_id","tissue","sample_weeks","percent_ckit","treatment")

Rx_m_proj_sel <- Rx_m_proj[,c("mRNA_PC2","mouse_id","sample_weeks","percent_ckit")]


Rx_m_proj_sel <- Rx_m_proj_sel[Rx_m_proj_sel$sample_weeks <= 10,]
Rx_m_proj_sel <- Rx_m_proj_sel[Rx_m_proj_sel$sample_weeks >= -3,]
Rx_m_proj_sel$sample_weeks <- round(Rx_m_proj_sel$sample_weeks)
Rx_m_proj_sel$mouse_time <- paste(Rx_m_proj_sel$mouse_id, Rx_m_proj_sel$sample_weeks, sep = "_")



# color corelation
relation_result <- merge(input_data, Rx_m_proj_sel, by = "mouse_time")

model <- lm(V2 ~ percent_ckit, data = relation_result)
R2 <- summary(model)$r.squared


relation_result$color <- factor(
  relation_result$color,
  levels = c("pink","yellow","blue","red","black","brown","grey","turquoise","green")
)

relation_result_norm <- relation_result %>%
  group_by(color) %>%
  mutate(
    mean_V2_norm = (V2 - min(V2, na.rm = TRUE)) /
      (max(V2, na.rm = TRUE) - min(V2, na.rm = TRUE)),
  ) %>%
  ungroup()

ggplot(relation_result_norm[relation_result_norm$color == "blue"|
                              relation_result_norm$color == "yellow"|
                              relation_result_norm$color == "turquoise",],
       aes(x = mean_V2_norm,
           y = percent_ckit,
           color = color)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 3) +
  coord_cartesian(ylim = c(0, NA)) +
  scale_color_manual(values = color_map_B)+
  labs(
    x = "Normalized mean_V1 (0–1)",
    y = "Normalized mean_ckit (0–1)",
    color = "Color"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    axis.text.x = element_text(size = 25),
    axis.text.y = element_text(size = 25),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  )


# ----  Figure 5A  -----------------------------------

mRNA_ids <- mRNA_module_df$gene_id
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
# load("mart.RData")

mRNA_results <- getBM(
  attributes = c("ensembl_gene_id",
                 "chromosome_name",
                 "start_position",
                 "end_position",
                 "strand",
                 "gene_biotype",
                 "external_gene_name"),
  filters = "ensembl_gene_id",
  values = mRNA_ids,
  mart = mart
)

mRNA_locations_result <- merge(mRNA_results,mRNA_module_df, by.x = "ensembl_gene_id", by.y = "gene_id" )
mRNA_locations_result <- mRNA_locations_result[!duplicated(mRNA_locations_result), ]

color_percentages <- mRNA_locations_result %>%
  group_by(chromosome_name, colors) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(chromosome_name) %>%
  mutate(percentage = count / sum(count) * 100)

color_percentages$chromosome_name <- factor(color_percentages$chromosome_name, 
                                            levels = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,"X","Y"))

color_percentages <- color_percentages[
  color_percentages$chromosome_name %in%
    as.character(c(1:19, "X", "Y")),
]

ggplot(color_percentages, aes(x = chromosome_name, y = count, fill = colors)) +
  geom_bar(stat = "identity") +
  labs(title = "100% Stacked Bar Chart", y = "Count") +
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = "white"),
    panel.border = element_rect(fill="transparent"),  
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 12))+
  scale_fill_manual(values = color_map_B)





# =============================================================================
# SECTION 6 — mRNA Figures (Figure 4B, D, F, H, 5B)
# =============================================================================

mi_data_ess <- df_mRNA_miRNA[grepl(mouse_list, df_mRNA_miRNA$mouse_id_y),]
mi_ess <- mi_data_ess[, c("mouse_id_y", "timepoint_x",'time_rx','timepoint_y','treatment_y','library_id_y',"miRNA_PC1")]

mi_ess$new_group <- ifelse(mi_ess$timepoint_y == "END",
                           "END",
                           paste0("TIME", round(mi_ess$time_rx)))
meta_Rx_mi <- mi_ess[,c('library_id_y','new_group')]
meta_Rx_mi$new_group <- gsub("TIME-(\\d+)", "TIME_N\\1", meta_Rx_mi$new_group)
meta_Rx_mi$new_group[meta_Rx_mi$new_group == 'END'] <- 'TLast'

rownames(mi_count) <- mi_count$miRNA
names(mi_count) <- gsub("_seqcluster$", "", names(mi_count))
mi_Rx_count <- mi_count[,meta_Rx_mi$library_id_y]

dds <- DESeqDataSetFromMatrix(mi_Rx_count,
                              meta_Rx_mi,
                              design = ~new_group)

dds <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds)
wpn_vsd <- getVarianceStabilizedData(dds)
rv_wpn <- rowVars(wpn_vsd)
q75_wpn <- quantile( rowVars(wpn_vsd), .75)  # <= original
expr_normalized <- wpn_vsd[ rv_wpn > q75_wpn, ]

# partial projection
{
  int_miR <- intersect(rownames(mi_Rx_cpm_rowmean),rownames(expr_normalized))
  mi.treat <- t(mi_Rx_cpm_rowmean[int_miR,])  %*% mi_bV[int_miR,c(1:9)]
  mi.treat <- as.data.frame(mi.treat)
  mi.treat$library_id_y <- rownames(mi.treat)
  proj_result <- merge(mi.treat, mi_ess, by = 'library_id_y')
  
  proj_result$mouse_id <- factor(proj_result$mouse_id , levels = c("4443", "4436", "4433", "4419"
                                                                   , "4329", "4324", "4321", 
                                                                   "4535", "4506"))
  plot_data <- proj_result[!is.na(proj_result$mouse_id), ]
  plot_sub <- plot_data[plot_data$mouse_id =='4329',]
  ggplot(plot_data, aes(x = time_rx, y = -V1, color = mouse_id)) +
    geom_point(size = 4) + 
    geom_line(size = 1) +
    scale_x_continuous(limits = c(-5, 20)) +
    theme(
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = "white"),
      panel.border = element_rect(fill="transparent"), 
      axis.title.x = element_text(size = 25),
      axis.title.y = element_text(size = 25), 
      axis.text.x = element_text(size = 25),
      axis.text.y = element_text(size = 25),
      
      # Increase legend text size (if you have a legend)
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20)
    ) +
    scale_color_manual(values = c("4443" = "red", "4436" = "orange", "4433" = 'yellow', "4419" = "lightgreen"
                                  , "4329" = "green", "4324" = "yellowgreen", "4321" = "lightskyblue", 
                                  "4535" = "dodgerblue", "4506" = "mediumblue")) +
    labs(title = "miRNA mouse trajectories", x = "time", y = "pc1") 
}

expr_normalized_df <- data.frame(expr_normalized) %>%
  mutate(
    Gene_id = row.names(expr_normalized)
  ) %>%
  pivot_longer(-Gene_id)
input_mat = t(expr_normalized)
allowWGCNAThreads() 
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

sft = pickSoftThreshold(
  input_mat,             # <= Input data
  blockSize = 30,
  powerVector = powers,
  verbose = 5
)

#plot the result
{
  par(mfrow = c(1,2));
  cex1 = 0.9;
  
  plot(sft$fitIndices[, 1],
       -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       xlab = "Soft Threshold (power)",
       ylab = "Scale Free Topology Model Fit, signed R^2",
       main = paste("Scale independence")
  )
  text(sft$fitIndices[, 1],
       -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       labels = powers, cex = cex1, col = "red"
  )
  abline(h = 0.90, col = "red")
  plot(sft$fitIndices[, 1],
       sft$fitIndices[, 5],
       xlab = "Soft Threshold (power)",
       ylab = "Mean Connectivity",
       type = "n",
       main = paste("Mean connectivity")
  )
  text(sft$fitIndices[, 1],
       sft$fitIndices[, 5],
       labels = powers,
       cex = cex1, col = "red")
}

picked_power = 7
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)

netwk <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)
cor <- temp_cor     # Return cor function to original namespace

# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
  
)
miR_module_df <- module_df

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$library_id_y = row.names(MEs0)
MEs0 <- merge(MEs0,meta_Rx_mi, by = 'library_id_y')
MEs0_f <- subset(MEs0, select = -library_id_y)
mME = MEs0_f %>%
  pivot_longer(-new_group) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME$new_group <- factor(mME$new_group, 
                        levels = c("TIME_N5", "TIME_N4", "TIME_N3",
                                   'TIME_N2','TIME_N1',"TIME0",
                                   "TIME2","TIME3","TIME4","TIME5",
                                   "TIME6","TIME7","TIME8","TIME9",
                                   "TIME10","TIME12","TIME14",'TLast'))
mME %>% ggplot(., aes(x=new_group, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-time Relationships", y = "Modules", fill="corr")

# ----  Figure 4B  -----------------------------------

inter_all_mir <- intersect(miR_module_df$gene_id,rownames(mi_bV))

sub_mi_bV <- as.data.frame(mi_bV[inter_all_mir,])
sub_mi_bV$gene_id <- rownames(sub_mi_bV)
group_mi_bV <- merge(sub_mi_bV,miR_module_df,by = "gene_id")


group_mi_bV$colors <- factor(
  group_mi_bV$colors,
  levels = c("red","turquoise","grey","yellow","brown","green","blue")
  
)

color_map_A <- c(
  "red"       = "red",  # soft coral-red
  "turquoise" = "#17BECF",  # medium mint
  "grey"      = "#6A3D9A",  # cool lavender-grey
  "yellow"    = "#FFD700",  # warm light yellow
  "brown"     = "#FF7F00",  # light tan-brown
  "green"     = "#addd8e",  # soft sage green
  "blue"      = "blue"   # cornflower blue
)


ggplot(group_mi_bV, aes(x = colors, y = -V1, fill = colors)) +
  geom_boxplot() + 
  labs(title = "", y = "Aerage Loading Value", x = "Co-expressed Group") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = "white"),
    panel.border = element_rect(fill="transparent"),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 12)) +
  scale_fill_manual(values = color_map_A)


# ----  Figure 4D  -----------------------------------

colors_to_process <- unique(miR_module_df$colors)
proj_results_list <- list()
for (color in colors_to_process){
  colored_group <- miR_module_df[miR_module_df$colors == color, ]
  inter_miR <- intersect(rownames(mi_Rx_cpm_rowmean),colored_group$gene_id)
  mi.treat <- t(mi_Rx_cpm_rowmean[inter_miR,])  %*% mi_bV[inter_miR,c(1:9)]
  
  mi.treat <- as.data.frame(mi.treat)
  mi.treat$library_id_y <- rownames(mi.treat)
  proj_result <- merge(mi.treat, mi_ess, by = 'library_id_y')
  
  proj_result$mouse_id <- factor(proj_result$mouse_id , levels = c("4443", "4436", "4433", "4419"
                                                                   , "4329", "4324", "4321", 
                                                                   "4535", "4506"))
  plot_data <- proj_result[!is.na(proj_result$mouse_id), ]
  plot_data$color <- color
  proj_results_list[[color]] <- plot_data
}
final_combined_table <- do.call(rbind, proj_results_list)
final_combined_table$mouse_color <- paste(final_combined_table$mouse_id, final_combined_table$color, sep = "_")
miR_final_WGCNA <- final_combined_table

input_data <- final_combined_table
input_data$time_rx <- round(input_data$time_rx)
input_data <- input_data[input_data$time_rx <= 10,]
input_data <- input_data[input_data$time_rx >= -3,]
result <- input_data %>%
  group_by(color, time_rx) %>%
  summarise(mean_V1 = mean(V1), .groups = 'drop')

ggplot(result, aes(x = time_rx, y = -mean_V1, color = color)) +
  geom_point(size = 4)+
  geom_line(size = 1)  +
  scale_x_continuous(limits = c(-3, 10), breaks = seq(-3, 10, by = 1)) +
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = "white"),
    panel.border = element_rect(fill="transparent"), 
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25), 
    axis.text.x = element_text(size = 25),
    axis.text.y = element_text(size = 25),
    
    # Increase legend text size (if you have a legend)
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  )  +
  labs(title = "miRNA mouse trajectories", x = "time", y = "pc1")  +
  scale_color_manual(values = color_map_A, breaks = names(color_map_A))


# ----  Figure 4F  -----------------------------------

miR_input <- miR_final_WGCNA
miR_input$time_rx <- round(miR_input$time_rx)
miR_input <- miR_input[miR_input$time_rx <= 10,]
miR_input <- miR_input[miR_input$time_rx >= -3,]
miR_input <- miR_input %>%
  group_by(color, time_rx) %>%
  summarise(mean_V1 = mean(V1), .groups = 'drop')
miR_input$mean_V1 <- miR_input$mean_V1 - min(miR_input$mean_V1) # shift everything to positive


miR_result_percet <- miR_input %>%
  group_by(time_rx) %>%
  mutate(Percent = mean_V1 / sum(mean_V1) * 100)


miR_result_percet$color <- factor(
  miR_result_percet$color,
  levels = rev(c("red","turquoise","grey","yellow","brown","green","blue"))
  
)
# Plot 100% stacked bar chart
ggplot(miR_result_percet, aes(x = factor(time_rx), y = Percent, fill = color)) +
  geom_bar(stat = "identity") +
  labs(title = "100% Stacked Bar Chart", y = "Percentage (%)", x = "Category") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = "white"),
    panel.border = element_rect(fill="transparent"), 
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25), 
    axis.text.x = element_text(size = 25),
    axis.text.y = element_text(size = 25),
    
    # Increase legend text size (if you have a legend)
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  ) +
  scale_fill_manual(values = color_map_A)


# ----  Figure 4H  -----------------------------------
input_data <- final_combined_table

input_data <- input_data[input_data$time_rx <= 10,]
input_data <- input_data[input_data$time_rx >= -3,]
input_data$mouse_time <- paste(input_data$mouse_id, input_data$time_rx, sep = "_")
input_data$time_rx <- round(input_data$time_rx)

mi_data_ess <- df_mRNA_miRNA[grepl(mouse_list, df_mRNA_miRNA$mouse_id_y),]
mi_Rx_ess <- mi_data_ess[, c("mouse_id_y",'library_id_y',"tissue_y", "time_rx","cohort_x",'treatment_y',"percent_ckit_x")]
colnames(mi_Rx_ess) <- c("mouse_id","library_id","tissue","sample_weeks","project","treatment","percent_ckit")


# Rx expression matrix
rownames(mi_count) <- mi_count$miRNA
names(mi_count) <- gsub("_seqcluster$", "", names(mi_count))
mi_Rx_count <- mi_count[,mi_Rx_ess$library_id]

# cpm
mi_Rx_count <- na.omit(mi_Rx_count)
mi_Rx_cpm <- cpm(mi_Rx_count, log = T)
inner_miRNA <- intersect(rownames(mi_Rx_cpm),rownames(mi_bV))

# minus 18 mean
mi_Rx_cpm_rowmean <- sweep(mi_Rx_cpm[inner_miRNA,], 1, rowmean18_miRNA[inner_miRNA], FUN="-")

# proj 18 bv
mi.treat <- t(mi_Rx_cpm_rowmean)  %*% mi_bV[inner_miRNA,c(1:9)]

# create result table
mi.treat <- as.data.frame(mi.treat)
mi.treat$library_id <- rownames(mi.treat)
proj_result <- merge(mi.treat, mi_Rx_ess, by = 'library_id')

# rename for standard output
Rx_mi_proj <- proj_result[,c("library_id","V1","mouse_id","tissue","sample_weeks","percent_ckit","treatment")]
colnames(Rx_mi_proj) <- c("library_id","miRNA_PC1","mouse_id","tissue","sample_weeks","percent_ckit","treatment")

Rx_mi_proj_sel <- Rx_mi_proj[,c("miRNA_PC1","mouse_id","sample_weeks","percent_ckit")]
Rx_mi_proj_sel <- Rx_mi_proj_sel[Rx_mi_proj_sel$sample_weeks <= 10,]
Rx_mi_proj_sel <- Rx_mi_proj_sel[Rx_mi_proj_sel$sample_weeks >= -3,]
Rx_mi_proj_sel$mouse_time <- paste(Rx_mi_proj_sel$mouse_id, Rx_mi_proj_sel$sample_weeks, sep = "_")
Rx_mi_proj_sel$sample_weeks <- round(Rx_mi_proj_sel$sample_weeks)

input_data$merge_col <- paste(input_data$mouse_id, input_data$time_rx, sep = "_")
Rx_mi_proj_sel$merge_col <- paste(Rx_mi_proj_sel$mouse_id, Rx_mi_proj_sel$sample_weeks, sep = "_")
relation_result <- merge(
  input_data,
  Rx_mi_proj_sel,
  by = "merge_col"
)
relation_result_norm <- relation_result %>%
  group_by(color) %>%
  mutate(
    mean_V1_norm = (V1 - min(V1, na.rm = TRUE)) /
      (max(V1, na.rm = TRUE) - min(V1, na.rm = TRUE))
  ) %>%
  ungroup()

ggplot(relation_result_norm[relation_result_norm$color == "blue"|
                              relation_result_norm$color == "red",],
       aes(x = mean_V1_norm,
           y = percent_ckit,
           color = color))  +
  geom_smooth(method = "lm", se = FALSE, linewidth = 3) +
  scale_color_manual(
    values = c(
      blue = "blue",
      brown = "brown",
      green = "green",
      purple = "purple",
      red = "red",
      turquoise = "turquoise",
      yellow = "yellow"
    )
  ) +
  labs(
    x = "Normalized mean_V1 (0–1)",
    y = "Normalized mean_ckit (0–1)",
    color = "Color"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    axis.text.x = element_text(size = 25),
    axis.text.y = element_text(size = 25),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  )



# ----  Figure 5B  -----------------------------------

clean_miR <- sub("-[35]p$", "", miR_module_df$gene_id)
attributes <- c("mirbase_id", "chromosome_name", "start_position")
# Get the data
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
#load("mart.RData")

miRNA_locations <- getBM(attributes = attributes,
                         filters = "mirbase_id",
                         values = clean_miR,
                         mart = ensembl)
miRNA_locations$gene_id <- gsub("mir", "miR", miRNA_locations$mirbase_id)

module_df_clean <- miR_module_df
module_df_clean$gene_id <- sub("-[35]p$", "", miR_module_df$gene_id)
miRNA_locations_result <- merge(miRNA_locations,module_df_clean, by = "gene_id" )
miRNA_locations_result <- miRNA_locations_result[!duplicated(miRNA_locations_result), ]

color_percentages <- miRNA_locations_result %>%
  group_by(chromosome_name, colors) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(chromosome_name) %>%
  mutate(percentage = count / sum(count) * 100)

color_percentages$chromosome_name <- factor(color_percentages$chromosome_name, 
                                            levels = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,"X","Y"))

ggplot(color_percentages, aes(x = chromosome_name, y = count, fill = colors)) +
  geom_bar(stat = "identity") +
  labs(title = "100% Stacked Bar Chart", y = "Count") +
  scale_fill_manual(values = color_map_A) +
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = "white"),
    panel.border = element_rect(fill="transparent"),  
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 12))
































