
# THP-1 macrophage secretome (PXD071481) — WT vs SRGN KO, M0/M1/M2

PATH_PG   <- "proteinGroups_B.txt"
PATH_SDRF <- "sdrf_FIXED.tsv"
OUT_DIR   <- "results"

FDR_THRESHOLD <- 0.05
FC_THRESHOLD  <- 1.0
MIN_PEPTIDES  <- 2

# 1.1 PACKAGES

cran_pkgs <- c("data.table", "tidyverse", "pheatmap", "FactoMineR",
               "factoextra", "scales", "ggrepel", "patchwork", "RColorBrewer",
               "viridis", "BiocManager", "dendextend")
new_cran <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly = TRUE)]
if (length(new_cran) > 0) install.packages(new_cran, repos = "https://cloud.r-project.org")

bioc_pkgs <- c("limma", "DEqMS", "clusterProfiler", "org.Hs.eg.db",
               "enrichplot", "AnnotationDbi", "ComplexHeatmap", "circlize")
new_bioc <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly = TRUE)]
if (length(new_bioc) > 0) BiocManager::install(new_bioc, ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(limma)
  library(DEqMS)
  library(pheatmap)
  library(FactoMineR)
  library(factoextra)
  library(scales)
  library(ggrepel)
  library(patchwork)
  library(RColorBrewer)
  library(viridis)
  library(dendextend)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(AnnotationDbi)
  library(ComplexHeatmap)
  library(circlize)
})

select <- dplyr::select

dir.create(file.path(OUT_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "tables"),  showWarnings = FALSE, recursive = TRUE)


# 1.2 LOAD FILES 

pg <- fread(PATH_PG, showProgress = FALSE, na.strings = c("", "NA"))
message(sprintf("  proteinGroups: %d rows × %d columns", nrow(pg), ncol(pg)))

sdrf_design <- NULL

if (file.exists(PATH_SDRF)) {
  sdrf_raw <- fread(PATH_SDRF, showProgress = FALSE)
  message(sprintf("  SDRF: %d rows × %d columns", nrow(sdrf_raw), ncol(sdrf_raw)))

  has_fixed_cols <- all(c("sample", "genotype", "polstate") %in% names(sdrf_raw))

  if (has_fixed_cols) {
    sdrf_design <- sdrf_raw %>%
      as_tibble() %>%
      dplyr::mutate(condition = paste(genotype, polstate, sep = "_")) %>%
      dplyr::select(sample, genotype, polstate, condition, everything())
    message("  Design loaded from sdrf_FIXED.tsv (direct column read).")
  } else {
    sdrf_norm <- sdrf_raw %>%
      as_tibble() %>%
      rename_with(~ gsub("\\[|\\]", "", gsub(" ", "_", tolower(.))))

    geno_col  <- intersect(c("factor_value_genotype", "characteristics_genotype"), names(sdrf_norm))[1]
    phase_col <- intersect(c("factor_value_phase", "characteristics_phase"), names(sdrf_norm))[1]

    if (!is.na(geno_col) && !is.na(phase_col)) {
      sdrf_design <- sdrf_norm %>%
        dplyr::mutate(
          genotype  = .data[[geno_col]],
          polstate  = .data[[phase_col]],
          condition = paste(genotype, polstate, sep = "_")
        ) %>%
        dplyr::group_by(genotype, polstate) %>%
        dplyr::mutate(rep = dplyr::row_number()) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(sample = paste(genotype, polstate, rep, sep = "_")) %>%
        dplyr::select(sample, genotype, polstate, condition, everything())
      message("  Design parsed from legacy SDRF factor value columns.")
    } else {
      message("  WARNING: SDRF columns not recognized — falling back to LFQ name reconstruction.")
    }
  }
}

if (is.null(sdrf_design)) {
  lfq_cols_all <- grep("^LFQ intensity ", names(pg), value = TRUE)
  sample_names <- trimws(sub("^LFQ intensity ", "", lfq_cols_all))
  sdrf_design <- tibble(
    sample    = sample_names,
    genotype  = sub("_(M[012])_.*", "", sample_names) %>% sub("_(M[012])$", "", .),
    polstate  = regmatches(sample_names, regexpr("M[012]", sample_names)),
    condition = paste(
      sub("_(M[012])_.*", "", sample_names) %>% sub("_(M[012])$", "", .),
      regmatches(sample_names, regexpr("M[012]", sample_names)),
      sep = "_"
    ),
    replicate = seq_along(sample_names)
  )
  message("  Design reconstructed:")
  print(sdrf_design %>% dplyr::select(sample, genotype, polstate, condition))
}

lfq_cols_check <- trimws(sub("^LFQ intensity ", "",
                              grep("^LFQ intensity ", names(pg), value = TRUE)))
missing_samples <- setdiff(lfq_cols_check, sdrf_design$sample)
extra_samples   <- setdiff(sdrf_design$sample, lfq_cols_check)

if (length(missing_samples) > 0 || length(extra_samples) > 0) {
  stop(sprintf(
    "SDRF ↔ LFQ mismatch!\n  In LFQ but not SDRF: %s\n  In SDRF but not LFQ: %s",
    paste(missing_samples, collapse = ", "),
    paste(extra_samples,   collapse = ", ")
  ))
}
message(sprintf("  Sample alignment check: OK (%d samples)", nrow(sdrf_design)))

# 1.3 QC FILTERING

rev_col  <- intersect(c("Reverse", "Decoy"), names(pg))[1]
cont_col <- "Potential contaminant"
site_col <- "Only identified by site"

n_before <- nrow(pg)

pg_clean <- pg[
  (is.na(pg[[rev_col]])  | pg[[rev_col]]  != "+") &
  (is.na(pg[[cont_col]]) | pg[[cont_col]] != "+") &
  (is.na(pg[[site_col]]) | pg[[site_col]] != "+")
]
pg_clean <- pg_clean[pg_clean[["Razor + unique peptides"]] >= MIN_PEPTIDES]

message(sprintf("  Before QC: %d proteins", n_before))
message(sprintf("  After QC:  %d proteins (removed %d)", nrow(pg_clean), n_before - nrow(pg_clean)))

lfq_cols <- grep("^LFQ intensity ", names(pg_clean), value = TRUE)
sample_names <- trimws(sub("^LFQ intensity ", "", lfq_cols))

for (col in lfq_cols) set(pg_clean, j = col, value = as.double(pg_clean[[col]]))

lfq_mat <- as.matrix(pg_clean[, ..lfq_cols])
rownames(lfq_mat) <- pg_clean[["Protein IDs"]]
colnames(lfq_mat) <- sample_names

msms_cols <- grep("^MS/MS count ", names(pg_clean), value = TRUE)
msms_mat  <- as.matrix(pg_clean[, ..msms_cols])
rownames(msms_mat) <- pg_clean[["Protein IDs"]]
colnames(msms_mat) <- trimws(sub("^MS/MS count ", "", msms_cols))

pep_cols <- grep("^Razor \\+ unique peptides ", names(pg_clean), value = TRUE)
if (length(pep_cols) == 0) pep_cols <- grep("^Razor.unique peptides ", names(pg_clean), value = TRUE)
pep_mat <- as.matrix(pg_clean[, ..pep_cols])
rownames(pep_mat) <- pg_clean[["Protein IDs"]]
colnames(pep_mat) <- trimws(sub("^Razor \\+ unique peptides ", "", pep_cols))

message(sprintf("  LFQ matrix: %d proteins × %d samples", nrow(lfq_mat), ncol(lfq_mat)))

# 1.4 LOG2 / NORMALIZE / IMPUTE 

lfq_mat[lfq_mat == 0] <- NA
log2_mat <- log2(lfq_mat)

n_valid_total <- rowSums(!is.na(log2_mat))
all_na_prots  <- names(n_valid_total[n_valid_total == 0])
if (length(all_na_prots) > 0) {
  message(sprintf("  Removing %d proteins with zero valid values in all samples",
                  length(all_na_prots)))
  log2_mat <- log2_mat[!rownames(log2_mat) %in% all_na_prots, ]
  lfq_mat  <- lfq_mat[!rownames(lfq_mat) %in% all_na_prots, ]
} else {
  message("  No proteins with all-NA across all samples — OK.")
}

global_median  <- median(log2_mat, na.rm = TRUE)
sample_medians <- apply(log2_mat, 2, median, na.rm = TRUE)
norm_mat <- sweep(log2_mat, 2, sample_medians - global_median, FUN = "-")

message(sprintf("  Median range before normalization: %.2f – %.2f",
                min(sample_medians), max(sample_medians)))
message(sprintf("  Median range after normalization:  %.2f – %.2f",
                min(apply(norm_mat, 2, median, na.rm = TRUE)),
                max(apply(norm_mat, 2, median, na.rm = TRUE))))

set.seed(42)
conditions   <- sdrf_design$condition
unique_conds <- unique(conditions)

valid_per_cond <- sapply(unique_conds, function(cond) {
  cols <- which(sdrf_design$condition == cond)
  rowSums(!is.na(norm_mat[, cols, drop = FALSE]))
})
rownames(valid_per_cond) <- rownames(norm_mat)

n_cond_present <- rowSums(valid_per_cond > 0)
n_cond_total   <- length(unique_conds)
mec_prots      <- rownames(valid_per_cond)[n_cond_present < n_cond_total]
mar_prots      <- rownames(valid_per_cond)[n_cond_present == n_cond_total]

message(sprintf("  Missing value classification:"))
message(sprintf("    MEC proteins (missing in >=1 entire condition): %d", length(mec_prots)))
message(sprintf("    MAR/MNAR-scattered proteins (present in all conditions): %d", length(mar_prots)))

mec_mask <- matrix(FALSE, nrow = nrow(norm_mat), ncol = ncol(norm_mat),
                   dimnames = dimnames(norm_mat))
for (cond in unique_conds) {
  cols <- which(sdrf_design$condition == cond)
  prots_missing_in_cond <- rownames(valid_per_cond)[valid_per_cond[, cond] == 0]
  for (p in prots_missing_in_cond) {
    mec_mask[p, cols] <- TRUE
  }
}
mec_mask <- mec_mask & is.na(norm_mat)

n_mec_cells <- sum(mec_mask)
n_mar_cells <- sum(is.na(norm_mat)) - n_mec_cells
message(sprintf("    MEC missing cells: %d (%.1f%% of all missing)",
                n_mec_cells, n_mec_cells / sum(is.na(norm_mat)) * 100))
message(sprintf("    MAR/MNAR-scattered missing cells: %d (%.1f%% of all missing)",
                n_mar_cells, n_mar_cells / sum(is.na(norm_mat)) * 100))

impute_minprob <- function(mat, q = 0.01, width = 0.3) {
  apply(mat, 2, function(x) {
    missing_idx <- is.na(x)
    if (sum(!missing_idx) < 10) return(x)
    obs    <- x[!missing_idx]
    mu_low <- quantile(obs, probs = q, na.rm = TRUE)
    sd_low <- mad(obs, na.rm = TRUE) * width
    x[missing_idx] <- rnorm(sum(missing_idx), mean = mu_low, sd = sd_low)
    x
  })
}

imp_mat <- impute_minprob(norm_mat)
rownames(imp_mat) <- rownames(norm_mat)

if (n_mec_cells > 0) {
  all_observed <- norm_mat[!is.na(norm_mat)]
  obs_min      <- min(all_observed, na.rm = TRUE)
  obs_q01      <- quantile(all_observed, 0.01, na.rm = TRUE)
  obs_mad      <- mad(all_observed, na.rm = TRUE)
  floor_value  <- min(obs_min, obs_q01 - 2.0 * obs_mad) - 0.5
  imp_mat[mec_mask] <- floor_value
  message(sprintf("    MEC floor value: %.2f", floor_value))
}

message(sprintf("  Missing before imputation: %d (%.1f%%)",
                sum(is.na(norm_mat)), mean(is.na(norm_mat)) * 100))
message(sprintf("  Missing after imputation:  %d", sum(is.na(imp_mat))))

remaining_prots <- rownames(imp_mat)
pg_clean <- pg_clean[pg_clean[["Protein IDs"]] %in% remaining_prots]
msms_mat <- msms_mat[rownames(msms_mat) %in% remaining_prots, ]
pep_mat  <- pep_mat[rownames(pep_mat) %in% remaining_prots, ]
message(sprintf("  Proteins after all-NA removal: %d", nrow(imp_mat)))

sdrf_design <- sdrf_design[match(colnames(imp_mat), sdrf_design$sample), ]
stopifnot(all(colnames(imp_mat) == sdrf_design$sample))
message("  SDRF ↔ matrix column order: aligned and verified.")

write.csv(as.data.frame(imp_mat),
          file.path(OUT_DIR, "tables", "01_normalized_imputed_matrix.csv"))

#  1.5 COLOR PALETTES & ANNOTATION 

cond_pal <- c(KO_M0 = "#1B7837", KO_M1 = "#E31A1C", KO_M2 = "#1F78B4",
              WT_M0 = "#7FBF7B", WT_M1 = "#FB9A99", WT_M2 = "#87CEFA")

col_ann <- sdrf_design %>%
  dplyr::select(sample, genotype, polstate) %>%
  column_to_rownames("sample") %>%
  as.data.frame()

ann_colors <- list(
  genotype  = c(KO = "#2166AC", WT = "#D6604D"),
  polstate  = c(M0 = "#4DAF4A", M1 = "#FD8D3C", M2 = "#800026")
)

polstate_bright <- c(M0 = "#4DAF4A", M1 = "#E31A1C", M2 = "#1F78B4")

# 1.6 QC PLOTS 

# 1.6.1 Boxplots before/after normalization 
mat_to_long <- function(mat, label) {
  as.data.frame(mat) %>%
    rownames_to_column("protein") %>%
    pivot_longer(-protein, names_to = "sample", values_to = "log2_LFQ") %>%
    dplyr::mutate(dataset = label) %>%
    dplyr::left_join(sdrf_design %>%
                       dplyr::select(sample, condition, genotype, polstate),
                     by = "sample")
}

long_raw  <- mat_to_long(log2_mat, "Before normalization")
long_norm <- mat_to_long(norm_mat, "After normalization")
long_both <- bind_rows(long_raw, long_norm) %>%
  dplyr::mutate(dataset = factor(dataset,
                                 levels = c("Before normalization",
                                            "After normalization")))

p_box <- ggplot(long_both %>% dplyr::filter(!is.na(log2_LFQ)),
                aes(x = sample, y = log2_LFQ, fill = condition)) +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.3, linewidth = 0.4) +
  scale_fill_manual(values = cond_pal) +
  facet_wrap(~ dataset, nrow = 2, scales = "free_y") +
  labs(title    = "LFQ intensity distributions — before and after median normalization",
       subtitle = "Each box = one sample | Boxes should align after normalization",
       x = NULL, y = expression(log[2]*"(LFQ intensity)"), fill = "Condition") +
  theme_bw(base_size = 11) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 8),
        strip.background = element_rect(fill = "grey92"),
        legend.position  = "right")

svg(file.path(OUT_DIR, "figures", "01_boxplots_norm.svg"), width = 12, height = 8)
print(p_box)
dev.off()
message("  Saved: 01_boxplots_norm.svg")

# 1.6.2 Missing value heatmap
miss_mat <- ifelse(is.na(norm_mat), 0, 1)
row_miss_order <- order(rowSums(miss_mat == 0), decreasing = TRUE)

miss_col_order <- order(match(sdrf_design$polstate, c("M0", "M1", "M2")),
                         match(sdrf_design$genotype, c("KO", "WT")))
miss_mat_ordered <- miss_mat[row_miss_order, miss_col_order]
miss_col_ann <- col_ann[miss_col_order, , drop = FALSE]

svg(file.path(OUT_DIR, "figures", "02_missing_value_heatmap.svg"),
    width = 10, height = 12)
pheatmap::pheatmap(miss_mat_ordered,
         color             = c("lightcoral", "black"),
         legend            = TRUE,
         annotation_col    = miss_col_ann,
         annotation_colors = ann_colors,
         cluster_cols      = FALSE,
         cluster_rows      = FALSE,
         show_rownames     = FALSE,
         fontsize_col      = 8,
         main              = "Missing value pattern (before imputation)\nBlack = Observed | Light red = Missing",
         border_color      = NA)
dev.off()
message("  Saved: 02_missing_value_heatmap.svg")

# 1.6.4 CV/SD per condition 
cv_list <- lapply(unique(sdrf_design$condition), function(cond) {
  samps   <- sdrf_design$sample[sdrf_design$condition == cond]
  sub_mat <- norm_mat[, samps, drop = FALSE]
  cv <- apply(sub_mat, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) return(NA_real_)
    sd(x)
  })
  tibble(protein = rownames(norm_mat), cv = cv, condition = cond)
})
cv_df <- bind_rows(cv_list) %>% dplyr::filter(!is.na(cv))

p_cv <- ggplot(cv_df, aes(x = cv, fill = condition, colour = condition)) +
  geom_density(alpha = 0.25, linewidth = 0.8) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values   = cond_pal) +
  scale_colour_manual(values = cond_pal) +
  coord_cartesian(xlim = c(0, quantile(cv_df$cv, 0.99, na.rm = TRUE))) +
  labs(title    = expression("Replicate SD per condition ("*log[2]*" scale)"),
       subtitle = "Dashed line = SD 0.5 | Lower SD = better replicate reproducibility",
       x = expression("SD between replicates ("*log[2]*" LFQ)"), y = "Density",
       fill = "Condition", colour = "Condition") +
  theme_bw(base_size = 13) +
  theme(legend.position = "right")

svg(file.path(OUT_DIR, "figures", "04_cv_per_condition.svg"), width = 8, height = 5)
print(p_cv)
dev.off()

cv_summary <- cv_df %>%
  dplyr::group_by(condition) %>%
  dplyr::summarise(median_sd    = median(cv, na.rm = TRUE),
                   pct_below_05 = mean(cv < 0.5, na.rm = TRUE) * 100,
                   .groups = "drop")
message("  Replicate SD summary:"); print(cv_summary)
write.csv(cv_summary, file.path(OUT_DIR, "tables", "02_cv_summary.csv"),
          row.names = FALSE)
message("  Saved: 04_cv_per_condition.svg")

# 1.7 EXPLORATORY ANALYSIS
pca_mat <- imp_mat
message(sprintf("  PCA: %d proteins (imputed matrix) × %d samples",
                nrow(pca_mat), ncol(pca_mat)))

pca_res <- prcomp(t(pca_mat), center = TRUE, scale. = TRUE)
var_exp <- summary(pca_res)$importance[2, 1:4] * 100

pca_df <- as_tibble(pca_res$x[, 1:4], rownames = "sample") %>%
  dplyr::left_join(sdrf_design %>%
                     dplyr::select(sample, condition, genotype, polstate),
                   by = "sample")

# 1.7.1 PCA 
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2,
                            colour = polstate, shape = genotype, label = sample)) +
  geom_point(size = 4, alpha = 0.9, stroke = 1.2) +
  geom_text_repel(size = 2.2, max.overlaps = 20, show.legend = FALSE,
                  segment.colour = "grey60", segment.size = 0.3) +
  scale_colour_manual(values = polstate_bright,
                      labels = c(M0 = "M0 (unstimulated)",
                                 M1 = "M1 (pro-inflammatory)",
                                 M2 = "M2 (anti-inflammatory)")) +
  scale_shape_manual(values = c(KO = 16, WT = 17),
                     labels = c(KO = "SRGN KO", WT = "Wild type")) +
  tryCatch(
    stat_ellipse(aes(group = interaction(genotype, polstate)),
                 type = "t", level = 0.68,
                 linetype = "dashed", linewidth = 0.4, alpha = 0.4),
    error = function(e) geom_blank()
  ) +
  labs(title    = expression("PCA — "*log[2]*" LFQ intensities (imputed, all proteins, z-scored)"),
       subtitle = sprintf("%d proteins | Color = polarization | Shape = genotype", nrow(pca_mat)),
       x = sprintf("PC1 (%.1f%% variance)", var_exp[1]),
       y = sprintf("PC2 (%.1f%% variance)", var_exp[2]),
       colour = "Polarization", shape = "Genotype") +
  theme_bw(base_size = 13) +
  theme(legend.position = "right",
        legend.text = element_text(size = 9))

svg(file.path(OUT_DIR, "figures", "05_pca_pc12.svg"), width = 9, height = 6)
print(p_pca)
dev.off()
message("  Saved: 05_pca_pc12.svg")

# 1.7.2 Hierarchical clustering dendrogram 
hc <- hclust(dist(t(pca_mat), method = "euclidean"), method = "ward.D2")
pol_vec        <- setNames(sdrf_design$polstate, sdrf_design$sample)
pol_dend_pal   <- c(M0 = "#4DAF4A", M1 = "#E31A1C", M2 = "#1F78B4")
branch_colors  <- pol_dend_pal[pol_vec[hc$labels]]

svg(file.path(OUT_DIR, "figures", "06_hierarchical_clustering.svg"),
    width = 10, height = 5)
par(mar = c(8, 4, 3, 1))
dend <- as.dendrogram(hc)
labels_colors(dend) <- branch_colors[order.dendrogram(dend)]
plot(dend, main = "Hierarchical clustering — Ward's D2 on Euclidean distance",
     ylab = "Distance")
legend("topright",
       legend = c("M0 (unstimulated)", "M1 (pro-inflammatory)", "M2 (anti-inflammatory)"),
       fill = pol_dend_pal, title = "Polarization", cex = 0.8, bty = "n")
dev.off()
message("  Saved: 06_hierarchical_clustering.svg")

# 2. DIFFERENTIAL ABUNDANCE ANALYSIS

group <- factor(sdrf_design$condition,
                levels = c("KO_M0", "KO_M1", "KO_M2",
                           "WT_M0", "WT_M1", "WT_M2"))

design_mat <- model.matrix(~ 0 + group)
colnames(design_mat) <- levels(group)
rownames(design_mat) <- colnames(imp_mat)

contrast_mat <- makeContrasts(
  KO_M1_vs_M0 = KO_M1 - KO_M0,
  KO_M2_vs_M0 = KO_M2 - KO_M0,
  KO_M2_vs_M1 = KO_M2 - KO_M1,
  WT_M1_vs_M0 = WT_M1 - WT_M0,
  WT_M2_vs_M0 = WT_M2 - WT_M0,
  WT_M2_vs_M1 = WT_M2 - WT_M1,
  KO_vs_WT_M0 = KO_M0 - WT_M0,
  KO_vs_WT_M1 = KO_M1 - WT_M1,
  KO_vs_WT_M2 = KO_M2 - WT_M2,
  levels = design_mat
)

contrast_names <- colnames(contrast_mat)

fit  <- lmFit(imp_mat, design_mat)
fit2 <- contrasts.fit(fit, contrast_mat)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

limma_results <- lapply(contrast_names, function(ct) {
  topTable(fit2, coef = ct, number = Inf, adjust.method = "BH",
           sort.by = "none") %>%
    rownames_to_column("protein") %>%
    as_tibble() %>%
    dplyr::rename(log2FC = logFC, pval = P.Value, FDR = adj.P.Val) %>%
    dplyr::mutate(contrast    = ct,
                  significant = FDR < FDR_THRESHOLD & abs(log2FC) >= FC_THRESHOLD,
                  direction   = dplyr::case_when(
                    significant & log2FC > 0 ~ "Up",
                    significant & log2FC < 0 ~ "Down",
                    TRUE                      ~ "NS"
                  ))
})
names(limma_results) <- contrast_names

limma_summary <- purrr::map_dfr(contrast_names, function(ct) {
  r <- limma_results[[ct]]
  tibble(contrast = ct,
         n_tested = sum(!is.na(r$pval)),
         n_sig    = sum(r$significant, na.rm = TRUE),
         n_up     = sum(r$direction == "Up",   na.rm = TRUE),
         n_down   = sum(r$direction == "Down", na.rm = TRUE))
})
message("  limma summary:"); print(limma_summary)
write.csv(limma_summary, file.path(OUT_DIR, "tables", "03_limma_summary.csv"),
          row.names = FALSE)

for (ct in contrast_names) {
  write.csv(limma_results[[ct]],
            file.path(OUT_DIR, "tables", sprintf("04_limma_%s.csv", ct)),
            row.names = FALSE)
}

# 2.4 DEqMS 

global_pep <- as.integer(pg_clean[["Razor + unique peptides"]])
names(global_pep) <- pg_clean[["Protein IDs"]]
global_pep[is.na(global_pep) | global_pep < 1L] <- 1L

fit2$count <- global_pep[rownames(fit2$coefficients)]
fit2$count[is.na(fit2$count)] <- 1L

fit_deqms <- spectraCounteBayes(fit2)


deqms_results <- lapply(contrast_names, function(ct) {
  tt <- topTable(fit_deqms, coef = ct, number = Inf,
                 adjust.method = "BH", sort.by = "none") %>%
    rownames_to_column("protein") %>%
    as_tibble() %>%
    dplyr::rename(log2FC = logFC, pval = P.Value, FDR = adj.P.Val)

  if (!is.null(fit_deqms$sca.p.value) && ct %in% colnames(fit_deqms$sca.p.value)) {
    sca_pval <- fit_deqms$sca.p.value[tt$protein, ct]
    sca_fdr  <- p.adjust(sca_pval, method = "BH")
    tt <- tt %>% dplyr::mutate(sca.pval = sca_pval, sca.FDR = sca_fdr)
  } else if ("sca.P.Value" %in% names(tt) && "sca.adj.pval" %in% names(tt)) {
    tt <- tt %>% dplyr::rename(sca.pval = sca.P.Value, sca.FDR = sca.adj.pval)
  } else {
    tt <- tt %>% dplyr::mutate(sca.pval = pval, sca.FDR = FDR)
  }

  tt %>%
    dplyr::mutate(contrast    = ct,
                  significant = sca.FDR < FDR_THRESHOLD & abs(log2FC) >= FC_THRESHOLD,
                  direction   = dplyr::case_when(
                    significant & log2FC > 0 ~ "Up",
                    significant & log2FC < 0 ~ "Down",
                    TRUE                      ~ "NS"
                  ))
})
names(deqms_results) <- contrast_names

deqms_summary <- purrr::map_dfr(contrast_names, function(ct) {
  r <- deqms_results[[ct]]
  tibble(contrast = ct,
         n_tested = sum(!is.na(r$sca.pval)),
         n_sig    = sum(r$significant, na.rm = TRUE),
         n_up     = sum(r$direction == "Up",   na.rm = TRUE),
         n_down   = sum(r$direction == "Down", na.rm = TRUE))
})
message("  DEqMS summary:"); print(deqms_summary)
write.csv(deqms_summary, file.path(OUT_DIR, "tables", "05_deqms_summary.csv"),
          row.names = FALSE)

for (ct in contrast_names) {
  write.csv(deqms_results[[ct]],
            file.path(OUT_DIR, "tables", sprintf("06_deqms_%s.csv", ct)),
            row.names = FALSE)
}

# DEqMS variance plot
tryCatch({
  var_df <- tibble(
    count    = fit_deqms$count,
    sigma2   = fit_deqms$sigma^2,
    sca_var  = if (!is.null(fit_deqms$sca.prior)) fit_deqms$sca.prior else NA_real_
  ) %>% dplyr::filter(!is.na(count), count > 0, !is.na(sigma2), sigma2 > 0)

  p_var <- ggplot(var_df, aes(x = log2(count), y = log2(sigma2))) +
    geom_point(alpha = 0.25, size = 0.8, colour = "steelblue") +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE, colour = "#D6604D",
                linewidth = 1, alpha = 0.2, span = 0.75) +
    labs(title    = "DEqMS: protein variance vs spectral count",
         subtitle = "Loess fit (red) = protein-specific prior variance | Negative slope = DEqMS working",
         x = expression(log[2]*"(Razor + unique peptides)"),
         y = expression(log[2]*"(residual variance "*sigma^2*")")) +
    theme_bw(base_size = 12)

  svg(file.path(OUT_DIR, "figures", "07_deqms_variance_plot.svg"), width = 7, height = 5)
  print(p_var)
  dev.off()
  message("  Saved: 07_deqms_variance_plot.svg")
}, error = function(e) {
  message(sprintf("  WARNING: DEqMS variance plot failed — %s", conditionMessage(e)))
})

# 3. VISUALIZATION 

sig_pal <- c("Up"   = "#D6604D",
             "Down" = "#2166AC",
             "NS"   = "grey72")

gene_map <- pg_clean %>%
  as_tibble() %>%
  dplyr::transmute(
    protein = `Protein IDs`,
    gene    = dplyr::coalesce(dplyr::na_if(`Gene names`, ""), `Protein IDs`),
    gene    = sub(";.*", "", gene)
  )

add_gene <- function(df) dplyr::left_join(df, gene_map, by = "protein")

# 3.1 Volcano plots 

make_volcano <- function(res, contrast_label, top_n = 8) {
  res <- add_gene(res) %>%
    dplyr::filter(!is.na(sca.FDR)) %>%
    dplyr::mutate(neg_log10_FDR = -log10(sca.FDR + 1e-300))

  y_cap <- max(10, ceiling(quantile(res$neg_log10_FDR, 0.98, na.rm = TRUE)))
  res <- res %>%
    dplyr::mutate(neg_log10_FDR_capped = pmin(neg_log10_FDR, y_cap))

  top_prots <- res %>%
    dplyr::filter(significant) %>%
    dplyr::arrange(sca.FDR) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(protein)

  res <- res %>%
    dplyr::mutate(label = dplyr::if_else(protein %in% top_prots, gene, NA_character_))

  n_up   <- sum(res$direction == "Up",   na.rm = TRUE)
  n_down <- sum(res$direction == "Down", na.rm = TRUE)

  if (nrow(res) == 0) return(NULL)

  x_range <- max(abs(res$log2FC), na.rm = TRUE)

  ggplot(res, aes(x = log2FC, y = neg_log10_FDR_capped,
                  colour = direction, size = direction)) +
    geom_point(alpha = 0.65) +
    geom_text_repel(
      aes(label = label),
      size              = 2.8,
      max.overlaps      = 12,
      min.segment.length = Inf,
      force             = 8,
      force_pull        = 0.5,
      box.padding       = 1.0,
      point.padding     = 0.5,
      segment.colour    = NA,
      segment.size      = 0,
      show.legend       = FALSE,
      na.rm             = TRUE,
      direction         = "both",
      seed              = 42
    ) +
    geom_vline(xintercept = c(-FC_THRESHOLD, FC_THRESHOLD),
               linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_hline(yintercept = -log10(FDR_THRESHOLD),
               linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    scale_colour_manual(values = sig_pal) +
    scale_size_manual(values  = c(Up = 1.8, Down = 1.8, NS = 0.6),
                      guide = "none") +
    
    annotate("text",
             x = x_range * 0.9, y = y_cap * 0.92,
             label = sprintf("Up: %d", n_up),
             colour = "#D6604D", size = 3.5, hjust = 1, vjust = 1) +
    annotate("text",
             x = -x_range * 0.9, y = y_cap * 0.92,
             label = sprintf("Down: %d", n_down),
             colour = "#2166AC", size = 3.5, hjust = 0, vjust = 1) +
    coord_cartesian(ylim = c(0, y_cap * 1.08)) +
    labs(title    = sprintf("Volcano: %s", gsub("_", " ", contrast_label)),
         subtitle = sprintf("DEqMS | FDR < %.2f & |log2FC| >= %.1f | n sig = %d", FDR_THRESHOLD, FC_THRESHOLD, n_up + n_down),
         x = expression(log[2]*" Fold Change"),
         y = expression(-log[10](FDR)),
         colour = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top",
          plot.subtitle   = element_text(size = 9, colour = "grey40"))
}

for (ct in contrast_names) {
  p <- make_volcano(deqms_results[[ct]], ct)
  if (!is.null(p)) {
    svg(file.path(OUT_DIR, "figures", sprintf("08_volcano_%s.svg", ct)), width = 7, height = 6)
    print(p)
    dev.off()
  }
}
message(sprintf("  Saved: 08_volcano_<contrast>.svg for all %d contrasts",
                length(contrast_names)))

# 3.2 Significant proteins dotplot 

dot_df <- deqms_summary %>%
  tidyr::pivot_longer(c(n_up, n_down), names_to = "direction", values_to = "n") %>%
  dplyr::mutate(
    direction = dplyr::recode(direction, n_up = "Up", n_down = "Down"),
    n_signed  = dplyr::if_else(direction == "Down", -n, n),
    contrast  = factor(contrast, levels = rev(contrast_names))
  )

p_dot <- ggplot(dot_df, aes(x = n_signed, y = contrast, fill = direction)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = abs(n),
                hjust = dplyr::if_else(direction == "Up", -0.2, 1.2)),
            size = 3.5, colour = "grey20") +
  geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.5) +
  scale_fill_manual(values = c(Up = "#D6604D", Down = "#2166AC")) +
  scale_x_continuous(labels = function(x) abs(x),
                     expand = expansion(mult = 0.15)) +
  labs(title    = "Significant proteins per contrast (DEqMS)",
       subtitle = sprintf("FDR < %.2f & |log2FC| >= %.1f", FDR_THRESHOLD, FC_THRESHOLD),
       x = "Number of significant proteins",
       y = NULL, fill = "Direction") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top",
        panel.grid.major.y = element_blank())

svg(file.path(OUT_DIR, "figures", "09_significant_proteins_dotplot.svg"), width = 7, height = 6)
print(p_dot)
dev.off()
message("  Saved: 09_significant_proteins_dotplot.svg")

# 3.3 Z-score heatmap of TOP significant proteins

sig_proteins <- unique(unlist(lapply(deqms_results, function(r) {
  r$protein[r$significant]
})))
sig_proteins <- intersect(sig_proteins, rownames(imp_mat))
message(sprintf("  Union of significant proteins: %d", length(sig_proteins)))

top_n_zscore <- 50
all_deqms <- bind_rows(deqms_results) %>%
  dplyr::filter(protein %in% sig_proteins) %>%
  dplyr::group_by(protein) %>%
  dplyr::summarise(min_fdr = min(sca.FDR, na.rm = TRUE),
                   max_abs_fc = max(abs(log2FC), na.rm = TRUE),
                   .groups = "drop") %>%
  dplyr::arrange(min_fdr, desc(max_abs_fc))

top_sig_proteins <- head(all_deqms$protein, top_n_zscore)
message(sprintf("  Top %d significant proteins for z-score heatmap", length(top_sig_proteins)))

if (length(top_sig_proteins) >= 5) {
  sig_mat <- imp_mat[top_sig_proteins, , drop = FALSE]
  sig_z <- t(scale(t(sig_mat)))

  col_order <- order(sdrf_design$condition)
  sig_z <- sig_z[, col_order]
  col_ann_ordered <- col_ann[col_order, , drop = FALSE]

  sig_z[sig_z > 3]  <- 3
  sig_z[sig_z < -3] <- -3

  sig_gene_map <- gene_map %>%
    dplyr::filter(protein %in% top_sig_proteins) %>%
    dplyr::distinct(protein, .keep_all = TRUE)
  row_labels <- sig_gene_map$gene[match(rownames(sig_z), sig_gene_map$protein)]

  sig_z_t <- t(sig_z)
  row_ann_ordered <- data.frame(
    genotype = sdrf_design$genotype[col_order],
    polstate = sdrf_design$polstate[col_order],
    row.names = colnames(sig_z)
  )

  svg(file.path(OUT_DIR, "figures", "10_heatmap_significant_zscore.svg"),
      width = max(14, length(top_sig_proteins) * 0.3 + 4), height = 6)
  pheatmap::pheatmap(sig_z_t,
           color             = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
           breaks            = seq(-3, 3, length.out = 101),
           annotation_row    = row_ann_ordered,
           annotation_colors = ann_colors,
           cluster_rows      = FALSE,
           cluster_cols      = TRUE,
           show_colnames     = TRUE,
           labels_col        = row_labels,
           fontsize_col      = 6,
           fontsize_row      = 8,
           main              = sprintf("Z-score heatmap: top %d most significant proteins",
                                       length(top_sig_proteins)),
           border_color      = NA)
  dev.off()
  message("  Saved: 10_heatmap_significant_zscore.svg")
}

# 4. BIOLOGICAL INTERPRETATION - GSEA 

all_proteins <- rownames(imp_mat)
uniprot_ids  <- sub(";.*", "", all_proteins)

uniprot_to_entrez <- suppressMessages(
  AnnotationDbi::select(org.Hs.eg.db,
                        keys    = uniprot_ids,
                        columns = c("ENTREZID", "SYMBOL", "GENENAME"),
                        keytype = "UNIPROT")
)

uniprot_to_entrez <- uniprot_to_entrez %>%
  as_tibble() %>%
  dplyr::rename(uniprot = UNIPROT, entrez = ENTREZID,
                symbol = SYMBOL, genename = GENENAME) %>%
  dplyr::group_by(uniprot) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

annot_table <- tibble(protein = all_proteins, uniprot = uniprot_ids) %>%
  dplyr::left_join(uniprot_to_entrez, by = "uniprot") %>%
  dplyr::left_join(gene_map, by = "protein")

n_mapped <- sum(!is.na(annot_table$entrez))
message(sprintf("  Mapped %d / %d proteins to Entrez IDs (%.1f%%)",
                n_mapped, nrow(annot_table),
                n_mapped / nrow(annot_table) * 100))

write.csv(annot_table, file.path(OUT_DIR, "tables", "07_protein_annotation.csv"),
          row.names = FALSE)

background_entrez <- annot_table %>%
  dplyr::filter(!is.na(entrez)) %>%
  dplyr::pull(entrez) %>%
  unique()

# 4.2 GSEA PER CONTRAST 

build_gene_list <- function(deqms_res, annot) {
  df <- deqms_res %>%
    dplyr::left_join(annot %>% dplyr::select(protein, entrez), by = "protein") %>%
    dplyr::filter(!is.na(entrez), !is.na(sca.pval), !is.na(log2FC)) %>%
    dplyr::mutate(rank_stat = -log10(sca.pval) * sign(log2FC)) %>%
    dplyr::distinct(entrez, .keep_all = TRUE)

  geneList <- df$rank_stat
  names(geneList) <- as.character(df$entrez)
  geneList <- sort(geneList, decreasing = TRUE)
  geneList
}

gsea_results <- list()

for (ct in contrast_names) {
  message(sprintf("  Running GSEA for: %s", ct))
  geneList <- build_gene_list(deqms_results[[ct]], annot_table)

  if (length(geneList) < 10) {
    message(sprintf("    Too few ranked genes (%d) — GSEA skipped", length(geneList)))
    next
  }

  gsea_results[[ct]] <- list()

  tryCatch({
    gse_go <- gseGO(geneList     = geneList,
                    OrgDb        = org.Hs.eg.db,
                    ont          = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    verbose      = FALSE,
                    seed         = 42)
    gsea_results[[ct]]$GO_BP <- gse_go
    n_sig <- sum(gse_go@result$p.adjust < 0.05, na.rm = TRUE)
    message(sprintf("    GO-BP GSEA: %d significant terms", n_sig))
  }, error = function(e) {
    message(sprintf("    GO-BP GSEA failed: %s", e$message))
  })

  tryCatch({
    gse_kegg <- gseKEGG(geneList     = geneList,
                        organism     = "hsa",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        verbose      = FALSE,
                        seed         = 42)
    gsea_results[[ct]]$KEGG <- gse_kegg
    n_sig <- sum(gse_kegg@result$p.adjust < 0.05, na.rm = TRUE)
    message(sprintf("    KEGG GSEA: %d significant terms", n_sig))
  }, error = function(e) {
    message(sprintf("    KEGG GSEA failed: %s", e$message))
  })
}

# 4.3 GSEA PLOTS — merged GO+KEGG, disease-filtered, p.adjust scale

disease_blacklist <- c(
  "Chagas disease", "Leishmaniasis", "African trypanosomiasis",
  "Amoebiasis", "Malaria", "Tuberculosis", "Hepatitis B", "Hepatitis C",
  "Influenza A", "Measles", "Pertussis", "Legionellosis",
  "Staphylococcus aureus infection", "Rheumatoid arthritis",
  "Allograft rejection", "Graft-versus-host disease",
  "Inflammatory bowel disease", "Systemic lupus erythematosus",
  "Prion diseases", "Viral myocarditis", "Toxoplasmosis",
  "Pathways in cancer", "Transcriptional misregulation in cancer",
  "Alcoholic liver disease"
)

for (ct in names(gsea_results)) {
  for (db in c("GO_BP", "KEGG")) {
    gse_obj <- gsea_results[[ct]][[db]]
    if (is.null(gse_obj)) next
    safe_ct <- gsub("[^A-Za-z0-9_]", "_", ct)
    db_lower <- tolower(db)
    tryCatch({
      write.csv(gse_obj@result,
                file.path(OUT_DIR, "tables", sprintf("08_gsea_%s_%s_results.csv", safe_ct, db_lower)),
                row.names = FALSE)
    }, error = function(e) NULL)
  }
}

# Select top contrasts for dotplots 
contrast_scores <- deqms_summary %>%
  dplyr::mutate(
    n_gsea_go = sapply(contrast, function(ct) {
      obj <- gsea_results[[ct]]$GO_BP
      if (is.null(obj)) return(0)
      sum(obj@result$p.adjust < 0.05, na.rm = TRUE)
    }),
    n_gsea_kegg = sapply(contrast, function(ct) {
      obj <- gsea_results[[ct]]$KEGG
      if (is.null(obj)) return(0)
      sum(obj@result$p.adjust < 0.05, na.rm = TRUE)
    }),
    total_enrichment = n_gsea_go + n_gsea_kegg
  ) %>% dplyr::arrange(desc(n_sig + total_enrichment))

message("  Contrast ranking for GSEA plots:")
print(contrast_scores %>% dplyr::select(contrast, n_sig, n_gsea_go, n_gsea_kegg))

# 4.3 COMBINED GSEA DOTPLOT 

all_gsea_terms <- list()
for (ct in names(gsea_results)) {
  for (db in c("GO_BP", "KEGG")) {
    gse_obj <- gsea_results[[ct]][[db]]
    if (is.null(gse_obj)) next
    sig_res <- gse_obj@result[!is.na(gse_obj@result$p.adjust) &
                               gse_obj@result$p.adjust < 0.05, ]
    if (nrow(sig_res) == 0) next
    # Filter disease terms
    sig_res <- sig_res %>% dplyr::filter(!Description %in% disease_blacklist)
    if (nrow(sig_res) == 0) next
    db_label <- ifelse(db == "GO_BP", "GO", "KEGG")
    sig_res <- sig_res %>%
      dplyr::mutate(source = db_label, contrast = ct) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = 15)
    all_gsea_terms[[paste(ct, db, sep = "_")]] <- sig_res
  }
}

if (length(all_gsea_terms) > 0) {
  combined_gsea_df <- bind_rows(all_gsea_terms) %>%
    dplyr::mutate(
      geno_label = dplyr::case_when(
        grepl("^KO_vs_WT", contrast) ~ "KO vs WT",
        grepl("^KO_", contrast)       ~ "KO",
        grepl("^WT_", contrast)       ~ "WT",
        TRUE ~ contrast
      ),
      pol_label = dplyr::case_when(
        grepl("M1_vs_M0$", contrast) ~ "M1 vs M0",
        grepl("M2_vs_M0$", contrast) ~ "M2 vs M0",
        grepl("M2_vs_M1$", contrast) ~ "M2 vs M1",
        grepl("M0$", contrast)       ~ "M0",
        grepl("M1$", contrast)       ~ "M1",
        grepl("M2$", contrast)       ~ "M2",
        TRUE ~ contrast
      ),
      contrast_label = paste0(gsub("_", " ", contrast), " [", source, "]"),
      term_id = paste0(Description, " [", source, "]")
    )

  term_freq <- combined_gsea_df %>%
    dplyr::count(term_id, sort = TRUE) %>%
    dplyr::rename(n_contrasts = n)


  top_terms <- term_freq %>%
    dplyr::arrange(dplyr::desc(n_contrasts)) %>%
    dplyr::slice_head(n = 30)

  plot_df <- combined_gsea_df %>%
    dplyr::filter(term_id %in% top_terms$term_id) %>%
    dplyr::left_join(term_freq, by = "term_id") %>%
    dplyr::mutate(
      term_id = factor(term_id, levels = rev(top_terms$term_id)),
      contrast_label = factor(contrast_label,
                              levels = sort(unique(combined_gsea_df$contrast_label)))
    )

  plot_df <- plot_df %>% dplyr::filter(contrast_label %in% levels(contrast_label))

  if (nrow(plot_df) > 0) {
    tryCatch({
      p_combined <- ggplot(plot_df, aes(x = contrast_label, y = term_id,
                                         colour = p.adjust, size = setSize)) +
        geom_point() +
        scale_colour_gradient2(low = "#D6604D", mid = "#FDC086", high = "#2166AC",
                               midpoint = 0.025,
                               name = "p.adjust", limits = c(0, 0.05),
                               oob = scales::squish) +
        scale_size_continuous(name = "Set size", range = c(3, 8)) +
        labs(title = "GSEA: all significant terms across contrasts",
             subtitle = "Merged GO-BP + KEGG | Disease-artifact terms removed | Size = gene set size",
             x = NULL, y = NULL) +
        theme_bw(base_size = 10) +
        theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 8),
              axis.text.y  = element_text(size = 7),
              plot.title   = element_text(size = 11),
              plot.subtitle = element_text(size = 8, colour = "grey40"),
              legend.position = "right",
              panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3))

      svg(file.path(OUT_DIR, "figures", "11_gsea_merged_all.svg"),
          width = 10, height = max(6, nrow(top_terms) * 0.32 + 2))
      print(p_combined)
      dev.off()
      message(sprintf("  Saved: 11_gsea_merged_all.svg (%d terms × %d contrasts)",
                      nrow(top_terms), length(unique(plot_df$contrast_label))))
    }, error = function(e) {
      message(sprintf("  Combined GSEA dotplot failed: %s", e$message))
    })
  }
} else {
  message("  No significant non-disease GSEA terms across any contrast — combined dotplot skipped.")
}


# 4.5 PER-PATHWAY PROTEIN EXPRESSION HEATMAPS 

specific_kegg_pathways <- c(
  "Cytokine-cytokine receptor interaction",
  "Intestinal immune network for IgA production",
  "Complement and coagulation cascades",
  "NOD-like receptor signaling pathway",
  "IL-17 signaling pathway"
)

priority_go_pathways <- c(
  "response to external biotic stimulus",
  "response to other organism",
  "response to biotic stimulus",
  "defense response to other organism",
  "defense response to symbiont"
)


all_pathway_ids <- list()
for (ct in names(gsea_results)) {
  for (db in c("GO_BP", "KEGG")) {
    gse_obj <- gsea_results[[ct]][[db]]
    if (is.null(gse_obj)) next
    sig_res <- gse_obj@result[!is.na(gse_obj@result$p.adjust) &
                               gse_obj@result$p.adjust < 0.05, ]
    if (nrow(sig_res) == 0) next
    sig_res <- sig_res %>% dplyr::filter(!Description %in% disease_blacklist)
    if (nrow(sig_res) == 0) next
    for (i in seq_len(nrow(sig_res))) {
      tid <- sig_res$ID[i]
      tdesc <- sig_res$Description[i]
      if (is.null(all_pathway_ids[[tid]])) {
        all_pathway_ids[[tid]] <- tdesc
      }
    }
  }
}


key_paths <- list()
for (tdesc in specific_kegg_pathways) {
  matching <- names(all_pathway_ids)[all_pathway_ids == tdesc]
  if (length(matching) > 0) {
    key_paths[[matching[1]]] <- tdesc
  } else {
    message(sprintf("    Pathway not found in GSEA results: %s — skipped.", tdesc))
  }
}
for (tdesc in priority_go_pathways) {
  matching <- names(all_pathway_ids)[all_pathway_ids == tdesc]
  if (length(matching) > 0) {
    key_paths[[matching[1]]] <- tdesc
  }
}

message(sprintf("  Selected %d pathways for protein heatmaps:", length(key_paths)))
for (tid in names(key_paths)) message(sprintf("    %s: %s", tid, key_paths[[tid]]))

for (tid in names(key_paths)) {
  tdesc <- key_paths[[tid]]
  member_entrez <- NULL

  for (ct in names(gsea_results)) {
    for (db in c("GO_BP", "KEGG")) {
      gse_obj <- gsea_results[[ct]][[db]]
      if (is.null(gse_obj)) next
      match_row <- gse_obj@result[gse_obj@result$ID == tid, ]
      if (nrow(match_row) > 0 && !is.na(match_row$p.adjust[1]) && match_row$p.adjust[1] < 0.05) {
        core_genes <- strsplit(as.character(match_row$core_enrichment[1]), "/")[[1]]
        member_entrez <- c(member_entrez, core_genes)
      }
    }
  }
  member_entrez <- unique(member_entrez)

  if (length(member_entrez) < 3) {
    message(sprintf("    Too few core enrichment genes (%d) for %s — skipped.",
                    length(member_entrez), tdesc))
    next
  }

  member_prots <- annot_table %>%
    dplyr::filter(as.character(entrez) %in% member_entrez, !is.na(entrez)) %>%
    dplyr::distinct(protein, .keep_all = TRUE) %>%
    dplyr::filter(protein %in% rownames(imp_mat))

  if (nrow(member_prots) < 3) next

  prot_ids   <- member_prots$protein
  prot_genes <- member_prots$gene

  path_mat <- imp_mat[prot_ids, , drop = FALSE]
  path_z   <- t(scale(t(path_mat)))
  path_z[path_z > 3]  <- 3
  path_z[path_z < -3] <- -3

  col_ord <- order(sdrf_design$condition)
  path_z <- path_z[, col_ord]
  col_ann_ord <- col_ann[col_ord, , drop = FALSE]

  safe_desc <- gsub("[^A-Za-z0-9_]", "_", substr(tdesc, 1, 50))

  svg(file.path(OUT_DIR, "figures",
                sprintf("13_protein_heatmap_%s.svg", safe_desc)),
      width = 10, height = max(5, length(prot_ids) * 0.25 + 3))
  pheatmap::pheatmap(path_z,
           color             = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
           breaks            = seq(-3, 3, length.out = 101),
           annotation_col    = col_ann_ord,
           annotation_colors = ann_colors,
           cluster_rows      = TRUE,
           cluster_cols      = FALSE,
           show_rownames     = TRUE,
           labels_row        = prot_genes,
           fontsize_row      = 7,
           fontsize_col      = 8,
           main              = sprintf("Protein expression: %s (core enrichment genes)", tdesc),
           border_color      = NA)
  dev.off()
  message(sprintf("  Saved: 13_protein_heatmap_%s.svg (%d proteins)", safe_desc, length(prot_ids)))
}


# 4.7 FUNCTIONAL-GROUP HEATMAPS 

# Collect union of all significant proteins across all contrasts
sig_prots_all <- unique(unlist(lapply(deqms_results, function(r) {
  r$protein[r$significant]
})))
sig_prots_all <- intersect(sig_prots_all, rownames(imp_mat))
message(sprintf("  Union of significant proteins: %d", length(sig_prots_all)))

# Map sig proteins to Entrez for ORA
sig_entrez <- annot_table %>%
  dplyr::filter(protein %in% sig_prots_all, !is.na(entrez)) %>%
  dplyr::distinct(entrez, .keep_all = TRUE)

if (nrow(sig_entrez) < 10) {
  message("  Too few annotated significant proteins for functional grouping — skipping.")
} else {
  tryCatch({
    ego_sig <- enrichGO(gene         = as.character(sig_entrez$entrez),
                        universe     = as.character(background_entrez),
                        OrgDb        = org.Hs.eg.db,
                        ont          = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2,
                        readable     = TRUE)

    if (is.null(ego_sig) || nrow(ego_sig@result[ego_sig@result$p.adjust < 0.05, ]) == 0) {
      message("  No significant GO-BP terms from ORA — trying relaxed cutoff (p < 0.1)...")
      ego_sig <- enrichGO(gene         = as.character(sig_entrez$entrez),
                          universe     = as.character(background_entrez),
                          OrgDb        = org.Hs.eg.db,
                          ont          = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.1,
                          qvalueCutoff = 0.3,
                          readable     = TRUE)
    }

    sig_go_res <- ego_sig@result[!is.na(ego_sig@result$p.adjust) & ego_sig@result$p.adjust < 0.1, ]

    if (nrow(sig_go_res) == 0) {
      message("  Still no GO-BP terms — using manual functional groups from GSEA results.")
      manual_groups <- list()
      for (ct in names(gsea_results)) {
        for (db in c("GO_BP", "KEGG")) {
          gse_obj <- gsea_results[[ct]][[db]]
          if (is.null(gse_obj)) next
          sig_r <- gse_obj@result[!is.na(gse_obj@result$p.adjust) & gse_obj@result$p.adjust < 0.05, ]
          if (nrow(sig_r) == 0) next
          for (i in seq_len(min(nrow(sig_r), 10))) {
            tid   <- sig_r$ID[i]
            tdesc <- sig_r$Description[i]
            core  <- strsplit(as.character(sig_r$core_enrichment[i]), "/")[[1]]
            if (is.null(manual_groups[[tid]])) {
              manual_groups[[tid]] <- list(desc = tdesc, entrez = core)
            } else {
              manual_groups[[tid]]$entrez <- union(manual_groups[[tid]]$entrez, core)
            }
          }
        }
      }
      sig_go_res <- data.frame(
        ID          = sapply(manual_groups, function(x) x$desc),
        Description = sapply(manual_groups, function(x) x$desc),
        GeneRatio   = NA_character_,
        BgRatio     = NA_character_,
        pvalue      = 0.01,
        p.adjust    = 0.05,
        qvalue      = 0.05,
        geneID      = sapply(manual_groups, function(x) paste(x$entrez, collapse = "/")),
        Count        = sapply(manual_groups, function(x) length(x$entrez)),
        stringsAsFactors = FALSE
      )
    }

    message(sprintf("  Found %d GO-BP functional groups", nrow(sig_go_res)))

    
    prot_to_group <- list()

    
    sig_go_res <- sig_go_res %>% dplyr::arrange(p.adjust)

    for (i in seq_len(nrow(sig_go_res))) {
      term_genes <- strsplit(as.character(sig_go_res$geneID[i]), "/")[[1]]
      term_desc  <- sig_go_res$Description[i]
      term_id    <- sig_go_res$ID[i]

      term_prots <- annot_table %>%
        dplyr::filter(symbol %in% term_genes | as.character(entrez) %in% term_genes) %>%
        dplyr::filter(protein %in% sig_prots_all) %>%
        dplyr::distinct(protein, .keep_all = TRUE) %>%
        dplyr::pull(protein)

      for (p in term_prots) {
        if (is.null(prot_to_group[[p]])) {
          prot_to_group[[p]] <- list(group_id = term_id, group_desc = term_desc)
        }
      }
    }

    # Proteins not assigned to any group → "Other"
    unassigned <- setdiff(sig_prots_all, names(prot_to_group))
    for (p in unassigned) {
      prot_to_group[[p]] <- list(group_id = "OTHER", group_desc = "Other")
    }

    # Build assignment table
    group_assign <- purrr::map_dfr(names(prot_to_group), function(p) {
      tibble(protein = p,
             group_id = prot_to_group[[p]]$group_id,
             group_desc = prot_to_group[[p]]$group_desc)
    }) %>%
      dplyr::left_join(gene_map, by = "protein") %>%
      dplyr::distinct(protein, .keep_all = TRUE)

    # Count proteins per group, keep groups with >= 3 proteins
    group_sizes <- group_assign %>%
      dplyr::count(group_id, group_desc, sort = TRUE, name = "n_prots")

    # Merge small groups (< 3 proteins) into "Other"
    small_groups <- group_sizes$group_id[group_sizes$n_prots < 3]
    group_assign <- group_assign %>%
      dplyr::mutate(
        group_id   = dplyr::if_else(group_id %in% small_groups, "OTHER", group_id),
        group_desc = dplyr::if_else(group_id %in% small_groups, "Other", group_desc)
      )

    # Re-count
    group_sizes <- group_assign %>%
      dplyr::count(group_id, group_desc, sort = TRUE, name = "n_prots")

   
    group_sizes <- group_sizes %>%
      dplyr::mutate(short_desc = ifelse(nchar(group_desc) > 40,
                                        paste0(substr(group_desc, 1, 37), "..."),
                                        group_desc))

    message("  Functional groups:")
    for (i in seq_len(nrow(group_sizes))) {
      message(sprintf("    %s (%d proteins): %s",
                      group_sizes$group_id[i], group_sizes$n_prots[i],
                      group_sizes$short_desc[i]))
    }

    
    group_order <- group_sizes %>%
      dplyr::arrange(dplyr::desc(n_prots)) %>%
      dplyr::filter(group_id != "OTHER") %>%
      dplyr::pull(group_id)
    if ("OTHER" %in% group_sizes$group_id) group_order <- c(group_order, "OTHER")

    
    short_desc_lookup <- setNames(group_sizes$short_desc, group_sizes$group_id)

    group_assign <- group_assign %>%
      dplyr::mutate(group_id = factor(group_id, levels = group_order)) %>%
      dplyr::arrange(group_id, gene) %>%
      dplyr::mutate(group_short = short_desc_lookup[as.character(group_id)])

    # For each group, label only the top 5 proteins by min FDR across contrasts
    all_deqms_for_rank <- bind_rows(deqms_results) %>%
      dplyr::filter(protein %in% sig_prots_all) %>%
      dplyr::group_by(protein) %>%
      dplyr::summarise(min_fdr = min(sca.FDR, na.rm = TRUE), .groups = "drop")

    group_assign <- group_assign %>%
      dplyr::left_join(all_deqms_for_rank, by = "protein")

  
    group_assign <- group_assign %>%
      dplyr::group_by(group_id) %>%
      dplyr::mutate(is_key = dplyr::row_number() <= 5) %>%
      dplyr::ungroup()

    group_assign_key <- group_assign %>%
      dplyr::filter(is_key)

    message(sprintf("  Key proteins for heatmaps: %d (from %d total)",
                    nrow(group_assign_key), nrow(group_assign)))

    # Build z-score matrix for sig proteins, averaged per condition
    sig_imp <- imp_mat[sig_prots_all, , drop = FALSE]

    # Average per condition — use explicit matrix construction
    unique_conds <- sort(unique(sdrf_design$condition))
    cond_means <- matrix(NA_real_, nrow = length(sig_prots_all), ncol = length(unique_conds),
                         dimnames = list(sig_prots_all, unique_conds))
    for (cond in unique_conds) {
      cols <- which(sdrf_design$condition == cond)
      cond_means[, cond] <- rowMeans(sig_imp[, cols, drop = FALSE], na.rm = TRUE)
    }


    group_assign_key <- group_assign_key %>% dplyr::filter(protein %in% rownames(cond_means))

    all_groups <- unique(group_assign_key$group_short)
    n_groups <- length(all_groups)
    ug_pal <- setNames(
      c(brewer.pal(min(n_groups, 12), "Set3")[1:min(n_groups, 12)],
        rep("grey90", max(0, n_groups - 12))),
      all_groups
    )
    if ("Other" %in% names(ug_pal)) ug_pal["Other"] <- "grey85"

    # HEATMAP: POLARIZATION EFFECT 

    pol_col_order <- c("KO_M0", "KO_M1", "KO_M2", "WT_M0", "WT_M1", "WT_M2")
    pol_col_order <- intersect(pol_col_order, unique(sdrf_design$condition))

    pol_mat <- cond_means[group_assign_key$protein, pol_col_order, drop = FALSE]
    pol_z <- t(scale(t(pol_mat)))
    pol_z[pol_z > 3]  <- 3
    pol_z[pol_z < -3] <- -3

    # Order proteins by group
    pol_prot_ann <- data.frame(
      group = group_assign_key$group_short[match(rownames(pol_z), group_assign_key$protein)],
      row.names = rownames(pol_z)
    )
    prot_order_pol <- order(pol_prot_ann$group, -rowMeans(pol_z, na.rm = TRUE))
    pol_z <- pol_z[prot_order_pol, , drop = FALSE]
    pol_prot_ann <- pol_prot_ann[prot_order_pol, , drop = FALSE]

    pol_col_labels <- group_assign_key$gene[match(colnames(t(pol_z)), group_assign_key$protein)]

    pol_z_t <- t(pol_z)

    pol_row_ann <- data.frame(
      genotype = c(rep("KO", 3), rep("WT", 3))[1:length(pol_col_order)],
      polstate = rep(c("M0", "M1", "M2"), 2)[1:length(pol_col_order)],
      row.names = pol_col_order
    )

    pol_col_ann <- data.frame(
      Group = pol_prot_ann$group,
      row.names = colnames(pol_z_t)
    )

  
    pol_row_split <- factor(pol_row_ann$genotype, levels = c("KO", "WT"))

    pol_col_split <- factor(pol_col_ann$Group, levels = unique(pol_prot_ann$group))

    n_key_prots_pol <- ncol(pol_z_t)

    svg(file.path(OUT_DIR, "figures", "16_heatmap_polarization_effect_functional_groups.svg"),
        width = max(16, n_key_prots_pol * 0.22 + 6), height = 7)

    ComplexHeatmap::Heatmap(pol_z_t,
      name                          = "Z-score",
      col                           = circlize::colorRamp2(c(-3, 0, 3),
                                                           c("#2166AC", "white", "#D6604D")),
      row_split                     = pol_row_split,
      row_title                     = c("KO (SRGN knockout)", "WT (wild type)"),
      row_title_gp                  = grid::gpar(fontsize = 11, fontface = "bold"),
      row_title_rot                 = 0,
      column_split                  = pol_col_split,
      column_title_gp               = grid::gpar(fontsize = 7),
      column_title_rot              = 45,
      cluster_rows                  = FALSE,
      cluster_columns               = TRUE,
      show_column_names             = TRUE,
      column_labels                 = pol_col_labels,
      column_names_gp               = grid::gpar(fontsize = 6),
      row_names_gp                  = grid::gpar(fontsize = 9),
      left_annotation               = ComplexHeatmap::rowAnnotation(
        df = pol_row_ann[, c("polstate"), drop = FALSE],
        col = list(polstate = c(M0 = "#4DAF4A", M1 = "#FD8D3C", M2 = "#800026")),
        annotation_name_side = "top",
        annotation_legend_param = list(polstate = list(title = "Polarization"))
      ),
      top_annotation                = ComplexHeatmap::HeatmapAnnotation(
        df = data.frame(Group = pol_col_ann$Group),
        col = list(Group = ug_pal)
        ),
      heatmap_legend_param          = list(title = "Z-score", direction = "horizontal",
                                           legend_height = unit(4, "cm")),
      border                        = TRUE,
      use_raster                    = TRUE,
      raster_quality                = 3,
      row_gap                       = unit(3, "mm"),
      column_gap                    = unit(2, "mm")
    ) |> ComplexHeatmap::draw(
      main_heatmap = "Z-score",
      column_title = "Polarization effect on secretome (key proteins per group)\nM0 → M1 → M2 within each genotype (mean per condition, z-scored)",
      heatmap_legend_side = "bottom",
      annotation_legend_side = "right"
    )
    dev.off()
    message("  Saved: 16_heatmap_polarization_effect_functional_groups.svg")

    write.csv(group_assign,
              file.path(OUT_DIR, "tables", "09_functional_group_assignments.csv"),
              row.names = FALSE)
    message("  Saved: 09_functional_group_assignments.csv")

  }, error = function(e) {
    message(sprintf("  Functional-group heatmaps failed: %s", e$message))
    print(traceback())
  })
}
# 6. SESSION INFO

message("")
message("EVERYTHING THAT HAS TRANSPIRED HAS DONE SO ACCORDING TO MY DESIGN")
message(sprintf("All outputs written to: %s/", OUT_DIR))
message("")
message("SESSION INFO")
print(sessionInfo())
