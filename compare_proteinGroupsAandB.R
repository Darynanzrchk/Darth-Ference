required_pkgs <- c("tidyverse", "data.table", "ggrepel", "patchwork", "scales")
new_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(new_pkgs) > 0) install.packages(new_pkgs, repos = "https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(ggrepel)
  library(patchwork); library(scales)
})

if ("MASS" %in% loadedNamespaces()) detach("package:MASS", unload = TRUE)

PATH_A  <- "proteinGroups_A.txt"
PATH_B  <- "proteinGroups_B.txt"
OUT_DIR <- "compareProteinGroups"
FIG_DIR <- file.path(OUT_DIR, "figures")
TAB_DIR <- file.path(OUT_DIR, "tables")

FDR_THRESHOLD <- 0.05
FC_THRESHOLD  <- 1.0
MIN_VALID     <- 2

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)

save_svg <- function(plot, path, width = 7, height = 6) {
  svg(path, width = width, height = height)
  print(plot)
  dev.off()
}

parse_cond <- function(s) {
  vapply(s, function(x) {
    m <- regmatches(x, regexpr("(KO|WT)_(M[012])", x))
    if (length(m) == 1L) m else NA_character_
  }, character(1L), USE.NAMES = FALSE)
}

apply_qc_filters <- function(dt, label) {
  flag_cols <- c("Reverse", "Decoy", "Potential contaminant", "Only identified by site")
  present <- intersect(flag_cols, names(dt))
  before  <- nrow(dt)
  for (col in present) {
    v <- dt[[col]]
    bad <- (!is.na(v)) & (v == "+" | v == TRUE | v == "true" | v == "1")
    dt <- dt[!bad]
  }
  message(sprintf("  [%s] QC: removed %d, retained %d", label, before - nrow(dt), nrow(dt)))
  dt
}

coerce_intensity_cols <- function(dt) {
  for (col in grep("^(LFQ intensity|Intensity |iBAQ |MS/MS count)", names(dt), value = TRUE))
    set(dt, j = col, value = as.double(dt[[col]]))
  dt
}

extract_lfq_matrix <- function(dt, lfq_cols, id_col = "leading_id") {
  m <- as.matrix(dt[, ..lfq_cols]); storage.mode(m) <- "double"
  m[m == 0] <- NA; m <- log2(m)
  rownames(m) <- dt[[id_col]]; colnames(m) <- sub("^LFQ intensity ", "", lfq_cols)
  m
}

run_ttest_per_protein <- function(mA, mB, min_valid = MIN_VALID) {
  proteins <- rownames(mA)
  stopifnot(identical(proteins, rownames(mB)))
  map_dfr(proteins, function(prot) {
    xA <- as.numeric(mA[prot, ]); xB <- as.numeric(mB[prot, ])
    nA <- sum(!is.na(xA)); nB <- sum(!is.na(xB))
    log2FC <- mean(xB, na.rm = TRUE) - mean(xA, na.rm = TRUE)
    p_val <- if (nA >= min_valid && nB >= min_valid)
      tryCatch(t.test(xB, xA, var.equal = FALSE)$p.value, error = function(e) NA_real_)
    else NA_real_
    tibble(protein = prot, mean_log2_A = mean(xA, na.rm = TRUE),
           mean_log2_B = mean(xB, na.rm = TRUE), log2FC = log2FC,
           n_valid_A = nA, n_valid_B = nB, p_value = p_val)
  }) %>%
    mutate(FDR = p.adjust(p_value, method = "BH"),
           significant = !is.na(FDR) & FDR < FDR_THRESHOLD & abs(log2FC) >= FC_THRESHOLD,
           direction = case_when(
             significant & log2FC > 0 ~ "Higher in B (reanalysis)",
             significant & log2FC < 0 ~ "Lower in B (reanalysis)",
             TRUE ~ "Not significant"))
}

# ---- Load & QC ----

dtA <- fread(PATH_A, showProgress = FALSE, na.strings = c("", "NA"))
dtA <- apply_qc_filters(dtA, "A"); dtA <- coerce_intensity_cols(dtA)
dtA[, leading_id := sub(";.*", "", `Protein IDs`)]

dtB <- fread(PATH_B, showProgress = FALSE, na.strings = c("", "NA"))
dtB <- apply_qc_filters(dtB, "B"); dtB <- coerce_intensity_cols(dtB)
dtB[, leading_id := sub(";.*", "", `Protein IDs`)]

# ---- Structure audit ----

col_type_patterns <- c(
  "LFQ intensity" = "^LFQ intensity ", "Intensity" = "^Intensity ",
  "iBAQ" = "^iBAQ ", "MS/MS count" = "^MS/MS count ",
  "Peptide counts" = "^Peptide counts", "Identification type" = "^Identification type ",
  "Sequence coverage" = "^Sequence coverage", "Unique peptides" = "^Unique peptides ",
  "Razor peptides" = "^Razor \\+ unique peptides ", "Peptides" = "^Peptides ")

ctA <- sapply(col_type_patterns, function(p) sum(grepl(p, names(dtA))))
ctB <- sapply(col_type_patterns, function(p) sum(grepl(p, names(dtB))))
structural_audit <- tibble(column_type = names(col_type_patterns),
                           A_original = as.integer(ctA), B_reanalysis = as.integer(ctB),
                           difference = as.integer(ctB) - as.integer(ctA))

lfq_colsA <- grep("^LFQ intensity ", names(dtA), value = TRUE)
lfq_colsB <- grep("^LFQ intensity ", names(dtB), value = TRUE)
sampA <- sub("^LFQ intensity ", "", lfq_colsA)
sampB <- sub("^LFQ intensity ", "", lfq_colsB)
condA <- parse_cond(sampA); condB <- parse_cond(sampB)
shared_conds <- sort(intersect(unique(na.omit(condA)), unique(na.omit(condB))))

sample_map <- bind_rows(
  tibble(dataset = "A (original)", sample = sampA, condition = condA),
  tibble(dataset = "B (reanalysis)", sample = sampB, condition = condB))

write_csv(structural_audit, file.path(TAB_DIR, "01_structural_audit.csv"))
write_csv(sample_map, file.path(TAB_DIR, "01_sample_condition_map.csv"))

# ---- Protein overlap ----

idsA <- dtA$leading_id; idsB <- dtB$leading_id
ids_common <- intersect(idsA, idsB)
ids_only_A <- setdiff(idsA, idsB); ids_only_B <- setdiff(idsB, idsA)
jaccard <- length(ids_common) / length(union(idsA, idsB))

gene_map <- dtA[, .(leading_id, `Gene names`)] %>%
  as_tibble() %>% rename(protein = leading_id, gene = `Gene names`) %>%
  mutate(gene = coalesce(na_if(gene, ""), protein))

overlap_summary <- tibble(
  category = c("A total", "B total", "Common", "Unique to A", "Unique to B", "Jaccard"),
  count = c(length(idsA), length(idsB), length(ids_common),
            length(ids_only_A), length(ids_only_B), round(jaccard, 4)))

write_csv(overlap_summary, file.path(TAB_DIR, "02_overlap_summary.csv"))
write_csv(tibble(leading_id = ids_only_A,
                 gene = dtA$`Gene names`[match(ids_only_A, dtA$leading_id)]),
          file.path(TAB_DIR, "02_proteins_unique_A.csv"))
write_csv(tibble(leading_id = ids_only_B,
                 gene = dtB$`Gene names`[match(ids_only_B, dtB$leading_id)]),
          file.path(TAB_DIR, "02_proteins_unique_B.csv"))

overlap_bar <- tibble(
  category = factor(c("A only\n(original)", "Common\n(both)", "B only\n(reanalysis)"),
                    levels = c("A only\n(original)", "Common\n(both)", "B only\n(reanalysis)")),
  n = c(length(ids_only_A), length(ids_common), length(ids_only_B)),
  fill = c("#2166AC", "#6BAED6", "#D6604D"))

p_overlap <- ggplot(overlap_bar, aes(x = category, y = n, fill = category)) +
  geom_col(width = 0.55, show.legend = FALSE) +
  geom_text(aes(label = comma(n)), vjust = -0.5, size = 4.5, fontface = "bold") +
  scale_fill_manual(values = setNames(overlap_bar$fill, overlap_bar$category)) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.14))) +
  labs(title = "Protein identification overlap between MaxQuant runs",
       subtitle = sprintf("Jaccard = %.3f | %d common / %d total", jaccard,
                          length(ids_common), length(union(idsA, idsB))),
       x = NULL, y = "Number of protein groups") +
  theme_bw(base_size = 13) +
  theme(plot.subtitle = element_text(size = 10, colour = "grey40"),
        panel.grid.major.x = element_blank())

save_svg(p_overlap, file.path(FIG_DIR, "03_protein_overlap_bar.svg"), 6, 5)

# ---- Shared proteins ----

dfA_c <- dtA[leading_id %in% ids_common][order(leading_id)]
dfB_c <- dtB[leading_id %in% ids_common][order(leading_id)]
stopifnot(identical(dfA_c$leading_id, dfB_c$leading_id))

matA <- extract_lfq_matrix(dfA_c, lfq_colsA)
matB <- extract_lfq_matrix(dfB_c, lfq_colsB)

# ---- Missing values ----

missing_per_cond <- bind_rows(lapply(shared_conds, function(cond) {
  iA <- which(condA == cond); iB <- which(condB == cond)
  tibble(condition = cond,
         A_missing_pct = round(mean(is.na(matA[, iA])) * 100, 1),
         B_missing_pct = round(mean(is.na(matB[, iB])) * 100, 1),
         delta_pct = round(mean(is.na(matB[, iB])) * 100 - mean(is.na(matA[, iA])) * 100, 1))
}))

overall_missing <- tibble(
  dataset = c("A (original)", "B (reanalysis)"),
  pct_missing = round(c(mean(is.na(matA)), mean(is.na(matB))) * 100, 2),
  median_log2_lfq = round(c(median(matA, na.rm = TRUE), median(matB, na.rm = TRUE)), 3))

write_csv(overall_missing, file.path(TAB_DIR, "04_missing_value_summary.csv"))
write_csv(missing_per_cond, file.path(TAB_DIR, "04_missing_per_condition.csv"))

# ---- LFQ distributions ----

pal2 <- c("A (original)" = "#2166AC", "B (reanalysis)" = "#D6604D")

long_all <- bind_rows(
  as_tibble(matA, rownames = "protein") %>% pivot_longer(-protein, names_to = "sample", values_to = "log2_LFQ") %>%
    mutate(dataset = "A (original)", condition = parse_cond(sample)),
  as_tibble(matB, rownames = "protein") %>% pivot_longer(-protein, names_to = "sample", values_to = "log2_LFQ") %>%
    mutate(dataset = "B (reanalysis)", condition = parse_cond(sample))
) %>% filter(!is.na(log2_LFQ))

p_density <- ggplot(long_all, aes(x = log2_LFQ, colour = dataset, fill = dataset)) +
  geom_density(alpha = 0.20, linewidth = 0.9) +
  scale_colour_manual(values = pal2) + scale_fill_manual(values = pal2) +
  labs(title = "Log2 LFQ intensity distributions — common proteins",
       x = expression(log[2](LFQ~intensity)), y = "Density", colour = NULL, fill = NULL) +
  theme_bw(base_size = 13) + theme(legend.position = "top")

p_box <- ggplot(long_all, aes(x = sample, y = log2_LFQ, fill = dataset)) +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.3, linewidth = 0.4) +
  scale_fill_manual(values = pal2) + facet_wrap(~ dataset, scales = "free_x", nrow = 2) +
  labs(x = "Sample", y = expression(log[2](LFQ~intensity)), fill = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "none", strip.background = element_rect(fill = "grey92"))

save_svg(p_density / p_box + plot_layout(heights = c(1, 1.6)) + plot_annotation(tag_levels = "A"),
         file.path(FIG_DIR, "05_lfq_distributions.svg"), 11, 9)

# ---- Statistical comparison ----

stats_global <- run_ttest_per_protein(matA, matB)
n_tested <- sum(!is.na(stats_global$p_value))
n_sig <- sum(stats_global$significant, na.rm = TRUE)
n_up <- sum(stats_global$significant & stats_global$log2FC > 0, na.rm = TRUE)
n_down <- sum(stats_global$significant & stats_global$log2FC < 0, na.rm = TRUE)

write_csv(stats_global, file.path(TAB_DIR, "06a_global_statistics.csv"))

stats_percond <- map_dfr(shared_conds, function(cond) {
  iA <- which(condA == cond); iB <- which(condB == cond)
  run_ttest_per_protein(matA[, iA, drop = FALSE], matB[, iB, drop = FALSE]) %>%
    mutate(condition = cond, .before = 1)
})

percond_summary <- stats_percond %>%
  group_by(condition) %>%
  summarise(n_tested = sum(!is.na(p_value)), n_sig = sum(significant, na.rm = TRUE),
            n_up = sum(significant & log2FC > 0, na.rm = TRUE),
            n_down = sum(significant & log2FC < 0, na.rm = TRUE), .groups = "drop")

write_csv(stats_percond, file.path(TAB_DIR, "06b_per_condition_statistics.csv"))
write_csv(percond_summary, file.path(TAB_DIR, "06b_per_condition_summary.csv"))

# ---- Per-condition volcano ----

volcano_colours <- c("Higher in B (reanalysis)" = "#D6604D",
                     "Lower in B (reanalysis)"  = "#2166AC",
                     "Not significant"           = "grey72")

top_percond <- stats_percond %>% filter(significant) %>%
  group_by(condition) %>% slice_min(FDR, n = 5) %>% ungroup() %>% dplyr::select(condition, protein)

plot_percond <- stats_percond %>% filter(!is.na(FDR)) %>%
  mutate(neg_log10_FDR = -log10(FDR), condition = factor(condition, levels = shared_conds)) %>%
  left_join(gene_map, by = "protein") %>%
  left_join(top_percond %>% mutate(is_top = TRUE), by = c("condition", "protein")) %>%
  mutate(label = if_else(!is.na(is_top), gene, NA_character_))

p_volcano_percond <- ggplot(plot_percond,
    aes(x = log2FC, y = neg_log10_FDR, colour = direction, size = direction)) +
  geom_point(alpha = 0.55) +
  geom_text_repel(aes(label = label), size = 2.6, max.overlaps = 15,
                  segment.colour = "grey50", segment.size = 0.3,
                  show.legend = FALSE, na.rm = TRUE) +
  geom_vline(xintercept = c(-FC_THRESHOLD, FC_THRESHOLD), linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = -log10(FDR_THRESHOLD), linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  scale_colour_manual(values = volcano_colours) +
  scale_size_manual(values = c("Higher in B (reanalysis)" = 1.6,
                               "Lower in B (reanalysis)" = 1.6, "Not significant" = 0.6)) +
  facet_wrap(~ condition, ncol = 3) +
  labs(title = "Per-condition volcano: B (reanalysis) vs A (original)",
       subtitle = sprintf("FDR<%.2f & |log2FC|>=%.1f | n=3 vs n=3 (low power)",
                          FDR_THRESHOLD, FC_THRESHOLD),
       x = expression(log[2]~Fold~Change~(B/A)),
       y = expression(-log[10](FDR)), colour = NULL) +
  guides(size = "none") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", strip.background = element_rect(fill = "grey92"),
        plot.subtitle = element_text(size = 9, colour = "grey40"))

save_svg(p_volcano_percond, file.path(FIG_DIR, "07b_volcano_per_condition.svg"), 11, 7)

# ---- Significant proteins table ----

sig_table <- stats_global %>% filter(significant) %>% arrange(FDR) %>%
  left_join(gene_map, by = "protein") %>%
  mutate(across(where(is.double), ~ round(.x, 4))) %>%
  dplyr::select(protein, gene, log2FC, p_value, FDR, direction,
                mean_log2_A, mean_log2_B, n_valid_A, n_valid_B)

write_csv(sig_table, file.path(TAB_DIR, "07_significant_proteins_global.csv"))

# ---- Correlation & PCA ----

meanA <- rowMeans(matA, na.rm = TRUE); meanB <- rowMeans(matB, na.rm = TRUE)
valid_idx <- is.finite(meanA) & is.finite(meanB)
corr_pear <- cor(meanA[valid_idx], meanB[valid_idx], method = "pearson")
corr_spear <- cor(meanA[valid_idx], meanB[valid_idx], method = "spearman")

corr_df <- tibble(protein = names(meanA)[valid_idx],
                  mean_log2_A = meanA[valid_idx], mean_log2_B = meanB[valid_idx]) %>%
  left_join(gene_map, by = "protein") %>%
  mutate(highlight = if_else(protein %in% (stats_global %>% filter(significant) %>% pull(protein)),
                             "Significant", "Not significant"))

p_corr <- ggplot(corr_df, aes(x = mean_log2_A, y = mean_log2_B)) +
  geom_point(data = filter(corr_df, highlight == "Not significant"),
             colour = "grey65", size = 0.8, alpha = 0.5) +
  geom_point(data = filter(corr_df, highlight == "Significant"),
             colour = "#D6604D", size = 1.8, alpha = 0.85) +
  geom_abline(slope = 1, intercept = 0, colour = "black", linetype = "dashed", linewidth = 0.7) +
  geom_smooth(method = "lm", se = TRUE, colour = "#2166AC", linewidth = 0.8, alpha = 0.12) +
  annotate("text", x = min(corr_df$mean_log2_A) + 0.3, y = max(corr_df$mean_log2_B) - 0.5,
           label = sprintf("Pearson r = %.4f\nSpearman rho = %.4f\nn = %d",
                           corr_pear, corr_spear, sum(valid_idx)), hjust = 0, size = 4) +
  labs(title = "Quantitative correlation: A vs B (common proteins)",
       subtitle = "Mean log2 LFQ per protein | red = significant (FDR<0.05, |log2FC|>=1)",
       x = expression(Mean~log[2](LFQ)~" — A (original)"),
       y = expression(Mean~log[2](LFQ)~" — B (reanalysis)")) +
  theme_bw(base_size = 13) + theme(plot.subtitle = element_text(size = 9, colour = "grey40"))

save_svg(p_corr, file.path(FIG_DIR, "08a_correlation_scatter.svg"), 6.5, 6)

matA_pca <- matA; colnames(matA_pca) <- paste0(colnames(matA_pca), "__A")
matB_pca <- matB; colnames(matB_pca) <- paste0(colnames(matB_pca), "__B")
mat_combined <- cbind(matA_pca, matB_pca)
complete_rows <- complete.cases(mat_combined)
mat_pca_in <- mat_combined[complete_rows, ]

if (nrow(mat_pca_in) >= 10) {
  pca_res <- prcomp(t(mat_pca_in), center = TRUE, scale. = TRUE)
  var_exp <- summary(pca_res)$importance[2, 1:2] * 100
  pca_df <- as_tibble(pca_res$x[, 1:2], rownames = "sample_id") %>%
    mutate(run = if_else(str_ends(sample_id, "__A"), "A (original)", "B (reanalysis)"),
           raw_name = str_remove(sample_id, "__(A|B)$"), condition = parse_cond(raw_name))

  cond_colours <- c(KO_M0 = "#1F78B4", KO_M1 = "#33A02C", KO_M2 = "#E31A1C",
                    WT_M0 = "#A6CEE3", WT_M1 = "#B2DF8A", WT_M2 = "#FB9A99")

  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = condition, shape = run, label = raw_name)) +
    geom_point(size = 3.5, alpha = 0.9, stroke = 0.8) +
    geom_text_repel(size = 2.5, max.overlaps = 20, show.legend = FALSE,
                    segment.colour = "grey60", segment.size = 0.3) +
    scale_colour_manual(values = cond_colours) +
    scale_shape_manual(values = c("A (original)" = 16, "B (reanalysis)" = 17)) +
    labs(title = "PCA: A + B samples combined",
         subtitle = sprintf("%d complete-case proteins | circles=A, triangles=B", nrow(mat_pca_in)),
         x = sprintf("PC1 (%.1f%% variance)", var_exp[1]),
         y = sprintf("PC2 (%.1f%% variance)", var_exp[2]),
         colour = "Condition", shape = "Run") +
    theme_bw(base_size = 13) +
    theme(legend.position = "right", plot.subtitle = element_text(size = 9, colour = "grey40"))

  save_svg(p_pca, file.path(FIG_DIR, "08b_pca_plot.svg"), 9, 7)
}

# ---- CV classification & annotated volcano ----

cv_A <- apply(matA, 1, function(x) sd(x, na.rm = TRUE) / abs(mean(x, na.rm = TRUE)))
cv_B <- apply(matB, 1, function(x) sd(x, na.rm = TRUE) / abs(mean(x, na.rm = TRUE)))
cv_table <- tibble(protein = names(cv_A), cv_A = cv_A, cv_B = cv_B, cv_ratio = cv_B / cv_A)

annotated <- stats_global %>%
  left_join(cv_table, by = "protein") %>%
  left_join(gene_map, by = "protein") %>%
  mutate(result_class = case_when(
    !significant                                                  ~ "Not significant",
    significant & log2FC > 0 & cv_A < 0.15 & cv_B < 0.15        ~ "Higher in B — biological",
    significant & log2FC > 0 & (cv_ratio > 2 | cv_ratio < 0.5)  ~ "Higher in B — processing",
    significant & log2FC < 0 & cv_A < 0.15 & cv_B < 0.15        ~ "Lower in B — biological",
    significant & log2FC < 0 & (cv_ratio > 2 | cv_ratio < 0.5)  ~ "Lower in B — processing",
    significant                                                   ~ "Ambiguous",
    TRUE                                                          ~ "Not significant"))

class_summary <- annotated %>% count(result_class, name = "n_proteins") %>% arrange(desc(n_proteins))
write_csv(annotated, file.path(TAB_DIR, "09_annotated_results.csv"))
write_csv(class_summary, file.path(TAB_DIR, "09_class_summary.csv"))

class_colours <- c("Higher in B — biological" = "#B2182B", "Higher in B — processing" = "#EF8A62",
                   "Lower in B — biological"  = "#2166AC", "Lower in B — processing"  = "#67A9CF",
                   "Ambiguous" = "#FEC44F", "Not significant" = "grey75")
class_sizes <- c("Higher in B — biological" = 2.2, "Higher in B — processing" = 2.0,
                 "Lower in B — biological" = 2.2, "Lower in B — processing" = 2.0,
                 "Ambiguous" = 2.0, "Not significant" = 0.7)

top_annot <- annotated %>% filter(significant) %>% arrange(FDR) %>% slice_head(n = 20) %>% pull(protein)

plot_annot <- annotated %>% filter(!is.na(FDR)) %>%
  mutate(neg_log10_FDR = -log10(FDR),
         result_class = factor(result_class,
                               levels = c("Higher in B — biological", "Higher in B — processing",
                                          "Lower in B — biological", "Lower in B — processing",
                                          "Ambiguous", "Not significant")),
         label = if_else(protein %in% top_annot, coalesce(gene, protein), NA_character_))

p_annot <- ggplot(plot_annot,
    aes(x = log2FC, y = neg_log10_FDR, colour = result_class, size = result_class)) +
  geom_point(alpha = 0.70) +
  geom_text_repel(aes(label = label), size = 2.8, max.overlaps = 25,
                  segment.colour = "grey50", segment.size = 0.3,
                  show.legend = FALSE, na.rm = TRUE) +
  geom_vline(xintercept = c(-FC_THRESHOLD, FC_THRESHOLD), linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = -log10(FDR_THRESHOLD), linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  scale_colour_manual(values = class_colours, guide = guide_legend(override.aes = list(size = 3))) +
  scale_size_manual(values = class_sizes) +
  labs(title = "Annotated volcano: B (reanalysis) vs A (original)",
       subtitle = sprintf("CV classification | FDR<%.2f & |log2FC|>=%.1f | %d significant",
                          FDR_THRESHOLD, FC_THRESHOLD, n_sig),
       x = expression(log[2]~Fold~Change~(B/A)),
       y = expression(-log[10](FDR)), colour = NULL) +
  guides(size = "none") +
  theme_bw(base_size = 13) +
  theme(legend.position = "right", plot.subtitle = element_text(size = 9, colour = "grey40"))

save_svg(p_annot, file.path(FIG_DIR, "09_volcano_annotated.svg"), 9, 7)

# ---- CV scatter ----

cv_plot_data <- annotated %>%
  filter(!is.na(cv_A) & !is.na(cv_B) & is.finite(cv_A) & is.finite(cv_B)) %>%
  mutate(result_class = factor(result_class,
                               levels = c("Higher in B — biological", "Higher in B — processing",
                                          "Lower in B — biological", "Lower in B — processing",
                                          "Ambiguous", "Not significant")))

cv_max <- min(max(c(cv_plot_data$cv_A, cv_plot_data$cv_B), na.rm = TRUE), 1.0)
cv_ref <- tibble(cv_A = seq(0.01, cv_max, length.out = 200)) %>%
  mutate(cv_B_2x = 2 * cv_A, cv_B_05x = 0.5 * cv_A)

p_cv <- ggplot(cv_plot_data, aes(x = cv_A, y = cv_B, colour = result_class)) +
  geom_point(alpha = 0.5, size = 0.9) +
  geom_abline(slope = 1, intercept = 0, colour = "grey30", linewidth = 0.6) +
  geom_line(data = cv_ref, aes(x = cv_A, y = cv_B_2x), inherit.aes = FALSE,
            linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_line(data = cv_ref, aes(x = cv_A, y = cv_B_05x), inherit.aes = FALSE,
            linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = 0.15, linetype = "dotted", colour = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = 0.15, linetype = "dotted", colour = "grey50", linewidth = 0.4) +
  annotate("text", x = cv_max * 0.95, y = cv_max * 0.95,
           label = "CV_A = CV_B", hjust = 1, size = 3, colour = "grey30") +
  annotate("text", x = cv_max * 0.95, y = min(cv_max, 2 * cv_max * 0.95),
           label = "CV_B = 2 x CV_A", hjust = 1, size = 3, colour = "grey50") +
  scale_colour_manual(values = class_colours,
                      guide = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  coord_cartesian(xlim = c(0, cv_max), ylim = c(0, cv_max)) +
  labs(title = "CV comparison: A (original) vs B (reanalysis)",
       subtitle = sprintf("Dotted = CV=0.15 | Dashed = CV ratio 2x/0.5x | n = %d", nrow(cv_plot_data)),
       x = "CV — A (original)", y = "CV — B (reanalysis)", colour = NULL) +
  theme_bw(base_size = 13) +
  theme(legend.position = "right", plot.subtitle = element_text(size = 9, colour = "grey40"))

save_svg(p_cv, file.path(FIG_DIR, "09b_cv_scatter.svg"), 8, 7)

message(sprintf("\n[DONE] %d figures in %s | %d tables in %s",
               length(list.files(FIG_DIR)), FIG_DIR,
               length(list.files(TAB_DIR)), TAB_DIR))
