# PXD002582 Reanalysis of mouse BMDM TMT proteomics


PATH_PG   <- "proteinGroups2.txt"
PATH_SDRF <- "SDRF_onlyProteom_CSF.tsv"
OUT_DIR   <- "results"

FDR_THRESHOLD <- 0.05
FC_STRICT     <- 1.0        # log2FC >= 1   ( 2x fold change )
FC_LENIENT    <- log2(1.5)  # log2FC >= 0.585 ( 1.5x fold change )
MIN_PEPTIDES  <- 2
SEED          <- 42

# ---- PACKAGES ----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

BiocManager::install(c("sva", "STRINGdb", "GSVA"), ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(limma)
  library(DEqMS)
  library(sva)
  library(pheatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(viridis)
  library(RColorBrewer)
  library(dendextend)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(ReactomePA)
  library(enrichplot)
  library(AnnotationDbi)
  library(httr)
  library(jsonlite)
  library(STRINGdb)
  library(GSVA)
  library(ggraph)
  library(tidygraph)
  library(igraph)
})

select <- dplyr::select
filter <- dplyr::filter
slice  <- dplyr::slice

dir.create(file.path(OUT_DIR, "figures"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(OUT_DIR, "tables"),  showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(OUT_DIR, "tmp"),     showWarnings=FALSE, recursive=TRUE)

set.seed(SEED)

stage <- function(n, title) {
  bar <- paste(rep("#", 78), collapse="")
  message("\n", bar, "\n#### STAGE ", n, "   ", title, "\n", bar)
}

COND_PAL <- c(GM = "#75A025", GM_LPS = "#FF9400", M = "#0279EE", M_LPS = "#E9ED4C")
SET_PAL  <- c(Set1 = "#000000", Set2 = "#FD9BED", Set3 = "#75A025")
DIR_PAL  <- c(Up = "#FF9400", Down = "#0279EE", NS = "grey72")
LOC_PAL  <- c(UniProt = "#0279EE", GO = "#FF9400", Unresolved = "grey72")
SECRET_PAL <- c(Cytokine = "#FF9400", Chemokine = "#E9ED4C",
                `ECM/remodelling` = "#75A025", `Growth factor` = "#0279EE",
                Other = "grey72")

save_svg <- function(filename, plot_expr, width=8, height=6) {
  path <- file.path(OUT_DIR, "figures", filename)
  svg(path, width=width, height=height)
  on.exit(dev.off())
  print(plot_expr)
}

# Label top significant proteins
volcano_labels <- function(df, n_label = 20) {
  df |>
    mutate(score = ifelse(direction_strict != "NS",
                          (-log10(pmax(sca.adj.pval, 1e-300))) * abs(logFC),
                          NA_real_)) |>
    mutate(label = ifelse(!is.na(score) & rank(-score, ties.method="first") <= n_label &
                          !is.na(gene_symbol) & gene_symbol != "",
                          gene_symbol, NA_character_)) |>
    pull(label)
}

# ---- STAGE 1 QUALITY CONTROL ----
pg <- fread(PATH_PG, showProgress=FALSE, na.strings=c("","NA"))
N_PROT_RAW_GLOBAL <- nrow(pg)

sdrf <- fread(PATH_SDRF, showProgress=FALSE)

## QC filter on protein groups
n_before <- nrow(pg)
pg <- pg[(is.na(Decoy)                  | Decoy                  != "+") &
        (is.na(`Potential contaminant`)| `Potential contaminant`!= "+") &
        (is.na(`Only identified by site`)|`Only identified by site` != "+") &
        (`Razor + unique peptides`     >= MIN_PEPTIDES)]
message(sprintf("  QC: %d -> %d proteins after decoy/contaminant/site/peptide filter",
                n_before, nrow(pg)))
N_PROT_QC_GLOBAL <- nrow(pg)

## Sample matrix from per-set active channels
design_tbl <- bind_rows(
  tibble(set="Set1", channel=1, condition="GM",     rep=1),
  tibble(set="Set1", channel=2, condition="GM_LPS", rep=1),
  tibble(set="Set1", channel=5, condition="M",      rep=1),
  tibble(set="Set1", channel=6, condition="M_LPS",  rep=1),
  tibble(set="Set2", channel=3, condition="GM",     rep=2),
  tibble(set="Set2", channel=4, condition="GM_LPS", rep=2),
  tibble(set="Set2", channel=5, condition="M",      rep=2),
  tibble(set="Set2", channel=6, condition="M_LPS",  rep=2),
  tibble(set="Set3", channel=3, condition="GM",     rep=3),
  tibble(set="Set3", channel=4, condition="GM_LPS", rep=3),
  tibble(set="Set3", channel=5, condition="M",      rep=3),
  tibble(set="Set3", channel=6, condition="M_LPS",  rep=3)
) |>
  mutate(
    sample = paste(condition, paste0("rep", rep), sep="_"),
    column = paste0("Reporter intensity corrected ", channel, " ", set)
  )

write.csv(design_tbl, file.path(OUT_DIR, "tables", "01_design_summary.csv"),
          row.names=FALSE)

stopifnot(all(design_tbl$column %in% names(pg)))

## Validate empty TMT channels
empty_channels <- list(
  Set1 = c("Reporter intensity corrected 3 Set1",
           "Reporter intensity corrected 4 Set1"),
  Set2 = c("Reporter intensity corrected 1 Set2",
           "Reporter intensity corrected 2 Set2"),
  Set3 = c("Reporter intensity corrected 1 Set3",
           "Reporter intensity corrected 2 Set3")
)
for (set_name in names(empty_channels)) {
  for (col in empty_channels[[set_name]]) {
    if (!col %in% names(pg)) next
    pct_nonzero <- mean(pg[[col]] > 0, na.rm = TRUE)
    if (pct_nonzero > 0.05)
      warning(sprintf(
        "TMT channel '%s' expected empty but %.1f%% of rows are non-zero. Check channel layout.",
        col, pct_nonzero * 100
      ))
  }
}

mat <- as.matrix(pg[, design_tbl$column, with=FALSE])
storage.mode(mat) <- "double"
colnames(mat) <- design_tbl$sample
rownames(mat) <- pg$`Protein IDs`
mat[mat == 0] <- NA

cond_factor  <- design_tbl$condition
unique_conds <- unique(cond_factor)

## Prefilter: require >=2 valid values in >=1 condition
valid_per_cond <- sapply(unique_conds, function(cc) {
  cols <- which(cond_factor == cc)
  rowSums(!is.na(mat[, cols, drop=FALSE]))
})
keep <- rowSums(valid_per_cond >= 2) >= 1
n_pref <- nrow(mat)
mat <- mat[keep, , drop=FALSE]
pg  <- pg[keep, ]
N_PROT_PREFILT_GLOBAL <- nrow(mat)

## Log2 + median normalization
log2_mat_raw <- log2(mat)
sample_meds  <- apply(log2_mat_raw, 2, median, na.rm=TRUE)
global_med   <- median(log2_mat_raw, na.rm=TRUE)
log2_mat     <- sweep(log2_mat_raw, 2, sample_meds - global_med, FUN="-")

## Missingness structure
valid_per_cond_post <- sapply(unique_conds, function(cc) {
  cols <- which(cond_factor == cc)
  rowSums(!is.na(log2_mat[, cols, drop=FALSE]))
})

mec_mask <- matrix(FALSE, nrow=nrow(log2_mat), ncol=ncol(log2_mat),
                   dimnames=dimnames(log2_mat))
for (cc in unique_conds) {
  cols <- which(cond_factor == cc)
  prot_missing_in_cond <- valid_per_cond_post[, cc] == 0
  mec_mask[prot_missing_in_cond, cols] <- TRUE
}
mec_mask <- mec_mask & is.na(log2_mat)

n_mec_cells <- sum(mec_mask)

## Imputation
impute_minprob <- function(x, q=0.01, width=0.3) {
  obs <- x[!is.na(x)]
  if (length(obs) < 5) return(x)
  mu <- as.numeric(quantile(obs, q, na.rm=TRUE))
  sd <- mad(obs, na.rm=TRUE) * width
  if (!is.finite(sd) || sd == 0) sd <- 0.1
  na_idx <- is.na(x)
  x[na_idx] <- rnorm(sum(na_idx), mean=mu, sd=sd)
  x
}

set.seed(SEED)
imp_mat <- apply(log2_mat, 2, impute_minprob)
rownames(imp_mat) <- rownames(log2_mat)

if (n_mec_cells > 0) {
  all_obs <- log2_mat[!is.na(log2_mat)]
  floor_value <- min(all_obs) - 0.5
  imp_mat[mec_mask] <- floor_value
}

write.csv(as.data.frame(imp_mat),
          file.path(OUT_DIR, "tables", "02_normalized_imputed_matrix.csv"))

## ComBat-adjusted matrix (visualization only)
combat_mat <- tryCatch({
  mod <- model.matrix(~ condition, data=data.frame(condition=cond_factor))
  ComBat(dat=imp_mat, batch=design_tbl$set, mod=mod, par.prior=TRUE,
         prior.plots=FALSE)
}, error = function(e) {
  imp_mat
})
write.csv(as.data.frame(combat_mat),
          file.path(OUT_DIR, "tables", "03_combat_adjusted_matrix.csv"))

## QC figures
ann_df <- design_tbl |>
  as.data.frame() |>
  column_to_rownames("sample") |>
  select(set, condition)
ann_colors <- list(set=SET_PAL, condition=COND_PAL)

### FIGURE 01   Missingness vs intensity

miss_intensity_df <- tibble(
  protein     = rownames(log2_mat),
  mean_log2   = rowMeans(log2_mat, na.rm=TRUE),
  n_missing   = rowSums(is.na(log2_mat)),
  pct_missing = rowSums(is.na(log2_mat)) / ncol(log2_mat) * 100
) |>
  drop_na(mean_log2)

p01 <- ggplot(miss_intensity_df, aes(mean_log2, n_missing)) +
  geom_jitter(aes(colour=pct_missing), alpha=0.5, size=0.8,
              height=0.15, width=0) +
  geom_smooth(method="loess", se=TRUE, colour="#000000", linewidth=0.7,
              fill="grey80", span=0.4) +
  scale_colour_viridis_c(
    option="plasma", direction=-1,
    guide=guide_colourbar(
      barwidth  = unit(0.4, "cm"),
      barheight = unit(4,   "cm"),
      title.position = "top",
      title.hjust    = 0.5
    )
  ) +
  scale_y_continuous(breaks=0:max(miss_intensity_df$n_missing)) +
  labs(title="Missingness vs mean log2 intensity",
       subtitle=sprintf("%d proteins; LOESS overlay; low-intensity proteins miss more often (left-censored)",
                        nrow(miss_intensity_df)),
       x=expression("Mean "*log[2]*" reporter intensity (observed values only)"),
       y="Number of missing samples (out of 12)",
       colour="% missing") +
  theme_bw(11) +
  theme(panel.grid.minor=element_blank(),
        legend.margin=margin(0, 0, 0, 4))
save_svg("01_missingness_vs_intensity.svg", p01, 10, 6)   # width 9->10

### FIGURE 02   Boxplots: Raw vs Normalized vs Imputed
long_raw  <- as_tibble(log2_mat_raw, rownames="protein") |>
  pivot_longer(-protein, names_to="sample", values_to="log2") |>
  mutate(stage="1. Raw (log2)") |> drop_na()
long_norm <- as_tibble(log2_mat, rownames="protein") |>
  pivot_longer(-protein, names_to="sample", values_to="log2") |>
  mutate(stage="2. Median-normalized") |> drop_na()
long_imp  <- as_tibble(imp_mat, rownames="protein") |>
  pivot_longer(-protein, names_to="sample", values_to="log2") |>
  mutate(stage="3. Normalized + imputed")
long_all  <- bind_rows(long_raw, long_norm, long_imp) |>
  left_join(design_tbl |> select(sample, set, condition), by="sample") |>
  mutate(stage = factor(stage,
                        levels=c("1. Raw (log2)",
                                 "2. Median-normalized",
                                 "3. Normalized + imputed")))

p02 <- ggplot(long_all, aes(sample, log2, fill=condition)) +
  geom_boxplot(outlier.size=0.3, outlier.alpha=0.3, linewidth=0.4) +
  scale_fill_manual(values=COND_PAL) +
  facet_wrap(~stage, nrow=3, scales="free_y") +
  labs(title="Intensity distributions: raw -> normalized -> imputed",
       x=NULL, y=expression(log[2]*" reporter intensity"), fill="Condition") +
  theme_bw(11) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=8),
        strip.background=element_rect(fill="grey92"))
save_svg("02_boxplots_raw_norm_imputed.svg", p02, 11, 9)

### FIGURE 03   Sample correlation heatmap
cor_mat <- cor(imp_mat, method="spearman")
svg(file.path(OUT_DIR, "figures", "03_sample_correlation_heatmap.svg"),
    width=8, height=7)
pheatmap(cor_mat,
         color = colorRampPalette(c("#0279EE","white","#FF9400"))(100),
         annotation_col = ann_df, annotation_row = ann_df,
         annotation_colors = ann_colors,
         display_numbers = TRUE, number_format = "%.2f",
         fontsize_number = 7, fontsize = 9,
         main = "Sample-sample Spearman correlation (normalized + imputed)",
         border_color = NA)
dev.off()

### FIGURES 04-06   PCA: raw / normalized / ComBat
do_pca_plot <- function(M, title, filename, subtitle=NULL) {
  pc <- prcomp(t(M), center=TRUE, scale.=TRUE)
  ve <- summary(pc)$importance[2, 1:2] * 100
  df <- as_tibble(pc$x[, 1:2], rownames="sample") |>
    left_join(design_tbl |> select(sample, set, condition), by="sample")
  p <- ggplot(df, aes(PC1, PC2, colour=condition, shape=set, label=sample)) +
    geom_point(size=4, stroke=1.2) +
    geom_text_repel(size=2.6, max.overlaps=15, show.legend=FALSE) +
    scale_colour_manual(values=COND_PAL) +
    scale_shape_manual(values=c(Set1=16, Set2=17, Set3=15)) +
    labs(title=title, subtitle=subtitle,
         x=sprintf("PC1 (%.1f%%)", ve[1]),
         y=sprintf("PC2 (%.1f%%)", ve[2])) +
    theme_bw(11)
  save_svg(filename, p, 8, 6)
  pc
}

naive_fill <- function(M) {
  out <- M
  for (j in seq_len(ncol(out))) {
    cm <- mean(out[, j], na.rm=TRUE)
    out[is.na(out[, j]), j] <- cm
  }
  out[apply(out, 1, var) > 0, , drop=FALSE]
}

pca_raw    <- do_pca_plot(naive_fill(log2_mat_raw),
                          "PCA (raw log2, NA mean-filled)",
                          "04_pca_raw.svg",
                          "Per-set median offset is large; raw data still reflects batch and biology mixed together")
pca_norm   <- do_pca_plot(imp_mat,
                          "PCA (normalized + imputed)",
                          "05_pca_normalized.svg",
                          "After median normalization; batch (set) still partially separates samples on PC1")
pca_combat <- do_pca_plot(combat_mat,
                          "PCA (ComBat-adjusted, visualization only)",
                          "06_pca_combat.svg",
                          "After ComBat: batch effect removed; biological clustering by condition becomes visible")

vp <- bind_rows(
  tibble(stage="Before ComBat",
         PC=paste0("PC",1:4),
         r2_set =sapply(1:4, function(k) summary(lm(pca_norm$x[,k]   ~ design_tbl$set))$r.squared),
         r2_cond=sapply(1:4, function(k) summary(lm(pca_norm$x[,k]   ~ design_tbl$condition))$r.squared)),
  tibble(stage="After ComBat",
         PC=paste0("PC",1:4),
         r2_set =sapply(1:4, function(k) summary(lm(pca_combat$x[,k] ~ design_tbl$set))$r.squared),
         r2_cond=sapply(1:4, function(k) summary(lm(pca_combat$x[,k] ~ design_tbl$condition))$r.squared))
)
print(vp)

# ---- LIMMA + DEqMS FIT ----

stage("DE fit", "Limma + DEqMS fit (design ~0 + condition + set)")

group <- factor(design_tbl$condition, levels = c("GM","GM_LPS","M","M_LPS"))
set_f <- factor(design_tbl$set,       levels = c("Set1","Set2","Set3"))
design_mat <- model.matrix(~ 0 + group + set_f)
colnames(design_mat) <- c("GM","GM_LPS","M","M_LPS","Set2","Set3")

contrasts_mat <- makeContrasts(
  M_vs_GM           = M - GM,
  GM_LPS_vs_GM      = GM_LPS - GM,
  M_LPS_vs_M        = M_LPS - M,
  M_LPS_vs_GM_LPS   = M_LPS - GM_LPS,
  interaction       = (M_LPS - M) - (GM_LPS - GM),
  levels = design_mat
)

fit  <- lmFit(imp_mat, design_mat)
fit2 <- contrasts.fit(fit, contrasts_mat)
fit2 <- eBayes(fit2, trend=TRUE, robust=TRUE)
fit2$count <- pg$`Razor + unique peptides`
fit2$count[is.na(fit2$count) | fit2$count < 1] <- 1L
fit2 <- spectraCounteBayes(fit2)

## Annotation
first_gene <- ifelse(is.na(pg$`Gene names`) | pg$`Gene names` == "",
                     NA_character_,
                     vapply(strsplit(pg$`Gene names`, ";"), `[`, character(1), 1))
.valid_genes <- unique(first_gene[!is.na(first_gene) & first_gene != ""])
sym2entrez <- AnnotationDbi::select(org.Mm.eg.db, keys=.valid_genes,
                                    columns=c("ENTREZID"), keytype="SYMBOL")
sym2entrez <- sym2entrez[!duplicated(sym2entrez$SYMBOL) & !is.na(sym2entrez$ENTREZID), ]
sym2entrez_v <- setNames(sym2entrez$ENTREZID, sym2entrez$SYMBOL)

annot_df <- tibble(
  protein_ids           = pg$`Protein IDs`,
  uniprot_lead          = vapply(strsplit(pg$`Protein IDs`, ";"), `[`, character(1), 1),
  gene_symbol           = first_gene,
  entrez_id             = unname(sym2entrez_v[first_gene]),
  razor_unique_peptides = pg$`Razor + unique peptides`,
  protein_names         = pg$`Protein names`
)

## Build DE tables (DEqMS, dual thresholds)
deqms_table <- function(fit_obj, coef, annot) {
  tt <- outputResult(fit_obj, coef = coef)
  tt$protein_ids <- rownames(tt)
  tt_tbl <- as_tibble(tt) |>
    left_join(annot, by="protein_ids") |>
    mutate(
      direction_strict = case_when(
        sca.adj.pval < FDR_THRESHOLD & logFC >=  FC_STRICT ~ "Up",
        sca.adj.pval < FDR_THRESHOLD & logFC <= -FC_STRICT ~ "Down",
        TRUE ~ "NS"),
      direction_lenient = case_when(
        sca.adj.pval < FDR_THRESHOLD & logFC >=  FC_LENIENT ~ "Up",
        sca.adj.pval < FDR_THRESHOLD & logFC <= -FC_LENIENT ~ "Down",
        TRUE ~ "NS"),
      sig_strict  = direction_strict  != "NS",
      sig_lenient = direction_lenient != "NS"
    ) |>
    arrange(sca.adj.pval)
  priority_cols <- c("protein_ids","uniprot_lead","gene_symbol","entrez_id",
                     "protein_names","razor_unique_peptides",
                     "logFC","AveExpr","sca.t","sca.P.Value","sca.adj.pval",
                     "direction_strict","sig_strict","direction_lenient","sig_lenient")
  remaining <- setdiff(names(tt_tbl), priority_cols)
  tt_tbl[, c(priority_cols, remaining)]
}

# Run contrasts
de_M_vs_GM         <- deqms_table(fit2, coef="M_vs_GM",          annot=annot_df)
de_GM_LPS_vs_GM    <- deqms_table(fit2, coef="GM_LPS_vs_GM",     annot=annot_df)
de_M_LPS_vs_M      <- deqms_table(fit2, coef="M_LPS_vs_M",       annot=annot_df)
de_M_LPS_vs_GM_LPS <- deqms_table(fit2, coef="M_LPS_vs_GM_LPS",  annot=annot_df)
de_interaction     <- deqms_table(fit2, coef="interaction",      annot=annot_df)

write.csv(de_M_vs_GM,         file.path(OUT_DIR, "tables", "04_deqms_M_vs_GM.csv"),         row.names=FALSE)
write.csv(de_GM_LPS_vs_GM,    file.path(OUT_DIR, "tables", "05_deqms_GM_LPS_vs_GM.csv"),    row.names=FALSE)
write.csv(de_M_LPS_vs_M,      file.path(OUT_DIR, "tables", "06_deqms_M_LPS_vs_M.csv"),      row.names=FALSE)
write.csv(de_M_LPS_vs_GM_LPS, file.path(OUT_DIR, "tables", "07_deqms_M_LPS_vs_GM_LPS.csv"), row.names=FALSE)
write.csv(de_interaction,     file.path(OUT_DIR, "tables", "08_deqms_interaction.csv"),     row.names=FALSE)

## Summary table
count_dir <- function(de, name) {
  tibble(
    contrast     = name,
    n_total      = nrow(de),
    up_strict    = sum(de$direction_strict  == "Up"),
    down_strict  = sum(de$direction_strict  == "Down"),
    up_lenient   = sum(de$direction_lenient == "Up"),
    down_lenient = sum(de$direction_lenient == "Down")
  )
}
dual_summary <- bind_rows(
  count_dir(de_M_vs_GM,         "M_vs_GM"),
  count_dir(de_GM_LPS_vs_GM,    "GM_LPS_vs_GM"),
  count_dir(de_M_LPS_vs_M,      "M_LPS_vs_M"),
  count_dir(de_M_LPS_vs_GM_LPS, "M_LPS_vs_GM_LPS"),
  count_dir(de_interaction,     "interaction")
)
write.csv(dual_summary, file.path(OUT_DIR, "tables", "16_dual_threshold_summary.csv"),
          row.names=FALSE)

# ---- STAGE 2   M vs GM ----

n_up      <- sum(de_M_vs_GM$direction_strict  == "Up")
n_down    <- sum(de_M_vs_GM$direction_strict  == "Down")
n_up_len  <- sum(de_M_vs_GM$direction_lenient == "Up")
n_dn_len  <- sum(de_M_vs_GM$direction_lenient == "Down")

## FIGURE 07   Volcano M vs GM
volcano_df_07 <- de_M_vs_GM |>
  mutate(neg_log10_p = -log10(pmax(sca.adj.pval, 1e-300)))
volcano_df_07$label <- volcano_labels(volcano_df_07, n_label = 20)

p07 <- ggplot(volcano_df_07, aes(logFC, neg_log10_p, colour=direction_strict)) +
  geom_point(alpha=0.7, size=1.4) +
  geom_text_repel(aes(label=label), size=3, max.overlaps=25,
                  box.padding=0.3, show.legend=FALSE, na.rm=TRUE) +
  geom_vline(xintercept=c(-FC_STRICT, FC_STRICT), linetype="dashed", colour="grey50") +
  geom_vline(xintercept=c(-FC_LENIENT, FC_LENIENT), linetype="dotted", colour="grey60") +
  geom_hline(yintercept=-log10(FDR_THRESHOLD), linetype="dashed", colour="grey50") +
  scale_colour_manual(values=DIR_PAL) +
  labs(title="Volcano: M-CSF vs GM-CSF (unstimulated)",
       subtitle=sprintf("Strict: Up=%d, Down=%d (|log2FC|>=1)   Lenient: Up=%d, Down=%d (|log2FC|>=0.585)",
                        n_up, n_down, n_up_len, n_dn_len),
       x=expression(log[2]*" fold change (M  -  GM)"),
       y=expression(-log[10]*" sca.adj.pval"),
       colour=NULL) +
  theme_bw(11)
save_svg("07_volcano_M_vs_GM.svg", p07, 8, 6)

## FIGURE 08   DE heatmap (top 25 by score, M vs GM)
sample_order_stage2 <- design_tbl |>
  arrange(factor(condition, levels=c("GM","M","GM_LPS","M_LPS")), rep) |>
  pull(sample)

hm08_pool <- de_M_vs_GM |>
  filter(sig_strict) |>
  mutate(score = (-log10(pmax(sca.adj.pval, 1e-300))) * abs(logFC))

if (nrow(hm08_pool) < 10) {
  hm08_pool <- de_M_vs_GM |>
    filter(sig_lenient) |>
    mutate(score = (-log10(pmax(sca.adj.pval, 1e-300))) * abs(logFC))
  message("  Fig 08: fewer than 10 strict-sig proteins; using lenient threshold")
}

hm08_up   <- hm08_pool |> filter(direction_strict == "Up"   | direction_lenient == "Up")   |> arrange(desc(score)) |> head(13)
hm08_down <- hm08_pool |> filter(direction_strict == "Down" | direction_lenient == "Down") |> arrange(desc(score)) |> head(13)
top25_M_vs_GM <- bind_rows(hm08_up, hm08_down) |>
  distinct(protein_ids, .keep_all=TRUE) |>
  arrange(desc(score)) |>
  head(25)

if (nrow(top25_M_vs_GM) >= 2) {
  hm_mat <- combat_mat[top25_M_vs_GM$protein_ids, sample_order_stage2, drop=FALSE]
  rn <- ifelse(!is.na(top25_M_vs_GM$gene_symbol) & top25_M_vs_GM$gene_symbol != "",
               top25_M_vs_GM$gene_symbol, top25_M_vs_GM$uniprot_lead)
  rownames(hm_mat) <- make.unique(rn)
  z <- t(scale(t(hm_mat)))
  ann_hm <- ann_df[sample_order_stage2, , drop=FALSE]


  dir_annot <- data.frame(
    direction = ifelse(top25_M_vs_GM$logFC > 0, "Up", "Down"),
    row.names = rownames(hm_mat)
  )
  ann_colors_08 <- c(ann_colors, list(direction=c(Up="#FF9400", Down="#0279EE")))

  svg(file.path(OUT_DIR, "figures", "08_de_heatmap_M_vs_GM.svg"),
      width=10, height=9)
  print(pheatmap(z,
           color = colorRampPalette(c("#0279EE","white","#FF9400"))(100),
           breaks = seq(-2, 2, length.out=101),
           annotation_col = ann_hm,
           annotation_row = dir_annot,
           annotation_colors = ann_colors_08,
           cluster_rows = TRUE, cluster_cols = FALSE,
           clustering_distance_rows = "correlation",
           show_rownames = TRUE, fontsize_row = 9, fontsize_col = 9,
           main = sprintf("Top %d DEPs: M-CSF vs GM-CSF (score = -log10(FDR) x |log2FC|)",
                          nrow(top25_M_vs_GM)),
           border_color = NA))
  dev.off()
}

# ---- STAGE 3   FUNCTIONAL INTERPRETATION (GSEA + STRING) ----

## GSEA function

run_gsea <- function(de_tt, contrast_id) {
  rk_df <- de_tt |>
    filter(!is.na(entrez_id)) |>
    mutate(rank_stat = sca.t) |>
    group_by(entrez_id) |>
    summarise(rank_stat = rank_stat[which.max(abs(rank_stat))], .groups="drop") |>
    arrange(desc(rank_stat))
  if (nrow(rk_df) < 50) {
    return(tibble())
  }
  rk <- setNames(rk_df$rank_stat, rk_df$entrez_id)
  rk <- rk[!duplicated(names(rk))]

  set.seed(SEED)  
  gg <- tryCatch(gseGO(geneList=rk, OrgDb=org.Mm.eg.db, ont="BP",
                        keyType="ENTREZID", pAdjustMethod="BH",
                        pvalueCutoff=0.05, verbose=FALSE),
                 error=function(e) NULL)
  set.seed(SEED)
  gk <- tryCatch(gseKEGG(geneList=rk, organism="mmu", pAdjustMethod="BH",
                          pvalueCutoff=0.05, verbose=FALSE),
                 error=function(e) NULL)
  set.seed(SEED)
  gr <- tryCatch(gsePathway(geneList=rk, organism="mouse", pAdjustMethod="BH",
                              pvalueCutoff=0.05, verbose=FALSE),
                 error=function(e) NULL)

  out <- list()
  for (nm in c("GO","KEGG","Reactome")) {
    obj <- switch(nm, GO=gg, KEGG=gk, Reactome=gr)
    if (!is.null(obj) && nrow(as.data.frame(obj)) > 0) {
      df <- as.data.frame(obj) |>
        mutate(source=nm, contrast=contrast_id) |>
        select(contrast, source, ID, Description, NES, pvalue, p.adjust,
               core_enrichment, setSize) |>
        as_tibble()
      out[[nm]] <- df
    }
  }
  if (length(out) == 0) return(tibble())
  bind_rows(out)
}

gsea_all <- bind_rows(
  run_gsea(de_M_vs_GM,         "M_vs_GM"),
  run_gsea(de_GM_LPS_vs_GM,    "GM_LPS_vs_GM"),
  run_gsea(de_M_LPS_vs_M,      "M_LPS_vs_M"),
  run_gsea(de_M_LPS_vs_GM_LPS, "M_LPS_vs_GM_LPS"),
  run_gsea(de_interaction,     "interaction")
)
if (nrow(gsea_all) == 0) {
  warning("All GSEA runs returned 0 terms.")
} else {
  write.csv(gsea_all,
            file.path(OUT_DIR, "tables", "09_gsea_combined_all_contrasts.csv"),
            row.names=FALSE)
}

## FIGURE 09   GSEA NES heatmap
if (nrow(gsea_all) > 0) {
  contrast_levels <- c("M_vs_GM","GM_LPS_vs_GM","M_LPS_vs_M",
                       "M_LPS_vs_GM_LPS","interaction")

  sig_terms <- gsea_all |>
    filter(p.adjust < 0.05) |>
    distinct(source, ID, Description)

  nes_per_term <- gsea_all |>
    semi_join(sig_terms, by=c("source","ID","Description")) |>
    group_by(source, ID, Description) |>
    summarise(max_abs_NES = max(abs(NES), na.rm=TRUE), .groups="drop")

  top_per_source <- nes_per_term |>
    group_by(source) |>
    arrange(desc(max_abs_NES)) |>
    slice_head(n=10) |>
    ungroup() |>
    mutate(Description = str_trunc(Description, 55))

  if (nrow(top_per_source) > 0) {
    full_grid <- expand_grid(
      contrast = factor(contrast_levels, levels=contrast_levels),
      top_per_source |> select(source, ID, Description)
    )

    nes_long <- gsea_all |>
      mutate(Description = str_trunc(Description, 55)) |>
      semi_join(top_per_source, by=c("source","ID","Description")) |>
      mutate(contrast = factor(contrast, levels=contrast_levels),
             sig = p.adjust < 0.05) |>
      select(contrast, source, ID, Description, NES, p.adjust, sig)

    nes_long <- full_grid |>
      left_join(nes_long, by=c("contrast","source","ID","Description"))

    row_order_df <- nes_long |>
      filter(contrast == "M_vs_GM") |>
      arrange(source, NES) |>
      mutate(row_rank = row_number()) |>
      select(source, ID, Description, row_rank)
    nes_long <- nes_long |>
      left_join(row_order_df, by=c("source","ID","Description")) |>
      arrange(source, row_rank) |>
      mutate(Description = factor(Description, levels=unique(Description)))

    p09 <- ggplot(nes_long,
                  aes(contrast, Description, fill=NES)) +
      geom_tile(colour="white", linewidth=0.3) +
      geom_text(aes(label=ifelse(!is.na(sig) & sig & !is.na(NES),
                                  sprintf("%.1f", NES), "")),
                size=2.6, colour="grey20") +
      facet_grid(source ~ ., scales="free_y", space="free_y", switch="y") +
      scale_fill_gradient2(
        low="#0279EE", mid="white", high="#FF9400",
        midpoint=0, na.value="grey90",
        limits=c(-3, 3), oob=scales::squish,
        guide=guide_colourbar(
          title          = "NES",
          barwidth       = unit(0.5, "cm"),
          barheight      = unit(5,   "cm"),
          title.position = "top",
          title.hjust    = 0.5
        )
      ) +
      labs(title="GSEA NES across contrasts",
           subtitle="Top 10 terms per source by max |NES|; only significant cells labelled (p.adj<0.05); grey = not enriched",
           x=NULL, y=NULL, fill="NES") +
      theme_bw(10) +
      theme(
        axis.text.x      = element_text(angle=30, hjust=1),
        strip.placement  = "outside",
        strip.background = element_rect(fill="grey92", colour="grey50"),
        strip.text.y.left= element_text(angle=90, size=11, face="bold", colour="black"),
        legend.position  = "right",
        legend.box       = "vertical",
        plot.margin      = margin(5, 25, 5, 5, "mm")   # B7: right margin for legend
      )
    save_svg("09_gsea_nes_heatmap_multicontrast.svg", p09, 13, 14)   # width 11->13
  }
}

## FIGURE 10   STRING network (M vs GM)
sig_M_vs_GM <- de_M_vs_GM |>
  filter(sig_strict, !is.na(gene_symbol), gene_symbol != "") |>
  arrange(sca.adj.pval)

if (nrow(sig_M_vs_GM) >= 5) {

  string_dir <- file.path(OUT_DIR, "tmp", "stringdb")
  dir.create(string_dir, showWarnings=FALSE, recursive=TRUE)

  string_db <- tryCatch(
    STRINGdb$new(version="12.0", species=10090,
                  score_threshold=700, input_directory=string_dir),
    error=function(e) NULL)

  if (!is.null(string_db)) {
    sig_M_vs_GM_capped <- head(sig_M_vs_GM, 400)
    map_in <- sig_M_vs_GM_capped |>
      select(gene_symbol, logFC, sca.adj.pval, direction_strict) |>
      as.data.frame()

    map_res <- tryCatch(
      string_db$map(map_in, "gene_symbol", removeUnmappedRows=TRUE),
      error=function(e) NULL)

    if (!is.null(map_res) && nrow(map_res) > 0) {

      edges <- tryCatch(
        string_db$get_interactions(map_res$STRING_id),
        error=function(e) NULL)

      if (!is.null(edges) && nrow(edges) > 0) {
        id2sym <- setNames(map_res$gene_symbol, map_res$STRING_id)
        edges$from_sym <- id2sym[edges$from]
        edges$to_sym   <- id2sym[edges$to]
        edges <- edges |>
          filter(!is.na(from_sym), !is.na(to_sym), from_sym != to_sym)

        nodes_df <- map_res |>
          mutate(label = gene_symbol,
                 dir   = direction_strict) |>
          select(label, STRING_id, gene_symbol, logFC, sca.adj.pval, dir) |>
          distinct(label, .keep_all = TRUE)

        edges_unique <- edges |>
          select(from_sym, to_sym, combined_score) |>
          distinct(from_sym, to_sym, .keep_all=TRUE) |>
          filter(from_sym %in% nodes_df$label, to_sym %in% nodes_df$label)

        if (nrow(edges_unique) > 0) {
          write.csv(edges_unique, file.path(OUT_DIR, "tables", "10_string_edges.csv"),
                    row.names=FALSE)
          write.csv(nodes_df,     file.path(OUT_DIR, "tables", "10_string_nodes.csv"),
                    row.names=FALSE)

          g <- igraph::graph_from_data_frame(
                  d=edges_unique |> dplyr::rename(from=from_sym, to=to_sym),
                  vertices=nodes_df |> dplyr::rename(name=label),
                  directed=FALSE)
          deg <- igraph::degree(g)
          g_lcc <- igraph::induced_subgraph(g, which(deg > 0))
          if (igraph::vcount(g_lcc) >= 5) {
            tg <- tidygraph::as_tbl_graph(g_lcc) |>
              tidygraph::activate(nodes) |>
              mutate(degree = tidygraph::centrality_degree())

            node_tbl <- as_tibble(tg) |>
              mutate(label_show = ifelse(rank(-abs(logFC),
                                              ties.method="first") <= 30,
                                          gene_symbol, NA_character_))
            tg <- tg |>
              tidygraph::activate(nodes) |>
              mutate(label_show = node_tbl$label_show)

            set.seed(SEED)
            p10 <- ggraph(tg, layout="fr") +
              ggraph::geom_edge_link(aes(width=combined_score/1000),
                                      alpha=0.25, colour="grey30") +
              ggraph::scale_edge_width(range=c(0.2, 1.4)) +
              geom_node_point(aes(size=degree, fill=dir),
                              shape=21, colour="grey20", stroke=0.4) +
              scale_fill_manual(values=DIR_PAL) +
              geom_node_text(aes(label=label_show),
                              size=2.8, repel=TRUE, max.overlaps=40,
                              colour="black") +
              labs(title="STRING interaction network: M vs GM significant proteins",
                   subtitle=sprintf("STRING v12, mouse, confidence >= 700;  %d nodes, %d edges",
                                    igraph::vcount(g_lcc), igraph::ecount(g_lcc)),
                   size="degree", fill="direction",
                   edge_width="STRING score") +
              theme_graph(base_family="sans")
            save_svg("10_string_network_M_vs_GM.svg", p10, 11, 10)
          }
        }
      }
    }
  }
}

# ---- STAGE 4-5   LPS RESPONSE + INTERACTION ----

n_int_up   <- sum(de_interaction$direction_strict  == "Up")
n_int_dn   <- sum(de_interaction$direction_strict  == "Down")
n_int_up_l <- sum(de_interaction$direction_lenient == "Up")
n_int_dn_l <- sum(de_interaction$direction_lenient == "Down")

## FIGURE 11   Volcano: interaction (M_LPS vs M) - (GM_LPS vs GM)
vol11 <- de_interaction |>
  mutate(
    neg_log10_p = -log10(pmax(sca.P.Value, 1e-300)),
    score       = abs(logFC) * neg_log10_p,
    dir = case_when(
      sca.P.Value < 0.05 & logFC >=  FC_LENIENT ~ "Up",
      sca.P.Value < 0.05 & logFC <= -FC_LENIENT ~ "Down",
      TRUE ~ "NS"
    ),
    dir = factor(dir, levels=c("Up","Down","NS"))
  )

top15_ids <- vol11 |>
  filter(sca.P.Value < 0.05) |>
  arrange(desc(score)) |>
  head(15) |>
  pull(protein_ids)

vol11 <- vol11 |>
  mutate(label = ifelse(protein_ids %in% top15_ids,
                         coalesce(gene_symbol, uniprot_lead),
                         NA_character_))

n_p01 <- sum(vol11$sca.P.Value < 0.01, na.rm=TRUE)

p11 <- ggplot(vol11, aes(logFC, neg_log10_p, colour=dir)) +
  geom_point(alpha=0.65, size=1.4) +
  geom_text_repel(aes(label=label), size=3, max.overlaps=30,
                  box.padding=0.35, show.legend=FALSE, na.rm=TRUE,
                  colour="black", min.segment.length=0) +
  geom_vline(xintercept=c(-FC_LENIENT, FC_LENIENT),
             linetype="dotted", colour="grey60") +
  geom_hline(yintercept=-log10(0.05), linetype="dotted", colour="grey60") +
  geom_hline(yintercept=-log10(0.01), linetype="dashed", colour="grey40") +
  annotate("text", x=Inf, y=-log10(0.01),
           label="raw P = 0.01", hjust=1.05, vjust=-0.5,
           size=3, colour="grey40") +
  scale_colour_manual(values=DIR_PAL, drop=FALSE) +
  labs(title="Interaction contrast: (M_LPS - M) - (GM_LPS - GM)",
       subtitle=sprintf("No proteins pass FDR<0.05 (min adj.P = %.2f); top-15 suggestive hits shown by |logFC|x-log10(P). raw P<0.01: n=%d.",
                        min(de_interaction$sca.adj.pval, na.rm=TRUE), n_p01),
       x=expression(log[2]*" interaction effect"),
       y=expression(-log[10]*"(raw P-value)"),
       colour="raw P<0.05 &\n|FC|>=1.5") +
  theme_bw(11) +
  theme(legend.position="right")
save_svg("11_volcano_interaction.svg", p11, 9, 6.5)

# ---- STAGE 7   PHENOTYPE MARKERS (SURFACE / IMMUNE / SECRETOME SIGNATURES) ----

marker_panel <- bind_rows(
  tibble(category="MHC-II / antigen presentation",
         gene=c("H2-Aa","H2-Ab1","H2-Eb1","H2-DMa","H2-DMb1","H2-DMb2",
                "H2-Oa","H2-Ob","Cd74","Ciita")),
  tibble(category="Antigen processing",
         gene=c("Tap1","Tap2","Tapbp","Psmb8","Psmb9","Psmb10")),
  tibble(category="Phagocytosis / scavenger / endocytic",
         gene=c("Mrc1","Mrc2","Msr1","Marco","Cd36","Cd163","Stab1","Stab2",
                "Mertk","Axl","Tyro3","Tfrc","Scarb1","Scarb2","Lrp1","Lrp2",
                "Cd68","Cd9")),
  tibble(category="Pattern recognition (TLR / CLR)",
         gene=c("Tlr2","Tlr3","Tlr4","Tlr7","Tlr9","Cd14","Ly96","Clec7a","Clec10a")),
  tibble(category="Phagosome / lysosome",
         gene=c("Lyz1","Lyz2","Ctsb","Ctsd","Ctsl","Lamp1","Lamp2")),
  tibble(category="M1 (classical)",
         gene=c("Nos2","Tnf","Il6","Il1b","Il12b","Cxcl9","Cxcl10","Ccl2",
                "Cd80","Cd86","Stat1","Irf5")),
  tibble(category="M2 (alternative)",
         gene=c("Arg1","Il10","Retnla","Chil3","Klf4","Stat6","Irf4","Pparg","Mgl2")),
  tibble(category="Cell-surface ID",
         gene=c("Itgam","Adgre1","Csf1r","Csf2ra","Csf2rb","Fcgr1","Fcgr2b",
                "Fcgr3","Fcgr4","Ly6c1","Ly6c2","Cx3cr1"))
)

marker_panel <- marker_panel |>
  left_join(annot_df |>
              select(gene_symbol, protein_ids, uniprot_lead, entrez_id, protein_names) |>
              filter(!is.na(gene_symbol)),
            by=c("gene"="gene_symbol")) |>
  mutate(detected = !is.na(protein_ids))

per_cond_mean <- function(prot_id) {
  if (is.na(prot_id) || !(prot_id %in% rownames(combat_mat))) {
    return(c(GM=NA_real_, GM_LPS=NA_real_, M=NA_real_, M_LPS=NA_real_))
  }
  v <- combat_mat[prot_id, ]
  c(GM     = mean(v[design_tbl$condition=="GM"]),
    GM_LPS = mean(v[design_tbl$condition=="GM_LPS"]),
    M      = mean(v[design_tbl$condition=="M"]),
    M_LPS  = mean(v[design_tbl$condition=="M_LPS"]))
}
means_mat <- do.call(rbind, lapply(marker_panel$protein_ids, per_cond_mean))
colnames(means_mat) <- paste0("mean_", colnames(means_mat))
marker_panel <- bind_cols(marker_panel, as_tibble(means_mat))

attach_logfc <- function(panel, de_tt, suffix) {
  sub <- de_tt |> select(protein_ids, logFC, sca.adj.pval)
  colnames(sub) <- c("protein_ids", paste0("logFC_", suffix), paste0("fdr_", suffix))
  left_join(panel, sub, by="protein_ids")
}
marker_panel <- marker_panel |>
  attach_logfc(de_M_vs_GM,         "M_vs_GM") |>
  attach_logfc(de_GM_LPS_vs_GM,    "GM_LPS_vs_GM") |>
  attach_logfc(de_M_LPS_vs_M,      "M_LPS_vs_M") |>
  attach_logfc(de_M_LPS_vs_GM_LPS, "M_LPS_vs_GM_LPS")

fc_cols  <- c("logFC_M_vs_GM","logFC_GM_LPS_vs_GM","logFC_M_LPS_vs_M","logFC_M_LPS_vs_GM_LPS")
fdr_cols <- c("fdr_M_vs_GM","fdr_GM_LPS_vs_GM","fdr_M_LPS_vs_M","fdr_M_LPS_vs_GM_LPS")
marker_panel$max_abs_logFC <- pmax(abs(marker_panel[[fc_cols[1]]]),
                                   abs(marker_panel[[fc_cols[2]]]),
                                   abs(marker_panel[[fc_cols[3]]]),
                                   abs(marker_panel[[fc_cols[4]]]),
                                   na.rm=TRUE)
marker_panel$min_fdr <- suppressWarnings(
  pmin(marker_panel[[fdr_cols[1]]],
       marker_panel[[fdr_cols[2]]],
       marker_panel[[fdr_cols[3]]],
       marker_panel[[fdr_cols[4]]],
       na.rm=TRUE))
marker_panel$min_fdr[!is.finite(marker_panel$min_fdr)] <- NA_real_
marker_panel$max_abs_logFC[!is.finite(marker_panel$max_abs_logFC)] <- NA_real_
marker_panel$sig_lenient_any <- with(marker_panel,
  detected & !is.na(min_fdr) & min_fdr < FDR_THRESHOLD & max_abs_logFC >= FC_LENIENT)

write.csv(marker_panel, file.path(OUT_DIR, "tables", "11_marker_panel_detection.csv"),
          row.names=FALSE)

## FIGURE 13   Marker heatmap (significant only)
sig_markers <- marker_panel |> filter(sig_lenient_any)
if (nrow(sig_markers) >= 2) {
  sample_order_h <- design_tbl |>
    arrange(factor(condition, levels=c("GM","GM_LPS","M","M_LPS")), rep) |>
    pull(sample)

  M_markers <- combat_mat[sig_markers$protein_ids, sample_order_h, drop=FALSE]
  z <- t(scale(t(M_markers)))
  rn <- ifelse(!is.na(sig_markers$gene) & sig_markers$gene != "",
               sig_markers$gene, sig_markers$uniprot_lead)
  rownames(z) <- make.unique(rn)

  cats <- unique(marker_panel$category)
  cat_pal <- setNames(colorRampPalette(brewer.pal(8, "Set2"))(length(cats)), cats)

  row_annot <- data.frame(category = sig_markers$category,
                          row.names = rownames(z))
  ann_col   <- ann_df[sample_order_h, , drop=FALSE]
  ann_colors_h <- list(set=SET_PAL, condition=COND_PAL, category=cat_pal)

  svg(file.path(OUT_DIR, "figures", "13_marker_heatmap_phenotype.svg"),
      width=10, height=8)
  print(pheatmap(z,
           color = colorRampPalette(c("#0279EE","white","#FF9400"))(100),
           breaks = seq(-2, 2, length.out=101),
           annotation_col = ann_col,
           annotation_row = row_annot,
           annotation_colors = ann_colors_h,
           cluster_rows = TRUE, cluster_cols = FALSE,
           clustering_distance_rows = "correlation",
           show_rownames = TRUE, fontsize_row = 9, fontsize_col = 9,
           main = sprintf("Significant macrophage markers (%d, FDR<0.05 & |log2FC|>=0.585)",
                          nrow(sig_markers)),
           border_color = NA))
  dev.off()
}

# ---- SUBCELLULAR LOCATION INVENTORY (UniProt + GO fallback) ----

uniprot_cache_file <- file.path(OUT_DIR, "tmp", "uniprot_cache.rds")
cache_env <- new.env(parent=emptyenv())
cache_env$cache <- if (file.exists(uniprot_cache_file))
  readRDS(uniprot_cache_file) else list()

classify_location <- function(locs) {
  if (length(locs) == 0) return(NA_character_)
  ls <- tolower(paste(locs, collapse=" | "))
  if (grepl("plasma membrane|cell membrane", ls)) return("Plasma membrane")
  if (grepl("secreted|extracellular", ls))         return("Secreted")
  if (grepl("membrane", ls))                       return("Cell membrane (other)")
  return("Intracellular")
}

go_fallback <- function(entrez_id) {
  if (is.na(entrez_id)) return(NA_character_)
  go_terms <- tryCatch({
    AnnotationDbi::select(org.Mm.eg.db, keys=as.character(entrez_id),
                          columns="GOALL", keytype="ENTREZID")$GOALL
  }, error=function(e) NULL)
  if (is.null(go_terms)) return(NA_character_)
  if (any(c("GO:0005886","GO:0009986","GO:0009897","GO:0004888") %in% go_terms))
    return("Plasma membrane")
  if ("GO:0005615" %in% go_terms || "GO:0005576" %in% go_terms) return("Secreted")
  if ("GO:0016020" %in% go_terms) return("Cell membrane (other)")
  return("Intracellular")
}

query_uniprot_batch <- function(accessions) {
  accs <- unique(na.omit(gsub("-[0-9]+$", "", accessions)))
  if (length(accs) == 0) return(list())
  query <- paste(sprintf("accession:%s", accs), collapse = " OR ")
  url   <- paste0(
    "https://rest.uniprot.org/uniprotkb/search?",
    "query=", URLencode(query, reserved = TRUE),
    "&fields=accession,cc_subcellular_location",
    "&format=json&size=", min(length(accs), 500)   # B4: capped at 500
  )
  resp <- tryCatch(GET(url, timeout(60)), error = function(e) NULL)
  if (is.null(resp) || status_code(resp) != 200) return(list())
  js <- tryCatch(content(resp, "parsed", encoding = "UTF-8"), error = function(e) NULL)
  if (is.null(js) || is.null(js$results)) return(list())
  out <- list()
  for (entry in js$results) {
    acc  <- entry$primaryAccession
    locs <- character()
    for (c1 in entry$comments) {
      if (!is.null(c1$commentType) && c1$commentType == "SUBCELLULAR LOCATION")
        for (loc in c1$subcellularLocations)
          if (!is.null(loc$location$value)) locs <- c(locs, loc$location$value)
    }
    out[[acc]] <- list(locations = locs, status = "ok")
  }
  out
}

dir.create(file.path(OUT_DIR, "tmp"), showWarnings=FALSE, recursive=TRUE)
inventory <- annot_df |>
  mutate(location_uniprot = NA_character_, location_evidence = NA_character_,
         location = NA_character_)
n_ok <- 0; n_fail <- 0
batch_size <- 100
accs_all   <- gsub("-[0-9]+$", "", inventory$uniprot_lead)
batches    <- split(seq_len(nrow(inventory)),
                    ceiling(seq_len(nrow(inventory)) / batch_size))
for (b in seq_along(batches)) {
  idx  <- batches[[b]]
  uncached_idx <- idx[!accs_all[idx] %in% names(cache_env$cache)]
  if (length(uncached_idx) > 0) {
    batch_result <- query_uniprot_batch(accs_all[uncached_idx])
    for (acc in names(batch_result))
      cache_env$cache[[acc]] <- batch_result[[acc]]
    Sys.sleep(0.5)
  }
  for (i in idx) {
    acc <- accs_all[i]
    if (is.na(acc) || acc == "") next
    obj <- cache_env$cache[[acc]]
    if (is.null(obj) || obj$status != "ok") {
      n_fail <- n_fail + 1
    } else {
      n_ok <- n_ok + 1
      inventory$location_uniprot[i] <- paste(obj$locations, collapse = " | ")
      inventory$location[i]         <- classify_location(obj$locations)
      inventory$location_evidence[i] <- "UniProt"
    }
  }
  if (b %% 5 == 0) saveRDS(cache_env$cache, uniprot_cache_file)
}
saveRDS(cache_env$cache, uniprot_cache_file)
need_fallback <- is.na(inventory$location)
for (i in which(need_fallback)) {
  loc <- go_fallback(inventory$entrez_id[i])
  if (!is.na(loc)) {
    inventory$location[i]         <- loc
    inventory$location_evidence[i] <- "GO"
  }
}
inventory$location[is.na(inventory$location)] <- "Unknown"
inventory$location_evidence[is.na(inventory$location_evidence)] <- "Unresolved"

inventory <- inventory |>
  left_join(de_M_vs_GM         |> select(protein_ids, logFC_M_vs_GM=logFC,
                                          fdr_M_vs_GM=sca.adj.pval),
            by="protein_ids") |>
  left_join(de_GM_LPS_vs_GM    |> select(protein_ids, logFC_GM_LPS_vs_GM=logFC,
                                          fdr_GM_LPS_vs_GM=sca.adj.pval),
            by="protein_ids") |>
  left_join(de_M_LPS_vs_M      |> select(protein_ids, logFC_M_LPS_vs_M=logFC,
                                          fdr_M_LPS_vs_M=sca.adj.pval),
            by="protein_ids") |>
  left_join(de_M_LPS_vs_GM_LPS |> select(protein_ids, logFC_M_LPS_vs_GM_LPS=logFC,
                                          fdr_M_LPS_vs_GM_LPS=sca.adj.pval),
            by="protein_ids")
inv_fc <- c("logFC_M_vs_GM","logFC_GM_LPS_vs_GM",
            "logFC_M_LPS_vs_M","logFC_M_LPS_vs_GM_LPS")
inventory$max_abs_logFC <- suppressWarnings(
  apply(inventory[, inv_fc], 1, function(r) max(abs(r), na.rm=TRUE)))
inventory$max_abs_logFC[!is.finite(inventory$max_abs_logFC)] <- NA_real_
mean_per_cond <- t(sapply(inventory$protein_ids, function(pid) {
  if (!(pid %in% rownames(combat_mat))) return(rep(NA_real_, 4))
  v <- combat_mat[pid, ]
  c(mean(v[design_tbl$condition=="GM"]),
    mean(v[design_tbl$condition=="GM_LPS"]),
    mean(v[design_tbl$condition=="M"]),
    mean(v[design_tbl$condition=="M_LPS"]))
}))
colnames(mean_per_cond) <- c("mean_GM","mean_GM_LPS","mean_M","mean_M_LPS")
inventory <- bind_cols(inventory, as_tibble(mean_per_cond))

write.csv(inventory, file.path(OUT_DIR, "tables", "12_surface_receptor_inventory.csv"),
          row.names=FALSE)

## FIGURE 14   Surface receptor heatmap — presentation-ready
pm_sig <- inventory |>
  filter(location == "Plasma membrane", is.finite(max_abs_logFC)) |>
  filter(
    (fdr_M_vs_GM         < FDR_THRESHOLD & abs(logFC_M_vs_GM)         >= FC_LENIENT) |
    (fdr_GM_LPS_vs_GM    < FDR_THRESHOLD & abs(logFC_GM_LPS_vs_GM)    >= FC_LENIENT) |
    (fdr_M_LPS_vs_M      < FDR_THRESHOLD & abs(logFC_M_LPS_vs_M)      >= FC_LENIENT) |
    (fdr_M_LPS_vs_GM_LPS < FDR_THRESHOLD & abs(logFC_M_LPS_vs_GM_LPS) >= FC_LENIENT)
  )

# Fallback to FDR<0.1 if <10 proteins pass
if (nrow(pm_sig) < 10) {
  message("  Fig 14: fewer than 10 PM proteins at FDR<0.05; relaxing to FDR<0.1")
  pm_sig <- inventory |>
    filter(location == "Plasma membrane", is.finite(max_abs_logFC)) |>
    filter(
      (fdr_M_vs_GM         < 0.1 & abs(logFC_M_vs_GM)         >= FC_LENIENT) |
      (fdr_GM_LPS_vs_GM    < 0.1 & abs(logFC_GM_LPS_vs_GM)    >= FC_LENIENT) |
      (fdr_M_LPS_vs_M      < 0.1 & abs(logFC_M_LPS_vs_M)      >= FC_LENIENT) |
      (fdr_M_LPS_vs_GM_LPS < 0.1 & abs(logFC_M_LPS_vs_GM_LPS) >= FC_LENIENT)
    )
}

pm <- pm_sig |> arrange(desc(max_abs_logFC)) |> head(20)

if (nrow(pm) >= 2) {
  sample_order_h <- design_tbl |>
    arrange(factor(condition, levels=c("GM","GM_LPS","M","M_LPS")), rep) |>
    pull(sample)
  M_pm <- combat_mat[pm$protein_ids, sample_order_h, drop=FALSE]
  z_pm <- t(scale(t(M_pm)))
  rn_pm <- ifelse(!is.na(pm$gene_symbol) & pm$gene_symbol != "",
                  pm$gene_symbol, pm$uniprot_lead)
  rownames(z_pm) <- make.unique(rn_pm)

  row_annot_pm <- data.frame(evidence = pm$location_evidence,
                             row.names = rownames(z_pm))
  ann_col_pm   <- ann_df[sample_order_h, , drop=FALSE]
  ann_colors_pm <- list(set=SET_PAL, condition=COND_PAL, evidence=LOC_PAL)

  svg(file.path(OUT_DIR, "figures", "14_surface_receptor_heatmap.svg"),
      width=10, height=8)
  print(pheatmap(z_pm,
           color = colorRampPalette(c("#0279EE","white","#FF9400"))(100),
           breaks = seq(-2, 2, length.out=101),
           annotation_col = ann_col_pm,
           annotation_row = row_annot_pm,
           annotation_colors = ann_colors_pm,
           cluster_rows = TRUE, cluster_cols = FALSE,
           clustering_distance_rows = "correlation",
           show_rownames = TRUE, fontsize_row = 10, fontsize_col = 9,
           main = sprintf("Top %d significant plasma-membrane receptors (FDR<0.05, |log2FC|>=0.585)",
                          nrow(pm)),
           border_color = NA))
  dev.off()
}

# ---- SECRETOME ----

curated_secretome <- tibble::tribble(
  ~gene,        ~secretome_class,
  "Il1a",       "Cytokine",
  "Il1b",       "Cytokine",
  "Il6",        "Cytokine",
  "Il10",       "Cytokine",
  "Il12a",      "Cytokine",
  "Il12b",      "Cytokine",
  "Il18",       "Cytokine",
  "Tnf",        "Cytokine",
  "Ifnb1",      "Cytokine",
  "Ccl2",       "Chemokine",
  "Ccl3",       "Chemokine",
  "Ccl4",       "Chemokine",
  "Ccl5",       "Chemokine",
  "Ccl12",      "Chemokine",
  "Ccl17",      "Chemokine",
  "Ccl22",      "Chemokine",
  "Cxcl1",      "Chemokine",
  "Cxcl2",      "Chemokine",
  "Cxcl9",      "Chemokine",
  "Cxcl10",     "Chemokine",
  "Cxcl11",     "Chemokine",
  "Mmp2",       "ECM/remodelling",
  "Mmp9",       "ECM/remodelling",
  "Mmp12",      "ECM/remodelling",
  "Mmp13",      "ECM/remodelling",
  "Timp1",      "ECM/remodelling",
  "Timp2",      "ECM/remodelling",
  "Serpinb2",   "ECM/remodelling",
  "Serpine1",   "ECM/remodelling",
  "Serping1",   "ECM/remodelling",
  "Vegfa",      "Growth factor",
  "Tgfb1",      "Growth factor",
  "Csf1",       "Growth factor",
  "Csf2",       "Growth factor"
)

secretome <- inventory |>
  filter(location == "Secreted" | gene_symbol %in% curated_secretome$gene) |>
  left_join(curated_secretome, by=c("gene_symbol"="gene")) |>
  mutate(secretome_class = ifelse(is.na(secretome_class), "Other", secretome_class)) |>
  filter(!duplicated(protein_ids))

write.csv(secretome,
          file.path(OUT_DIR, "tables", "13_secretome_inventory.csv"),
          row.names=FALSE)

## FIGURE 15   Secretome heatmap 

is_sig_any <- function(df) {
  (df$fdr_M_vs_GM         < FDR_THRESHOLD & abs(df$logFC_M_vs_GM)         >= FC_LENIENT) |
  (df$fdr_GM_LPS_vs_GM    < FDR_THRESHOLD & abs(df$logFC_GM_LPS_vs_GM)    >= FC_LENIENT) |
  (df$fdr_M_LPS_vs_M      < FDR_THRESHOLD & abs(df$logFC_M_LPS_vs_M)      >= FC_LENIENT) |
  (df$fdr_M_LPS_vs_GM_LPS < FDR_THRESHOLD & abs(df$logFC_M_LPS_vs_GM_LPS) >= FC_LENIENT)
}

sec_curated_sig <- secretome |>
  filter(gene_symbol %in% curated_secretome$gene,
         is.finite(max_abs_logFC),
         is_sig_any(pick(everything()))) |>
  arrange(desc(max_abs_logFC))

if (nrow(sec_curated_sig) >= 15) {
  sec_show <- head(sec_curated_sig, 20)
} else {
  extra <- secretome |>
    filter(!protein_ids %in% sec_curated_sig$protein_ids,
           location == "Secreted",
           is.finite(max_abs_logFC)) |>
    arrange(desc(max_abs_logFC)) |>
    head(20 - nrow(sec_curated_sig))
  sec_show <- bind_rows(sec_curated_sig, extra) |>
    arrange(desc(max_abs_logFC)) |>
    head(20)
}

if (nrow(sec_show) >= 2) {
  sample_order_h <- design_tbl |>
    arrange(factor(condition, levels=c("GM","GM_LPS","M","M_LPS")), rep) |>
    pull(sample)
  M_sec <- combat_mat[sec_show$protein_ids, sample_order_h, drop=FALSE]
  z_sec <- t(scale(t(M_sec)))
  rn_sec <- ifelse(!is.na(sec_show$gene_symbol) & sec_show$gene_symbol != "",
                   sec_show$gene_symbol, sec_show$uniprot_lead)
  rownames(z_sec) <- make.unique(rn_sec)

  row_annot_sec <- data.frame(class = sec_show$secretome_class,
                              row.names = rownames(z_sec))
  ann_col_sec   <- ann_df[sample_order_h, , drop=FALSE]
  ann_colors_sec <- list(set=SET_PAL, condition=COND_PAL, class=SECRET_PAL)

  svg(file.path(OUT_DIR, "figures", "15_secretome_heatmap.svg"),
      width=10, height=8)
  print(pheatmap(z_sec,
           color = colorRampPalette(c("#0279EE","white","#FF9400"))(100),
           breaks = seq(-2, 2, length.out=101),
           annotation_col = ann_col_sec,
           annotation_row = row_annot_sec,
           annotation_colors = ann_colors_sec,
           cluster_rows = TRUE, cluster_cols = FALSE,
           clustering_distance_rows = "correlation",
           show_rownames = TRUE, fontsize_row = 10, fontsize_col = 9,
           main = sprintf("Secretome: curated + significant proteins (%d, FDR<0.05 & |log2FC|>=0.585)",
                          nrow(sec_show)),
           border_color = NA))
  dev.off()
}

# ---- FIGURE 16   GSVA pathway-activity heatmap ----

get_kegg_entrez_cached <- NULL  

kegg_cache_file <- file.path(OUT_DIR, "tmp", "kegg_mmu_paths.rds")
if (file.exists(kegg_cache_file)) {
  kegg_path_df <- readRDS(kegg_cache_file)
} else {
  kegg_path_df <- tryCatch(clusterProfiler::download_KEGG("mmu")$KEGGPATHID2EXTID,
                            error=function(e) NULL)
  if (!is.null(kegg_path_df)) saveRDS(kegg_path_df, kegg_cache_file)
}
get_kegg_entrez_cached <- function(pid) {
  if (is.null(kegg_path_df)) return(character(0))
  ids <- kegg_path_df$to[kegg_path_df$from == pid]
  unique(gsub("^mmu:", "", ids))
}

get_reactome_entrez <- function(rid) {
  tryCatch({
    res <- AnnotationDbi::select(reactome.db::reactome.db,
                                  keys=rid, keytype="PATHID",
                                  columns="ENTREZID")
    unique(res$ENTREZID)
  }, error=function(e) character(0))
}

sym_to_entrez <- function(symbols) {
  if (length(symbols) == 0) return(character(0))
  res <- tryCatch({
    AnnotationDbi::select(org.Mm.eg.db, keys=symbols, keytype="SYMBOL",
                          columns="ENTREZID")
  }, error=function(e) NULL)
  if (is.null(res)) return(character(0))
  unique(na.omit(res$ENTREZID))
}

pathway_defs <- list(
  Glycolysis              = list(source="KEGG",     id="mmu00010", curated=NULL),
  OXPHOS                  = list(source="KEGG",     id="mmu00190", curated=NULL),
  `NF-kB signaling`       = list(source="KEGG",     id="mmu04064", curated=NULL),
  `IFN-alpha/beta`        = list(source="Reactome", id="R-MMU-909733", curated=NULL),
  `MHC-II presentation`   = list(source="KEGG",     id="mmu04612", curated=NULL),
  Phagocytosis            = list(source="Reactome", id="R-MMU-2168880", curated=NULL),
  Lysosome                = list(source="KEGG",     id="mmu04142", curated=NULL),
  `ECM remodelling`       = list(source="Reactome", id="R-MMU-1474228", curated=NULL),
  `Fatty-acid beta-ox`    = list(source="KEGG",     id="mmu00071", curated=NULL),
  `M1 cytokines`          = list(source="Curated",  id=NA,
                                  curated=c("Tnf","Il6","Il1b","Il12b","Nos2","Cxcl9","Cxcl10","Ccl2")),
  `M2 cytokines`          = list(source="Curated",  id=NA,
                                  curated=c("Il10","Arg1","Retnla","Chil3","Mgl2","Tgfb1"))
)

gs_list <- list()
for (nm in names(pathway_defs)) {
  d <- pathway_defs[[nm]]
  if (d$source == "KEGG")          g <- get_kegg_entrez_cached(d$id)
  else if (d$source == "Reactome") g <- get_reactome_entrez(d$id)
  else if (d$source == "Curated")  g <- sym_to_entrez(d$curated)
  else                              g <- character(0)
  gs_list[[nm]] <- as.character(g)
  message(sprintf("    %-25s  size=%d  source=%s",
                  nm, length(gs_list[[nm]]), d$source))
}

prot_entrez <- annot_df |>
  filter(!is.na(entrez_id), protein_ids %in% rownames(combat_mat)) |>
  select(protein_ids, entrez_id)
detected_entrez <- unique(prot_entrez$entrez_id)
gs_list <- lapply(gs_list, function(g) intersect(g, detected_entrez))
gs_keep <- sapply(gs_list, length) >= 5
if (any(!gs_keep)) {
  message(sprintf("  Dropping %d pathways with <5 detected members: %s",
                  sum(!gs_keep), paste(names(gs_list)[!gs_keep], collapse=", ")))
}
gs_list <- gs_list[gs_keep]

ent2pid <- prot_entrez |>
  group_by(entrez_id) |>
  summarise(
    protein_ids = {
      pids <- pull(pick(protein_ids), protein_ids)
      pids[which.max(rowMeans(combat_mat[pids, , drop=FALSE], na.rm=TRUE))]
    },
    .groups="drop"
  )
ent_mat <- combat_mat[ent2pid$protein_ids, , drop=FALSE]
rownames(ent_mat) <- ent2pid$entrez_id

gsva_scores <- tryCatch({
  if (utils::packageVersion("GSVA") >= "1.50.0") {
    par <- GSVA::gsvaParam(exprData=ent_mat, geneSets=gs_list, kcdf="Gaussian")
    GSVA::gsva(par, verbose=FALSE)
  } else {
    GSVA::gsva(ent_mat, gs_list, method="gsva", kcdf="Gaussian", verbose=FALSE)
  }
}, error=function(e) {
  NULL
})

if (!is.null(gsva_scores)) {
  gsva_long <- as.data.frame(gsva_scores) |>
    tibble::rownames_to_column("pathway") |>
    tidyr::pivot_longer(-pathway, names_to="sample", values_to="gsva_score") |>
    left_join(design_tbl |> select(sample, condition, set, rep), by="sample")
  write.csv(gsva_long, file.path(OUT_DIR, "tables", "14_pathway_gsva_scores.csv"),
            row.names=FALSE)

  sample_order_g <- design_tbl |>
    arrange(factor(condition, levels=c("GM","GM_LPS","M","M_LPS")), rep) |>
    pull(sample)
  gsva_show    <- gsva_scores[, sample_order_g, drop=FALSE]
  ann_col_g    <- ann_df[sample_order_g, , drop=FALSE]
  ann_colors_g <- list(set=SET_PAL, condition=COND_PAL)

  svg(file.path(OUT_DIR, "figures", "16_pathway_activity_heatmap.svg"),
      width=10, height=max(6, 0.45 * nrow(gsva_show) + 3))
  print(pheatmap(gsva_show,
           color = colorRampPalette(c("#0279EE","white","#FF9400"))(100),
           breaks = seq(-1, 1, length.out=101),
           annotation_col = ann_col_g,
           annotation_colors = ann_colors_g,
           cluster_rows = TRUE, cluster_cols = FALSE,
           show_rownames = TRUE, fontsize_row = 10, fontsize_col = 9,
           main = sprintf("Pathway activity (GSVA on ComBat-adjusted log2 intensity, %d pathways)",
                          nrow(gsva_show)),
           border_color = NA))
  dev.off()

  pathway_le <- tibble(pathway=character(), entrez_id=character(),
                       gene_symbol=character())
  ent_to_sym <- annot_df |> select(entrez_id, gene_symbol) |>
    filter(!is.na(entrez_id)) |> distinct()
  for (pw in names(gs_list)) {
    sym <- ent_to_sym |> filter(entrez_id %in% gs_list[[pw]]) |> pull(gene_symbol)
    if (length(sym) == 0) next
    pathway_le <- bind_rows(pathway_le,
      tibble(pathway = pw,
             entrez_id = gs_list[[pw]][gs_list[[pw]] %in% ent_to_sym$entrez_id],
             gene_symbol = ent_to_sym$gene_symbol[match(
               gs_list[[pw]][gs_list[[pw]] %in% ent_to_sym$entrez_id],
               ent_to_sym$entrez_id)]))
  }
  write.csv(pathway_le, file.path(OUT_DIR, "tables", "15_pathway_leading_edge.csv"),
            row.names=FALSE)
}


# ---- FIGURE 17   Pathway-driver heatmap ----

pathway_drivers <- tibble::tribble(
  ~pathway_group,           ~gene,
  "Glycolysis",             "Eno1",
  "Glycolysis",             "Aldoc",
  "Glycolysis",             "Ldhb",
  "Glycolysis",             "Pfkp",
  "OXPHOS / ETC",           "Sod2",
  "OXPHOS / ETC",           "Ndufa4",
  "NF-kB / TLR regulators", "Irak3",
  "IFN-alpha/beta",         "Ifit3",
  "IFN-alpha/beta",         "Isg15",
  "IFN-alpha/beta",         "Ifitm3",
  "MHC-II presentation",    "H2-Aa",
  "MHC-II presentation",    "H2-Ab1",
  "New M-CSF markers",      "Pik3ap1",
  "New M-CSF markers",      "Tmem176a",
  "New M-CSF markers",      "Tmem176b"
)
pathway_group_levels <- c("Glycolysis", "OXPHOS / ETC",
                          "NF-kB / TLR regulators", "IFN-alpha/beta",
                          "MHC-II presentation", "New M-CSF markers")

driver_tbl <- pathway_drivers |>
  left_join(annot_df |> select(gene_symbol, protein_ids) |>
              filter(!is.na(gene_symbol)) |>
              distinct(gene_symbol, .keep_all = TRUE),
            by = c("gene" = "gene_symbol")) |>
  left_join(de_M_vs_GM |> select(protein_ids,
                                  logFC_M_vs_GM = logFC,
                                  fdr_M_vs_GM   = sca.adj.pval),
            by = "protein_ids") |>
  left_join(de_GM_LPS_vs_GM |> select(protein_ids,
                                       logFC_GM_LPS_vs_GM = logFC,
                                       fdr_GM_LPS_vs_GM   = sca.adj.pval),
            by = "protein_ids") |>
  left_join(de_M_LPS_vs_M |> select(protein_ids,
                                     logFC_M_LPS_vs_M = logFC,
                                     fdr_M_LPS_vs_M   = sca.adj.pval),
            by = "protein_ids") |>
  left_join(de_M_LPS_vs_GM_LPS |> select(protein_ids,
                                          logFC_M_LPS_vs_GM_LPS = logFC,
                                          fdr_M_LPS_vs_GM_LPS   = sca.adj.pval),
            by = "protein_ids")

driver_tbl <- driver_tbl |>
  mutate(
    detected = !is.na(protein_ids) & protein_ids %in% rownames(combat_mat),
    sig_any = detected &
      ((!is.na(fdr_M_vs_GM)         & fdr_M_vs_GM         < FDR_THRESHOLD &
        !is.na(logFC_M_vs_GM)       & abs(logFC_M_vs_GM)         >= FC_LENIENT) |
       (!is.na(fdr_GM_LPS_vs_GM)    & fdr_GM_LPS_vs_GM    < FDR_THRESHOLD &
        !is.na(logFC_GM_LPS_vs_GM)  & abs(logFC_GM_LPS_vs_GM)    >= FC_LENIENT) |
       (!is.na(fdr_M_LPS_vs_M)      & fdr_M_LPS_vs_M      < FDR_THRESHOLD &
        !is.na(logFC_M_LPS_vs_M)    & abs(logFC_M_LPS_vs_M)      >= FC_LENIENT) |
       (!is.na(fdr_M_LPS_vs_GM_LPS) & fdr_M_LPS_vs_GM_LPS < FDR_THRESHOLD &
        !is.na(logFC_M_LPS_vs_GM_LPS)&abs(logFC_M_LPS_vs_GM_LPS) >= FC_LENIENT)),
    reason = case_when(
      !detected ~ "not_detected",
      !sig_any  ~ "not_significant",
      TRUE      ~ "kept"
    )
  )

# Log exclusions
excluded_drivers <- driver_tbl |> filter(reason != "kept")
if (nrow(excluded_drivers) > 0) {
  message(sprintf("  Fig 17: excluded %d protein(s) — %s",
                  nrow(excluded_drivers),
                  paste(sprintf("%s (%s)", excluded_drivers$gene,
                                excluded_drivers$reason),
                        collapse = ", ")))
}
write.csv(
  driver_tbl |> filter(reason != "kept") |>
    select(gene, pathway_group, reason,
           logFC_M_vs_GM,         fdr_M_vs_GM,
           logFC_GM_LPS_vs_GM,    fdr_GM_LPS_vs_GM,
           logFC_M_LPS_vs_M,      fdr_M_LPS_vs_M,
           logFC_M_LPS_vs_GM_LPS, fdr_M_LPS_vs_GM_LPS),
  file.path(OUT_DIR, "tables", "17_pathway_drivers_excluded.csv"),
  row.names = FALSE
)

driver_keep <- driver_tbl |> filter(reason == "kept")

if (nrow(driver_keep) >= 2) {
  driver_keep <- driver_keep |>
    mutate(pathway_group = factor(pathway_group, levels = pathway_group_levels)) |>
    arrange(pathway_group, gene)

  sample_order_17 <- design_tbl |>
    arrange(factor(condition, levels = c("GM","GM_LPS","M","M_LPS")), rep) |>
    pull(sample)

  M_drv <- combat_mat[driver_keep$protein_ids, sample_order_17, drop = FALSE]
  rownames(M_drv) <- make.unique(driver_keep$gene)
  z_drv <- t(scale(t(M_drv)))

  ann_col_17 <- ann_df[sample_order_17, , drop = FALSE]
  col_ha <- ComplexHeatmap::HeatmapAnnotation(
    condition = ann_col_17$condition,
    set       = ann_col_17$set,
    col       = list(condition = COND_PAL, set = SET_PAL),
    annotation_name_side = "left",
    annotation_name_gp = grid::gpar(fontsize = 9),
    show_legend = TRUE
  )

  col_fun <- circlize::colorRamp2(
    seq(-2, 2, length.out = 5),
    c("#0279EE", "#7BB6F4", "white", "#FFC97A", "#FF9400")
  )

  ht <- ComplexHeatmap::Heatmap(
    z_drv,
    name             = "z-score",
    col              = col_fun,
    row_split        = driver_keep$pathway_group,
    row_title_rot    = 0,
    row_title_gp     = grid::gpar(fontsize = 10, fontface = "bold"),
    row_gap          = grid::unit(2, "mm"),
    cluster_rows     = TRUE,
    cluster_row_slices = FALSE,
    cluster_columns  = FALSE,
    show_row_names   = TRUE,
    row_names_gp     = grid::gpar(fontsize = 10),
    show_column_names= TRUE,
    column_names_gp  = grid::gpar(fontsize = 9),
    column_names_rot = 45,
    top_annotation   = col_ha,
    border           = FALSE,
    column_title     = sprintf(
      "Pathway-driver proteins (%d significant in >=1 contrast, FDR<%.2f & |log2FC|>=%.2f)",
      nrow(driver_keep), FDR_THRESHOLD, round(FC_LENIENT, 3)),
    column_title_gp  = grid::gpar(fontsize = 11, fontface = "bold"),
    heatmap_legend_param = list(
      title = "z-score (row)",
      at = c(-2, -1, 0, 1, 2),
      legend_height = grid::unit(3.5, "cm")
    )
  )

  svg(file.path(OUT_DIR, "figures", "17_pathway_driver_heatmap.svg"),
      width = 10, height = 8)
  ComplexHeatmap::draw(ht, merge_legend = TRUE,
                       heatmap_legend_side = "right",
                       annotation_legend_side = "right")
  dev.off()
} else {
  message("  Fig 17: <2 proteins passed filter; heatmap skipped.")
}

# ---- FINALIZATION (sessionInfo + parameters.json + run banner) ----

si_path <- file.path(OUT_DIR, "sessionInfo.txt")
writeLines(capture.output(sessionInfo()), si_path)

params <- list(
  pipeline         = "PXD002582 polished reanalysis",
  pipeline_version = "v22-corrected",
  run_date         = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  R_version        = paste0(R.version$major, ".", R.version$minor),
  seed             = SEED,
  thresholds = list(
    FDR           = FDR_THRESHOLD,
    log2FC_strict = FC_STRICT,
    log2FC_lenient= FC_LENIENT,
    min_peptides  = MIN_PEPTIDES
  ),
  imputation = list(
    method        = "MinProb + MEC floor",
    minprob_q     = 0.01,
    minprob_width = 0.3,
    mec_floor     = "min(observed) - 0.5"
  ),
  model = list(
    formula          = "~ 0 + condition + set",
    batch_correction = "ComBat (visualisation only)",
    de_engine        = "limma + DEqMS (spectraCounteBayes, count = razor+unique peptides)"
  ),
  contrasts = c("M_vs_GM","GM_LPS_vs_GM","M_LPS_vs_M","M_LPS_vs_GM_LPS","interaction"),
  inputs = list(
    proteinGroups = normalizePath(PATH_PG,   mustWork=FALSE),
    sdrf          = normalizePath(PATH_SDRF, mustWork=FALSE)
  ),
  data_counts = list(
    n_proteins_raw          = N_PROT_RAW_GLOBAL,
    n_proteins_post_qc      = N_PROT_QC_GLOBAL,
    n_proteins_prefiltered  = N_PROT_PREFILT_GLOBAL,
    n_samples               = ncol(combat_mat)
  ),
  string = list(species=10090, confidence_min=700, version="12.0"),
  gsva   = list(pathways_used = if (exists("gs_list")) names(gs_list) else NULL),
  changes_vs_v1 = list(
    B2 = "deqms_table: removed duplicate columns from everything()",
    B3 = "run_gsea: set.seed(SEED) added before gseGO",
    B4 = "query_uniprot_batch: capped &size= at min(n,500)",
    B5 = "ent2pid: explicit pull(protein_ids) before matrix indexing",
    B6 = "Fig01: added guide_colourbar() with explicit barwidth/barheight; width 9->10",
    B7 = "Fig09: added plot.margin for legend; explicit guide_colourbar(); width 11->13",
    Fig08 = "Redesigned: top 25 by score=(-log10FDR x |logFC|); fontsize 9; direction annotation",
    Fig14 = "Redesigned: sig in >=1 contrast (FDR<0.05), top 20; fontsize 10; 10x8",
    Fig15 = "Redesigned: curated+sig priority, fill to 20; fontsize 10; 10x8",
    Fig17 = "NEW: pathway-driver heatmap from Part IV report - 15 curated proteins in 6 pathway groups (Glycolysis, OXPHOS, NF-kB regulators, IFN, MHC-II, new M-CSF markers); ComplexHeatmap row_split; FDR<0.05 & |log2FC|>=0.585 in >=1 contrast; 10x8 SVG; exclusions logged to tables/17_pathway_drivers_excluded.csv",
    Q1 = "Removed unused rename from dplyr masking",
    Q4 = "Removed unused universe_entrez",
    Q5 = "Removed dead functions kegg_genes() and reactome_genes()"
  )
)
writeLines(jsonlite::toJSON(params, pretty=TRUE, auto_unbox=TRUE, na="null"),
           file.path(OUT_DIR, "parameters.json"))

## Run summary banner
cat("\n", strrep("=", 78), "\n", sep="")
cat("PXD002582 corrected pipeline (v22)   RUN COMPLETE\n")
cat(strrep("=", 78), "\n", sep="")
cat(sprintf("Date:                  %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("Output dir:            %s\n", normalizePath(OUT_DIR, mustWork=FALSE)))
cat(sprintf("Proteins (post-QC):    %d\n", N_PROT_QC_GLOBAL))
cat(sprintf("Proteins (analysed):   %d\n", nrow(combat_mat)))
cat(sprintf("Samples:               %d\n", ncol(combat_mat)))
cat("\nDual-threshold summary (Up | Down):\n")
print(dual_summary, row.names=FALSE)
cat("\nFigures written to:    results/figures/\n")
cat("Tables  written to:    results/tables/\n")
cat("Parameters & sessionInfo: results/parameters.json, results/sessionInfo.txt\n")
cat(strrep("=", 78), "\n\n", sep="")
message("The Dark Side Of The Force Is A Pathway To Many Abilities Some Consider To Be Unnatural")
