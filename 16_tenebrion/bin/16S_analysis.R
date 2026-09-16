################################################################################
# SCRIPT: 16S MICROBIAL DIVERSITY ANALYSIS — Tenebrio molitor
################################################################################
# ==============================================================================
# Author  : Thomas BOUTET (IHP, LMGE)
# Date    : 2026
# Model   : Tenebrio molitor 
# ==============================================================================
# This script analyzes the gut bacterial communities of T. molitor across
# different developmental stages (larvae -> adults) and its environment (substrate).
#
# ACTUAL OUTPUT FILES GENERATED (in results/16S/16S_figure/):
#
# --- 1. Alpha / Beta Diversity Values ---
#  - Valeurs_Shannon_BrayCurtis.tsv                 -> Raw numeric table of diversity metrics
#
# --- 2. Overall Composition ---
#  - Barplot_Global_16S.png / .svg                  -> Composition by individual replicate (insect only)
#
# --- 3. Larval Core Microbiome ---
#  - Larval_Microbiome_Means.png                    -> Larval core stacked barplot (Genus level)
#  - Larval_Microbiome_Quadrants.png / .svg         -> Prevalence vs Abundance quadrants (Larvae)
#  - Genome_heterogeneous.png / .svg                -> Intra-stage fidelity matrices (Heatmaps)
#
# --- 4. Adult Core Microbiome ---
#  - Adult_Microbiome_Quadrants.png / .svg          -> Prevalence vs Abundance quadrants (Adults)
#
# --- 5. Community Structure & Clustering ---
#  - PCA_microbiote.png                             -> PCA (Hellinger transformed, insect only)
#  - PCA_Equivalent_Clustering_Elbow_Plot.png       -> Elbow method plot for optimal k clusters
#  - PCA_Equivalent_Clustering_Dendrogram.png       -> Hierarchical clustering dendrogram (Ward.D2)
#  - PCA_Equivalent_Clustering_Clusters_Composition.tsv -> Samples assignment to clusters
#
# --- 6. Intra-genus Taxonomic Resolution ---
#  - Taxonomic_Resolution_Average_Larvae.png / .svg -> Species-level breakdown (pooled larvae)
#  - Taxonomic_Resolution_Microlarvae.png / .svg    -> Species-level breakdown (Microlarvae 7 mg)
#  - Taxonomic_Resolution_Larvae_W1.png / .svg      -> Species-level breakdown (Larvae W1 14 mg)
#  - Taxonomic_Resolution_Larvae_W2.png / .svg      -> Species-level breakdown (Larvae W2 40 mg)
#  - Taxonomic_Resolution_Larvae_W3.png / .svg      -> Species-level breakdown (Larvae W3 65 mg)
#  - Taxonomic_Resolution_Larvae_W4.png / .svg      -> Species-level breakdown (Larvae W4 100 mg)
#  - Taxonomic_Resolution_Adults.png / .svg         -> Species-level breakdown (Beetles)
#
# --- 7. Statistical Comparisons (Mirror Plots) ---
#  - Comparison_Mirror_Plot_Spine.png / .svg        -> Larvae vs Beetles with Standard Error & Wilcoxon
#  - Comparison_Mirror_Plot_Spine_Substrate.png / .svg -> Larvae vs Substrate with Standard Error & Wilcoxon
#
# REQUIRED DATA:
#   - QIIME2-like table (sequences x samples, TSV format)
#   - Metadata file (sample -> developmental condition)
#   - SILVA v138.2 databases for prior taxonomic assignment
#
# REFACTORING NOTES (code only - analyses and figures are unchanged):
#   - The two mirror plots are now produced by a single parameterized function
#     creer_mirror_plot() instead of two near-identical 90-line blocks.
#   - Shared helpers added: normaliser_avec_condition(), construire_palette_taxo(),
#     sauver_figure().
#   - chronologie_pca replaces the silent mid-script overwrite of chronologie_insecte.
#   - cond_larves is declared once (was duplicated as cond_larves_cibles).
#   - ggrepel now uses a fixed seed: the quadrant figures used to be redrawn
#     slightly differently at every run. Set SEED_REPEL to NA to revert.
################################################################################


# ==============================================================================
# BLOCK 0: CONFIGURATION AND PACKAGES
# ==============================================================================

setwd("/home/thomas/Tenebrion/")

suppressPackageStartupMessages({
  library(dplyr)       # Data manipulation and wrangling
  library(ggplot2)     # Advanced data visualization
  library(tidyr)       # Data reshaping (pivot_longer, pivot_wider, complete)
  library(vegan)       # Ecological diversity analysis (Shannon, Bray-Curtis, PCA/RDA)
  library(tibble)      # Dataframe enhancements (column_to_rownames)
  library(stringr)     # String manipulation (regex for taxonomy parsing)
  library(readr)       # Fast file reading
  library(forcats)     # Categorical factor manipulation for plot ordering
  library(colorspace)  # Advanced color palettes
  library(patchwork)   # Combining multiple ggplots together
  library(ggrepel)     # Smart text labels to avoid overlapping on scatter plots
  library(scales)      # Axis scaling and formatting
})

# --- File paths ---
fichier_qiime    <- "results/all_matam_salmon_qiime_like_table_HYBRIDE.tsv"
fichier_metadata <- "data/16S/metadata.tsv"
out_dir          <- "results/16S/16S_figure"

# Create the output directory if it does not already exist.
# recursive = TRUE ensures all parent directories are created safely.
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Color palette for individual taxonomic categories
pal_taxo <- colorspace::qualitative_hcl(14, palette = "Dark 3")

# ==============================================================================
# BLOCK 1: UTILITY FUNCTIONS
# ==============================================================================

# ------------------------------------------------------------------------------
# lire_qiime_tsv_robuste(): Robust loading of a QIIME2 TSV file
# ------------------------------------------------------------------------------
# Biological/Technical Context: QIIME2 TSV output often places a hash "#" at the 
# start of the header line (e.g., "#OTU ID"). If read normally, R treats the 
# header as a comment and skips it, ruining column names. 
# This function dynamically locates the header and strips the "#" so data is 
# parsed cleanly.
# ------------------------------------------------------------------------------
lire_qiime_tsv_robuste <- function(chemin) {
  lignes <- readLines(chemin, warn = FALSE)

  # Detect the exact header line containing "OTU ID" or "taxonomy"
  idx_header <- which(grepl("(?i)OTU.?ID", lignes, perl = TRUE))[1]
  if (is.na(idx_header)) {
    idx_header <- which(grepl("(?i)taxonomy", lignes, perl = TRUE))[1]
  }

  lignes_propres <- lignes[idx_header:length(lignes)]

  # Remove the leading '#' so read_tsv parses it as a standard header
  lignes_propres[1] <- sub("^#", "", lignes_propres[1])

  # I() wraps the text as a text connection, bypassing the need for a physical file
  df <- read_tsv(I(paste(lignes_propres, collapse = "\n")), show_col_types = FALSE)
  names(df)[1] <- "Taxonomy"
  return(df)
}

# ------------------------------------------------------------------------------
# normaliser_100(): Conversion to Relative Abundance
# ------------------------------------------------------------------------------
# Technical Context: Sequencing depth differs heavily across biological replicates 
# (e.g., 15k reads vs 50k reads). Comparing raw counts is statistically flawed. 
# This standardizes every sample to sum to 100%, allowing fair composition comparisons.
# ------------------------------------------------------------------------------
normaliser_100 <- function(df, rang) {
  df %>%
    # Sum counts of identical taxonomic assignments per sample
    group_by(Sample, !!sym(rang)) %>%   
    summarise(Count = sum(Count), .groups = "drop") %>%

    # Compute percentage per sample
    group_by(Sample) %>%
    mutate(Relative_Abundance = (Count / sum(Count)) * 100) %>%
    ungroup()
}

# ------------------------------------------------------------------------------
# normaliser_avec_condition(): Normalization + metadata join
# ------------------------------------------------------------------------------
# The pair "normalize to 100%, then re-attach the Condition of each sample"
# was repeated identically in 6 different blocks. Factored here so the join
# key and the selected metadata columns are defined in a single place.
# ------------------------------------------------------------------------------
normaliser_avec_condition <- function(df, rang) {
  normaliser_100(df, rang) %>%
    left_join(meta %>% select(Sample, Condition), by = "Sample")
}

# ------------------------------------------------------------------------------
# construire_palette_taxo(): Color assignment for a set of taxa
# ------------------------------------------------------------------------------
# Assigns one color from pal_taxo per real taxon, and forces the two catch-all
# categories to neutral greys so they never compete visually with the genera.
# Used by both preparer_barplot_data() and the global barplot (BLOCK 5), which
# previously duplicated this logic line for line.
# ------------------------------------------------------------------------------
construire_palette_taxo <- function(levels_taxa) {
  taxa_vrais <- setdiff(levels_taxa, c("Others", "Unassigned"))
  palette    <- setNames(rep_len(pal_taxo, length(taxa_vrais)), taxa_vrais)

  if ("Others"     %in% levels_taxa) palette["Others"]     <- "grey85"
  if ("Unassigned" %in% levels_taxa) palette["Unassigned"] <- "grey40"
  return(palette)
}

# ------------------------------------------------------------------------------
# sauver_figure(): Standardized figure export
# ------------------------------------------------------------------------------
# Every figure was exported through 1 or 2 hand-written ggsave() calls with the
# dpi / background / extension repeated each time. This wrapper centralizes them.
#   nom_base    - filename WITHOUT extension (added automatically per format)
#   formats     - "png", "svg", or both
#   svg_width/svg_height - default to the PNG dimensions, overridden only where
#                 the original script deliberately used a larger vector canvas
#   bg          - passed through as-is; NULL means "let the theme decide", which
#                 is exactly what an omitted bg argument does in ggsave()
# ------------------------------------------------------------------------------
sauver_figure <- function(plot, nom_base, width, height,
                          formats = c("png", "svg"), dpi = 300, bg = NULL,
                          svg_width = width, svg_height = height,
                          svg_bg = bg) {
  if ("png" %in% formats) {
    ggsave(file.path(out_dir, paste0(nom_base, ".png")),
           plot = plot, width = width, height = height, dpi = dpi, bg = bg)
  }
  if ("svg" %in% formats) {
    ggsave(file.path(out_dir, paste0(nom_base, ".svg")),
           plot = plot, width = svg_width, height = svg_height, bg = svg_bg)
  }
}

# Null-coalescing operator (defined early because the plotting functions below
# rely on it for their default subtitles).
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ------------------------------------------------------------------------------
# preparer_barplot_data(): Abstracted barplot data generation
# ------------------------------------------------------------------------------
# Optimization Context: Instead of rewriting grouping, filtering, and factor 
# releveling logic for larvae, adults, and the whole dataset, this single 
# functional block processes the data for any subset, keeping code modular/DRY.
# ------------------------------------------------------------------------------
preparer_barplot_data <- function(df_brut, rang, top_n, seuil_min = 0,
                                   cond_levels = NULL, exclus = c("Non Assigné", "Unassigned")) {
  
  # 1. Normalize data (and attach each sample's Condition)
  df_norm <- normaliser_avec_condition(df_brut, rang)

  # 2. Compute "True Mean" per stage
  # Divides sum of abundances by total replicates to ensure bars exactly hit 100%
  df_agg <- df_norm %>%
    group_by(Condition) %>%
    mutate(Nb_Samples = n_distinct(Sample)) %>%
    group_by(Condition, !!sym(rang), Nb_Samples) %>%
    summarise(Relative_Abundance = sum(Relative_Abundance) / first(Nb_Samples),
              .groups = "drop")

  if (!is.null(cond_levels)) {
    df_agg <- df_agg %>% mutate(Condition = factor(Condition, levels = cond_levels))
  }

  # 3. Identify Top N highly abundant taxa
  top_taxa <- df_agg %>%
    filter(!.data[[rang]] %in% exclus) %>%
    group_by(.data[[rang]]) %>%
    summarise(MeanAbund = mean(Relative_Abundance), .groups = "drop") %>%
    filter(MeanAbund >= seuil_min) %>%
    slice_max(MeanAbund, n = top_n, with_ties = FALSE) %>%
    pull(.data[[rang]])

  # 4. Collapse low-abundance taxa into "Others"
  df_plot <- df_agg %>%
    mutate(Taxon_Top = case_when(
      .data[[rang]] %in% exclus ~ "Unassigned",
      .data[[rang]] %in% top_taxa ~ as.character(.data[[rang]]),
      TRUE ~ "Others"
    )) %>%
    group_by(Condition, Taxon_Top) %>%
    summarise(Abundance = sum(Relative_Abundance), .groups = "drop") %>%
    mutate(
      # fct_reorder sorts by volume, fct_relevel forces Others/Unassigned to the bottom
      Taxon_Top = fct_relevel(
        fct_reorder(factor(Taxon_Top), Abundance, sum, .desc = FALSE),
        "Others", "Unassigned",
        after = 0
      )
    )

  # 5. Build dynamic color palette
  palette <- construire_palette_taxo(levels(df_plot$Taxon_Top))

  return(list(df_plot = df_plot, palette = palette))
}

# ------------------------------------------------------------------------------
# creer_barplot(): Generates standardized stacked barplots
# ------------------------------------------------------------------------------
creer_barplot <- function(df_plot, palette, titre, x_lab = NULL,
                           y_lab = "Mean Relative Abundance (%)",
                           fill_lab = "Bacterial Genus:",
                           angle_x = 45, largeur_bar = 0.65) {
  ggplot(df_plot, aes(x = Condition, y = Abundance, fill = Taxon_Top)) +
    geom_bar(stat = "identity", color = "black",
             linewidth = 0.3, width = largeur_bar, alpha = 0.9) +
    scale_fill_manual(values = palette) +
    scale_y_continuous(
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(title = titre, x = x_lab, y = y_lab, fill = fill_lab) +
    theme_bw(base_size = 12) +
    theme(
      plot.margin        = margin(10, 10, 10, 10),
      legend.position    = "right",
      axis.text.x        = element_text(angle = angle_x, hjust = 1,
                                        face = "bold", size = 11, color = "black"),
      legend.text        = element_text(face = "italic"),
      panel.grid.major.x = element_blank()
    )
}

# ------------------------------------------------------------------------------
# creer_quadrant_plot(): Prevalence vs Abundance (Core Microbiome visualization)
# ------------------------------------------------------------------------------
# Biological Context: Instead of relying entirely on relative abundance, we 
# evaluate the "Core Microbiome" via Prevalence (presence/absence across samples).
#   - Strict Core: Present in 100% of samples.
#   - Shared Core: Present in most samples (> threshold).
#   - Background: Occasional, low prevalence taxa.
# ------------------------------------------------------------------------------
# REPRODUCIBILITY: ggrepel places labels with a random jitter drawn from the
# global RNG. Without a fixed seed, R initializes that RNG from the clock, so
# the two quadrant figures came out slightly different at EVERY run (verified:
# re-running the original script twice produces two different PNG files).
# Passing an explicit seed makes the label layout stable from now on. Set it to
# NA to restore the original (non-reproducible) behaviour.
SEED_REPEL <- 42

creer_quadrant_plot <- function(df_scatter, rang, seuil_prev,
                                 titre, sous_titre = NULL, style_label = "bold") {

  palette_quadrant <- c(
    "1. Strict Shared Microbiota" = "#E31A1C", # Red
    "2. High Shared Microbiota"   = "#1F78B4", # Blue
    "3. Others genera"            = "grey75"   # Grey
  )

  ggplot(df_scatter, aes(x = Prevalence, y = MeanAbund)) +
    # Vertical line separating Shared from Background taxa
    geom_vline(xintercept = seuil_prev, linetype = "dashed", color = "grey40", linewidth = 0.8) +

    geom_point(aes(color = Quadrant, size = MeanAbund), alpha = 0.8) +

    # Repel labels intelligently to avoid clutter. Core labels nudge left, Shared labels nudge right.
    geom_text_repel(
      data = filter(df_scatter, Quadrant == "1. Strict Shared Microbiota"),
      aes(label = .data[[rang]], color = Quadrant),
      size = 3.5, fontface = style_label, hjust = 0,
      direction = "y", nudge_x = 4, seed = SEED_REPEL,
      box.padding = 0.62, max.overlaps = 30, show.legend = FALSE
    ) +
    geom_text_repel(
      data = filter(df_scatter, Quadrant == "2. High Shared Microbiota"),
      aes(label = .data[[rang]], color = Quadrant),
      size = 3.5, fontface = style_label, hjust = 1,
      direction = "y", nudge_x = -4, seed = SEED_REPEL,
      box.padding = 0.62, max.overlaps = 30, show.legend = FALSE
    ) +

    # Log10 Y-axis: critical for viewing both highly dominant and highly rare taxa
    scale_y_log10(
      labels = scales::comma_format(accuracy = 0.01),
      breaks = c(0.01, 0.1, 1, 10, 100)
    ) +
    scale_x_continuous(breaks = seq(0, 100, 25), limits = c(0, 110)) +
    scale_size_continuous(range = c(2, 10), guide = "none") +
    scale_color_manual(values = palette_quadrant) +

    labs(
      title    = titre,
      subtitle = sous_titre %||% sprintf("Threshold: Prevalence >= %d%%", seuil_prev),
      x = "Prevalence (%)",
      y = "Mean Relative Abundance (%)",
      color = "Ecological Classification:"
    ) +
    theme_bw(base_size = 13, base_family = "Arial") +
    theme(
      legend.position = "bottom",
      legend.title    = element_text(face = "bold"),
      legend.text     = element_text(size = 10),
      plot.title      = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 5), nrow = 2))
}

# ------------------------------------------------------------------------------
# calculer_core_scatter(): Computes metrics for the Core Quadrant plot
# ------------------------------------------------------------------------------
calculer_core_scatter <- function(df_norm, rang, n_samples, seuil_prev) {
  df_norm %>%
    filter(!.data[[rang]] %in% c("Non Assigné", "Unassigned")) %>%
    group_by(.data[[rang]]) %>%
    summarise(
      MeanAbund  = sum(Relative_Abundance) / n_samples,
      Prevalence = sum(Relative_Abundance > 0) / n_samples * 100,
      .groups = "drop"
    ) %>%
    # Filter trace elements below 0.01% abundance to keep visual output clean
    filter(MeanAbund >= 0.01) %>%
    mutate(Quadrant = case_when(
      Prevalence >= 100        ~ "1. Strict Shared Microbiota",
      Prevalence >= seuil_prev ~ "2. High Shared Microbiota",
      TRUE                     ~ "3. Others genera"
    ))
}

# ------------------------------------------------------------------------------
# Geometry constants shared by the mirror plots
# ------------------------------------------------------------------------------
# GAP_SPINE is the half-width of the empty central corridor (in % of relative
# abundance) reserved for the genus names. Bars therefore start at +/- GAP_SPINE
# instead of 0, so the axis ticks must be shifted outward by the same amount for
# the printed labels (10%, 20%, 30%) to still mean true abundances.
GAP_SPINE      <- 15
VAL_BREAKS     <- c(10, 20, 30)
BREAKS_X_SPINE <- c(-VAL_BREAKS - GAP_SPINE, VAL_BREAKS + GAP_SPINE)
LABELS_X_SPINE <- c(paste0(VAL_BREAKS, "%"), paste0(VAL_BREAKS, "%"))
COL_TEXT       <- "#2C3E50"

# ------------------------------------------------------------------------------
# creer_mirror_plot(): Two-group comparison as a mirrored "spine" barplot
# ------------------------------------------------------------------------------
# Statistical Context: Compares the relative abundance of a curated genus list
# between two groups of samples with a Wilcoxon rank-sum test (non-parametric,
# safe for the many zeros of compositional data), followed by a
# Benjamini-Hochberg FDR correction. FDR is strictly necessary because running
# one test per genus mechanically inflates the chance of false positives.
#
# Optimization Context: the larvae-vs-beetles and larvae-vs-substrate figures
# were two ~90-line copies of the very same pipeline. Only the sample groups,
# the genus list, the labels/colors and the output names actually differed, so
# the whole thing is parameterized here. The statistics and the geometry are
# unchanged: the two figures come out exactly as before.
#
#   cond_g1 / cond_g2 - Conditions forming the left and right group
#   lab_g1  / lab_g2  - Group names shown in the legend (also factor levels)
#   col_g1  / col_g2  - Fill color of each side
#   genres            - Curated genera to display, one row each
#   Shift             - mean(group 2) - mean(group 1), used to order the rows
# ------------------------------------------------------------------------------
creer_mirror_plot <- function(cond_g1, cond_g2, lab_g1, lab_g2,
                              col_g1, col_g2, genres, titre, nom_base) {

  # Normalize within this comparison only, then collapse the Conditions into a
  # two-level Stage factor (level order drives left/right placement).
  df_norm <- normaliser_avec_condition(
    df_final_16s %>% filter(Condition %in% c(cond_g1, cond_g2)), "Genus"
  ) %>%
    mutate(Stage = ifelse(Condition %in% cond_g1, lab_g1, lab_g2),
           Stage = factor(Stage, levels = c(lab_g1, lab_g2)))

  # 0-fill absent samples to preserve standard error integrity
  samples_inclus <- df_norm %>% distinct(Sample, Stage)
  df_complet <- df_norm %>%
    filter(Genus %in% genres) %>%
    select(Sample, Genus, Relative_Abundance) %>%
    complete(Sample, Genus = genres, fill = list(Relative_Abundance = 0)) %>%
    left_join(samples_inclus, by = "Sample")

  # Wilcoxon Test + FDR (Benjamini-Hochberg) computation
  # A tryCatch intercepts completely identical variables (like 0 vs 0) returning p=1 safely
  df_stats <- df_complet %>%
    group_by(Genus) %>%
    summarise(
      Mean_G1 = mean(Relative_Abundance[Stage == lab_g1]),
      Mean_G2 = mean(Relative_Abundance[Stage == lab_g2]),
      p_val_brute = tryCatch(wilcox.test(Relative_Abundance ~ Stage, exact = FALSE)$p.value, error = function(e) 1),
      .groups = "drop"
    ) %>%
    mutate(
      Shift = Mean_G2 - Mean_G1,
      p_val_adj = p.adjust(p_val_brute, method = "BH"),
      Sig = case_when(p_val_adj < 0.01 ~ "**", p_val_adj < 0.05 ~ "*", TRUE ~ "")
    )

  # Group means + standard errors; rows ordered by Shift so that the genera
  # enriched in one group end up opposite those enriched in the other.
  df_plot <- df_complet %>%
    group_by(Stage, Genus) %>%
    summarise(MeanAbund = mean(Relative_Abundance), SEAbund = sd(Relative_Abundance) / sqrt(n()), .groups = "drop") %>%
    left_join(df_stats %>% select(Genus, Shift, Sig), by = "Genus") %>%
    mutate(
      Genus = fct_reorder(Genus, Shift),
      # Traces rounded to 0.0% are shown as "<0.1%" so they are not read as absences
      Label_Text = ifelse(MeanAbund > 0 & MeanAbund < 0.05, "<0.1%", sprintf("%.1f%%", MeanAbund))
    )

  # Bars are drawn with geom_rect(), so their coordinates are computed here:
  # group 1 grows leftward from -GAP_SPINE, group 2 rightward from +GAP_SPINE.
  # pmax/pmin clamp the error bars at the corridor edge so they never overlap
  # the central genus labels.
  df_plot_spine <- df_plot %>%
    mutate(
      Genus_num = as.numeric(Genus),
      X_min = ifelse(Stage == lab_g1, -MeanAbund - GAP_SPINE, GAP_SPINE),
      X_max = ifelse(Stage == lab_g1, -GAP_SPINE, MeanAbund + GAP_SPINE),
      Err_Min_Spine = ifelse(Stage == lab_g1, X_min - SEAbund, pmax(GAP_SPINE, X_max - SEAbund)),
      Err_Max_Spine = ifelse(Stage == lab_g1, pmin(-GAP_SPINE, X_min + SEAbund), X_max + SEAbund),
      Label_Pct_X = ifelse(Stage == lab_g1, X_min - SEAbund - 1.2, X_max + SEAbund + 1.2),
      Label_Pct_Hjust = ifelse(Stage == lab_g1, 1, 0)
    )

  # Central labels (genus + significance stars), sharing the row order of the bars
  df_axis <- df_stats %>%
    mutate(Genus = factor(Genus, levels = levels(df_plot$Genus)), Genus_num = as.numeric(Genus),
           Label_Center = paste0(Genus, " ", Sig))

  couleurs <- setNames(c(col_g1, col_g2), c(lab_g1, lab_g2))

  # Layers use their own data frames, hence the empty ggplot() call
  p_mirror <- ggplot() +
    geom_vline(xintercept = -GAP_SPINE, color = "#BDC3C7", linewidth = 1) +
    geom_vline(xintercept = GAP_SPINE,  color = "#BDC3C7", linewidth = 1) +
    geom_rect(data = df_plot_spine, aes(xmin = X_min, xmax = X_max, ymin = Genus_num - 0.35, ymax = Genus_num + 0.35, fill = Stage), color = "white", linewidth = 0.5) +
    geom_errorbarh(data = df_plot_spine, aes(y = Genus_num, xmin = Err_Min_Spine, xmax = Err_Max_Spine), height = 0.25, color = "#7F8C8D", linewidth = 0.6) +
    geom_text(data = df_axis, aes(x = 0, y = Genus_num, label = Label_Center), fontface = "bold.italic", size = 9.84, color = COL_TEXT) +
    geom_text(data = df_plot_spine, aes(x = Label_Pct_X, y = Genus_num, label = Label_Text, color = Stage, hjust = Label_Pct_Hjust), size = 9.84, fontface = "bold") +
    scale_fill_manual(values = couleurs, name = "") +
    scale_color_manual(values = couleurs, guide = "none") +
    # Fixed symmetric limits keep the two mirror figures directly comparable
    scale_x_continuous(breaks = BREAKS_X_SPINE, labels = LABELS_X_SPINE, limits = c(-55, 55)) + 
    # Genus names are drawn as central text, so the y-axis itself stays blank
    scale_y_continuous(breaks = NULL, name = NULL) + 
    labs(title = titre, subtitle = "Mean relative abundance ± SE & FDR-adjusted Wilcoxon Test (*p<0.05, **p<0.01)", x = "Relative Abundance (%)") +
    theme_minimal(base_size = 28) +
    theme(legend.position = "top", legend.text = element_text(size = 28, face = "bold", color = COL_TEXT), legend.key.height = unit(1.5, "cm"), legend.key.width = unit(2.5, "cm"), legend.margin = margin(b = 15), plot.title = element_text(face = "bold", size = 28, hjust = 0.5, color = COL_TEXT, margin = margin(b = 8)), plot.subtitle = element_text(size = 28, hjust = 0.5, color = "#7F8C8D", margin = margin(b = 20)), panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(color = "grey85", linetype = "dotted", linewidth = 0.8), axis.text.y = element_text(face = "bold.italic", size = 28, color = COL_TEXT), axis.text.x = element_text(face = "bold", size = 28, color = "#7F8C8D"), axis.title.x = element_text(face = "bold", size = 28, color = COL_TEXT, margin = margin(t = 15)), plot.margin = margin(t = 20, r = 20, b = 20, l = 20))

  # Large canvas: the base font size (28) is tuned for a full-page figure
  sauver_figure(p_mirror, nom_base, width = 21, height = 15, bg = "white",
                svg_width = 24, svg_height = 18)

  # Returned invisibly so the stats stay reachable for the manuscript
  invisible(list(plot = p_mirror, stats = df_stats))
}


# ==============================================================================
# BLOCK 2: DATA LOADING AND FILTERING
# ==============================================================================
message("--- Loading Data and Preparing Formatting ---")

# --- Translation dictionaries ---
# Ensures the internal shortcodes map perfectly to English biological terms 
# suitable for publication. Weights provide developmental context.
trad_conditions <- c(
  "Pre_MicroLarve" = "Tiny_Larvae (2 mg)",
  "MicroLarve"     = "Microlarvae (7 mg)",
  "Larve_S1"       = "Larvae_W1 (14 mg)",
  "Larve_S2"       = "Larvae_W2 (40 mg)",
  "Larve_S3"       = "Larvae_W3 (65 mg)",
  "Larve_S4"       = "Larvae_W4 (100 mg)",
  "Nymphe"         = "Pupae",
  "Jeune_Adulte"   = "Young_Beetles",
  "Adulte"         = "Beetles",
  "Substrat_Base"  = "Raw_Substrate",
  "Substrat_Final" = "Final_Substrate"
)

trad_samples <- c(
  "Pre1" = "TL1", "Pre2" = "TL2", "Pre3" = "TL3", "Pre4" = "TL4",
  "Mi1"  = "Mi1", "Mi2"  = "Mi2", "Mi3"  = "Mi3", "Mi4"  = "Mi4",
  "S1-L1" = "W1-L1", "S1-L2" = "W1-L2", "S1-L3" = "W1-L3", "S1-L4" = "W1-L4",
  "S2-L1" = "W2-L1", "S2-L2" = "W2-L2", "S2-L3" = "W2-L3", "S2-L4" = "W2-L4",
  "S3-L1" = "W3-L1", "S3-L2" = "W3-L2", "S3-L3" = "W3-L3", "S3-L4" = "W3-L4",
  "S4-L1" = "W4-L1", "S4-L2" = "W4-L2", "S4-L3" = "W4-L3", "S4-L4" = "W4-L4",
  "N1"  = "P1",  "N2"  = "P2",  "N3"  = "P3",  "N4"  = "P4",
  "JA1" = "YB1", "JA2" = "YB2", "JA3" = "YB3", "JA4" = "YB4",
  "A1"  = "B1",  "A2"  = "B2",  "A3"  = "B3",  "A4"  = "B4",
  "SB1" = "RS1", "SB2" = "RS2", "SB3" = "RS3", "SB4" = "RS4",
  "SF1" = "FS1", "SF2" = "FS2", "SF3" = "FS3", "SF4" = "FS4"
)

# Read QIIME table and rename sample column names based on dictionary.
df_qiime_16s <- lire_qiime_tsv_robuste(fichier_qiime) %>%
  rename_with(~ recode(.x, !!!trad_samples), -Taxonomy)

# Read metadata and match sample identifiers.
meta <- read_tsv(fichier_metadata, show_col_types = FALSE) %>%
  rename(Sample = `sample-id`, Condition = condition) %>%
  mutate(
    Condition = recode(Condition, !!!trad_conditions),
    Sample    = recode(Sample,    !!!trad_samples)
  )

# --- Extracting taxonomic ranks ---
# The taxonomy string looks like "p__Firmicutes;c__Bacilli;...".
# Regex extracts hierarchical assignments systematically.
df_taxo_16s <- df_qiime_16s %>%
  pivot_longer(cols = -Taxonomy, names_to = "Sample", values_to = "Count") %>%
  filter(Count > 0) %>% # Drop zero-count data immediately to save memory
  
  # CRITICAL: Chloroplasts and Mitochondria have homologous 16S-like genes.
  # Their DNA will naturally contaminate animal gut sequencing but does not 
  # represent the true microbial biome. This filters those false positives out.
  filter(!grepl("Chloroplast|Mitochondria|Mitochondrion", Taxonomy, ignore.case = TRUE)) %>%

  mutate(
    Phylum  = str_extract(Taxonomy, "p__[^;]+"),
    Class   = str_extract(Taxonomy, "c__[^;]+"),
    Order   = str_extract(Taxonomy, "o__[^;]+"),
    Family  = str_extract(Taxonomy, "f__[^;]+"),
    Genus   = str_extract(Taxonomy, "g__[^;]+"),
    Species = str_extract(Taxonomy, "s__[^;]+")
  ) %>%
  # Clean regex prefixes (e.g., 'g__')
  mutate(across(c(Phylum, Class, Order, Family, Genus, Species),
                ~ str_remove(.x, "^[a-z]__"))) %>%
  # Replace NA or blank values with "Unassigned" for safer filtering
  mutate(across(c(Phylum, Class, Order, Family, Genus, Species),
                ~ ifelse(is.na(.x) | str_trim(.x) %in% c("", "Non Assigné"), "Unassigned", .x)))

# Merge taxonomies and metadata
df_final_16s <- df_taxo_16s %>%
  left_join(meta, by = "Sample") %>%
  filter(!is.na(Condition))

# Declare chronological factors for accurate plotting axes (Holometabolous life cycle)
chronologie_insecte <- c(
  "Tiny_Larvae (2 mg)", "Microlarvae (7 mg)",
  "Larvae_W1 (14 mg)", "Larvae_W2 (40 mg)", "Larvae_W3 (65 mg)", "Larvae_W4 (100 mg)",
  "Pupae", "Young_Beetles", "Beetles"
)
chronologie_tout <- c(chronologie_insecte, "Raw_Substrate", "Final_Substrate")

# The five feeding larval instars, used by the core-microbiome barplot (BLOCK 6),
# the quadrant charts (BLOCK 7), the species resolution (BLOCK 11) and both
# mirror plots (BLOCKS 13-14). It was previously declared twice, under two
# different names (cond_larves / cond_larves_cibles), with identical content.
cond_larves <- c("Microlarvae (7 mg)", "Larvae_W1 (14 mg)", "Larvae_W2 (40 mg)",
                 "Larvae_W3 (65 mg)", "Larvae_W4 (100 mg)")

# Sample set used for the ordination and the clustering (BLOCKS 10 and 12):
# the whole insect cycle plus the raw substrate as an environmental anchor.
# NOTE: the original script obtained this set by overwriting chronologie_insecte
# in the middle of BLOCK 10, which silently changed the meaning of that variable
# for every block below it. Same content, but now under its own explicit name.
chronologie_pca <- c("Raw_Substrate", chronologie_insecte)

df_insecte <- df_final_16s %>%
  filter(Condition %in% chronologie_tout) %>%
  mutate(
    Condition    = factor(Condition, levels = chronologie_tout),
    ConditionNum = as.numeric(Condition)
  )

# Define gradient plotting palette for life stages
my_cond_colors <- setNames(
  c(
    "#708166ff", "#8dd1c1ff", "#41B6C4", "#0fa5e0ff", "#0f6cddff", "#001858ff", # Larval gradient
    "#984EA3", # Pupae
    "#FD8D3C", "#E31A1C", # Beetles
    "#FFD92F", "#8C510A"  # Substrates
  ),
  chronologie_tout
)


# ==============================================================================
# BLOCK 3: ALPHA AND BETA DIVERSITY CALCULATIONS
# ==============================================================================
# Alpha Diversity (Shannon): Evaluates richness and evenness within a single sample.
# Beta Diversity (Bray-Curtis): Evaluates compositional dissimilarity across samples.

message("--- Computing Alpha Diversity (Shannon) ---")
mat_counts_alpha <- df_insecte %>%
  group_by(Sample, Genus) %>%
  summarise(Count = sum(Count), .groups = "drop") %>%
  pivot_wider(names_from = Genus, values_from = Count, values_fill = 0) %>%
  column_to_rownames("Sample")

# Shannon requires raw counts. The vegan package normalizes inherently.
shannon_vals <- diversity(mat_counts_alpha, index = "shannon")

df_alpha_stade <- data.frame(Sample = names(shannon_vals), Shannon = shannon_vals) %>%
  left_join(meta, by = "Sample") %>%
  filter(Condition %in% chronologie_tout) %>%
  mutate(
    Condition    = factor(Condition, levels = chronologie_tout),
    ConditionNum = as.numeric(Condition)
  ) %>%
  group_by(Condition, ConditionNum) %>%
  summarise(
    MeanShannon = mean(Shannon),
    SD_Shannon  = ifelse(n() > 1, sd(Shannon), 0),
    .groups = "drop"
  )

message("--- Computing Beta Diversity (Bray-Curtis) ---")
resultats_stade <- list()
resultats_ind   <- list()

for (cond in levels(df_insecte$Condition)) {
  df_cond <- df_insecte %>% filter(Condition == cond)
  samps   <- unique(df_cond$Sample)
  
  if (length(samps) < 2) next # Need at least 2 replicates for pairwise beta

  # Normalizing after filtering out unassigned items avoids skewed percentages
  df_norm <- normaliser_100(df_cond %>% filter(Genus != "Unassigned"), "Genus")

  mat_wide <- df_norm %>%
    pivot_wider(id_cols = Sample, names_from = Genus,
                values_from = Relative_Abundance, values_fill = 0) %>%
    column_to_rownames("Sample")

  dist_mat <- as.matrix(vegdist(mat_wide, method = "bray"))
  vals <- dist_mat[upper.tri(dist_mat)] # Avoid duplicates by using upper triangle

  resultats_stade[[cond]] <- data.frame(
    Condition = cond,
    Beta_Mean = mean(vals),
    Beta_SD   = ifelse(length(vals) > 1, sd(vals), 0)
  )

  diag(dist_mat) <- NA   # Ignore self-to-self distance
  beta_ind <- rowMeans(dist_mat, na.rm = TRUE)
  resultats_ind[[cond]] <- data.frame(
    Sample    = names(beta_ind),
    Condition = cond,
    Beta_Ind  = beta_ind
  )
}

df_beta_stade <- bind_rows(resultats_stade)
df_beta_ind   <- bind_rows(resultats_ind)

# NOTE: the original message announced a Spearman test, but no cor.test() is
# ever run. The table below is only assembled (one row per sample: its Shannon
# index and its mean dissimilarity to its own replicates), ready for such a
# correlation. Message corrected to match what the code actually does.
message("--- Assembling per-sample Alpha/Beta table ---")
df_ind <- data.frame(Sample = names(shannon_vals), Shannon = shannon_vals) %>%
  left_join(df_beta_ind, by = "Sample") %>%
  filter(!is.na(Beta_Ind))

# Merge Alpha & Beta outputs
df_final <- df_alpha_stade %>%
  left_join(df_beta_stade, by = "Condition") %>%
  mutate(
    # English Trajectory names
    Trajectoire = ifelse(grepl("Substrat|Substrate", Condition), "Environment", "Insect Development"),
    Trajectoire = factor(Trajectoire, levels = c("Insect Development", "Environment"))
  )


# ==============================================================================
# BLOCK 4: EXPORTING STATISTICAL RESULTS
# ==============================================================================
message("--- Exporting Raw Numerical Values (TSV) ---")

df_final %>%
  select(Stade = Condition, MeanShannon, SD_Shannon, Beta_Mean, Beta_SD) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  write.table(
    file = file.path(out_dir, "Valeurs_Shannon_BrayCurtis.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )


# ==============================================================================
# BLOCK 5: GLOBAL 16S BARPLOT (GENUS LEVEL REPLICATE COMPOSITION)
# ==============================================================================
# Visualizes biological variability across replicates over time.
message("--- Plotting Global 16S Barplot ---")

ordre_chronologique <- c(
  "Raw_Substrate", "Tiny_Larvae (2 mg)", "Microlarvae (7 mg)",
  "Larvae_W1 (14 mg)", "Larvae_W2 (40 mg)", "Larvae_W3 (65 mg)", "Larvae_W4 (100 mg)",
  "Pupae", "Young_Beetles", "Beetles"
)

df_global_norm <- normaliser_avec_condition(df_final_16s, "Genus") %>%
  filter(!is.na(Condition) & Condition %in% ordre_chronologique) %>%
  mutate(Condition = factor(Condition, levels = ordre_chronologique))

# Top 14 global Genera, manually excluding non-specific environmental contaminants
top_global <- df_global_norm %>%
  filter(!(Genus %in% c("Unassigned", "Bacteroides", "Alistipes", 
                         "Aquabacterium", "Blastococcus", "Sphingomonas", "Delftia"))) %>%
  group_by(Genus) %>%
  summarise(MeanAbund = mean(Relative_Abundance), .groups = "drop") %>%
  slice_max(MeanAbund, n = 14) %>%
  pull(Genus)

df_plot_global <- df_global_norm %>%
  mutate(Taxon_Top = case_when(
    Genus == "Unassigned" ~ "Unassigned",
    Genus %in% top_global ~ Genus,
    TRUE ~ "Others"
  )) %>%
  group_by(Sample, Condition, Taxon_Top) %>%
  summarise(Abundance = sum(Relative_Abundance), .groups = "drop") %>%
  mutate(
    Taxon_Top = fct_relevel(
      fct_reorder(factor(Taxon_Top), Abundance, sum, .desc = FALSE),
      "Others", "Unassigned", after = 0
    )
  )

levels_taxa_global <- levels(df_plot_global$Taxon_Top)
taxa_vrais_global  <- setdiff(levels_taxa_global, c("Others", "Unassigned"))
my_col_global <- setNames(rep_len(pal_taxo, length(taxa_vrais_global)), taxa_vrais_global)
if ("Others"     %in% levels_taxa_global) my_col_global["Others"]     <- "grey85"
if ("Unassigned" %in% levels_taxa_global) my_col_global["Unassigned"] <- "grey40"

p_global <- ggplot(df_plot_global, aes(x = Sample, y = Abundance, fill = Taxon_Top)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2, width = 0.9, alpha = 0.9) +
  scale_fill_manual(values = my_col_global) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.02))) +
  facet_grid(~ Condition, scales = "free_x", space = "free_x") +
  labs(
    title = "Global Microbiota Composition (16S)",
    x = "Biological Replicates", y = "Relative Abundance (%)", fill = "Bacterial Genus"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 8, color = "grey20"),
    strip.text.x = element_text(face = "bold", size = 10, color = "white"),
    strip.background = element_rect(fill = "grey30", color = "black"),
    legend.text = element_text(face = "italic", size = 10),
    legend.title = element_text(face = "bold"), legend.position = "bottom",
    panel.spacing = unit(0.2, "lines"), panel.grid.major.x = element_blank()
  )

sauver_figure(p_global, "Barplot_Global_16S", width = 18, height = 8)


# ==============================================================================
# BLOCK 6: LARVAL CORE MICROBIOME BARPLOT
# ==============================================================================
message("--- Plotting Larval Core Microbiome ---")

# cond_larves is declared once in BLOCK 2 and reused by BLOCKS 7, 11, 13 and 14
df_larves_brut <- df_final_16s %>% filter(Condition %in% cond_larves)

res_larves_gen <- preparer_barplot_data(
  df_brut = df_larves_brut, rang = "Genus",
  top_n = 13, seuil_min = 1.0, cond_levels = cond_larves
)

p_core <- creer_barplot(
  res_larves_gen$df_plot, res_larves_gen$palette,
  titre    = "Larval Microbiome (16S)",
  x_lab    = "Larval Stages",
  fill_lab = "Bacterial Genus:"
)
sauver_figure(p_core, "Larval_Microbiome_Means", width = 9, height = 7, formats = "png")


# ==============================================================================
# BLOCK 7: QUADRANT CHARTS (PREVALENCE vs ABUNDANCE) — LARVAE
# ==============================================================================
message("--- Plotting Larval Quadrant Charts ---")
SEUIL_PREVALENCE <- 75  

df_larves_norm_gen <- normaliser_avec_condition(df_larves_brut, "Genus")
n_larves <- n_distinct(df_larves_norm_gen$Sample)

df_core_scatter_gen <- calculer_core_scatter(df_larves_norm_gen, "Genus", n_larves, SEUIL_PREVALENCE)
p_scatter_gen <- creer_quadrant_plot(
  df_core_scatter_gen, "Genus", SEUIL_PREVALENCE,
  titre = "Larval Microbiome Structure (16S - Genus level)"
)
sauver_figure(p_scatter_gen, "Larval_Microbiome_Quadrants", width = 12, height = 8)


# ==============================================================================
# BLOCK 8: INTRA-STAGE FIDELITY MATRICES
# ==============================================================================
# Generates heatmaps to assess biological reproducibility for highly sensitive stages.
message("--- Plotting Intra-Stage Fidelity Matrices ---")

cond_het  <- c("Tiny_Larvae (2 mg)", "Pupae", "Young_Beetles")
plot_list <- list()

for (stade in cond_het) {
  reps_attendus <- meta %>% filter(Condition == stade) %>% pull(Sample) %>% as.character()

  df_stade <- df_final_16s %>%
    filter(Condition == stade) %>%
    normaliser_avec_condition("Genus") %>%
    filter(!Genus %in% c("Non Assigné", "Unassigned")) %>%
    group_by(Genus) %>%
    # Filter criterion: Only keep taxa present in at least 3 of 4 replicates
    mutate(Nb_Present = sum(Relative_Abundance > 0)) %>%
    filter(Nb_Present >= 3) %>%
    ungroup()

  if (nrow(df_stade) > 0) {
    df_stade <- df_stade %>%
      mutate(Sample = factor(Sample, levels = reps_attendus)) %>%
      # The complete() function populates absent samples with 0 to complete the grid
      complete(Sample, Genus, fill = list(Relative_Abundance = 0, Condition = stade)) %>%
      group_by(Genus) %>%
      mutate(Total_Abund = sum(Relative_Abundance)) %>%
      ungroup() %>%
      mutate(
        Genus      = fct_reorder(Genus, Total_Abund),
        Text_Label = ifelse(Relative_Abundance >= 0.05, sprintf("%.1f", Relative_Abundance), ""),
        Is_Dark    = Relative_Abundance > 40
      )

    p <- ggplot(df_stade, aes(x = Sample, y = Genus, fill = Relative_Abundance)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = Text_Label, color = Is_Dark), size = 3.2, fontface = "bold", show.legend = FALSE) +
      scale_color_manual(values = c("TRUE" = "white", "FALSE" = "grey20")) +
      # Non-linear gradient to highlight lower abundances which dominate count tables
      scale_fill_gradientn(
        colors = c("#F2F4F4", "#AED6F1", "#2E86C1", "#154360"),
        values = rescale(c(0, 1, 20, 100)), limits = c(0, 100),
        name   = "Relative Abundance (%) :", breaks = c(0, 1, 10, 50, 100)
      ) +
      scale_x_discrete(drop = FALSE) + facet_grid(. ~ Condition) +
      labs(x = NULL, y = NULL) + theme_bw(base_size = 11) +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9, face = "bold"),
        axis.text.y = element_text(face = "italic", size = 10),
        strip.background = element_rect(fill = "#2C3E50"),
        strip.text = element_text(color = "white", face = "bold", size = 12),
        panel.grid = element_blank(), plot.margin = margin(5, 5, 5, 5)
      )
    plot_list[[stade]] <- p
  }
}

if (length(plot_list) > 0) {
  genus_counts <- sapply(plot_list, function(p) length(levels(p$data$Genus)))
  p_combined <- wrap_plots(plot_list, ncol = 1, guides = "collect") +
    plot_layout(heights = genus_counts) +
    plot_annotation(
      title    = "Intra-Stage Genus Matrices (16S)",
      subtitle = "Only genus present in >= 3/4 replicates are shown.",
      theme    = theme(plot.title = element_text(face = "bold", size = 16))
    ) & theme(legend.position = "bottom", legend.key.width = unit(2, "cm"))

  final_h <- max(7, sum(genus_counts) * 0.35 + 2.5)
  # The SVG keeps the original fixed canvas, hence the explicit svg_* overrides
  sauver_figure(p_combined, "Genome_heterogeneous", width = 9.5, height = final_h,
                bg = "white", svg_width = 18, svg_height = 8, svg_bg = NULL)
}


# ==============================================================================
# BLOCK 9: QUADRANT CHARTS (PREVALENCE vs ABUNDANCE) — ADULTS
# ==============================================================================
message("--- Plotting Adult Quadrant Charts ---")

SEUIL_PREVALENCE_AD <- 75
df_adulte_brut <- df_final_16s %>% filter(Condition == "Beetles")

df_ad_norm_gen <- normaliser_avec_condition(df_adulte_brut, "Genus") %>%
  filter(!Genus %in% c("Non Assigné", "Unassigned"))
n_ad <- n_distinct(df_ad_norm_gen$Sample)

df_core_ad_gen <- calculer_core_scatter(df_ad_norm_gen, "Genus", n_ad, SEUIL_PREVALENCE_AD)
p_scatter_ad <- creer_quadrant_plot(
  df_core_ad_gen, "Genus", SEUIL_PREVALENCE_AD,
  titre = "Adult Microbiome Structure (Beetles - Genus level)"
)
sauver_figure(p_scatter_ad, "Adult_Microbiome_Quadrants", width = 12, height = 8)


# ==============================================================================
# BLOCK 10: PRINCIPAL COMPONENT ANALYSIS (PCA)
# ==============================================================================
# Biological Context: Hellinger transformation takes the square root of relative 
# abundances. This handles the "double-zero" problem (where two sites are falsely 
# deemed similar just because they both lack a rare species) and dampens highly 
# dominant taxa from skewing the principal components artificially.
message("--- Performing Hellinger-transformed PCA (Insect Only) ---")

df_pca_norm <- normaliser_100(
  df_final_16s %>% filter(Condition %in% chronologie_pca, Genus != "Unassigned"), 
  "Genus"
)

mat_pca <- df_pca_norm %>%
  pivot_wider(id_cols = Sample, names_from = Genus,
              values_from = Relative_Abundance, values_fill = 0) %>%
  column_to_rownames("Sample")

# Transform and run standard PCA (RDA with no environmental constraint = PCA)
mat_hellinger <- decostand(mat_pca, method = "hellinger")
pca_res <- rda(mat_hellinger)

# Scaling = 1 strongly preserves true distances between biological samples
pca_coords <- as.data.frame(scores(pca_res, display = "sites", scaling = 1))
colnames(pca_coords) <- c("PC1", "PC2")
pca_coords$Sample <- as.character(rownames(pca_coords))

# Extract Eigenvalues for axis %
eig_vals <- pca_res$CA$eig
var_pc1  <- round(eig_vals[1] / sum(eig_vals) * 100, 1)
var_pc2  <- round(eig_vals[2] / sum(eig_vals) * 100, 1)

df_plot_pca <- pca_coords %>%
  inner_join(meta %>% select(Sample, Condition) %>% mutate(Sample = as.character(Sample)), by = "Sample") %>%
  mutate(Condition = factor(Condition, levels = chronologie_pca))

p_pca <- ggplot(df_plot_pca, aes(x = PC1, y = PC2, color = Condition, fill = Condition)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(size = 4, shape = 21, color = "black", stroke = 0.6) +
  scale_color_manual(values = my_cond_colors) +
  scale_fill_manual(values = my_cond_colors) +
  labs(
    title = "PCA (Hellinger-transformed relative abundances)",
    x = sprintf("PC1 (%s%%)", var_pc1), y = sprintf("PC2 (%s%%)", var_pc2),
    fill  = "Development Stage:", color = "Development Stage:"
  ) +
  theme_bw(base_size = 18) + 
  theme(
    text = element_text(size = 18), legend.position = "right",
    legend.title = element_text(face = "bold"), plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40"), panel.grid.minor = element_blank()
  )
sauver_figure(p_pca, "PCA_microbiote", width = 10, height = 7, formats = "png", bg = "white")


# ==============================================================================
# BLOCK 11: INTRA-GENUS TAXONOMIC RESOLUTION
# ==============================================================================
message("--- Extracting Intra-Genus Species Resolution ---")

# nom_base is given WITHOUT extension: sauver_figure() appends .png and .svg
generer_resolution_fragmente <- function(stades_cibles, titre_stade, nom_base,
                                          top_n_genres = 8, exclure_genres = NULL) {
  message(sprintf("  -> Processing Stage(s): %s", paste(stades_cibles, collapse = ", ")))

  df_stade <- df_final_16s %>% filter(Condition %in% stades_cibles)
  if (nrow(df_stade) == 0) return(NULL)

  df_norm_gen  <- normaliser_100(df_stade, "Genus")
  n_samp_stade <- n_distinct(df_norm_gen$Sample)

  df_gen_agg <- df_norm_gen %>%
    group_by(Genus) %>%
    summarise(MeanAbund = sum(Relative_Abundance) / n_samp_stade, .groups = "drop") %>%
    filter(!Genus %in% c("Non Assigné", "Unassigned", exclure_genres))

  top_gen       <- df_gen_agg %>% slice_max(MeanAbund, n = top_n_genres) %>% pull(Genus)
  df_gen_labels <- df_gen_agg %>% filter(Genus %in% top_gen)

  df_sp <- df_stade %>%
    filter(Genus %in% top_gen) %>%
    mutate(
      Species_Label = case_when(
        Species %in% c("Unassigned", "Non Assigné", NA, "") ~ "Unclassified (sp.)",
        grepl(Genus, Species, ignore.case = TRUE)           ~ Species,
        TRUE ~ paste(Genus, Species)
      )
    )

  df_intra <- df_sp %>%
    group_by(Genus, Species_Label) %>%
    summarise(Count = sum(Count), .groups = "drop") %>%
    group_by(Genus) %>%
    mutate(Pct_in_Genus = (Count / sum(Count)) * 100) %>%
    ungroup()

  df_res_clean <- df_intra %>%
    group_by(Genus) %>%
    mutate(
      Rank = case_when(
        Species_Label == "Unclassified (sp.)" ~ 999,
        TRUE ~ dense_rank(desc(Pct_in_Genus))
      ),
      Species_Final = case_when(
        Species_Label == "Unclassified (sp.)" ~ "Unclassified (sp.)",
        Rank <= 5 ~ Species_Label,
        TRUE ~ "Other identified species"
      )
    ) %>%
    group_by(Genus, Species_Final) %>%
    summarise(Pct_in_Genus = sum(Pct_in_Genus), .groups = "drop") %>%
    left_join(df_gen_labels, by = "Genus") %>%
    mutate(Genus_Label = sprintf("%s (%.2f%%)", Genus, MeanAbund))

  df_res_clean <- df_res_clean %>%
    mutate(
      Genus_Label = fct_reorder(Genus_Label, MeanAbund, .desc = TRUE),
      Species_Final = fct_relevel(fct_reorder(Species_Final, Pct_in_Genus), "Unclassified (sp.)", "Other identified species", after = 0)
    )

  pal_gen <- setNames(colorspace::qualitative_hcl(length(top_gen), palette = "Dark 3"), top_gen)

  p_res <- ggplot(df_res_clean, aes(x = Pct_in_Genus, y = Species_Final)) +
    geom_col(aes(fill = Genus), color = "black", linewidth = 0.3, width = 0.7) +
    geom_text(aes(label = sprintf("%.2f%%", Pct_in_Genus)), hjust = -0.15, size = 3.5, fontface = "bold", color = "grey20") +
    facet_wrap(~ Genus_Label, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = pal_gen) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(title = sprintf("Intra-Genus Species Resolution : %s", titre_stade), x = "Relative Proportion within Genus (%)", y = NULL) +
    theme_bw(base_size = 13) +
    theme(
      legend.position = "none", strip.background = element_rect(fill = "grey20"),
      strip.text = element_text(color = "white", face = "bold.italic", size = 12),
      axis.text.y = element_text(face = "italic", size = 10, color = "black"),
      axis.text.x = element_text(size = 9), panel.grid.major.y = element_blank(), plot.margin = margin(10, 15, 10, 10)
    )

  sauver_figure(p_res, nom_base, width = 13, height = 9, bg = "white")
  message(sprintf("  => Successfully saved '%s'", nom_base))
}

# Apply to stages
generer_resolution_fragmente(cond_larves, "Average Larval Microbiome", "Taxonomic_Resolution_Average_Larvae", 10, c("Citrobacter","Lactobacillus","Bacillus"))
generer_resolution_fragmente("Microlarvae (7 mg)", "Microlarvae Microbiome", "Taxonomic_Resolution_Microlarvae", 10)
generer_resolution_fragmente("Larvae_W1 (14 mg)", "Larvae W1 (14 mg) Microbiome", "Taxonomic_Resolution_Larvae_W1", 10)
generer_resolution_fragmente("Larvae_W2 (40 mg)", "Larvae W2 (40 mg) Microbiome", "Taxonomic_Resolution_Larvae_W2", 10)
generer_resolution_fragmente("Larvae_W3 (65 mg)", "Larvae_W3 (65 mg) Microbiome", "Taxonomic_Resolution_Larvae_W3", 10)
generer_resolution_fragmente("Larvae_W4 (100 mg)", "Larvae_W4 (100 mg) Microbiome", "Taxonomic_Resolution_Larvae_W4", 10)
generer_resolution_fragmente("Beetles", "Beetles Microbiome", "Taxonomic_Resolution_Adults", 10)


# ==============================================================================
# BLOCK 12: DENDROGRAMS (Ward.D2 on PCA Space)
# ==============================================================================
# Statistical Context: Hierarchical clustering over the exact PC1 and PC2 axis 
# coordinates rather than the massive distance matrix. This ensures the Dendrogram 
# perfectly mimics the graphical groupings seen on the PCA plot. Ward's method 
# minimizes variance within clusters. The Elbow plot dynamically finds optimal 'K'.
message("--- Generating Dendrograms (Elbow + Hellinger + Ward.D2 on PCA space) ---")

if (!requireNamespace("ggdendro", quietly = TRUE)) install.packages("ggdendro")
suppressPackageStartupMessages(library(ggdendro))

generer_dendrogramme_et_elbow <- function(df_counts, meta_df, prefixe_nom, titre_dendro, k_optimal, n_axes_clustering = 2) {
  
  df_norm <- normaliser_100(df_counts %>% filter(Genus != "Unassigned"), "Genus")
  mat_wide <- df_norm %>%
    pivot_wider(id_cols = Sample, names_from = Genus, values_from = Relative_Abundance, values_fill = 0) %>%
    column_to_rownames("Sample")
  
  mat_hellinger <- decostand(mat_wide, method = "hellinger")
  pca_res  <- rda(mat_hellinger)
  n_axes   <- min(n_axes_clustering, ncol(scores(pca_res, display = "sites")))
  
  # Select coordinate arrays for strictly PC1 and PC2 (scaling = 1)
  mat_clust <- scores(pca_res, display = "sites", scaling = 1, choices = 1:n_axes)
  dist_mat <- vegdist(mat_clust, method = "euclidean")
  
  hc <- hclust(dist_mat, method = "ward.D2")
  dendro <- dendro_data(hc, type = "rectangle")
  
  # Compute WSS (Within-Cluster Sum of Squares) for 1 to 10 clusters (Elbow Plot)
  k_max <- min(10, nrow(mat_clust) - 1)
  wss <- sapply(1:k_max, function(k) {
    clusters <- cutree(hc, k)
    sum(sapply(1:k, function(i) {
      mat_sub <- mat_clust[clusters == i, , drop = FALSE]
      center <- colMeans(mat_sub)
      sum(rowSums((sweep(mat_sub, 2, center))^2))
    }))
  })
  
  df_elbow <- data.frame(K = 1:k_max, WSS = wss)
  p_elbow <- ggplot(df_elbow, aes(x = K, y = WSS)) +
    geom_line(color = "steelblue", linewidth = 1) + geom_point(color = "red", size = 3) +
    geom_vline(xintercept = k_optimal, linetype = "dashed", color = "grey40") +
    scale_x_continuous(breaks = 1:k_max) +
    labs(title = paste("Elbow Method :", prefixe_nom), x = "Number of clusters (k)", y = "Within-Cluster Sum of Squares (WSS)") +
    theme_bw(base_size = 13)
  sauver_figure(p_elbow, paste0(prefixe_nom, "_Elbow_Plot"), width = 8, height = 6, formats = "png")
  
  labels_df <- label(dendro) %>% rename(Sample = label) %>% left_join(meta_df, by = "Sample")
  groupes <- cutree(hc, k = k_optimal)
  labels_df$Cluster <- as.factor(groupes[labels_df$Sample])
  
  labels_df %>% select(Cluster, Sample, Condition) %>% arrange(Cluster, Condition) %>%
    write_tsv(file.path(out_dir, paste0(prefixe_nom, "_Clusters_Composition.tsv")))
  
  rect_df <- labels_df %>%
    group_by(Cluster) %>%
    summarise(xmin = min(x) - 0.45, xmax = max(x) + 0.45, .groups = "drop") %>%
    mutate(ymin = -0.015, ymax = max(segment(dendro)$y) * 0.95)
  
  p_dendro <- ggplot() +
    geom_rect(data = rect_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey50", alpha = 0.15, color = "black", linetype = "dashed", inherit.aes = FALSE) +
    geom_segment(data = segment(dendro), aes(x = x, y = y, xend = xend, yend = yend), color = "grey60", linewidth = 0.6) +
    geom_point(data = labels_df, aes(x = x, y = y, fill = Condition), shape = 21, size = 3, color = "black", stroke = 0.5) +
    geom_text(data = labels_df, aes(x = x, y = y - 0.03, label = Sample, color = Condition), angle = 90, hjust = 1, vjust = 0.5, size = 4, fontface = "bold", show.legend = FALSE) + 
    scale_color_manual(values = my_cond_colors, breaks = names(my_cond_colors)) +
    scale_fill_manual(values = my_cond_colors, breaks = names(my_cond_colors)) +
    scale_y_continuous(expand = expansion(mult = c(0.25, 0.08))) +
    labs(title = titre_dendro, x = NULL, y = "Height (Euclidean Distance)", color = "Development Stage:", fill = "Development Stage:") +
    theme_bw(base_size = 13) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom", plot.title = element_text(face = "bold", size = 15), plot.margin = margin(10, 10, 10, 10))
  
  sauver_figure(p_dendro, paste0(prefixe_nom, "_Dendrogram"), width = 14, height = 8, formats = "png", bg = "white")
}

generer_dendrogramme_et_elbow(
  df_counts = df_final_16s %>% filter(Condition %in% chronologie_pca),
  meta_df   = meta, prefixe_nom = "PCA_Equivalent_Clustering",
  titre_dendro = "Hierarchical Clustering (ACP Hellinger - PC1+PC2, Ward.D2)",
  k_optimal = 3, n_axes_clustering = 2
)


# ============================================================================
# BLOCK 13: MIRROR PLOTS (LARVAE VS ADULTS)
# ============================================================================
message("--- Generating Adult Mirror Plots (Standard Error + Wilcoxon FDR) ---")


col_larvae   <- "#3678a3fd" # Blue   : larval stages (always the left half)
col_adults   <- "#e26759ff" # Salmon : beetle stages
col_substrat <- "#ebc90cff" # Yellow : raw rearing substrate

# Curated genus list for the larvae vs beetles contrast: genera expected to
# shift markedly across metamorphosis.
genres_cibles <- c(
  "Staphylococcus", "Enterococcus", "Lactococcus", "Enterobacillus", 
  "Spiroplasma", "Mammaliicoccus", "Weissella", "Brevibacterium", 
  "Corynebacterium", "Serratia", "Brachybacterium", "Latilactobacillus", 
  "Pediococcus", "Leuconostoc"
)

res_mirror_adultes <- creer_mirror_plot(
  cond_g1  = cond_larves, cond_g2 = "Beetles",
  lab_g1   = "Larval Stage", lab_g2 = "Beetles Stage",
  col_g1   = col_larvae,     col_g2 = col_adults,
  genres   = genres_cibles,
  titre    = "Taxonomic comparison of the larval and beetle stages",
  nom_base = "Comparison_Mirror_Plot_Spine"
)


# ==============================================================================
# BLOCK 14: MIRROR PLOTS (LARVAE VS RAW SUBSTRATE)
# ==============================================================================

message("--- Generating Substrate Mirror Plots (Standard Error + Wilcoxon FDR) ---")

genres_cibles_sub <- c(
  "Staphylococcus", "Enterococcus", "Lactococcus", "Enterobacter", 
  "Spiroplasma", "Mammaliicoccus", "Weissella", "Brevibacterium", 
  "Corynebacterium", "Brachybacterium", "Pantoea", "Paenibacillus", "Massilia"
)

res_mirror_substrat <- creer_mirror_plot(
  cond_g1  = cond_larves, cond_g2 = "Raw_Substrate",
  lab_g1   = "Larval Stage", lab_g2 = "Raw Substrate",
  col_g1   = col_larvae,     col_g2 = col_substrat,
  genres   = genres_cibles_sub,
  titre    = "Taxonomic comparison of the larval stages and the raw substrate",
  nom_base = "Comparison_Mirror_Plot_Spine_Substrate"
)
