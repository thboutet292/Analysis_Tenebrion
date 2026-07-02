################################################################################
# SCRIPT: 16S MICROBIAL DIVERSITY ANALYSIS — Tenebrio molitor
################################################################################
#
# This script analyzes the gut bacterial communities of T. molitor across
# the different developmental stages (larvae -> adults + substrate)
#
# OUTPUT FILES (in results/16S/alpha_beta/):
#
# --- Alpha / Beta Diversity ---
#  1. Valeurs_Shannon_BrayCurtis.tsv                  -> Raw numeric table (appendix)
#  2. shannon_vs_braycurtis.png / .svg                -> Alpha/beta curves (insect only)
#
# --- Overall composition ---
#  3. Barplot_Global_16S.png / .svg                   -> Composition by replicate (insect only)
#  4. Matrice_BrayCurtis_Tous_Replicats.png           -> Dissimilarity heatmap (all replicates)
#  5. Matrice_BrayCurtis_Insecte_Seul.png             -> Dissimilarity heatmap (insect only)
#
# --- Master figure ---
#  6. Master_Figure_substrats.png / .svg              -> Combined figure (insect + environment)
#  7. Master_Figure.png / .svg                        -> Combined figure (insect only)
#
# --- Larval core microbiome ---
#  8. Larval_Microbiome_Means.png                     -> Larval core barplot (Genus)
#  9. Larval_Microbiome_Means_Species.png             -> Larval core barplot (Species)
# 10. Larval_Microbiome_Quadrants.png                 -> Prevalence/Abundance quadrants (Genus)
# 11. Larval_Microbiome_Quadrants_Species.png         -> Prevalence/Abundance quadrants (Species)
# 12. Genome_heterogeneous.png                        -> Intra-stage fidelity matrices
#
# --- Adult core microbiome ---
# 13. Adult_Microbiome_Means.png                      -> Adult barplot (Genus)
# 14. Adult_Microbiome_Donut.png                      -> Adult donut chart (Genus)
# 15. Adult_Microbiome_Means_Species.png              -> Adult barplot (Species)
# 16. Adult_Microbiome_Donut_Species.png              -> Adult donut chart (Species)
# 17. Adult_Microbiome_Quadrants.png                  -> Adult quadrants (Genus)
# 18. Adult_Microbiome_Quadrants_Species.png          -> Adult quadrants (Species)
#
# --- Community structure ---
# 19. PCA_microbiote_insecte_seul.png                 -> PCA without ellipses (insect only)
# 20. PCA_microbiote_insecte_seul_ellipses.png        -> PCA with 95% ellipses (insect only)
# 21. Correlation_Bray_vs_Shannon.png / .svg          -> Alpha vs beta correlation
#
# --- Intra-genus taxonomic resolution ---
# 22. Taxonomic_Resolution_Average_Larvae.png         -> Species-level resolution (pooled larvae)
# 23. Taxonomic_Resolution_Microlarvae.png            -> Species-level resolution (Microlarvae 7 mg)
# 24. Taxonomic_Resolution_Larvae_S1.png              -> Species-level resolution (Larvae S1 14 mg)
# 25. Taxonomic_Resolution_Larvae_S2.png              -> Species-level resolution (Larvae S2 40 mg)
# 26. Taxonomic_Resolution_Larvae_S3.png              -> Species-level resolution (Larvae S3 65 mg)
# 27. Taxonomic_Resolution_Larvae_S4.png              -> Species-level resolution (Larvae S4 100 mg)
# 28. Taxonomic_Resolution_Adults.png                 -> Species-level resolution (adults)
#
# --- Dendrograms ---
# 29. Dendrogramme_Global_Epure.png                   -> Hierarchical clustering (all samples)
# 30. Dendrogramme_Insecte_Epure.png                  -> Hierarchical clustering (insect only)
#
# REQUIRED DATA:
#   - QIIME2-like table (sequences x samples, TSV format)
#   - Metadata file (sample -> developmental condition)
#   - SILVA v138.2 databases for prior taxonomic assignment
################################################################################


# ==============================================================================
# BLOCK 0: CONFIGURATION AND PACKAGES
# ==============================================================================

setwd("/home/thomas/Tenebrion/")

suppressPackageStartupMessages({
  library(dplyr)       
  library(ggplot2)    
  library(tidyr)       
  library(vegan)       
  library(tibble)      
  library(stringr)     
  library(readr)       
  library(forcats)     
  library(colorspace)  
  library(patchwork)   
  library(ggrepel)     
  library(scales)
  library(readxl)
})

# File paths
fichier_qiime    <- "results/all_matam_salmon_qiime_like_table_HYBRIDE.tsv"
fichier_metadata <- "data/16S/metadata.tsv"
out_dir          <- "results/16S/alpha_beta_test"

# Creates the output folder if it does not already exist.
# recursive = TRUE: also creates missing parent folders if needed
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Color palette for taxa 
pal_taxo <- colorspace::qualitative_hcl(14, palette = "Dark 3")

# Diversity curve colors 
COL_SHANNON <- "#1F78B4"          # alpha diversity
COL_BRAY    <- "#E31A1C"          # beta diversity


# ==============================================================================
# BLOCK 1: UTILITY FUNCTIONS
# ==============================================================================
# ------------------------------------------------------------------------------
# lire_qiime_tsv_robuste(): Robust reading of a QIIME2 count file
# ------------------------------------------------------------------------------
# The QIIME2 format has a quirk: the header line starts with "#"
# (e.g. "#OTU ID\tSample1\tSample2..."). This function automatically detects
# and cleans up that header so R can read the file correctly,
# regardless of the exact format version.
#
lire_qiime_tsv_robuste <- function(chemin) {
  lignes <- readLines(chemin, warn = FALSE)

  # Header line detection: we look for "OTU ID" or "taxonomy"
  idx_header <- which(grepl("(?i)OTU.?ID", lignes, perl = TRUE))[1]
  if (is.na(idx_header)) {
    idx_header <- which(grepl("(?i)taxonomy", lignes, perl = TRUE))[1]
  }

  lignes_propres <- lignes[idx_header:length(lignes)]

  # Removal of the leading "#": read_tsv() would otherwise interpret it as a comment
  lignes_propres[1] <- sub("^#", "", lignes_propres[1])

  # I() wraps the text so read_tsv() reads it as if it were
  # an actual file on disk (text connection)
  df <- read_tsv(I(paste(lignes_propres, collapse = "\n")), show_col_types = FALSE)
  names(df)[1] <- "Taxonomy"
  return(df)
}


# ------------------------------------------------------------------------------
# harmoniser_taxonomie(): Harmonization of bacterial genus names
# ------------------------------------------------------------------------------
# Biological problem: 16S databases (SILVA, Greengenes, NCBI) have
# undergone several major taxonomic revisions in recent years.
# Older genera have been split, renamed, or merged.
#
# Examples:
#   - Escherichia and Shigella are phylogenetically indistinguishable by 16S
#   - The genus Lactobacillus was split into ~25 distinct genera in SILVA 138
#   - Bacillus was split into ~10 genera (Cytobacillus, Priestia, etc.)
#
# This function applies a correction dictionary to ensure consistent
# names across the different SILVA versions and publications.
#
# Arguments:
#   tab: data.frame containing at least the Genus and Family columns
# Returns:
#   The same data.frame with harmonized genus names
harmoniser_taxonomie <- function(tab) {
  tab %>%
    mutate(
      Genus = case_when(
        # --- Genera indistinguishable by 16S ---
        # Escherichia and Shigella share >97% identity in 16S
        Genus %in% c("Escherichia", "Shigella",
                     "Escherichia-Shigella", "Escherichia/Shigella") ~ "Escherichia-Shigella",

        # --- Genus complexes revised in SILVA 138 ---
        # Lactobacillus: split into ~25 genera (Ligilactobacillus, Lacticaseibacillus...)
        grepl("lactobacillus", Genus, ignore.case = TRUE) ~ "Lactobacillus",

        # Mycobacterium: major 2018 revision (Gupta et al.)
        Genus %in% c("Mycobacterium", "Mycolicibacterium", "Mycolicibacter",
                     "Mycobacteroides", "Mycolicibacillus") ~ "Mycobacterium",

        # Bacillus: split into ~10 genera in SILVA 138
        Genus %in% c("Bacillus", "Cytobacillus", "Mesobacillus", "Neobacillus",
                     "Peribacillus", "Alkalihalobacillus", "Litchfieldia",
                     "Metabacillus", "Priestia") ~ "Bacillus",

        # Pseudomonas: a few satellite genera reintegrated
        Genus %in% c("Pseudomonas", "Halopseudomonas", "Neopseudomonas") ~ "Pseudomonas",

        # Burkholderia: heavily reshuffled genus in the 2010s
        Genus %in% c("Burkholderia", "Cupriavidus", "Paraburkholderia",
                     "Caballeronia", "Trinickia") ~ "Burkholderia",

        # Mycoplasma: major 2019 taxonomic revision (Gupta et al.)
        Genus %in% c("Mycoplasma", "Mesomycoplasma", "Metamycoplasma",
                     "Mycoplasmopsis", "Mycoplasmoides", "Malacoplasma",
                     "Ureaplasma", "Williamsoniiplasma") ~ "Mycoplasma",

        # Pantoea and related genera (Enterobacteriaceae)
        Genus %in% c("Pantoea", "Erwinia", "Tatumella") ~ "Pantoea",

        # Enterobacter complex: cryptic genera from recent revisions
        Genus %in% c("Mixta", "Dickeya", "Leclercia", "Mangrovibacter",
                     "Enterobacillus", "Jejubacter", "Pluralibacter",
                     "Klebsiella", "Raoultella") ~ "Enterobacter",

        # Individual reclassifications
        Genus == "Vagococcus"  ~ "Enterococcus",  # Phylogenetically nested
        Genus == "Enhydrobacter" ~ "Moraxella",   # 2012 reclassification

        # Everything else: keep the original name
        TRUE ~ Genus
      )
    ) %>%
    mutate(
      Family = case_when(
        # If Family is missing but Genus is known,
        # create an informative name so the taxonomic information is not lost.
        # This keeps a trace of the position in the tree of life.
        Family == "Unassigned" & Genus != "Unassigned" ~ paste0("Famille_de_", Genus),
        TRUE ~ Family
      )
    )
}


# ------------------------------------------------------------------------------
# normaliser_100(): Conversion of raw counts to relative abundances (%)
# ------------------------------------------------------------------------------
# Technical problem: Sequencing depth varies between samples.
# A sample with 10,000 reads and one with 50,000 reads are not
# directly comparable using raw values.
# Solution: Each sample is converted to a percentage (sum = 100%),
# which allows comparing compositions independently of the total volume
#
normaliser_100 <- function(df, rang) {
  df %>%
    # Step 1: Sum reads for each taxon per sample.
    # Several rows can point to the same genus (different sequences
    # but the same taxonomic assignment) here they are merged
    group_by(Sample, !!sym(rang)) %>%   # sym() converts the string into an R symbol
    summarise(Count = sum(Count), .groups = "drop") %>%

    # Step 2: Compute the percentage of each taxon within its sample
    group_by(Sample) %>%
    mutate(Relative_Abundance = (Count / sum(Count)) * 100) %>%
    ungroup()
}


# ------------------------------------------------------------------------------
# forcer_nomenclature_binomiale(): Building the full Genus + Species name
# ------------------------------------------------------------------------------
# In microbiology, binomial nomenclature convention requires that a
# species name always include the genus (e.g. "Escherichia coli" and not "coli").
# SILVA sometimes stores the species without the genus -> this function ensures
# species-name consistency across the whole script.
#
# Arguments:
#   df: data.frame containing the Genus and Species columns
# Returns:
#   The same data.frame with the corrected Species column
forcer_nomenclature_binomiale <- function(df) {
  df %>%
    mutate(Species = case_when(
      # Case 1: Unassigned species -> state it explicitly
      Species %in% c("Non Assigné", "Unassigned") ~ "Unassigned",
      # Case 2: Unassigned genus -> a binomial name cannot be built
      Genus   %in% c("Non Assigné", "Unassigned") ~ Species,
      # Case 3: The species already contains the genus (e.g. "Lactobacillus acidophilus") -> OK
      grepl(Genus, Species, ignore.case = TRUE)   ~ Species,
      # Case 4: The species stands alone (e.g. "acidophilus") -> prefix it with the genus
      TRUE ~ paste(Genus, Species)
    ))
}


# ------------------------------------------------------------------------------
# preparer_barplot_data(): Full preparation pipeline for a barplot
# ------------------------------------------------------------------------------
# This pattern is used identically for larvae, adults, and the global dataset.
# Factoring it into one function avoids 6 nearly identical code blocks.
#
# Arguments:
#   df_brut   : raw data filtered for the desired stage
#   rang      : "Genus" or "Species"
#   top_n     : number of taxa to display individually (e.g. 13)
#   seuil_min : minimum abundance threshold to enter the Top (e.g. 1.0%)
#   cond_levels : vector of conditions in the desired chronological order
#   exclus    : vector of genera/species to exclude from the Top (e.g. "Unassigned")
# Returns:
#   List with $df_plot (data for ggplot) and $palette (colors)
preparer_barplot_data <- function(df_brut, rang, top_n, seuil_min = 0,
                                   cond_levels = NULL, exclus = c("Non Assigné", "Unassigned")) {
  # --- 1. Normalization to 100% per sample ---
  df_norm <- normaliser_100(df_brut, rang) %>%
    left_join(meta %>% select(Sample, Condition), by = "Sample")

  # --- 2. Computing the true mean per stage ---
  # The "true mean" = sum of all relative abundances / number of replicates
  # This guarantees that the sum of the bars = exactly 100% (no averaging artifact)
  df_agg <- df_norm %>%
    group_by(Condition) %>%
    mutate(Nb_Samples = n_distinct(Sample)) %>%
    group_by(Condition, !!sym(rang), Nb_Samples) %>%
    summarise(Relative_Abundance = sum(Relative_Abundance) / first(Nb_Samples),
              .groups = "drop")

  if (!is.null(cond_levels)) {
    df_agg <- df_agg %>% mutate(Condition = factor(Condition, levels = cond_levels))
  }

  # --- 3. Identifying the Top N most abundant taxa ---
  top_taxa <- df_agg %>%
    filter(!.data[[rang]] %in% exclus) %>%
    group_by(.data[[rang]]) %>%
    summarise(MeanAbund = mean(Relative_Abundance), .groups = "drop") %>%
    filter(MeanAbund >= seuil_min) %>%
    slice_max(MeanAbund, n = top_n, with_ties = FALSE) %>%
    pull(.data[[rang]])

  # --- 4. Grouping rare taxa into "Others" ---
  # Taxa outside the Top N are grouped to avoid an unreadable legend.
  # "Others" and "Unassigned" are forced to the very bottom of the bars.
  df_plot <- df_agg %>%
    mutate(Taxon_Top = case_when(
      .data[[rang]] %in% exclus ~ "Unassigned",
      .data[[rang]] %in% top_taxa ~ as.character(.data[[rang]]),
      TRUE ~ "Others"
    )) %>%
    group_by(Condition, Taxon_Top) %>%
    summarise(Abundance = sum(Relative_Abundance), .groups = "drop") %>%
    mutate(
      # fct_reorder() sorts by total abundance (most abundant on top)
      # fct_relevel() then forces "Others" and "Unassigned" to the very bottom
      Taxon_Top = fct_relevel(
        fct_reorder(factor(Taxon_Top), Abundance, sum, .desc = FALSE),
        "Others", "Unassigned",
        after = 0
      )
    )

  # --- 5. Building the color palette ---
  # Real taxa receive qualitative colors (pal_taxo)
  # Generic categories receive neutral greys so they are not confused
  # with real taxa in the legend
  levels_taxa <- levels(df_plot$Taxon_Top)
  taxa_vrais  <- setdiff(levels_taxa, c("Others", "Unassigned"))

  palette <- setNames(rep_len(pal_taxo, length(taxa_vrais)), taxa_vrais)
  if ("Others"     %in% levels_taxa) palette["Others"]     <- "grey85"
  if ("Unassigned" %in% levels_taxa) palette["Unassigned"] <- "grey40"

  return(list(df_plot = df_plot, palette = palette))
}


# ------------------------------------------------------------------------------
# creer_barplot(): Creating a standardized stacked barplot
# ------------------------------------------------------------------------------
# Arguments:
#   df_plot     : output of preparer_barplot_data()$df_plot
#   palette     : output of preparer_barplot_data()$palette
#   titre, x_lab, y_lab, fill_lab : chart texts
#   angle_x     : angle of the X-axis labels (default = 45°)
#   largeur_bar : bar width (default = 0.65)
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
# creer_quadrant_plot(): Prevalence vs Abundance chart (Core Microbiome)
# ------------------------------------------------------------------------------
# This type of chart classifies each bacterium into 4 quadrants based on:
#   - Its PREVALENCE: in what % of samples is it present?
#   - Its mean ABUNDANCE: what % of sequences does it represent?
#
# Ecological classification:
#   Q1 (red)    : High abundance + High prevalence -> Strict core microbiome
#   Q2 (blue)   : Low abundance + High prevalence -> Satellite core microbiome
#   Q3 (orange) : High abundance + Low prevalence -> Transient bacterium
#   Q4 (grey)   : Low abundance + Low prevalence -> Background noise
#
# Arguments:
#   df_scatter : data.frame with columns Prevalence, MeanAbund, Quadrant, [rang]
#   rang       : column used for labels ("Genus" or "Species")
#   seuil_prev, seuil_abond : threshold lines to draw
#   titre, sous_titre : chart texts
#   style_label : "bold" for genera, "italic" for species
creer_quadrant_plot <- function(df_scatter, rang, seuil_prev, seuil_abond,
                                 titre, sous_titre = NULL, style_label = "bold") {

  # Semantic colors: red = important, blue = common but rare,
  # orange = occasional but abundant, grey = negligible
  palette_quadrant <- c(
    "1. High Abundance Shared Microbiota" = "#E31A1C",
    "2. Low Abundance Shared Microbiota"  = "#1F78B4",
    "3. Transient"                        = "#FF7F00",
    "4. Background Noise"                 = "grey75"
  )

  ggplot(df_scatter, aes(x = Prevalence, y = MeanAbund)) +
    # Threshold lines (the "boundaries" between quadrants)
    geom_vline(xintercept = seuil_prev,  linetype = "dashed", color = "grey40", linewidth = 0.8) +
    geom_hline(yintercept = seuil_abond, linetype = "dashed", color = "grey40", linewidth = 0.8) +

    # Points whose size is proportional to abundance
    geom_point(aes(color = Quadrant, size = MeanAbund), alpha = 0.8) +

    # Smart labels (ggrepel automatically avoids overlaps)
    # Background noise is not labeled to keep the chart readable
    geom_text_repel(
      data = filter(df_scatter, Quadrant != "4. Background Noise"),
      aes(label = .data[[rang]], color = Quadrant),
      size = 3.5, fontface = style_label,
      box.padding = 0.62, max.overlaps = 30, show.legend = FALSE
    ) +

    # Log-scale Y axis: allows viewing both rare taxa
    # (0.01%) and abundant ones (>10%) on the same chart
    scale_y_log10(
      labels = scales::comma_format(accuracy = 0.01),
      breaks = c(0.01, 0.1, 1, 10, 100)
    ) +
    scale_x_continuous(breaks = seq(0, 100, 25), limits = c(0, 105)) +
    scale_size_continuous(range = c(2, 10), guide = "none") +
    scale_color_manual(values = palette_quadrant) +

    labs(
      title    = titre,
      subtitle = sous_titre %||% sprintf(
        "Thresholds: Prevalence >= %d%% | Abundance >= %.1f%%",
        seuil_prev, seuil_abond
      ),
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

# "Null-coalescing" operator: returns the 2nd argument if the 1st is NULL
# (used in creer_quadrant_plot for the default subtitle)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ------------------------------------------------------------------------------
# calculer_core_scatter(): Computing Prevalence + Mean Abundance metrics
# ------------------------------------------------------------------------------
# For each taxon, we compute:
#   - MeanAbund   : mean relative abundance across ALL samples of the group
#   - Prevalence  : % of samples in which this taxon is detected (abundance > 0)
#
# Arguments:
#   df_norm   : data normalized to 100% (output of normaliser_100)
#   rang      : "Genus" or "Species"
#   n_samples : total number of samples in the group
#   seuil_prev, seuil_abond : thresholds for quadrant classification
# Returns:
#   data.frame with columns: [rang], MeanAbund, Prevalence, Quadrant
calculer_core_scatter <- function(df_norm, rang, n_samples, seuil_prev, seuil_abond) {
  df_norm %>%
    filter(!.data[[rang]] %in% c("Non Assigné", "Unassigned")) %>%
    group_by(.data[[rang]]) %>%
    summarise(
      # Mean abundance: sum across all samples / total number of samples
      # (samples where the bacterium is absent contribute 0 to the numerator)
      MeanAbund  = sum(Relative_Abundance) / n_samples,
      # Prevalence: proportion of samples where abundance is strictly > 0
      Prevalence = sum(Relative_Abundance > 0) / n_samples * 100,
      .groups = "drop"
    ) %>%
    # Remove near-absent taxa (<0.01%) to keep the chart light
    filter(MeanAbund >= 0.01) %>%
    mutate(Quadrant = case_when(
      Prevalence >= seuil_prev & MeanAbund >= seuil_abond ~ "1. High Abundance Shared Microbiota",
      Prevalence >= seuil_prev & MeanAbund <  seuil_abond ~ "2. Low Abundance Shared Microbiota",
      Prevalence <  seuil_prev & MeanAbund >= seuil_abond ~ "3. Transient",
      Prevalence <  seuil_prev & MeanAbund <  seuil_abond ~ "4. Background Noise"
    ))
}


# ==============================================================================
# BLOCK 2: DATA LOADING AND FILTERING
# ==============================================================================
message("--- Chargement des données ---")

# --- Translation dictionaries ---
# MATAM and the metadata use short French codes.
# They are translated into full English names with weights (mg) so that the
# charts are directly usable for a publication.
trad_conditions <- c(
  "Pre_MicroLarve" = "Tiny_Larvae (2 mg)",
  "MicroLarve"     = "Microlarvae (7 mg)",
  "Larve_S1"       = "Larvae_S1 (14 mg)",
  "Larve_S2"       = "Larvae_S2 (40 mg)",
  "Larve_S3"       = "Larvae_S3 (65 mg)",
  "Larve_S4"       = "Larvae_S4 (100 mg)",
  "Nymphe"         = "Pupae",
  "Jeune_Adulte"   = "Young_Beetles",
  "Adulte"         = "Beetles",
  "Substrat_Base"  = "Raw_Substrate",
  "Substrat_Final" = "Final_Substrate"
)

trad_samples <- c(
  "Pre1" = "TL1", "Pre2" = "TL2", "Pre3" = "TL3", "Pre4" = "TL4",
  "Mi1"  = "Mi1", "Mi2"  = "Mi2", "Mi3"  = "Mi3", "Mi4"  = "Mi4",
  "S1-L1" = "S1-L1", "S1-L2" = "S1-L2", "S1-L3" = "S1-L3", "S1-L4" = "S1-L4",
  "S2-L1" = "S2-L1", "S2-L2" = "S2-L2", "S2-L3" = "S2-L3", "S2-L4" = "S2-L4",
  "S3-L1" = "S3-L1", "S3-L2" = "S3-L2", "S3-L3" = "S3-L3", "S3-L4" = "S3-L4",
  "S4-L1" = "S4-L1", "S4-L2" = "S4-L2", "S4-L3" = "S4-L3", "S4-L4" = "S4-L4",
  "N1"  = "P1",  "N2"  = "P2",  "N3"  = "P3",  "N4"  = "P4",
  "JA1" = "YB1", "JA2" = "YB2", "JA3" = "YB3", "JA4" = "YB4",
  "A1"  = "B1",  "A2"  = "B2",  "A3"  = "B3",  "A4"  = "B4",
  "SB1" = "RB1", "SB2" = "RB2", "SB3" = "RB3", "SB4" = "RB4",
  "SF1" = "FS1", "SF2" = "FS2", "SF3" = "FS3", "SF4" = "FS4"
)

# --- Reading the count table ---
# rename_with() renames the sample columns using the trad_samples dictionary
# -Taxonomy: the Taxonomy column is excluded from renaming
df_qiime_16s <- lire_qiime_tsv_robuste(fichier_qiime) %>%
  rename_with(~ recode(.x, !!!trad_samples), -Taxonomy)

# --- Reading the metadata ---
meta <- read_tsv(fichier_metadata, show_col_types = FALSE) %>%
  rename(Sample = `sample-id`, Condition = condition) %>%
  mutate(
    Condition = recode(Condition, !!!trad_conditions),
    Sample    = recode(Sample,    !!!trad_samples)
  )

# --- Extracting taxonomic ranks via regex ---
# The Taxonomy column contains strings of the form:
#   "d__Bacteria;p__Firmicutes;c__Bacilli;...;g__Lactobacillus;s__acidophilus"
# str_extract() captures everything following the prefix up to the next ";"
df_taxo_16s <- df_qiime_16s %>%
  # Switch from wide format (1 column/sample) to long format (1 row/taxon x sample)
  pivot_longer(cols = -Taxonomy, names_to = "Sample", values_to = "Count") %>%

  # Optimization: keep only rows with at least 1 detected read
  filter(Count > 0) %>%

  # Crucial biological filtering: chloroplasts and mitochondria carry
  # a 16S rRNA gene homologous to bacteria. They systematically contaminate
  # plant or animal DNA extractions and do not represent the insect's actual
  # gut microbiome -> they are removed.
  filter(!grepl("Chloroplast|Mitochondria|Mitochondrion", Taxonomy, ignore.case = TRUE)) %>%

  # Extracting each taxonomic rank via its SILVA prefix
  mutate(
    Phylum  = str_extract(Taxonomy, "p__[^;]+"),
    Class   = str_extract(Taxonomy, "c__[^;]+"),
    Order   = str_extract(Taxonomy, "o__[^;]+"),
    Family  = str_extract(Taxonomy, "f__[^;]+"),
    Genus   = str_extract(Taxonomy, "g__[^;]+"),
    Species = str_extract(Taxonomy, "s__[^;]+")
  ) %>%

  # Cleaning up prefixes ("g__", "p__", etc.)
  mutate(across(c(Phylum, Class, Order, Family, Genus, Species),
                ~ str_remove(.x, "^[a-z]__"))) %>%

  # Replacing NA with "Unassigned": preferable to NA for filtering
  mutate(across(c(Phylum, Class, Order, Family, Genus, Species),
                ~ ifelse(is.na(.x) | str_trim(.x) == "", "Unassigned", .x)))

# --- Assembling the final table ---
df_final_16s <- harmoniser_taxonomie(df_taxo_16s) %>%
  left_join(meta, by = "Sample") %>%
  filter(!is.na(Condition))   # Discard samples without metadata

# --- Chronological order of T. molitor development ---
# THIS ORDER IS CRUCIAL: it governs the X axis of all time-course charts.
# T. molitor is a holometabolous insect (complete metamorphosis):
#   Egg -> Larvae (several stages) -> Pupa -> Adult
chronologie_insecte <- c(
  "Tiny_Larvae (2 mg)", "Microlarvae (7 mg)",
  "Larvae_S1 (14 mg)", "Larvae_S2 (40 mg)", "Larvae_S3 (65 mg)", "Larvae_S4 (100 mg)",
  "Pupae", "Young_Beetles", "Beetles"
)

# Extended version including the feed substrates (environment)
chronologie_tout <- c(chronologie_insecte, "Raw_Substrate", "Final_Substrate")

# Table filtered on all conditions (insect + substrates)
# NOTE: JA1 and N3 (pathological samples) are intentionally kept.
df_insecte <- df_final_16s %>%
  filter(Condition %in% chronologie_tout) %>%
  mutate(
    Condition    = factor(Condition, levels = chronologie_tout),
    # ConditionNum: numeric version of the factor for ggplot2's continuous X axis
    ConditionNum = as.numeric(Condition)
  )

# Color palette by developmental stage
# Gradient from light green (small larvae) to dark blue (older larvae) -> purple (pupa) -> red (adult)
my_cond_colors <- setNames(
  c(
    "#d8f5c7ff", "#8dd1c1ff", "#41B6C4", "#0fa5e0ff", "#0f6cddff", "#001858ff",
    "#984EA3",
    "#FD8D3C", "#E31A1C",
    "#FFD92F", "#8C510A"
  ),
  chronologie_tout
)


# ==============================================================================
# BLOCK 3: ALPHA AND BETA DIVERSITY CALCULATIONS
# ==============================================================================

# ==============================================================================
# A. ALPHA DIVERSITY: SHANNON INDEX
# ==============================================================================
# ALPHA DIVERSITY measures diversity WITHIN each sample.
# The Shannon index (H') accounts for both the number of species
# (richness) and their evenness (regularity of the distribution).
# High H' = diverse and balanced community
# Low H'  = community dominated by a few taxa
message("--- Calcul de l'alpha-diversité (Shannon) ---")

# Building the Samples x Genera matrix (format required by vegan)
mat_counts_alpha <- df_insecte %>%
  group_by(Sample, Genus) %>%
  summarise(Count = sum(Count), .groups = "drop") %>%
  pivot_wider(names_from = Genus, values_from = Count, values_fill = 0) %>%
  column_to_rownames("Sample")

# Computing Shannon on RAW COUNTS (not the normalized abundances)
# vegan's diversity() function handles normalization internally
shannon_vals <- diversity(mat_counts_alpha, index = "shannon")

# Summary by developmental stage (mean +/- standard deviation)
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
    SD_Shannon  = ifelse(n() > 1, sd(Shannon), 0),   # SD = 0 if only 1 replicate
    .groups = "drop"
  )


# ==============================================================================
# B. BETA DIVERSITY: BRAY-CURTIS DISTANCE
# ==============================================================================
# BETA DIVERSITY measures the DISSIMILARITY between samples.
# The Bray-Curtis distance ranges from 0 (identical communities) to 1
# (no shared species). Here it is computed BETWEEN REPLICATES of the same
# stage to measure within-group heterogeneity: if the replicates of a stage
# are very different from one another, the within-group beta diversity will be high.
message("--- Calcul de la bêta-diversité (Bray-Curtis) ---")

resultats_stade <- list()   # Aggregated results per stage (mean +/- SD)
resultats_ind   <- list()   # Results per individual replicate

for (cond in levels(df_insecte$Condition)) {
  df_cond <- df_insecte %>% filter(Condition == cond)
  samps   <- unique(df_cond$Sample)

  # At least 2 samples are needed to compute a pairwise distance
  if (length(samps) < 2) next

  # Normalization to 100% ("Unassigned" entries are removed before
  # normalization so they do not bias the relative composition)
  df_norm <- normaliser_100(df_cond %>% filter(Genus != "Unassigned"), "Genus")

  # Conversion to a Samples x Genera matrix (vegan format)
  mat_wide <- df_norm %>%
    pivot_wider(id_cols = Sample, names_from = Genus,
                values_from = Relative_Abundance, values_fill = 0) %>%
    column_to_rownames("Sample")

  # Computing the pairwise distance matrix (all replicates, pairwise)
  dist_mat <- as.matrix(vegdist(mat_wide, method = "bray"))

  # upper.tri(): only the upper triangle is kept to avoid duplicates
  vals <- dist_mat[upper.tri(dist_mat)]

  resultats_stade[[cond]] <- data.frame(
    Condition = cond,
    Beta_Mean = mean(vals),
    Beta_SD   = ifelse(length(vals) > 1, sd(vals), 0)
  )

  # For each replicate: mean distance to the other replicates of the same stage
  diag(dist_mat) <- NA   # The diagonal = 0 (distance to itself) -> excluded
  beta_ind <- rowMeans(dist_mat, na.rm = TRUE)
  resultats_ind[[cond]] <- data.frame(
    Sample    = names(beta_ind),
    Condition = cond,
    Beta_Ind  = beta_ind
  )
}

df_beta_stade <- bind_rows(resultats_stade)
df_beta_ind   <- bind_rows(resultats_ind)


# ==============================================================================
# C. SPEARMAN CORRELATION TEST: SHANNON vs BRAY-CURTIS
# ==============================================================================
# We test whether the most diverse samples (high Shannon) are also
# the most homogeneous between replicates (low Bray-Curtis), or the opposite.
# The Spearman test is non-parametric -> suited to non-normal distributions
# such as microbiome data.
message("--- Test de corrélation de Spearman ---")

# Individual table: 1 row per replicate with both Shannon AND Bray-Curtis
df_ind <- data.frame(Sample = names(shannon_vals), Shannon = shannon_vals) %>%
  left_join(df_beta_ind, by = "Sample") %>%
  filter(!is.na(Beta_Ind))

# Test 1: WITHOUT the substrates (insect only)
df_ind_sans   <- df_ind %>% filter(!grepl("Substrat|Substrate", Condition))
test_stat_sans <- cor.test(df_ind_sans$Shannon, df_ind_sans$Beta_Ind, method = "spearman", exact = FALSE)
rho_sans  <- round(test_stat_sans$estimate, 2)
p_val_sans <- paste0("= ", format(test_stat_sans$p.value, scientific = TRUE, digits = 3))

# Test 2: WITH the substrates (full dataset)
test_stat_avec <- cor.test(df_ind$Shannon, df_ind$Beta_Ind, method = "spearman", exact = FALSE)
rho_avec  <- round(test_stat_avec$estimate, 2)
p_val_avec <- paste0("= ", format(test_stat_avec$p.value, scientific = TRUE, digits = 3))

# Merging alpha + beta by stage for the charts
df_final <- df_alpha_stade %>%
  left_join(df_beta_stade, by = "Condition") %>%
  mutate(
    # Trajectoire: distinguishes insect vs substrates for chart 2's facet
    Trajectoire = ifelse(grepl("Substrat|Substrate", Condition), "Environnement", "Développement de l'Insecte"),
    Trajectoire = factor(Trajectoire, levels = c("Développement de l'Insecte", "Environnement"))
  )


# ==============================================================================
# BLOCK 4: EXPORTING RESULTS AS TSV
# ==============================================================================
# A summary table with the raw numeric values is exported.
# Useful for: traceability, future re-analysis, or supplementary material.
message("--- Export des valeurs numériques ---")

df_final %>%
  select(Stade = Condition, MeanShannon, SD_Shannon, Beta_Mean, Beta_SD) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  write.table(
    file = file.path(out_dir, "Valeurs_Shannon_BrayCurtis.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
message("=> Fichier TSV des valeurs exporté.")


# ==============================================================================
# BLOCK 5: DUAL-AXIS CHART (SHANNON + BRAY-CURTIS)
# ==============================================================================
# This chart overlays two metrics on two separate Y axes.
# The coefficient (coeff) trick allows aligning both curves
# on the same visual scale while keeping independent axes.
#
# Coefficient = max(Shannon) / max(Bray-Curtis)
# -> All Bray-Curtis values are multiplied by this coefficient
#   before being plotted on the Shannon Y axis.
# -> The right Y axis is rescaled by the inverse (division) to display
#   the true Bray-Curtis values.
message("--- Génération des Graphiques Alpha/Bêta ---")

# Internal function to create a dual-axis chart
creer_graphique_double_axe <- function(df_data, coeff, chronologie_labels,
                                        rho, p_val, titre, titre_x,
                                        df_facette = NULL) {
  # Position of the statistical annotation (bottom-left corner)
  df_label <- data.frame(
    ConditionNum = 1,
    y_pos = max(df_data$MeanShannon, na.rm = TRUE) * 0.1,
    label = paste0("rho = ", rho, "\np ", p_val)
  )
  if (!is.null(df_facette)) df_label$Trajectoire <- df_facette

  p <- ggplot(df_data, aes(x = ConditionNum)) +

    # Error bars (standard deviation) for both metrics
    geom_errorbar(aes(ymin = (Beta_Mean - Beta_SD) * coeff,
                      ymax = (Beta_Mean + Beta_SD) * coeff),
                  color = COL_BRAY, width = 0.2, linewidth = 0.8) +
    geom_errorbar(aes(ymin = MeanShannon - SD_Shannon,
                      ymax = MeanShannon + SD_Shannon),
                  color = COL_SHANNON, width = 0.2, linewidth = 0.8) +

    # Bray-Curtis curve (dashed + triangles)
    geom_line(aes(y = Beta_Mean * coeff, color = "Heterogeneity (Bray-Curtis)"),
              linewidth = 1.5, linetype = "dashed") +
    geom_point(aes(y = Beta_Mean * coeff, color = "Heterogeneity (Bray-Curtis)"),
               size = 4, shape = 17) +

    # Shannon curve (solid + circles)
    geom_line(aes(y = MeanShannon, color = "Diversity (Shannon)"),
              linewidth = 1.5) +
    geom_point(aes(y = MeanShannon, color = "Diversity (Shannon)"),
               size = 5, shape = 16) +

    # Spearman correlation test annotation
    geom_label(data = df_label,
               aes(x = ConditionNum, y = y_pos, label = label),
               fill = "white", fontface = "bold", size = 3.8,
               hjust = 0, color = COL_BRAY, inherit.aes = FALSE) +

    # X axis: developmental stages
    scale_x_continuous(breaks = seq_along(chronologie_labels),
                       labels = chronologie_labels) +

    # Dual Y axis: Shannon (left) and rescaled Bray-Curtis (right)
    scale_y_continuous(
      name = "Shannon index",
      sec.axis = sec_axis(~ . / coeff, name = "Bray-Curtis dissimilarity")
    ) +

    scale_color_manual(values = c(
      "Diversity (Shannon)"          = COL_SHANNON,
      "Heterogeneity (Bray-Curtis)"  = COL_BRAY
    )) +

    labs(title = titre, x = titre_x, color = NULL) +
    theme_bw(base_size = 13) +
    theme(
      axis.text.x        = element_text(angle = 50, hjust = 1, face = "bold"),
      axis.title.y.left  = element_text(color = COL_SHANNON, face = "bold"),
      axis.title.y.right = element_text(color = COL_BRAY, face = "bold"),
      legend.position    = "bottom",
      panel.grid.minor   = element_blank()
    )

  return(p)
}

# --- Chart 1: Insect only (without substrates) ---
df_sans    <- df_final %>% filter(Trajectoire == "Développement de l'Insecte")
coeff_sans <- max(df_sans$MeanShannon, na.rm = TRUE) / max(df_sans$Beta_Mean, na.rm = TRUE)

p_dyn_sans <- creer_graphique_double_axe(
  df_data           = df_sans,
  coeff             = coeff_sans,
  chronologie_labels = chronologie_insecte,
  rho    = rho_sans,
  p_val  = p_val_sans,
  titre  = "Alpha vs Beta diversity ",
  titre_x = "Stages of development"
)
ggsave(file.path(out_dir, "shannon_vs_braycurtis.png"),
       plot = p_dyn_sans, width = 12, height = 8, dpi = 300)

ggsave(filename = file.path(out_dir, "shannon_vs_braycurtis.svg"),
       plot = p_dyn_sans, 
       width = 18, 
       height = 8)

# --- Chart 2: Insect + Substrates (with facets) ---
coeff_avec <- max(df_final$MeanShannon, na.rm = TRUE) / max(df_final$Beta_Mean, na.rm = TRUE)

p_dyn_avec <- creer_graphique_double_axe(
  df_data           = df_final,
  coeff             = coeff_avec,
  chronologie_labels = chronologie_tout,
  rho    = rho_avec,
  p_val  = p_val_avec,
  titre  = "Dynamique du Microbiote (Insecte et Environnement)",
  titre_x = "Stades de développement et Environnement",
  df_facette = factor("Développement de l'Insecte",
                      levels = c("Développement de l'Insecte", "Environnement"))
) +
  facet_grid(~ Trajectoire, scales = "free_x", space = "free_x") +
  theme(
    strip.text.x     = element_text(face = "bold", size = 11, color = "white"),
    strip.background = element_rect(fill = "grey30", color = "black")
  )

#ggsave(file.path(out_dir, "shannon_vs_braycurtis_substrats.png"),
#       plot = p_dyn_avec, width = 12, height = 8, dpi = 300)


# ==============================================================================
# BLOCK 6: GLOBAL 16S BARPLOT — REPLICATE COMPOSITION (GENUS LEVEL)
# ==============================================================================
# This chart displays the bacterial genus composition of EACH replicate,
# organized by development stage. It allows evaluating reproducibility between
# replicates and community evolution across development.
message("--- Global 16S Barplot ---")

# Define chronological order excluding substrates
ordre_chronologique <- c(
  "Tiny_Larvae (2 mg)", "Microlarvae (7 mg)",
  "Larvae_S1 (14 mg)", "Larvae_S2 (40 mg)", "Larvae_S3 (65 mg)", "Larvae_S4 (100 mg)",
  "Pupae", "Young_Beetles", "Beetles"
)

# Normalize data and filter out substrates immediately
df_global_norm <- normaliser_100(df_final_16s, "Genus") %>%
  left_join(meta %>% select(Sample, Condition), by = "Sample") %>%
  filter(!is.na(Condition)) %>%
  filter(Condition %in% ordre_chronologique) %>% # Exclude substrates
  mutate(Condition = factor(Condition, levels = ordre_chronologique))

# Identify Top 14 genera globally
top_global <- df_global_norm %>%
  filter(!(Genus %in% c("Unassigned", "Bacteroides", "Alistipes", 
                         "Aquabacterium", "Blastococcus", "Sphingomonas", "Delftia", "Listeria"))) %>%
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

# Color palette
levels_taxa_global <- levels(df_plot_global$Taxon_Top)
taxa_vrais_global  <- setdiff(levels_taxa_global, c("Others", "Unassigned"))
my_col_global <- setNames(rep_len(pal_taxo, length(taxa_vrais_global)), taxa_vrais_global)
if ("Others"     %in% levels_taxa_global) my_col_global["Others"]     <- "grey85"
if ("Unassigned" %in% levels_taxa_global) my_col_global["Unassigned"] <- "grey40"

p_global <- ggplot(df_plot_global, aes(x = Sample, y = Abundance, fill = Taxon_Top)) +
  geom_bar(stat = "identity", color = "black",
           linewidth = 0.2, width = 0.9, alpha = 0.9) +
  scale_fill_manual(values = my_col_global) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.02))
  ) +
  facet_grid(~ Condition, scales = "free_x", space = "free_x") +
  labs(
    title = "Global Microbiota Composition (16S)",
    x = "Biological Replicates", 
    y = "Relative Abundance (%)",
    fill = "Bacterial Genus"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, face = "bold",
                                    size = 8, color = "grey20"),
    strip.text.x     = element_text(face = "bold", size = 10, color = "white"),
    strip.background = element_rect(fill = "grey30", color = "black"),
    legend.text      = element_text(face = "italic", size = 10),
    legend.title     = element_text(face = "bold"),
    legend.position  = "bottom",
    panel.spacing    = unit(0.2, "lines"),
    panel.grid.major.x = element_blank()
  )

ggsave(filename = file.path(out_dir, "Barplot_Global_16S.png"),
       plot = p_global, width = 18, height = 8, dpi = 300)

ggsave(filename = file.path(out_dir, "Barplot_Global_16S.svg"),
       plot = p_global, 
       width = 18, 
       height = 8)

# ==============================================================================
# BLOCK 7: BRAY-CURTIS DISSIMILARITY HEATMAP (ALL REPLICATES)
# ==============================================================================
# This matrix shows the Bray-Curtis dissimilarity between EACH pair of
# replicates (across all conditions). Diagonal blocks (green cells)
# indicate that the replicates of a given stage are similar to one another.
# Off-diagonal green blocks would indicate stages sharing a
# common microbial composition.
message("--- Heatmap Bray-Curtis globale ---")

# Internal function to create a distance heatmap
creer_heatmap_bray <- function(df_dist_long, ordre, titre, largeur = 12, hauteur = 10,
                                nom_fichier) {
  df_plot <- df_dist_long %>%
    filter(Sample1 %in% ordre & Sample2 %in% ordre) %>%
    mutate(
      Sample1 = factor(Sample1, levels = ordre),
      # Sample2's order is reversed so the diagonal runs
      # from top-left to bottom-right (standard matrix convention)
      Sample2 = factor(Sample2, levels = rev(ordre))
    )

  p <- ggplot(df_plot, aes(x = Sample1, y = Sample2, fill = Distance)) +
    geom_tile(color = "white", linewidth = 0.3) +
    # Semantic gradient: Green (0 = identical) -> Yellow -> Red (1 = different)
    scale_fill_gradientn(
      colors = c("#2e7d32", "#fff59d", "#c62828"),
      limits = c(0, 1),
      name   = "Bray-Curtis\ndissimilarity"
    ) +
    labs(title = titre, x = "Biological Replicates", y = "Biological Replicates") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                  size = 9, face = "bold"),
      axis.text.y  = element_text(size = 9, face = "bold"),
      panel.grid   = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      plot.title   = element_text(face = "bold", size = 14)
    )

  ggsave(filename = file.path(out_dir, nom_fichier),
         plot = p, width = largeur, height = hauteur, dpi = 300)
}

# Preparing the full distance matrix
df_matrice <- df_final_16s %>%
  filter(Genus != "Unassigned") %>%
  group_by(Sample, Genus) %>%
  summarise(Count = sum(Count), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(Relative_Abundance = (Count / sum(Count)) * 100) %>%
  ungroup()

mat_bray_complete  <- df_matrice %>%
  pivot_wider(id_cols = Sample, names_from = Genus,
              values_from = Relative_Abundance, values_fill = 0) %>%
  column_to_rownames("Sample")

dist_mat_complete <- as.matrix(vegdist(mat_bray_complete, method = "bray"))

# Switching to long format for ggplot2
df_dist_long <- as.data.frame(dist_mat_complete) %>%
  rownames_to_column("Sample1") %>%
  pivot_longer(cols = -Sample1, names_to = "Sample2", values_to = "Distance")

# Chronological ordering of replicates
ordre_echantillons <- df_global_norm %>%
  distinct(Sample, Condition) %>%
  arrange(Condition, Sample) %>%
  pull(Sample)

# Heatmap 1: All replicates (insect + substrates)
creer_heatmap_bray(
  df_dist_long, ordre_echantillons,
  titre = "Inter-replicates dissimilarity matrix (Bray-Curtis)",
  nom_fichier = "Matrice_BrayCurtis_Tous_Replicats.png"
)

# Heatmap 2: Insect only (without substrates)
samples_insecte <- meta %>% filter(Condition %in% chronologie_insecte) %>% pull(Sample)
ordre_insecte   <- ordre_echantillons[ordre_echantillons %in% samples_insecte]

creer_heatmap_bray(
  df_dist_long, ordre_insecte,
  titre = "Inter-replicates dissimilarity matrix (Insect only)",
  largeur = 10, hauteur = 8,
  nom_fichier = "Matrice_BrayCurtis_Insecte_Seul.png"
)


# ==============================================================================
# BLOCK 8: MASTER FIGURE
# ==============================================================================
# This figure combines into a single chart:
#   - A header with the stage and trajectory names
#   - The Shannon + Bray-Curtis curves (middle)
#   - The taxonomic composition barplot (bottom)
# All three panels share the same finely-tuned X axis, aligned to the column
# of each replicate. Assembly is done using the patchwork package.
message("--- Génération de la Figure Maîtresse ---")

generer_master_ultime <- function(stades_filtre, nom_fichier) {

  # --- 1. FINE-GRAINED POSITIONING OF REPLICATES ---
  # A numeric X position is computed for each replicate, with wider gaps
  # between the major categories (insect vs environment)
  df_samples <- meta %>%
    filter(Condition %in% stades_filtre) %>%
    mutate(Condition = factor(Condition, levels = stades_filtre)) %>%
    arrange(Condition, Sample) %>%
    mutate(
      Trajectoire = ifelse(grepl("Substrat|Substrate", Condition),
                           "Environment", "Insect development"),
      Trajectoire = factor(Trajectoire, levels = c("Insect development", "Environment"))
    )

  df_samples$x_pos <- 0
  cur_x <- 1
  last_cond <- df_samples$Condition[1]

  for(i in 1:nrow(df_samples)) {
    if(df_samples$Condition[i] != last_cond) {
      # Wide gap (1.5) between trajectories, normal gap (0.8) between stages
      if(df_samples$Trajectoire[i] != df_samples$Trajectoire[i-1]) {
        cur_x <- cur_x + 1.5
      } else {
        cur_x <- cur_x + 0.8
      }
    }
    df_samples$x_pos[i] <- cur_x
    cur_x     <- cur_x + 1
    last_cond <- df_samples$Condition[i]
  }

  # --- 2. COMPUTING COLUMNS AND BANDS ---
  df_conds <- df_samples %>%
    group_by(Condition, Trajectoire) %>%
    summarise(x_mid = mean(x_pos),
              x_min = min(x_pos) - 0.48,
              x_max = max(x_pos) + 0.48, .groups = "drop")

  # Gaps between condition groups (grey separator zones)
  df_gaps <- data.frame(
    xmin = head(df_conds$x_max, -1),
    xmax = tail(df_conds$x_min, -1)
  )

  # Trajectory bands (span several conditions)
  df_traj <- df_conds %>%
    group_by(Trajectoire) %>%
    summarise(x_min = min(x_min), x_max = max(x_max),
              x_mid = mean(c(min(x_min), max(x_max))), .groups = "drop")

  LIMITES_X <- c(min(df_conds$x_min) - 0.1, max(df_conds$x_max) + 0.1)

  # --- 3. CORRELATION STATISTICS ON THIS SUBSET ---
  df_stats_loc <- df_ind %>%
    filter(Sample %in% df_samples$Sample & !is.na(Shannon) & !is.na(Beta_Ind))
  test_loc <- cor.test(df_stats_loc$Shannon, df_stats_loc$Beta_Ind,
                        method = "spearman", exact = FALSE)
  rho_txt <- round(test_loc$estimate, 2)
  p_txt   <- format(test_loc$p.value, scientific = TRUE, digits = 3)

  df_label <- data.frame(
    x = df_conds$x_mid[1],
    y = max(df_final$MeanShannon[df_final$Condition %in% stades_filtre], na.rm = TRUE) * 0.1,
    label = paste0("rho = ", rho_txt, "\np = ", p_txt)
  )

  coeff_loc <- max(df_final$MeanShannon[df_final$Condition %in% stades_filtre], na.rm = TRUE) /
               max(df_final$Beta_Mean[df_final$Condition %in% stades_filtre], na.rm = TRUE)

  df_top <- df_final %>%
    filter(Condition %in% stades_filtre) %>%
    select(-any_of("Trajectoire")) %>%
    left_join(df_conds, by = "Condition")

  df_bot <- df_plot_global %>%
    filter(Sample %in% df_samples$Sample) %>%
    select(-any_of("Trajectoire")) %>%
    left_join(df_samples, by = c("Sample", "Condition"))

  # --- 4. HEADER PANEL (Trajectory + Condition bands) ---
  p_labels <- ggplot() +
    geom_rect(data = df_gaps,
              aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 2),
              fill = "grey92", inherit.aes = FALSE) +
    geom_rect(data = df_traj,
              aes(xmin = x_min, xmax = x_max, ymin = 1, ymax = 2),
              fill = "#2C3E50", color = "black", linewidth = 0.3) +
    geom_text(data = df_traj, aes(x = x_mid, y = 1.5, label = Trajectoire),
              color = "white", fontface = "bold", size = 4) +
    geom_rect(data = df_conds,
              aes(xmin = x_min, xmax = x_max, ymin = 0, ymax = 1),
              fill = "grey35", color = "black", linewidth = 0.3) +
    geom_text(data = df_conds, aes(x = x_mid, y = 0.5, label = Condition),
              color = "white", fontface = "bold", size = 3) +
    scale_x_continuous(limits = LIMITES_X, expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 2),   expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(t = 10, r = 10, b = 2, l = 10))

  # --- 5. MIDDLE PANEL (Shannon + Bray-Curtis curves) ---
  p_top_plot <- ggplot(df_top, aes(x = x_mid)) +
    geom_rect(data = df_gaps, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "grey96", inherit.aes = FALSE) +
    geom_vline(data = df_conds, aes(xintercept = x_min),
               color = "grey50", linetype = "dashed", linewidth = 0.5) +
    geom_vline(data = df_conds, aes(xintercept = x_max),
               color = "grey50", linetype = "dashed", linewidth = 0.5) +
    geom_errorbar(aes(ymin = (Beta_Mean - Beta_SD) * coeff_loc,
                      ymax = (Beta_Mean + Beta_SD) * coeff_loc),
                  color = COL_BRAY, width = 0.5) +
    geom_errorbar(aes(ymin = MeanShannon - SD_Shannon,
                      ymax = MeanShannon + SD_Shannon),
                  color = COL_SHANNON, width = 0.5) +
    geom_line(aes(y = Beta_Mean * coeff_loc, group = Trajectoire,
                  color = "Heterogeneity (Bray-Curtis)",
                  linetype = "Heterogeneity (Bray-Curtis)"), linewidth = 1.3) +
    geom_line(aes(y = MeanShannon, group = Trajectoire,
                  color = "Diversity (Shannon)",
                  linetype = "Diversity (Shannon)"), linewidth = 1.3) +
    geom_point(aes(y = Beta_Mean * coeff_loc,
                   color = "Heterogeneity (Bray-Curtis)",
                   shape = "Heterogeneity (Bray-Curtis)"), size = 4) +
    geom_point(aes(y = MeanShannon,
                   color = "Diversity (Shannon)",
                   shape = "Diversity (Shannon)"), size = 5) +
    geom_label(data = df_label, aes(x = x, y = y, label = label),
               fill = "white", fontface = "bold", size = 3.5,
               color = COL_BRAY, hjust = 0, inherit.aes = FALSE) +
    scale_x_continuous(limits = LIMITES_X, expand = c(0, 0)) +
    scale_y_continuous(
      name = "Shannon index",
      sec.axis = sec_axis(~ . / coeff_loc, name = "Bray-Curtis dissimilarity")
    ) +
    scale_color_manual(
      values = c("Diversity (Shannon)" = COL_SHANNON,
                 "Heterogeneity (Bray-Curtis)" = COL_BRAY), name = "") +
    scale_linetype_manual(
      values = c("Diversity (Shannon)" = "solid",
                 "Heterogeneity (Bray-Curtis)" = "dashed"), name = "") +
    scale_shape_manual(
      values = c("Diversity (Shannon)" = 16,
                 "Heterogeneity (Bray-Curtis)" = 17), name = "") +
    theme_bw(base_size = 13) +
    theme(
      axis.title.x = element_blank(), axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y.left  = element_text(color = COL_SHANNON, face = "bold"),
      axis.title.y.right = element_text(color = COL_BRAY, face = "bold"),
      legend.position = "top", panel.grid = element_blank(),
      plot.margin = margin(t = 0, r = 10, b = 0, l = 10)
    )

  # --- 6. BOTTOM PANEL (Taxonomic composition barplots) ---
  p_bot_plot <- ggplot(df_bot, aes(x = x_pos, y = Abundance)) +
    geom_rect(data = df_gaps, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "grey96", inherit.aes = FALSE) +
    geom_vline(data = df_conds, aes(xintercept = x_min),
               color = "grey50", linetype = "dashed", linewidth = 0.5) +
    geom_vline(data = df_conds, aes(xintercept = x_max),
               color = "grey50", linetype = "dashed", linewidth = 0.5) +
    geom_bar(aes(fill = Taxon_Top), stat = "identity",
             color = "black", linewidth = 0.15, width = 0.9) +
    scale_fill_manual(values = my_col_global, name = "Bacterial Genus :") +
    scale_x_continuous(breaks = df_samples$x_pos, labels = df_samples$Sample,
                       limits = LIMITES_X, expand = c(0, 0)) +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       name = "Abundance (%)", expand = c(0, 0)) +
    theme_bw(base_size = 13) +
    theme(
      axis.text.x    = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                    face = "bold", size = 8),
      axis.title.x   = element_blank(),
      panel.grid     = element_blank(),
      legend.position = "bottom",
      legend.text    = element_text(face = "italic"),
      plot.margin    = margin(t = 0, r = 10, b = 10, l = 10)
    )

  # --- 7. ASSEMBLY WITH PATCHWORK ---
  # The three panels are stacked vertically.
  # heights = c(0.12, 1, 1.25): the header is very compact,
  # the barplot slightly taller than the curves.
  final_p <- (p_labels / p_top_plot / p_bot_plot) +
    plot_layout(heights = c(0.12, 1, 1.25), guides = "collect") &
    theme(legend.position = "bottom", legend.box = "vertical")

  ggsave(file.path(out_dir, nom_fichier),
         final_p, width = 18, height = 13, dpi = 600, bg = "white")
  nom_svg <- str_replace(nom_fichier, "\\.png$", ".svg")
  ggsave(file.path(out_dir, nom_svg),
         final_p, width = 18, height = 13, bg = "white")
  message(sprintf("  => Master figure '%s' sauvegardée.", nom_fichier))
}

generer_master_ultime(chronologie_tout,    "Master_Figure_substrats.png")
generer_master_ultime(chronologie_insecte, "Master_Figure.png")

# ==============================================================================
# BLOCK 9: LARVAL CORE MICROBIOME (GENUS + SPECIES BARPLOTS)
# ==============================================================================
# The "core microbiome" of a group of stages is the set of bacteria
# regularly and abundantly present at those stages.
# It is represented here as a mean per stage (barplot with X axis = stage).
message("--- Core Microbiote Larvaire ---")

# Larval stages of interest (very young pre-larvae and pupae are excluded)
cond_larves <- c("Microlarvae (7 mg)", "Larvae_S1 (14 mg)", "Larvae_S2 (40 mg)",
                  "Larvae_S3 (65 mg)", "Larvae_S4 (100 mg)")
df_larves_brut <- df_final_16s %>% filter(Condition %in% cond_larves)

# --- Part A: Genus level ---
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
ggsave(file.path(out_dir, "Larval_Microbiome_Means.png"),
       plot = p_core, width = 9, height = 7, dpi = 300)

# --- Part B: Species level ---
df_larves_sp <- forcer_nomenclature_binomiale(df_larves_brut)
res_larves_sp <- preparer_barplot_data(
  df_brut = df_larves_sp, rang = "Species",
  top_n = 13, seuil_min = 1.0, cond_levels = cond_larves
)
p_core_sp <- creer_barplot(
  res_larves_sp$df_plot, res_larves_sp$palette,
  titre    = "Larval Microbiome (16S - Species level)",
  x_lab    = "Larval Stages",
  fill_lab = "Bacterial Species:"
)
ggsave(file.path(out_dir, "Larval_Microbiome_Means_Species.png"),
       plot = p_core_sp, width = 9.5, height = 7, dpi = 300)
message("=> Larval Microbiome barplots générés.")


# ==============================================================================
# BLOCK 10: QUADRANT CHARTS (PREVALENCE vs ABUNDANCE) — LARVAE
# ==============================================================================
# This chart places each genus (or species) in a two-dimensional plane:
#   - X axis: Prevalence (% of larvae in which this genus is detected)
#   - Y axis: Mean relative abundance (log scale)
# The threshold lines divide the space into 4 ecological quadrants.
message("--- Graphiques Quadrants Larves ---")

SEUIL_PREVALENCE <- 80   # Present in >= 80% of larvae -> "shared"
SEUIL_ABONDANCE  <- 1.0  # Represents >= 1% of the community -> "abundant"

# --- Genus ---
df_larves_norm_gen <- normaliser_100(df_larves_brut, "Genus") %>%
  left_join(meta %>% select(Sample, Condition), by = "Sample")
n_larves <- n_distinct(df_larves_norm_gen$Sample)

df_core_scatter_gen <- calculer_core_scatter(
  df_larves_norm_gen, "Genus", n_larves, SEUIL_PREVALENCE, SEUIL_ABONDANCE
)
p_scatter_gen <- creer_quadrant_plot(
  df_core_scatter_gen, "Genus", SEUIL_PREVALENCE, SEUIL_ABONDANCE,
  titre = "Larval Microbiome Structure (16S - Genus level)"
)
ggsave(file.path(out_dir, "Larval_Microbiome_Quadrants.png"),
       plot = p_scatter_gen, width = 11, height = 8, dpi = 300)
ggsave(file.path(out_dir, "Larval_Microbiome_Quadrants.svg"),
       plot = p_scatter_gen, width = 11, height = 8)
# --- Species ---
df_larves_sp_quad <- forcer_nomenclature_binomiale(df_larves_brut)
df_larves_norm_sp_quad <- normaliser_100(df_larves_sp_quad, "Species") %>%
  left_join(meta %>% select(Sample, Condition), by = "Sample") %>%
  filter(!Species %in% c("Non Assigné", "Unassigned"))
n_larves_sp <- n_distinct(df_larves_norm_sp_quad$Sample)

df_core_scatter_sp <- calculer_core_scatter(
  df_larves_norm_sp_quad, "Species", n_larves_sp, SEUIL_PREVALENCE, SEUIL_ABONDANCE
)
p_scatter_sp <- creer_quadrant_plot(
  df_core_scatter_sp, "Species", SEUIL_PREVALENCE, SEUIL_ABONDANCE,
  titre = "Larval Microbiome Structure (16S - Species level)",
  style_label = "italic"
)
ggsave(file.path(out_dir, "Larval_Microbiome_Quadrants_Species.png"),
       plot = p_scatter_sp, width = 11, height = 8, dpi = 300)


# ==============================================================================
# BLOCK 11: INTRA-STAGE FIDELITY MATRICES
# ==============================================================================
# For the biologically "critical" stages (very young larvae, pupae,
# young adults), a heatmap is generated showing the abundance of each genus
# in each of the 4 replicates.
# Only genera present in >= 3 out of 4 replicates are kept (fidelity).
# This allows assessing biological reproducibility at each stage.
message("--- Matrices de fidélité intra-stade ---")

cond_het  <- c("Tiny_Larvae (2 mg)", "Pupae", "Young_Beetles")
plot_list <- list()

for (stade in cond_het) {

  reps_attendus <- meta %>% filter(Condition == stade) %>% pull(Sample) %>% as.character()

  df_stade <- df_final_16s %>%
    filter(Condition == stade) %>%
    normaliser_100("Genus") %>%
    left_join(meta %>% select(Sample, Condition), by = "Sample") %>%
    filter(!Genus %in% c("Non Assigné", "Unassigned")) %>%
    group_by(Genus) %>%
    # Fidelity criterion: present in at least 3 out of 4 replicates
    mutate(Nb_Present = sum(Relative_Abundance > 0)) %>%
    filter(Nb_Present >= 3) %>%
    ungroup()

  if (nrow(df_stade) > 0) {

    df_stade <- df_stade %>%
      mutate(Sample = factor(Sample, levels = reps_attendus)) %>%
      # complete() generates the missing cells (genus absent from a replicate) with 0
      # -> they are displayed as light grey on the heatmap
      complete(Sample, Genus, fill = list(Relative_Abundance = 0, Condition = stade)) %>%
      group_by(Genus) %>%
      mutate(Total_Abund = sum(Relative_Abundance)) %>%
      ungroup() %>%
      mutate(
        Genus      = fct_reorder(Genus, Total_Abund),
        # Text is only displayed if abundance >= 0.05% (avoids visual clutter)
        Text_Label = ifelse(Relative_Abundance >= 0.05,
                            sprintf("%.1f", Relative_Abundance), ""),
        # White text on dark background (>40%), black on light background
        Is_Dark    = Relative_Abundance > 40
      )

    p <- ggplot(df_stade, aes(x = Sample, y = Genus, fill = Relative_Abundance)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = Text_Label, color = Is_Dark),
                size = 3.2, fontface = "bold", show.legend = FALSE) +
      scale_color_manual(values = c("TRUE" = "white", "FALSE" = "grey20")) +
      # Non-linear gradient: 0=grey, 1%=light blue, >20%=dark blue
      # rescale() allows "zooming in" on the 0-20% range where most data lie
      scale_fill_gradientn(
        colors = c("#F2F4F4", "#AED6F1", "#2E86C1", "#154360"),
        values = rescale(c(0, 1, 20, 100)),
        limits = c(0, 100),
        name   = "Relative Abundance (%) :",
        breaks = c(0, 1, 10, 50, 100)
      ) +
      scale_x_discrete(drop = FALSE) +
      facet_grid(. ~ Condition) +
      labs(x = NULL, y = NULL) +
      theme_bw(base_size = 11) +
      theme(
        axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                        size = 9, face = "bold"),
        axis.text.y      = element_text(face = "italic", size = 10),
        strip.background = element_rect(fill = "#2C3E50"),
        strip.text       = element_text(color = "white", face = "bold", size = 12),
        panel.grid       = element_blank(),
        plot.margin      = margin(5, 5, 5, 5)
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
    ) &
    theme(legend.position = "bottom", legend.key.width = unit(2, "cm"))

  final_h <- max(7, sum(genus_counts) * 0.35 + 2.5)
  ggsave(filename = file.path(out_dir, "Genome_heterogeneous.png"),
         plot = p_combined, width = 9.5, height = final_h, dpi = 300, bg = "white")
  ggsave(filename = file.path(out_dir, "Genome_heterogeneous.svg"),
       plot = p_combined, 
       width = 18, 
       height = 8)

}


# ==============================================================================
# BLOCK 12: ADULT CORE MICROBIOME (BARPLOTS + DONUTS, GENUS + SPECIES)
# ==============================================================================
# For adults (Beetles), two types of visualization are generated:
#   - Stacked barplot (same as the larval core)
#   - Donut chart (pie chart with a central hole): particularly
#     suited for a single stage since it shows proportions in a circular form
message("--- Core Microbiote Adulte ---")

cond_adulte    <- c("Beetles")
df_adulte_brut <- df_final_16s %>% filter(Condition %in% cond_adulte)

# Internal function to create a donut chart
creer_donut <- function(df_plot, palette, titre) {
  df_donut <- df_plot %>%
    arrange(Taxon_Top) %>%
    mutate(
      ymax          = cumsum(Abundance),
      ymin          = lag(ymax, default = 0),
      labelPosition = (ymax + ymin) / 2,   # Center of each sector
      label         = paste0(Taxon_Top, "\n", round(Abundance, 1), "%")
    )

  ggplot(df_donut, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2.5, fill = Taxon_Top)) +
    geom_rect(color = "white", linewidth = 0.5) +
    # Labels are only shown for sectors > 2% (avoids overlaps)
    geom_text(data = filter(df_donut, Abundance > 2.0),
              aes(x = 4.8, y = labelPosition, label = label),
              size = 3.2, fontface = "bold", lineheight = 0.8) +
    # coord_polar turns the rectangular chart into a donut
    coord_polar(theta = "y") +
    xlim(c(1, 5.5)) +
    scale_fill_manual(values = palette) +
    labs(title = titre) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.text     = element_text(face = "italic", size = 10),
      legend.title    = element_text(face = "bold"),
      plot.title      = element_text(face = "bold", hjust = 0.5, size = 14,
                                     margin = margin(b = 5)),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin     = margin(10, 10, 10, 10)
    )
}

# --- Genus part ---
res_adulte_gen <- preparer_barplot_data(
  df_brut = df_adulte_brut, rang = "Genus",
  top_n = 12, seuil_min = 1.0, cond_levels = cond_adulte
)
p_core_adulte <- creer_barplot(
  res_adulte_gen$df_plot, res_adulte_gen$palette,
  titre = "Adult Microbiome (16S)", x_lab = NULL,
  fill_lab = "Bacterial Genus:", angle_x = 0, largeur_bar = 0.45
) +
  theme(axis.text.x = element_text(angle = 0, face = "bold", size = 12))

p_donut_adulte <- creer_donut(
  res_adulte_gen$df_plot, res_adulte_gen$palette,
  titre = "Adult Stage Microbiote Diversity  (16S)"
)
ggsave(file.path(out_dir, "Adult_Microbiome_Means.png"),
       plot = p_core_adulte, width = 6.5, height = 7, dpi = 300)
ggsave(file.path(out_dir, "Adult_Microbiome_Donut.png"),
       plot = p_donut_adulte, width = 9.5, height = 7, dpi = 300)

# --- Species part ---
df_adulte_sp <- forcer_nomenclature_binomiale(df_adulte_brut)
res_adulte_sp <- preparer_barplot_data(
  df_brut = df_adulte_sp, rang = "Species",
  top_n = 12, seuil_min = 1.0, cond_levels = cond_adulte
)
p_core_adulte_sp <- creer_barplot(
  res_adulte_sp$df_plot, res_adulte_sp$palette,
  titre = "Adult Microbiome (16S — Species level)", x_lab = NULL,
  fill_lab = "Bacterial Species:", angle_x = 0, largeur_bar = 0.45
) +
  theme(axis.text.x = element_text(angle = 0, face = "bold", size = 12))

p_donut_adulte_sp <- creer_donut(
  res_adulte_sp$df_plot, res_adulte_sp$palette,
  titre = "Adult Stage Microbiote Diversity (16S — Species level)"
)
ggsave(file.path(out_dir, "Adult_Microbiome_Means_Species.png"),
       plot = p_core_adulte_sp, width = 7.5, height = 7, dpi = 300)
ggsave(file.path(out_dir, "Adult_Microbiome_Donut_Species.png"),
       plot = p_donut_adulte_sp, width = 10.5, height = 7, dpi = 300)


# ==============================================================================
# BLOCK 13: QUADRANT CHARTS (PREVALENCE vs ABUNDANCE) — ADULTS
# ==============================================================================
message("--- Graphiques Quadrants Adultes ---")

# Thresholds for adults are slightly less strict than for larvae:
# there are fewer replicates, and adult biological diversity is narrower.
SEUIL_PREVALENCE_AD <- 75   # Present in >= 3 out of 4 replicates = 75%
SEUIL_ABONDANCE_AD  <- 1.0  # >= 1% mean relative abundance

# --- Genus ---
df_ad_norm_gen <- normaliser_100(df_adulte_brut, "Genus") %>%
  left_join(meta %>% select(Sample, Condition), by = "Sample") %>%
  filter(!Genus %in% c("Non Assigné", "Unassigned"))
n_ad <- n_distinct(df_ad_norm_gen$Sample)

df_core_ad_gen <- calculer_core_scatter(
  df_ad_norm_gen, "Genus", n_ad, SEUIL_PREVALENCE_AD, SEUIL_ABONDANCE_AD
)
p_scatter_ad <- creer_quadrant_plot(
  df_core_ad_gen, "Genus", SEUIL_PREVALENCE_AD, SEUIL_ABONDANCE_AD,
  titre = "Adult Microbiome Structure (Beetles - Genus level)"
)
ggsave(file.path(out_dir, "Adult_Microbiome_Quadrants.png"),
       plot = p_scatter_ad, width = 11, height = 8, dpi = 300)

ggsave(file.path(out_dir, "Adult_Microbiome_Quadrants.svg"),
       plot = p_scatter_ad, width = 11, height = 8)

# --- Species ---
df_adulte_sp_quad <- forcer_nomenclature_binomiale(df_adulte_brut)
df_ad_norm_sp <- normaliser_100(df_adulte_sp_quad, "Species") %>%
  left_join(meta %>% select(Sample, Condition), by = "Sample") %>%
  filter(!Species %in% c("Non Assigné", "Unassigned"))
n_ad_sp <- n_distinct(df_ad_norm_sp$Sample)

df_core_ad_sp <- calculer_core_scatter(
  df_ad_norm_sp, "Species", n_ad_sp, SEUIL_PREVALENCE_AD, SEUIL_ABONDANCE_AD
)
p_scatter_ad_sp <- creer_quadrant_plot(
  df_core_ad_sp, "Species", SEUIL_PREVALENCE_AD, SEUIL_ABONDANCE_AD,
  titre = "Adult Microbiome Structure (Beetles - Species level)",
  style_label = "italic"
)
ggsave(file.path(out_dir, "Adult_Microbiome_Quadrants_Species.png"),
       plot = p_scatter_ad_sp, width = 11, height = 8, dpi = 300)


# ==============================================================================
# BLOCK 14: PCA (WITHOUT SUBSTRATES)
# ==============================================================================
# PCA summarizes the structure of the microbial community in 2 dimensions.
# Samples close together on the PCA have similar microbial compositions.
# Separate groups of points indicate distinct communities.
#
# The HELLINGER TRANSFORMATION is essential before running a PCA on count
# data:
#   - Reduces the weight of dominant taxa (which would otherwise swamp the signal)
#   - Mitigates the "double zero problem" (two samples lacking a rare taxon
#     should not necessarily be considered similar because of that)
message("--- ACP Classique (Hellinger) - Insecte Uniquement ---")

# 1. Strict definition of the insect stages (without substrates)
chronologie_insecte <- c("Tiny_Larvae (2 mg)", "Microlarvae (7 mg)", 
                         "Larvae_S1 (14 mg)", "Larvae_S2 (40 mg)", 
                         "Larvae_S3 (65 mg)", "Larvae_S4 (100 mg)", 
                         "Pupae", "Young_Beetles", "Beetles")

# 2. Filtering to keep ONLY the insect AND remove "Unassigned" genera
df_pca_norm <- normaliser_100(
  df_final_16s %>% filter(Condition %in% chronologie_insecte, Genus != "Unassigned"), 
  "Genus"
)

mat_pca <- df_pca_norm %>%
  pivot_wider(id_cols = Sample, names_from = Genus,
              values_from = Relative_Abundance, values_fill = 0) %>%
  column_to_rownames("Sample")

# Hellinger transformation: takes the square root of relative abundances
mat_hellinger <- decostand(mat_pca, method = "hellinger")

# rda() without an environmental variable = standard PCA (covariance matrix)
pca_res <- rda(mat_hellinger)

# scaling = 1: preserves distances between samples (recommended in ecology)
pca_coords <- as.data.frame(scores(pca_res, display = "sites", scaling = 1))
colnames(pca_coords) <- c("PC1", "PC2")
pca_coords$Sample <- as.character(rownames(pca_coords))

# Variance explained by each axis (in %)
eig_vals <- pca_res$CA$eig
var_pc1  <- round(eig_vals[1] / sum(eig_vals) * 100, 1)
var_pc2  <- round(eig_vals[2] / sum(eig_vals) * 100, 1)

# Merging with metadata and forcing the factor order (chronologie_insecte)
df_plot_pca <- pca_coords %>%
  inner_join(meta %>% select(Sample, Condition) %>% mutate(Sample = as.character(Sample)),
             by = "Sample") %>%
  mutate(Condition = factor(Condition, levels = chronologie_insecte))

p_pca <- ggplot(df_plot_pca, aes(x = PC1, y = PC2, color = Condition, fill = Condition)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  # shape = 21: circle with a black outline (makes the points more visible)
  geom_point(size = 4, shape = 21, color = "black", stroke = 0.6) +
  scale_color_manual(values = my_cond_colors) +
  scale_fill_manual(values = my_cond_colors) +
  labs(
    title    = "PCA Hellinger-transformed relative abundances",
    x = sprintf("PC1 (%s%%)", var_pc1),
    y = sprintf("PC2 (%s%%)", var_pc2),
    fill  = "Development Stage:",
    color = "Development Stage:"
  ) +
  theme_bw(base_size = 18) + 
  theme(
    text             = element_text(size = 18), 
    legend.position  = "right",
    legend.title     = element_text(face = "bold"), 
    legend.text      = element_text(),              
    plot.title       = element_text(face = "bold"), 
    plot.subtitle    = element_text(color = "grey40"), 
    axis.title       = element_text(),              
    panel.grid.minor = element_blank()
  )
ggsave(file.path(out_dir, "PCA_microbiote_insecte_seul.png"),
       plot = p_pca, width = 10, height = 7, dpi = 300, bg = "white")
# ------------------------------------------------------------------------------
# CHART 2: PCA WITH 95% CONFIDENCE ELLIPSES
# ------------------------------------------------------------------------------
p_pca_ellipses <- ggplot(df_plot_pca, aes(x = PC1, y = PC2, color = Condition, fill = Condition)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  
  # 1. The ellipse fill (semi-transparent polygon)
  stat_ellipse(geom = "polygon", alpha = 0.15, level = 0.95, color = NA) +
  
  # 2. The ellipse outline (bolder line)
  stat_ellipse(geom = "path", linewidth = 0.8, level = 0.95) +
  
  # 3. The points displayed on top of the ellipses
  geom_point(size = 4, shape = 21, color = "black", stroke = 0.6) +
  
  scale_color_manual(values = my_cond_colors) +
  scale_fill_manual(values = my_cond_colors) +
  labs(
    title    = "PCA Hellinger-transformed relative abundances",
    x = sprintf("PC1 (%s%%)", var_pc1),
    y = sprintf("PC2 (%s%%)", var_pc2),
    fill  = "Development Stage:",
    color = "Development Stage:"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position  = "right",
    legend.title     = element_text(face = "bold"),
    plot.title       = element_text(face = "bold", size = 16),
    plot.subtitle    = element_text(color = "grey40", size = 12),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir, "PCA_microbiote_insecte_seul_ellipses.png"), plot = p_pca_ellipses, width = 10, height = 7, dpi = 300, bg = "white")

message("=> ACP avec et sans ellipses générées avec succès !")
# ==============================================================================
# BLOCK 15: ALPHA vs BETA CORRELATION (REPLICATES + CENTROIDS)
# ==============================================================================
# This chart tests and visualizes the relationship between:
#   - WITHIN-SAMPLE diversity (Shannon): on the Y axis
#   - HETEROGENEITY between replicates (mean Bray-Curtis): on the X axis
message("--- Corrélation Alpha vs Bêta ---")

# --- 1. Chart preparation ---
label_cor_ins <- sprintf("Spearman's rho = %s\np-value %s", rho_sans, p_val_sans)

# The X axis is capped at 1.0 (Bray-Curtis <= 1 by definition)
max_bray    <- max(df_ind_sans$Beta_Ind, na.rm = TRUE)
max_x_limit <- ifelse(max_bray > 0.9, 1.0, max_bray + 0.05)

p_corr_ins <- ggplot() +
  # Linear regression line (on individual replicates)
  geom_smooth(data = df_ind_sans, aes(x = Beta_Ind, y = Shannon),
              method = "lm", color = "black", linetype = "dashed",
              alpha = 0.15, linewidth = 0.8, fullrange = FALSE) +
              
  # Individual replicates (small points of FIXED size = 3)
  geom_point(data = df_ind_sans,
             aes(x = Beta_Ind, y = Shannon, fill = Condition),
             shape = 21, size = 3, color = "black", stroke = 0.4, alpha = 0.4,
             position = position_jitter(width = 0.005, height = 0.005)) +
             
  # PER-STAGE CENTROIDS (Fixed size, larger to stand out)
  geom_point(data = df_sans,
             aes(x = Beta_Mean, y = MeanShannon, fill = Condition),
             shape = 21, size = 7, color = "black", stroke = 1.2, alpha = 1) +
             
  # Color scales and legend
  scale_fill_manual(values = my_cond_colors, breaks = chronologie_insecte, name = "Development Stage:") +
  scale_color_manual(values = my_cond_colors, guide = "none") + # Disables the duplicate legend
  
  scale_x_continuous(limits = c(NA, max_x_limit), breaks = scales::pretty_breaks(n = 6)) +
  
  # Spearman test annotation
  annotate("text", x = -Inf, y = -Inf,
           label = label_cor_ins,
           hjust = -0.1, vjust = -0.5,
           fontface = "bold", color = "#C0392B", size = 5) +
           
  labs(
    title = "Relationship between inter-individual instability and microbial complexity",
    x = "Intra-group heterogeneity (Mean Bray-Curtis distance)",
    y = "Individual Alpha diversity (Shannon Index)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position  = "right",
    legend.title     = element_text(face = "bold"),
    plot.title       = element_text(face = "bold", size = 15),
    plot.subtitle    = element_text(color = "grey40", size = 12, margin = margin(b = 15)),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", linewidth = 1)
  ) +
  # Trick: forces the point size in the legend so they are clearly visible
  # and fully opaque (bypassing the alpha = 0.4 of the small points)
  guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1)))

ggsave(file.path(out_dir, "Correlation_Bray_vs_Shannon.png"),
       plot = p_corr_ins, width = 12, height = 7, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "Correlation_Bray_vs_Shannon.svg"),
       plot = p_corr_ins, width = 12, height = 7, bg = "white")

message("=> Graphique de corrélation Alpha/Bêta généré.")
# ==============================================================================
# BLOCK 16: INTRA-GENUS TAXONOMIC RESOLUTION (Fragmented Horizontal Barplots)
# ==============================================================================
# For the most abundant genera, this shows WHICH SPECIES compose them
# and in what proportion. Each genus is displayed in its own facet.
# The X axis = % of that genus's reads assigned to each species (sum = 100% per genus).
# It is a "non-stacked" chart (bars not overlapping): each species
# is an independent bar, which makes comparison between species easier.
message("--- Résolution Taxonomique Intra-Genre ---")

# ---> NEW: Added the 'exclure_genres' argument
generer_resolution_fragmente <- function(stades_cibles, titre_stade, nom_fichier,
                                          top_n_genres = 8, exclure_genres = NULL) {
  message(sprintf("  -> Stade(s) : %s", paste(stades_cibles, collapse = ", ")))

  df_stade <- df_final_16s %>% filter(Condition %in% stades_cibles)
  if (nrow(df_stade) == 0) return(NULL)

  # --- A. Exact mean abundance per genus (for the facet label) ---
  df_norm_gen  <- normaliser_100(df_stade, "Genus")
  n_samp_stade <- n_distinct(df_norm_gen$Sample)

  df_gen_agg <- df_norm_gen %>%
    group_by(Genus) %>%
    summarise(MeanAbund = sum(Relative_Abundance) / n_samp_stade, .groups = "drop") %>%
    # ---> NEW: The filter excludes unassigned entries AND the genera you chose
    filter(!Genus %in% c("Non Assigné", "Unassigned", exclure_genres))

  top_gen       <- df_gen_agg %>% slice_max(MeanAbund, n = top_n_genres) %>% pull(Genus)
  df_gen_labels <- df_gen_agg %>% filter(Genus %in% top_gen)

  # --- B. Intra-genus proportions by species ---
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
    # % computed over the total reads OF THAT GENUS (not the whole sample)
    mutate(Pct_in_Genus = (Count / sum(Count)) * 100) %>%
    ungroup()

  # --- C. Simplification (Top 5 species per genus + "Other") ---
  df_res_clean <- df_intra %>%
    group_by(Genus) %>%
    mutate(
      Rank = case_when(
        Species_Label == "Unclassified (sp.)" ~ 999,  # Forced to the bottom
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
    # The facet label includes the genus's overall abundance
    mutate(Genus_Label = sprintf("%s (%.2f%%)", Genus, MeanAbund))

  # --- D. Aesthetic ordering ---
  df_res_clean <- df_res_clean %>%
    mutate(
      Genus_Label = fct_reorder(Genus_Label, MeanAbund, .desc = TRUE),
      # fct_relevel forces the generic categories to the bottom of the chart
      Species_Final = fct_relevel(
        fct_reorder(Species_Final, Pct_in_Genus),
        "Unclassified (sp.)", "Other identified species",
        after = 0
      )
    )

  pal_gen <- setNames(colorspace::qualitative_hcl(length(top_gen), palette = "Dark 3"), top_gen)

# --- E. Fragmented chart (horizontal bars, 1 bar = 1 species) ---
  p_res <- ggplot(df_res_clean, aes(x = Pct_in_Genus, y = Species_Final)) +
    geom_col(aes(fill = Genus), color = "black", linewidth = 0.3, width = 0.7) +
    geom_text(aes(label = sprintf("%.2f%%", Pct_in_Genus)),
              hjust = -0.15, size = 3.5, fontface = "bold", color = "grey20") +
    facet_wrap(~ Genus_Label, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = pal_gen) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(
      title = sprintf("Intra-Genus Species Resolution : %s", titre_stade),
      x = "Relative Proportion within Genus (%)",
      y = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(
      legend.position  = "none",
      strip.background = element_rect(fill = "grey20"),
      strip.text       = element_text(color = "white", face = "bold.italic", size = 12),
      axis.text.y      = element_text(face = "italic", size = 10, color = "black"),
      axis.text.x      = element_text(size = 9),
      panel.grid.major.y = element_blank(),
      plot.margin      = margin(10, 15, 10, 10)
    )

  # ============================================================================
  # MULTI-FORMAT SAVING (PNG + SVG)
  # ============================================================================
  # 1. Saving in PNG format (Classic)
  ggsave(file.path(out_dir, nom_fichier), plot = p_res, width = 13, height = 9, dpi = 300, bg = "white")
  
  # 2. Replacing the .png extension with .svg via a regular expression
  nom_fichier_svg <- sub("\\.png$", ".svg", nom_fichier)
  
  # 3. Saving in SVG format (Vector)
  # Note: SVG format does not need the 'dpi' argument since it has infinite resolution
  ggsave(file.path(out_dir, nom_fichier_svg), plot = p_res, width = 13, height = 9, bg = "white")
  
  message(sprintf("  => '%s' (et format .svg) sauvegardés avec succès.", nom_fichier))
}

# ==============================================================================
# CHART EXECUTIONS
# ==============================================================================

# 1. For larvae
generer_resolution_fragmente(
  stades_cibles = cond_larves,
  titre_stade   = "Average Larval Microbiome",
  nom_fichier   = "Taxonomic_Resolution_Average_Larvae.png",
  top_n_genres  = 10,
  exclure_genres = c("Listeria","Citrobacter","Lactobacillus","Bacillus")
)

generer_resolution_fragmente(
  stades_cibles = c("Microlarvae (7 mg)"),
  titre_stade   = "Microlarvae Microbiome",
  nom_fichier   = "Taxonomic_Resolution_Microlarvae.png",
  top_n_genres  = 10
)

generer_resolution_fragmente(
  stades_cibles = c("Larvae_S1 (14 mg)"),
  titre_stade   = "Larvae S1 (14 mg) Microbiome",
  nom_fichier   = "Taxonomic_Resolution_Larvae_S1.png",
  top_n_genres  = 10
)

generer_resolution_fragmente(
  stades_cibles = c("Larvae_S2 (40 mg)"),
  titre_stade   = "Larvae S2 (40 mg) Microbiome",
  nom_fichier   = "Taxonomic_Resolution_Larvae_S2.png",
  top_n_genres  = 10
)

generer_resolution_fragmente(
  stades_cibles = c("Larvae_S3 (65 mg)"),
  titre_stade   = "Larvae_S3 (65 mg) Microbiome",
  nom_fichier   = "Taxonomic_Resolution_Larvae_S3.png",
  top_n_genres  = 10
)

generer_resolution_fragmente(
  stades_cibles = c("Larvae_S4 (100 mg)"),
  titre_stade   = "Larvae_S4 (100 mg) Microbiome",
  nom_fichier   = "Taxonomic_Resolution_Larvae_S4.png",
  top_n_genres  = 10
)

generer_resolution_fragmente(
  stades_cibles  = c("Beetles"),
  titre_stade    = "Adult Microbiome (Beetles)",
  nom_fichier    = "Taxonomic_Resolution_Adults.png",
  top_n_genres   = 10
)

# ==============================================================================
# BLOCK 17: DENDROGRAMS (ELBOW METHOD + CLUSTERS)
# ==============================================================================
message("--- Génération des Dendrogrammes (Elbow + Hellinger) ---")

if (!requireNamespace("ggdendro", quietly = TRUE)) install.packages("ggdendro")
suppressPackageStartupMessages(library(ggdendro))

# Main function
generer_dendrogramme_et_elbow <- function(df_counts, meta_df, prefixe_nom, titre_dendro, k_optimal) {
  
  # 1. Normalization and matrix preparation
  df_norm <- normaliser_100(df_counts %>% filter(Genus != "Unassigned"), "Genus")
  
  mat_wide <- df_norm %>%
    pivot_wider(id_cols = Sample, names_from = Genus, values_from = Relative_Abundance, values_fill = 0) %>%
    column_to_rownames("Sample")
  
  # 2. Hellinger transformation and matrix
  mat_hellinger <- decostand(mat_wide, method = "hellinger")
  dist_mat <- vegdist(mat_hellinger, method = "euclidean")
  hc <- hclust(dist_mat, method = "ward.D2")
  dendro <- dendro_data(hc, type = "rectangle")
  
  # --- PART A: ELBOW COMPUTATION AND CHART ---
  # Computing the WSS (Within-Cluster Sum of Squares) for 1 to 10 clusters
  k_max <- min(10, nrow(mat_hellinger) - 1)
  wss <- sapply(1:k_max, function(k) {
    clusters <- cutree(hc, k)
    sum(sapply(1:k, function(i) {
      mat_sub <- mat_hellinger[clusters == i, , drop = FALSE]
      center <- colMeans(mat_sub)
      sum(rowSums((sweep(mat_sub, 2, center))^2))
    }))
  })
  
  df_elbow <- data.frame(K = 1:k_max, WSS = wss)
  
  p_elbow <- ggplot(df_elbow, aes(x = K, y = WSS)) +
    geom_line(color = "steelblue", linewidth = 1) +
    geom_point(color = "red", size = 3) +
    # Dashed line showing the K value you chose as optimal
    geom_vline(xintercept = k_optimal, linetype = "dashed", color = "grey40") +
    scale_x_continuous(breaks = 1:k_max) +
    labs(
      title = paste("Elbow Method :", prefixe_nom),
      x = "Number of clusters (k)", 
      y = "Within-Cluster Sum of Squares (WSS)"
    ) +
    theme_bw(base_size = 13)
    
  ggsave(file.path(out_dir, paste0(prefixe_nom, "_Elbow_Plot.png")), plot = p_elbow, width = 8, height = 6, dpi = 300)
  
  # --- PART B: DENDROGRAM WITH CLUSTERS ---
  labels_df <- label(dendro) %>%
    rename(Sample = label) %>%
    left_join(meta_df, by = "Sample")
    
  # The tree is cut at k_optimal
  groupes <- cutree(hc, k = k_optimal)
  labels_df$Cluster <- as.factor(groupes[labels_df$Sample])
  
  # Exporting the cluster composition to a table
  nom_tsv <- paste0(prefixe_nom, "_Clusters_Composition.tsv")
  labels_df %>%
    select(Cluster, Sample, Condition) %>%
    arrange(Cluster, Condition) %>%
    write_tsv(file.path(out_dir, nom_tsv))
  
  # Computing the boxes (rectangles) framing the clusters
  rect_df <- labels_df %>%
    group_by(Cluster) %>%
    summarise(xmin = min(x) - 0.45, xmax = max(x) + 0.45, .groups = "drop") %>%
    mutate(ymin = -0.015, ymax = max(segment(dendro)$y) * 0.95)
  
  # Creating the Dendrogram chart
  p_dendro <- ggplot() +
    # Grey cluster boxes
    geom_rect(data = rect_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "grey50", alpha = 0.15, color = "black", linetype = "dashed", inherit.aes = FALSE) +
    
    # Tree branches
    geom_segment(data = segment(dendro), aes(x = x, y = y, xend = xend, yend = yend), 
                 color = "grey60", linewidth = 0.6) +
                 
    # Terminal points (colored by larval stage)
    geom_point(data = labels_df, aes(x = x, y = y, fill = Condition), 
               shape = 21, size = 3, color = "black", stroke = 0.5) +
               
    # Sample text labels
    geom_text(data = labels_df, aes(x = x, y = y - 0.03, label = Sample, color = Condition),
              angle = 90, hjust = 1, vjust = 0.5, size = 4, fontface = "bold", show.legend = FALSE) + 
              
    scale_color_manual(values = my_cond_colors, breaks = names(my_cond_colors)) +
    scale_fill_manual(values = my_cond_colors, breaks = names(my_cond_colors)) +
    scale_y_continuous(expand = expansion(mult = c(0.25, 0.08))) +
    labs(
      title = titre_dendro,
      x = NULL, y = "Height (Euclidean Distance)", color = "Development Stage:", fill = "Development Stage:"
    ) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      legend.position = "bottom", plot.title = element_text(face = "bold", size = 15),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  ggsave(file.path(out_dir, paste0(prefixe_nom, "_Dendrogram.png")), plot = p_dendro, width = 14, height = 8, dpi = 300, bg = "white")
  message(sprintf("  => Elbow & Dendrogramme '%s' sauvegardés.", prefixe_nom))
}

# --- EXECUTION ---

generer_dendrogramme_et_elbow(
  df_counts = df_final_16s %>% filter(Condition %in% chronologie_tout),
  meta_df   = meta,
  prefixe_nom = "Global",
  titre_dendro = "Global Hierarchical Clustering (All Samples)",
  k_optimal = 4 
)

generer_dendrogramme_et_elbow(
  df_counts = df_final_16s %>% filter(Condition %in% chronologie_insecte),
  meta_df   = meta,
  prefixe_nom = "Insect_Only",
  titre_dendro = "Insect Hierarchical Clustering (Substrates Excluded)",
  k_optimal = 3 
)



# ==============================================================================
# BLOCK 18: COMMUNITY STRUCTURE VIA BRAY-CURTIS (INSECT ONLY)
# ==============================================================================
# This block explores community structuring using the Bray-Curtis distance
# (taxonomic dissimilarity), via two complementary algorithms:
#   1. PCoA (Metric MDS): Attempts to project the exact distance in 2D.
#   2. NMDS (Non-Metric MDS): Preserves the rank order of distances,
#      very robust for microbial ecology.
# Substrates are excluded to focus on host dynamics.
# ==============================================================================
message("--- Generating PCoA and NMDS (Bray-Curtis / Insect Only) ---")

# 1. Strict filtering (Insect only) and normalization
# 'chronologie_insecte' is reused since it naturally excludes substrates
df_bray_norm <- normaliser_100(
  df_final_16s %>% filter(Condition %in% chronologie_insecte, Genus != "Unassigned"),
  "Genus"
)

# 2. Switching to a wide matrix for vegan
mat_bray_multiv <- df_bray_norm %>%
  pivot_wider(id_cols = Sample, names_from = Genus, values_from = Relative_Abundance, values_fill = 0) %>%
  column_to_rownames("Sample")

# 3. Computing the actual Bray-Curtis distance matrix
dist_bray_multiv <- vegdist(mat_bray_multiv, method = "bray")


# ------------------------------------------------------------------------------
# PART A: PCoA (METRIC MDS / "The Bray-Curtis PCA")
# ------------------------------------------------------------------------------
# The cmdscale() function performs the PCoA. eig = TRUE lets us recover the variance.
pcoa_res <- cmdscale(dist_bray_multiv, k = 2, eig = TRUE)

# Extracting the coordinates
pcoa_coords <- as.data.frame(pcoa_res$points)
colnames(pcoa_coords) <- c("PCoA1", "PCoA2")
pcoa_coords$Sample <- rownames(pcoa_coords)

# Computing the percentage of variance explained by the axes
pcoa_eig <- pcoa_res$eig
var_pcoa1 <- round(pcoa_eig[1] / sum(pcoa_eig[pcoa_eig > 0]) * 100, 1)
var_pcoa2 <- round(pcoa_eig[2] / sum(pcoa_eig[pcoa_eig > 0]) * 100, 1)

# Joining with metadata for ggplot
df_plot_pcoa <- pcoa_coords %>%
  inner_join(meta %>% select(Sample, Condition), by = "Sample") %>%
  mutate(Condition = factor(Condition, levels = chronologie_insecte))

# Creating the PCoA chart
p_pcoa <- ggplot(df_plot_pcoa, aes(x = PCoA1, y = PCoA2, color = Condition, fill = Condition)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  
  geom_point(size = 4, shape = 21, color = "black", stroke = 0.6) +
  
  scale_color_manual(values = my_cond_colors) +
  scale_fill_manual(values = my_cond_colors) +
  
  labs(
    title = "Principal Coordinates Analysis (PCoA)",
    subtitle = "Based on Bray-Curtis dissimilarity (Insect only)",
    x = sprintf("PCoA1 (%s%%)", var_pcoa1),
    y = sprintf("PCoA2 (%s%%)", var_pcoa2),
    fill = "Development Stage:",
    color = "Development Stage:"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir, "PCoA_BrayCurtis_Insecte.png"), plot = p_pcoa, width = 10, height = 7, dpi = 300, bg = "white")


# ------------------------------------------------------------------------------
# PART B: NMDS (NON-METRIC MDS)
# ------------------------------------------------------------------------------
# NMDS is an iterative process. A 'seed' is fixed so that the chart
# is 100% identical if you rerun the script tomorrow.
set.seed(42) 

# trymax = 100 forces the algorithm to actively search for the best convergence
nmds_res <- metaMDS(mat_bray_multiv, distance = "bray", k = 2, trymax = 100, trace = FALSE)

# Extracting the coordinates
nmds_coords <- as.data.frame(scores(nmds_res, display = "sites"))
nmds_coords$Sample <- rownames(nmds_coords)

df_plot_nmds <- nmds_coords %>%
  inner_join(meta %>% select(Sample, Condition), by = "Sample") %>%
  mutate(Condition = factor(Condition, levels = chronologie_insecte))

# 'Stress' is THE vital NMDS metric.
# < 0.2 = Good | < 0.1 = Excellent | < 0.05 = Perfect
stress_val <- round(nmds_res$stress, 3)

# Creating the NMDS chart
p_nmds <- ggplot(df_plot_nmds, aes(x = NMDS1, y = NMDS2, color = Condition, fill = Condition)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  
  geom_point(size = 4, shape = 21, color = "black", stroke = 0.6) +
  
  scale_color_manual(values = my_cond_colors) +
  scale_fill_manual(values = my_cond_colors) +
  
  labs(
    title = "Non-Metric Multidimensional Scaling (NMDS)",
    subtitle = sprintf("Bray-Curtis dissimilarity (Insect only) | Stress: %s", stress_val),
    # Note: Unlike PCA/PCoA, NMDS axes have no "percentage" of variance
    x = "NMDS Axis 1",
    y = "NMDS Axis 2",
    fill = "Development Stage:",
    color = "Development Stage:"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir, "NMDS_BrayCurtis_Insecte.png"), plot = p_nmds, width = 10, height = 7, dpi = 300, bg = "white")
message("=> PCoA et NMDS (Insecte seul) générés avec succès.")

# ==============================================================================
# BLOCK 20: GLOBAL LARVAL CORE MICROBIOME (BARPLOT + DONUT, GENUS + SPECIES)
# ==============================================================================
# Goal: Get an "overall mean" view of the larval phase (Microlarvae to S4)
# so it can be visually compared to the Adult stage "identity card".
# Method: Strict manual computation. Relative abundances are summed across the
# 20 larval samples, then divided by 20.
# ==============================================================================
message("--- Generating Global Larval Microbiome (Barplot & Donut) ---")

# 1. Strict selection of the 20 larval samples
cond_larves_cibles <- c(
  "Microlarvae (7 mg)", "Larvae_S1 (14 mg)", 
  "Larvae_S2 (40 mg)", "Larvae_S3 (65 mg)", "Larvae_S4 (100 mg)"
)

df_larves_global_brut <- df_final_16s %>% 
  filter(Condition %in% cond_larves_cibles)

# ------------------------------------------------------------------------------
# PART A: GENUS LEVEL
# ------------------------------------------------------------------------------
# 1. Normalization per sample
df_norm_gen_glob <- normaliser_100(df_larves_global_brut, "Genus")

# 2. Manual computation of the overall mean (Sum / 20 replicates)
n_replicats <- n_distinct(df_norm_gen_glob$Sample) # Automatically finds 20

df_glob_moy_gen <- df_norm_gen_glob %>%
  group_by(Genus) %>%
  summarise(Relative_Abundance = sum(Relative_Abundance) / n_replicats, .groups = "drop") %>%
  mutate(Condition = "Global Larval Stage") # Creating the single condition
  
# 3. Top 12 and grouping into "Others"
top_12_gen <- df_glob_moy_gen %>%
  filter(!Genus %in% c("Non Assigné", "Unassigned")) %>%
  slice_max(Relative_Abundance, n = 12) %>%
  pull(Genus)
  
df_plot_glob_gen <- df_glob_moy_gen %>%
  mutate(Taxon_Top = case_when(
    Genus %in% c("Non Assigné", "Unassigned") ~ "Unassigned",
    Genus %in% top_12_gen ~ Genus,
    TRUE ~ "Others"
  )) %>%
  group_by(Condition, Taxon_Top) %>%
  summarise(Abundance = sum(Relative_Abundance), .groups = "drop") %>%
  mutate(
    Taxon_Top = fct_relevel(fct_reorder(factor(Taxon_Top), Abundance, sum, .desc = FALSE), "Others", "Unassigned", after = 0)
  )
  
# 4. Manual palette creation (to match the Top 12 exactly)
taxa_vrais_glob_gen <- setdiff(levels(df_plot_glob_gen$Taxon_Top), c("Others", "Unassigned"))
pal_glob_gen <- setNames(rep_len(pal_taxo, length(taxa_vrais_glob_gen)), taxa_vrais_glob_gen)
if ("Others" %in% levels(df_plot_glob_gen$Taxon_Top)) pal_glob_gen["Others"] <- "grey85"
if ("Unassigned" %in% levels(df_plot_glob_gen$Taxon_Top)) pal_glob_gen["Unassigned"] <- "grey40"

# 5. Charts
p_core_larves_glob <- creer_barplot(
  df_plot_glob_gen, pal_glob_gen,
  titre = "Global Larval Microbiome (16S)", x_lab = NULL,
  fill_lab = "Bacterial Genus:", angle_x = 0, largeur_bar = 0.45
) + theme(axis.text.x = element_text(angle = 0, face = "bold", size = 12))

p_donut_larves_glob <- creer_donut(
  df_plot_glob_gen, pal_glob_gen,
  titre = "Global Larval Stage Identity Card (16S)"
)

ggsave(file.path(out_dir, "Global_Larval_Microbiome_Means.png"), plot = p_core_larves_glob, width = 6.5, height = 7, dpi = 300)
ggsave(file.path(out_dir, "Global_Larval_Microbiome_Donut.png"), plot = p_donut_larves_glob, width = 9.5, height = 7, dpi = 300)


ggsave(filename = file.path(out_dir, "Global_Larval_Microbiome_Means.svg"),
       plot = p_core_larves_glob, 
       width = 12, 
       height = 14)

ggsave(filename = file.path(out_dir, "Global_Larval_Microbiome_Donuts.svg"),
       plot = p_donut_larves_glob, 
       width = 12, 
       height = 14)
# ------------------------------------------------------------------------------
# PART B: SPECIES LEVEL
# ------------------------------------------------------------------------------
# 1. Binomial nomenclature and normalization
df_larves_global_brut_sp <- forcer_nomenclature_binomiale(df_larves_global_brut)
df_norm_sp_glob <- normaliser_100(df_larves_global_brut_sp, "Species")

# 2. Manual computation of the overall mean
df_glob_moy_sp <- df_norm_sp_glob %>%
  group_by(Species) %>%
  summarise(Relative_Abundance = sum(Relative_Abundance) / n_replicats, .groups = "drop") %>%
  mutate(Condition = "Global Larval Stage")
  
# 3. Top 12
top_12_sp <- df_glob_moy_sp %>%
  filter(!Species %in% c("Non Assigné", "Unassigned")) %>%
  slice_max(Relative_Abundance, n = 12) %>%
  pull(Species)
  
df_plot_glob_sp <- df_glob_moy_sp %>%
  mutate(Taxon_Top = case_when(
    Species %in% c("Non Assigné", "Unassigned") ~ "Unassigned",
    Species %in% top_12_sp ~ Species,
    TRUE ~ "Others"
  )) %>%
  group_by(Condition, Taxon_Top) %>%
  summarise(Abundance = sum(Relative_Abundance), .groups = "drop") %>%
  mutate(
    Taxon_Top = fct_relevel(fct_reorder(factor(Taxon_Top), Abundance, sum, .desc = FALSE), "Others", "Unassigned", after = 0)
  )


# ==============================================================================
# BLOCK 25: GENUS COMPARISON (LARVAE VS ADULTS) - CUSTOM DUMBBELL PLOT
# ==============================================================================
message("--- Génération du Dumbbell Plot (Échelle non-linéaire) ---")

# 1. Selecting the target stages
cond_larves_cibles <- c("Microlarvae (7 mg)", "Larvae_S1 (14 mg)", "Larvae_S2 (40 mg)", "Larvae_S3 (65 mg)", "Larvae_S4 (100 mg)")
cond_adulte        <- c("Beetles")

# 2. List of manually chosen genera
genres_cibles <- c(
  "Staphylococcus", "Enterococcus", "Lactococcus", "Enterobacter", 
  "Spiroplasma", "Mammaliicoccus", "Weissella", "Brevibacterium", 
  "Corynebacterium", "Kocuria", "Gordonia", "Serratia", 
  "Citrobacter", "Pectobacterium"
)

# 3. Using the global table so the 100% baseline is not distorted
df_sp_raw <- df_final_16s %>%
  filter(Condition %in% c(cond_larves_cibles, cond_adulte))

# 4. Exact normalization (100% computed BEFORE any genus filtering)
df_norm_sp <- normaliser_100(df_sp_raw, "Genus") %>%
  left_join(meta %>% select(Sample, Condition), by = "Sample") %>%
  mutate(Stage = ifelse(Condition %in% cond_larves_cibles, "Larval Stage", "Adult Stage"))

# 5. Computing the weighted OVERALL MEAN
n_larves  <- n_distinct(df_norm_sp %>% filter(Stage == "Larval Stage") %>% pull(Sample))
n_adultes <- n_distinct(df_norm_sp %>% filter(Stage == "Adult Stage") %>% pull(Sample))

df_gen_moy <- df_norm_sp %>%
  # FILTERING HERE ON YOUR CUSTOM LIST
  filter(Genus %in% genres_cibles) %>%
  group_by(Stage, Genus) %>%
  summarise(
    MeanAbund = ifelse(first(Stage) == "Larval Stage", 
                       sum(Relative_Abundance) / n_larves, 
                       sum(Relative_Abundance) / n_adultes),
    .groups = "drop"
  )

# 6. Wide formatting for the Dumbbell Plot
df_dumbbell <- df_gen_moy %>%
  pivot_wider(names_from = Stage, values_from = MeanAbund, values_fill = 0) %>%
  # The difference (Shift) is computed to order the Y axis
  mutate(Shift = `Adult Stage` - `Larval Stage`) %>%
  arrange(Shift) %>%
  mutate(Genus = factor(Genus, levels = Genus))

# 7. Creating the Dumbbell chart (with sqrt transformation)
p_dumbbell <- ggplot(df_dumbbell) +
  geom_segment(aes(y = Genus, yend = Genus, 
                   x = `Larval Stage`, xend = `Adult Stage`), 
               color = "grey60", linewidth = 1.5) +
  geom_point(aes(y = Genus, x = `Larval Stage`, color = "Larval Stage"), size = 5) +
  geom_point(aes(y = Genus, x = `Adult Stage`, color = "Adult Stage"), size = 5) +
  
  scale_color_manual(values = c("Larval Stage" = "#2980B9", "Adult Stage" = "#E74C3C"), 
                     name = "Development Stage:") +
  
  # ---> CHANGE HERE: X axis transformation <---
  scale_x_continuous(
    trans = "sqrt",  # Mathematically transforms the scale (stretches small values)
    breaks = c(0, 0.5, 1, 2.5, 5, 10, 20, 30), # Forces the tick marks we want to see
    labels = function(x) paste0(x, "%")
  ) +
  
  labs(
    title = "Taxonomic comparison of the larval and adult stages",
    x = "Mean Relative Abundance (%)",
    y = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.x = element_line(color = "grey85", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(face = "italic", size = 12, color = "black"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16)
  )

ggsave(file.path(out_dir, "Comparison_Selected_Genus_Dumbbell.png"), plot = p_dumbbell, width = 10, height = 8, dpi = 300)
message("=> Dumbbell Plot (Échelle Racine Carrée) généré avec succès.")


# ==============================================================================
# BLOCK 26: MIRROR PLOTS (WITH CONFIDENCE SPREAD / STANDARD ERROR) - FIXED
# ==============================================================================
message("--- Génération des Mirror Plots (Classique et Colonne Vertébrale + SE) ---")

# 1. Data selection
cond_larves_cibles <- c("Microlarvae (7 mg)", "Larvae_S1 (14 mg)", "Larvae_S2 (40 mg)", "Larvae_S3 (65 mg)", "Larvae_S4 (100 mg)")
cond_adulte        <- c("Beetles")

genres_cibles <- c(
  "Staphylococcus", "Enterococcus", "Lactococcus", "Enterobacter", 
  "Spiroplasma", "Mammaliicoccus", "Weissella", "Brevibacterium", 
  "Corynebacterium", "Serratia", 
  "Citrobacter", "Brachybacterium", "Lactobacillus", "Pediococcus"
)

# 2. Retrieval and normalization
df_sp_raw <- df_final_16s %>% filter(Condition %in% c(cond_larves_cibles, cond_adulte))

df_norm_sp <- normaliser_100(df_sp_raw, "Genus") %>%
  left_join(meta %>% select(Sample, Condition), by = "Sample") %>%
mutate(Stage = ifelse(Condition %in% cond_larves_cibles, "Larval Stage", "Beetles Stage"),
       Stage = factor(Stage, levels = c("Larval Stage", "Beetles Stage")))

# 3. Filling in the zeros (essential for the statistical test and error to be valid)
samples_inclus <- df_norm_sp %>% distinct(Sample, Stage)

df_complet <- df_norm_sp %>%
  filter(Genus %in% genres_cibles) %>%
  select(Sample, Genus, Relative_Abundance) %>%
  complete(Sample, Genus = genres_cibles, fill = list(Relative_Abundance = 0)) %>%
  left_join(samples_inclus, by = "Sample")

# 4. COMPUTING MEANS AND STATISTICAL TEST (Wilcoxon + FDR)
df_stats <- df_complet %>%
  group_by(Genus) %>%
  summarise(
    Mean_Larvae = mean(Relative_Abundance[Stage == "Larval Stage"]),
    Mean_Adults = mean(Relative_Abundance[Stage == "Beetles Stage"]),
    p_val_brute = tryCatch(wilcox.test(Relative_Abundance ~ Stage, exact = FALSE)$p.value, 
                           error = function(e) 1),
    .groups = "drop"
  ) %>%
  mutate(
    Shift = Mean_Adults - Mean_Larvae,
    p_val_adj = p.adjust(p_val_brute, method = "BH"),
    Sig = case_when(
      p_val_adj < 0.01  ~ "**",
      p_val_adj < 0.05  ~ "*",
      TRUE ~ "" 
    )
  )

# 5. Formatting the data for the chart (Including Standard Error computation)
df_plot <- df_complet %>%
  group_by(Stage, Genus) %>%
  summarise(
    MeanAbund = mean(Relative_Abundance),
    SEAbund   = sd(Relative_Abundance) / sqrt(n()), # Computing the Standard Error
    .groups = "drop"
  ) %>%
  left_join(df_stats %>% select(Genus, Shift, Sig), by = "Genus") %>%
  mutate(
    Plot_Mean = ifelse(Stage == "Larval Stage", -MeanAbund, MeanAbund),
    Genus = fct_reorder(Genus, Shift),
    Label_Text = ifelse(MeanAbund > 0 & MeanAbund < 0.05, "<0.1%", sprintf("%.1f%%", MeanAbund)),
    
    Err_Min_Classic = ifelse(Stage == "Larval Stage", -MeanAbund - SEAbund, pmax(0, MeanAbund - SEAbund)),
    Err_Max_Classic = ifelse(Stage == "Larval Stage", -pmax(0, MeanAbund - SEAbund), MeanAbund + SEAbund),
    Label_X_Classic = ifelse(Stage == "Larval Stage", -MeanAbund - SEAbund - 1.2, MeanAbund + SEAbund + 1.2)
  )

# Global settings
col_larvae <- "#3678a3fd" 
col_adults <- "#e26759ff" 
col_text   <- "#2C3E50" 

# ==============================================================================
# FIGURE A: THE CLASSIC PREMIUM MIRROR PLOT
# ==============================================================================
breaks_x_std <- seq(-30, 30, by = 10)
labels_x_std <- paste0(abs(breaks_x_std), "%")

p_mirror_classic <- ggplot(df_plot, aes(y = Genus)) +
  geom_vline(xintercept = 0, color = col_text, linewidth = 1) +
  geom_col(aes(x = Plot_Mean, fill = Stage), width = 0.75, color = "white", linewidth = 0.5) +
  
  geom_errorbarh(aes(xmin = Err_Min_Classic, xmax = Err_Max_Classic), 
                 height = 0.25, color = "#7F8C8D", linewidth = 0.6) +
  
  geom_text(aes(x = Label_X_Classic, label = Label_Text, color = Stage,
                hjust = ifelse(Stage == "Larval Stage", 1, 0)), 
            size = 3.8, fontface = "bold") +
  
  geom_label(data = df_stats %>% filter(Sig != ""), 
             aes(x = 0, y = Genus, label = Sig),
             fill = alpha("white", 0.85), color = col_text, label.size = 0, 
             size = 6, fontface = "bold", label.padding = unit(0.15, "lines")) +
  
  scale_fill_manual(values = c("Larval Stage" = col_larvae, "Beetles Stage" = col_adults), name = "") +
  scale_color_manual(values = c("Larval Stage" = col_larvae, "Beetles Stage" = col_adults), guide = "none") +
  scale_x_continuous(breaks = breaks_x_std, labels = labels_x_std, limits = c(-38, 38)) + 
  
  labs(title = "Taxonomic comparison of the larval and adult stages",
       subtitle = "Mean relative abundance ± SE & FDR-adjusted Wilcoxon Test (*p<0.05, **p<0.01)",
       x = "Relative Abundance (%)", y = NULL) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "top", legend.text = element_text(size = 14, face = "bold", color = col_text),
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5, color = col_text, margin = margin(b = 8)),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "#7F8C8D", margin = margin(b = 20)),
    panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey85", linetype = "dotted", linewidth = 0.8),
    axis.text.y = element_text(face = "bold.italic", size = 14, color = col_text),
    axis.text.x = element_text(face = "bold", size = 12, color = "#7F8C8D"),
    axis.title.x = element_text(face = "bold", size = 14, color = col_text, margin = margin(t = 15)),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  )

ggsave(file.path(out_dir, "Comparison_Mirror_Plot_Classic.png"), plot = p_mirror_classic, width = 12, height = 9, dpi = 600, bg = "white")

# ==============================================================================
# FIGURE B: THE "CENTRAL SPINE" MIRROR PLOT
# ==============================================================================
gap <- 15 
val_breaks <- c(10, 20, 30)
breaks_x_spine <- c(-val_breaks - gap, val_breaks + gap)
labels_x_spine <- c(paste0(val_breaks, "%"), paste0(val_breaks, "%"))

df_plot_spine <- df_plot %>%
  mutate(
    Genus_num = as.numeric(Genus),
    X_min = ifelse(Stage == "Larval Stage", -MeanAbund - gap, gap),
    X_max = ifelse(Stage == "Larval Stage", -gap, MeanAbund + gap),
    
    # THE FIX IS HERE: pmin() instead of min()
    Err_Min_Spine = ifelse(Stage == "Larval Stage", X_min - SEAbund, pmax(gap, X_max - SEAbund)),
    Err_Max_Spine = ifelse(Stage == "Larval Stage", pmin(-gap, X_min + SEAbund), X_max + SEAbund),
    
    Label_Pct_X = ifelse(Stage == "Larval Stage", X_min - SEAbund - 1.2, X_max + SEAbund + 1.2),
    Label_Pct_Hjust = ifelse(Stage == "Larval Stage", 1, 0)
  )

df_axis <- df_stats %>%
  mutate(Genus = factor(Genus, levels = levels(df_plot$Genus)), Genus_num = as.numeric(Genus),
         Label_Center = paste0(Genus, " ", Sig))

p_mirror_spine <- ggplot() +
  geom_vline(xintercept = -gap, color = "#BDC3C7", linewidth = 1) +
  geom_vline(xintercept = gap,  color = "#BDC3C7", linewidth = 1) +
  geom_rect(data = df_plot_spine, aes(xmin = X_min, xmax = X_max, ymin = Genus_num - 0.35, ymax = Genus_num + 0.35, fill = Stage), 
            color = "white", linewidth = 0.5) +
  
  geom_errorbarh(data = df_plot_spine, aes(y = Genus_num, xmin = Err_Min_Spine, xmax = Err_Max_Spine), 
                 height = 0.25, color = "#7F8C8D", linewidth = 0.6) +
                 
# Central labels (genus names) — was size = 12
geom_text(data = df_axis, aes(x = 0, y = Genus_num, label = Label_Center), 
          fontface = "bold.italic", size = 9.84, color = col_text) +

# Percentages — was size = 10
geom_text(data = df_plot_spine, aes(x = Label_Pct_X, y = Genus_num, label = Label_Text, 
          color = Stage, hjust = Label_Pct_Hjust), 
          size = 9.84, fontface = "bold") +
  
  scale_fill_manual(values = c("Larval Stage" = col_larvae, "Beetles Stage" = col_adults), name = "") +
  scale_color_manual(values = c("Larval Stage" = col_larvae, "Beetles Stage" = col_adults), guide = "none") +
  scale_x_continuous(breaks = breaks_x_spine, labels = labels_x_spine, limits = c(-55, 55)) + 
  scale_y_continuous(breaks = NULL, name = NULL) + 
  
  labs(title = "Taxonomic comparison of the larval and adult stages",
       subtitle = "Mean relative abundance ± SE & FDR-adjusted Wilcoxon Test (*p<0.05, **p<0.01)",
       x = "Relative Abundance (%)") +
       
theme_minimal(base_size = 28) +
theme(
    legend.position   = "top", 
    legend.text       = element_text(size = 28, face = "bold", color = col_text),
    legend.key.height = unit(1.5, "cm"),
    legend.key.width  = unit(2.5, "cm"),
    legend.spacing.x  = unit(0.5, "cm"),
    legend.margin     = margin(b = 15),

    plot.title    = element_text(face = "bold", size = 28, hjust = 0.5, color = col_text, margin = margin(b = 8)),
    plot.subtitle = element_text(size = 28, hjust = 0.5, color = "#7F8C8D", margin = margin(b = 20)),
    
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey85", linetype = "dotted", linewidth = 0.8),
    
    axis.text.y  = element_text(face = "bold.italic", size = 28, color = col_text),
    axis.text.x  = element_text(face = "bold", size = 28, color = "#7F8C8D"),
    axis.title.x = element_text(face = "bold", size = 28, color = col_text, margin = margin(t = 15)),
    plot.margin  = margin(t = 20, r = 20, b = 20, l = 20)
)

ggsave(file.path(out_dir, "Comparison_Mirror_Plot_Spine.png"), plot = p_mirror_spine, width = 21, height = 15, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "Comparison_Mirror_Plot_Spine.svg"), 
       plot = p_mirror_spine, width = 24, height = 18, bg = "white")
message("=> Mirror Plots (SE + Wilcoxon-FDR) générés avec succès.")
