################################################################################
# SCRIPT : ANALYSE DE LA DIVERSITÉ MICROBIENNE 16S — Tenebrio molitor
################################################################################
#
# Ce script analyse les communautés bactériennes intestinales de T. molitor
# à travers les différents stades de développement (larves → adultes + substrat)
#
# GRAPHIQUES PRODUITS :
#
#  1. Valeurs_Shannon_BrayCurtis.tsv    → Tableau numérique brut (annexe)
#  2. shannon_vs_braycurtis_.png        → Courbes alpha/bêta (insecte seul)
#  3. shannon_vs_braycurtis_substrats.png → Courbes alpha/bêta (+ environnement)
#  4. Barplot_Global_16S.png            → Composition par réplicat (tous stades)
#  5. Matrice_BrayCurtis_Tous_Replicats.png → Heatmap similarité (tous)
#  6. Matrice_BrayCurtis_Insecte_Seul.png   → Heatmap similarité (insecte)
#  7. Master_Figure.png / _substrats.png    → Figure maîtresse combinée
#  8. Core_Larval_Microbiome_Means.png      → Core larvaire (Genre)
#  9. Core_Larval_Microbiome_Means_Species.png → Core larvaire (Espèce)
# 10. Core_Microbiome_Quadrants.png         → Quadrants Prévalence/Abondance (larves, Genre)
# 11. Core_Microbiome_Quadrants_Species.png → Quadrants (larves, Espèce)
# 12. Core_genome_heterogeneous.png         → Matrices de fidélité intra-stade
# 13. Core_Adult_Microbiome_Means.png       → Core adulte (Genre)
# 14. Core_Adult_Microbiome_Donut.png       → Donut adulte (Genre)
# 15. Core_Adult_Microbiome_Means_Species.png → Core adulte (Espèce)
# 16. Core_Adult_Microbiome_Donut_Species.png → Donut adulte (Espèce)
# 17. Core_Adult_Microbiome_Quadrants.png   → Quadrants adultes (Genre)
# 18. Core_Adult_Microbiome_Quadrants_Species.png → Quadrants adultes (Espèce)
# 19. PCA_microbiote.png                    → ACP (structure des communautés)
# 20. Correlation_Bray_vs_Shannon.png/.svg  → Corrélation alpha vs bêta
# 21. Taxonomic_Resolution_Average_Larvae.png → Résolution intra-genre (larves)
# 22. Taxonomic_Resolution_Adults.png         → Résolution intra-genre (adultes)
#
# DONNÉES REQUISES :
#   - Table QIIME2-like (séquences × échantillons, format TSV)
#   - Fichier de métadonnées (échantillon → condition de développement)
#   - Bases SILVA v138.2 pour l'assignation taxonomique préalable
################################################################################


# ==============================================================================
# BLOC 0 : CONFIGURATION ET PACKAGES
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
})

# Chemins des fichiers 
fichier_qiime    <- "results/all_matam_salmon_qiime_like_table_HYBRIDE.tsv"
fichier_metadata <- "data/16S/metadata.tsv"
out_dir          <- "results/16S/alpha_beta"

# Crée le dossier de sortie s'il n'existe pas encore.
# recursive = TRUE : crée aussi les dossiers parents manquants si nécessaire
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Palette de couleurs pour les taxons 
pal_taxo <- colorspace::qualitative_hcl(14, palette = "Dark 3")

# Couleurs des courbes diversité 
COL_SHANNON <- "#1F78B4"          # diversité alpha
COL_BRAY    <- "#E31A1C"          # diversité beta


# ==============================================================================
# BLOC 1 : FONCTIONS UTILITAIRES
# ==============================================================================
# ------------------------------------------------------------------------------
# lire_qiime_tsv_robuste() : Lecture robuste d'un fichier de comptages QIIME2
# ------------------------------------------------------------------------------
# Le format QIIME2 a une particularité : la ligne d'en-tête commence par "#"
# (ex: "#OTU ID\tSample1\tSample2..."). Cette fonction détecte et nettoie
# automatiquement cet en-tête pour que R puisse lire le fichier correctement,
# quelle que soit la version exacte du format.
#
lire_qiime_tsv_robuste <- function(chemin) {
  lignes <- readLines(chemin, warn = FALSE)

  # Détection de la ligne d'en-tête : on cherche "OTU ID" ou "taxonomy"
  idx_header <- which(grepl("(?i)OTU.?ID", lignes, perl = TRUE))[1]
  if (is.na(idx_header)) {
    idx_header <- which(grepl("(?i)taxonomy", lignes, perl = TRUE))[1]
  }

  lignes_propres <- lignes[idx_header:length(lignes)]

  # Suppression du "#" initial : read_tsv() l'interpréterait comme un commentaire
  lignes_propres[1] <- sub("^#", "", lignes_propres[1])

  # I() encapsule le texte pour que read_tsv() le lise comme s'il s'agissait
  # d'un vrai fichier sur le disque (connexion textuelle)
  df <- read_tsv(I(paste(lignes_propres, collapse = "\n")), show_col_types = FALSE)
  names(df)[1] <- "Taxonomy"
  return(df)
}


# ------------------------------------------------------------------------------
# harmoniser_taxonomie() : Harmonisation des noms de genres bactériens
# ------------------------------------------------------------------------------
# Problème biologique : Les bases de données 16S (SILVA, Greengenes, NCBI) ont
# subi plusieurs révisions taxonomiques majeures ces dernières années.
# Des genres anciens ont été scindés, renommés ou regroupés.
#
# Exemples :
#   - Escherichia et Shigella sont phylogénétiquement indiscernables en 16S
#   - Le genre Lactobacillus a été divisé en ~25 genres distincts dans SILVA 138
#   - Bacillus a été scindé en ~10 genres (Cytobacillus, Priestia, etc.)
#
# Cette fonction applique un dictionnaire de corrections pour garantir des
# noms cohérents entre les différentes versions de SILVA et publications.
#
# Arguments :
#   tab : data.frame contenant au minimum les colonnes Genus et Family
# Retourne :
#   Le même data.frame avec les noms de genres harmonisés
harmoniser_taxonomie <- function(tab) {
  tab %>%
    mutate(
      Genus = case_when(
        # --- Genres indiscernables par le 16S ---
        # Escherichia et Shigella partagent >97% d'identité en 16S
        Genus %in% c("Escherichia", "Shigella",
                     "Escherichia-Shigella", "Escherichia/Shigella") ~ "Escherichia-Shigella",

        # --- Complexes de genres révisés dans SILVA 138 ---
        # Lactobacillus : scindé en ~25 genres (Ligilactobacillus, Lacticaseibacillus...)
        grepl("lactobacillus", Genus, ignore.case = TRUE) ~ "Lactobacillus",

        # Mycobacterium : révision majeure de 2018 (Gupta et al.)
        Genus %in% c("Mycobacterium", "Mycolicibacterium", "Mycolicibacter",
                     "Mycobacteroides", "Mycolicibacillus") ~ "Mycobacterium",

        # Bacillus : scindé en ~10 genres dans SILVA 138
        Genus %in% c("Bacillus", "Cytobacillus", "Mesobacillus", "Neobacillus",
                     "Peribacillus", "Alkalihalobacillus", "Litchfieldia",
                     "Metabacillus", "Priestia") ~ "Bacillus",

        # Pseudomonas : quelques genres satellites réintégrés
        Genus %in% c("Pseudomonas", "Halopseudomonas", "Neopseudomonas") ~ "Pseudomonas",

        # Burkholderia : genre très remanié dans les années 2010
        Genus %in% c("Burkholderia", "Cupriavidus", "Paraburkholderia",
                     "Caballeronia", "Trinickia") ~ "Burkholderia",

        # Mycoplasma : révision taxonomique majeure de 2019 (Gupta et al.)
        Genus %in% c("Mycoplasma", "Mesomycoplasma", "Metamycoplasma",
                     "Mycoplasmopsis", "Mycoplasmoides", "Malacoplasma",
                     "Ureaplasma", "Williamsoniiplasma") ~ "Mycoplasma",

        # Pantoea et genres apparentés (Enterobacteriaceae)
        Genus %in% c("Pantoea", "Erwinia", "Tatumella") ~ "Pantoea",

        # Complexe Enterobacter : genres cryptiques issus de révisions récentes
        Genus %in% c("Mixta", "Dickeya", "Leclercia", "Mangrovibacter",
                     "Enterobacillus", "Jejubacter", "Pluralibacter",
                     "Klebsiella", "Raoultella") ~ "Enterobacter",

        # Reclassifications individuelles
        Genus == "Vagococcus"  ~ "Enterococcus",  # Phylogénétiquement imbriqué
        Genus == "Enhydrobacter" ~ "Moraxella",   # Reclassification 2012

        # Tout le reste : conserver le nom d'origine
        TRUE ~ Genus
      )
    ) %>%
    mutate(
      Family = case_when(
        # Si la Famille est manquante mais le Genre est connu,
        # on crée un nom informatif pour ne pas perdre l'information taxonomique.
        # Cela permet de garder une trace de la position dans l'arbre du vivant.
        Family == "Unassigned" & Genus != "Unassigned" ~ paste0("Famille_de_", Genus),
        TRUE ~ Family
      )
    )
}


# ------------------------------------------------------------------------------
# normaliser_100() : Conversion des comptages bruts en abondances relatives (%)
# ------------------------------------------------------------------------------
# Problème technique : La profondeur de séquençage varie entre échantillons.
# Un échantillon avec 10 000 reads et un avec 50 000 reads ne sont pas
# directement comparables en valeurs brutes.
# Solution : On convertit chaque échantillon en pourcentage (somme = 100%),
# ce qui permet de comparer des compositions, indépendamment du volume total
#
normaliser_100 <- function(df, rang) {
  df %>%
    # Étape 1 : Additionner les reads de chaque taxon par échantillon.
    # Plusieurs lignes peuvent pointer vers le même genre (séquences différentes
    # mais même assignation taxonomique) ici fusionner
    group_by(Sample, !!sym(rang)) %>%   # sym() convertit la chaîne en symbole R
    summarise(Count = sum(Count), .groups = "drop") %>%

    # Étape 2 : Calculer le pourcentage de chaque taxon dans son échantillon
    group_by(Sample) %>%
    mutate(Relative_Abundance = (Count / sum(Count)) * 100) %>%
    ungroup()
}


# ------------------------------------------------------------------------------
# forcer_nomenclature_binomiale() : Construction du nom complet Genre + Espèce
# ------------------------------------------------------------------------------
# En microbiologie, la convention de nomenclature binomiale exige que le nom
# d'une espèce inclue toujours le genre (ex: "Escherichia coli" et non "coli").
# SILVA stocke parfois l'espèce sans le genre → cette fonction assure
# l'uniformité des noms d'espèces à travers tout le script.
#
# Arguments :
#   df : data.frame contenant les colonnes Genus et Species
# Retourne :
#   Le même data.frame avec la colonne Species corrigée
forcer_nomenclature_binomiale <- function(df) {
  df %>%
    mutate(Species = case_when(
      # Cas 1 : Espèce non assignée → on le dit explicitement
      Species %in% c("Non Assigné", "Unassigned") ~ "Unassigned",
      # Cas 2 : Genre non assigné → on ne peut pas construire un nom binomial
      Genus   %in% c("Non Assigné", "Unassigned") ~ Species,
      # Cas 3 : L'espèce contient déjà le genre (ex: "Lactobacillus acidophilus") → OK
      grepl(Genus, Species, ignore.case = TRUE)   ~ Species,
      # Cas 4 : L'espèce est seule (ex: "acidophilus") → on préfixe avec le genre
      TRUE ~ paste(Genus, Species)
    ))
}


# ------------------------------------------------------------------------------
# preparer_barplot_data() : Pipeline complet de préparation pour un barplot
# ------------------------------------------------------------------------------
# Ce pattern est utilisé identiquement pour les larves, les adultes, le global.
# Le mutualiser en une fonction évite 6 blocs de code quasi-identiques.
#
# Arguments :
#   df_brut   : données brutes filtrées sur le stade voulu
#   rang      : "Genus" ou "Species"
#   top_n     : nombre de taxons à afficher individuellement (ex: 13)
#   seuil_min : seuil d'abondance minimale pour entrer dans le Top (ex: 1.0%)
#   cond_levels : vecteur des conditions dans l'ordre chronologique voulu
#   exclus    : vecteur de genres/espèces à exclure du Top (ex: "Unassigned")
# Retourne :
#   Liste avec $df_plot (données pour ggplot) et $palette (couleurs)
preparer_barplot_data <- function(df_brut, rang, top_n, seuil_min = 0,
                                   cond_levels = NULL, exclus = c("Non Assigné", "Unassigned")) {
  # --- 1. Normalisation à 100% par échantillon ---
  df_norm <- normaliser_100(df_brut, rang) %>%
    left_join(meta %>% select(Sample, Condition), by = "Sample")

  # --- 2. Calcul de la vraie moyenne par stade ---
  # La "vraie moyenne" = somme de toutes les abondances relatives / nb de réplicats
  # Cela garantit que la somme des barres = exactement 100% (pas d'artefact de moyenne)
  df_agg <- df_norm %>%
    group_by(Condition) %>%
    mutate(Nb_Samples = n_distinct(Sample)) %>%
    group_by(Condition, !!sym(rang), Nb_Samples) %>%
    summarise(Relative_Abundance = sum(Relative_Abundance) / first(Nb_Samples),
              .groups = "drop")

  if (!is.null(cond_levels)) {
    df_agg <- df_agg %>% mutate(Condition = factor(Condition, levels = cond_levels))
  }

  # --- 3. Identification du Top N des taxons les plus abondants ---
  top_taxa <- df_agg %>%
    filter(!.data[[rang]] %in% exclus) %>%
    group_by(.data[[rang]]) %>%
    summarise(MeanAbund = mean(Relative_Abundance), .groups = "drop") %>%
    filter(MeanAbund >= seuil_min) %>%
    slice_max(MeanAbund, n = top_n, with_ties = FALSE) %>%
    pull(.data[[rang]])

  # --- 4. Regroupement des taxons rares en "Others" ---
  # Les taxons hors Top N sont regroupés pour éviter une légende illisible.
  # "Others" et "Unassigned" sont forcés tout en bas des barres.
  df_plot <- df_agg %>%
    mutate(Taxon_Top = case_when(
      .data[[rang]] %in% exclus ~ "Unassigned",
      .data[[rang]] %in% top_taxa ~ as.character(.data[[rang]]),
      TRUE ~ "Others"
    )) %>%
    group_by(Condition, Taxon_Top) %>%
    summarise(Abundance = sum(Relative_Abundance), .groups = "drop") %>%
    mutate(
      # fct_reorder() trie par abondance totale (les plus abondants en haut)
      # fct_relevel() force ensuite "Others" et "Unassigned" tout en bas
      Taxon_Top = fct_relevel(
        fct_reorder(factor(Taxon_Top), Abundance, sum, .desc = FALSE),
        "Others", "Unassigned",
        after = 0
      )
    )

  # --- 5. Construction de la palette de couleurs ---
  # Les vrais taxons reçoivent des couleurs qualitatives (pal_taxo)
  # Les catégories génériques reçoivent des gris neutres pour ne pas les
  # confondre avec de vrais taxons dans la légende
  levels_taxa <- levels(df_plot$Taxon_Top)
  taxa_vrais  <- setdiff(levels_taxa, c("Others", "Unassigned"))

  palette <- setNames(rep_len(pal_taxo, length(taxa_vrais)), taxa_vrais)
  if ("Others"     %in% levels_taxa) palette["Others"]     <- "grey85"
  if ("Unassigned" %in% levels_taxa) palette["Unassigned"] <- "grey40"

  return(list(df_plot = df_plot, palette = palette))
}


# ------------------------------------------------------------------------------
# creer_barplot() : Création d'un barplot empilé standardisé
# ------------------------------------------------------------------------------
# Arguments :
#   df_plot     : sortie de preparer_barplot_data()$df_plot
#   palette     : sortie de preparer_barplot_data()$palette
#   titre, x_lab, y_lab, fill_lab : textes du graphique
#   angle_x     : angle des étiquettes de l'axe X (défaut = 45°)
#   largeur_bar : largeur des barres (défaut = 0.65)
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
# creer_quadrant_plot() : Graphique Prévalence vs Abondance (Core Microbiome)
# ------------------------------------------------------------------------------
# Ce type de graphique classe chaque bactérie en 4 quadrants selon :
#   - Sa PRÉVALENCE : dans quel % des échantillons est-elle présente ?
#   - Son ABONDANCE moyenne : quel % des séquences représente-t-elle ?
#
# Classification écologique :
#   Q1 (rouge)   : Haute abondance + Haute prévalence → Core microbiome strict
#   Q2 (bleu)    : Faible abondance + Haute prévalence → Core microbiome satellite
#   Q3 (orange)  : Haute abondance + Faible prévalence → Bactérie transitoire
#   Q4 (gris)    : Faible abondance + Faible prévalence → Bruit de fond
#
# Arguments :
#   df_scatter : data.frame avec colonnes Prevalence, MeanAbund, Quadrant, [rang]
#   rang       : colonne à utiliser pour les étiquettes ("Genus" ou "Species")
#   seuil_prev, seuil_abond : lignes de seuils à tracer
#   titre, sous_titre : textes du graphique
#   style_label : "bold" pour les genres, "italic" pour les espèces
creer_quadrant_plot <- function(df_scatter, rang, seuil_prev, seuil_abond,
                                 titre, sous_titre = NULL, style_label = "bold") {

  # Couleurs sémantiques : rouge = important, bleu = commun mais rare,
  # orange = occasionnel mais abondant, gris = négligeable
  palette_quadrant <- c(
    "1. High Abundance Shared Microbiota" = "#E31A1C",
    "2. Low Abundance Shared Microbiota"  = "#1F78B4",
    "3. Transient"                        = "#FF7F00",
    "4. Background Noise"                 = "grey75"
  )

  ggplot(df_scatter, aes(x = Prevalence, y = MeanAbund)) +
    # Lignes de seuils (les "frontières" entre les quadrants)
    geom_vline(xintercept = seuil_prev,  linetype = "dashed", color = "grey40", linewidth = 0.8) +
    geom_hline(yintercept = seuil_abond, linetype = "dashed", color = "grey40", linewidth = 0.8) +

    # Points dont la taille est proportionnelle à l'abondance
    geom_point(aes(color = Quadrant, size = MeanAbund), alpha = 0.8) +

    # Étiquettes intelligentes (ggrepel évite les chevauchements automatiquement)
    # On n'étiquette pas le bruit de fond pour alléger le graphique
    geom_text_repel(
      data = filter(df_scatter, Quadrant != "4. Background Noise"),
      aes(label = .data[[rang]], color = Quadrant),
      size = 3.5, fontface = style_label,
      box.padding = 0.6, max.overlaps = 20, show.legend = FALSE
    ) +

    # Axe Y en échelle logarithmique : permet de voir à la fois les taxons rares
    # (0.01%) et abondants (>10%) sur le même graphique
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
    theme_bw(base_size = 13) +
    theme(
      legend.position = "bottom",
      legend.title    = element_text(face = "bold"),
      legend.text     = element_text(size = 10),
      plot.title      = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 5), nrow = 2))
}

# Opérateur "null-coalesce" : renvoie le 2e argument si le 1er est NULL
# (utilisé dans creer_quadrant_plot pour le sous-titre par défaut)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ------------------------------------------------------------------------------
# calculer_core_scatter() : Calcul des métriques Prévalence + Abondance moyenne
# ------------------------------------------------------------------------------
# Pour chaque taxon, on calcule :
#   - MeanAbund   : abondance relative moyenne sur TOUS les échantillons du groupe
#   - Prevalence  : % d'échantillons où ce taxon est détecté (abondance > 0)
#
# Arguments :
#   df_norm   : données normalisées à 100% (sortie de normaliser_100)
#   rang      : "Genus" ou "Species"
#   n_samples : nombre total d'échantillons dans le groupe
#   seuil_prev, seuil_abond : seuils pour la classification en quadrants
# Retourne :
#   data.frame avec colonnes : [rang], MeanAbund, Prevalence, Quadrant
calculer_core_scatter <- function(df_norm, rang, n_samples, seuil_prev, seuil_abond) {
  df_norm %>%
    filter(!.data[[rang]] %in% c("Non Assigné", "Unassigned")) %>%
    group_by(.data[[rang]]) %>%
    summarise(
      # Abondance moyenne : somme sur tous les échantillons / nb total d'échantillons
      # (les échantillons où la bactérie est absente contribuent 0 au numérateur)
      MeanAbund  = sum(Relative_Abundance) / n_samples,
      # Prévalence : proportion d'échantillons où l'abondance est strictement > 0
      Prevalence = sum(Relative_Abundance > 0) / n_samples * 100,
      .groups = "drop"
    ) %>%
    # On retire les taxons quasi-absents (<0.01%) pour alléger le graphique
    filter(MeanAbund >= 0.01) %>%
    mutate(Quadrant = case_when(
      Prevalence >= seuil_prev & MeanAbund >= seuil_abond ~ "1. High Abundance Shared Microbiota",
      Prevalence >= seuil_prev & MeanAbund <  seuil_abond ~ "2. Low Abundance Shared Microbiota",
      Prevalence <  seuil_prev & MeanAbund >= seuil_abond ~ "3. Transient",
      Prevalence <  seuil_prev & MeanAbund <  seuil_abond ~ "4. Background Noise"
    ))
}


# ==============================================================================
# BLOC 2 : CHARGEMENT ET FILTRAGE DES DONNÉES
# ==============================================================================
message("--- Chargement des données ---")

# --- Dictionnaires de traduction ---
# MATAM et les métadonnées utilisent des codes courts en français.
# On les traduit en noms anglais complets avec les poids (mg) pour que les
# graphiques soient directement exploitables pour une publication.
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

# --- Lecture de la table de comptages ---
# rename_with() renomme les colonnes échantillons selon le dictionnaire trad_samples
# -Taxonomy : la colonne Taxonomy est exclue du renommage
df_qiime_16s <- lire_qiime_tsv_robuste(fichier_qiime) %>%
  rename_with(~ recode(.x, !!!trad_samples), -Taxonomy)

# --- Lecture des métadonnées ---
meta <- read_tsv(fichier_metadata, show_col_types = FALSE) %>%
  rename(Sample = `sample-id`, Condition = condition) %>%
  mutate(
    Condition = recode(Condition, !!!trad_conditions),
    Sample    = recode(Sample,    !!!trad_samples)
  )

# --- Extraction des rangs taxonomiques par regex ---
# La colonne Taxonomy contient des chaînes de la forme :
#   "d__Bacteria;p__Firmicutes;c__Bacilli;...;g__Lactobacillus;s__acidophilus"
# str_extract() capture tout ce qui suit le préfixe jusqu'au prochain ";"
df_taxo_16s <- df_qiime_16s %>%
  # Passage du format large (1 colonne/échantillon) au format long (1 ligne/taxon×échantillon)
  pivot_longer(cols = -Taxonomy, names_to = "Sample", values_to = "Count") %>%

  # Optimisation : on ne garde que les lignes avec au moins 1 read détecté
  filter(Count > 0) %>%

  # Filtrage biologique crucial : les chloroplastes et mitochondries possèdent
  # un gène ARNr 16S homologue aux bactéries. Ils contaminent systématiquement
  # les extractions ADN végétales ou animales et ne représentent pas le
  # microbiote intestinal réel de l'insecte → on les retire.
  filter(!grepl("Chloroplast|Mitochondria|Mitochondrion", Taxonomy, ignore.case = TRUE)) %>%

  # Extraction de chaque rang taxonomique par son préfixe SILVA
  mutate(
    Phylum  = str_extract(Taxonomy, "p__[^;]+"),
    Class   = str_extract(Taxonomy, "c__[^;]+"),
    Order   = str_extract(Taxonomy, "o__[^;]+"),
    Family  = str_extract(Taxonomy, "f__[^;]+"),
    Genus   = str_extract(Taxonomy, "g__[^;]+"),
    Species = str_extract(Taxonomy, "s__[^;]+")
  ) %>%

  # Nettoyage des préfixes ("g__", "p__", etc.)
  mutate(across(c(Phylum, Class, Order, Family, Genus, Species),
                ~ str_remove(.x, "^[a-z]__"))) %>%

  # Remplacement des NA par "Unassigned" : préférable à NA pour les filtres
  mutate(across(c(Phylum, Class, Order, Family, Genus, Species),
                ~ ifelse(is.na(.x) | str_trim(.x) == "", "Unassigned", .x)))

# --- Assemblage du tableau final ---
df_final_16s <- harmoniser_taxonomie(df_taxo_16s) %>%
  left_join(meta, by = "Sample") %>%
  filter(!is.na(Condition))   # Écarter les échantillons sans métadonnées

# --- Ordre chronologique du développement de T. molitor ---
# CET ORDRE EST CRUCIAL : il gouverne l'axe X de tous les graphiques temporels.
# T. molitor est un insecte à métamorphose complète (holométabole) :
#   Œuf → Larves (plusieurs stades) → Nymphe → Adulte
chronologie_insecte <- c(
  "Tiny_Larvae (2 mg)", "Microlarvae (7 mg)",
  "Larvae_S1 (14 mg)", "Larvae_S2 (40 mg)", "Larvae_S3 (65 mg)", "Larvae_S4 (100 mg)",
  "Pupae", "Young_Beetles", "Beetles"
)

# Version étendue incluant les substrats alimentaires (environnement)
chronologie_tout <- c(chronologie_insecte, "Raw_Substrate", "Final_Substrate")

# Tableau filtré sur toutes les conditions (insecte + substrats)
# NOTE : JA1 et N3 (échantillons pathologiques) sont conservés intentionnellement.
df_insecte <- df_final_16s %>%
  filter(Condition %in% chronologie_tout) %>%
  mutate(
    Condition    = factor(Condition, levels = chronologie_tout),
    # ConditionNum : version numérique du facteur pour l'axe X continu de ggplot2
    ConditionNum = as.numeric(Condition)
  )

# Palette de couleurs par stade de développement
# Dégradé du vert clair (petites larves) au bleu foncé (larves âgées) → violet (nymphe) → rouge (adulte)
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
# BLOC 3 : CALCULS ALPHA ET BÊTA-DIVERSITÉ
# ==============================================================================

# ==============================================================================
# A. ALPHA-DIVERSITÉ : INDICE DE SHANNON
# ==============================================================================
# L'ALPHA-DIVERSITÉ mesure la diversité AU SEIN de chaque échantillon.
# L'indice de Shannon (H') tient compte à la fois du nombre d'espèces
# (richesse) et de leur équitabilité (régularité de la distribution).
# H' élevé = communauté diverse et équilibrée
# H' faible = communauté dominée par quelques taxons
message("--- Calcul de l'alpha-diversité (Shannon) ---")

# Construction de la matrice Échantillons × Genres (format requis par vegan)
mat_counts_alpha <- df_insecte %>%
  group_by(Sample, Genus) %>%
  summarise(Count = sum(Count), .groups = "drop") %>%
  pivot_wider(names_from = Genus, values_from = Count, values_fill = 0) %>%
  column_to_rownames("Sample")

# Calcul de Shannon sur les COMPTAGES BRUTS (pas les abondances normalisées)
# La fonction diversity() de vegan gère la normalisation en interne
shannon_vals <- diversity(mat_counts_alpha, index = "shannon")

# Résumé par stade de développement (moyenne ± écart-type)
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
    SD_Shannon  = ifelse(n() > 1, sd(Shannon), 0),   # SD = 0 si 1 seul réplicat
    .groups = "drop"
  )


# ==============================================================================
# B. BÊTA-DIVERSITÉ : DISTANCE DE BRAY-CURTIS
# ==============================================================================
# La BÊTA-DIVERSITÉ mesure la DISSIMILARITÉ entre échantillons.
# La distance de Bray-Curtis varie de 0 (communautés identiques) à 1
# (pas d'espèces en commun). On la calcule ici ENTRE RÉPLICATS d'un même stade
# pour mesurer l'hétérogénéité intra-groupe : si les réplicats d'un stade
# sont très différents entre eux, la bêta-diversité intra-groupe sera élevée.
message("--- Calcul de la bêta-diversité (Bray-Curtis) ---")

resultats_stade <- list()   # Résultats agrégés par stade (moyenne ± SD)
resultats_ind   <- list()   # Résultats par réplicat individuel

for (cond in levels(df_insecte$Condition)) {
  df_cond <- df_insecte %>% filter(Condition == cond)
  samps   <- unique(df_cond$Sample)

  # Il faut au minimum 2 échantillons pour calculer une distance pairwise
  if (length(samps) < 2) next

  # Normalisation à 100% (les "Unassigned" sont retirés avant normalisation
  # pour ne pas biaiser la composition relative)
  df_norm <- normaliser_100(df_cond %>% filter(Genus != "Unassigned"), "Genus")

  # Conversion en matrice Échantillons × Genres (format vegan)
  mat_wide <- df_norm %>%
    pivot_wider(id_cols = Sample, names_from = Genus,
                values_from = Relative_Abundance, values_fill = 0) %>%
    column_to_rownames("Sample")

  # Calcul de la matrice de distances pairwise (tous réplicats 2 à 2)
  dist_mat <- as.matrix(vegdist(mat_wide, method = "bray"))

  # upper.tri() : on ne prend que le triangle supérieur pour éviter les doublons
  vals <- dist_mat[upper.tri(dist_mat)]

  resultats_stade[[cond]] <- data.frame(
    Condition = cond,
    Beta_Mean = mean(vals),
    Beta_SD   = ifelse(length(vals) > 1, sd(vals), 0)
  )

  # Pour chaque réplicat : distance moyenne aux autres réplicats du même stade
  diag(dist_mat) <- NA   # La diagonale = 0 (distance à soi-même) → on l'exclut
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
# C. TEST DE CORRÉLATION DE SPEARMAN : SHANNON vs BRAY-CURTIS
# ==============================================================================
# On teste si les échantillons les plus divers (Shannon élevé) sont aussi
# les plus homogènes entre réplicats (Bray-Curtis faible) ou l'inverse.
# Le test de Spearman est non-paramétrique → adapté à des distributions
# non-normales comme les données de microbiome.
message("--- Test de corrélation de Spearman ---")

# Tableau individuel : 1 ligne par réplicat avec Shannon ET Bray-Curtis
df_ind <- data.frame(Sample = names(shannon_vals), Shannon = shannon_vals) %>%
  left_join(df_beta_ind, by = "Sample") %>%
  filter(!is.na(Beta_Ind))

# Test 1 : SANS les substrats (insecte seul)
df_ind_sans   <- df_ind %>% filter(!grepl("Substrat|Substrate", Condition))
test_stat_sans <- cor.test(df_ind_sans$Shannon, df_ind_sans$Beta_Ind, method = "spearman", exact = FALSE)
rho_sans  <- round(test_stat_sans$estimate, 2)
p_val_sans <- paste0("= ", format(test_stat_sans$p.value, scientific = TRUE, digits = 3))

# Test 2 : AVEC les substrats (tout le jeu de données)
test_stat_avec <- cor.test(df_ind$Shannon, df_ind$Beta_Ind, method = "spearman", exact = FALSE)
rho_avec  <- round(test_stat_avec$estimate, 2)
p_val_avec <- paste0("= ", format(test_stat_avec$p.value, scientific = TRUE, digits = 3))

# Fusion alpha + bêta par stade pour les graphiques
df_final <- df_alpha_stade %>%
  left_join(df_beta_stade, by = "Condition") %>%
  mutate(
    # Trajectoire : différencie insecte et substrats pour la facette du graphique 2
    Trajectoire = ifelse(grepl("Substrat|Substrate", Condition), "Environnement", "Développement de l'Insecte"),
    Trajectoire = factor(Trajectoire, levels = c("Développement de l'Insecte", "Environnement"))
  )


# ==============================================================================
# BLOC 4 : EXPORT DES RÉSULTATS EN FORMAT TSV
# ==============================================================================
# On exporte un tableau synthétique avec les valeurs numériques brutes.
# Utile pour : traçabilité, réanalyse future, ou fourniture en annexe.
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
# BLOC 5 : GRAPHIQUE DOUBLE-AXE (SHANNON + BRAY-CURTIS)
# ==============================================================================
# Ce graphique superpose deux métriques sur deux axes Y distincts.
# L'astuce du coefficient (coeff) permet d'aligner les deux courbes
# sur la même échelle visuelle tout en conservant des axes indépendants.
#
# Coefficient = max(Shannon) / max(Bray-Curtis)
# → Toutes les valeurs Bray-Curtis sont multipliées par ce coefficient
#   avant d'être tracées sur l'axe Y de Shannon.
# → L'axe Y droit est re-scalé par l'inverse (division) pour afficher
#   les vraies valeurs Bray-Curtis.
message("--- Génération des Graphiques Alpha/Bêta ---")

# Fonction interne pour créer un graphique double-axe
creer_graphique_double_axe <- function(df_data, coeff, chronologie_labels,
                                        rho, p_val, titre, titre_x,
                                        df_facette = NULL) {
  # Position de l'annotation statistique (coin bas-gauche)
  df_label <- data.frame(
    ConditionNum = 1,
    y_pos = max(df_data$MeanShannon, na.rm = TRUE) * 0.1,
    label = paste0("rho = ", rho, "\np ", p_val)
  )
  if (!is.null(df_facette)) df_label$Trajectoire <- df_facette

  p <- ggplot(df_data, aes(x = ConditionNum)) +

    # Barres d'erreur (écart-type) pour les deux métriques
    geom_errorbar(aes(ymin = (Beta_Mean - Beta_SD) * coeff,
                      ymax = (Beta_Mean + Beta_SD) * coeff),
                  color = COL_BRAY, width = 0.2, linewidth = 0.8) +
    geom_errorbar(aes(ymin = MeanShannon - SD_Shannon,
                      ymax = MeanShannon + SD_Shannon),
                  color = COL_SHANNON, width = 0.2, linewidth = 0.8) +

    # Courbe Bray-Curtis (pointillée + triangles)
    geom_line(aes(y = Beta_Mean * coeff, color = "Heterogeneity (Bray-Curtis)"),
              linewidth = 1.5, linetype = "dashed") +
    geom_point(aes(y = Beta_Mean * coeff, color = "Heterogeneity (Bray-Curtis)"),
               size = 4, shape = 17) +

    # Courbe Shannon (continue + ronds)
    geom_line(aes(y = MeanShannon, color = "Diversity (Shannon)"),
              linewidth = 1.5) +
    geom_point(aes(y = MeanShannon, color = "Diversity (Shannon)"),
               size = 5, shape = 16) +

    # Annotation du test de corrélation de Spearman
    geom_label(data = df_label,
               aes(x = ConditionNum, y = y_pos, label = label),
               fill = "white", fontface = "bold", size = 3.8,
               hjust = 0, color = COL_BRAY, inherit.aes = FALSE) +

    # Axe X : stades de développement
    scale_x_continuous(breaks = seq_along(chronologie_labels),
                       labels = chronologie_labels) +

    # Double axe Y : Shannon (gauche) et Bray-Curtis re-scalé (droite)
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
      axis.text.x        = element_text(angle = 35, hjust = 1, face = "bold"),
      axis.title.y.left  = element_text(color = COL_SHANNON, face = "bold"),
      axis.title.y.right = element_text(color = COL_BRAY, face = "bold"),
      legend.position    = "bottom",
      panel.grid.minor   = element_blank()
    )

  return(p)
}

# --- Graphique 1 : Insecte seul (sans substrats) ---
df_sans    <- df_final %>% filter(Trajectoire == "Développement de l'Insecte")
coeff_sans <- max(df_sans$MeanShannon, na.rm = TRUE) / max(df_sans$Beta_Mean, na.rm = TRUE)

p_dyn_sans <- creer_graphique_double_axe(
  df_data           = df_sans,
  coeff             = coeff_sans,
  chronologie_labels = chronologie_insecte,
  rho    = rho_sans,
  p_val  = p_val_sans,
  titre  = "Dynamique du Microbiote",
  titre_x = "Stades de développement"
)
ggsave(file.path(out_dir, "shannon_vs_braycurtis_.png"),
       plot = p_dyn_sans, width = 12, height = 8, dpi = 300)

# --- Graphique 2 : Insecte + Substrats (avec facettes) ---
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

ggsave(file.path(out_dir, "shannon_vs_braycurtis_substrats.png"),
       plot = p_dyn_avec, width = 12, height = 8, dpi = 300)


# ==============================================================================
# BLOC 6 : BARPLOT GLOBAL 16S — COMPOSITION PAR RÉPLICAT (NIVEAU GENRE)
# ==============================================================================
# Ce graphique montre la composition en genres bactériens de CHAQUE réplicat,
# organisée par condition. Il permet d'évaluer la reproductibilité entre
# réplicats et de voir l'évolution de la communauté à travers le développement.
message("--- Barplot Global 16S ---")

ordre_chronologique <- chronologie_tout  # Alias pour la cohérence des noms

# Normalisation à 100% sur l'ensemble du jeu de données
df_global_norm <- normaliser_100(df_final_16s, "Genus") %>%
  left_join(meta %>% select(Sample, Condition), by = "Sample") %>%
  filter(!is.na(Condition)) %>%
  mutate(Condition = factor(Condition, levels = ordre_chronologique))

# Identification du Top 14 des genres les plus abondants globalement.
# Certains genres écologiquement non pertinents pour T. molitor sont exclus.
top_global <- df_global_norm %>%
  filter(!(Genus %in% c("Unassigned", "Bacteroides", "Alistipes",
                         "Aquabacterium", "Blastococcus"))) %>%
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

# Palette du barplot global
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
  # facet_grid crée un panneau par condition
  # scales = "free_x" + space = "free_x" : chaque panneau a sa propre largeur
  # proportionnelle au nombre de réplicats (évite des barres de largeurs inégales)
  facet_grid(~ Condition, scales = "free_x", space = "free_x") +
  labs(
    title = "Catalogue Global du Microbiote Intestinal (16S)",
    x = "Réplicats Biologiques", y = "Abondance Relative (%)",
    fill = "Genre Bactérien"
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


# ==============================================================================
# BLOC 7 : HEATMAP DE DISSIMILARITÉ BRAY-CURTIS (TOUS RÉPLICATS)
# ==============================================================================
# Cette matrice montre la dissimilarité Bray-Curtis entre CHAQUE paire de
# réplicats (toutes conditions confondues). Les blocs diagonaux (cases vertes)
# indiquent que les réplicats d'un même stade sont similaires entre eux.
# Des blocs hors-diagonale verts indiqueraient des stades partageant une
# composition microbienne commune.
message("--- Heatmap Bray-Curtis globale ---")

# Fonction interne pour créer une heatmap de distance
creer_heatmap_bray <- function(df_dist_long, ordre, titre, largeur = 12, hauteur = 10,
                                nom_fichier) {
  df_plot <- df_dist_long %>%
    filter(Sample1 %in% ordre & Sample2 %in% ordre) %>%
    mutate(
      Sample1 = factor(Sample1, levels = ordre),
      # On inverse l'ordre de Sample2 pour que la diagonale soit orientée
      # de haut-gauche vers bas-droite (convention standard des matrices)
      Sample2 = factor(Sample2, levels = rev(ordre))
    )

  p <- ggplot(df_plot, aes(x = Sample1, y = Sample2, fill = Distance)) +
    geom_tile(color = "white", linewidth = 0.3) +
    # Dégradé sémantique : Vert (0 = identique) → Jaune → Rouge (1 = différent)
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

# Préparation de la matrice de distance complète
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

# Passage en format long pour ggplot2
df_dist_long <- as.data.frame(dist_mat_complete) %>%
  rownames_to_column("Sample1") %>%
  pivot_longer(cols = -Sample1, names_to = "Sample2", values_to = "Distance")

# Ordonnancement chronologique des réplicats
ordre_echantillons <- df_global_norm %>%
  distinct(Sample, Condition) %>%
  arrange(Condition, Sample) %>%
  pull(Sample)

# Heatmap 1 : Tous les réplicats (insecte + substrats)
creer_heatmap_bray(
  df_dist_long, ordre_echantillons,
  titre = "Inter-replicates dissimilarity matrix (Bray-Curtis)",
  nom_fichier = "Matrice_BrayCurtis_Tous_Replicats.png"
)

# Heatmap 2 : Insecte seul (sans substrats)
samples_insecte <- meta %>% filter(Condition %in% chronologie_insecte) %>% pull(Sample)
ordre_insecte   <- ordre_echantillons[ordre_echantillons %in% samples_insecte]

creer_heatmap_bray(
  df_dist_long, ordre_insecte,
  titre = "Inter-replicates dissimilarity matrix (Insect only)",
  largeur = 10, hauteur = 8,
  nom_fichier = "Matrice_BrayCurtis_Insecte_Seul.png"
)


# ==============================================================================
# BLOC 8 : FIGURE MAÎTRESSE (MASTER FIGURE)
# ==============================================================================
# Cette figure combine en un seul graphique :
#   - Un en-tête avec les noms des stades et des trajectoires
#   - Les courbes Shannon + Bray-Curtis (milieu)
#   - Le barplot de composition taxonomique (bas)
# Les trois panneaux partagent le même axe X millimétré, aligné à la colonne
# de chaque réplicat. L'assemblage est réalisé avec le package patchwork.
message("--- Génération de la Figure Maîtresse ---")

generer_master_ultime <- function(stades_filtre, nom_fichier) {

  # --- 1. POSITIONNEMENT MILLIMÉTRÉ DES RÉPLICATS ---
  # On calcule une position X numérique pour chaque réplicat, avec des espaces
  # plus larges entre les grandes catégories (insecte vs environnement)
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
      # Espace large (1.5) entre trajectoires, normal (0.8) entre stades
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

  # --- 2. CALCUL DES COLONNES ET BANDEAUX ---
  df_conds <- df_samples %>%
    group_by(Condition, Trajectoire) %>%
    summarise(x_mid = mean(x_pos),
              x_min = min(x_pos) - 0.48,
              x_max = max(x_pos) + 0.48, .groups = "drop")

  # Espaces entre groupes de conditions (zones grises de séparation)
  df_gaps <- data.frame(
    xmin = head(df_conds$x_max, -1),
    xmax = tail(df_conds$x_min, -1)
  )

  # Bandeaux des trajectoires (s'étendent sur plusieurs conditions)
  df_traj <- df_conds %>%
    group_by(Trajectoire) %>%
    summarise(x_min = min(x_min), x_max = max(x_max),
              x_mid = mean(c(min(x_min), max(x_max))), .groups = "drop")

  LIMITES_X <- c(min(df_conds$x_min) - 0.1, max(df_conds$x_max) + 0.1)

  # --- 3. STATISTIQUES DE CORRÉLATION SUR CE SOUS-ENSEMBLE ---
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

  # --- 4. PANNEAU EN-TÊTE (Bandeaux Trajectoire + Condition) ---
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

  # --- 5. PANNEAU DU MILIEU (Courbes Shannon + Bray-Curtis) ---
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

  # --- 6. PANNEAU DU BAS (Barplots de composition taxonomique) ---
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

  # --- 7. ASSEMBLAGE AVEC PATCHWORK ---
  # Les trois panneaux sont empilés verticalement.
  # heights = c(0.12, 1, 1.25) : l'en-tête est très compact,
  # le barplot légèrement plus grand que les courbes.
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
# BLOC 9 : CORE MICROBIOTE LARVAIRE (BARPLOTS GENRE + ESPÈCE)
# ==============================================================================
# Le "core microbiome" d'un groupe de stades est l'ensemble des bactéries
# régulièrement et abondamment présentes à ces stades.
# On le représente ici comme une moyenne par stade (barplot sur l'axe X = stade).
message("--- Core Microbiote Larvaire ---")

# Stades larvaires d'intérêt (on exclut les pré-larves trop jeunes et les nymphes)
cond_larves <- c("Microlarvae (7 mg)", "Larvae_S1 (14 mg)", "Larvae_S2 (40 mg)",
                  "Larvae_S3 (65 mg)", "Larvae_S4 (100 mg)")
df_larves_brut <- df_final_16s %>% filter(Condition %in% cond_larves)

# --- Partie A : Niveau Genre ---
res_larves_gen <- preparer_barplot_data(
  df_brut = df_larves_brut, rang = "Genus",
  top_n = 13, seuil_min = 1.0, cond_levels = cond_larves
)
p_core <- creer_barplot(
  res_larves_gen$df_plot, res_larves_gen$palette,
  titre    = "Core Larval Microbiome (16S)",
  x_lab    = "Larval Stages",
  fill_lab = "Bacterial Genus:"
)
ggsave(file.path(out_dir, "Core_Larval_Microbiome_Means.png"),
       plot = p_core, width = 9, height = 7, dpi = 300)

# --- Partie B : Niveau Espèce ---
df_larves_sp <- forcer_nomenclature_binomiale(df_larves_brut)
res_larves_sp <- preparer_barplot_data(
  df_brut = df_larves_sp, rang = "Species",
  top_n = 13, seuil_min = 1.0, cond_levels = cond_larves
)
p_core_sp <- creer_barplot(
  res_larves_sp$df_plot, res_larves_sp$palette,
  titre    = "Core Larval Microbiome (16S - Species level)",
  x_lab    = "Larval Stages",
  fill_lab = "Bacterial Species:"
)
ggsave(file.path(out_dir, "Core_Larval_Microbiome_Means_Species.png"),
       plot = p_core_sp, width = 9.5, height = 7, dpi = 300)
message("=> Core Larval Microbiome barplots générés.")


# ==============================================================================
# BLOC 10 : GRAPHIQUES QUADRANTS (PRÉVALENCE vs ABONDANCE) — LARVES
# ==============================================================================
# Ce graphique place chaque genre (ou espèce) dans un plan à deux dimensions :
#   - Axe X : Prévalence (% des larves où ce genre est détecté)
#   - Axe Y : Abondance relative moyenne (en échelle log)
# Les lignes de seuils divisent l'espace en 4 quadrants écologiques.
message("--- Graphiques Quadrants Larves ---")

SEUIL_PREVALENCE <- 80   # Présent dans ≥ 80% des larves → "partagé"
SEUIL_ABONDANCE  <- 1.0  # Représente ≥ 1% de la communauté → "abondant"

# --- Genre ---
df_larves_norm_gen <- normaliser_100(df_larves_brut, "Genus") %>%
  left_join(meta %>% select(Sample, Condition), by = "Sample")
n_larves <- n_distinct(df_larves_norm_gen$Sample)

df_core_scatter_gen <- calculer_core_scatter(
  df_larves_norm_gen, "Genus", n_larves, SEUIL_PREVALENCE, SEUIL_ABONDANCE
)
p_scatter_gen <- creer_quadrant_plot(
  df_core_scatter_gen, "Genus", SEUIL_PREVALENCE, SEUIL_ABONDANCE,
  titre = "Larval Core Microbiome Structure (16S - Genus level)"
)
ggsave(file.path(out_dir, "Core_Microbiome_Quadrants.png"),
       plot = p_scatter_gen, width = 11, height = 8, dpi = 300)

# --- Espèce ---
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
ggsave(file.path(out_dir, "Core_Microbiome_Quadrants_Species.png"),
       plot = p_scatter_sp, width = 11, height = 8, dpi = 300)


# ==============================================================================
# BLOC 11 : MATRICES DE FIDÉLITÉ INTRA-STADE
# ==============================================================================
# Pour les stades biologiquement "critiques" (très jeunes larves, nymphes,
# jeunes adultes), on génère une heatmap montrant l'abondance de chaque genre
# dans chacun des 4 réplicats.
# On ne retient que les genres présents dans ≥ 3 réplicats sur 4 (fidélité).
# Cela permet d'évaluer la reproductibilité biologique à chaque stade.
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
    # Critère de fidélité : présent dans au moins 3 réplicats sur 4
    mutate(Nb_Present = sum(Relative_Abundance > 0)) %>%
    filter(Nb_Present >= 3) %>%
    ungroup()

  if (nrow(df_stade) > 0) {

    df_stade <- df_stade %>%
      mutate(Sample = factor(Sample, levels = reps_attendus)) %>%
      # complete() génère les cases manquantes (genre absent d'un réplicat) avec 0
      # → elles s'affichent en gris clair sur la heatmap
      complete(Sample, Genus, fill = list(Relative_Abundance = 0, Condition = stade)) %>%
      group_by(Genus) %>%
      mutate(Total_Abund = sum(Relative_Abundance)) %>%
      ungroup() %>%
      mutate(
        Genus      = fct_reorder(Genus, Total_Abund),
        # N'affiche le texte que si l'abondance est ≥ 0.05% (évite la surcharge visuelle)
        Text_Label = ifelse(Relative_Abundance >= 0.05,
                            sprintf("%.1f", Relative_Abundance), ""),
        # Texte blanc sur fond foncé (>40%), noir sur fond clair
        Is_Dark    = Relative_Abundance > 40
      )

    p <- ggplot(df_stade, aes(x = Sample, y = Genus, fill = Relative_Abundance)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = Text_Label, color = Is_Dark),
                size = 3.2, fontface = "bold", show.legend = FALSE) +
      scale_color_manual(values = c("TRUE" = "white", "FALSE" = "grey20")) +
      # Gradient non-linéaire : 0=gris, 1%=bleu clair, >20%=bleu foncé
      # rescale() permet de "zoomer" sur la plage 0-20% où sont la plupart des données
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
  ggsave(filename = file.path(out_dir, "Core_genome_heterogeneous.png"),
         plot = p_combined, width = 9.5, height = final_h, dpi = 300, bg = "white")
}


# ==============================================================================
# BLOC 12 : CORE MICROBIOTE ADULTE (BARPLOTS + DONUTS, GENRE + ESPÈCE)
# ==============================================================================
# Pour les adultes (Beetles), on génère deux types de visualisations :
#   - Barplot empilé (identique au core larvaire)
#   - Graphique en donut (camembert avec trou central) : particulièrement
#     adapté pour un seul stade car il montre les proportions sous forme circulaire
message("--- Core Microbiote Adulte ---")

cond_adulte    <- c("Beetles")
df_adulte_brut <- df_final_16s %>% filter(Condition %in% cond_adulte)

# Fonction interne pour créer un donut chart
creer_donut <- function(df_plot, palette, titre) {
  df_donut <- df_plot %>%
    arrange(Taxon_Top) %>%
    mutate(
      ymax          = cumsum(Abundance),
      ymin          = lag(ymax, default = 0),
      labelPosition = (ymax + ymin) / 2,   # Centre de chaque secteur
      label         = paste0(Taxon_Top, "\n", round(Abundance, 1), "%")
    )

  ggplot(df_donut, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2.5, fill = Taxon_Top)) +
    geom_rect(color = "white", linewidth = 0.5) +
    # N'affiche les étiquettes que pour les secteurs > 2% (évite les chevauchements)
    geom_text(data = filter(df_donut, Abundance > 2.0),
              aes(x = 4.8, y = labelPosition, label = label),
              size = 3.2, fontface = "bold", lineheight = 0.8) +
    # coord_polar transforme le graphique rectangulaire en donut
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

# --- Partie Genre ---
res_adulte_gen <- preparer_barplot_data(
  df_brut = df_adulte_brut, rang = "Genus",
  top_n = 12, seuil_min = 1.0, cond_levels = cond_adulte
)
p_core_adulte <- creer_barplot(
  res_adulte_gen$df_plot, res_adulte_gen$palette,
  titre = "Core Adult Microbiome (16S)", x_lab = NULL,
  fill_lab = "Bacterial Genus:", angle_x = 0, largeur_bar = 0.45
) +
  theme(axis.text.x = element_text(angle = 0, face = "bold", size = 12))

p_donut_adulte <- creer_donut(
  res_adulte_gen$df_plot, res_adulte_gen$palette,
  titre = "Adult Stage Identity Card (16S)"
)
ggsave(file.path(out_dir, "Core_Adult_Microbiome_Means.png"),
       plot = p_core_adulte, width = 6.5, height = 7, dpi = 300)
ggsave(file.path(out_dir, "Core_Adult_Microbiome_Donut.png"),
       plot = p_donut_adulte, width = 9.5, height = 7, dpi = 300)

# --- Partie Espèce ---
df_adulte_sp <- forcer_nomenclature_binomiale(df_adulte_brut)
res_adulte_sp <- preparer_barplot_data(
  df_brut = df_adulte_sp, rang = "Species",
  top_n = 12, seuil_min = 1.0, cond_levels = cond_adulte
)
p_core_adulte_sp <- creer_barplot(
  res_adulte_sp$df_plot, res_adulte_sp$palette,
  titre = "Core Adult Microbiome (16S — Species level)", x_lab = NULL,
  fill_lab = "Bacterial Species:", angle_x = 0, largeur_bar = 0.45
) +
  theme(axis.text.x = element_text(angle = 0, face = "bold", size = 12))

p_donut_adulte_sp <- creer_donut(
  res_adulte_sp$df_plot, res_adulte_sp$palette,
  titre = "Adult Stage Identity Card (16S — Species level)"
)
ggsave(file.path(out_dir, "Core_Adult_Microbiome_Means_Species.png"),
       plot = p_core_adulte_sp, width = 7.5, height = 7, dpi = 300)
ggsave(file.path(out_dir, "Core_Adult_Microbiome_Donut_Species.png"),
       plot = p_donut_adulte_sp, width = 10.5, height = 7, dpi = 300)


# ==============================================================================
# BLOC 13 : GRAPHIQUES QUADRANTS (PRÉVALENCE vs ABONDANCE) — ADULTES
# ==============================================================================
message("--- Graphiques Quadrants Adultes ---")

# Les seuils pour les adultes sont légèrement moins stricts que pour les larves :
# on a moins de réplicats, et la diversité biologique adulte est plus restreinte.
SEUIL_PREVALENCE_AD <- 75   # Présent dans ≥ 3 réplicats sur 4 = 75%
SEUIL_ABONDANCE_AD  <- 1.0  # ≥ 1% d'abondance relative moyenne

# --- Genre ---
df_ad_norm_gen <- normaliser_100(df_adulte_brut, "Genus") %>%
  left_join(meta %>% select(Sample, Condition), by = "Sample") %>%
  filter(!Genus %in% c("Non Assigné", "Unassigned"))
n_ad <- n_distinct(df_ad_norm_gen$Sample)

df_core_ad_gen <- calculer_core_scatter(
  df_ad_norm_gen, "Genus", n_ad, SEUIL_PREVALENCE_AD, SEUIL_ABONDANCE_AD
)
p_scatter_ad <- creer_quadrant_plot(
  df_core_ad_gen, "Genus", SEUIL_PREVALENCE_AD, SEUIL_ABONDANCE_AD,
  titre = "Adult Core Microbiome Structure (Beetles - Genus level)"
)
ggsave(file.path(out_dir, "Core_Adult_Microbiome_Quadrants.png"),
       plot = p_scatter_ad, width = 11, height = 8, dpi = 300)

# --- Espèce ---
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
  titre = "Adult Core Microbiome Structure (Beetles - Species level)",
  style_label = "italic"
)
ggsave(file.path(out_dir, "Core_Adult_Microbiome_Quadrants_Species.png"),
       plot = p_scatter_ad_sp, width = 11, height = 8, dpi = 300)


# ==============================================================================
# BLOC 14 : ACP (ANALYSE EN COMPOSANTES PRINCIPALES)
# ==============================================================================
# L'ACP résume la structure de la communauté microbienne en 2 dimensions.
# Des échantillons proches sur l'ACP ont des compositions microbiennes similaires.
# Des groupes de points séparés indiquent des communautés distinctes.
#
# La TRANSFORMATION DE HELLINGER est indispensable avant une ACP sur des données
# de comptages :
#   - Réduit le poids des taxons dominants (qui sinon écraseraient tout le signal)
#   - Atténue le "problème du double zéro" (deux échantillons sans un taxon rare
#     ne doivent pas être considérés comme similaires pour autant)
message("--- ACP Classique (Hellinger) ---")

df_pca_norm <- normaliser_100(df_final_16s %>% filter(Genus != "Unassigned"), "Genus")

mat_pca <- df_pca_norm %>%
  pivot_wider(id_cols = Sample, names_from = Genus,
              values_from = Relative_Abundance, values_fill = 0) %>%
  column_to_rownames("Sample")

# Transformation de Hellinger : prend la racine carrée des abondances relatives
mat_hellinger <- decostand(mat_pca, method = "hellinger")

# rda() sans variable environnementale = ACP standard (matrice de covariance)
pca_res <- rda(mat_hellinger)

# scaling = 1 : préserve les distances entre échantillons (recommandé en écologie)
pca_coords <- as.data.frame(scores(pca_res, display = "sites", scaling = 1))
colnames(pca_coords) <- c("PC1", "PC2")
pca_coords$Sample <- as.character(rownames(pca_coords))

# Variance expliquée par chaque axe (en %)
eig_vals <- pca_res$CA$eig
var_pc1  <- round(eig_vals[1] / sum(eig_vals) * 100, 1)
var_pc2  <- round(eig_vals[2] / sum(eig_vals) * 100, 1)

df_plot_pca <- pca_coords %>%
  inner_join(meta %>% select(Sample, Condition) %>% mutate(Sample = as.character(Sample)),
             by = "Sample") %>%
  mutate(Condition = factor(Condition, levels = chronologie_tout))

p_pca <- ggplot(df_plot_pca, aes(x = PC1, y = PC2, color = Condition, fill = Condition)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  # shape = 21 : cercle avec un contour noir (rend les points mieux visibles)
  geom_point(size = 4, shape = 21, color = "black", stroke = 0.6) +
  scale_color_manual(values = my_cond_colors) +
  scale_fill_manual(values = my_cond_colors) +
  labs(
    title    = "Microbial Community Structure",
    subtitle = "Hellinger-transformed relative abundances",
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

ggsave(file.path(out_dir, "PCA_microbiote.png"),
       plot = p_pca, width = 10, height = 7, dpi = 300, bg = "white")


# ==============================================================================
# BLOC 15 : CORRÉLATION ALPHA vs BÊTA (RÉPLICATS + CENTROÏDES)
# ==============================================================================
# Ce graphique teste et visualise la relation entre :
#   - La diversité INTRA-ÉCHANTILLON (Shannon) : sur l'axe Y
#   - L'HÉTÉROGÉNÉITÉ entre réplicats (Bray-Curtis moyen) : sur l'axe X
#
# Deux couches sont superposées :
#   - Petits points transparents : les réplicats individuels
#   - Grands points opaques : les centroïdes (moyennes exactes par stade)
# La droite de régression linéaire est calculée sur les réplicats individuels.
message("--- Corrélation Alpha vs Bêta ---")

label_cor_ins <- sprintf("Spearman's rho = %s\np-value %s", rho_sans, p_val_sans)

# On borne l'axe X à 1.0 maximum (Bray-Curtis ≤ 1 par définition)
max_bray    <- max(df_ind_sans$Beta_Ind, na.rm = TRUE)
max_x_limit <- ifelse(max_bray > 0.9, 1.0, max_bray + 0.05)

p_corr_ins <- ggplot() +
  # Droite de régression linéaire (sur les réplicats individuels)
  geom_smooth(data = df_ind_sans, aes(x = Beta_Ind, y = Shannon),
              method = "lm", color = "black", linetype = "dashed",
              alpha = 0.15, linewidth = 0.8, fullrange = FALSE) +
  # Réplicats individuels (petits points semi-transparents)
  geom_point(data = df_ind_sans,
             aes(x = Beta_Ind, y = Shannon, fill = Condition),
             shape = 21, size = 3, color = "black", stroke = 0.4, alpha = 0.4,
             position = position_jitter(width = 0.005, height = 0.005)) +
  # Centroïdes par stade (grands points opaques)
  geom_point(data = df_sans,
             aes(x = Beta_Mean, y = MeanShannon, fill = Condition),
             shape = 21, size = 6, color = "black", stroke = 1.2, alpha = 1) +
  # Étiquettes des centroïdes
  geom_text_repel(data = df_sans,
                  aes(x = Beta_Mean, y = MeanShannon, label = Condition, color = Condition),
                  fontface = "bold", size = 4,
                  box.padding = 0.8, point.padding = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = my_cond_colors) +
  scale_color_manual(values = my_cond_colors) +
  scale_x_continuous(limits = c(NA, max_x_limit), breaks = scales::pretty_breaks(n = 6)) +
  # Annotation du test de Spearman (coin bas-gauche)
  annotate("text", x = -Inf, y = -Inf,
           label = label_cor_ins,
           hjust = -0.1, vjust = -0.5,
           fontface = "bold", color = "#C0392B", size = 5) +
  labs(
    title    = "Relationship between inter-individual instability and microbial complexity",
    subtitle = "Large points = Stage exact mean | Small points = Individual replicates",
    x = "Intra-group heterogeneity (Mean Bray-Curtis distance)",
    y = "Individual Alpha diversity (Shannon Index)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position  = "none",
    plot.title       = element_text(face = "bold", size = 15),
    plot.subtitle    = element_text(color = "grey40", size = 12, margin = margin(b = 15)),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", linewidth = 1)
  )

ggsave(file.path(out_dir, "Correlation_Bray_vs_Shannon.png"),
       plot = p_corr_ins, width = 11, height = 7, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "Correlation_Bray_vs_Shannon.svg"),
       plot = p_corr_ins, width = 11, height = 7, bg = "white")
message("=> Graphique de corrélation Alpha/Bêta généré.")


# ==============================================================================
# BLOC 16 : RÉSOLUTION TAXONOMIQUE INTRA-GENRE (Barplots Horizontaux Fragmentés)
# ==============================================================================
# Pour les genres les plus abondants, on montre QUELLES ESPÈCES les composent
# et dans quelle proportion. Chaque genre est représenté dans sa propre facette.
# L'axe X = % des reads de ce genre attribués à chaque espèce (somme = 100% par genre).
# C'est un graphique "non-empilé" (barres non superposées) : chaque espèce
# est une barre indépendante, ce qui facilite la comparaison entre espèces.
message("--- Résolution Taxonomique Intra-Genre ---")

generer_resolution_fragmente <- function(stades_cibles, titre_stade, nom_fichier,
                                          top_n_genres = 8) {
  message(sprintf("  -> Stade(s) : %s", paste(stades_cibles, collapse = ", ")))

  df_stade <- df_final_16s %>% filter(Condition %in% stades_cibles)
  if (nrow(df_stade) == 0) return(NULL)

  # --- A. Abondance moyenne exacte par genre (pour l'étiquette de la facette) ---
  df_norm_gen  <- normaliser_100(df_stade, "Genus")
  n_samp_stade <- n_distinct(df_norm_gen$Sample)

  df_gen_agg <- df_norm_gen %>%
    group_by(Genus) %>%
    summarise(MeanAbund = sum(Relative_Abundance) / n_samp_stade, .groups = "drop") %>%
    filter(!Genus %in% c("Non Assigné", "Unassigned"))

  top_gen       <- df_gen_agg %>% slice_max(MeanAbund, n = top_n_genres) %>% pull(Genus)
  df_gen_labels <- df_gen_agg %>% filter(Genus %in% top_gen)

  # --- B. Proportions intra-genre par espèce ---
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
    # % calculé sur le total des reads DE CE GENRE (pas de l'échantillon entier)
    mutate(Pct_in_Genus = (Count / sum(Count)) * 100) %>%
    ungroup()

  # --- C. Simplification (Top 5 espèces par genre + "Other") ---
  df_res_clean <- df_intra %>%
    group_by(Genus) %>%
    mutate(
      Rank = case_when(
        Species_Label == "Unclassified (sp.)" ~ 999,  # Forcé en bas
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
    # L'étiquette de la facette inclut l'abondance globale du genre
    mutate(Genus_Label = sprintf("%s (%.2f%%)", Genus, MeanAbund))

  # --- D. Ordonnancement esthétique ---
  df_res_clean <- df_res_clean %>%
    mutate(
      Genus_Label = fct_reorder(Genus_Label, MeanAbund, .desc = TRUE),
      # fct_relevel force les catégories génériques en bas du graphique
      Species_Final = fct_relevel(
        fct_reorder(Species_Final, Pct_in_Genus),
        "Unclassified (sp.)", "Other identified species",
        after = 0
      )
    )

  pal_gen <- setNames(colorspace::qualitative_hcl(length(top_gen), palette = "Dark 3"), top_gen)

  # --- E. Graphique fragmenté (barres horizontales, 1 barre = 1 espèce) ---
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

  ggsave(file.path(out_dir, nom_fichier), plot = p_res, width = 13, height = 9, dpi = 300)
  message(sprintf("  => '%s' sauvegardé.", nom_fichier))
}

generer_resolution_fragmente(
  stades_cibles = cond_larves,
  titre_stade   = "Average Larval Microbiome",
  nom_fichier   = "Taxonomic_Resolution_Average_Larvae.png",
  top_n_genres  = 10
)

generer_resolution_fragmente(
  stades_cibles = c("Beetles"),
  titre_stade   = "Adult Microbiome (Beetles)",
  nom_fichier   = "Taxonomic_Resolution_Adults.png",
  top_n_genres  = 10
)


# ==============================================================================
# FIN DU SCRIPT
# ==============================================================================
message("\n=== ANALYSE DIVERSITÉ 16S TERMINÉE AVEC SUCCÈS ===")
message(sprintf("Résultats sauvegardés dans : %s", out_dir))
