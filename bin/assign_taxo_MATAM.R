# ==============================================================================
# SCRIPT DE PRODUCTION : Assignation Taxonomique DADA2 (Double Export)
# OBJECTIF : Générer le dico PICRUSt2 (.rds) ET la table de diversité (.tsv)
# VERSION  : Optimisée multi-cœurs (mclapply) & Formatage QIIME2 natif
#
# Flux général du script :
#   Fichier TSV (MATAM)
#       └─> Extraction des séquences uniques
#               └─> assignTaxonomy()  →  rang jusqu'au Genre (classifieur Bayésien)
#               └─> addSpecies()      →  rang Espèce (identité 100 % SILVA)
#                       └─> Export 1 : dictionnaire RDS  (pour PICRUSt2)
#                       └─> Export 2 : table TSV format QIIME2 (pour l'analyse de diversité)
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. CHARGEMENT DES BIBLIOTHÈQUES
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)   
  library(tidyr)   
  library(dada2)    
  library(readr)   
  library(parallel)
})

# ------------------------------------------------------------------------------
# 1. PARAMÈTRES ET CHEMINS
# ------------------------------------------------------------------------------
setwd("/home/thomas/Tenebrion/")

fichier_entree  <- "data/16S/MATAM/all_matam_sequences_compiled.tsv"

# SILVA v138.2 est la référence standard pour le 16S.
fichier_genus   <- "data/16S/SILVA/silva_nr99_v138.2_toGenus_trainset.fa.gz"
fichier_species <- "data/16S/SILVA/silva_v138.2_assignSpecies.fa.gz"

fichier_sortie_rds   <- "data/16S/tax_info_sans_chimere.rds"
fichier_sortie_qiime <- "results/all_matam_salmon_qiime_like_table_counts_wSpecies.tsv"

# ------------------------------------------------------------------------------
# 2. EXTRACTION DES SÉQUENCES UNIQUES
# ------------------------------------------------------------------------------
message("[1/5] Lecture des données brutes...")


df_complet <- read_tsv(fichier_entree, show_col_types = FALSE)

# On déduplique les séquences : DADA2 n'a besoin de classifier chaque
# séquence qu'une seule fois, même si elle apparaît dans plusieurs échantillons.
# Cela réduit considérablement la charge de calcul.
seqs_uniques <- unique(df_complet$Sequence)
n_total <- length(seqs_uniques)

# df_complet est conservé intact pour la jointure ultérieure (étape 5),
# afin de pouvoir relier chaque séquence à ses comptes par échantillon.
df <- df_complet

# ------------------------------------------------------------------------------
# 3. ASSIGNATION DADA2 (PARALLÉLISÉE)
# ------------------------------------------------------------------------------
message("[2/5] Lancement de l'assignation taxonomique...")

# 3A : Assignation jusqu'au Genre (Bayésien) 
# assignTaxonomy() utilise le classifieur Bayésien naïf de Wang et al. (2007),
# réimplémenté dans DADA2.
# - minBoot = 80 : seuil de bootstrap minimal pour valider un rang taxonomique.
#   En dessous de 80, le rang est laissé à NA (recommandation SILVA/DADA2).
# - multithread = TRUE : utilise tous les cœurs disponibles via la
#   parallélisation interne de DADA2 (OpenMP).
message("      -> A : Genre (Bayésien, multithread natif)...")
tax_info <- assignTaxonomy(seqs_uniques, fichier_genus, minBoot = 80, multithread = TRUE)

# 3B : Ajout de l'Espèce 
# addSpecies() compare chaque séquence à la base d'espèces SILVA par alignement
# exact (100 % d'identité)
message("      -> B : Espèce (Parallélisation forcée par mclapply)...")
chunk_size <- 2000
indices <- seq(1, n_total, by = chunk_size)
n_cores <- max(1, detectCores() - 3)
message(sprintf("         Utilisation de %d cœurs pour le traitement parallèle.", n_cores))

tax_list <- mclapply(seq_along(indices), function(i) {
  start <- indices[i]
  end <- min(start + chunk_size - 1, n_total)  # Dernier indice du paquet (sans déborder)

  # addSpecies() : cherche une correspondance à 100 % d'identité dans SILVA.
  chunk_res <- addSpecies(tax_info[start:end, , drop=FALSE], fichier_species)
  return(chunk_res)
}, mc.cores = n_cores)

tax_final <- do.call(rbind, tax_list)

# ------------------------------------------------------------------------------
# 4. EXPORT N°1 : LE DICTIONNAIRE (.rds)
# ------------------------------------------------------------------------------
# Ce fichier RDS sera chargé directement par PICRUSt2 pour la prédiction
# fonctionnelle (métagénomique implicite). Il contient la correspondance
# séquence → taxonomie complète.
message("[3/5] Sauvegarde du dictionnaire pour PICRUSt2...")

# Création du répertoire de sortie si inexistant
if (!dir.exists(dirname(fichier_sortie_rds))) dir.create(dirname(fichier_sortie_rds), recursive = TRUE)
saveRDS(tax_final, fichier_sortie_rds)

# ------------------------------------------------------------------------------
# 5. EXPORT N°2 : LA TABLE FORMAT QIIME2 (PRÉSERVATION BIOLOGIQUE)
# ------------------------------------------------------------------------------
message("[4/5] Formatage de la taxonomie façon QIIME2...")

# 5A : Conversion de la matrice DADA2 en data.frame
tax_df <- as.data.frame(tax_final)
tax_df$Sequence <- rownames(tax_final)
rownames(tax_df) <- NULL  # Suppression des rownames devenus redondants

# 5B : Jointure avec les données brutes
df_joint <- df %>% left_join(tax_df, by = "Sequence")

# 5C : Filtrage biologique minimal

df_filtered <- df_joint %>% filter(!is.na(Kingdom))           # Les séquences sans Règne (Kingdom = NA) donc sans taxonomie sont supprimées

# 5D : Construction de la chaîne taxonomique au format QIIME2
# Convention QIIME2 / Silva : chaque rang est préfixé par son initiale et "__"
df_formatted <- df_filtered %>%
  mutate(
    d = ifelse(!is.na(Kingdom), paste0("d__", Kingdom), "__"),  # Domain/Kingdom
    p = ifelse(!is.na(Phylum),  paste0("p__", Phylum),  "__"),  # Phylum
    c = ifelse(!is.na(Class),   paste0("c__", Class),   "__"),  # Class
    o = ifelse(!is.na(Order),   paste0("o__", Order),   "__"),  # Order
    f = ifelse(!is.na(Family),  paste0("f__", Family),  "__"),  # Family
    g = ifelse(!is.na(Genus),   paste0("g__", Genus),   "__"),  # Genus
    s = ifelse(!is.na(Species), paste0("s__", Species), "__"),  # Species

    # Concaténation finale séparée par ";" → chaîne taxonomique complète
    Taxonomy = paste(d, p, c, o, f, g, s, sep = ";")
  ) %>%
  # Nettoyage : suppression des colonnes intermédiaires d, p, c, o, f, g, s
  select(-d, -p, -c, -o, -f, -g, -s)

# --- 5E : Pivotement en table d'abondance ---
message("[5/5] Création de la matrice et sauvegarde TSV...")

df_qiime <- df_formatted %>%
  # Agrégation : somme des lectures par combinaison (Taxonomie × Échantillon).
  # Plusieurs séquences différentes peuvent partager la même assignation
  # taxonomique → leurs reads sont additionnés.
  group_by(Taxonomy, Sample) %>%
  summarise(Count = sum(NumReads, na.rm = TRUE), .groups = "drop") %>%

  # Pivotement : on passe de la forme longue (1 ligne = 1 taxon×échantillon)
  # à la forme large (1 ligne = 1 taxon, 1 colonne par échantillon).
  # Les combinaisons absentes sont remplies par 0 (taxon absent de l'échantillon).
  pivot_wider(names_from = Sample, values_from = Count, values_fill = list(Count = 0)) %>%

  # Renommage de la colonne Taxonomy en "#OTU ID" : convention QIIME2 /
  # format BIOM, attendu par la plupart des outils d'analyse de diversité.
  rename(`#OTU ID` = Taxonomy)

# Création du répertoire de sortie si nécessaire, puis écriture du TSV
if (!dir.exists(dirname(fichier_sortie_qiime))) dir.create(dirname(fichier_sortie_qiime), recursive = TRUE)
write_tsv(df_qiime, file = fichier_sortie_qiime)

message("\n=== ASSIGNATION DE PRODUCTION TERMINÉE AVEC SUCCÈS ===")
message("1. Fichier PICRUSt2 prêt : ", fichier_sortie_rds)
message("2. Fichier Diversité prêt : ", fichier_sortie_qiime)