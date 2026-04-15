# Pipeline Métagénomique 16S — *Tenebrio molitor*

Analyse bioinformatique des communautés bactériennes associées aux différents stades de développement de *Tenebrio molitor* par séquençage du 16S full-length et Shotgun.

---

## Objectif

Traiter des données de séquençage brutes pour aboutir à une quantification précise des taxons bactériens, en utilisant une approche **hybride d'assemblage de novo** ciblée sur le gène de l'ARNr 16S.

Contrairement aux approches ASV standards, ce pipeline reconstruit des séquences 16S quasi-complètes (scaffolds) par assemblage de novo via MATAM, permettant une résolution taxonomique bien supérieure. Les abondances sont ensuite quantifiées de manière probabiliste par Salmon via l'algorithme Expectation-Maximization (EM).

L'ensemble du pipeline est conçu pour s'exécuter sur un **cluster de calcul HPC** via le gestionnaire de tâches **SLURM**, avec parallélisation par job arrays.

---

## Outils et versions

| Outil | Version | Méthode d'appel | Rôle |
|---|---|---|---|
| FastQC | 0.11.7 | Module SLURM | Contrôle qualité des lectures brutes (scores Phred, GC%, adaptateurs) |
| MultiQC | 1.7 | Module SLURM | Agrégation des rapports FastQC en un tableau de bord HTML unique |
| Python | 3.7.1 | Module SLURM | Nettoyage et parsing des lectures en flux continu |
| GCC | 4.8.4 / 8.1.0 | Modules SLURM | Compilateur C/C++ (dépendances système) |
| Singularity | cluster | Binaire système | Isolation et exécution des environnements MATAM et SortMeRNA |
| SortMeRNA | 2.1b | Image Singularity (Biocontainers) | Indexation de la base de données de référence |
| USEARCH | 9.2.64 | Module SLURM | Déréplication ultrarapide des séquences |
| MATAM | 1.6.1 / 1.6.2* | Images Singularity | Assemblage de novo ciblé sur l'ARNr 16S |
| Salmon | 1.10.2 | Conda 23.3.1 | Quantification par pseudo-alignement et algorithme EM |

*Note : le script de téléchargement pointe vers l'image MATAM 1.6.2, tandis que le script de production appelle l'image 1.6.1.*

---

## Architecture du projet

```
16_tenebrion/
├── bin/                        # Scripts de soumission SLURM (cœur du pipeline)
│   ├── 16S_fastqc.slurm        # Étape 1 — Contrôle qualité initial
│   ├── 16S_multiqc.slurm       # Étape 2 — Agrégation des rapports QC
│   ├── pull_MATAM_sif.slurm    # Étape 3 — Initialisation des environnements
│   └── 16S_MATAM.slurm         # Étape 4 — Pipeline hybride de production
├── containers/                 # Images Singularity (.sif)
├── data/
│   └── raw/                    # Fichiers bruts (*_R1.fastq.gz, *_R2.fastq.gz)
├── log/                        # Fichiers de log SLURM (.out / .err)
├── resources/                  # Base de données SILVA + environnements Conda
└── results/
    ├── qc/                     # Rapports FastQC et MultiQC
    └── PRODUCTION_HYBRID/      # Scaffolds MATAM + table d'abondance Salmon
```

---

## Lancer le pipeline

Les scripts doivent être soumis dans l'ordre suivant depuis la racine du projet :

```bash
sbatch bin/16S_fastqc.slurm
sbatch bin/16S_multiqc.slurm
sbatch bin/pull_MATAM_sif.slurm
sbatch bin/16S_MATAM.slurm
```

Chaque étape doit être complétée avant de soumettre la suivante. Les étapes 1 et 4 s'exécutent en mode SLURM Array.

---

## Description des scripts

---

### Étape 1 — `16S_fastqc.slurm` : Contrôle qualité initial

**Objectif :** Évaluer la qualité intrinsèque des données de séquençage brutes avant tout traitement.

Le script exploite les **SLURM Arrays** pour traiter chaque fichier `*_R1.fastq.gz` en parallèle. Il détecte automatiquement le nombre d'échantillons présents dans `data/raw/`. Les fichiers temporaires sont écrits dans un dossier dédié sur le `/storage/scratch` du cluster afin d'éviter les conflits d'écriture concurrents entre les tâches de l'array et de ne pas saturer les I/O du système de fichiers réseau principal.

---

### Étape 2 — `16S_multiqc.slurm` : Agrégation du contrôle qualité

**Objectif :** Centraliser l'interprétation qualité de l'ensemble du jeu de données en un seul rapport.

MultiQC parse l'ensemble des rapports FastQC générés à l'étape précédente et produit un **rapport HTML interactif unique**. Ce rapport permet d'identifier visuellement en un seul coup d'œil les échantillons dont le comportement s'écarte de la norme (outliers), qui pourraient nécessiter un traitement particulier avant l'assemblage.

---

### Étape 3 — `pull_MATAM_sif.slurm` : Initialisation de l'environnement

**Objectif :** Préparer toute l'infrastructure logicielle requise par le pipeline d'assemblage hybride.

Ce script réalise trois opérations successives :

**Construction de l'image Singularity MATAM.** La conversion depuis une image Docker (Biocontainers) vers un fichier `.sif` encapsule l'intégralité des dépendances lourdes de MATAM (Python 2/3, assembleur SGA) sans polluer l'environnement hôte du cluster.

**Téléchargement de la base SILVA SSURef NR95.** La récupération s'effectue via un mécanisme de fallback : plusieurs URLs sont testées successivement pour garantir le téléchargement même en cas d'indisponibilité d'un miroir.

**Indexation adaptative.** L'indexation de la base de référence nécessite le binaire `indexdb_rna` de SortMeRNA. Le script vérifie d'abord si ce binaire est disponible dans l'image MATAM principale. Dans le cas contraire, il télécharge dynamiquement un conteneur SortMeRNA dédié pour exécuter cette opération.

---

### Étape 4 — `16S_MATAM.slurm` : Pipeline hybride de production

**Objectif :** Transformer les séquences brutes validées en une table d'abondance taxonomique. S'exécute en mode SLURM Array, un job par échantillon.

Le script enchaîne quatre sous-étapes :

**Nettoyage strict (Python).** Un script Python implémenté à la volée lit les séquences en flux continu (*streaming*), sans chargement intégral en mémoire — stratégie adaptée aux gros volumes. Les 15 premiers nucléotides de chaque lecture sont retirés (région souvent bruitée par les amorces), et les séquences trop courtes sont éliminées.

**Déréplication (USEARCH).** Les lectures strictement identiques sont fusionnées en entités uniques. Cette étape agit comme une forte compression de l'information et réduit considérablement l'espace de recherche et la complexité temporelle pour l'assembleur qui suit.

**Assemblage de novo (MATAM).** Contrairement aux approches ASV standards, MATAM utilise les short-reads et la base SILVA pour reconstruire des séquences d'ARNr 16S quasi-complètes (scaffolds). Cette reconstruction permet d'atteindre une résolution taxonomique bien plus fine, au prix d'un calcul plus intensif.

**Quantification probabiliste (Salmon).** Une fois le catalogue de séquences 16S assemblé, Salmon l'indexe et aligne virtuellement les lectures initiales dessus. L'algorithme **Expectation-Maximization (EM)** résout l'ambiguïté des reads s'alignant avec une probabilité équivalente sur plusieurs taxons proches, produisant une table d'abondance finale robuste.

---

## Schéma du pipeline

```
Données brutes (FASTQ)
        |
        v
[1] FastQC              Contrôle qualité par échantillon (SLURM Array)
        |
        v
[2] MultiQC             Rapport QC agrégé (HTML interactif)
        |
        v
[3] pull_MATAM_sif      Image Singularity MATAM + base SILVA indexée
        |
        v
[4a] Nettoyage          Clipping amorces + filtrage longueur (Python streaming)
        |
        v
[4b] Déréplication      Compression des lectures identiques (USEARCH)
        |
        v
[4c] Assemblage         Reconstruction de scaffolds 16S quasi-complets (MATAM)
        |
        v
[4d] Quantification     Table d'abondance taxonomique (Salmon + EM)
```

---

## Prérequis

- Cluster HPC avec gestionnaire de tâches **SLURM**
- **Singularity / Apptainer** disponible en tant que binaire système
- Modules disponibles sur le cluster : `fastqc/0.11.7`, `MultiQC/1.7`, `python/3.7.1`, `gcc/4.8.4`, `gcc/8.1.0`, `usearch/9.2.64`
- Environnement **Conda** (version 23.3.1) incluant Salmon 1.10.2
- Accès internet depuis les noeuds de calcul pour le téléchargement de SILVA et des images Singularity

---

## Auteur

Thomas BOUTET — Projet Ténébrion, analyse métagénomique 16S de *Tenebrio molitor*
