# 16S Metagenomic Pipeline — *Tenebrio molitor*

Bioinformatic analysis of bacterial communities associated with the different developmental stages of *Tenebrio molitor* using full-length 16S and Shotgun sequencing.

---

## Objective

Process raw sequencing data to produce accurate bacterial taxon quantification, using a **hybrid de novo assembly approach** targeting the 16S rRNA gene.

Unlike standard ASV-based approaches, this pipeline reconstructs near-complete 16S sequences (scaffolds) through de novo assembly with MATAM, enabling far superior taxonomic resolution. Abundances are then quantified probabilistically by Salmon using the Expectation-Maximization (EM) algorithm, and each sequence is assigned a taxonomy down to species level via DADA2 and the SILVA reference database. Finally, a comprehensive diversity analysis is performed in R to characterise microbial communities across developmental stages.

The entire pipeline is designed to run on an **HPC computing cluster** via the **SLURM** job scheduler, with parallelisation through job arrays. The final statistical analysis step runs locally.

---

## Tools and versions

| Tool | Version | Invocation method | Role |
|---|---|---|---|
| FastQC | 0.11.7 | SLURM module | Quality control of raw reads (Phred scores, GC%, adapters) |
| MultiQC | 1.7 | SLURM module | Aggregation of FastQC reports into a single interactive HTML dashboard |
| Python | 3.7.1 | SLURM module | Read cleaning and streaming parsing |
| GCC | 4.8.4 / 8.1.0 | SLURM modules | C/C++ compiler (system dependencies) |
| Singularity | cluster | System binary | Isolation and execution of MATAM and SortMeRNA environments |
| SortMeRNA | 2.1b | Singularity image (Biocontainers) | Reference database indexing |
| USEARCH | 9.2.64 | SLURM module | Ultra-fast sequence dereplication |
| MATAM | 1.6.1 / 1.6.2* | Singularity images | De novo assembly targeting 16S rRNA |
| Salmon | 1.10.2 | Conda 23.3.1 | Quantification via pseudo-alignment and EM algorithm |
| R / DADA2 | — | Local R script | Taxonomic assignment (Bayesian classifier + exact identity) |
| R / ggplot2, vegan, etc. | — | Local R script | Alpha/beta diversity analysis and visualisation |
| SILVA | v138.2 | Local reference database | Taxonomic reference for 16S assignment |

---

## Project architecture

```
16_tenebrion/
├── bin/                        # SLURM submission scripts
│   ├── 16S_fastqc.slurm        # Step 1 — Initial quality control
│   ├── 16S_multiqc.slurm       # Step 2 — QC report aggregation
│   ├── pull_MATAM_sif.slurm    # Step 3 — Install MATAM via Singularity
│   ├── 16S_MATAM.slurm         # Step 4 — Hybrid production pipeline
│   ├── assign_taxo_MATAM.R     # Step 5 — Taxonomic assignment (DADA2/SILVA)
│   └── 16S_analysis.R          # Step 6 — Diversity analysis and visualisation
│
├── containers/                 # Singularity images (.sif)
├── data/
│   ├── raw/                    # Raw files (*_R1.fastq.gz, *_R2.fastq.gz)
│   └── 16S/
│       ├── MATAM/              # Compiled sequence table (Salmon output)
│       ├── SILVA/              # SILVA v138.2 reference databases
│       ├── metadata.tsv        # Sample metadata (stage → condition)
│       └── tax_info_sans_chimere.rds   # Taxonomic dictionary (Step 5 output)
├── log/                        # SLURM log files (.out / .err)
├── resources/                  # Conda environments
└── results/
    ├── qc/                     # FastQC and MultiQC reports
    ├── PRODUCTION_HYBRID/      # MATAM scaffolds + Salmon abundance table
    ├── 16S/
    │   └── alpha_beta/         # Diversity plots (Step 6 output)
    └── all_matam_salmon_qiime_like_table_counts_wSpecies.tsv  # Final table (Step 5 output)
```

---

## Running the pipeline

Scripts must be submitted in the following order from the project root:

```bash
sbatch bin/16S_fastqc.slurm
sbatch bin/16S_multiqc.slurm
sbatch bin/pull_MATAM_sif.slurm
sbatch bin/16S_MATAM.slurm
Rscript bin/assign_taxo_MATAM.R
Rscript bin/16S_analysis.R
```

Each step must be completed before submitting the next. Steps 1 and 4 run as SLURM Arrays. Step 6 runs locally (off-cluster).

---

## Script descriptions

---

### Step 1 — `16S_fastqc.slurm`: Initial quality control

**Objective:** Assess the intrinsic quality of raw sequencing data before any processing.

The script leverages **SLURM Arrays** to process each `*_R1.fastq.gz` file in parallel. It automatically detects the number of samples present in `data/raw/`. Temporary files are written to a dedicated folder on the cluster's `/storage/scratch` to avoid concurrent write conflicts between array tasks and to prevent saturating the main network filesystem I/O.

---

### Step 2 — `16S_multiqc.slurm`: Quality control aggregation

**Objective:** Centralise the quality assessment of the entire dataset into a single report.

MultiQC parses all FastQC reports generated in the previous step and produces a **single interactive HTML report**. This report allows visual identification at a glance of samples whose behaviour deviates from the norm (outliers), which may require special treatment before assembly.

---

### Step 3 — `pull_MATAM_sif.slurm`: Environment initialisation

**Objective:** Prepare the entire software infrastructure required by the hybrid assembly pipeline.

This script performs three successive operations:

**Building the MATAM Singularity image.** Conversion from a Docker image (Biocontainers) to a `.sif` file encapsulates all of MATAM's heavy dependencies (Python 2/3, SGA assembler) without polluting the cluster host environment.

**Downloading the SILVA SSURef NR95 database.** Retrieval uses a fallback mechanism: multiple URLs are tested in succession to guarantee the download even if a mirror is unavailable.

**Adaptive indexing.** Indexing the reference database requires the `indexdb_rna` binary from SortMeRNA. The script first checks whether this binary is available in the main MATAM image. If not, it dynamically downloads a dedicated SortMeRNA container to run this operation.

---

### Step 4 — `16S_MATAM.slurm`: Hybrid production pipeline

**Objective:** Transform validated raw sequences into a taxonomic abundance table. Runs as a SLURM Array, one job per sample.

The script chains four sub-steps:

**Strict cleaning (Python).** An on-the-fly Python script reads sequences as a continuous stream (*streaming*), without loading the entire file into memory — a strategy suited to large volumes. The first 15 nucleotides of each read are trimmed (a region often noisy due to primers), and sequences that are too short are discarded.

**Dereplication (USEARCH).** Strictly identical reads are merged into unique entities. This step acts as a strong information compression and considerably reduces the search space and time complexity for the downstream assembler.

**De novo assembly (MATAM).** Unlike standard ASV approaches, MATAM uses short reads and the SILVA database to reconstruct near-complete 16S rRNA sequences (scaffolds). This reconstruction achieves much finer taxonomic resolution.

**Probabilistic quantification (Salmon).** Once the 16S sequence catalogue is assembled, Salmon indexes it and virtually aligns the original reads against it. The **Expectation-Maximization (EM)** algorithm resolves the ambiguity of reads mapping with equivalent probability to multiple closely related taxa, producing a robust final abundance table.

---

### Step 5 — `assign_taxo_MATAM.R`: Taxonomic assignment

**Objective:** Assign a complete taxonomy (down to species level) to each sequence assembled by MATAM, and produce an abundance table in standard QIIME2 format.

The script takes as input the compiled sequence table produced by Salmon (`all_matam_sequences_compiled.tsv`) and chains five sub-steps:

**Sequence deduplication.** Only unique sequences are extracted before assignment, to avoid classifying the same sequence multiple times when it appears in several samples. The full table is retained for the final join.

**Genus assignment (Bayesian classifier, `assignTaxonomy`).** The DADA2 naïve Bayesian classifier (Wang et al. 2007) is applied against the SILVA v138.2 database with a minimum bootstrap threshold of 80 %. Ranks below this threshold remain `NA` rather than forcing an uncertain assignment. DADA2's internal parallelisation (OpenMP) is enabled.

**Species assignment (`addSpecies`, manually parallelised).** `addSpecies()` searches for 100 % identity matches in the SILVA database. As this operation is inherently sequential, it is manually parallelised in chunks of 2,000 sequences via `mclapply()`, reserving 3 cores for system stability.

**Minimal biological filtering.** Only sequences with no assignment at the Kingdom rank (`Kingdom = NA`) are removed — these most likely correspond to chimeras or artefacts. Unresolved lower ranks (`NA`) are retained as they represent valid biological information.

**Dual export.** The script produces two files:
- `tax_info_sans_chimere.rds`: sequence → complete taxonomy dictionary, in RDS format (for PICRUSt2).
- `all_matam_salmon_qiime_like_table_counts_wSpecies.tsv`: pivoted abundance table (rows = taxa in `d__`;`p__`;...;`s__` format, columns = samples), compliant with the QIIME2/BIOM convention.

---

### Step 6 — `16S_analysis.R`: Diversity analysis and visualisation

**Objective:** Characterise the gut microbial communities of *T. molitor* across developmental stages (larvae → adults) and produce the full set of analysis figures.

The script takes as input the QIIME2-like table produced at Step 5 and a metadata file (`data/16S/metadata.tsv`), and runs through 16 analysis blocks:

**Data loading and harmonisation (Blocks 1–2).** The count file is read robustly (handling the `#OTU ID` header specific to QIIME2 format). Taxonomic ranks are extracted by regex from SILVA strings (`d__`, `p__`, ..., `s__`). A taxonomic revision dictionary (`harmoniser_taxonomie()`) corrects genera reclassified in SILVA 138 (splits of *Lactobacillus*, *Bacillus*, *Mycobacterium*, *Burkholderia*, etc.) to ensure cross-sample consistency. Plastids (chloroplasts, mitochondria) are filtered out. Binomial nomenclature is standardised.

**Alpha and beta diversity (Blocks 3–4).** The Shannon index is computed per sample on 100 %-normalised counts. Bray-Curtis dissimilarity is calculated via the `vegan` library (sample × sample matrix) to measure community divergence between replicates within each stage. Both metrics are displayed jointly on time-course curves following the developmental chronology.

**Taxonomic compositions (Blocks 5–8).** Stacked barplots at Genus and Species level are generated for larval stages, adult stages, and the full dataset (insect + environmental substrates). Rare taxa are grouped into an "Others" category to keep figures readable.

**Similarity matrices (Block 9).** Bray-Curtis heatmaps are produced for all replicates (insect + substrate) and for the insect-only subset, enabling visual identification of stage-level clustering.

**Master figures (Blocks 10–11).** Alpha/beta curves and barplots are combined via `patchwork` into multi-panel publication-ready figures, for both the insect-only and the full dataset including substrates.

**Core microbiome (Blocks 12–15).** For larvae and adults, the core microbiome is characterised at two taxonomic ranks (Genus and Species) using two complementary approaches: mean abundance barplots per stage, and Prevalence/Abundance quadrant plots (ecological classification into *Strict core*, *Satellite*, *Transient* and *Background noise*). Donut charts summarise adult core composition.

**Alpha/beta correlation (Block 15 bis).** A Spearman test evaluates the relationship between intra-group heterogeneity (mean Bray-Curtis distance) and individual alpha diversity (Shannon), with stage centroids and individual replicates represented.

**Intra-genus taxonomic resolution (Block 16).** For the most abundant genera in larval and adult stages, fragmented horizontal barplots show the intra-genus breakdown into species, enabling assessment of the taxonomic resolution achieved by MATAM/SILVA.

**Output figures** (saved to `results/16S/alpha_beta/`):

| File | Content |
|---|---|
| `Valeurs_Shannon_BrayCurtis.tsv` | Raw numerical table of alpha/beta indices |
| `shannon_vs_braycurtis_.png` | Alpha/beta curves (insect only) |
| `shannon_vs_braycurtis_substrats.png` | Alpha/beta curves (insect + substrates) |
| `Barplot_Global_16S.png` | Taxonomic composition per replicate (all stages) |
| `Matrice_BrayCurtis_Tous_Replicats.png` | Bray-Curtis heatmap (all replicates) |
| `Matrice_BrayCurtis_Insecte_Seul.png` | Bray-Curtis heatmap (insect only) |
| `Master_Figure.png` / `_substrats.png` | Multi-panel master figures |
| `Core_Larval_Microbiome_Means.png` | Mean larval core microbiome (Genus) |
| `Core_Larval_Microbiome_Means_Species.png` | Mean larval core microbiome (Species) |
| `Core_Microbiome_Quadrants.png` | Prevalence/Abundance quadrants, larvae (Genus) |
| `Core_Microbiome_Quadrants_Species.png` | Prevalence/Abundance quadrants, larvae (Species) |
| `Core_genome_heterogeneous.png` | Intra-stage fidelity matrices |
| `Core_Adult_Microbiome_Means.png` | Mean adult core microbiome (Genus) |
| `Core_Adult_Microbiome_Donut.png` | Adult donut chart (Genus) |
| `Core_Adult_Microbiome_Means_Species.png` | Mean adult core microbiome (Species) |
| `Core_Adult_Microbiome_Donut_Species.png` | Adult donut chart (Species) |
| `Core_Adult_Microbiome_Quadrants.png` | Prevalence/Abundance quadrants, adults (Genus) |
| `Core_Adult_Microbiome_Quadrants_Species.png` | Prevalence/Abundance quadrants, adults (Species) |
| `PCA_microbiote.png` | PCA of community structures |
| `Correlation_Bray_vs_Shannon.png` / `.svg` | Alpha vs beta correlation (Spearman) |
| `Taxonomic_Resolution_Average_Larvae.png` | Intra-genus resolution, larvae (top 10 genera) |
| `Taxonomic_Resolution_Adults.png` | Intra-genus resolution, adults (top 10 genera) |

---

## Pipeline diagram

```
Raw data (FASTQ)
        |
        v
[1] FastQC              Per-sample quality control (SLURM Array)
        |
        v
[2] MultiQC             Aggregated QC report (interactive HTML)
        |
        v
[3] pull_MATAM_sif      MATAM Singularity image + indexed SILVA database
        |
        v
[4a] Cleaning           Primer clipping + length filtering (Python streaming)
        |
        v
[4b] Dereplication      Compression of identical reads (USEARCH)
        |
        v
[4c] Assembly           Reconstruction of near-complete 16S scaffolds (MATAM)
        |
        v
[4d] Quantification     Per-sequence abundance table (Salmon + EM)
        |
        v
[5a] Genus assignment   DADA2 Bayesian classifier × SILVA v138.2 (minBoot=80)
        |
        v
[5b] Species assignment addSpecies() 100% identity × SILVA (mclapply, chunks 2000)
        |
        v
[5c] Export             .rds dictionary + QIIME2 .tsv table
        |
        v
[6a] Harmonisation      Taxonomic corrections + plastid filtering
        |
        v
[6b] α/β Diversity      Shannon per sample + Bray-Curtis between replicates
        |
        v
[6c] Visualisation      Barplots, heatmaps, quadrants, PCA, correlations (22 figures)
```

---

## Prerequisites

- HPC cluster with **SLURM** job scheduler
- **Singularity / Apptainer** available as a system binary
- Modules available on the cluster: `fastqc/0.11.7`, `MultiQC/1.7`, `python/3.7.1`, `gcc/4.8.4`, `gcc/8.1.0`, `usearch/9.2.64`
- **Conda** environment (version 23.3.1) including Salmon 1.10.2
- **R** with packages: `dada2`, `dplyr`, `tidyr`, `readr`, `parallel` (Step 5) and `ggplot2`, `vegan`, `tibble`, `stringr`, `forcats`, `colorspace`, `patchwork`, `ggrepel`, `scales` (Step 6)
- SILVA v138.2 reference databases (`toGenus_trainset.fa.gz` and `assignSpecies.fa.gz`) in `data/16S/SILVA/`
- Metadata file (`data/16S/metadata.tsv`) with columns `sample-id` and `condition`
- Internet access from compute nodes for downloading SILVA and Singularity images

---

## Author

Thomas BOUTET — Ténébrion Project, 16S metagenomic analysis of *Tenebrio molitor*
