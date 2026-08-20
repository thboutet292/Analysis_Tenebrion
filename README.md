## Table of Contents

- [16S rRNA Pipeline](#16s-metagenomic-pipeline--tenebrio-molitor)
- [Shotgun Pipeline](#shotgun-metagenomic-pipeline--tenebrio-molitor)

---

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

**Strict cleaning (Python).** An on-the-fly Python script reads sequences as a continuous stream (*streaming*), without loading the entire file into memory — a strategy suited to large volumes. The first 15 nucleotides of each read are trimmed (a region often noisy due to primers), and sequences that are too short are discarded (>50 nt).

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
| `shannon_vs_braycurtis.png` / `.svg` | Alpha/beta curves (insect only) |
| `Barplot_Global_16S.png` / `.svg` | Taxonomic composition per replicate (insect only) |
| `Matrice_BrayCurtis_Tous_Replicats.png` | Bray-Curtis heatmap (all replicates) |
| `Matrice_BrayCurtis_Insecte_Seul.png` | Bray-Curtis heatmap (insect only) |
| `Master_Figure_substrats.png` / `.svg` | Multi-panel master figure (insect + substrates) |
| `Master_Figure.png` / `.svg` | Multi-panel master figure (insect only) |
| `Larval_Microbiome_Means.png` | Mean larval core microbiome (Genus) |
| `Larval_Microbiome_Means_Species.png` | Mean larval core microbiome (Species) |
| `Larval_Microbiome_Quadrants.png` | Prevalence/Abundance quadrants, larvae (Genus) |
| `Larval_Microbiome_Quadrants_Species.png` | Prevalence/Abundance quadrants, larvae (Species) |
| `Genome_heterogeneous.png` | Intra-stage fidelity matrices |
| `Adult_Microbiome_Means.png` | Mean adult core microbiome (Genus) |
| `Adult_Microbiome_Donut.png` | Adult donut chart (Genus) |
| `Adult_Microbiome_Means_Species.png` | Mean adult core microbiome (Species) |
| `Adult_Microbiome_Donut_Species.png` | Adult donut chart (Species) |
| `Adult_Microbiome_Quadrants.png` | Prevalence/Abundance quadrants, adults (Genus) |
| `Adult_Microbiome_Quadrants_Species.png` | Prevalence/Abundance quadrants, adults (Species) |
| `PCA_microbiote_insecte_seul.png` | PCA of community structures (insect only) |
| `PCA_microbiote_insecte_seul_ellipses.png` | PCA with 95% confidence ellipses (insect only) |
| `Correlation_Bray_vs_Shannon.png` / `.svg` | Alpha vs beta correlation (Spearman) |
| `Taxonomic_Resolution_Average_Larvae.png` | Intra-genus resolution, all larval stages pooled |
| `Taxonomic_Resolution_Microlarvae.png` | Intra-genus resolution, Microlarvae (7 mg) |
| `Taxonomic_Resolution_Larvae_S1.png` | Intra-genus resolution, Larvae S1 (14 mg) |
| `Taxonomic_Resolution_Larvae_S2.png` | Intra-genus resolution, Larvae S2 (40 mg) |
| `Taxonomic_Resolution_Larvae_S3.png` | Intra-genus resolution, Larvae S3 (65 mg) |
| `Taxonomic_Resolution_Larvae_S4.png` | Intra-genus resolution, Larvae S4 (100 mg) |
| `Taxonomic_Resolution_Adults.png` | Intra-genus resolution, adults |
| `Dendrogramme_Global_Epure.png` | Hierarchical clustering dendrogram (all samples) |
| `Dendrogramme_Insecte_Epure.png` | Hierarchical clustering dendrogram (insect only) |

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
[6c] Visualisation      Barplots, heatmaps, quadrants, PCA, correlations
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

# Shotgun Metagenomic Pipeline — *Tenebrio molitor*

Bioinformatic analysis of bacterial communities associated with the different developmental stages of *Tenebrio molitor* using Shotgun metagenomic sequencing.

---

## Objective

Process raw shotgun sequencing data to characterise bacterial communities through two complementary strategies:

- A **read-based approach** (Kraken2) for rapid, whole-community taxonomic profiling directly from cleaned reads, enabling fast comparisons across developmental stages.
- An **assembly-based approach** (assembly → binning → QC → taxonomy → functional annotation) to reconstruct Metagenome-Assembled Genomes (MAGs) and access functional gene content, complementing the taxonomic resolution obtained from the 16S pipeline.

Before any downstream analysis, host-derived reads (*Tenebrio molitor* genome) are systematically removed to avoid contamination of the microbial signal.

The entire pipeline is designed to run on an **HPC computing cluster** via the **SLURM** job scheduler, with parallelisation through job arrays. Raw and intermediate data are stored on an **S3 remote** (`s3_uca:lmge-tenebrion`), accessed via `rclone`.

---

## Tools and versions

| Tool | Version | Invocation method | Role |
|---|---|---|---|
| rclone | 1.55.1 | SLURM module | Data transfer to/from S3 storage |
| Bowtie2 | 2.3.4.3 | SLURM module | Host read alignment (decontamination) |
| samtools | 1.16.1 | SLURM module | BAM filtering, sorting, FASTQ extraction |
| Apptainer/Singularity | cluster | System binary | Execution of the Kraken2/Bracken, CoverM, binning, CheckM2, GTDB-Tk and DRAM containers |
| Kraken2 | *(per container image)* | Apptainer image | K-mer-based taxonomic assignment of paired reads |
| Bracken | *(per container image)* | Apptainer image | Bayesian re-estimation of species-level abundances from Kraken2 reports |
| NCBI Datasets CLI | v2 (linux-amd64) | Downloaded binary | Genome retrieval for the custom insect-focused Kraken2 addon database |
| MEGAHIT | 1.2.9 | SLURM module | De novo metagenomic assembly, run independently per sample (`--presets meta-sensitive`) |
| MMseqs2 | 13-45111 | SLURM module | Contig catalog construction via linear-time clustering (dereplication) |
| CoverM | 0.7.0 | Apptainer image | Read mapping (minimap2-sr) and coverage/abundance quantification against the non-redundant catalog |
| MetaBAT2 | 2.15 | Apptainer image | Genome binning via TNF composition + differential coverage (`jgi_summarize_bam_contig_depths`) |
| SemiBin2 | 2.1.0 | Apptainer image | Genome binning via self-supervised deep learning (`single_easy_bin`) |
| CONCOCT | 1.1.0 | Apptainer image | Genome binning via Gaussian mixture clustering on chunked contigs |
| MaxBin2 | 2.2.7 | Apptainer image | Genome binning via Expectation-Maximization on tetranucleotide frequency + abundance (`run_MaxBin.pl`) |
| CheckM2 | 1.0.2 | Apptainer image | MAG quality assessment (completeness/contamination) via Machine Learning (Diamond BlastP) |
| GTDB-Tk | 2.3.2 (DB release R214) | Apptainer image | Taxonomic assignment of MAGs (marker gene placement + FastANI) |
| DRAM | 1.4.6 | Apptainer image | Functional annotation of MAGs (KEGG, dbCAN/CAZymes, MEROPS, Pfam, VOG) via `annotate` + `distill` |

---

## Project architecture

```
shotgun_tenebrion/
├── bin/                                     # SLURM submission scripts
│   ├── host_decontamination_bowtie2.slurm   # Step 1 — Host read removal
│   ├── build_kraken2_db.slurm               # Step 2a — Download the Kraken2 PlusPFP database
│   ├── build_kraken_custom.slurm            # Step 2b — Build the custom insect addon database
│   ├── shotgun_kraken2_bracken.slurm        # Step 2c — Kraken2/Bracken profiling (PlusPFP only)
│   ├── shotgun_kraken2_bracken_addon.slurm  # Step 2d — Kraken2/Bracken profiling (PlusPFP + Addon, dual)
│   ├── shotgun_megahit_full.slurm           # Step 3 — Per-sample de novo assembly (MEGAHIT)
│   ├── pull_coverm.slurm                    # Step 4a — Install the CoverM Apptainer image
│   ├── shotgun_coverm.slurm                 # Step 4b — Non-redundant catalog + coverage profiling
│   ├── pull_binners.slurm                   # Step 5a — Install the MetaBAT2 + CONCOCT + SemiBin2 Apptainer images
│   ├── shotgun_metabat2.slurm               # Step 5b — Binning with MetaBAT2
│   ├── shotgun_semibin2.slurm               # Step 5c — Binning with SemiBin2
│   ├── shotgun_concoct.slurm                # Step 5d — Binning with CONCOCT
│   ├── shotgun_maxbin2.slurm                # Step 5e — Binning with MaxBin2
│   ├── pull_checkm2.slurm                   # Step 6a — Install the CheckM2 Apptainer image + ML database
│   ├── build_checkm2_db.slurm               # Step 6b — (one-off) Upload the CheckM2 database to S3
│   ├── shotgun_checkm2.slurm                # Step 6c — Quality control of the 3 binners' MAGs (CheckM2)
│   ├── pull_gtdbtk.slurm                    # Step 7a — Install the GTDB-Tk Apptainer image
│   ├── build_gtdbtk.slurm                   # Step 7b — (one-off) Upload the GTDB-Tk R214 database to S3
│   ├── shotgun_gtdbtk.slurm                 # Step 7c — Taxonomic assignment of the 3 binners' MAGs (GTDB-Tk)
│   ├── setup_dram_db.slurm                  # Step 8a — (maintenance) Rebuild the viral + peptidase DRAM sub-databases on S3
│   ├── shotgun_DRAM.slurm                   # Step 8b — Functional annotation of one MAG (any binner), one job per bin
│   └── ...                                  # Step 9 — MAG abundance quantification
│
├── containers/                        # Singularity/Apptainer images (.sif)
├── resources/                        # Bowtie2 index, reference genome, databases
├── data/
│   └── raw/                          # Raw files (mirrored on S3: shotgun_bruts/)
├── log/                               # SLURM log files (.out / .err)
└── results/
    ├── shotgun_cleaned/               # Host-decontaminated FASTQ + stats (S3: shotgun_cleaned/)
    ├── kraken2_pluspfp/               # Kraken2 reports vs. PlusPFP database
    ├── bracken_pluspfp/               # Bracken species-level abundance tables (PlusPFP)
    ├── kraken2_addon/                 # Kraken2 reports vs. custom Insect Addon database
    ├── bracken_addon/                 # Bracken species-level abundance tables (Addon)
    ├── assembly_megahit/              # Per-sample contigs + MEGAHIT logs
    ├── global_catalog/                # Non-redundant contig catalog (MMseqs2)
    ├── coverm/                        # Per-sample coverage tables + cached BAMs
    ├── metabat2_bins/                 # MAGs from MetaBAT2 + global depth matrix
    ├── semibin2_bins/                 # MAGs from SemiBin2
    ├── concoct_bins/                  # MAGs from CONCOCT
    ├── checkm2_qc/                    # Completeness/contamination reports, one subfolder per binner
    │   ├── metabat2/                  # quality_report.tsv for MetaBAT2 MAGs
    │   ├── semibin2/                  # quality_report.tsv for SemiBin2 MAGs
    │   └── concoct/                   # quality_report.tsv for CONCOCT MAGs
    ├── filtered_bins_metabat2/        # MetaBAT2 MAGs passing the CheckM2 quality filter
    ├── filtered_bins_semibin2/        # SemiBin2 MAGs passing the CheckM2 quality filter
    ├── filtered_bins_concoct/         # CONCOCT MAGs passing the CheckM2 quality filter
    ├── gtdbtk_taxo_metabat2/          # GTDB-Tk taxonomy assignments for MetaBAT2 MAGs
    ├── gtdbtk_taxo_semibin2/          # GTDB-Tk taxonomy assignments for SemiBin2 MAGs
    ├── gtdbtk_taxo_concoct/           # GTDB-Tk taxonomy assignments for CONCOCT MAGs
    ├── dram_metabat2_raw/             # DRAM annotate output, one subfolder per MetaBAT2 MAG
    ├── dram_metabat2_distill/         # DRAM distill output, one subfolder per MetaBAT2 MAG
    ├── dram_semibin2_raw/             # DRAM annotate output, one subfolder per SemiBin2 MAG
    ├── dram_semibin2_distill/         # DRAM distill output, one subfolder per SemiBin2 MAG
    ├── dram_concoct_raw/              # DRAM annotate output, one subfolder per CONCOCT MAG
    └── dram_concoct_distill/          # DRAM distill output, one subfolder per CONCOCT MAG
```

> Note: unlike the 16S pipeline (fully local `data/` and `results/`), the shotgun pipeline uses an **S3 bucket** (`lmge-tenebrion`) as the primary storage for raw and cleaned reads, accessed via `rclone`.

---

## Running the pipeline

```bash
sbatch bin/host_decontamination_bowtie2.slurm
sbatch bin/build_kraken2_db.slurm               # once — builds the PlusPFP database on S3
sbatch bin/build_kraken_custom.slurm            # once — builds the insect addon database on S3
sbatch bin/shotgun_kraken2_bracken_addon.slurm  # per-sample array — dual-database profiling
sbatch bin/shotgun_megahit_full.slurm           # per-sample array — de novo assembly (MEGAHIT)
sbatch bin/pull_coverm.slurm                    # once — pulls the CoverM Apptainer image
sbatch bin/shotgun_coverm.slurm                 # per-sample array — catalog build (leader task) + coverage profiling
sbatch bin/pull_binners.slurm                   # once — pulls the MetaBAT2 + CONCOCT + SemiBin2 Apptainer images
sbatch bin/shotgun_metabat2.slurm               # single job — binning with MetaBAT2
sbatch bin/shotgun_semibin2.slurm               # single job — binning with SemiBin2
sbatch bin/shotgun_concoct.slurm                # single job — binning with CONCOCT
sbatch bin/shotgun_maxbin2.slurm                # single job — binning with MaxBin2
sbatch bin/pull_checkm2.slurm                   # once — pulls the CheckM2 image + downloads the ML database locally
sbatch bin/build_checkm2_db.slurm               # once — re-uploads the CheckM2 database to S3, for cluster-wide reuse
sbatch bin/shotgun_checkm2.slurm                # single job — QC of all 3 binners' MAGs (CheckM2)
sbatch bin/pull_gtdbtk.slurm                    # once — pulls the GTDB-Tk Apptainer image
sbatch bin/build_gtdbtk.slurm                   # once — uploads the GTDB-Tk R214 database to S3
sbatch bin/shotgun_gtdbtk.slurm                 # single job — taxonomic assignment of all 3 binners' MAGs (GTDB-Tk)
sbatch bin/setup_dram_db.slurm                  # maintenance — (re)builds the viral + peptidase DRAM sub-databases on S3

# one sbatch per MAG, for each binner (see loop example below)
sbatch bin/shotgun_DRAM.slurm metabat2 bin.42
sbatch bin/shotgun_DRAM.slurm semibin2 SemiBin_125
sbatch bin/shotgun_DRAM.slurm concoct   71
```

`shotgun_kraken2_bracken.slurm` (PlusPFP-only) is the earlier, single-database version of the profiling step; `shotgun_kraken2_bracken_addon.slurm` supersedes it once the insect addon database is available, since it reproduces the PlusPFP run and adds the addon pass in the same job.

Step 1 is a **self-submitting SLURM Array**: launched once without `--array`, it auto-detects the sample list on S3 and resubmits itself with the correct array range. Unlike Steps 1–4, the three binning tools in Step 5, the CheckM2 QC in Step 6, and GTDB-Tk in Step 7 are **single jobs, not arrays** — they each process the whole catalog / all bins together in one run, since these steps inherently work across samples rather than per sample. **Step 8 (DRAM) breaks that pattern again**: `DRAM.py annotate` only processes one genome per invocation, so `shotgun_DRAM.slurm` is launched once per MAG (one `sbatch` call per bin, looped manually or via a small wrapper), rather than as a single job or a SLURM array over all bins at once.

---

## Script descriptions

---

### Step 1 — `host_decontamination_bowtie2.slurm`: Host read removal

**Objective:** Remove *Tenebrio molitor* host reads from raw shotgun sequencing data before any taxonomic or functional analysis, keeping only microbiome-derived reads.

**Dynamic array submission.** The script first lists all `*_R1.fastq.gz` samples available on the S3 bucket (`shotgun_bruts/`) via `rclone ls`. When launched without a SLURM array task ID, it counts the samples and resubmits itself as a job array (`--array=0-N-1`), one task per sample.

**Safe, concurrency-proof index construction.** Since all array tasks start simultaneously, the first task to acquire an atomic lock (`mkdir` on a dedicated lock folder) downloads the *T. molitor* reference genome (`GCA_963966145.1_icTenMoli1.1`) from S3 and builds the Bowtie2 index. Other tasks detect the lock and wait (polling every 30 s) until the index is available, avoiding redundant downloads/builds and race conditions.

**Local scratch execution.** Each array task works in a dedicated folder under `/storage/scratch/${USER}/`, automatically cleaned up on exit (`trap cleanup EXIT`), to avoid I/O bottlenecks and leftover files on the shared filesystem.

**Host mapping (Bowtie2).** Paired-end reads are downloaded from S3 and aligned against the indexed host genome using `bowtie2 --very-sensitive`, piped directly into a BAM file via `samtools view`.

**Microbiome extraction (samtools).** Read pairs where *both* mates are unmapped (`-f 12`) are kept — i.e. strictly host-free pairs — and exported as clean paired FASTQ files (`samtools fastq`).

**Decontamination report.** For each sample, the script computes the total read count, host read count/percentage, and microbiome read count/percentage, saved as a `<sample>_decontam_report.txt` table.

**Export to S3.** Clean FASTQ files and the per-sample report are pushed back to S3, under `shotgun_cleaned/fastq/` and `shotgun_cleaned/stats/` respectively.

---

### Step 2 — Kraken2 / Bracken: Read-based taxonomic profiling

**Objective:** Assign a taxonomy to every decontaminated read pair and estimate species-level relative abundances directly from reads (no assembly), enabling fast whole-community comparisons across developmental stages.

**Database construction — `build_kraken2_db.slurm`.** Downloads the pre-built Kraken2 **PlusPFP** database (Standard + Protozoa + Fungi + Plant, `k2_pluspfp_20240112`, ~75 GB) from the official Kraken2 index mirror, extracts it on `/storage/scratch` to avoid saturating `/home`, and uploads the resulting `.k2d` files to S3 (`ref/kraken2_pluspfp/`).

**Custom database construction — `build_kraken_custom.slurm`.** Builds a **custom "insect addon" database** targeting taxa poorly represented in PlusPFP but expected in the *T. molitor* microbiome: entomopathogenic fungi (*Beauveria*, *Metarhizium*, *Cordyceps*, *Ophiocordyceps*, *Lecanicillium*...), insect-associated yeasts, insect viruses (Baculoviridae, Iridovirus, Alphaentomopoxvirus...), gut protozoan parasites (*Gregarina*, *Nosema*, *Crithidia*, *Leptomonas*...) and archaea (*Methanobrevibacter* spp.). Genomes are retrieved via the **NCBI Datasets CLI**, and TaxIDs are injected directly into FASTA headers (`kraken:taxid|...`) to **bypass** the standard NCBI accession2taxid mapping — a lightweight taxonomy (`nodes.dmp`/`names.dmp` only, with empty accession2taxid placeholders) is enough for `kraken2-build` to compile the database, considerably speeding up construction. The database is built and cleaned with `kraken2-build`, then uploaded to S3 (`ref/kraken2_insect_addon/`). Both database scripts run Kraken2 tools through the same Apptainer image already built for the 16S pipeline (`16_tenebrion/containers/kraken2_bracken.sif`).

**Per-sample profiling, single database — `shotgun_kraken2_bracken.slurm`.** Initial version of the profiling step, run as a **SLURM Array** (one task per sample, self-submitting like the decontamination script). The PlusPFP database is synced once from S3 to a shared scratch location, guarded by a `flock`-based lock and a `.sync_done` flag file so concurrent array tasks don't duplicate the transfer. Kraken2 assigns paired, gzip-compressed reads (`--paired --gzip-compressed --use-names`) against the database, then Bracken re-estimates species-level abundances (`-l S`, assumed read length 150 bp) from the resulting Kraken2 report.

**Per-sample profiling, dual database — `shotgun_kraken2_bracken_addon.slurm`.** Extended version run once the insect addon database is available: profiles each sample against **both databases** (PlusPFP and Insect Addon) within the same job, each with its own scratch sync/lock, producing two independent Kraken2/Bracken report pairs per sample. The addon Bracken step is wrapped in `|| true` since `bracken-build` may not yet have been run on every addon database revision, allowing the job to complete even if that particular re-estimation fails.

**Outputs**, copied to `results/`:

| Folder | Content |
|---|---|
| `kraken2_pluspfp/` | Kraken2 reports (`.kreport`) vs. the PlusPFP database |
| `bracken_pluspfp/` | Bracken species-level abundance tables and re-estimated reports (PlusPFP) |
| `kraken2_addon/` | Kraken2 reports vs. the custom Insect Addon database |
| `bracken_addon/` | Bracken species-level abundance tables and re-estimated reports (Addon) |

---

### Step 3 — `shotgun_megahit_full.slurm`: Per-sample de novo assembly

**Objective:** Reconstruct contigs from decontaminated reads via de Bruijn graph assembly, independently for each sample (no co-assembly across samples).

**Dynamic array submission.** As in Step 1, the script lists all `*_clean_R1.fastq.gz` samples available on S3 (`shotgun_cleaned/fastq/`) via `rclone ls`. When launched without a SLURM array task ID, it counts the samples and resubmits itself as a job array (`--array=0-N-1`), one task per sample.

**Local scratch execution.** Each array task works in a dedicated folder under `/storage/scratch/${USER}/`, automatically cleaned up on exit (`trap cleanup EXIT`).

**Full-file retrieval.** Unlike the on-the-fly streaming used in the 16S cleaning step, the entire cleaned FASTQ pair for the sample is downloaded from S3 before assembly.

**Assembly (MEGAHIT).** Reads are assembled with the `meta-sensitive` preset, tuned for complex, low-abundance metagenomic communities at the cost of longer runtime. Each sample is assembled **independently** — reads from different samples are never pooled, so this is a per-sample assembly strategy rather than a co-assembly.

**Export to project directory.** The resulting contigs (`<sample>_megahit.contigs.fa`) and the MEGAHIT log are copied to `results/assembly_megahit/`.

| Folder | Content |
|---|---|
| `assembly_megahit/` | Per-sample contigs (`<sample>_megahit.contigs.fa`) + MEGAHIT logs |

---

### Step 4 — `pull_coverm.slurm` + `shotgun_coverm.slurm`: Contig catalog & coverage profiling

**Objective:** Merge all per-sample MEGAHIT assemblies into a single non-redundant contig catalog, then quantify the coverage of every sample's reads against that catalog — producing the cross-sample differential coverage signal that downstream binning tools rely on.

**Image installation — `pull_coverm.slurm`.** Pulls the CoverM 0.7.0 Apptainer image from the Biocontainers registry, with the Apptainer cache/tmp redirected to scratch to avoid filling `/home`.

**Dynamic array submission — `shotgun_coverm.slurm`.** As in the previous steps, the script lists all cleaned samples on S3 and resubmits itself as a SLURM array (one task per sample) when launched without an array task ID.

**Catalog construction (mutex/leader task).** Before any per-sample work, exactly one array task builds the global catalog, guarded by an atomic `mkdir`-based lock — other tasks poll every 30 s until a `.catalog_done.lock` flag file appears:
- All per-sample MEGAHIT contigs are concatenated, with each contig header prefixed by its sample name (e.g. `>k141_1` → `>AR_E40_k141_1`) to prevent ID collisions across samples.
- Contigs shorter than 1,500 bp are discarded via an explicit `awk` filter on the merged FASTA — a deliberate cutoff applied before clustering, not an automatic behaviour of MMseqs2.
- The filtered set is dereplicated with `mmseqs easy-linclust` (`--min-seq-id 0.95 -c 0.90 --cov-mode 1`): contigs at least 95% identical over 90% of the longer sequence's length are collapsed, keeping one representative per cluster. This produces `non_redundant_catalog.fasta`.

**Per-sample coverage profiling (array task).** Once the catalog is available, each task downloads its sample's full cleaned read pair from S3 into an isolated scratch folder (`trap cleanup EXIT`) and runs CoverM (`coverm contig`) inside the Apptainer image: reads are mapped onto the catalog with `minimap2-sr`, filtered at 95% read identity, and summarised as `trimmed_mean` and `count` per contig. BAM files are cached alongside the per-sample coverage table for reuse by downstream steps.

| Folder | Content |
|---|---|
| `global_catalog/` | Non-redundant contig catalog (`non_redundant_catalog.fasta`) |
| `coverm/` | Per-sample coverage tables (`<sample>_coverage.tsv`) + cached BAMs (`bams/`) |

---

### Step 5 — Binning: MetaBAT2, MaxBin2, SemiBin2 & CONCOCT

**Objective:** From the non-redundant contig catalog (Step 4) and its per-sample coverage BAMs, recover draft genomes (MAGs) using three independent binning strategies run in parallel, to be cross-compared and quality-checked downstream (Step 6).

**Image installation — `pull_binners.slurm`.** Pulls all four binning Apptainer images in a single job: MetaBAT2 2.15, MaxBin2 2.2.7, CONCOCT 1.1.0, and SemiBin2 2.1.0, all from Biocontainers. The Apptainer cache/tmp is redirected to scratch and removed once every image is downloaded and version-tested.

**MetaBAT2 — `shotgun_metabat2.slurm`.** A single SLURM job (no array) processing the whole catalog and all per-sample BAMs at once:
- *Depth matrix.* `jgi_summarize_bam_contig_depths` combines every sorted BAM from `results/coverm/bams/` into one global coverage matrix (`global_depth.txt`); BAMs are listed alphabetically first so the column order stays consistent across reruns.
- *Catalog re-filter.* The catalog is filtered again to ≥1,500 bp (`catalog_min1500.fasta`) immediately before binning — a defensive step, since the catalog from Step 4 is already filtered at the same threshold, but MetaBAT2 requires the guarantee at read time and won't run on a file containing shorter contigs.
- *Binning.* `metabat2` clusters contigs by TNF composition + differential coverage, with the minimum contig size (`-m 1500`) matching the pre-filter.
- *Output:* one `bin.<N>.fa` file per recovered MAG in `results/metabat2_bins/`.

**SemiBin2 — `shotgun_semibin2.slurm`.** A single SLURM job using the self-supervised deep-learning binning mode (`single_easy_bin`):
- *Pre-flight check.* Confirms every BAM in `results/coverm/bams/` has a matching `.bai` index (mandatory for SemiBin2), aborting with an explicit error otherwise.
- *Catalog re-filter.* Same ≥1,500 bp defensive filter as MetaBAT2, cached (`catalog_min1500.fasta`) so it isn't recomputed on reruns.
- *Binning.* `SemiBin2 single_easy_bin` jointly learns from k-mer composition and multi-sample coverage via a self-supervised neural network, instead of the manual TNF/coverage heuristics used by MetaBAT2 — generally better suited to separating closely related strains.
- *Output:* bins are copied from SemiBin2's `output_bins/` into a standardised `results/semibin2_bins/bins_fasta/` folder.

**CONCOCT — `shotgun_concoct.slurm`.** A single SLURM job implementing CONCOCT's full multi-step workflow:
- *Pre-flight check & catalog re-filter:* same BAM index check and ≥1,500 bp filter as above.
- *Contig chunking* (`cut_up_fasta.py`, 10 kb windows, `--merge_last`): CONCOCT clusters fixed-size fragments rather than full-length contigs, to stabilise its Gaussian mixture model.
- *Coverage table* (`concoct_coverage_table.py`) computed on the chunked contigs against all BAMs.
- *Clustering* (`concoct`) via Gaussian mixture models on composition + coverage.
- *Fragment re-merging* (`merge_cutup_clustering.py`) reassembles original contigs from their chunk-level cluster assignments.
- *Bin extraction* (`extract_fasta_bins.py`) writes one FASTA per cluster; outputs are renamed from `.fasta` to `.fa` for consistency with MetaBAT2/CheckM2 conventions.
- *Output:* `results/concoct_bins/bins_fasta/`.

**MaxBin2 — `shotgun_maxbin2.slurm`.** A single SLURM job (no array) that reuses MetaBAT2's depth-calculation tool before running MaxBin2 itself, since MaxBin2 has no built-in coverage estimator:
- *Pre-flight check.* Confirms every BAM in `results/coverm/bams/` has a matching `.bai` index, aborting with an explicit error otherwise.
- *Catalog re-filter.* Same ≥1,500 bp defensive filter as the other binners, cached as `catalog_min1500.fasta`.
- *Depth matrix.* `jgi_summarize_bam_contig_depths` (called from the MetaBAT2 image) computes a global depth matrix (`depth_matrix.txt`) across all BAMs.
- *Abundance files.* The depth matrix is split into one 2-column (contig, depth) file per sample under `abundance_files/`, listed in `abund_list.txt` — the format MaxBin2 expects via `-abund_list`.
- *Binning.* `run_MaxBin.pl` clusters contigs using an Expectation-Maximization algorithm on tetranucleotide frequency composition and multi-sample abundance.
- *Output:* bins (`maxbin2_out.<NNN>.fasta`) are renamed to `.fa` and copied into `results/maxbin2_bins/bins_fasta/`, standardised like the other three binners.

**Outputs:**

| Folder | Content |
|---|---|
| `metabat2_bins/` | MAGs from MetaBAT2 (`bin.<N>.fa`) + global depth matrix |
| `semibin2_bins/bins_fasta/` | MAGs from SemiBin2 |
| `concoct_bins/bins_fasta/` | MAGs from CONCOCT |
| `maxbin2_bins/bins_fasta/` | MAGs from MaxBin2 + depth matrix & per-sample abundance files |

---

### Step 6 — Bin quality control: CheckM2

**Objective:** Assess the completeness and contamination of every MAG recovered by the three binners (Step 5), producing one comparable quality report per tool, used in Step 7 to decide which MAGs are worth a taxonomic assignment.

**Image installation — `pull_checkm2.slurm`.** Pulls the CheckM2 1.0.2 Apptainer image from Biocontainers, then downloads CheckM2's Machine Learning database (`checkm2 database --download`, ~3 GB Diamond `.dmnd` file) directly into a local `databases/CheckM2_database/` folder.

**Database upload — `build_checkm2_db.slurm`.** A one-off companion script: re-downloads the same database into scratch via the Apptainer image, then uploads it to S3 (`ref/checkm2_db/`) so any compute node can retrieve it on demand, rather than relying on the local copy pulled by `pull_checkm2.slurm`.

**QC across all three binners — `shotgun_checkm2.slurm`.** A single SLURM job (no array) that evaluates all three bin sets in one run:
- *Database sync.* The `.dmnd` database is rclone-copied once from S3 into scratch (`trap cleanup EXIT`) and reused for all three CheckM2 runs below — avoiding three redundant ~3 GB downloads.
- *Per-binner loop.* The script iterates over the three bin directories (`metabat2_bins/`, `semibin2_bins/bins_fasta/`, `concoct_bins/bins_fasta/`). For each, it counts the `.fa` bins present; if a binner produced none (e.g. not yet run), it logs a warning and skips that tool rather than failing the whole job.
- *Prediction (Diamond BlastP + ML).* For each non-empty bin set, `checkm2 predict` searches conserved marker genes via Diamond against the database, then estimates completeness and contamination through CheckM2's gradient-boosting models — without relying on lineage-specific marker sets, unlike the older CheckM1 approach.
- *Per-tool output.* Each binner gets its own `quality_report.tsv` in a dedicated subfolder, keeping the three tools' results directly comparable side-by-side.

**Outputs:**

| Folder | Content |
|---|---|
| `checkm2_qc/metabat2/` | `quality_report.tsv` — completeness/contamination for MetaBAT2 MAGs |
| `checkm2_qc/semibin2/` | `quality_report.tsv` — completeness/contamination for SemiBin2 MAGs |
| `checkm2_qc/concoct/` | `quality_report.tsv` — completeness/contamination for CONCOCT MAGs |

---

### Step 7 — Taxonomic assignment: GTDB-Tk

**Objective:** Assign a taxonomy to the MAGs recovered by each of the three binners, using the completeness/contamination scores from Step 6 to decide which MAGs are passed to GTDB-Tk.

**Image installation — `pull_gtdbtk.slurm`.** Pulls the GTDB-Tk 2.3.2 Apptainer image from the official `ecogenomic/gtdbtk` Docker repository.

**Database upload — `build_gtdbtk.slurm`.** A one-off script: downloads the GTDB-Tk **R214** reference package (~85 GB, marker gene HMMs, reference tree, taxonomy and FastANI genome set) directly from the official GTDB mirror, extracts it, and uploads it to S3 (`ref/gtdbtk_r214/`) for reuse across cluster nodes.

**Taxonomic assignment across all three binners — `shotgun_gtdbtk.slurm`.** A single SLURM job (no array) that processes all three bin sets sequentially, reusing one shared function (`run_gtdbtk_for_binner`) for MetaBAT2, SemiBin2 and CONCOCT:
- *Database sync (once).* The R214 reference package is rclone-copied once from S3 into scratch at the start of the job (`trap cleanup EXIT`) and reused for all three GTDB-Tk runs — avoiding three redundant ~85 GB downloads.
- *Quality-based pre-filter.* For each binner, the corresponding `checkm2_qc/<binner>/quality_report.tsv` (Step 6) is parsed to select which MAGs get a taxonomy attempt. The current threshold is **deliberately permissive** (Completeness ≥ 10%, Contamination < 100%) rather than the stricter MIMAG Medium-Quality standard (Comp ≥ 50% / Cont < 10%) — the goal here is exploratory: get a taxonomic signal for as many recovered bins as possible, including partial ones, rather than only the highest-confidence MAGs. Passing bins are copied into `filtered_bins_<binner>/`.
- *Classification (GTDB-Tk `classify_wf`).* Each filtered bin set is run through `gtdbtk classify_wf --skip_ani_screen`, which places genomes in the GTDB reference tree via concatenated marker gene alignment and resolves species-level calls with FastANI where the tree placement alone isn't conclusive.
- *Sequential execution.* The three runs happen one after another rather than in parallel: GTDB-Tk needs ~150 GB RAM and saturates the allocated CPUs, so running binners concurrently would exceed the job's 250 GB memory budget.
- *Per-binner tolerance.* If a `quality_report.tsv` or bin directory is missing for a given binner (e.g. Step 5/6 not yet run for it), the script logs a warning and moves on to the next binner instead of aborting the whole job.

**Outputs:**

| Folder | Content |
|---|---|
| `filtered_bins_metabat2/` | MetaBAT2 MAGs passing the quality pre-filter |
| `filtered_bins_semibin2/` | SemiBin2 MAGs passing the quality pre-filter |
| `filtered_bins_concoct/` | CONCOCT MAGs passing the quality pre-filter |
| `gtdbtk_taxo_metabat2/` | GTDB-Tk classification output for MetaBAT2 MAGs |
| `gtdbtk_taxo_semibin2/` | GTDB-Tk classification output for SemiBin2 MAGs |
| `gtdbtk_taxo_concoct/` | GTDB-Tk classification output for CONCOCT MAGs |

---

### Step 8 — Functional annotation: DRAM

**Objective:** Annotate the metabolic and functional gene content of each MAG that passed the Step 7 quality filter — central databases (KEGG, dbCAN/CAZymes, MEROPS, Pfam) plus viral (VOG) markers — then distill the raw annotation into pathway-level completeness scores.

**Database maintenance — `setup_dram_db.slurm`.** A maintenance/patch script, not part of the standard per-run pipeline: rebuilds specifically the **viral** and **peptidase** DRAM sub-databases (downloading `viral.merged.protein.faa.gz` and `merops_peptidases_nr.faa` from S3, then compiling them with `DRAM-setup.py prepare_databases --select_db`) and re-uploads the result to S3 (`ref/dram_db/`). Run this only when those two sub-databases specifically need rebuilding — the rest of the DRAM database (KEGG, dbCAN, Pfam...) is assumed already built and present on S3.

**Per-MAG annotation — `shotgun_DRAM.slurm`.** A single script parameterised by binner and bin name (`sbatch bin/shotgun_DRAM.slurm <metabat2|semibin2|concoct> <bin_name>`), merging what used to be three near-identical scripts (one per binner). Unlike Steps 5–7, this step stays **one SLURM job per MAG** rather than a single job or an array over everything at once, since `DRAM.py annotate` only processes one genome per invocation — jobs are submitted in a loop, one `sbatch` call per bin (example loop included at the bottom of the script):
- *Database sync.* The DRAM database is rclone-copied from S3 into scratch, excluding the raw `database_files/` and unused `kofam_profiles/profiles/` subfolders to save transfer time.
- *Config injection.* A dummy CONFIG file is bound over DRAM's internal package config (`mag_annotator/CONFIG`) inside the container, then populated via `DRAM-setup.py set_database_locations`, pointing at whichever database files were actually found on scratch (each location argument is only added if the corresponding file exists, so missing sub-databases — e.g. viral/peptidase before their first build — don't break the run).
- *Runtime patch.* `annotate_bins.py` is copied out of the container and patched with `sed` to guard a Pfam-hits access that otherwise crashes when Pfam annotation is absent for a given gene — the patched copy is bound back over the original inside the container at runtime.
- *Annotation (`DRAM.py annotate`).* Predicts genes (Prodigal) and searches them against all configured databases (MMseqs2/HMMER) for the selected MAG (`--min_contig_size 1000`).
- *Distillation (`DRAM.py distill`).* Cross-references the raw gene hits against DRAM's metabolic pathway maps to produce per-pathway completeness scores, incorporating rRNA/tRNA calls from the annotate step when available.
- *Input source.* Bins are read from `filtered_bins_<binner>/` (Step 7's quality-filtered sets) — consistent across all three binners.

**Outputs:**

| Folder | Content |
|---|---|
| `dram_metabat2_raw/<bin>/` | Raw DRAM annotation (`annotations.tsv`, rRNA/tRNA calls) — MetaBAT2 |
| `dram_metabat2_distill/<bin>/` | Distilled pathway completeness — MetaBAT2 |
| `dram_semibin2_raw/<bin>/` | Raw DRAM annotation — SemiBin2 |
| `dram_semibin2_distill/<bin>/` | Distilled pathway completeness — SemiBin2 |
| `dram_concoct_raw/<bin>/` | Raw DRAM annotation — CONCOCT |
| `dram_concoct_distill/<bin>/` | Distilled pathway completeness — CONCOCT |

---

## Pipeline diagram

```
Raw shotgun data (FASTQ, on S3)
        |
        v
[1] Host decontamination     Bowtie2 alignment vs T. molitor genome + samtools filtering (-f 12)
        |
        +-----------------------------------------------------+
        |                                                     |
        v                                                     v
[2] Read-based profiling                                [3] Assembly       Per-sample de novo assembly (MEGAHIT, meta-sensitive)
    PlusPFP DB   <--[2a] build_kraken2_db                     |
    Insect Addon <--[2b] build_kraken_custom                  |
        |                                                     v
        v                                             [4a] Catalog construction
    [2c/2d] Kraken2 + Bracken (per sample, dual-DB)            Concat + rename + length filter (>=1500bp)
        |                                                      + MMseqs2 dereplication (95% id / 90% cov)
        v                                                     |
Taxonomic abundance tables                                    v
(PlusPFP + Addon)                                    [4b] Coverage profiling
                                                            Per-sample read mapping (CoverM, minimap2-sr)
                                                              |
                                        +---------------------+---------------------+
                                        |                     |                     |
                                        v                     v                     v
                                  [5] MetaBAT2           [5] SemiBin2           [5] CONCOCT
                                  TNF + jgi depth        Self-supervised DL     Chunking + GMM
                                        |                     |                     |
                                        +---------------------+---------------------+
                                                              |
                                                              v
                                                    [6] CheckM2 (completeness/contamination)
                                                        one quality_report.tsv per binner
                                                              |
                                                              v
                                                    [7] GTDB-Tk (quality pre-filter + classify_wf)
                                                        one taxonomy per binner, sequential runs
                                                              |
                                                              |                                        
                                                              v                                    
                                                    [8] DRAM (annotate + distill)            
            
```

---

## Prerequisites

- HPC cluster with **SLURM** job scheduler
- **rclone** configured with an S3 remote (`s3_uca`) pointing to the `lmge-tenebrion` bucket
- **Apptainer/Singularity** available as a system binary, with the following images in `containers/`: `kraken2_bracken.sif` (shared with the 16S pipeline, `16_tenebrion/containers/`), `coverm_v0.7.0.sif`, `metabat2_v2.15.sif`, `concoct_1.1.0.sif`, `semibin_2.1.0.sif`, `checkm2_v1.0.2.sif`, `gtdbtk_v2.3.2.sif`, `dram_v1.4.6.sif`
- CheckM2 ML database (~3 GB `.dmnd` file) hosted on S3 (`ref/checkm2_db/`), built once via `build_checkm2_db.slurm`
- GTDB-Tk R214 reference package (~85 GB) hosted on S3 (`ref/gtdbtk_r214/`), built once via `build_gtdbtk.slurm`
- DRAM database (KEGG, dbCAN, MEROPS, Pfam, VOG, description DB) hosted on S3 (`ref/dram_db/`); the viral and peptidase sub-databases specifically can be rebuilt via `setup_dram_db.slurm`
- Modules available on the cluster: `rclone/1.55.1`, `bowtie2/2.3.4.3`, `gcc/8.1.0`, `samtools/1.16.1`, `megahit/1.2.9`, `MMseqs2/13-45111`
- Reference genome of *Tenebrio molitor* (`GCA_963966145.1_icTenMoli1.1`) hosted on S3 (`ref/genome/`)
- Kraken2 databases hosted on S3: PlusPFP (`ref/kraken2_pluspfp/`) and the custom Insect Addon (`ref/kraken2_insect_addon/`)
- **NCBI Datasets CLI** (fetched automatically by `build_kraken_custom.slurm`) for the addon database genomes
- Internet/S3 access from compute nodes

---

## Author

Thomas BOUTET — Tenebrion Project
