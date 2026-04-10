# Forest Soil Microbiome Analysis
Comparative analysis of microbial communities in tropical and temperate forest soils, using 16S rRNA amplicon data from the [Earth Microbiome Project (EMP)](https://doi.org/10.1038/nature24621).

## Portfolio Context

This is Project 3 in a series of soil biogeochemistry portfolio projects I am building as part of my transition from nanotechnology into soil science.

- **Project 1:** [Interactive CO₂ emissions dashboard (R Shiny)](https://vinmk-shiny.shinyapps.io/co2_viz_0/)
- **Project 2:** [Forest soil respiration modelling under climate warming (Python)](https://github.com/vinmk-git/Soil-Respiration-Model)
- **Project 3:** This project (Python, R)
- **Project 4:** [Arbuscular mycorrhiza meta-analysis](https://github.com/vinmk-git/mycorrhiza)

More at [vinilegueuno.com/portfolio](https://vinilegueuno.com/portfolio) and [github.com/vinmk-git](https://github.com/vinmk-git).

## Overview

The central question: do tropical and temperate forest soils host measurably different microbial communities — and if so, how?

The fix_biom.ipynb Jupyter notebook loads the EMP Release 1 biom dataset, filters for forest soils, then extracts a phylogenetic tree and taxonomy table, and finally creates separate OTU (observed taxonomic units) tables and metadata files for temperate and tropical forests.

The R pipeline then loads OTU tables, taxonomy, and a phylogenetic tree, combines tropical and temperate soil datasets into a single phyloseq object, and runs alpha diversity, beta diversity, and differential abundance analyses.

## Data

All source data is from the EMP Release 1 dataset [1]

Publicly available at: https://zenodo.org/records/890000

## Data preparation (Python)

Before running the R pipeline, the raw EMP data needs to be extracted and reformatted. The `fix_biom.ipynb` notebook handles this: it downloads the four required EMP Release 1 files from Zenodo, loads the full BIOM table, and filters the metadata to forest soil samples only — excluding boreal forests. It then splits the remaining samples into tropical and temperate subsets based on the `env_biome` field, filters the BIOM table to match, and exports the OTU count tables (`tropical.tsv`, `temperate.tsv`), sample metadata (`tropical_metadata.tsv`, `temperate_metadata.tsv`), and a taxonomy table (`taxonomy.tsv`) as tab-separated files ready for phyloseq. The phylogenetic tree (`silva_123.97_otus.tre`) is extracted directly from the `emp_observation_info_cr_silva.tar.gz` archive. Run this notebook once before running the R pipeline.

| Output Files | Description |
|---|---|
| `tropical.tsv` | OTU count table for tropical forest soil samples |
| `temperate.tsv` | OTU count table for temperate forest soil samples |
| `tropical_metadata.tsv` | Sample metadata for tropical samples |
| `temperate_metadata.tsv` | Sample metadata for temperate samples |
| `taxonomy.tsv` | Taxonomy table for all OTUs (SILVA 123 reference) |
| `silva_123.97_otus.tre` | Phylogenetic tree (from `emp_observation_info_cr_silva.tar.gz`) |

## Dependencies

```r
install.packages(c("tidyverse", "ggplot2", "car", "ggpubr"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("phyloseq", "vegan", "DESeq2", "microbiome", "ggtree"))
```

## R Pipeline

### 1. Data loading and filtering
OTU tables are filtered to remove samples with fewer than 558 reads (chosen based on the distribution of per-sample read counts across both datasets) and taxa with zero total counts across remaining samples. Tropical and temperate tables are then merged using a full outer join, filling absent OTUs with zero.

### 2. Phyloseq object construction
A combined phyloseq object is built from the merged OTU table, taxonomy, sample metadata (with a `biome` column indicating Tropical/Temperate), and the EMP phylogenetic tree. Phyloseq prunes the full EMP tree automatically to the OTUs present in the dataset — no manual pruning needed.

### 3. Rarefaction
Data is rarefied to even depth (`sample.size = 10727`) — the minimum sample depth after filtering — using `rarefy_even_depth()`. The effect of rarefaction is validated by comparing Bray-Curtis dissimilarities before and after.

### 4. Alpha diversity
Richness (observed OTUs, Chao1) and evenness (Shannon, Inverse Simpson) are calculated on the rarefied OTU table. I tested normality with Shapiro-Wilk and homoscedasticity with Levene's test — assumptions were violated across the board, so group differences are tested with unpaired Wilcoxon tests.

| Metric | p-value |
|---|---|
| Observed OTUs | 0.0006 |
| Shannon Index | 1.3 × 10⁻⁵ |
| Inverse Simpson | 8.5 × 10⁻⁵ |
| Chao1 | 0.031 |

All four metrics show significant differences between biomes.

### 5. Taxonomic composition
OTUs are agglomerated to phylum level and plotted as relative abundance stacked bar charts, showing the top 15 most abundant phyla across biomes.

### 6. Beta diversity
Three complementary beta diversity methods are applied to the rarefied dataset:

- **Unweighted UniFrac** (PCoA) — presence/absence based, phylogeny-aware
- **Weighted UniFrac** (PCoA) — abundance-weighted, phylogeny-aware
- **Bray-Curtis NMDS** — abundance-based, non-phylogenetic

Group separation is tested with PERMANOVA (`adonis2`) on weighted UniFrac distances. Homogeneity of dispersion is checked with `betadisper`/`permutest`.

### 7. Differential abundance
Differential abundance between biomes is assessed using DESeq2 on the **unrarefied** count table with a Wald test.
DESeq2 handles its own normalisation internally.

## Outputs

| File | Description |
|---|---|
| `AlphaDiversity.png` | 2×2 panel of alpha diversity violin/boxplots |
| `top15_plot.png` | Relative abundance bar chart by biome (top 15 phyla) |
| `unwt_unifrac.png` | PCoA plot — unweighted UniFrac |
| `wt_unifrac.png` | PCoA plot — weighted UniFrac |
| `nmds_plot.png` | NMDS plot — Bray-Curtis |

## Notes on methods

A few decisions worth flagging explicitly:

- **Rarefaction vs. DESeq2 normalisation:** rarefaction is used for alpha and beta diversity; DESeq2 uses the unrarefied table for its own normalisation. This is intentional and consistent with current best practice.
- **558 read threshold:** chosen empirically from the per-sample read count distribution. Samples below this threshold were outliers with unusually low sequencing depth.
- **Full EMP tree:** the full SILVA tree (~100,000 OTUs) is passed to phyloseq, which prunes it to the OTUs in the dataset. This avoids any manual filtering errors.

## References

[1] Thompson, L.R., Sanders, J.G., McDonald, D., et al. (2017). A communal catalogue reveals Earth's multiscale microbial diversity. *Nature*, 551, 457–463. https://doi.org/10.1038/nature24621
