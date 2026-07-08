# RPS27L_MASLD_bioinformatics

This repository contains the bioinformatics analysis code used to evaluate
metabolic-aging-associated hepatocyte states and RPS27L expression in public
human MASLD/MASH transcriptomic datasets.

The repository is structured for public release and Zenodo archival. It does
not include raw sequencing matrices, processed Seurat objects, local working
paths, or working notes. Public data accessions and required input files are
listed in `data_availability/`.

## Public Datasets

| Dataset | Use in analysis | Data type |
|---|---|---|
| GSE289173 | Primary discovery cohort for hepatocyte metabolic-aging score, MA-high hepatocyte proportion, RPS27L hepatocyte pseudobulk expression and donor-internal MA-high versus MA-low analyses | Human liver snRNA-seq |
| GSE202379 | Histology activity validation using SAF activity groups | Human liver single-cell/nucleus RNA-seq |
| GSE213621 | External bulk-liver validation of RPS27L across fibrosis groups | Human liver bulk RNA-seq |

## Repository Contents

| Path | Description |
|---|---|
| `code/` | Numbered R scripts for reproducing the final analyses |
| `env/` | R session and package version records |
| `data_availability/` | Public data sources, expected input files, sample inclusion/exclusion rules and gene sets |
| `docs/` | Method overview, output descriptions and Zenodo archiving checklist |

## Running the Analysis

Run scripts from the repository root after placing the required public input
files under `data/raw/` as described in `data_availability/input_files.md`.
The scripts use relative paths by default and write outputs to `results/`.

```bash
Rscript code/01_GSE289173_metabolic_aging_and_RPS27L.R
Rscript code/02_GSE289173_MA_high_low_edgeR_GSEA.R
Rscript code/03_GSE202379_SAF_validation.R
Rscript code/04_GSE213621_bulk_RPS27L_validation.R
Rscript code/05_make_publication_tables_and_figures.R
```

Alternatively, set a repository root explicitly:

```bash
export RPS27L_REPO_ROOT=/path/to/RPS27L_MASLD_bioinformatics
```

## Script Order and Purpose

| Script | Purpose | Main inputs | Main outputs |
|---|---|---|---|
| `00_setup_paths_and_functions.R` | Shared paths, package loading, plotting theme, statistics and helper functions | Repository root and installed R packages | Shared functions sourced by downstream scripts |
| `01_GSE289173_metabolic_aging_and_RPS27L.R` | Compute metabolic-aging module scores, global MA-high calls, donor summaries and hepatocyte RPS27L pseudobulk expression for GSE289173 | GSE289173 GEO series matrix, hepatocyte covariates, sample matrices/barcodes/features, metabolic-aging gene sets | `results/tables/GSE289173_*` |
| `02_GSE289173_MA_high_low_edgeR_GSEA.R` | Perform donor-internal MA-high versus MA-low pseudobulk differential analysis and Hallmark/Reactome GSEA | Outputs from script 01 and matched raw matrices | edgeR DE table, GSEA table, selected pathway table |
| `03_GSE202379_SAF_validation.R` | Validate MA-high hepatocyte enrichment across SAF activity groups | GSE202379 donor-level hepatocyte MA-high summary with SAF metadata | SAF activity source tables and statistics |
| `04_GSE213621_bulk_RPS27L_validation.R` | Analyze bulk-liver RPS27L expression across fibrosis groups | GSE213621 FPKM matrix and GEO series matrix | RPS27L sample-level table, group statistics and pairwise tests |
| `05_make_publication_tables_and_figures.R` | Assemble final source tables and publication-style plots from prior outputs | Outputs from scripts 01-04 | `results/figures/` and `results/source_data/` |

The candidate metabolic-aging gene-set table is kept in
`data_availability/gene_sets/metabolic_aging_gene_sets.tsv`. `RPS27L` is
excluded in `01_GSE289173_metabolic_aging_and_RPS27L.R` before MA score
calculation, so downstream RPS27L expression analyses are not part of the score
definition.

## Software Versions

The analysis was prepared with R 4.5.3. Package versions are recorded in:

- `env/sessionInfo.txt`
- `env/R_package_versions.tsv`
- `env/conda_environment.yml`

## Data Availability

Raw public datasets should be downloaded from GEO or the corresponding public
repository linked by the GEO records. This repository intentionally excludes raw
matrices and large processed objects. See:

- `data_availability/public_datasets.md`
- `data_availability/input_files.md`
- `data_availability/sample_inclusion_exclusion.md`

## License

Code is released under the MIT License. Public datasets should be used according
to the terms of their original repositories and publications.
