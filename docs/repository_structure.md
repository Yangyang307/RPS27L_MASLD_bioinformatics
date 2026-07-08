# Repository Structure

```text
RPS27L_MASLD_bioinformatics/
├── README.md
├── LICENSE
├── CITATION.cff
├── code/
│   ├── 00_setup_paths_and_functions.R
│   ├── 01_GSE289173_metabolic_aging_and_RPS27L.R
│   ├── 02_GSE289173_MA_high_low_edgeR_GSEA.R
│   ├── 03_GSE202379_SAF_validation.R
│   ├── 04_GSE213621_bulk_RPS27L_validation.R
│   └── 05_make_publication_tables_and_figures.R
├── env/
│   ├── sessionInfo.txt
│   ├── R_package_versions.tsv
│   └── conda_environment.yml
├── data_availability/
│   ├── public_datasets.md
│   ├── input_files.md
│   ├── sample_inclusion_exclusion.md
│   └── gene_sets/
│       └── metabolic_aging_gene_sets.tsv
└── docs/
    ├── methods_overview.md
    ├── output_files.md
    ├── repository_structure.md
    └── zenodo_archiving_checklist.md
```

`data/raw/`, `data/processed/` and `results/` are generated or user-supplied
runtime folders and are excluded from version control.

