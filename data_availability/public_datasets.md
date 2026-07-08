# Public Datasets

This analysis uses public, de-identified transcriptomic datasets. Raw matrices
and large processed objects are not included in this repository.

| Accession | Dataset role | Data type | Use in this repository |
|---|---|---|---|
| GSE289173 | Primary discovery dataset | Human liver snRNA-seq | Hepatocyte metabolic-aging scoring, MA-high hepatocyte proportion, RPS27L hepatocyte pseudobulk expression, donor-internal MA-high versus MA-low edgeR/GSEA |
| GSE202379 | Histology validation dataset | Human liver single-cell/nucleus RNA-seq | SAF activity validation of MA-high hepatocyte proportions |
| GSE213621 | External validation dataset | Human liver bulk RNA-seq | RPS27L expression across external fibrosis groups |

The repository expects public input files to be placed under `data/raw/` before
running the scripts. The raw data should be downloaded from GEO or the public
repository linked by each GEO accession.

The metabolic-aging gene set used by the scripts is included in
`data_availability/gene_sets/metabolic_aging_gene_sets.tsv`.
