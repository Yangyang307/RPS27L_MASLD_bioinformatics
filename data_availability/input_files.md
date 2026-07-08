# Expected Input Files

Place public input files under the following paths before running the analysis.
The paths are relative to the repository root.

## GSE289173

```text
data/raw/GSE289173/GSE289173_series_matrix.txt.gz
data/raw/GSE289173/covariates.hepatocyte.txt.gz
data/raw/GSE289173/features/GSM*_features.tsv.gz
data/raw/GSE289173/barcodes/GSM*_barcodes.tsv.gz
data/raw/GSE289173/matrix/GSM*_matrix.mtx.gz
```

`01_GSE289173_metabolic_aging_and_RPS27L.R` parses the GEO series matrix to map
sample titles, GEO accessions, disease groups and matrix/barcode file names.
If files are missing locally, downloads are blocked by default. Set
`RPS27L_ALLOW_DOWNLOAD=1` only when scripted downloading from GEO is intended.

## GSE202379

```text
data/processed/GSE202379_donor_level_summary.tsv
```

or:

```text
data/raw/GSE202379/donor_level_summary.tsv
```

Required columns:

```text
patient_id
disease_clean
saf_score_raw
ma_high_pct
metabolic_aging_score
fibrosis_stage
```

The GSE202379 table is a donor-level summary used for SAF validation. It should
be generated from the public GSE202379 single-cell/nucleus data using the same
MA-high definition as the primary analysis.

## GSE213621

```text
data/raw/GSE213621/GSE213621_FPKMs_allsamples.txt.gz
data/raw/GSE213621/GSE213621_series_matrix.txt.gz
```

`04_GSE213621_bulk_RPS27L_validation.R` extracts the RPS27L FPKM row and parses
the GEO series matrix fibrotic-stage metadata.
