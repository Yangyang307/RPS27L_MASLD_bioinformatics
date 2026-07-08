# Methods Overview

## Metabolic-Aging Score

For GSE289173 hepatocytes, genes in the metabolic-aging gene-set table were
grouped into stress/senescence modules and hepatocyte functional modules. Counts
were normalized per cell, log transformed and summarized as module scores.
The candidate gene-set table was retained as listed in
`data_availability/gene_sets/metabolic_aging_gene_sets.tsv`. Before score
calculation, RPS27L was excluded from the scoring gene set to avoid circular
association testing between RPS27L expression and the metabolic-aging score.
Module scores were z-transformed across hepatocytes. The final metabolic-aging
score was calculated as:

```text
stress/senescence score - hepatocyte function score
```

Global MA-high hepatocytes were defined as the top 25% of all GSE289173
hepatocytes by metabolic-aging score. Donor-internal MA-high and MA-low
hepatocytes were defined as the top and bottom 25% of hepatocytes within each
donor.

## RPS27L Pseudobulk

For GSE289173, RPS27L hepatocyte expression was computed by summing raw counts
across matched hepatocytes per donor and converting to log1p CPM. Donor-internal
MA-high and MA-low pseudobulk counts were used for paired RPS27L comparison and
edgeR differential expression.

## MA-High Versus MA-Low Differential Analysis

Raw counts from donor-internal MA-high and MA-low hepatocytes were aggregated
into pseudobulk samples. edgeR quasi-likelihood testing was run with donor as a
blocking term:

```text
~ donor + MA_state
```

The MA-high coefficient was used for differential expression ranking. Hallmark
and Reactome GSEA were run with genes ranked by edgeR log fold change.

## SAF Validation

GSE202379 donor summaries were parsed by SAF activity score. Main Fig. 1d
compares A0-1 versus A2-4 after excluding healthy/control samples. Supplementary
Fig. 1f assigns healthy/control samples to the corresponding parsed A-score
group.

## External Bulk Validation

For GSE213621, RPS27L FPKM values were extracted from the public bulk RNA-seq
matrix and transformed as log2(FPKM + 1). Fibrotic-stage metadata were parsed
from the GEO series matrix and summarized as Control, F0-1, F2 and F3-4.
