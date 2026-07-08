# Sample Inclusion and Exclusion

## GSE289173

GSE289173 was used as the primary human liver snRNA-seq discovery dataset.

| Analysis layer | Included samples |
|---|---:|
| GEO metadata parsed from series matrix | 48 |
| Final donor-level hepatocyte panel | 44 |

Final donor-level panel by disease group:

| Disease group | Donors |
|---|---:|
| No-MASLD | 23 |
| MASL | 4 |
| eMASH | 7 |
| aMASH | 10 |

Four samples present in GEO metadata were not included in the final donor-level
panel because matched hepatocyte cell-score/pseudobulk inputs were incomplete:
`FL511`, `FL531`, `FL284` and `FL793`.

MA-high hepatocytes were defined globally as the top 25% of hepatocytes by
metabolic-aging score in GSE289173. Donor-internal MA-high and MA-low cells were
defined as the top and bottom 25% of hepatocytes within each donor.

## GSE202379

GSE202379 was used as the SAF activity validation dataset.

Main Fig. 1d uses SAF activity A0-1 versus A2-4 after excluding
healthy/control samples:

| SAF activity group | Samples |
|---|---:|
| A0-1 | 7 |
| A2-4 | 31 |

The two healthy/control samples with structured SAF scores were excluded from
the main Fig. 1d analysis:

| Patient ID | SAF score | Parsed A group |
|---|---|---|
| P30 | S0A1F0 | A1 |
| P98 | S1A1F1 | A1 |

Supplementary Fig. 1f assigns healthy/control samples to the corresponding
parsed A-score group instead of displaying a separate Healthy group:

| Full SAF activity group | Samples |
|---|---:|
| A0 | 1 |
| A1 | 8 |
| A2 | 4 |
| A3 | 18 |
| A4 | 9 |

## GSE213621

GSE213621 was used as an external bulk-liver validation dataset for RPS27L.

| Fibrosis group | Samples |
|---|---:|
| Control | 69 |
| F0-1 | 97 |
| F2 | 107 |
| F3-4 | 95 |

The final external bulk validation includes 368 samples with matched expression
and GEO fibrosis metadata.
