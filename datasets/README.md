# Datasets for the BRCA MPL IBI-DT Pilot

This folder contains the TCGA/Xena-derived inputs needed to audit and rerun the
BRCA MPL pilot, except for the very large raw expression matrix.

## Included in Git

`tcga_brca_inputs/` contains:

- `brca_deg_outputs/brca_DEG_binary_matrix.csv.gz`
- `brca_deg_outputs/brca_normal_baseline_by_gene.csv`
- `brca_deg_outputs/brca_DEG_gene_summary.csv`
- `brca_deg_outputs/brca_selected_samples.csv`
- `brca_deg_outputs/brca_DEG_run_summary.txt`
- `brca_gistic_thresholded_by_genes.gz`
- `mc3.v0.2.8.PUBLIC.nonsilentGene.xena.gz`
- `TCGA_phenotype_denseDataOnlyDownload.tsv.gz`
- `TCGA.PANCAN.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz`
- `04_make_brca_deg_binary_matrix.py`
- `DOWNLOADS.md`

## Large Raw Expression Matrix

The raw TCGA/Xena expression file used to construct the BRCA DEG matrix is:

```text
tcga_RSEM_Hugo_norm_count.gz
```

It is about 892 MB, so it is too large for normal GitHub git storage. It is
uploaded as a GitHub Release asset when available.

The pilot output already includes the derived BRCA DEG matrix and MPL-specific
audit files, so reviewers can inspect the current MPL run without downloading
the full raw expression matrix.
