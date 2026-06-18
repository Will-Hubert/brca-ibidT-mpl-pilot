# TCGA Pan-Cancer Downloads

Downloaded from UCSC Xena on 2026-05-29.

## RNA-seq datasets

| Teacher label | Local file | Source |
|---|---|---|
| Batch effects normalized mRNA data | `EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz` | `https://pancanatlas.xenahubs.net/download/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz` |
| TOIL RSEM expected_count | `tcga_gene_expected_count.gz` | `https://toil.xenahubs.net/download/tcga_gene_expected_count.gz` |
| TOIL RSEM fpkm | `tcga_RSEM_gene_fpkm.gz` | `https://toil.xenahubs.net/download/tcga_RSEM_gene_fpkm.gz` |
| TOIL RSEM norm_count | `tcga_RSEM_Hugo_norm_count.gz` | `https://toil.xenahubs.net/download/tcga_RSEM_Hugo_norm_count.gz` |
| TOIL RSEM tpm | `tcga_RSEM_gene_tpm.gz` | `https://toil.xenahubs.net/download/tcga_RSEM_gene_tpm.gz` |

The current BRCA IBI-DT pilot uses `tcga_RSEM_Hugo_norm_count.gz`. The other
large RNA-seq files are archived under `archive_not_used_for_current_ibidT/`.

## Genomic alteration datasets

| Purpose | Local file | Source |
|---|---|---|
| GISTIC2 thresholded gene-level copy number | `TCGA.PANCAN.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz` | UCSC Xena Pan-Cancer Atlas Hub, `gene-level copy number (gistic2)` |
| MC3 gene-level non-silent mutation 0/1 calls | `mc3.v0.2.8.PUBLIC.nonsilentGene.xena.gz` | UCSC Xena Pan-Cancer Atlas Hub, `Gene level non-silent mutation` |
| BRCA GISTIC2 subset used for faster runs | `brca_gistic_thresholded_by_genes.gz` | Extracted BRCA subset of the thresholded GISTIC2 file; verified against the full file by gene and TCGA sample barcode. |

## Metadata files

| Purpose | Local file | Source |
|---|---|---|
| Pan-Cancer phenotype/sample annotation | `TCGA_phenotype_denseDataOnlyDownload.tsv.gz` | `https://pancanatlas.xenahubs.net/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz` |
| Pan-Cancer curated survival/clinical table | `Survival_SupplementalTable_S1_20171025_xena_sp` | `https://pancanatlas.xenahubs.net/download/Survival_SupplementalTable_S1_20171025_xena_sp` |
| Toil TCGA/TARGET/GTEx phenotype annotation | `TcgaTargetGTEX_phenotype.txt.gz` | `https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz` |

All `.gz` files passed `gzip -t` after download.
