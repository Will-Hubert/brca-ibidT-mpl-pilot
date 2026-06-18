#!/usr/bin/env python3
"""
Create a BRCA tumor-by-gene DEG binary matrix following the TCI preprocessing idea.

Default input is the UCSC Xena TCGA RSEM Hugo norm_count matrix. These Xena
RSEM values are already log2-based, so the script does not apply log2(x + 1)
again unless --force-log2p is supplied.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_Z_CUTOFF_FOR_TWO_SIDED_P_005 = 2.807033768343811


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a BRCA DEG 0/1 matrix using normal-sample Gaussian baseline."
    )
    parser.add_argument(
        "--expression",
        default="tcga_RSEM_Hugo_norm_count.gz",
        help="Gene-by-sample expression matrix, gzipped TSV. First column must be gene/sample.",
    )
    parser.add_argument(
        "--phenotype",
        default="TCGA_phenotype_denseDataOnlyDownload.tsv.gz",
        help="TCGA phenotype file with sample_type_id and _primary_disease.",
    )
    parser.add_argument(
        "--cancer",
        default="breast invasive carcinoma",
        help="Cancer name in phenotype _primary_disease.",
    )
    parser.add_argument(
        "--outdir",
        default="brca_deg_outputs",
        help="Output directory.",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=1000,
        help="Number of genes to process per chunk.",
    )
    parser.add_argument(
        "--z-cutoff",
        type=float,
        default=DEFAULT_Z_CUTOFF_FOR_TWO_SIDED_P_005,
        help="Absolute z-score cutoff. Default is two-sided Gaussian p <= 0.005.",
    )
    parser.add_argument(
        "--sd-floor",
        type=float,
        default=1e-8,
        help="Minimum normal-sample standard deviation used to avoid division by zero.",
    )
    parser.add_argument(
        "--force-log2p",
        action="store_true",
        help="Apply log2(x + 1) to expression values before DEG calling.",
    )
    return parser.parse_args()


def read_header(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return handle.readline().rstrip("\n").split("\t")


def sample_type_from_barcode(sample: str) -> str:
    parts = sample.split("-")
    if len(parts) < 4:
        return ""
    return parts[3][:2]


def select_brca_samples(phenotype_path: Path, expression_samples: set[str], cancer: str):
    pheno = pd.read_csv(phenotype_path, sep="\t", compression="infer", dtype=str)
    required = {"sample", "sample_type_id", "_primary_disease"}
    missing = required.difference(pheno.columns)
    if missing:
        raise ValueError(f"Phenotype file is missing required columns: {sorted(missing)}")

    pheno = pheno[pheno["sample"].isin(expression_samples)].copy()
    cancer_mask = pheno["_primary_disease"].str.lower().eq(cancer.lower())
    cancer_pheno = pheno[cancer_mask].copy()

    tumor = cancer_pheno[cancer_pheno["sample_type_id"].eq("01")]["sample"].tolist()
    normal = cancer_pheno[cancer_pheno["sample_type_id"].eq("11")]["sample"].tolist()

    if not tumor:
        raise ValueError(f"No primary tumor samples found for cancer={cancer!r}.")
    if len(normal) < 2:
        raise ValueError(
            f"Need at least 2 normal samples to estimate sd; found {len(normal)} for {cancer!r}."
        )

    return tumor, normal, cancer_pheno


def main() -> None:
    args = parse_args()
    expression_path = Path(args.expression)
    phenotype_path = Path(args.phenotype)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    header = read_header(expression_path)
    gene_col = header[0]
    expression_samples = header[1:]
    expression_sample_set = set(expression_samples)

    tumor_samples, normal_samples, cancer_pheno = select_brca_samples(
        phenotype_path, expression_sample_set, args.cancer
    )
    selected_samples = normal_samples + tumor_samples
    usecols = [gene_col] + selected_samples

    sample_table = cancer_pheno[
        cancer_pheno["sample"].isin(selected_samples)
    ][["sample", "sample_type_id", "sample_type", "_primary_disease"]].copy()
    sample_table.to_csv(outdir / "brca_selected_samples.csv", index=False)

    deg_chunks = []
    baseline_chunks = []
    gene_summaries = []

    reader = pd.read_csv(
        expression_path,
        sep="\t",
        compression="infer",
        usecols=usecols,
        chunksize=args.chunksize,
    )

    for chunk_index, chunk in enumerate(reader, start=1):
        genes = chunk[gene_col].astype(str)
        values = chunk[selected_samples].apply(pd.to_numeric, errors="coerce")
        if args.force_log2p:
            values = np.log2(values + 1.0)

        normal_values = values[normal_samples]
        tumor_values = values[tumor_samples]

        normal_mean = normal_values.mean(axis=1, skipna=True)
        normal_sd_raw = normal_values.std(axis=1, skipna=True, ddof=1)
        normal_sd = normal_sd_raw.mask(normal_sd_raw < args.sd_floor, args.sd_floor)

        z = tumor_values.sub(normal_mean, axis=0).div(normal_sd, axis=0)
        deg = z.abs().ge(args.z_cutoff).astype("int8")
        deg.insert(0, "gene", genes.values)

        deg_t = deg.set_index("gene").T
        deg_t.index.name = "sample"
        deg_chunks.append(deg_t)

        baseline_chunks.append(
            pd.DataFrame(
                {
                    "gene": genes.values,
                    "normal_n": normal_values.notna().sum(axis=1).astype(int).values,
                    "normal_mean": normal_mean.values,
                    "normal_sd": normal_sd_raw.values,
                    "normal_sd_used": normal_sd.values,
                }
            )
        )

        gene_summaries.append(
            pd.DataFrame(
                {
                    "gene": genes.values,
                    "deg_tumor_count": deg.drop(columns=["gene"]).sum(axis=1).astype(int).values,
                    "deg_tumor_fraction": deg.drop(columns=["gene"]).mean(axis=1).values,
                }
            )
        )

        print(f"Processed chunk {chunk_index}: {len(chunk)} genes")

    deg_matrix = pd.concat(deg_chunks, axis=1)
    baseline = pd.concat(baseline_chunks, axis=0, ignore_index=True)
    gene_summary = pd.concat(gene_summaries, axis=0, ignore_index=True)

    matrix_path = outdir / "brca_DEG_binary_matrix.csv.gz"
    baseline_path = outdir / "brca_normal_baseline_by_gene.csv"
    gene_summary_path = outdir / "brca_DEG_gene_summary.csv"
    run_summary_path = outdir / "brca_DEG_run_summary.txt"

    deg_matrix.to_csv(matrix_path, compression="gzip")
    baseline.to_csv(baseline_path, index=False)
    gene_summary.to_csv(gene_summary_path, index=False)

    total_entries = int(deg_matrix.shape[0] * deg_matrix.shape[1])
    total_deg = int(deg_matrix.to_numpy(dtype=np.int8).sum())
    summary_lines = [
        "BRCA DEG binary matrix construction",
        f"cancer: {args.cancer}",
        f"expression: {expression_path}",
        f"phenotype: {phenotype_path}",
        f"log2p_applied_by_script: {args.force_log2p}",
        f"z_cutoff: {args.z_cutoff}",
        "p_value_rule: two-sided Gaussian p <= 0.005",
        f"normal_samples: {len(normal_samples)}",
        f"tumor_samples: {len(tumor_samples)}",
        f"genes: {deg_matrix.shape[1]}",
        f"matrix_shape: {deg_matrix.shape[0]} tumor samples x {deg_matrix.shape[1]} genes",
        f"total_DEG_ones: {total_deg}",
        f"overall_DEG_fraction: {total_deg / total_entries:.6f}",
        f"matrix_file: {matrix_path}",
        f"normal_baseline_file: {baseline_path}",
        f"gene_summary_file: {gene_summary_path}",
    ]
    run_summary_path.write_text("\n".join(summary_lines) + "\n")

    print("\n".join(summary_lines))


if __name__ == "__main__":
    main()
