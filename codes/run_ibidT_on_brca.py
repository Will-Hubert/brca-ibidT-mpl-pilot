#!/usr/bin/env python3
"""
Run the original IBI-DT tree search on the BRCA DEG matrix.

The public demo implementation in IBIDT_core.py runs one binary DEG phenotype
at a time. This wrapper prepares matched BRCA DEG and SGA inputs, then calls the
same functions from IBIDT_core.py.

Current SGA source:
  - TCGA BRCA GISTIC thresholded gene-level copy-number calls.
  - MC3 gene-level non-silent mutation calls from UCSC Xena.
  - SGA = 1 when GISTIC call is -2/2 or non-silent mutation is 1.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import math
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import torch

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import IBIDT_core as core


def parse_args() -> argparse.Namespace:
    repo_root = HERE.parent.parent
    tcga_dir = repo_root / "TCGA_PanCancer_RNAseq"
    parser = argparse.ArgumentParser(
        description="Prepare BRCA inputs and run IBI-DT for one DEG target gene."
    )
    parser.add_argument(
        "--deg-matrix",
        default=str(tcga_dir / "brca_deg_outputs" / "brca_DEG_binary_matrix.csv.gz"),
        help="Tumor-by-gene BRCA DEG binary matrix from 04_make_brca_deg_binary_matrix.py.",
    )
    parser.add_argument(
        "--gistic",
        default=str(tcga_dir / "brca_gistic_thresholded_by_genes.gz"),
        help="Gene-by-sample GISTIC thresholded copy-number file.",
    )
    parser.add_argument(
        "--mutation",
        default=str(tcga_dir / "mc3.v0.2.8.PUBLIC.nonsilentGene.xena.gz"),
        help="Gene-by-sample MC3 non-silent mutation 0/1 file.",
    )
    parser.add_argument(
        "--target-gene",
        default="MPL",
        help="One DEG gene to use as the binary phenotype.",
    )
    parser.add_argument(
        "--top-sga",
        type=int,
        default=500,
        help="Keep the top N most frequent combined SGA genes for the tree search.",
    )
    parser.add_argument(
        "--bootstraps",
        type=int,
        default=20,
        help="Number of bootstrap trees. Original demo uses 200.",
    )
    parser.add_argument(
        "--q",
        type=float,
        default=1.0,
        help="IBI-DT structure prior parameter. Original demo uses 1.",
    )
    parser.add_argument(
        "--outdir",
        default=str(HERE.parent / "brca_run_outputs_complete_sga"),
        help="Output directory.",
    )
    parser.add_argument(
        "--scna-only",
        action="store_true",
        help="Use only GISTIC copy-number calls and skip mutation calls.",
    )
    parser.add_argument(
        "--verbose-tree",
        action="store_true",
        help="Print the original per-node logs from the IBI-DT dt() function.",
    )
    return parser.parse_args()


def read_gene_list_from_deg(deg_path: Path) -> list[str]:
    return pd.read_csv(deg_path, nrows=0).columns.tolist()[1:]


def load_deg_target(deg_path: Path, target_gene: str) -> pd.DataFrame:
    header = read_gene_list_from_deg(deg_path)
    if target_gene not in header:
        raise ValueError(f"Target gene {target_gene!r} was not found in DEG matrix.")
    deg = pd.read_csv(deg_path, usecols=["sample", target_gene])
    deg = deg.rename(columns={target_gene: "pheno"})
    deg["sample"] = deg["sample"].astype(str)
    deg["pheno"] = pd.to_numeric(deg["pheno"], errors="coerce").fillna(0).astype("int8")
    return deg.set_index("sample")


def load_scna_sga(gistic_path: Path, samples: list[str]) -> pd.DataFrame:
    sample_set = set(samples)
    header = pd.read_csv(gistic_path, sep="\t", compression="infer", nrows=0).columns.tolist()
    gene_col = header[0]
    available_samples = [s for s in header[1:] if s in sample_set]
    if not available_samples:
        raise ValueError("No overlapping samples between DEG matrix and GISTIC file.")

    usecols = [gene_col] + available_samples
    gistic = pd.read_csv(gistic_path, sep="\t", compression="infer", usecols=usecols)
    gistic = gistic.rename(columns={gene_col: "gene"})
    gistic = gistic.drop_duplicates(subset=["gene"], keep="first")
    gistic = gistic.set_index("gene")
    gistic = gistic.apply(pd.to_numeric, errors="coerce").fillna(0).astype("int8")

    sga = gistic.isin([-2, 2]).astype("int8")
    return sga.loc[:, available_samples]


def load_mutation_sga(mutation_path: Path, samples: list[str]) -> pd.DataFrame:
    sample_set = set(samples)
    header = pd.read_csv(mutation_path, sep="\t", compression="infer", nrows=0).columns.tolist()
    gene_col = header[0]
    available_samples = [s for s in header[1:] if s in sample_set]
    if not available_samples:
        raise ValueError("No overlapping samples between DEG matrix and mutation file.")

    usecols = [gene_col] + available_samples
    mutation = pd.read_csv(mutation_path, sep="\t", compression="infer", usecols=usecols)
    mutation = mutation.rename(columns={gene_col: "gene"})
    mutation["gene"] = mutation["gene"].astype(str)
    mutation = mutation.set_index("gene")
    mutation = mutation.apply(pd.to_numeric, errors="coerce").fillna(0)
    mutation = mutation.ne(0).astype("int8")
    if mutation.index.has_duplicates:
        mutation = mutation.groupby(level=0).max()
    return mutation.loc[:, available_samples]


def combine_sga(
    scna_sga: pd.DataFrame,
    mutation_sga: pd.DataFrame | None,
    samples: list[str],
    top_sga: int,
) -> tuple[pd.DataFrame, dict[str, int]]:
    scna = scna_sga.loc[:, samples]
    if mutation_sga is None:
        combined = scna.copy()
        mutation_gene_count = 0
    else:
        mutation = mutation_sga.loc[:, samples]
        genes = scna.index.union(mutation.index)
        combined = (
            scna.reindex(index=genes, fill_value=0).astype("int8")
            | mutation.reindex(index=genes, fill_value=0).astype("int8")
        ).astype("int8")
        mutation_gene_count = int((mutation.sum(axis=1) > 0).sum())

    counts = combined.sum(axis=1).sort_values(ascending=False)
    counts = counts[counts > 0]
    if top_sga:
        counts = counts.head(top_sga)
    combined = combined.loc[counts.index, samples]
    stats = {
        "scna_genes_with_any_event": int((scna.sum(axis=1) > 0).sum()),
        "mutation_genes_with_any_event": mutation_gene_count,
        "combined_genes_with_any_event": int((combined.sum(axis=1) > 0).sum()),
        "combined_sga_ones_in_selected_matrix": int(combined.to_numpy(dtype=np.int8).sum()),
    }
    return combined, stats


def run_one_target(sga_deg: pd.DataFrame, bootstraps: int, q: float, verbose_tree: bool) -> pd.DataFrame:
    all_rows = []
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    clock_start = datetime.now()

    for qq in range(1, bootstraps + 1):
        start_time = time.time()
        boot = sga_deg.sample(frac=1, replace=True, random_state=qq).copy()
        boot.index = list(range(1, boot.shape[0] + 1))
        boot.insert(loc=1, column="A0", value=1)

        traits = torch.as_tensor(boot[["pheno"]].values.copy(), dtype=torch.float64, device=device)
        variants = torch.as_tensor(boot.iloc[:, 1:].values.T.copy(), dtype=torch.float64, device=device)
        var_ids = list(boot.columns[1:])

        core.traits = traits
        core.variants = variants
        core.varIDs = var_ids

        _, _, _, top_index = core.Finding_RootNode(traits, variants)
        root_var = var_ids[top_index[0]]

        sga_pheno = torch.as_tensor(boot["pheno"].values.copy(), dtype=torch.float64, device=device)
        sga_var = torch.as_tensor(boot[root_var].values.copy(), dtype=torch.float64, device=device)
        root_lgm, or_details = core.root_node_function(
            sga_pheno, sga_var, root_var, variants, q=q, device=device
        )

        if verbose_tree:
            rows = core.dt(
                root_var,
                -math.inf,
                0,
                2,
                [],
                root_lgm,
                or_details,
                collect=[],
                pointer=[],
                q=q,
                device=device,
            )
        else:
            with contextlib.redirect_stdout(io.StringIO()):
                rows = core.dt(
                    root_var,
                    -math.inf,
                    0,
                    2,
                    [],
                    root_lgm,
                    or_details,
                    collect=[],
                    pointer=[],
                    q=q,
                    device=device,
                )
        for row in rows:
            row.append(qq)
        all_rows.extend(rows)
        print(f"bootstrap {qq}/{bootstraps}: root={root_var}, seconds={time.time() - start_time:.2f}")

    tree_df = pd.DataFrame(
        all_rows,
        columns=["Marginal", "Node", "Level", "Branch", "OR_Details", "Additional_Stat", "Tree"],
    )
    print(f"total elapsed seconds: {(datetime.now() - clock_start).seconds}")
    return tree_df


def run_original_single_node_and_pathway_pair_analysis(outdir: Path) -> None:
    env = os.environ.copy()
    env.setdefault("JOBLIB_MULTIPROCESSING", "0")
    subprocess.run(
        [sys.executable, str(HERE / "SingleNode_and_PathwayPair_analysis.py")],
        cwd=outdir,
        env=env,
        check=True,
    )


def main() -> None:
    args = parse_args()
    deg_path = Path(args.deg_matrix)
    gistic_path = Path(args.gistic)
    mutation_path = Path(args.mutation)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    deg = load_deg_target(deg_path, args.target_gene)
    scna_sga = load_scna_sga(gistic_path, deg.index.tolist())
    mutation_sga = None if args.scna_only else load_mutation_sga(mutation_path, deg.index.tolist())

    common_sample_set = set(scna_sga.columns)
    if mutation_sga is not None:
        common_sample_set &= set(mutation_sga.columns)
    common_samples = [s for s in deg.index.tolist() if s in common_sample_set]
    if not common_samples:
        raise ValueError("No common samples after loading DEG and SGA sources.")

    deg = deg.loc[common_samples]
    sga_by_gene, sga_stats = combine_sga(scna_sga, mutation_sga, common_samples, args.top_sga)
    sga = sga_by_gene.T
    sga_deg = pd.concat([deg, sga], axis=1).fillna(0)
    sga_deg["pheno"] = sga_deg["pheno"].astype("int8")
    for col in sga_deg.columns[1:]:
        sga_deg[col] = sga_deg[col].astype("int8")

    original_deg_path = outdir / "DEG.csv"
    original_sga_path = outdir / "SGA.csv"
    original_trees_path = outdir / "IBIDT_trees.csv"
    original_single_path = outdir / "SingleNode_analysis.csv"
    original_pair_path = outdir / "PathwayPair_analysis.csv"
    input_deg_path = outdir / f"DEG_{args.target_gene}.csv"
    input_sga_path = outdir / ("SGA_SCNA_only.csv" if args.scna_only else "SGA_complete_SCNA_or_mutation.csv")
    trees_path = outdir / f"IBIDT_trees_{args.target_gene}.csv"
    single_path = outdir / f"SingleNode_analysis_{args.target_gene}.csv"
    pair_path = outdir / f"PathwayPair_analysis_{args.target_gene}.csv"
    summary_path = outdir / f"run_summary_{args.target_gene}.txt"

    deg.to_csv(original_deg_path)
    sga_by_gene.to_csv(original_sga_path)
    shutil.copyfile(original_deg_path, input_deg_path)
    shutil.copyfile(original_sga_path, input_sga_path)

    print(
        f"Running IBI-DT target={args.target_gene}, samples={len(common_samples)}, "
        f"SGA genes={sga.shape[1]}, DEG positives={int(deg['pheno'].sum())}, "
        f"bootstraps={args.bootstraps}"
    )
    tree_df = run_one_target(sga_deg, args.bootstraps, args.q, args.verbose_tree)
    tree_df.to_csv(original_trees_path)
    shutil.copyfile(original_trees_path, trees_path)

    run_original_single_node_and_pathway_pair_analysis(outdir)
    shutil.copyfile(original_single_path, single_path)
    shutil.copyfile(original_pair_path, pair_path)

    summary = [
        "BRCA IBI-DT run",
        f"target_gene: {args.target_gene}",
        "IBI_DT_algorithm_source: original IBIDT_core.py",
        "single_node_and_pathway_pair_source: original SingleNode_and_PathwayPair_analysis.py",
        "SGA_source: TCGA BRCA GISTIC thresholded copy-number"
        + (" only" if args.scna_only else " + MC3 non-silent gene-level mutation"),
        "SGA_rule: GISTIC -2 or 2 -> 1; MC3 non-silent mutation 1 -> 1; combined SGA = SCNA OR mutation",
        f"mutation_events_included: {not args.scna_only}",
        f"gistic_file: {gistic_path}",
        f"mutation_file: {'not used' if args.scna_only else mutation_path}",
        f"samples_used: {len(common_samples)}",
        f"DEG_samples_available: {load_deg_target(deg_path, args.target_gene).shape[0]}",
        f"SCNA_samples_overlap_DEG: {len(scna_sga.columns)}",
        f"mutation_samples_overlap_DEG: {'not used' if mutation_sga is None else len(mutation_sga.columns)}",
        f"SGA_genes_used_top_by_frequency: {sga.shape[1]}",
        f"SCNA_genes_with_any_event_before_top_filter: {sga_stats['scna_genes_with_any_event']}",
        f"mutation_genes_with_any_event_before_top_filter: {sga_stats['mutation_genes_with_any_event']}",
        f"combined_SGA_genes_with_any_event_after_top_filter: {sga_stats['combined_genes_with_any_event']}",
        f"combined_SGA_ones_in_selected_matrix: {sga_stats['combined_sga_ones_in_selected_matrix']}",
        f"DEG_positive_samples: {int(deg['pheno'].sum())}",
        f"DEG_positive_fraction: {float(deg['pheno'].mean()):.6f}",
        f"bootstraps: {args.bootstraps}",
        f"original_DEG_file: {original_deg_path}",
        f"original_SGA_file: {original_sga_path}",
        f"original_trees_file: {original_trees_path}",
        f"original_single_node_file: {original_single_path}",
        f"original_pathway_pair_file: {original_pair_path}",
        f"trees_file: {trees_path}",
        f"single_node_file: {single_path}",
        f"pathway_pair_file: {pair_path}",
        f"prepared_DEG_file: {input_deg_path}",
        f"prepared_SGA_file: {input_sga_path}",
    ]
    summary_path.write_text("\n".join(summary) + "\n")
    print("\n".join(summary))


if __name__ == "__main__":
    main()
