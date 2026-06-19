#!/usr/bin/env python3
"""Create debugging artifacts for the BRCA MPL IBI-DT pilot run."""

from __future__ import annotations

import argparse
import contextlib
import io
import math
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import IBIDT_core as core
import run_ibidT_on_brca as brca


def parse_args() -> argparse.Namespace:
    tcga_dir = REPO_ROOT / "TCGA_PanCancer_RNAseq"
    parser = argparse.ArgumentParser(description="Audit BRCA MPL IBI-DT output files.")
    parser.add_argument("--target-gene", default="MPL")
    parser.add_argument("--bootstraps", type=int, default=20)
    parser.add_argument("--outdir", default=str(HERE.parent / "brca_run_outputs_complete_sga"))
    parser.add_argument("--top-sga", type=int, default=500)
    parser.add_argument(
        "--deg-matrix",
        default=str(tcga_dir / "brca_deg_outputs" / "brca_DEG_binary_matrix.csv.gz"),
    )
    parser.add_argument(
        "--gistic",
        default=str(tcga_dir / "brca_gistic_thresholded_by_genes.gz"),
    )
    parser.add_argument(
        "--mutation",
        default=str(tcga_dir / "mc3.v0.2.8.PUBLIC.nonsilentGene.xena.gz"),
    )
    parser.add_argument(
        "--expression",
        default=str(tcga_dir / "tcga_RSEM_Hugo_norm_count.gz"),
    )
    parser.add_argument(
        "--selected-samples",
        default=str(tcga_dir / "brca_deg_outputs" / "brca_selected_samples.csv"),
    )
    parser.add_argument(
        "--baseline",
        default=str(tcga_dir / "brca_deg_outputs" / "brca_normal_baseline_by_gene.csv"),
    )
    parser.add_argument("--stability-seeds", default="101,202,303")
    parser.add_argument("--skip-stability", action="store_true")
    return parser.parse_args()


def beta_pairs_for_tree(tree_group: pd.DataFrame) -> set[str]:
    """Reproduce the pair construction in SingleNode_and_PathwayPair_analysis.py."""
    tree_temp = tree_group[tree_group["Node"].astype(str) != "0"]
    z = np.array(tree_temp)
    all_leaf_node = [i for i in np.arange(len(z)) if z[i][1] != "0"]
    pairs: set[str] = set()
    for line_number in all_leaf_node[1:]:
        current_gene = z[line_number][1]
        direction = z[line_number][3]
        next_search = z[line_number][2] - 1
        next_search_index = line_number
        numbers = list(np.arange(next_search_index))
        numbers.reverse()
        for i in numbers:
            if z[i][1] != "0" and z[i][2] == next_search:
                prev_gene = z[i][1]
                prev_gene_dir = prev_gene + "_" + str(direction)
                pairs.add(prev_gene_dir + "___" + current_gene)
                next_search = next_search - 1
                direction = z[i][3]
    return pairs


def write_alpha_beta_logs(tree_df: pd.DataFrame, audit_dir: Path, bootstraps: int) -> None:
    alpha_rows = []
    beta_rows = []
    alpha_records = []
    beta_records = []

    for tree_id in sorted(tree_df["Tree"].unique()):
        group = tree_df[tree_df["Tree"] == tree_id]
        genes = sorted(group[group["Node"].astype(str) != "0"]["Node"].astype(str).unique())
        pairs = sorted(beta_pairs_for_tree(group))

        alpha_rows.append(
            {"Tree": tree_id, "unique_gene_count": len(genes), "genes": "|".join(genes)}
        )
        beta_rows.append(
            {"Tree": tree_id, "unique_pair_count": len(pairs), "pairs": "|".join(pairs)}
        )
        alpha_records.extend({"Tree": tree_id, "SGA": gene, "contribution": 1} for gene in genes)
        beta_records.extend({"Tree": tree_id, "Gene_pair": pair, "contribution": 1} for pair in pairs)

    alpha_per_tree = pd.DataFrame(alpha_rows)
    beta_per_tree = pd.DataFrame(beta_rows)
    alpha_per_tree.to_csv(audit_dir / "alpha_per_tree_unique_genes.csv", index=False)
    beta_per_tree.to_csv(audit_dir / "beta_per_tree_unique_pairs.csv", index=False)

    tree_cols = [f"tree_{i}" for i in range(1, bootstraps + 1)]
    if alpha_records:
        alpha_long = pd.DataFrame(alpha_records)
        alpha_matrix = (
            alpha_long.pivot_table(index="SGA", columns="Tree", values="contribution", fill_value=0)
            .reindex(columns=list(range(1, bootstraps + 1)), fill_value=0)
            .astype(int)
        )
        alpha_matrix.columns = tree_cols
        alpha_matrix["alpha_total"] = alpha_matrix[tree_cols].sum(axis=1)
        alpha_matrix = alpha_matrix.sort_values("alpha_total", ascending=False)
    else:
        alpha_matrix = pd.DataFrame(columns=tree_cols + ["alpha_total"])
    alpha_matrix.to_csv(audit_dir / "alpha_gene_by_tree_matrix.csv")

    if beta_records:
        beta_long = pd.DataFrame(beta_records)
        beta_matrix = (
            beta_long.pivot_table(index="Gene_pair", columns="Tree", values="contribution", fill_value=0)
            .reindex(columns=list(range(1, bootstraps + 1)), fill_value=0)
            .astype(int)
        )
        beta_matrix.columns = tree_cols
        beta_matrix["beta_total"] = beta_matrix[tree_cols].sum(axis=1)
        beta_matrix = beta_matrix.sort_values("beta_total", ascending=False)
    else:
        beta_matrix = pd.DataFrame(columns=tree_cols + ["beta_total"])
    beta_matrix.to_csv(audit_dir / "beta_pair_by_tree_matrix.csv")

    checks = [
        "alpha_beta_invariant_checks",
        f"bootstraps: {bootstraps}",
        f"max_alpha_total: {int(alpha_matrix['alpha_total'].max()) if not alpha_matrix.empty else 0}",
        f"max_beta_total: {int(beta_matrix['beta_total'].max()) if not beta_matrix.empty else 0}",
        f"alpha_cell_max: {int(alpha_matrix[tree_cols].max().max()) if not alpha_matrix.empty else 0}",
        f"beta_cell_max: {int(beta_matrix[tree_cols].max().max()) if not beta_matrix.empty else 0}",
        f"alpha_total_le_bootstraps: {bool(alpha_matrix.empty or alpha_matrix['alpha_total'].max() <= bootstraps)}",
        f"beta_total_le_bootstraps: {bool(beta_matrix.empty or beta_matrix['beta_total'].max() <= bootstraps)}",
        f"alpha_each_tree_contributes_at_most_one: {bool(alpha_matrix.empty or alpha_matrix[tree_cols].max().max() <= 1)}",
        f"beta_each_tree_contributes_at_most_one: {bool(beta_matrix.empty or beta_matrix[tree_cols].max().max() <= 1)}",
    ]
    (audit_dir / "alpha_beta_invariant_checks.txt").write_text("\n".join(checks) + "\n")


def reconstruct_bootstrap_deg_counts(sga_deg: pd.DataFrame, bootstraps: int) -> pd.DataFrame:
    rows = []
    for tree_id in range(1, bootstraps + 1):
        boot = sga_deg.sample(frac=1, replace=True, random_state=tree_id)
        rows.append(
            {
                "Tree": tree_id,
                "bootstrap_DEG_0": int((boot["pheno"] == 0).sum()),
                "bootstrap_DEG_1": int((boot["pheno"] == 1).sum()),
            }
        )
    return pd.DataFrame(rows)


def write_tree_summaries(
    tree_df: pd.DataFrame,
    sga_deg: pd.DataFrame,
    audit_dir: Path,
    bootstraps: int,
    outdir: Path,
) -> None:
    class_counts_path = outdir / "bootstrap_class_counts.csv"
    if class_counts_path.exists():
        class_counts = pd.read_csv(class_counts_path)
        class_counts.to_csv(audit_dir / "bootstrap_class_counts.csv", index=False)
        deg_counts = class_counts.rename(
            columns={
                "tree_id": "Tree",
                "DEG_0_count": "bootstrap_DEG_0",
                "DEG_1_count": "bootstrap_DEG_1",
            }
        )
    else:
        deg_counts = reconstruct_bootstrap_deg_counts(sga_deg, bootstraps)
        deg_counts.rename(
            columns={
                "Tree": "tree_id",
                "bootstrap_DEG_0": "DEG_0_count",
                "bootstrap_DEG_1": "DEG_1_count",
            }
        ).assign(total_count=lambda df: df["DEG_0_count"] + df["DEG_1_count"]).to_csv(
            audit_dir / "bootstrap_class_counts.csv",
            index=False,
        )
    rows = []
    for tree_id, group in tree_df.groupby("Tree"):
        nonzero = group[group["Node"].astype(str) != "0"]
        root_rows = group[(group["Level"] == 1) & (group["Branch"] == 2)]
        rows.append(
            {
                "tree_id": tree_id,
                "num_rows_in_tree_output": int(group.shape[0]),
                "num_nonzero_nodes": int(nonzero.shape[0]),
                "unique_SGA_used": int(nonzero["Node"].astype(str).nunique()),
                "duplicated_nonzero_node_rows": int(
                    nonzero.shape[0] - nonzero["Node"].astype(str).nunique()
                ),
                "root_node": str(root_rows.iloc[0]["Node"]) if not root_rows.empty else "",
                "depth": int(group["Level"].max()) if not group.empty else 0,
                "terminal_zero_rows": int((group["Node"].astype(str) == "0").sum()),
            }
        )
    summary = pd.DataFrame(rows).merge(deg_counts, left_on="tree_id", right_on="Tree", how="left")
    summary = summary.drop(columns=["Tree"])
    summary.to_csv(audit_dir / "bootstrap_tree_summaries.csv", index=False)


def prepare_sga_deg(args: argparse.Namespace):
    deg = brca.load_deg_target(Path(args.deg_matrix), args.target_gene)
    scna_sga = brca.load_scna_sga(Path(args.gistic), deg.index.tolist())
    mutation_sga = brca.load_mutation_sga(Path(args.mutation), deg.index.tolist())
    common_samples = [
        s for s in deg.index.tolist() if s in set(scna_sga.columns) and s in set(mutation_sga.columns)
    ]
    deg = deg.loc[common_samples]
    sga_by_gene, _ = brca.combine_sga(scna_sga, mutation_sga, common_samples, args.top_sga)
    sga = sga_by_gene.T
    sga_deg = pd.concat([deg, sga], axis=1).fillna(0)
    sga_deg["pheno"] = sga_deg["pheno"].astype("int8")
    for col in sga_deg.columns[1:]:
        sga_deg[col] = sga_deg[col].astype("int8")
    return deg, sga_by_gene, sga_deg, scna_sga, mutation_sga, common_samples


def write_sample_alignment(
    audit_dir: Path, deg: pd.DataFrame, sga_by_gene: pd.DataFrame, sga_deg: pd.DataFrame
) -> None:
    x_samples = list(sga_by_gene.columns)
    y_samples = list(deg.index)
    x_set = set(x_samples)
    y_set = set(y_samples)
    diagnostics = [
        "sample_alignment_diagnostics",
        f"X_sample_count: {len(x_samples)}",
        f"Y_sample_count: {len(y_samples)}",
        f"intersection_count: {len(x_set & y_set)}",
        f"missing_in_X_count: {len(y_set - x_set)}",
        f"missing_in_Y_count: {len(x_set - y_set)}",
        f"sample_order_match: {x_samples == y_samples}",
        f"Y_NA_count: {int(deg.isna().sum().sum())}",
        f"X_NA_count: {int(sga_by_gene.isna().sum().sum())}",
        f"joined_matrix_NA_count: {int(sga_deg.isna().sum().sum())}",
        f"Y_non_binary_values: {sorted(set(deg['pheno'].dropna().unique()) - {0, 1})}",
        f"X_min: {int(sga_by_gene.min().min())}",
        f"X_max: {int(sga_by_gene.max().max())}",
    ]
    (audit_dir / "sample_alignment_diagnostics.txt").write_text("\n".join(diagnostics) + "\n")


def write_raw_sga_frequency(
    audit_dir: Path,
    scna_sga: pd.DataFrame,
    mutation_sga: pd.DataFrame,
    common_samples: list[str],
    selected_sga: pd.DataFrame,
) -> None:
    scna = scna_sga.loc[:, common_samples]
    mutation = mutation_sga.loc[:, common_samples]
    genes = scna.index.union(mutation.index)
    scna_counts = scna.reindex(index=genes, fill_value=0).sum(axis=1)
    mutation_counts = mutation.reindex(index=genes, fill_value=0).sum(axis=1)
    combined_counts = (
        scna.reindex(index=genes, fill_value=0).astype("int8")
        | mutation.reindex(index=genes, fill_value=0).astype("int8")
    ).sum(axis=1)
    n = len(common_samples)
    freq = pd.DataFrame(
        {
            "gene": genes,
            "scna_count": scna_counts.astype(int).values,
            "scna_frequency": (scna_counts / n).values,
            "mutation_count": mutation_counts.astype(int).values,
            "mutation_frequency": (mutation_counts / n).values,
            "combined_count": combined_counts.astype(int).values,
            "combined_frequency": (combined_counts / n).values,
            "selected_top500": [gene in set(selected_sga.index) for gene in genes],
        }
    )
    freq = freq.sort_values(["combined_count", "scna_count", "mutation_count"], ascending=False)
    freq.to_csv(audit_dir / "raw_sga_frequency_before_top500.csv", index=False)
    quantiles = freq[["scna_frequency", "mutation_frequency", "combined_frequency"]].quantile(
        [0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1.0]
    )
    quantiles.to_csv(audit_dir / "raw_sga_frequency_quantiles.csv")


def describe(values: pd.Series, prefix: str) -> dict[str, float]:
    values = pd.to_numeric(values, errors="coerce").dropna()
    return {
        f"{prefix}_n": int(values.shape[0]),
        f"{prefix}_mean": float(values.mean()),
        f"{prefix}_sd": float(values.std(ddof=1)),
        f"{prefix}_min": float(values.min()),
        f"{prefix}_q25": float(values.quantile(0.25)),
        f"{prefix}_median": float(values.median()),
        f"{prefix}_q75": float(values.quantile(0.75)),
        f"{prefix}_max": float(values.max()),
    }


def write_mpl_deg_audit(args: argparse.Namespace, audit_dir: Path, deg: pd.DataFrame) -> None:
    selected = pd.read_csv(args.selected_samples, dtype={"sample": str, "sample_type_id": str})
    normal_samples = selected[selected["sample_type_id"].astype(str).eq("11")]["sample"].tolist()
    tumor_samples = selected[selected["sample_type_id"].astype(str).eq("01")]["sample"].tolist()
    common_tumor_samples = [s for s in deg.index if s in set(tumor_samples)]

    header = pd.read_csv(args.expression, sep="\t", compression="infer", nrows=0).columns.tolist()
    gene_col = header[0]
    available = set(header[1:])
    normal_samples = [s for s in normal_samples if s in available]
    common_tumor_samples = [s for s in common_tumor_samples if s in available]
    expr = pd.read_csv(
        args.expression,
        sep="\t",
        compression="infer",
        usecols=[gene_col] + normal_samples + common_tumor_samples,
    )
    expr = expr[expr[gene_col].astype(str).eq(args.target_gene)]
    if expr.empty:
        raise ValueError(f"{args.target_gene} was not found in expression file.")
    values = pd.to_numeric(expr.iloc[0].drop(labels=[gene_col]), errors="coerce")

    baseline = pd.read_csv(args.baseline)
    base = baseline[baseline["gene"].astype(str).eq(args.target_gene)].iloc[0]
    normal_mean = float(base["normal_mean"])
    normal_sd_used = float(base["normal_sd_used"])
    z = (values.reindex(common_tumor_samples) - normal_mean) / normal_sd_used
    p = z.abs().map(lambda value: math.erfc(float(value) / math.sqrt(2.0)))
    gistic_header = pd.read_csv(args.gistic, sep="\t", compression="infer", nrows=0).columns.tolist()
    gistic_gene_col = gistic_header[0]
    gistic_samples = [s for s in common_tumor_samples if s in set(gistic_header[1:])]
    gistic = pd.read_csv(
        args.gistic,
        sep="\t",
        compression="infer",
        usecols=[gistic_gene_col] + gistic_samples,
    )
    gistic = gistic[gistic[gistic_gene_col].astype(str).eq(args.target_gene)]
    if gistic.empty:
        gistic_values = pd.Series(index=common_tumor_samples, dtype="float64")
    else:
        gistic_values = pd.to_numeric(gistic.iloc[0].drop(labels=[gistic_gene_col]), errors="coerce")

    audit = pd.DataFrame(
        {
            "sample": common_tumor_samples,
            "tumor_expression": values.reindex(common_tumor_samples).values,
            "normal_mean": normal_mean,
            "normal_sd_used": normal_sd_used,
            "z_score": z.values,
            "two_sided_gaussian_p": p.values,
            "DEG_label": deg.loc[common_tumor_samples, "pheno"].astype(int).values,
            "abs_z_ge_threshold": z.abs().ge(2.807033768343811).astype(int).values,
            "MPL_GISTIC": gistic_values.reindex(common_tumor_samples).values,
        }
    )
    audit["cis_SCNA_event"] = audit["MPL_GISTIC"].isin([-2, 2]).astype(int)
    audit["DEG_if_cis_SCNA_excluded"] = (
        audit["DEG_label"].eq(1) & audit["cis_SCNA_event"].eq(0)
    ).astype(int)
    audit.to_csv(audit_dir / "MPL_DEG_audit_per_sample.csv", index=False)
    pd.DataFrame(
        {"sample": normal_samples, "normal_expression": values.reindex(normal_samples).values}
    ).to_csv(audit_dir / "MPL_normal_expression_values.csv", index=False)

    summary = {
        "target_gene": args.target_gene,
        "threshold_rule": "abs(z) >= 2.807033768343811, equivalent to two-sided Gaussian p <= 0.005",
        "normal_mean_from_baseline": normal_mean,
        "normal_sd_used_from_baseline": normal_sd_used,
        "DEG_positive_common_tumors": int(audit["DEG_label"].sum()),
        "MPL_cis_SCNA_event_count": int(audit["cis_SCNA_event"].sum()),
        "DEG_positive_after_cis_SCNA_exclusion": int(audit["DEG_if_cis_SCNA_excluded"].sum()),
    }
    summary.update(describe(pd.Series(values.reindex(normal_samples).values), "normal_expression"))
    summary.update(describe(audit["tumor_expression"], "tumor_expression"))
    summary.update(describe(audit["z_score"], "z_score"))
    pd.Series(summary).to_csv(audit_dir / "MPL_DEG_audit_summary.csv", header=["value"])


def run_tree_search_with_seed(sga_deg: pd.DataFrame, bootstraps: int, seed: int) -> pd.DataFrame:
    all_rows = []
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    for qq in range(1, bootstraps + 1):
        random_state = seed * 1000 + qq
        boot = sga_deg.sample(frac=1, replace=True, random_state=random_state).copy()
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
            sga_pheno, sga_var, root_var, variants, q=1.0, device=device
        )
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
                q=1.0,
                device=device,
            )
        for row in rows:
            row.append(qq)
        all_rows.extend(rows)
    return pd.DataFrame(
        all_rows,
        columns=["Marginal", "Node", "Level", "Branch", "OR_Details", "Additional_Stat", "Tree"],
    )


def top_alpha_genes(tree_df: pd.DataFrame, n: int = 20) -> list[str]:
    records = []
    for tree_id in sorted(tree_df["Tree"].unique()):
        group = tree_df[(tree_df["Tree"] == tree_id) & (tree_df["Node"].astype(str) != "0")]
        records.extend(str(gene) for gene in group["Node"].astype(str).unique())
    counts = pd.Series(records).value_counts()
    return counts.head(n).index.tolist()


def write_stability(args: argparse.Namespace, audit_dir: Path, sga_deg: pd.DataFrame) -> None:
    seeds = [int(seed.strip()) for seed in args.stability_seeds.split(",") if seed.strip()]
    rows = []
    top_sets = {}
    start = time.time()
    for seed in seeds:
        tree_df = run_tree_search_with_seed(sga_deg, args.bootstraps, seed)
        genes = top_alpha_genes(tree_df, 20)
        top_sets[seed] = set(genes)
        rows.extend({"seed": seed, "rank": i + 1, "SGA": gene} for i, gene in enumerate(genes))
    pd.DataFrame(rows).to_csv(audit_dir / "bootstrap_stability_top20_by_seed.csv", index=False)

    overlap = pd.DataFrame(index=seeds, columns=seeds, dtype=int)
    for a in seeds:
        for b in seeds:
            overlap.loc[a, b] = len(top_sets[a] & top_sets[b])
    overlap.to_csv(audit_dir / "bootstrap_stability_top20_overlap_matrix.csv")
    (audit_dir / "bootstrap_stability_runtime.txt").write_text(
        f"seeds: {seeds}\nbootstraps_per_seed: {args.bootstraps}\nseconds: {time.time() - start:.2f}\n"
    )


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    audit_dir = outdir / "debug_artifacts"
    audit_dir.mkdir(parents=True, exist_ok=True)

    tree_df = pd.read_csv(outdir / "IBIDT_trees.csv", index_col=0)
    deg, selected_sga, sga_deg, scna_sga, mutation_sga, common_samples = prepare_sga_deg(args)

    write_alpha_beta_logs(tree_df, audit_dir, args.bootstraps)
    write_tree_summaries(tree_df, sga_deg, audit_dir, args.bootstraps, outdir)
    write_sample_alignment(audit_dir, deg, selected_sga, sga_deg)
    write_raw_sga_frequency(audit_dir, scna_sga, mutation_sga, common_samples, selected_sga)
    write_mpl_deg_audit(args, audit_dir, deg)
    if not args.skip_stability:
        write_stability(args, audit_dir, sga_deg)

    manifest = sorted(path.name for path in audit_dir.iterdir() if path.is_file())
    (audit_dir / "MANIFEST.txt").write_text("\n".join(manifest) + "\n")
    print(f"Wrote {len(manifest)} audit files to {audit_dir}")


if __name__ == "__main__":
    main()
