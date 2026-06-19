#!/usr/bin/env python3
"""Compare two MPL IBI-DT BRCA run directories."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare q=1.0 and q=0.75 MPL runs.")
    parser.add_argument("--q1-dir", default="brca_run_outputs_mpl_200_q1_stratified")
    parser.add_argument("--q075-dir", default="brca_run_outputs_mpl_200_q075_stratified")
    parser.add_argument("--out", default="mpl_q1_vs_q075_comparison.md")
    return parser.parse_args()


def read_checks(path: Path) -> dict[str, str]:
    checks = {}
    for line in (path / "debug_artifacts" / "alpha_beta_invariant_checks.txt").read_text().splitlines():
        if ": " in line:
            key, value = line.split(": ", 1)
            checks[key] = value
    return checks


def read_single(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path / "SingleNode_analysis_MPL.csv")
    return df.sort_values(
        ["Number of appearances as the best estimators", "SGA"],
        ascending=[False, True],
    )


def read_pairs(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path / "PathwayPair_analysis_MPL.csv")
    return df.sort_values(
        ["Number of appearances as partners in the same branch", "Gene_pair"],
        ascending=[False, True],
    )


def rank_lookup(single: pd.DataFrame, gene: str) -> str:
    matches = single[single["SGA"].astype(str).eq(gene)]
    if matches.empty:
        return "not in table"
    row = matches.iloc[0]
    return f"rank {int(row['Rank of number of appearances as the best estimators'])}, alpha {int(row['Number of appearances as the best estimators'])}"


def summarize_run(label: str, path: Path) -> dict[str, object]:
    single = read_single(path)
    pairs = read_pairs(path)
    trees = pd.read_csv(path / "debug_artifacts" / "bootstrap_tree_summaries.csv")
    checks = read_checks(path)
    roots = trees["root_node"].value_counts().reset_index()
    roots.columns = ["root_node", "count"]
    return {
        "label": label,
        "path": path,
        "single": single,
        "pairs": pairs,
        "trees": trees,
        "checks": checks,
        "top20_genes": single.head(20)["SGA"].astype(str).tolist(),
        "top20_pairs": pairs.head(20)["Gene_pair"].astype(str).tolist(),
        "avg_depth": float(trees["depth"].mean()),
        "max_depth": int(trees["depth"].max()),
        "avg_unique_sga": float(trees["unique_SGA_used"].mean()),
        "root_distribution": roots,
        "max_alpha": int(checks.get("max_alpha_total", 0)),
        "max_beta": int(checks.get("max_beta_total", 0)),
        "invariants_ok": all(
            checks.get(key) == "True"
            for key in [
                "alpha_total_le_bootstraps",
                "beta_total_le_bootstraps",
                "alpha_each_tree_contributes_at_most_one",
                "beta_each_tree_contributes_at_most_one",
            ]
        ),
    }


def markdown_table(df: pd.DataFrame, columns: list[str], n: int | None = None) -> str:
    if n is not None:
        df = df.head(n)
    df = df[columns].copy()
    rows = []
    rows.append("| " + " | ".join(columns) + " |")
    rows.append("| " + " | ".join("---" for _ in columns) + " |")
    for _, row in df.iterrows():
        values = [str(row[col]) for col in columns]
        rows.append("| " + " | ".join(values) + " |")
    return "\n".join(rows)


def main() -> None:
    args = parse_args()
    q1 = summarize_run("q=1.0", Path(args.q1_dir))
    q075 = summarize_run("q=0.75", Path(args.q075_dir))

    gene_overlap = sorted(set(q1["top20_genes"]) & set(q075["top20_genes"]))
    pair_overlap = sorted(set(q1["top20_pairs"]) & set(q075["top20_pairs"]))
    brca_genes = ["PIK3CA", "TP53", "GATA3", "CDH1", "ERBB2"]
    passenger_genes = ["TTN", "CSMD1", "CSMD3"]

    lines = [
        "# MPL q=1.0 vs q=0.75 Stratified Bootstrap Comparison",
        "",
        "Both runs use original `IBIDT_core.py` for tree generation and original `SingleNode_and_PathwayPair_analysis.py` for alpha/beta statistics.",
        "",
        "## Run Directories",
        "",
        f"- q=1.0: `{q1['path']}`",
        f"- q=0.75: `{q075['path']}`",
        "",
        "## Invariant Summary",
        "",
        "| metric | q=1.0 | q=0.75 |",
        "|---|---:|---:|",
        f"| max alpha | {q1['max_alpha']} | {q075['max_alpha']} |",
        f"| max beta | {q1['max_beta']} | {q075['max_beta']} |",
        f"| average tree depth | {q1['avg_depth']:.2f} | {q075['avg_depth']:.2f} |",
        f"| maximum tree depth | {q1['max_depth']} | {q075['max_depth']} |",
        f"| average unique SGA per tree | {q1['avg_unique_sga']:.2f} | {q075['avg_unique_sga']:.2f} |",
        f"| alpha/beta invariants OK | {q1['invariants_ok']} | {q075['invariants_ok']} |",
        "",
        f"Top-20 alpha gene overlap: {len(gene_overlap)} / 20",
        "",
        ", ".join(gene_overlap),
        "",
        f"Top-20 beta pair overlap: {len(pair_overlap)} / 20",
        "",
        "## Top 20 Alpha Genes: q=1.0",
        "",
        markdown_table(
            q1["single"],
            ["SGA", "Number of appearances as the best estimators", "Rank of number of appearances as the best estimators"],
            20,
        ),
        "",
        "## Top 20 Alpha Genes: q=0.75",
        "",
        markdown_table(
            q075["single"],
            ["SGA", "Number of appearances as the best estimators", "Rank of number of appearances as the best estimators"],
            20,
        ),
        "",
        "## Top 20 Beta Pairs: q=1.0",
        "",
        markdown_table(
            q1["pairs"],
            ["Gene_pair", "Number of appearances as partners in the same branch", "Rank of number of appearances as partners in the same branch", "Branch", "Start_gene", "Finish_gene"],
            20,
        ),
        "",
        "## Top 20 Beta Pairs: q=0.75",
        "",
        markdown_table(
            q075["pairs"],
            ["Gene_pair", "Number of appearances as partners in the same branch", "Rank of number of appearances as partners in the same branch", "Branch", "Start_gene", "Finish_gene"],
            20,
        ),
        "",
        "## Root Node Distribution: q=1.0",
        "",
        markdown_table(q1["root_distribution"], ["root_node", "count"], 25),
        "",
        "## Root Node Distribution: q=0.75",
        "",
        markdown_table(q075["root_distribution"], ["root_node", "count"], 25),
        "",
        "## BRCA Driver Stability",
        "",
        "| gene | q=1.0 | q=0.75 |",
        "|---|---|---|",
    ]
    for gene in brca_genes:
        lines.append(f"| {gene} | {rank_lookup(q1['single'], gene)} | {rank_lookup(q075['single'], gene)} |")

    lines.extend(
        [
            "",
            "## Passenger / Frequency-Sensitive Genes",
            "",
            "| gene | q=1.0 | q=0.75 |",
            "|---|---|---|",
        ]
    )
    for gene in passenger_genes:
        lines.append(f"| {gene} | {rank_lookup(q1['single'], gene)} | {rank_lookup(q075['single'], gene)} |")

    lines.extend(
        [
            "",
            "## Preliminary Interpretation",
            "",
            f"- q=0.75 has shallower trees than q=1.0: {q075['avg_depth']:.2f} vs {q1['avg_depth']:.2f} average depth.",
            f"- q=0.75 uses fewer unique SGAs per tree than q=1.0: {q075['avg_unique_sga']:.2f} vs {q1['avg_unique_sga']:.2f}.",
            f"- Top-20 alpha gene overlap is {len(gene_overlap)} / 20.",
            f"- q=1.0 invariants OK: {q1['invariants_ok']}; q=0.75 invariants OK: {q075['invariants_ok']}.",
        ]
    )

    Path(args.out).write_text("\n".join(lines) + "\n")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
