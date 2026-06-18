# IBI-DT
IBI-DT (Individualized Bayesian Inference - Decision Tree) algorithm advances the IBI methodology by integrating decision trees to examine complex traits at the subgroup level in cancer research. This approach offers a framework for exploring and understanding the interactions of genomic variants and their impacts on cancer.

# IBIDT_core
IBIDT_core provides the structure of 'n' decision trees. This included all the essential components required for an in-depth analysis of tree-based data. The key elements included are:


Marginal: Provides corresponding marginal.

Node: Represents the best node.

Level: Indicates the depth.

Branch: Details the diverging paths (direction).

OR_Details: Provides the odds ration and corresponding details.

Tree: Tree number

# SingleNode_and_PathwayPair_analysis

This module is an integral part of our framework, designed to utilize the tree information and provide key statistics.

α Statistic for Single Nodes: The file calculates the α statistic for single nodes and saves it in a file named SingleNode_analysis.csv. 

β Statistic for Pathway Pairs: In addition to single node analysis, this file also computes the β statistic for pathway pairs and saves it in PathwayPair_analysis.csv.

## BRCA MPL Pilot Run

This repository also contains a TCGA BRCA pilot run for one target DEG gene:

- Target DEG phenotype: `MPL`
- SGA definition: TCGA BRCA GISTIC gene-level copy-number alteration (`-2` or `2`) OR MC3 non-silent mutation
- Samples used after X/Y alignment: `779`
- SGA genes used: top `500` by combined alteration frequency
- Bootstrap trees: `20`
- Tree generation: original `codes/IBIDT_core.py`
- Single-node and pathway-pair statistics: original `codes/SingleNode_and_PathwayPair_analysis.py`

Main results are in `brca_run_outputs_complete_sga/`:

- `SingleNode_analysis_MPL.csv`
- `PathwayPair_analysis_MPL.csv`
- `IBIDT_trees_MPL.csv`
- `run_summary_MPL.txt`

Debugging and validation artifacts are in `brca_run_outputs_complete_sga/debug_artifacts/`:

- `alpha_gene_by_tree_matrix.csv`: gene-by-tree 0/1 alpha contribution matrix
- `beta_pair_by_tree_matrix.csv`: pair-by-tree 0/1 beta contribution matrix
- `alpha_beta_invariant_checks.txt`: verifies alpha/beta counts do not exceed bootstrap count
- `bootstrap_tree_summaries.csv`: one-row summary per bootstrap tree
- `sample_alignment_diagnostics.txt`: X/Y sample order, intersection, and NA checks
- `raw_sga_frequency_before_top500.csv`: SCNA/mutation/combined frequency before the top-500 filter
- `MPL_DEG_audit_per_sample.csv`: MPL expression, z-score, DEG label, and MPL cis-SCNA audit
- `bootstrap_stability_top20_by_seed.csv`: top-20 genes across 3 seed settings

Key invariant checks from the current run:

```text
bootstraps: 20
max_alpha_total: 20
max_beta_total: 12
alpha_cell_max: 1
beta_cell_max: 1
alpha_total_le_bootstraps: True
beta_total_le_bootstraps: True
alpha_each_tree_contributes_at_most_one: True
beta_each_tree_contributes_at_most_one: True
```

## What to Review First

For an external reviewer or a new ChatGPT session, start here:

1. `brca_run_outputs_complete_sga/run_summary_MPL.txt`
2. `brca_run_outputs_complete_sga/SingleNode_analysis_MPL.csv`
3. `brca_run_outputs_complete_sga/PathwayPair_analysis_MPL.csv`
4. `brca_run_outputs_complete_sga/debug_artifacts/alpha_beta_invariant_checks.txt`
5. `brca_run_outputs_complete_sga/debug_artifacts/alpha_gene_by_tree_matrix.csv`
6. `brca_run_outputs_complete_sga/debug_artifacts/beta_pair_by_tree_matrix.csv`
7. `brca_run_outputs_complete_sga/debug_artifacts/bootstrap_tree_summaries.csv`
8. `brca_run_outputs_complete_sga/debug_artifacts/sample_alignment_diagnostics.txt`
9. `brca_run_outputs_complete_sga/debug_artifacts/raw_sga_frequency_before_top500.csv`
10. `brca_run_outputs_complete_sga/debug_artifacts/MPL_DEG_audit_per_sample.csv`

Papers are in `papers/`.

TCGA input datasets are in `datasets/tcga_brca_inputs/`. The only omitted raw
input from git is `tcga_RSEM_Hugo_norm_count.gz` because it is about 892 MB; it
is intended to be attached as a GitHub Release asset.
