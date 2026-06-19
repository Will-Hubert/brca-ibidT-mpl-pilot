# Method Run Completion

## Canonical Configuration

- Target gene: `MPL`
- Samples: `779`
- SGA genes: top `500`
- Bootstraps: `200`
- q: `0.75`
- Seed: `101`
- Stratified bootstrap: `True`
- DEG 0 per tree: `553`
- DEG 1 per tree: `226`
- Tree generation: original `codes/IBIDT_core.py`
- Alpha/beta analysis: original `codes/SingleNode_and_PathwayPair_analysis.py`
- Output directory: `brca_run_outputs_mpl_200_q075_stratified/`

This is a completed method run. It confirms that the pipeline executes and that
the required technical invariants hold. It is not a final biological conclusion.

## Completion Checks

| Check | Status |
|---|---|
| 200 trees generated | PASS |
| stratified class counts | PASS |
| X/Y sample alignment | PASS |
| no NA values | PASS |
| binary X/Y inputs | PASS |
| alpha <= 200 | PASS |
| beta <= 200 | PASS |
| each tree contributes at most one count per gene | PASS |
| each tree contributes at most one count per pair | PASS |
| original IBIDT_core.py used | PASS |
| original SingleNode_and_PathwayPair_analysis.py used | PASS |

## Evidence Files

- `brca_run_outputs_mpl_200_q075_stratified/run_summary_MPL.txt`
- `brca_run_outputs_mpl_200_q075_stratified/bootstrap_class_counts.csv`
- `brca_run_outputs_mpl_200_q075_stratified/debug_artifacts/alpha_beta_invariant_checks.txt`
- `brca_run_outputs_mpl_200_q075_stratified/debug_artifacts/sample_alignment_diagnostics.txt`
- `brca_run_outputs_mpl_200_q075_stratified/debug_artifacts/bootstrap_tree_summaries.csv`
- `brca_run_outputs_mpl_200_q075_stratified/debug_artifacts/alpha_gene_by_tree_matrix.csv`
- `brca_run_outputs_mpl_200_q075_stratified/debug_artifacts/beta_pair_by_tree_matrix.csv`

Key invariant results:

```text
bootstraps: 200
max_alpha_total: 98
max_beta_total: 23
alpha_cell_max: 1
beta_cell_max: 1
alpha_total_le_bootstraps: True
beta_total_le_bootstraps: True
alpha_each_tree_contributes_at_most_one: True
beta_each_tree_contributes_at_most_one: True
```

Stratified bootstrap class counts:

```text
DEG_0_count = 553 for every tree
DEG_1_count = 226 for every tree
total_count = 779 for every tree
```

## Re-run Command

From the repository root:

```bash
.venv_brca/bin/python codes/run_ibidT_on_brca.py \
  --target-gene MPL \
  --bootstraps 200 \
  --q 0.75 \
  --stratified-bootstrap \
  --seed 101 \
  --top-sga 500 \
  --outdir brca_run_outputs_mpl_200_q075_stratified

.venv_brca/bin/python codes/audit_brca_mpl_run.py \
  --outdir brca_run_outputs_mpl_200_q075_stratified \
  --bootstraps 200 \
  --top-sga 500 \
  --skip-stability
```

## Software Versions

```text
Python: 3.12.13
pandas: 3.0.3
numpy: 2.4.6
torch: 2.12.0
scipy: 1.17.1
joblib: 1.5.3
```

## Scope Notes

- No permutation test was run.
- No additional target genes were run.
- No algorithm changes were made to `codes/IBIDT_core.py`.
- No alpha/beta algorithm changes were made to `codes/SingleNode_and_PathwayPair_analysis.py`.
- `q=1.0` is retained only as a parameter contrast in `brca_run_outputs_mpl_200_q1_stratified/`.
