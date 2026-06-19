# MPL q=1.0 vs q=0.75 Stratified Bootstrap Comparison

Both runs use original `IBIDT_core.py` for tree generation and original `SingleNode_and_PathwayPair_analysis.py` for alpha/beta statistics.

## Run Directories

- q=1.0: `brca_run_outputs_mpl_200_q1_stratified`
- q=0.75: `brca_run_outputs_mpl_200_q075_stratified`

## Invariant Summary

| metric | q=1.0 | q=0.75 |
|---|---:|---:|
| max alpha | 196 | 98 |
| max beta | 96 | 23 |
| average tree depth | 31.64 | 6.64 |
| maximum tree depth | 52 | 10 |
| average unique SGA per tree | 63.27 | 13.85 |
| alpha/beta invariants OK | True | True |

Top-20 alpha gene overlap: 13 / 20

AHCTF1, CDH1, CSMD1, ERBB2, GATA3, MYEOV, NOTCH2, OR2T2, PIK3CA, RYR2, TP53, TTN, ZFPM2

Top-20 beta pair overlap: 2 / 20

## Top 20 Alpha Genes: q=1.0

| SGA | Number of appearances as the best estimators | Rank of number of appearances as the best estimators |
| --- | --- | --- |
| PIK3CA | 196 | 1.0 |
| TP53 | 194 | 2.0 |
| GATA3 | 187 | 3.0 |
| TTN | 183 | 4.0 |
| CDH1 | 180 | 5.0 |
| MYEOV | 164 | 6.0 |
| CSMD1 | 163 | 8.0 |
| NOTCH2 | 163 | 8.0 |
| OR2T2 | 160 | 9.0 |
| ERBB2 | 151 | 10.0 |
| CSMD3 | 146 | 11.0 |
| RYR2 | 143 | 12.0 |
| FLG | 140 | 13.0 |
| TCHH | 139 | 14.0 |
| RIMS2 | 136 | 15.0 |
| AHCTF1 | 132 | 16.0 |
| HMCN1 | 131 | 17.0 |
| TRPS1 | 130 | 18.0 |
| OBSCN | 129 | 19.0 |
| ZFPM2 | 128 | 20.0 |

## Top 20 Alpha Genes: q=0.75

| SGA | Number of appearances as the best estimators | Rank of number of appearances as the best estimators |
| --- | --- | --- |
| MYEOV | 98 | 1.0 |
| TTN | 83 | 2.0 |
| GATA3 | 76 | 3.0 |
| TP53 | 73 | 4.0 |
| CSMD1 | 61 | 5.0 |
| PIK3CA | 59 | 6.0 |
| CDC42BPA | 57 | 7.0 |
| ZFPM2 | 56 | 8.0 |
| NOTCH2 | 50 | 9.0 |
| PLEKHA6 | 48 | 10.0 |
| PPFIA1 | 47 | 11.0 |
| DNAH14 | 43 | 12.0 |
| TRIM67 | 42 | 13.0 |
| USH2A | 41 | 14.0 |
| CDH1 | 38 | 16.0 |
| ERBB2 | 38 | 16.0 |
| SPTA1 | 36 | 17.0 |
| AHCTF1 | 34 | 18.0 |
| RYR2 | 34 | 18.0 |
| OR2T2 | 32 | 20.0 |

## Top 20 Beta Pairs: q=1.0

| Gene_pair | Number of appearances as partners in the same branch | Rank of number of appearances as partners in the same branch | Branch | Start_gene | Finish_gene |
| --- | --- | --- | --- | --- | --- |
| OR2T2_0___PIK3CA | 96 | 1.0 | 0 | OR2T2 | PIK3CA |
| MYEOV_0___PIK3CA | 86 | 2.0 | 0 | MYEOV | PIK3CA |
| CSMD1_0___PIK3CA | 85 | 4.0 | 0 | CSMD1 | PIK3CA |
| TP53_0___PIK3CA | 85 | 4.0 | 0 | TP53 | PIK3CA |
| ZFPM2_0___PIK3CA | 85 | 4.0 | 0 | ZFPM2 | PIK3CA |
| CSMD1_0___TP53 | 83 | 6.0 | 0 | CSMD1 | TP53 |
| MYEOV_0___TP53 | 82 | 8.0 | 0 | MYEOV | TP53 |
| TTN_0___PIK3CA | 82 | 8.0 | 0 | TTN | PIK3CA |
| GATA3_0___PIK3CA | 81 | 10.0 | 0 | GATA3 | PIK3CA |
| OR2T2_0___TP53 | 81 | 10.0 | 0 | OR2T2 | TP53 |
| MYEOV_0___CDH1 | 77 | 12.0 | 0 | MYEOV | CDH1 |
| OR2T2_0___GATA3 | 77 | 12.0 | 0 | OR2T2 | GATA3 |
| TTN_0___TP53 | 77 | 12.0 | 0 | TTN | TP53 |
| OR2T2_0___OR2T3 | 75 | 14.0 | 0 | OR2T2 | OR2T3 |
| CSMD1_0___GATA3 | 73 | 16.0 | 0 | CSMD1 | GATA3 |
| MYEOV_0___TTN | 73 | 16.0 | 0 | MYEOV | TTN |
| TTN_0___GATA3 | 73 | 16.0 | 0 | TTN | GATA3 |
| ZFPM2_0___CDH1 | 73 | 16.0 | 0 | ZFPM2 | CDH1 |
| MYEOV_0___GATA3 | 71 | 20.0 | 0 | MYEOV | GATA3 |
| TTN_0___CDH1 | 71 | 20.0 | 0 | TTN | CDH1 |

## Top 20 Beta Pairs: q=0.75

| Gene_pair | Number of appearances as partners in the same branch | Rank of number of appearances as partners in the same branch | Branch | Start_gene | Finish_gene |
| --- | --- | --- | --- | --- | --- |
| MYEOV_0___TTN | 23 | 1.0 | 0 | MYEOV | TTN |
| MYEOV_1___PPFIA1 | 19 | 2.0 | 1 | MYEOV | PPFIA1 |
| MYEOV_0___CSMD1 | 15 | 4.0 | 0 | MYEOV | CSMD1 |
| MYEOV_0___PPFIA1 | 15 | 4.0 | 0 | MYEOV | PPFIA1 |
| MYEOV_1___SPTA1 | 15 | 4.0 | 1 | MYEOV | SPTA1 |
| TTN_0___CSMD1 | 15 | 4.0 | 0 | TTN | CSMD1 |
| MYEOV_1___TP53 | 13 | 8.0 | 1 | MYEOV | TP53 |
| TTN_0___CDC42BPA | 13 | 8.0 | 0 | TTN | CDC42BPA |
| CDC42BPA_0___OR2T2 | 12 | 10.0 | 0 | CDC42BPA | OR2T2 |
| PPFIA1_0___TTN | 12 | 10.0 | 0 | PPFIA1 | TTN |
| TTN_0___MYEOV | 12 | 10.0 | 0 | TTN | MYEOV |
| MYEOV_0___AHCTF1 | 11 | 14.0 | 0 | MYEOV | AHCTF1 |
| MYEOV_0___CDC42BPA | 11 | 14.0 | 0 | MYEOV | CDC42BPA |
| MYEOV_0___DNAH14 | 11 | 14.0 | 0 | MYEOV | DNAH14 |
| MYEOV_1___GATA3 | 11 | 14.0 | 1 | MYEOV | GATA3 |
| ZFPM2_1___ABRA | 11 | 14.0 | 1 | ZFPM2 | ABRA |
| CDC42BPA_0___MYEOV | 10 | 20.0 | 0 | CDC42BPA | MYEOV |
| CDC42BPA_0___TTN | 10 | 20.0 | 0 | CDC42BPA | TTN |
| GATA3_0___PIK3CA | 10 | 20.0 | 0 | GATA3 | PIK3CA |
| MYEOV_1___DNAH14 | 10 | 20.0 | 1 | MYEOV | DNAH14 |

## Root Node Distribution: q=1.0

| root_node | count |
| --- | --- |
| MYEOV | 34 |
| CDC42BPA | 16 |
| ZFPM2 | 16 |
| TTN | 12 |
| PLEKHA6 | 9 |
| GATA3 | 9 |
| TRIM67 | 8 |
| RIMS2 | 6 |
| PIK3CA | 5 |
| KCNU1 | 4 |
| OBSCN | 4 |
| CNBD1 | 4 |
| DNAH14 | 4 |
| INTS8 | 4 |
| PDP1 | 4 |
| RUNX1T1 | 3 |
| ERBB2 | 3 |
| SPTA1 | 3 |
| PPFIA1 | 3 |
| DPY19L4 | 3 |
| CDH1 | 3 |
| RN7SL350P | 3 |
| COL22A1 | 3 |
| MIR548A3 | 2 |
| NECAB1 | 2 |

## Root Node Distribution: q=0.75

| root_node | count |
| --- | --- |
| MYEOV | 34 |
| CDC42BPA | 16 |
| ZFPM2 | 16 |
| TTN | 12 |
| PLEKHA6 | 9 |
| GATA3 | 9 |
| TRIM67 | 8 |
| RIMS2 | 6 |
| PIK3CA | 5 |
| KCNU1 | 4 |
| OBSCN | 4 |
| CNBD1 | 4 |
| DNAH14 | 4 |
| INTS8 | 4 |
| PDP1 | 4 |
| RUNX1T1 | 3 |
| ERBB2 | 3 |
| SPTA1 | 3 |
| PPFIA1 | 3 |
| DPY19L4 | 3 |
| CDH1 | 3 |
| RN7SL350P | 3 |
| COL22A1 | 3 |
| MIR548A3 | 2 |
| NECAB1 | 2 |

## BRCA Driver Stability

| gene | q=1.0 | q=0.75 |
|---|---|---|
| PIK3CA | rank 1, alpha 196 | rank 6, alpha 59 |
| TP53 | rank 2, alpha 194 | rank 4, alpha 73 |
| GATA3 | rank 3, alpha 187 | rank 3, alpha 76 |
| CDH1 | rank 5, alpha 180 | rank 16, alpha 38 |
| ERBB2 | rank 10, alpha 151 | rank 16, alpha 38 |

## Passenger / Frequency-Sensitive Genes

| gene | q=1.0 | q=0.75 |
|---|---|---|
| TTN | rank 4, alpha 183 | rank 2, alpha 83 |
| CSMD1 | rank 8, alpha 163 | rank 5, alpha 61 |
| CSMD3 | rank 11, alpha 146 | rank 88, alpha 9 |

## Preliminary Interpretation

- q=0.75 has shallower trees than q=1.0: 6.64 vs 31.64 average depth.
- q=0.75 uses fewer unique SGAs per tree than q=1.0: 13.85 vs 63.27.
- Top-20 alpha gene overlap is 13 / 20.
- q=1.0 invariants OK: True; q=0.75 invariants OK: True.
