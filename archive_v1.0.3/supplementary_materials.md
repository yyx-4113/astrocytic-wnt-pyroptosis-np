# Supplementary Materials

**Manuscript title:** Computational re-analysis of convergent Wnt/β-catenin and pyroptosis signaling in neuropathic pain

**Corresponding author:** Yang Yong-Xin (Yongxin Yang), ORCID 0009-0004-9698-6552

**Journal of submission:** *Molecular Pain* (SAGE; Original Research Article)

> These supplementary materials accompany the honest computational re-analysis.
> They contain NO Mendelian-randomization analyses (withdrawn) and NO WGCNA
> module-preservation statistics (never computed; the earlier Z-summary = 1.0
> file was forged by self-comparison and has been deleted). All tables below are
> derived from the genuine public datasets (GSE217932, GSE134003) and are
> reproducible from the version-tagged repository.

**Supplementary table files:**
- `Supplementary_Tables_S1-S3.docx` — the three tables below as formatted Word tables (submission file).
- `Supplementary_Tables_S1-S3.xlsx` — same three tables in Excel.

> Table numbering matches the manuscript exactly: S1 = WGCNA module assignments; S2 = descriptive drug three-tier classification; S3 = CellChat Wnt–Fzd communication probability matrix (SNI).

---

## Supplementary Table S1. WGCNA Module Assignments (GSE217932, exploratory)

Signed WGCNA on the single genuine rat bulk cohort (n = 12; n = 4 per primary contrast), top 5000 variable genes, soft-threshold β = 6, targeted minimum module size = 30. Module **preservation** across an independent cohort was NOT assessed (no independent bulk NP cohort exists), so no Langfelder et al. (2011) Z-summary was computed. Assignments are descriptive and not replicated.

| Module | Size (n genes) | Wnt / pyroptosis genes assigned |
|---|---|---|
| M5 | 1395 | — |
| M7 | 1140 | Fzd10, Tcf4 |
| M3 | 905 | Dvl3, Fzd3, Gsk3b, Wnt7a |
| M8 | 886 | Ccnd1, Ctnnb1, Dkk3, Dvl1 |
| M9 | 585 | Nfkb1 |
| M6 | 39 | — |
| M2 | 26 | — |
| M1 | 22 | — |
| M4 | 2 | — |
| **Total** | **5000** | 10 Wnt + 1 pyroptosis (Nfkb1) |

Full 5000-gene assignment: `results/MC5_WGCNA_Modules.csv` in the repository.

---

## Supplementary Table S2. Descriptive Drug Three-Tier Classification

Descriptive annotation of DrugBank-listed compounds against the Wnt/β-catenin and pyroptosis nodes. Classification is **descriptive only**: no statistical testing was performed because target coverage is incomplete; tiers reflect mechanism-based relevance, not efficacy.

| Drug | Class | FDA-approved | Tier | Target(s) | n targets |
|---|---|---|---|---|---|
| MCC950 | NLRP3i | No | 1 | NLRP3 | 1 |
| CY-09 | NLRP3i | No | 1 | NLRP3 | 1 |
| VX-765 | CASP1i | No | 1 | CASP1 | 1 |
| Disulfiram | GSDMDi | Yes | 1 | GSDMD,ALDH2 | 2 |
| XAV939 | Wnti | No | 1 | TNKS1,TNKS2 | 2 |
| IWR-1 | Wnti | No | 1 | AXIN2 | 1 |
| LGK-974 | PORCNi | No | 1 | PORCN | 1 |
| ICG-001 | CBPi | No | 1 | CTNNB1,CREBBP | 2 |
| PRI-724 | CBPi | No | 1 | CTNNB1,CREBBP | 2 |
| KYA1797K | WntD | No | 2 | CTNNB1 | 1 |
| LY2090314 | GSK3i | No | 2 | GSK3B | 1 |
| Ruxolitinib | JAKi | Yes | 2 | JAK1,JAK2 | 2 |
| Tofacitinib | JAKi | Yes | 2 | JAK1,JAK3 | 2 |
| Gabapentin | a2d | Yes | 3 | CACNA2D1,CACNA2D2 | 2 |
| Pregabalin | a2d | Yes | 3 | CACNA2D1,CACNA2D2 | 2 |
| Duloxetine | SNRI | Yes | 3 | SLC6A4,SLC6A2 | 2 |
| Venlafaxine | SNRI | Yes | 3 | SLC6A4,SLC6A2 | 2 |
| Amitriptyline | TCA | Yes | 3 | SLC6A4,SLC6A2,HRH1 | 3 |
| Tramadol | Opioid | Yes | 3 | OPRM1,SLC6A4 | 2 |
| Capsaicin | TRPV1 | Yes | 3 | TRPV1 | 1 |

---

## Supplementary Table S3. CellChat Wnt–Fzd Communication Probability Matrix (SNI condition, GSE134003)

CellChat aggregated communication probabilities (sender × receiver) across 14 annotated spinal-cord cell types + "Other", computed on the SNI (spared-nerve-injury) single-cell cohort. Wnt-family ligand → Frizzled (Fzd) receptor pairs (e.g., Wnt7a–Fzd3/Fzd10) are the subset driving Figure 5. Values rounded to 3 decimals.

| Sender \\ Receiver | Astrocytes | Microglia | Oligodendrocytes | Myelinating_Oligo | OPCs | Excitatory_Neurons | Inhibitory_Neurons | Endothelial | Pericytes | Vascular_SMC | Myelinating_Schwann | NM_Schwann | Macrophages | Leptomeningeal | Other |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Astrocytes | 0.785 | 0.054 | 0.002 | 0.005 | 0.028 | 0.015 | 0.062 | 0.142 | 0.192 | 0.167 | 0.000 | 0.161 | 0.000 | 0.661 | 0.000 |
| Microglia | 0.012 | 0.005 | 0.009 | 0.013 | 0.026 | 0.009 | 0.013 | 0.014 | 0.019 | 0.032 | 0.023 | 0.014 | 0.008 | 0.047 | 0.019 |
| Oligodendrocytes | 0.263 | 0.037 | 0.029 | 0.024 | 0.026 | 0.016 | 0.037 | 0.037 | 0.061 | 0.130 | 0.011 | 0.122 | 0.004 | 0.432 | 0.010 |
| Myelinating_Oligo | 0.368 | 0.051 | 0.039 | 0.031 | 0.032 | 0.019 | 0.048 | 0.049 | 0.080 | 0.175 | 0.006 | 0.170 | 0.003 | 0.602 | 0.007 |
| OPCs | 0.475 | 0.082 | 0.048 | 0.089 | 0.269 | 0.066 | 0.104 | 0.298 | 0.165 | 0.242 | 0.061 | 0.141 | 0.126 | 0.589 | 0.119 |
| Excitatory_Neurons | 0.199 | 0.042 | 0.044 | 0.055 | 0.101 | 0.044 | 0.069 | 0.194 | 0.128 | 0.208 | 0.084 | 0.115 | 0.041 | 0.396 | 0.083 |
| Inhibitory_Neurons | 0.351 | 0.110 | 0.131 | 0.174 | 0.348 | 0.131 | 0.189 | 0.319 | 0.263 | 0.489 | 0.264 | 0.255 | 0.174 | 0.895 | 0.268 |
| Endothelial | 0.109 | 0.043 | 0.068 | 0.077 | 0.107 | 0.065 | 0.092 | 0.117 | 0.140 | 0.250 | 0.163 | 0.130 | 0.061 | 0.426 | 0.141 |
| Pericytes | 0.030 | 0.015 | 0.026 | 0.030 | 0.045 | 0.027 | 0.036 | 0.038 | 0.047 | 0.102 | 0.065 | 0.058 | 0.019 | 0.179 | 0.061 |
| Vascular_SMC | 0.058 | 0.026 | 0.035 | 0.053 | 0.140 | 0.046 | 0.054 | 0.024 | 0.011 | 0.188 | 0.047 | 0.178 | 0.000 | 0.502 | 0.100 |
| Myelinating_Schwann | 0.015 | 0.011 | 0.015 | 0.016 | 0.020 | 0.014 | 0.019 | 0.028 | 0.030 | 0.053 | 0.031 | 0.035 | 0.010 | 0.105 | 0.027 |
| NM_Schwann | 0.036 | 0.124 | 0.094 | 0.073 | 0.059 | 0.060 | 0.081 | 0.402 | 0.204 | 0.253 | 0.000 | 0.314 | 0.000 | 0.866 | 0.007 |
| Macrophages | 0.076 | 0.005 | 0.000 | 0.000 | 0.003 | 0.001 | 0.006 | 0.014 | 0.019 | 0.016 | 0.000 | 0.016 | 0.000 | 0.064 | 0.000 |
| Leptomeningeal | 0.610 | 0.347 | 0.420 | 0.524 | 0.942 | 0.424 | 0.601 | 2.981 | 1.244 | 1.944 | 0.743 | 0.971 | 0.247 | 2.996 | 0.736 |
| Other | 0.122 | 0.031 | 0.030 | 0.029 | 0.038 | 0.029 | 0.045 | 0.412 | 0.122 | 0.203 | 0.027 | 0.093 | 0.009 | 0.314 | 0.042 |

---

**Data availability.** All processed data and analysis code are in the public, version-tagged repository (https://github.com/yyx-4113/astrocytic-wnt-pyroptosis-np). No "available on request" clause is used.
