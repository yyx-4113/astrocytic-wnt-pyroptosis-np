# MANIFEST — Reproducibility Package (astrocytic-wnt-pyroptosis-np)

Generated 2026-10-28. All checksums via md5 of the packaged files.
Repository: https://github.com/yyx-4113/astrocytic-wnt-pyroptosis-np (public, tag v1.0.0).

## data/ (genuine processed outputs; no MR / cross-dataset / GSE186237)
| File | MD5 |
|---|---|
| Figure10_summary.json | e306bbabf65b2f7d56aa306e546bad72 |
| Figure5_communication_SNI.csv | 65cb8b596d466b209a6b2dc875679041 |
| Figure6_deconv_GSE217932.csv | 70195d28d46e3339991c871f28c822f2 |
| Figure7_HubGenes.csv | af1c3d5427db5f1f26fea76f5d5acb62 |
| Figure7_summary.json | 3eb977ea0f1cca3bf523a3431970da6a |
| Figure8_DrugConnectivity.csv | 2b0ff8a6187bbfef02b4307bb94ca449 |
| Figure8_TargetEnrichment.csv | d5dd3a5987ec4fe9b4c0b6a59a8a5be7 |
| Figure8_summary.json | 573bd0ca96bb60d636d24b328392f4c5 |
| Figure9_Perturbation_OE_KD.csv | f2425a12e7daa078063386bd88faede7 |
| Figure9_RWR_Scores.csv | 097f9c52e9f7b6230261ff90f811a2a8 |
| Figure9_VirtualKO.csv | 99d7e1e9a5b5950b47f96578d582f9f8 |
| Figure9_summary.json | e416f6775e76a14c04fd05cab6fc2110 |
| GSE217932_DEGs.csv | a7b4df4f57a2145df1b503cd2dfb9e8c |
| MC2_Drug_ThreeTier.csv | 818e2ecf78d5ad01bb0b612b37f8bd94 |
| MC4_Pathway_Comparison.csv | 650767009dc2df676b1fb0975e1af394 |
| MC5_GeneModuleMembership.csv | fb4fca1595975a1a75601e08dc8da1f8 |
| MC5_WGCNA_Modules.csv | a08d0dff2beddc74f6669d2f88d55058 |
| MC5_WGCNA_Summary.json | 886f40626af066a074c6f354a7b2f278 |

## figures/ (12 genuine figures, Figure0–Figure11)
| File | MD5 |
|---|---|
| Figure0_Graphical_Abstract.png | d531191ca0e52f2d0a54fa0302057a2a |
| Figure10_Structural_Analysis.png | acc477b3ed57f754f9e20468a7f0c681 |
| Figure11_Reproducibility_Audit.png | 2a3f9ec37777c05d7233ced40907671e |
| Figure1_DEG_overview.png | 597420ee0b40a2e07acd2cf848a42cd0 |
| Figure2_Wnt_pathway.png | 34b266d213c25907140d1bcb6355c8bd |
| Figure3_Pyroptosis_Wnt_crosstalk.png | 08da2cc8336490ce22d6b40949155c7b |
| Figure4_scRNAseq_cell_atlas.png | cd0797d6fce006a8607669ab076c9b56 |
| Figure5_CellChat_Pseudotime.png | de21783dc96111cd0d4367e85e558507 |
| Figure6_Deconvolution.png | 3eb0bea115216f44037e2b42f97a89fa |
| Figure7_PPI_Hub_TCF4.png | 1a238b6daacf1e1a85d886ddd326b1ce |
| Figure8_Drug_ThreeTier_Classification.png | c1f08a92d8d420b7575e7c4700a8a261 |
| Figure9_GRN_VirtualKO.png | 007d4c57cb909f9cc489ddf1a1e6a9c3 |

## scripts/ (honest regeneration pipeline: rebuild_*.py)
| Script | Purpose |
|---|---|
| rebuild_01_gse186237_deg.py | Method check: reproduces GSE217932 DEG (proves r=1 identity) |
| rebuild_02_raw_compare.py | First raw-count identity check |
| rebuild_03_definitive_raw_compare.py | Definitive 328,624-cell identity proof (GSE186237=GSE217932) |
| rebuild_04_ref_audit.py | Crossref reference-DOI audit |
| rebuild_08_genuine_single_dataset.py | Regenerates genuine single-cohort DEG |
| rebuild_09_scRNA_astrocyte_pyroptosis.py | GSE134003 single-cell re-analysis |
| rebuild_10_fig4.py | Fig 4 single-cell atlas |
| rebuild_11_fig1.py | Fig 1 single-dataset DEG overview |
| rebuild_12_fig6.py | Fig 6 CIBERSORTx deconvolution (9.6%->9.9%) |
| rebuild_13_fig0.py | Fig 0 graphical abstract |
| rebuild_14_wgcna_modules.py | Honest WGCNA module assignment (no preservation) |
| rebuild_16_fig2_gsea_plot.py | Fig 2 genuine GSEA (NES ~= -1.6, negative) |
| rebuild_17_fig10_structural.py | Fig 10 structural analysis (literature Kd; 0/4 MD) |
| rebuild_18_fig11_audit.py | Fig 11 reproducibility self-audit matrix |
| rebuild_19_build_docx.py | Rebuilds manuscript_final_with_figures.docx |
| rebuild_20_assemble_repo.py | Assembles this repository |

## Integrity flags (resolved)
1. Single bulk cohort only (GSE217932, n=4/condition); no cross-dataset replication claimed.
2. GSE186237 is a byte-identical duplicate of GSE217932 (MD5 a7b4df4f57a2145df1b503cd2dfb9e8c) and is excluded.
3. Mendelian-randomization section withdrawn; no MR/IV results reported.
4. WGCNA module preservation (forged Z_summary=1.0) deleted; only descriptive module assignment retained.
5. Deconvolution astrocyte fraction 9.6%->9.9% (p=0.343, not significant); not cross-dataset.
