# Computational re-analysis of convergent Wnt/β-catenin and pyroptosis signaling in neuropathic pain, with single-cell resolution of cell-type specificity: a hypothesis-generating study

**Authors:** Yang Yong-Xin
**Affiliations:** Department of Anesthesiology, The Second Affiliated Hospital of Fujian University of Traditional Chinese Medicine, Fuzhou, Fujian, China
**Corresponding author:** Yang Yong-Xin, 960856791@qq.com
**Running title:** Wnt–Pyroptosis Signaling in Neuropathic Pain: A Computational Re-analysis
**Keywords:** Neuropathic pain, Wnt/β-catenin, Pyroptosis, Astrocyte, Drug repurposing, Single-cell RNA-seq, Computational biology, Hypothesis generation
**Figures:** 12 (Figures 0–11)  |  **Supplementary Tables:** 3  |  **References:** 26

---

## ABSTRACT

Neuropathic pain involves neuroinflammatory signalling, but whether Wnt/β-catenin and NLRP3-inflammasome-mediated pyroptosis converge in specific spinal-cord cell types is unresolved. We performed a computational, hypothesis-generating re-analysis of publicly available data: one rat spinal-cord bulk RNA-seq cohort (GSE217932; n=4 per condition, no independent replication) and one mouse spared-nerve-injury single-cell RNA-seq cohort (GSE134003), together with STRING protein–protein interactions, DrugBank annotations, and four Protein Data Bank structures. Bulk analysis showed no significant Wnt/β-catenin enrichment after correction (NES ≈ −1.6, nominal p ≈ 0.03). Single-cell analysis localized pyroptosis gene expression (Nlrp3, Gsdmd, Casp1) predominantly to microglia and macrophages, whereas astrocytes showed low pyroptosis expression without injury-associated change; astrocyte-specific Wnt–pyroptosis co-activation was not supported. Protein–protein interaction analysis (216 nodes, 189 edges) identified Il6, Ccl2, Wnt2, Nfkb1 and Rela as top hubs. A descriptive three-tier classification nominated nine Tier-1 drug candidates, including MCC950 (NLRP3, Kd=5 nM), disulfiram (GSDMD) and XAV939 (TNKS2, Kd=46 nM). Reported limitations include a single bulk cohort without replication, small samples, cross-species data and absent wet-lab validation; all findings are computational predictions requiring experimental confirmation, and no simulated data were used.

---

## 1. INTRODUCTION

Neuropathic pain (NP) affects 7–10% of the general population and represents a major unmet clinical need, with fewer than 50% of patients achieving adequate pain relief from existing pharmacotherapy (van Hecke et al., 2014; Finnerup et al., 2015). The molecular pathogenesis of NP involves complex neuroimmune interactions within the spinal dorsal horn, in which glial activation plays a central role in pain initiation and chronification.

Wnt/β-catenin signaling has been implicated in the pathogenesis of NP (Zhang et al., 2013; Itokazu et al., 2014), and NLRP3 inflammasome-mediated pyroptosis is a recognized neuroinflammatory mechanism (Bortolotti et al., 2018). Astrocytes, the most abundant glial cell type in the spinal cord, are well-established contributors to NP maintenance (Ji et al., 2019). Two previously independent lines of investigation, Wnt signaling in synaptic plasticity and pain (Zhang et al., 2013; Itokazu et al., 2014) and NLRP3 inflammasome-mediated pyroptosis in neuroinflammation (Bortolotti et al., 2018), have not been jointly examined in astrocytes. This leaves an open question: is Wnt–pyroptosis co-activation restricted to microglia, or does it also occur in astrocytes?

The present study addresses this question through parallel computational analyses integrating rodent transcriptomics (bulk and single-cell) with drug-target pharmacology and structural biology. A graphical overview of the study design is provided in Figure 0. This work is positioned as a hypothesis-generating computational exploration rather than a de novo discovery or a definitive mechanistic demonstration. We emphasize that all findings are computational predictions requiring independent experimental validation, and we explicitly report negative and null results alongside positive ones.

A distinctive feature of this study is the systematic reporting of negative and null results. Specifically, we report: Wnt and pyroptosis genes that were not differentially expressed; the limited statistical power of bulk differential expression (only 2 genes significant after correction at n=4 per group); molecular dynamics simulations that could not be completed (0/4); the absence of an independent replication cohort for the bulk analysis; and statistical limitations that preclude certain claims (no Fisher exact test for drug-target enrichment). We contend that such complete reporting strengthens scientific inference and provides a more accurate representation of the evidence than selective presentation of supportive results alone.

## 2. MATERIALS AND METHODS

### 2.1 Data Acquisition

**Bulk RNA-seq.** One rat neuropathic-pain transcriptomic dataset was analyzed. GSE217932 comprises spinal dorsal horn tissue from 12 male rats (4 sham, 4 chronic constriction injury [CCI], and 4 CCI plus the experimental compound THP) profiled by bulk RNA-seq. The primary differential-expression contrast uses the sham (n=4) versus CCI (n=4) groups; the CCI+THP arm was reserved for a separate intervention contrast and is not used for the primary sham-versus-CCI comparison. Expression values are log2(normalized count + 1). CCI is a well-established rodent model of neuropathic pain (Decosterd & Woolf, 2000). The n=4 per condition limitation is acknowledged; all bulk-derived results are exploratory. No second, independent bulk dataset is available in this study; claims of cross-dataset replication are therefore not made.

**Single-cell RNA-seq.** GSE134003 (mouse spinal cord after SNI, 10x Genomics platform; 18,855 cells across 14 annotated cell types, 19 animals) was analyzed at full-transcriptome resolution (14,700 genes). Because pre-filtering to a 140-gene Wnt/pyroptosis set before annotation can artifactually define the result, pyroptosis and Wnt module scores were computed across all genes and compared across annotated cell types; cell-type specificity was assessed at the animal level (SI: n=10 animals; SNI: n=9 animals), not by treating individual cells as independent replicates (Squair et al., 2021).

**Protein structures.** Four PDB structures were retrieved: 7PZC (NLRP3 NACHT domain), 6N9O (CTNNB1/β-catenin), 1JDH (TNKS2 tankyrase-2), and 3UA9 (GSDMD pore). All were processed with PDBFixer (Eastman et al., 2017) for missing-atom remediation. Molecular dynamics simulations were attempted using OpenMM (Eastman et al., 2017) but 0/4 structures completed due to computational resource limitations; binding energies reported herein are from published experimental literature.

**Gene annotation.** Rat transcript identifiers (Ensembl release 105) were mapped to gene symbols using the rat_gene_map.pkl dictionary (49,331 entries). Transcript-level expression values were aggregated to gene level by mean.

**Data reuse disclosure.** GSE217932 was used in three complementary analyses (differential expression, WGCNA, GRN). These are not independent replications; results from these analyses reflect different analytical perspectives on the same underlying data. The single-cell (GSE134003) and structural (PDB) datasets are independent of the bulk cohort but address different biological questions.

**Gel and blot images.** This study is entirely computational: no new experimental data were generated, and no gel electrophoresis or Western blot experiments were performed. Consequently, no uncropped gel or blot images exist for this manuscript. All figures present quantitative analyses of publicly available transcriptomic, structural, and pharmacological data. Any experimental validation proposed in the Discussion (e.g., ChIP-qPCR, Western blot of NLRP3/GSDMD, behavioral testing) is a future experimental direction and has not been conducted in the present work.

### 2.2 Differential Expression

Differentially expressed genes (DEGs) were identified using a threshold of |log2 fold change| > 0.5 and nominal p-value < 0.05, with Benjamini–Hochberg (BH) adjusted p-value reported. Because only a single bulk cohort is available (n=4 per condition), no cross-dataset aggregation (e.g., Robust Rank Aggregation) was performed; all bulk DEG results are exploratory. Wnt pathway genes were defined by the KEGG Wnt signaling pathway (rno04310, 167 genes), of which 32 were selected for focused analysis based on core signaling components (ligands, receptors, transducers, and transcriptional effectors); pyroptosis genes were curated from literature (16 genes: Nlrp3, Pycard, Casp1, Gsdmd, Il1b, Il18, Nlrp1, Nlrc4, Aim2, Casp4, Casp5, Gsdme, Il1a, Hmgb1, Tnf, Nfkb1). We explicitly report all 32 Wnt and 16 pyroptosis genes detected in expression data, including those not meeting the DEG threshold.

### 2.3 Functional Enrichment and Co-Expression Network Analysis

Gene Set Enrichment Analysis (GSEA) was performed against MSigDB Hallmark gene sets (v2023.2) using gene rankings by log2 fold change. Weighted Gene Co-expression Network Analysis (WGCNA; Langfelder & Horvath, 2008) was conducted with a signed network (Spearman correlation, soft-threshold power beta=6, targeted minimum module size of 30) on the top 5,000 most variable genes in GSE217932. Due to the limited sample size (n=12 total, n=4 per primary contrast), permutation-based p-values are not reported and WGCNA results are interpreted qualitatively and treated as exploratory. Module preservation across an independent dataset was not assessed because no independent bulk cohort is available in this study; the Langfelder et al. (2011) Z-summary statistic therefore could not be applied to a genuine test–reference split.

### 2.4 Single-Cell Analysis

Cell types were annotated using canonical markers: Astrocytes (Gfap, Aldh1l1), Microglia (Aif1, Cx3cr1), Oligodendrocytes (Mbp, Mog), Neurons (Rbfox3, Syt1), Endothelial cells (Cldn5), Pericytes (Pdgfrb), and others. CellChat v1.6 (Jin et al., 2021) was used for ligand–receptor interaction detection with the triMean method and 1,000 permutations. Given the pre-filtered gene set (140 Wnt/pyroptosis-related genes), CellChat analysis was restricted to Wnt-Fzd receptor–ligand pairs; cross-pathway quantitative comparisons are therefore not feasible. Pseudotime trajectory inference was performed with Monocle3 (Cao et al., 2019), using DDRTree dimensionality reduction with root selection at minimal pseudotime in homeostatic astrocytes.

### 2.5 Cell Type Deconvolution

CIBERSORTx (Newman et al., 2019) deconvolution was performed with a custom signature matrix derived from the annotated scRNA-seq data. Astrocyte proportions were compared between sham and CCI conditions within the single GSE217932 bulk cohort by paired analysis. No cross-dataset validation was performed because only one bulk cohort is available.

### 2.6 Protein–Protein Interaction Network

High-confidence protein–protein interactions (combined score >0.7) were retrieved from STRING v12.0 (Szklarczyk et al., 2023) for the differentially expressed genes. Hub genes were identified by multi-centrality scoring: degree, betweenness, closeness, and eigenvector centrality were computed using NetworkX, and genes ranking in the top 20th percentile for at least three of four metrics were classified as hubs. TCF4 transcription factor binding-site prediction used JASPAR 2022 (Castro-Mondragon et al., 2022), matrix MA0830.1, with a relative profile score threshold ≥85%. Putative TCF4 binding sites in Nlrp3 and Gsdmd promoter regions (−2,000 to +200 bp relative to transcription start site) were reported as hypothetical.

### 2.7 Drug–Target Mapping: Descriptive Three-Tier Classification

Drug–target annotations were obtained from DrugBank v5.1 (Wishart et al., 2018) and supplemented by manual curation of published literature. Drug identities were verified by PubChem CID lookup (20/20 confirmed). We identified 21 unique molecular targets across 20 drugs spanning three pharmacological categories. Fisher exact test and other statistical enrichment analyses for drug–target pathway associations were not performed. The 21 targets represent an estimated 10–20% of true DrugBank interactions for these 20 drugs, providing incomplete coverage. Testing pre-selected, pathway-biased targets against the same pathway gene sets would also constitute circular reasoning. Instead, we report a descriptive three-tier expert annotation intended to generate hypotheses for preclinical testing, not an evidence-based ranking.

- **Tier 1 (Direct pathway targeting, n=9):** Drugs whose documented molecular targets include core Wnt pathway proteins (CTNNB1, AXIN2, GSK3B, TNKS1/2, PORCN) or pyroptosis pathway proteins (NLRP3, GSDMD, CASP1). Examples: MCC950 (NLRP3 inhibitor, Kd=5 nM; Coll et al., 2019), disulfiram (GSDMD pore blocker, FDA-approved; Hu et al., 2020), XAV939 (TNKS2 inhibitor, Kd=46 nM; Huang et al., 2009), and six others.
- **Tier 2 (Indirect/modulatory targeting, n=4):** Drugs modulating upstream regulators or paralogous nodes without direct documented binding to the core Wnt/pyroptosis proteins.
- **Tier 3 (Exploratory/repositioning candidates, n=7):** Compounds with literature-reported anti-neuroinflammatory or anti-NP activity but without a clearly mapped Wnt/pyroptosis target.

### 2.8 Gene Co-Expression Network with Virtual Perturbation

An ARACNe-style gene regulatory network (Margolin et al., 2006) was inferred from GSE217932 using Spearman |rho| > 0.75 (n=12 samples; no permutation-based false-discovery-rate correction was applied because 31,125 pairwise tests at n=12 lack statistical power). Edges were treated as undirected, unsigned correlations; no causal direction was inferred. Virtual knockout/overexpression/knockdown (KO/OE/KD) was implemented as random-walk-with-restart (RWR) seed-weight perturbation (±2.0 SD) on this correlation network. The resulting scores are network-propagation outputs and do not represent measured gene-expression changes. All virtual-perturbation results are hypothesis-generating only and must not be read as biological conclusions.

### 2.9 Structural Analysis

Four PDB structures (7PZC, 6N9O, 1JDH, 3UA9) were prepared with PDBFixer (Eastman et al., 2017) and inspected. Molecular dynamics simulations were attempted with OpenMM (Eastman et al., 2017) but 0/4 completed; binding energies reported are from published experimental literature (Coll et al., 2019; Huang et al., 2009). No new MD trajectories or de novo binding-energy calculations were performed; the structural section is therefore descriptive and does not constitute structural validation of the drug priorities.

## 3. RESULTS

### 3.1 Single-Dataset Differential Expression of GSE217932

Differential expression between CCI and sham (n=4 each) in GSE217932 identified 1,698 genes at the exploratory nominal threshold (|log2FC|>0.5, nominal p<0.05), of which only 2 genes (1 up, 1 down) survived BH correction (padj<0.05). This weak post-correction signal reflects limited statistical power at n=4 per group and is reported transparently. Of 32 Wnt pathway genes tested, 8 were differentially expressed and 20 were not at the exploratory nominal threshold; for pyroptosis, 15 of 16 pathway genes were detected and 6 were differentially expressed. These counts are suggestive rather than definitive given the limited power and are reported alongside the null result that most pathway genes were not significant after correction.

### 3.2 GSEA Confirms Wnt/β-catenin Signaling Enrichment

GSEA against MSigDB Hallmark gene sets (v2023.2), using log2FC-ranked genes from GSE217932, did not show enrichment of the Wnt/β-catenin signaling hallmark in the CCI (injury) direction. A recomputation on the genuine ranking found the hallmark set modestly depleted (NES ≈ −1.6, nominal one-sided p ≈ 0.03). This difference does not survive correction for the number of Hallmark sets tested and is consistent with the weak bulk differential-expression signal (only 2 genes significant after BH correction). Core Wnt genes are shown as a heatmap for completeness. The absence of bulk Wnt/β-catenin enrichment indicates that any Wnt–pyroptosis signal is cell-type-restricted (see single-cell results) rather than a uniform bulk transcriptome effect. Because only a single bulk cohort is available, no cross-dataset GSEA comparison was performed.

### 3.3 WGCNA Reveals Co-Expression Modules (Exploratory, No Independent Replication)

WGCNA across the 5,000 most variable genes in GSE217932 identified co-expression modules; Wnt/pyroptosis genes distributed across several modules. Module assignments are reported descriptively. Module preservation across an independent dataset was not assessed (no independent bulk cohort available; the Langfelder et al., 2011 Z-summary requires a genuine test–reference split and could not be applied). WGCNA results are therefore exploratory and should not be interpreted as evidence of reproducible module structure.

### 3.4 Single-Cell Analysis Reveals Cell-Type-Specific Wnt and Pyroptosis Expression

Full-transcriptome re-analysis of GSE134003 (14,700 genes; 18,855 cells; 19 animals) showed that pyroptosis genes (Nlrp3, Pycard, Casp1, Gsdmd, Il1b, Tnf) are expressed predominantly in microglia and macrophages, with mean module scores of +0.70 and +0.52 respectively, whereas astrocytes showed the lowest pyroptosis expression among major cell types (module score −0.07, 25% detection rate). At the animal level (SI: n=10; SNI: n=9), astrocyte pyroptosis scores did not differ between SNI and SI (delta = −0.03, not significant), whereas several pyroptosis genes changed significantly within microglia. Wnt signaling genes were more broadly distributed across cell types. In short, the single-cell data support microglial/macrophage pyroptosis in this model but do not support astrocyte-specific pyroptosis co-activation; the originally hypothesized astrocyte mechanism remains untested. These displays were computed across all genes rather than a pre-filtered 140-gene set, addressing the artifact risk of pre-annotation gene restriction (see Limitations).

### 3.5 CellChat Detects Wnt-Fzd Signaling; Cross-Pathway Comparison Limited by Gene Set

CellChat detected Wnt-Fzd ligand–receptor pairs between astrocytes and neurons/microglia in the SNI condition. Cross-pathway dominance comparisons remain limited because the Wnt/pyroptosis-focused hypothesis constrains interpretation, and the relative contribution of Wnt versus other injury signals in spinal cord cannot be fully resolved from this dataset.

### 3.6 Deconvolution Within GSE217932

CIBERSORTx deconvolution of GSE217932 (sham vs CCI) estimated astrocyte proportions. Within this single cohort, the estimated astrocyte fraction changed from 9.6% (sham) to 9.9% (CCI) with nominal p=0.343; the change was not statistically significant. No cross-dataset correlation was computed because only one bulk cohort is available. An independent orthogonal measurement of astrocyte abundance was not performed; the deconvolution result is therefore a computational estimate requiring experimental confirmation.

### 3.7 PPI Network Identifies Inflammatory Hubs and TCF4 as a Transcriptional Link

PPI analysis (216 nodes, 189 edges) identified Il6, Ccl2, Wnt2, Nox4, Rela, and Nfkb1 as top hubs (multi-centrality ranking). JASPAR predicted TCF4 (MA0830.1) binding sites in Nlrp3 and Gsdmd promoter regions. These motif predictions are hypothetical regulatory links and do not establish promoter occupancy or regulation; they are reported as testable hypotheses.

### 3.8 Descriptive Drug Three-Tier Classification Identifies Candidates for Preclinical Evaluation

Descriptive three-tier drug classification identified nine Tier 1 candidates directly targeting Wnt/pyroptosis nodes: MCC950 (Kd=5 nM, NLRP3), disulfiram (FDA-approved GSDMD inhibitor), XAV939 (Kd=46 nM, TNKS2), and six others. This classification is a descriptive catalog, not a statistically validated ranking (see Methods 2.7 and Limitations).

### 3.9 Gene Regulatory Network Analysis Predicts Wnt-to-Pyroptosis Directionality (Hypothesis-Generating)

An ARACNe-style GRN from GSE217932 (Spearman |rho|>0.75, no permutation FDR) was used for RWR-based virtual perturbation. Tcf4 and Ctnnb1 were identified as high-impact hubs; RWR seed-weight perturbation (+2.0 SD) of Wnt/β-catenin nodes propagated most strongly to pyroptosis genes, with combined RWR scores of Il1b = 0.011, Gsdmd = 0.009, Nlrp3 = 0.007, Casp1 = 0.009, and Pycard = 0.011 (network-propagation outputs, not measured expression changes). These are in silico predictions, not causal inferences, and are reported as hypothesis-generating only (see Methods 2.8 and Limitations).

### 3.10 Structural Analysis: Literature Binding Energies, MD Not Completed

Four PDB structures were prepared and inspected: 7PZC (NLRP3 NACHT domain), 6N9O (CTNNB1/β-catenin), 1JDH (TNKS2 tankyrase-2), and 3UA9 (GSDMD pore). Published experimental binding data are reported: MCC950–NLRP3 Kd=5 nM (Coll et al., 2019); XAV939–TNKS2 Kd=46 nM (Huang et al., 2009). MD simulations were not completed for any of the four structures (0/4); no new binding-energy calculations were performed (structural details and literature binding affinities are provided in Supplementary Figure S1).

## 4. DISCUSSION

### 4.1 Cell-Type Specificity of Wnt–Pyroptosis Signaling: A Hypothesis-Generating Reassessment

Wnt/β-catenin signaling (Zhang et al., 2013; Itokazu et al., 2014) and NLRP3 inflammasome-mediated pyroptosis (Bortolotti et al., 2018) are each independently linked to NP, and both have been studied primarily in the microglial context. The present single-cell re-analysis found that pyroptosis signaling localizes to microglia and macrophages rather than astrocytes, and therefore does not support astrocyte-specific Wnt–pyroptosis co-activation in GSE134003. This is a hypothesis-generating reassessment: the originally proposed astrocyte mechanism remains experimentally untested and could be addressed by future astrocyte-restricted perturbation with a pyroptosis readout and pain endpoint. We distinguish reactive transcriptional states, pathway activation, and pyroptotic cell death, and we note that non-pyroptotic astrocyte effects on nociceptive transmission may also be relevant.

### 4.2 Drug Prioritization as a Translational Bridge

The descriptive three-tier classification identifies candidates (MCC950, disulfiram, XAV939, and others) for preclinical evaluation. We emphasize that this is a candidate catalog, not an evidence-based ranking: no statistical testing was performed, target coverage is incomplete (10–20%), and pharmacological reversal requires concordant molecular and functional outcomes. The classification should guide future experimental testing rather than inform clinical use.

### 4.3 Negative and Null Findings

A distinctive feature of this study is the systematic reporting of findings that did not support, or that qualified, our primary hypothesis. We consider this transparency essential for computational research, where selective reporting of only positive results can create a misleadingly coherent narrative that overstates the strength of evidence.

First, 20 of 28 detected Wnt pathway genes (71.4%) were not differentially expressed at the exploratory threshold, and only 2 genes (1 up, 1 down) survived multiple-testing correction in the bulk cohort, indicating that bulk Wnt/pyroptosis signal at n=4 per group is weak. Second, WGCNA module preservation could not be assessed because no independent bulk cohort was available; WGCNA results are exploratory. Third, molecular dynamics simulations were not completed for any of the four PDB structures (0/4), limiting structural characterization to published experimental binding energies. Fourth, the scRNA-seq analysis was restricted to a pre-filtered set of 140 Wnt/pyroptosis genes, precluding quantitative cross-pathway CellChat comparisons. Fifth, virtual-perturbation scores are network-propagation outputs, not measured expression changes, and are reported as hypothesis-generating only.

We contend that transparent reporting of these negative and null findings strengthens, rather than weakens, the credibility of this computational study. By providing readers with a complete picture of both the signals detected and the limitations encountered, we aim to facilitate realistic interpretation and guide future experimental validation efforts.

### 4.4 Limitations

Beyond the negative findings detailed in Section 4.3, we acknowledge the following additional limitations:

1. **Single bulk cohort, no independent replication.** This study analyzes one bulk RNA-seq cohort (GSE217932, n=4 per condition). No independent bulk dataset is available; claims of cross-dataset replication are not made. All bulk-derived results are exploratory.
2. **Small sample sizes.** The bulk cohort has n=4 per primary contrast. All bulk-derived results, including WGCNA modules (no permutation-based FDR) and GRN edges (no permutation-based FDR), should be interpreted as exploratory.
3. **Incomplete drug-target annotation.** Only 21 targets were annotated for 20 drugs (estimated 10–20% of true DrugBank interactions). The three-tier classification is a descriptive expert annotation; complete DrugBank API queries would enable unbiased statistical enrichment analyses.
4. **Data reuse.** GSE217932 was used in three complementary analyses (differential expression, WGCNA, GRN). These are not independent replications but different analytical perspectives on the same underlying data.
5. **Cross-species integration.** Rat bulk transcriptomic data, mouse scRNA-seq data, and human structural/pharmacology data were integrated. Cross-species extrapolation requires caution, particularly for drug target validation where species-specific pharmacokinetics may differ.
6. **scRNA-seq gene-set restriction (corrected by full-transcriptome re-analysis).** An earlier version of this analysis pre-filtered the single-cell input to 140 Wnt/pyroptosis-related genes before annotation, which can artifactually define cell-type specificity. The single-cell conclusions reported here were therefore recomputed across the full 14,700-gene transcriptome and assessed at the animal level (n=19 animals), which showed pyroptosis enrichment in microglia/macrophages rather than astrocytes. Cell-type identity and cross-pathway dominance should still be confirmed in an independent cohort.
7. **No simulated data were used.** All analyses use real, publicly available datasets. No synthetic or simulated data were generated. The scope, limitations, and reproducibility posture of each analytical module are summarized in a self-audit matrix (Supplementary Figure S2).

### 4.5 Future Directions

**Experimental validation priorities.** (1) Astrocyte-specific Ctnnb1 conditional knockout (GFAP-CreERT2; Ctnnb1 flox/flox) in the SNI model, with behavioral pain endpoints (von Frey, Hargreaves) and NLRP3/GSDMD expression measurement by qPCR and Western blot. (2) TCF4 ChIP-qPCR at Nlrp3 and Gsdmd promoters in primary spinal cord astrocytes treated with Wnt3a, to test the JASPAR-predicted binding sites identified in this study. (3) Preclinical evaluation of MCC950 (10 mg/kg, i.p.) and disulfiram (50 mg/kg, i.p.) in a rodent SNI model, prioritized from Tier 1 of our drug classification. (4) Orthogonal measurement (e.g., immunohistochemistry or flow cytometry) of astrocyte abundance to confirm the deconvolution estimate.

**Computational refinement priorities.** (1) Analysis of a second, genuinely independent bulk NP transcriptomic cohort to test reproducibility of the Wnt/pyroptosis signal. (2) WGCNA module preservation analysis in larger independent NP cohorts (n≥20 per condition) to overcome the instability observed at the current sample size. (3) Complete DrugBank REST API query for unbiased drug-target enrichment with full target annotation. (4) Re-analysis of the scRNA-seq data without the pre-filtering restriction, with sample-level QC and marker validation, to establish astrocyte identity and cross-pathway context. (5) Completion of molecular dynamics simulations for the four drug-target complexes to provide computational binding-energy validation, or replacement with molecular docking scores.

## ACKNOWLEDGMENTS

The author gratefully acknowledges the Gene Expression Omnibus (GEO), the STRING database, DrugBank, and the Protein Data Bank (PDB) for making their resources publicly available. This study would not have been possible without open data sharing.

The author used an AI-assisted language tool (a large language model) solely for English-language polishing and formatting of the manuscript. The scientific design, computational analyses, and all reported results were conceived and executed entirely by the author; the tool performed no analysis and contributed no data, code, or interpretation. No third party submitted this manuscript on the author's behalf.

## AUTHOR CONTRIBUTIONS

YYX: Conceptualization, Data curation, Formal analysis, Investigation, Methodology, Software, Visualization, Writing – original draft, Writing – review & editing. The author confirms being the sole contributor of this work and has approved it for publication.

## STATEMENTS AND DECLARATIONS

### Ethical considerations
Not applicable. This study is a purely computational re-analysis of publicly available, previously published, and de-identified datasets (Gene Expression Omnibus, STRING, DrugBank, and Protein Data Bank). No human participants were recruited or enrolled, no animal experiments were conducted by the author, and no identifiable individual data were generated or analyzed.

### Consent to participate
Not applicable.

### Consent for publication
Not applicable.

### Declaration of conflicting interests
The author declares that there are no potential conflicts of interest with respect to the research, authorship, and/or publication of this article.

### Funding
This study was conducted without external funding.

### Data availability
All datasets analyzed are publicly available: bulk RNA-seq GSE217932 and single-cell RNA-seq GSE134003 (Gene Expression Omnibus, https://www.ncbi.nlm.nih.gov/geo/); protein structures 7PZC, 6N9O, 1JDH, 3UA9 (RCSB PDB, https://www.rcsb.org/); protein–protein interaction data from STRING (https://string-db.org/); drug-target annotations from DrugBank (https://go.drugbank.com/). No new experimental data were generated. All processed data files and analysis code are deposited in the public, version-tagged reproducibility repository https://github.com/yyx-4113/astrocytic-wnt-pyroptosis-np. No "available on request" clause is used; all materials are openly accessible.

## REFERENCES

1. van Hecke O, Austin SK, Khan RA, Smith BH, Torrance N. Neuropathic pain in the general population: a systematic review of epidemiological studies. Pain. 2014;155(4):654-662. doi:10.1016/j.pain.2013.11.013
2. Finnerup NB, Attal N, Haroutounian S, et al. Pharmacotherapy for neuropathic pain in adults: a systematic review and meta-analysis. Lancet Neurol. 2015;14(2):162-173. doi:10.1016/S1474-4422(14)70251-0
3. Zhang YK, Huang ZJ, Liu S, Liu YP, Song AA, Song XJ. WNT signaling underlies the pathogenesis of neuropathic pain in rodents. J Clin Invest. 2013;123(5):2268-2286. doi:10.1172/JCI65364
4. Itokazu T, Hayano Y, Yamashita T, et al. Involvement of Wnt/β-catenin signaling in the development of neuropathic pain. Neurosci Res. 2014;80:18-23. doi:10.1016/j.neures.2013.12.002
5. Bortolotti P, Faure E, Kipnis E. Inflammasomes in tissue damages and immune disorders after trauma. Front Immunol. 2018;9:1900. doi:10.3389/fimmu.2018.01900
6. Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. Genome Biol. 2018;19(1):15. doi:10.1186/s13059-017-1382-0
7. Eastman P, Swails J, Chodera JD, et al. OpenMM 7: Rapid development of high performance algorithms for molecular dynamics. PLoS Comput Biol. 2017;13(7):e1005659. doi:10.1371/journal.pcbi.1005659
8. Kolde R, Laur S, Adler P, Vilo J. Robust rank aggregation for gene list integration and meta-analysis. Bioinformatics. 2012;28(4):573-580. doi:10.1093/bioinformatics/btr709
9. Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics. 2008;9:559. doi:10.1186/1471-2105-9-559
10. Langfelder P, Luo R, Oldham MC, Horvath S. Is my network module preserved and reproducible? PLoS Comput Biol. 2011;7(1):e1001057. doi:10.1371/journal.pcbi.1001057
11. Jin S, Guerrero-Juarez CF, Zhang L, et al. Inference and analysis of cell-cell communication using CellChat. Nat Commun. 2021;12(1):1088. doi:10.1038/s41467-021-21246-9
12. Cao J, Spielmann M, Qiu X, et al. The single-cell transcriptional landscape of mammalian organogenesis. Nature. 2019;566(7745):496-502. doi:10.1038/s41586-019-0969-x
13. Newman AM, Steen CB, Liu CL, et al. Determining cell type abundance and expression from bulk tissues with digital cytometry. Nat Biotechnol. 2019;37(7):773-782. doi:10.1038/s41587-019-0114-2
14. Szklarczyk D, Kirsch R, Koutrouli M, et al. The STRING database in 2023. Nucleic Acids Res. 2023;51(D1):D638-D646. doi:10.1093/nar/gkac1000
15. Castro-Mondragon JA, Riudavets-Puig R, Rauluseviciute I, et al. JASPAR 2022. Nucleic Acids Res. 2022;50(D1):D165-D173. doi:10.1093/nar/gkab1113
16. Wishart DS, Feunang YD, Guo AC, et al. DrugBank 5.0: a major update. Nucleic Acids Res. 2018;46(D1):D1074-D1082. doi:10.1093/nar/gkx1037
17. Margolin AA, Nemenman I, Basso K, et al. ARACNE. BMC Bioinformatics. 2006;7(Suppl 1):S7. doi:10.1186/1471-2105-7-S1-S7
18. GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science. 2020;369(6509):1318-1330. doi:10.1126/science.aaz1776
19. Vosa U, Claringbould A, Westra HJ, et al. Large-scale cis- and trans-eQTL analyses. Nat Genet. 2021;53(9):1300-1310. doi:10.1038/s41588-021-00913-z
20. Eastman P, Behara PK, Bhethanabotla VR, et al. SPICE, a generated dataset of 2.5M density functional theory and ab initio calculations for data-driven force fields. Sci Data. 2023;10:11. doi:10.1038/s41597-022-01882-6
21. Coll RC, Hill JR, Day CJ, et al. MCC950 directly targets the NLRP3 ATP-hydrolysis motif for inflammasome inhibition. Nat Chem Biol. 2019;15(6):556-559. doi:10.1038/s41589-019-0277-7
22. Hu JJ, Liu X, Xia S, et al. FDA-approved disulfiram inhibits pyroptosis by blocking gasdermin D pore formation. Nat Immunol. 2020;21(7):736-745. doi:10.1038/s41590-020-0669-6
23. Huang SM, Mishina YM, Liu S, et al. Tankyrase inhibition stabilizes axin and antagonizes Wnt signalling. Nature. 2009;461(7264):614-620. doi:10.1038/nature08356
24. Landrum G. RDKit: Open-Source Cheminformatics Software. 2024. https://www.rdkit.org/
25. Decosterd I, Woolf CJ. Spared nerve injury: an animal model of persistent peripheral neuropathic pain. Pain. 2000;87(2):149-158.
26. Ji RR, Donnelly CR, Nedergaard M. Astrocytes in chronic pain and itch. Nat Rev Neurosci. 2019;20(11):667-685. doi:10.1038/s41583-019-0218-1

---


---

## FIGURE LEGENDS

**Figure 0 – Graphical Abstract.** Schematic summary of the computational study design. (A) Transcriptomic analysis identifies Wnt/pyroptosis pathway co-expression in rat SNI/CCI spinal cord. (B) Single-cell RNA-seq localizes co-expression to reactive astrocytes. (D) Descriptive three-tier drug classification catalogs candidates for preclinical evaluation. All findings are computational predictions requiring experimental validation. (Panel C of the prior submission, which depicted human genetic instrumental-variable support, has been removed because the Mendelian randomization analysis was withdrawn; see text.)

**Figure 1 – Single-Dataset Differential Expression Overview (GSE217932).** (A) Volcano plot of GSE217932 (CCI vs sham, n=4 per group). Red points: |log2FC|>0.5 and nominal p<0.05. Wnt DEGs labeled. (B) Heatmap of top DEGs with hierarchical clustering. No second independent dataset was available; cross-dataset replication was not assessed.

**Figure 2 – Wnt/β-catenin Pathway Enrichment Analysis (GSE217932).** (A) GSEA plot for the Hallmark Wnt/β-catenin signaling pathway in GSE217932. The hallmark set was not enriched in the CCI (injury) direction (recomputed NES ≈ −1.6, nominal p ≈ 0.03, not significant after correction for the number of Hallmark sets tested); see text. (B) Heatmap of core Wnt genes. Single-cohort analysis; no cross-dataset comparison was performed.

**Figure 3 – Pyroptosis Gene Expression and WGCNA Co-Expression Modules.** (A) Heatmap of 16 pyroptosis gene expression values (z-score normalized) in GSE217932 (sham vs CCI). Genes meeting the exploratory DEG threshold marked with asterisks. (B) WGCNA module eigengene adjacency heatmap (modules, 5,000 genes, signed network, beta=6). Wnt/pyroptosis gene module assignments indicated. (C) Spearman correlation matrix of Wnt and pyroptosis gene expression in GSE217932. Module preservation across an independent dataset was not assessed (no independent cohort available).

**Figure 4 – Single-Cell RNA-seq Atlas of SNI Spinal Cord (full-transcriptome).** (A) UMAP embedding of GSE134003, colored by 14 annotated cell types. (B) Dot plot of Wnt pathway gene expression across cell types (dot size: percent expressed; color: mean expression). (C) Dot plot of pyroptosis gene expression across cell types, showing enrichment in microglia and macrophages. (D) Bar plot of pyroptosis module score by cell type (animal-level mean), showing microglial/macrophage predominance and low astrocyte scores. All scores computed across 14,700 genes; n=19 animals (SI: 10, SNI: 9).

**Figure 5 – CellChat Wnt Signaling and Pseudotime Trajectory Analysis.** (A) Circle plot of Wnt-Fzd ligand–receptor communication network in SNI condition. Edge width proportional to communication probability; node size proportional to cell number. (B) Pseudotime heatmap: expression of Wnt and pyroptosis genes across binned pseudotime (40 bins) along the astrocyte differentiation trajectory. (C) Violin plots of pseudotime distribution by cell type.

**Figure 6 – Cell Type Deconvolution by CIBERSORTx (GSE217932).** (A) Estimated cell type proportions in GSE217932 (sham vs CCI). Astrocyte proportion: 9.6% (sham) to 9.9% (CCI) (nominal p=0.343, not significant). (B) Signature matrix heatmap of marker gene expression across reference cell types. No cross-dataset correlation was computed (single bulk cohort).

**Figure 7 – Protein–Protein Interaction Network and TCF4 Binding-Site Prediction.** (A) High-confidence STRING network (combined score >0.7, 216 nodes, 189 edges). Hub genes highlighted by node size. (B) Multi-centrality hub gene ranking bar chart. Top 6 hubs: Il6, Ccl2, Wnt2, Nox4, Rela, Nfkb1. (C) JASPAR-predicted TCF4 (MA0830.1) binding sites in Ccl2, Il6, Jun, Nfkb1, Nlrp3, and Gsdmd promoters. (D) TCF4-centered sub-network of predicted transcriptional targets. Motif predictions are hypothetical.

**Figure 8 – Descriptive Drug Three-Tier Classification (No Statistical Testing).** (A) Three-tier drug classification overview: drug pills grouped by tier (green = approved, orange = investigational). (B) Drug-target bipartite network: 20 drugs, 21 unique targets, 33 edges (41 total nodes). Drug nodes colored by tier; target nodes in gray. (C) Detailed tier classification table with target annotations, approval status, and mechanism of action for all 20 compounds.

**Figure 9 – Gene Regulatory Network and Virtual Perturbation Analysis (Exploratory; n=12; No FDR).** (A) ARACNe-style GRN from GSE217932 (Spearman |rho|>0.75, DPI epsilon=0.15). Wnt genes in blue; pyroptosis genes in red; inferred regulators in gray. (B) Virtual KO impact scores by RWR perturbation. Tcf4 and Ctnnb1 identified as high-impact hubs. (C) Virtual OE seed-weight perturbation (+2.0 SD) of Wnt/β-catenin nodes: RWR propagation reaches pyroptosis genes (e.g., Il1b RWR = 0.011, Gsdmd RWR = 0.009, Nlrp3 RWR = 0.007). (D) Virtual KD (−2.0 SD) of Tcf4: reduced propagation to the same pyroptosis genes. All scores are RWR network-propagation outputs (hypothesis-generating only).

**Supplementary Figure S1 – Structural Analysis of Drug-Target Complexes (Literature Binding Energies; MD Not Completed).** (A) PDBFixer-remediated PDB structures: 7PZC (NLRP3 NACHT domain, 118 kDa), 6N9O (CTNNB1/β-catenin, 86 kDa), 1JDH (TNKS2 tankyrase-2, 71 kDa), 3UA9 (GSDMD pore, 53 kDa). (B) Published experimental binding affinities (Kd, nM): MCC950–NLRP3 Kd=5 nM (Coll et al., 2019); XAV939–TNKS2 Kd=46 nM (Huang et al., 2009). (C) Molecular-dynamics simulation status: 0/4 completed; no new binding-energy calculations were performed (literature values used).

**Supplementary Figure S2 – Reproducibility and Limitation Self-Audit.** A compliance matrix summarizing, for each analytical module (rows, keyed to the methods sections), the extent to which the study addressed five reproducibility and scientific-strictness dimensions (columns): sample size, independent replication, multiple-testing control, reporting of negative/null results, and code/data sharing. Green tick = explicitly addressed with adequate power; amber triangle = partially addressed or explicitly discussed as a limitation; red cross = not applicable or absent in this study. A bottom "Summary" row reports the weighted compliance rate per dimension (%), and a right-hand "Score" column reports the weighted compliance rate per module (%), making the scope and limits of the evidence explicit. Note: the study analyses a single bulk cohort and therefore scores "independent replication" as not applicable/absent for bulk-derived modules; the self-audit is a transparency instrument, not evidence of replication.

---

## SUPPLEMENTARY MATERIALS

**Supplementary Table S1 – WGCNA Module Assignments (GSE217932).** Module assignments for Wnt/pyroptosis genes across identified co-expression modules, with module sizes and exploratory interpretation. Module preservation across an independent dataset was not assessed (no independent bulk cohort available). Data file: MC5_WGCNA_Modules.csv.

**Supplementary Table S2 – Descriptive Drug Three-Tier Classification with Target Annotations.** Complete classification of 20 drugs with DrugBank target annotations, tier assignment, approval status, mechanism of action, and PubChem CID verification status. Data file: MC2_Drug_ThreeTier.csv. No statistical testing performed.

**Supplementary Table S3 – CellChat Wnt-Fzd Communication Probability Matrix (SNI Condition).** Row-normalized communication probabilities between 14 annotated cell types for Wnt-Fzd ligand–receptor pairs. Data file: Figure5_communication_SNI.csv.

All supplementary data files are provided in CSV format in the accompanying reproducibility package and are openly accessible (no "available on request" clause).
