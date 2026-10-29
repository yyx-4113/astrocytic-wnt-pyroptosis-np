# The Wnt/β-catenin transcriptional programme is not induced in neuropathic pain spinal cord: evidence from three independent rat cohorts

**Yongxin Yang**^1^

^1^ Department of Anesthesiology, The Second Affiliated Hospital of Fujian University of Traditional Chinese Medicine, Fuzhou, Fujian 350003, China

Correspondence: Yongxin Yang, 960856791@qq.com; ORCID 0009-0004-9698-6552

---

## Abstract

**Background.** Wnt/β-catenin signalling and pyroptosis have each been proposed as drivers of neuropathic pain, and their convergence in spinal astrocytes has been suggested. We asked whether the corresponding transcriptional programmes are detectable in spinal cord after nerve injury, and whether transcriptome-wide inference is achievable at the sample sizes on which these proposals rest.

**Results.** Three independent rat bulk RNA-seq cohorts were analysed through one pre-specified pipeline: GSE217932 (chronic constriction injury, 4 vs 4), GSE256472 (spared nerve injury, 3 vs 3) and GSE205906 (spared nerve injury, 3 vs 3). Transcripts were collapsed to gene level, tested by Welch t-test with Benjamini–Hochberg correction, summarised as 95% confidence intervals against a pre-specified ±0.5 log2 equivalence bound, and combined by random-effects meta-analysis; preranked gene set enrichment analysis (GSEA) gave a set-level readout. The assay was calibrated with *Atf3*, which rose 4.03 log2 (16.3-fold, *p* = 0.0044) and nevertheless reached only *q* = 0.466; no gene survived correction in the first two cohorts, and the single survivor in the third was zero-variance. Across 36 Wnt/β-catenin genes quantified in ≥2 cohorts, the largest pooled effect was +0.27 log2 (*Ccnd2*, 1.20-fold; 95% CI 0.00–0.53), 31 of 36 pooled intervals lay wholly within ±0.5, and GSEA gave NES −1.29 (*q* = 0.100), +1.02 (0.370) and +0.48 (0.931); the null held under four pre-processing variants (*p* = 0.10–0.25). Excluding the cohort whose own positive control failed, the largest induction was 1.23-fold. Pyroptosis transcripts were coordinately but modestly up-regulated in the chronic constriction cohort (NES +2.68, *q* < 0.001; 13 of 19 genes in the top 6% of the ranking, hypergeometric *p* = 2.4 × 10^−12^), whereas the rate-limiting nodes were unchanged: NLRP3 +0.01 (*p* = 0.97), NEK7 −0.01 (0.95), CASP1 +0.28 (0.11) and GSDMD +0.22 (0.06). This is a pattern of transcriptional priming, not of inflammasome execution, and it did not reproduce in the second cohort.

**Conclusions.** A Wnt/β-catenin transcriptional programme is not detectable in neuropathic pain spinal cord, and genome-wide false-discovery control is unattainable at n = 3–4. Pathway-activation claims in this literature therefore rest on nominal *p*-values, and require either larger samples or explicit equivalence framing.

**Keywords:** neuropathic pain; Wnt/β-catenin; pyroptosis; spinal cord; data reanalysis; equivalence testing; false discovery rate; sample size

---

## Background

Neuropathic pain arises from a lesion or disease of the somatosensory system and affects a substantial minority of patients after nerve injury, diabetes or chemotherapy [1]. Its treatment remains unsatisfactory, and the search for mechanism-based targets has generated a large literature built largely on rodent spinal cord and dorsal root ganglion tissue.

Wnt/β-catenin is one such candidate. Wnt ligands stabilise β-catenin, which partners with TCF/LEF factors to activate a transcriptional programme [2]. Several studies report that this programme is engaged after peripheral nerve injury: Wnt/β-catenin signalling was reported to contribute to the development of neuropathic pain [3], to up-regulate spinal VGLUT2 [4], to drive BDNF release from spinal microglia [5], and to mediate analgesic effects in chronic constriction injury through a CYP24A1-dependent route [6]. A parallel literature proposes pyroptosis — NLRP3 inflammasome assembly, caspase-1 activation and gasdermin D pore formation — as a mechanism of spinal cord injury and neuropathic pain [7–9], and a recent bioinformatic study nominated pyroptosis-related genes in a spared nerve injury model [10]. The two proposals have been combined into a single narrative in which Wnt/β-catenin and pyroptosis are co-activated in spinal astrocytes.

These conclusions rest almost entirely on bulk transcriptomic comparisons with three or four animals per group, in which differential expression is called at nominal *p* < 0.05. That inferential step is fragile for two well-characterised reasons. First, RNA-seq experiments of this size have low power for all but very large effects [11, 12], and false-discovery rate control across ~20,000 genes at n = 3–4 is close to unattainable [13]. Second, a non-significant result is not evidence of absence: without an equivalence bound, a null result cannot distinguish "no effect" from "not enough data" [14]. The two errors are complementary, and together they allow a literature in which both positive and negative claims are made without the statistical resolution to support either.

We therefore asked a narrower and answerable question. Using three independent rat spinal cord bulk RNA-seq cohorts generated in different laboratories, with two pain models (chronic constriction injury and spared nerve injury) and two quantification chemistries, and applying one identical pipeline to all three, we tested whether the Wnt/β-catenin transcriptional programme is detectable after nerve injury. We used (i) gene-level effect sizes with 95% confidence intervals against a pre-specified equivalence bound, (ii) random-effects meta-analysis, and (iii) preranked GSEA with fixed settings. We calibrated the assay with *Atf3*, a canonical injury-response transcription factor [15], so that a negative result could be attributed to biology rather than to insensitivity. We report two findings: the Wnt programme is not induced, and the transcriptome-wide false-discovery problem is severe enough to explain why the opposite conclusion was plausible in the first place.

---

## Methods

### Cohorts and pre-specified exclusions

Three rat (*Rattus norvegicus*) spinal cord bulk RNA-seq datasets were obtained from the Gene Expression Omnibus. Cohort C1 was GSE217932 (Luo et al. [16]), a chronic constriction injury (CCI) model [17]; four sham and four CCI libraries were used, and the four CCI + levo-tetrahydropalmatine libraries in the same series were excluded before any comparison. Cohort C2 was GSE256472, a spared nerve injury (SNI) model [18]; three control and three SNI libraries were used (no associated publication; the series description was used). Cohort C3 was GSE205906 (Chen et al. [19]), an SNI model; three sham (JA) and three SNI libraries were used, and the three paeoniflorin–liquiritin treated libraries were excluded. C1 and C2 provided raw counts; C3 provided processed FPKM values. Exclusions were applied at the level of sample columns and were identical for every analysis reported here.

A single-cell dataset, GSE134003 (Niehaus et al. [20]), was used only to quantify how many astrocytes and microglia the reference atlas contributes, and not for any expression comparison.

### Processing and gene-level aggregation

Counts (C1, C2) and FPKM values (C3) were transformed as log2(x + 1). Transcripts were collapsed to gene level by retaining, for each gene symbol, the transcript with the highest mean log2 expression across that cohort. This rule is deterministic and independent of the tested contrast; it avoids the selection bias inherent in choosing the transcript with the smallest *p*-value. All subsequent steps were identical across cohorts.

### Differential expression and equivalence testing

Each gene was tested with a two-sample Welch t-test on log2 values, which does not assume equal variances. The 95% confidence interval of the log2 fold change was constructed from the Welch standard error and Welch–Satterthwaite degrees of freedom. *P*-values were corrected across all genes tested in that cohort by the Benjamini–Hochberg procedure [21].

Equivalence was assessed against a bound of ±0.5 log2 (1.41-fold), fixed in the analysis code before any cohort was tested. A gene was recorded as equivalent to no change when its entire 95% confidence interval lay within ±0.5. This bound was chosen as the smallest effect that the field routinely interprets as pathway induction; it is a judgement about biological salience, not a statistical constant, and results are reported in a way that allows the reader to apply a different bound.

### Cross-cohort meta-analysis

For genes quantified in at least two cohorts, log2 fold changes were pooled by DerSimonian–Laird random-effects meta-analysis [22], with *I*² and τ² reported. A fixed-effect estimate was also computed and is used only in the sensitivity analysis described below.

### Gene set enrichment

Preranked GSEA [23] was run separately in each cohort. Genes were ranked by Welch t statistic (case versus control), with the statistic winsorised at ±50 to prevent near-zero-variance transcripts from occupying the top of the list. Genes below the cohort's own lower-quartile mean expression were removed before ranking; this filter is reported because an unfiltered t-based ranking is not numerically reproducible across re-runs, not because it changes any conclusion (see Results, sensitivity analyses). Enrichment scores were computed as the weighted Kolmogorov–Smirnov statistic with exponent *p* = 1. The null distribution was generated by 10,000 random gene sets of identical size drawn from the same ranked list; the ranking itself was never permuted, because phenotype permutation is not valid for preranked data with unequal group sizes. Normalised enrichment scores (NES) were obtained by dividing the observed score by the mean of the null scores on the matching sign side. Nominal *p* is one-sided in the observed direction, and *q*-values are Benjamini–Hochberg across the six gene sets tested in each cohort — not across MSigDB, and they should not be read as genome-wide Hallmark *q*-values.

Six gene sets were used. Two of them concern Wnt, and they are not interchangeable. Gene-level testing and equivalence assessment used a 37-gene curated Wnt/β-catenin panel spanning the core signal-transduction machinery (Frizzled receptors, the Axin–Apc–Gsk3b destruction complex, β-catenin itself and the Tcf/Lef factors), the ligands, and canonical transcriptional targets; the full list is in Additional file 1. Set-level testing used the HALLMARK_WNT_BETA_CATENIN collection as retrieved (38 symbols) [24]. The overlap between the two is only 14 genes, so they ask different questions: the curated panel tests whether the pathway's own components and its most responsive targets move, whereas the Hallmark set tests whether a target-gene programme is coordinately shifted. The remaining five were curated marker panels, not MSigDB collections: PYROPTOSIS_CORE (22 genes), INJURY_RESPONSE_POSCTRL (23 genes, the positive control), MICROGLIA_ACTIVATION (20), ASTROCYTE_REACTIVITY (19) and NFKB_CYTOKINE_SIGNALLING (22). Membership of every panel is given in Additional file 1.

Rank clustering was additionally tested without permutation. For a gene subset of size *k*, the number of members falling within the top 6% of the ranked list was compared with the hypergeometric expectation, which gives a check that does not depend on the permutation scheme.

### Sensitivity analyses

Five analyses probed the dependence of the conclusions on analysis choices. (i) Each of the 13 leading-edge pyroptosis genes was removed in turn and GSEA re-run (13 separate analyses). (ii) The 22-gene pyroptosis panel was split into an execution subset (NLRP3, NLRC4, AIM2, PYCARD, NEK7, CASP1, CASP4, GSDMD, GSDME; 9 genes) and a priming subset (NF-κB-inducible substrates and cytokines; 13 genes), and GSEA was run on each. (iii) GSEA was repeated under four pre-processing variants: the primary setting, log2 fold change instead of *t* as the ranking metric, expression filtering at the cohort median instead of the lower quartile, and no expression filter at all. (iv) Wnt effects were re-pooled using C1 and C2 only, excluding the cohort whose positive control failed. (v) For the Wnt panel the pooled estimate was also computed as a fixed-effect estimate, to separate a change in point estimate from a change in precision.

### Single-cell reference

GSE134003 was processed to cell-type labels as deposited. Cell numbers and the number of contributing animals were tabulated per cell type and condition. No differential expression was performed at single-cell level, because the astrocyte compartment proved too small to support it (Results).

### Software and reproducibility

Analyses were run in Python 3.11 with NumPy and SciPy. Scripts, gene panels, per-cohort result tables and the GSEA outputs are available at https://github.com/yyx-4113/astrocytic-wnt-pyroptosis-np (tag v1.1.0), which also carries a SHA-256 manifest for every data file and result table.

---

## Results

### The assay detects injury but correction removes everything, including the positive control

The three cohorts contributed 20 libraries in total (C1: 4 sham, 4 CCI; C2: 3 control, 3 SNI; C3: 3 sham, 3 SNI). Between 17,412 and 19,864 genes were testable per cohort after gene-level collapse (Table 1).

*Atf3*, the canonical injury-induced transcription factor, behaved as required of a positive control in the discovery cohort: +4.03 log2 (16.3-fold; 95% CI 2.07 to 5.98; *p* = 0.0044). It rose in the same direction in both other cohorts (+1.33, *p* = 0.0053; +1.03, *p* = 0.119), and the random-effects pooled estimate was +1.99 log2 (95% CI 0.64 to 3.33; *p* = 0.0039) with substantial heterogeneity (*I*² = 86%). The injury response is therefore present and detectable in this material.

Benjamini–Hochberg correction nonetheless removed it. *Atf3* reached *q* = 0.466 in C1, 0.700 in C2 and 0.761 in C3. Across the whole transcriptome, 1,371 of 19,009 genes reached nominal *p* < 0.05 in C1, 1,210 of 17,412 in C2 and 734 of 19,864 in C3; the numbers surviving *q* < 0.05 were 0, 0 and 1 (Table 1). The single survivor (LOC688442 in C3) is a zero-variance transcript with a fold change of 1.09 and a *t* statistic of 1.27 × 10^16^; it carries no biological information.

This is the first substantive result, and it is methodological rather than biological: at n = 3–4, a 16-fold change in a gene that is unambiguously regulated does not survive correction across ~19,000 tests. Any claim about a pathway in this literature that rests on nominal *p* is therefore a claim made below the detection threshold of the design.

### The Wnt/β-catenin transcriptional programme is not induced

Of the 37 genes in the curated Wnt panel, 36 were quantified in at least two cohorts (30 in all three) and could be pooled. The largest pooled effect in either direction was +0.27 log2 for *Ccnd2* (1.20-fold; 95% CI 0.00 to 0.53); its interval crosses the equivalence bound and it is the one gene in the panel for which equivalence cannot be claimed. *Jun* (+0.26), *Birc5* (+0.24) and *Cd44* (+0.19) followed. Thirty-one of the 36 pooled intervals lay wholly within ±0.5 log2, which excludes an induction of 1.41-fold or more for those genes (Fig. 1C).

The core pathway components were flat. *Ctnnb1* (β-catenin) itself was +0.08 log2 (pooled, three cohorts; *p* = 0.16), *Axin2* — the most responsive canonical target gene — was −0.02, *Lef1* −0.01, *Tcf4* −0.12 and *Gsk3b* −0.05. The largest decrease in the panel was *Wnt2* at −0.17, whose interval (−0.75 to 0.41) is too wide to support equivalence in either direction; no Wnt gene was equivalent to a decrease.

At the set level the result was the same. For HALLMARK_WNT_BETA_CATENIN (29 to 33 of its 38 symbols passing the expression filter, depending on cohort) GSEA gave NES −1.29 (*p* = 0.100, *q* = 0.100) in C1, +1.02 (*p* = 0.304, *q* = 0.370) in C2 and +0.48 (*p* = 0.931, *q* = 0.931) in C3 (Fig. 2A). None reached significance, and the sign reversed between cohorts. In C1 the negative score was produced by mild depletion of *AXIN2*, *TCF7*, *TCF7L1*, *LEF1*, *PAX2*, *WIF1* and *CXXC4* (leading edge 13 of 33 genes): the direction is depletion, not activation (Fig. 2B).

Because the Wnt null is the paper's main claim, its dependence on pre-processing was examined directly. Across four ranking and filtering variants the NES ranged from −1.29 to −1.01 in C1 with *p* between 0.10 and 0.25, and never approached significance in any variant (Additional file 2). The null is not an artefact of the *t*-statistic ranking, of the expression filter, or of winsorising.

### Excluding the cohort with a failed positive control does not change the Wnt result

C3 did not behave as a replication cohort. Its injury programme showed a *negative* enrichment score (INJURY_RESPONSE_POSCTRL NES −1.30, *q* = 0.018), and its glial markers were at or below sham levels (*Aif1* 0.78-fold, *C1qa* 0.89-fold, *Gfap* 0.97-fold). Its *Atf3* response was 2.27-fold and did not reach significance (*p* = 0.119). A cohort in which the canonical injury programme runs in the opposite direction cannot function as confirmation of anything, and we therefore re-pooled the Wnt panel using C1 and C2 only.

The point estimates barely moved: 31 genes were poolable, the largest induction was *Tcf3* at +0.30 log2 (1.23-fold), and the only gene whose point estimate exceeded 0.5 was *Wnt2*, which moved *down* (−0.50). The number of genes meeting the equivalence criterion fell from 31 of 36 to 20 of 31, but this reflects wider confidence intervals from having one fewer cohort, not a larger effect: no point estimate crossed the bound. This distinction matters, because a fall in the equivalence count is easily misread as instability of the result.

C3 is retained in the manuscript for two purposes only: as a demonstration that a public dataset can fail its own positive control, and as the source of the single zero-variance gene that passed correction. It is not counted as independent support.

### Pyroptosis transcripts are primed but the execution machinery is not engaged

The pyroptosis panel behaved differently from Wnt, and not in the direction we expected. In C1 the panel was strongly enriched (NES +2.68, *p* = 1 × 10^−4^, *q* < 0.001), with 13 of 19 quantified members in the leading edge, all within the top 844 of 14,257 ranked genes. Rank clustering was far beyond chance: 13 of 19 in the top 6% (expected 1.14; hypergeometric *p* = 2.4 × 10^−12^). The enrichment was not carried by one or two genes — removing each leading-edge gene in turn reduced NES from 2.68 to no lower than 2.60, and all 13 leave-one-out analyses remained significant (Additional file 2).

The composition of that signal is the informative part. The leading edge consisted of *Pycard*, *Nlrc4*, *Tlr4*, *Il6*, *Il1b*, *Gsdmd*, *Casp1*, *Txnip*, *Casp4*, *Il18*, *Rela*, *Aim2* and *Casp8*. *Nlrp3* was absent from it. Gene by gene, the rate-limiting nodes were unchanged: NLRP3 +0.01 log2 (*p* = 0.97), NEK7 −0.01 (0.95), CASP1 +0.28 (0.11), GSDMD +0.22 (0.06) and GSDME +0.06 (0.69), against +0.58 for *Pycard* (0.003), +1.73 for *Il1b* (0.063) and +1.98 for *Il6* (0.021).

We tested whether the set-level statistic could distinguish the two stages of the pathway by splitting the panel into an execution subset and a priming subset. It could not: NES was +2.34 for the execution subset and +2.12 for the priming subset, a difference with no statistical content. The set-level statistic therefore does not support a claim about which stage is engaged; the gene-level effect sizes do. Eight of nine execution genes nevertheless fell in the top 6% of the ranking (hypergeometric *p* = 1.4 × 10^−9^), which is consistent with coordinated, low-amplitude induction of the whole panel rather than selective activation of the terminal step.

The reading we adopt is priming without execution: the NF-κB-dependent inducible components are modestly and coordinately up-regulated, while the sensor (NLRP3) and the effectors that require post-translational activation (CASP1, GSDMD) are not. This is a statement about transcription, and it cannot exclude inflammasome assembly or gasdermin cleavage, both of which occur after and independently of transcription.

The pattern did not replicate. In C2 the pyroptosis panel gave NES +0.76 (*p* = 0.514) and in C3 −1.11 (*p* = 0.033, i.e. opposite in sign). The priming signal is therefore a property of the CCI cohort, not a general feature of neuropathic pain models in this dataset collection.

### The single-cell reference cannot localise the signal to astrocytes

GSE134003 contains 18,856 cells, of which 735 are microglia (385 sham, 350 SNI; 10 and 9 animals) and 40 are astrocytes (27 sham, 13 SNI; 9 and 5 animals). The astrocyte compartment averages 2.5 cells per animal in the injured group. Single-cell atlases of the injured spinal cord have documented glial heterogeneity and myelin pathology after injury [25, 26], but none of them resolves the astrocyte compartment at the depth required to attribute a bulk signal. No cell-type attribution of the pyroptosis signal is possible at this depth, and we make none. The claim that Wnt and pyroptosis converge in spinal astrocytes cannot be tested in the available single-cell material irrespective of what the bulk data show.

---

## Discussion

### Principal findings

In three independent rat spinal cord cohorts, the Wnt/β-catenin transcriptional programme was not induced after nerve injury. The largest pooled effect across 36 genes was 1.20-fold, 31 of 36 intervals excluded an effect of 1.41-fold or more, and set-level enrichment was non-significant in all three cohorts and under all four pre-processing variants. At the same time, a gene that is unambiguously induced — *Atf3*, 16-fold — failed correction, and no gene in two of three cohorts survived false-discovery control. The second finding explains the first: in this literature, a pathway can be reported as "activated" on evidence that would not survive the multiple-testing burden the same data impose.

### What a transcriptional null does and does not exclude

β-catenin activation is predominantly post-translational: destruction-complex disassembly, phosphorylation, nuclear translocation and co-factor recruitment all occur without changes in transcript abundance [2]. A transcriptomic null therefore cannot exclude Wnt signalling activity, and we do not claim that it does. What it excludes is the specific and commonly asserted proposition that nerve injury induces a Wnt/β-catenin *transcriptional programme* in spinal cord — the coordinated induction of canonical target genes such as *Axin2*, *Ccnd1*, *Lef1* and *Tcf7*. That proposition is testable in bulk tissue, and it fails: *Axin2* moved −0.02 log2, and the β-catenin transcript itself moved +0.08.

The distinction is not pedantic. Three of the studies that motivated this analysis report target-gene induction or pathway-dependence in tissue-level assays [3–6]; if the programme were engaged, at least the most responsive canonical targets should move. They do not, in any of three cohorts, at effect sizes the design can resolve.

### Priming without execution, and why the distinction survives

The pyroptosis result initially appeared to contradict the Wnt result — an enriched inflammatory panel alongside a flat Wnt panel — and it is worth being explicit about how it was handled. The set-level enrichment is real and robust: it survives leave-one-out at every leading-edge gene, it is present under all four pre-processing variants (NES 1.83 to 2.68), and it is confirmed by a permutation-free hypergeometric test.

It is also, on inspection, not about pyroptosis execution. The genes that carry the signal are the inducible substrates — *Il1b*, *Il6*, *Pycard*, *Tlr4* — while the sensor that initiates assembly (*Nlrp3*) is unchanged, and the effectors whose activity depends on cleavage (*Casp1*, *Gsdmd*) show effects of 1.16–1.21-fold that do not reach nominal significance. We attempted to formalise this by splitting the panel into execution and priming subsets, and report that the split failed to separate them statistically (NES +2.34 versus +2.12). The conclusion rests instead on the gene-level estimates and on the composition of the leading edge. Presenting the split is necessary precisely because it was the obvious way to support the claim, and it did not work.

The biological content is modest but coherent: after CCI, spinal cord tissue mounts a low-amplitude, coordinated NF-κB-driven inflammatory transcriptional response without induction of the inflammasome sensor. Whether this reflects priming that requires a second signal, or simply the macrophage/microglial infiltration that accompanies nerve injury, cannot be resolved here.

### Why the literature contains the opposite conclusion

The two errors documented above are mutually reinforcing. A study that tests 19,000 genes at n = 3–4 and reports nominal *p* < 0.05 will find roughly 1,000 "differentially expressed" genes by chance alone; among 20–40 genes of a candidate pathway, several will fall below 0.05 and the pathway will appear activated. Our C1 data illustrate this exactly: 1,371 genes at nominal *p* < 0.05, zero after correction, and a positive control with *q* = 0.47.

The consequence is not that individual studies are dishonest, but that the field's inferential standard is below its design's resolution. Two remedies are available and both are cheap: report effect sizes with confidence intervals against a stated equivalence bound, and report set-level statistics with the pre-processing settings that produced them. The first converts "not significant" into "excludes effects above 1.41-fold"; the second makes the analysis reproducible. Neither requires new experiments.

Published estimates of the required replicate number for RNA-seq [11] and of the consequences of underpowered neuroscience [12] have been available for a decade; the literature on neuropathic pain transcriptomics has largely not adopted them. Our contribution is to show what happens when they are applied to a specific, widely cited mechanistic claim.

### Limitations

First, all data are transcriptomic. No protein, phosphorylation or localisation data were generated, and the post-translational events that define Wnt and pyroptosis activity are invisible to this design. The conclusions are correspondingly restricted to transcription.

Second, the cohorts are small (n = 3–4 per group), and this is the object of study as much as a limitation. Power is sufficient to detect a 16-fold change and insufficient to detect a 1.4-fold change after correction; genes whose intervals exceed the equivalence bound are indeterminate rather than null.

Third, C3 failed its own positive control and its injury programme ran in the opposite direction to the other two cohorts. We report it, downgrade it, and do not count it as replication. Whether this reflects a different injury time point, a pre-activated sham baseline, or a processing difference cannot be determined from the deposited data.

Fourth, most panels are curated rather than drawn from MSigDB: the 37-gene Wnt panel used for gene-level testing, the pyroptosis panel and the four cell-type and signalling panels were assembled by hand, and only the set-level Wnt analysis uses a Hallmark collection. The *q*-values apply to the six sets tested, not to a genome-wide collection. A full Hallmark scan would give different, larger *q*-values and was not performed. Because the curated and Hallmark Wnt sets overlap at only 14 genes, the two analyses are not two views of one object, and we report them as separate lines of evidence.

Fifth, the equivalence bound of ±0.5 log2 is a convention, not a measurement. A reader who prefers ±0.3 would count fewer genes as equivalent; the direction of the conclusion would not change, because the point estimates are all below 0.3 as well.

Sixth, we analysed one contrast per dataset. Other contrasts, time points and spinal levels exist in the same series and were outside the pre-specified scope.

### Conclusions

The Wnt/β-catenin transcriptional programme is not induced in rat spinal cord after nerve injury, across three independent cohorts and four analysis variants. The pyroptosis literature's central observation — coordinated induction of inflammasome-associated transcripts — is reproducible in the CCI cohort, but it reflects priming of NF-κB-dependent substrates rather than engagement of the execution machinery, and it does not replicate in a second model. More consequentially, the design on which these claims rest cannot support transcriptome-wide inference at all: a 16-fold change in a canonical injury gene fails correction. Progress in this literature requires larger samples or explicit equivalence framing, and the second is available immediately.

---

## Declarations

**Ethics approval and consent to participate.** Not applicable. This study used only publicly available, de-identified rodent data deposited in the Gene Expression Omnibus; no animals were used and no human participants were involved.

**Consent for publication.** Not applicable.

**Availability of data and materials.** All datasets analysed are publicly available: GSE217932, GSE256472, GSE205906 and GSE134003 (Gene Expression Omnibus, https://www.ncbi.nlm.nih.gov/geo/). Analysis code, gene panels, per-cohort result tables and GSEA outputs are available at https://github.com/yyx-4113/astrocytic-wnt-pyroptosis-np under tag v1.1.0, together with a SHA-256 manifest of all files.

**Competing interests.** The author declares no competing interests.

**Funding.** This work received no specific grant from any funding agency in the public, commercial or not-for-profit sectors.

**Authors' contributions.** YY conceived the study, designed the analysis, performed all analyses, interpreted the results and wrote the manuscript.

**Acknowledgements.** Not applicable.

---

## Figure legends

**Fig. 1** Three independent rat spinal cord cohorts, one pipeline, and no induction of the Wnt/β-catenin transcriptional programme. (**A**) Pre-specified design: three bulk RNA-seq cohorts analysed with one pipeline; intervention arms (CCI + levo-tetrahydropalmatine, n = 4; paeoniflorin–liquiritin, n = 3) excluded before any comparison. (**B**) Positive control. *Atf3* rose 16.3-fold in the discovery cohort (nominal *p* = 0.0044), establishing that injury-induced transcription is detectable at this sample size; the same gene nevertheless failed Benjamini–Hochberg correction (*q* = 0.466), as it did in the two further cohorts (*q* = 0.700 and 0.761). (**C**) All 36 genes of the 37-gene curated Wnt/β-catenin panel that were quantified in at least two cohorts. Points show per-cohort log2 fold change (C1 red, C2 blue, C3 green); diamonds and bars show the DerSimonian–Laird pooled estimate and its 95% confidence interval. The grey band marks the pre-specified equivalence bound of ±0.5 log2 (1.41-fold); 31 of 36 pooled intervals lie wholly within it and the largest pooled effect is 0.27 log2. (**D**) Genes reaching nominal *p* < 0.05 versus *q* < 0.05 in each cohort. The single gene passing correction (C3) is a zero-variance transcript (LOC688442, 1.09-fold), so no gene shows a biologically interpretable effect after correction in any cohort. *Note*: excluding C3, the cohort whose positive control did not reach significance, widens the pooled intervals (20 of 31 remaining genes still equivalent) but changes no point estimate; the largest induction among the remaining genes is 1.23-fold.

**Fig. 2** Preranked gene set enrichment in three cohorts. (**A**) Normalised enrichment scores for the six gene sets, per cohort. Ranking metric: Welch *t* (case versus control), winsorised at ±50, genes below the cohort's lower-quartile mean expression removed; 10,000 gene-set permutations; *q* values are Benjamini–Hochberg across the six sets tested in each cohort. HALLMARK_WNT_BETA_CATENIN is non-significant and direction-inconsistent across all three cohorts. (**B**) Running enrichment score for the Wnt panel in the C1 cohort (NES −1.29, *p* = 0.100). The negative score is generated by mild depletion of *Axin2*, *Tcf7*, *Tcf7l1*, *Lef1*, *Pax2*, *Wif1* and *Cxxc4*, whose positions are marked; the direction is depletion, not activation.

---

## Additional files

**Additional file 1** Supplementary tables (XLSX). Table S1: per-cohort gene-level statistics for all analysed genes. Table S2: random-effects and fixed-effect pooled estimates for the Wnt, pyroptosis and positive-control panels. Table S3: GSEA output for all cohorts, sets and pre-processing variants.

**Additional file 2** Supplementary figures (PDF). Figure S1: dependence of the Wnt and pyroptosis enrichment scores on pre-processing settings. Figure S2: leave-one-out analysis of the 13 leading-edge pyroptosis genes.

---

## References

1. Colloca L, Ludman T, Bouhassira D, Baron R, Dickenson AH, Yarnitsky D, et al. Neuropathic pain. Nat Rev Dis Primers. 2017;3:17002. doi:10.1038/nrdp.2017.2
2. Nusse R, Clevers H. Wnt/β-Catenin Signaling, Disease, and Emerging Therapeutic Modalities. Cell. 2017;169(6):985-999. doi:10.1016/j.cell.2017.05.016
3. Itokazu T, Hayano Y, Takahashi R, Yamashita T. Involvement of Wnt/β-catenin signaling in the development of neuropathic pain. Neurosci Res. 2014;79:34-40. doi:10.1016/j.neures.2013.12.002
4. Zhang ZL, Yu G, Peng J, Wang HB, Li YL, Liang XN, et al. Wnt1/β-catenin signaling upregulates spinal VGLUT2 expression to control neuropathic pain in mice. Neuropharmacology. 2020;164:107869. doi:10.1016/j.neuropharm.2019.107869
5. Zhou X, Tao L, Zhao M, Wu S, Obeng E, Wang D, et al. Wnt/β-catenin signaling regulates brain-derived neurotrophic factor release from spinal microglia to mediate HIV(1) gp120-induced neuropathic pain. Mol Pain. 2020;16:1744806920922100. doi:10.1177/1744806920922100
6. Liao J, Liu J, Long G, Lv X. MiR-30b-5p attenuates neuropathic pain by the CYP24A1-Wnt/β-catenin signaling in CCI rats. Exp Brain Res. 2022;240(1):263-277. doi:10.1007/s00221-021-06253-y
7. Liu P, Cheng J, Ma S, Zhou J. Paeoniflorin attenuates chronic constriction injury-induced neuropathic pain by suppressing spinal NLRP3 inflammasome activation. Inflammopharmacology. 2020;28(6):1495-1508. doi:10.1007/s10787-020-00737-z
8. Zhang X, Miao Y, Li Z, Xu H, Niu Z. ACVR1 drives neuropathic pain by regulating NLRP3-Induced neuronal pyroptosis through the p38 and Smad1/5/8 pathways. Neuropharmacology. 2025;274:110469. doi:10.1016/j.neuropharm.2025.110469
9. Deng D, Ma L, Shen J, Huang L, Zhang T, Yang X, et al. PD-1 Attenuates Neuropathic Pain by Ameliorating NLRP3 Inflammasome-Mediated Microglia Pyroptosis. Mol Neurobiol. 2025;62(12):15478-15493. doi:10.1007/s12035-025-05225-5
10. Li W, Hu Z, Lin P, He L, Liu R. Identification of Pyroptosis-Related Genes in Male Rats with Spared Nerve Injury-Induced Neuropathic Pain. J Pain Res. 2025;18:7029-7041. doi:10.2147/JPR.S568630
11. Schurch NJ, Schofield P, Gierliński M, Cole C, Sherstnev A, Singh V, et al. How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use? RNA. 2016;22(6):839-51. doi:10.1261/rna.053959.115
12. Button KS, Ioannidis JP, Mokrysz C, Nosek BA, Flint J, Robinson ES, et al. Power failure: why small sample size undermines the reliability of neuroscience. Nat Rev Neurosci. 2013;14(5):365-76. doi:10.1038/nrn3475
13. Pawitan Y, Michiels S, Koscielny S, Gusnanto A, Ploner A. False discovery rate, sensitivity and sample size for microarray studies. Bioinformatics. 2005;21(13):3017-24. doi:10.1093/bioinformatics/bti448
14. Lakens D. Equivalence Tests: A Practical Primer for t Tests, Correlations, and Meta-Analyses. Soc Psychol Personal Sci. 2017;8(4):355-362. doi:10.1177/1948550617697177
15. Renthal W, Tochitsky I, Yang L, Cheng YC, Li E, Kawaguchi R, et al. Transcriptional Reprogramming of Distinct Peripheral Sensory Neuron Subtypes after Axonal Injury. Neuron. 2020;108(1):128-144.e9. doi:10.1016/j.neuron.2020.07.026
16. Luo Z, Zhang Z, Li P, Yi M, Luo A, Zeng H, et al. The analgesic effect and mechanism of the active components screening from Corydalis yanhusuo by P2X3 receptors. J Ethnopharmacol. 2025;337(Pt 2):118989. doi:10.1016/j.jep.2024.118989
17. Bennett GJ, Xie YK. A peripheral mononeuropathy in rat that produces disorders of pain sensation like those seen in man. Pain. 1988;33(1):87-107. doi:10.1016/0304-3959(88)90209-6
18. Decosterd I, Woolf CJ. Spared nerve injury: an animal model of persistent peripheral neuropathic pain. Pain. 2000;87(2):149-158. doi:10.1016/S0304-3959(00)00276-1
19. Chen YY, Feng LM, Xu DQ, Yue SJ, Fu RJ, Zhang MM, et al. Combination of paeoniflorin and liquiritin alleviates neuropathic pain by lipid metabolism and calcium signaling coordination. Front Pharmacol. 2022;13:944386. doi:10.3389/fphar.2022.944386
20. Niehaus JK, Taylor-Blake B, Loo L, Simon JM, Zylka MJ. Spinal macrophages resolve nociceptive hypersensitivity after peripheral injury. Neuron. 2021;109(8):1274-1282.e6. doi:10.1016/j.neuron.2021.02.018
21. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B Stat Methodol. 1995;57(1):289-300. doi:10.1111/j.2517-6161.1995.tb02031.x
22. DerSimonian R, Laird N. Meta-analysis in clinical trials. Control Clin Trials. 1986;7(3):177-88. doi:10.1016/0197-2456(86)90046-2
23. Subramanian A, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, Gillette MA, et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proc Natl Acad Sci U S A. 2005;102(43):15545-50. doi:10.1073/pnas.0506580102
24. Liberzon A, Birger C, Thorvaldsdóttir H, Ghandi M, Mesirov JP, Tamayo P. The Molecular Signatures Database (MSigDB) hallmark gene set collection. Cell Syst. 2015;1(6):417-425. doi:10.1016/j.cels.2015.12.004
25. Li D, Yang K, Li J, Xu X, Gong L, Yue S, et al. Single-cell sequencing reveals glial cell involvement in development of neuropathic pain via myelin sheath lesion formation in the spinal cord. J Neuroinflammation. 2024;21(1):213. doi:10.1186/s12974-024-03207-3
26. Milich LM, Choi JS, Ryan C, Cerqueira SR, Benavides S, Yahn SL, et al. Single-cell analysis of the cellular heterogeneity and interactions in the injured mouse spinal cord. J Exp Med. 2021;218(8):e20210040. doi:10.1084/jem.20210040

---

## Table 1

**Table 1** Cohort characteristics and transcriptome-wide results

| | C1 GSE217932 | C2 GSE256472 | C3 GSE205906 |
|---|---|---|---|
| Model | CCI | SNI | SNI |
| Groups (control vs injury) | 4 vs 4 | 3 vs 3 | 3 vs 3 |
| Quantification | raw counts | raw counts | processed FPKM |
| Excluded arms | CCI + THP (n = 4) | none | PLM (n = 3) |
| Genes tested | 19,009 | 17,412 | 19,864 |
| Nominal *p* < 0.05 | 1,371 | 1,210 | 734 |
| *q* < 0.05 | 0 | 0 | 1^a^ |
| Intervals within ±0.5 log2 | 9,149 (48.1%) | 3,373 (19.4%) | 8,902 (44.8%) |
| *Atf3* log2FC (95% CI) | +4.03 (2.07 to 5.98) | +1.33 (0.75 to 1.91) | +1.03 (−0.43 to 2.48) |
| *Atf3* *p* / *q* | 0.0044 / 0.466 | 0.0053 / 0.700 | 0.119 / 0.761 |
| Wnt NES (*q*) | −1.29 (0.100) | +1.02 (0.370) | +0.48 (0.931) |
| Pyroptosis NES (*q*) | +2.68 (<0.001) | +0.76 (0.514) | −1.11 (0.040) |

^a^ Zero-variance transcript (LOC688442), 1.09-fold change; not biologically interpretable.
