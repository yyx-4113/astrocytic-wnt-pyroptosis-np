# No induction of a Wnt/β-catenin transcriptional programme in bulk rat spinal cord after peripheral nerve injury: a three-cohort reanalysis with equivalence testing

**Yongxin Yang**^1^

^1^ Department of Anesthesiology, The Second Affiliated Hospital of Fujian University of Traditional Chinese Medicine, Fuzhou, Fujian 350003, China

Correspondence: Yongxin Yang, 960856791@qq.com; ORCID 0009-0004-9698-6552

---

## Abstract

**Background.** Wnt/β-catenin signalling and pyroptosis have each been proposed as drivers of neuropathic pain, and their convergence in spinal astrocytes has been suggested. These proposals rest largely on bulk transcriptomic comparisons of three or four animals per group called at nominal *p* < 0.05. We asked whether the corresponding transcriptional programmes are detectable in spinal cord after nerve injury, and what such a design can resolve.

**Results.** Three independent rat bulk RNA-seq cohorts were analysed through one pipeline fixed before any cohort result was read: GSE217932 (chronic constriction injury, 4 vs 4), GSE256472 (spared nerve injury, 3 vs 3) and GSE205906 (spared nerve injury, 3 vs 3), using equivalence testing against a ±0.5 log2 bound and random-effects pooling. The assay was calibrated with *Atf3*, which rose 4.03 log2 (16.3-fold, *p* = 0.0044) yet reached only *q* = 0.466; no gene survived correction in two cohorts, and the single survivor in the third was a near-zero-variance artefact. Three library-normalisation regimes, with and without independent filtering, left this unchanged (0–6 of ~19,000 genes), so it is a property of the design, not the pipeline: after correction the smallest detectable effect is 5.6-fold in the largest cohort. Across 36 Wnt/β-catenin genes the largest pooled effect was +0.27 log2 (1.20-fold); 31 of 36 pooled intervals lay within ±0.5 under the normal-quantile interval but only 18 of 36 under the Hartung–Knapp correction appropriate at k = 2–3, the remainder indeterminate. Pyroptosis transcripts were coordinately up-regulated in the chronic constriction cohort only (NES +2.68, *p* ≤ 1 × 10^−4^), with NLRP3 unchanged (+0.01, *p* = 0.97) and the effectors showing 1.16–1.21-fold changes whose intervals do not exclude induction.

**Conclusions.** A spinal-cord-wide Wnt/β-catenin transcriptional programme is not detectable in bulk rat tissue after nerve injury. Because all three cohorts sampled whole spinal segments, a dorsal-horn-restricted programme would be diluted into exactly this range and is not excluded. Transcriptome-wide inference at n = 3–4 lacks the power to resolve effects of the size routinely read as pathway induction.

**Keywords:** neuropathic pain; Wnt/β-catenin; pyroptosis; spinal cord; data reanalysis; equivalence testing; statistical power; sample size

---

## Background

Neuropathic pain arises from a lesion or disease of the somatosensory system and affects a substantial minority of patients after nerve injury, diabetes or chemotherapy [1]. Its treatment remains unsatisfactory, and the search for mechanism-based targets has generated a large literature built largely on rodent spinal cord and dorsal root ganglion tissue.

Wnt/β-catenin is one such candidate. Wnt ligands stabilise β-catenin, which partners with TCF/LEF factors to activate a transcriptional programme [2]. Several studies report that this programme is engaged after nerve injury, most of them in mice: Wnt/β-catenin signalling was reported to contribute to neuropathic pain after partial sciatic nerve ligation [3] and to up-regulate spinal VGLUT2 after spared nerve injury [4], both localised to the dorsal horn. In a centrally administered gp120 model of HIV-associated neuropathic pain the same pathway was reported to drive BDNF release from spinal microglia [5], and in rat chronic constriction injury a CYP24A1-dependent Wnt/β-catenin route was reported to mediate analgesia [6]. Only the last of these is in the species and model combination analysed here. A parallel literature proposes pyroptosis — NLRP3 inflammasome assembly, caspase-1 activation and gasdermin D pore formation — as a mechanism of neuropathic pain [7–9], and a recent study nominated pyroptosis-related genes in rat spared nerve injury [10]. The two proposals have been combined into a single narrative in which Wnt/β-catenin and pyroptosis are co-activated in spinal astrocytes.

To our knowledge this proposition has not previously been tested as a hypothesis able to fail. The supporting reports are pathway interpretations of nominal *p*-values [3–6]; none pre-states an equivalence bound, none reports interval estimates, and none analyses more than one cohort. Two inferential problems make such reports fragile, and they are distinct. First, experiments of this size have low power: at three replicates, current tools recover only 20–40% of truly differentially expressed genes, rising above 85% only for changes greater than fourfold [11]. This is a loss of sensitivity, not a failure of false-discovery control — the same tools do hold the false discovery rate near 5% at any replicate number [11] — and low power in small animal studies is a general problem, not one specific to transcriptomics [12]. Second, a non-significant result is not evidence of absence: without an equivalence bound, a null result cannot distinguish "no effect" from "not enough data" [13–14]. Low power and uncorrected nominal testing reinforce each other, and together they permit a literature in which both positive and negative claims are made without the resolution to support either.

We therefore asked a narrower question. Using three independent rat spinal cord bulk RNA-seq cohorts generated in different laboratories, with two nerve-injury models, and applying one identical pipeline to all three, we tested whether a Wnt/β-catenin transcriptional programme is detectable after nerve injury. We used gene-level effect sizes with confidence intervals against a pre-specified equivalence bound, random-effects pooling, and preranked gene set enrichment. We calibrated the assay with *Atf3*, a canonical injury-response transcription factor [15], so that a negative result could be attributed to biology rather than to insensitivity. This analysis was not formally pre-registered; the pipeline was fixed in code before any cohort result was read and is archived with a checksum manifest at tag v1.1.4.

---

## Methods

### Cohorts and pre-specified exclusions

Three rat (*Rattus norvegicus*) spinal cord bulk RNA-seq datasets were obtained from the Gene Expression Omnibus [16]. Cohort C1 was GSE217932 (Wu et al. [17]) [18], a chronic constriction injury (CCI) model [19]; four sham and four CCI libraries were used, and the four CCI + levo-tetrahydropalmatine libraries were excluded (GSM6730109–6730112). Cohort C2 was GSE256472 [20], a spared nerve injury (SNI) model [21]; three control and three SNI libraries were used, none excluded; no publication is associated with the series (GEO returns an empty PubMed list as of 2026-10-30) and the series description was used. Cohort C3 was GSE205906 (Chen et al. [22]) [23], an SNI model; three sham (JA) and three SNI libraries were used and the three paeoniflorin–liquiritin libraries were excluded (GSM6235298–6235300). C1 and C2 provided raw counts (the C2 deposit also contains FPKM columns; counts were used); C3 provided processed normalised abundance, recorded in GEO as RPKM.

Exclusions were applied at the level of sample columns and were identical for every analysis reported. Segmental coverage and post-injury interval were taken from the GEO records and are reported in Table 1; they are not uniform across cohorts, and this is reported rather than corrected for. All three cohorts sampled whole spinal cord segments; no deposited record describes dorsal-horn or dorsal-root-ganglion microdissection in any of the three series.

A mouse (*Mus musculus*) single-cell dataset, GSE134003 (Niehaus et al. [24]) [25], was used only to quantify how deeply astrocytes and microglia are sampled in current spinal atlases, not for any expression comparison.

### Processing and gene-level aggregation

Counts (C1, C2) and RPKM values (C3) were transformed as log2(x + 1). Transcripts were collapsed to gene level by retaining, for each gene symbol, the transcript with the highest mean log2 expression across that cohort. This rule is deterministic and independent of the tested contrast, and avoids the selection bias inherent in choosing the transcript with the smallest *p*-value. All subsequent steps were identical across cohorts.

### Differential expression, equivalence testing and pooling

Each gene was tested with a two-sample Welch *t*-test on log2 values. The 95% confidence interval of the log2 fold change was constructed from the Welch standard error and Welch–Satterthwaite degrees of freedom. *P*-values were corrected across all genes tested in that cohort by the Benjamini–Hochberg procedure [26].

Equivalence was assessed against a bound of ±0.5 log2 (1.41-fold), fixed in the analysis code before any cohort was tested. A gene was recorded as equivalent to no change when its entire 95% confidence interval lay within ±0.5. The bound is a judgement about biological salience, not a statistical constant, and results are reported so that a reader can apply a different one.

For genes quantified in at least two cohorts, log2 fold changes were pooled by DerSimonian–Laird random-effects meta-analysis [27], with *I*² and τ² reported. Because k = 2–3 here, the normal-quantile random-effects interval under-covers, since it treats τ² as if it were known; all pooled intervals were therefore also computed with the Hartung–Knapp–Sidik–Jonkman correction [28] (standard error rescaled by √(q/(k−1)), *t*_{k−1} quantile). Both are reported, and equivalence counts are given under both.

### Gene set enrichment

Preranked GSEA [29] was run separately in each cohort. Genes were ranked by Welch *t* (case versus control), winsorised at ±50 as a safeguard which in the event never bound after filtering (maximum |*t*| 21.8 across the three cohorts). Genes below the cohort's own lower-quartile mean expression were removed before ranking. This is a pre-processing choice and a researcher degree of freedom, so its effect is quantified rather than left implicit: the same sets were re-run with the filter at the cohort median and with no filter (Additional file 2). Enrichment scores were the weighted Kolmogorov–Smirnov statistic with exponent 1. The null distribution was 10,000 random gene sets of identical size drawn from the same ranked list under `numpy.random.default_rng(seed = 20261029)`; gene-set permutation is the only option available for a preranked list in the Broad implementation, and the reason the ranking is not permuted is that a preranked list carries no sample-level information — all three contrasts here are balanced. This null asks whether a set is more extreme than a random set of the same size on this ranking; it does not preserve gene–gene correlation and is anti-conservative as a test of phenotype association. Nominal *p* is one-sided in the observed direction and is bounded below by 1/(10,000 + 1) = 1 × 10^−4^, a bound attained by five of the six sets in C1. Normalised enrichment scores divide the observed score by the mean of the null on the matching sign side; that divisor is sign-specific and set-size-specific, so NES values are comparable only within the same cohort, sign and set size, and cross-cohort NES comparisons in this paper are descriptive. *Q*-values are Benjamini–Hochberg across the six sets tested in each cohort, not across MSigDB.

Six gene sets were used. Two concern Wnt and are not interchangeable. Gene-level testing and equivalence assessment used a 37-gene curated Wnt/β-catenin panel spanning core signal-transduction machinery, ligands and canonical transcriptional targets. Set-level testing used HALLMARK_WNT_BETA_CATENIN as retrieved (38 symbols) [30]. The two overlap at 14 genes only, so they ask different questions, and they are reported as separate lines of evidence. The remaining five were curated marker panels: PYROPTOSIS_CORE (22 genes), INJURY_RESPONSE_POSCTRL (23), MICROGLIA_ACTIVATION (20), ASTROCYTE_REACTIVITY (19) and NFKB_CYTOKINE_SIGNALLING (22). Membership is given in Additional file 1.

Rank clustering was additionally tested without permutation: for a subset of size *k*, the number of members within the top 6% of the ranked list was compared with the hypergeometric expectation. This is independent of the permutation scheme but not of the ranking or the filtering that produced it.

### Sensitivity analyses

Six analyses probed the dependence of the conclusions on analysis choices. (i) Each of the 13 leading-edge pyroptosis genes was removed in turn and GSEA re-run; all 13 were also removed at once. (ii) The 22-gene pyroptosis panel was split into disjoint execution (9 genes) and priming (13 genes defined, 10 measurable) subsets. (iii) GSEA was repeated under four pre-processing variants. (iv) The Wnt panel was re-pooled using C1 and C2 only, excluding the cohort whose positive control failed. (v) Library normalisation was varied: raw counts, DESeq2-style median-of-ratios size factors, and upper-quartile factors, each with and without independent filtering before Benjamini–Hochberg correction (C1 and C2 only; C3 is already normalised abundance). (vi) Fixed-effect pooling was computed alongside random-effects pooling.

### Power

Because the paper's method claim concerns resolution, the minimum detectable effect was computed from the data rather than assumed. For each cohort the median within-group standard deviation of log2 expression was taken across all genes, and the effect detectable at 80% power with α = 0.05 was computed for the median gene, both at nominal α and at the *p*-value threshold a gene must beat to survive Benjamini–Hochberg correction in that cohort.

### Single-cell reference

GSE134003 was tabulated from the deposited ClusterMetadata (66 subtypes; the 735 microglia are the sum of deposited subtypes Micro1–Micro3, the 40 astrocytes the single deposited "Astrocytes" label). Cell numbers and contributing animals were counted per type and condition. No differential expression was performed at single-cell level.

### Software and reproducibility

Analyses were run in Python 3.13.14 with NumPy 2.5.4 and SciPy 1.18.1; plots and the additional files used Matplotlib 3.11.1 and openpyxl 3.1.5. Exact versions are pinned in the repository's requirements.txt. HALLMARK_WNT_BETA_CATENIN is reproduced as retrieved from MSigDB (38 symbols), with its membership frozen in the repository's gene_panels.json; the curated panels are likewise versioned there. Scripts, gene panels, per-cohort result tables, GSEA outputs and a SHA-256 manifest of every file are available at [31] (tag v1.1.4).

---

## Results

### The assay detects injury, but the design lacks the resolution to declare it

The three cohorts contributed 20 libraries in total (C1: 4 sham, 4 CCI; C2: 3 control, 3 SNI; C3: 3 sham, 3 SNI). Between 17,412 and 19,864 genes were testable per cohort after gene-level collapse (Table 1).

*Atf3*, the canonical injury-induced transcription factor, behaved as required of a positive control in the discovery cohort: +4.03 log2 (16.3-fold; 95% CI 2.07 to 5.98; *p* = 0.0044). It rose in the same direction in both other cohorts (+1.33, *p* = 0.0053; +1.03, *p* = 0.119), and the pooled random-effects estimate was +1.99 log2 (95% CI 0.64 to 3.33; *p* = 0.0039) with substantial heterogeneity (*I*² = 86%). Injury-induced transcription is therefore present and detectable in this material.

Correction nonetheless removed it. *Atf3* reached *q* = 0.466 in C1, 0.700 in C2 and 0.761 in C3. Across the transcriptome, 1,371 of 19,009 genes reached nominal *p* < 0.05 in C1, 1,210 of 17,412 in C2 and 734 of 19,864 in C3; the numbers surviving *q* < 0.05 were 0, 0 and 1 (Table 1). The single survivor (LOC688442 in C3) has a fold change of 1.09 with a *t* statistic of 1.27 × 10^16^ — a near-zero-variance transcript carrying no biological information.

Two readings of this are possible and must be separated. One is that the pipeline is at fault: it applies log2(x + 1) to raw counts with no library-size normalisation and no independent filtering, both of which are standard and both of which increase power. It is not. Re-running C1 and C2 from raw counts under median-of-ratios and upper-quartile normalisation, with and without independent filtering, gave between 0 and 6 survivors out of ~19,000 genes in C1, 0 in C2 and 0–1 in C3, and *Atf3* remained non-significant under every regime (16.3- to 17.2-fold; *q* = 0.34–0.47) (Additional file 2). The other reading is that the design is at fault, and the power calculation supports it: for the median gene, the smallest effect detectable at 80% power is 1.38-fold at nominal *p* < 0.05 in C1 but 5.6-fold after correction, and 1.79-fold rising to 428-fold in C2 (Table 1). The limiting factor is sensitivity, not a failure of false-discovery control [11]. It follows that any claim about a pathway in this literature resting on nominal *p* is a claim made below the detection threshold of the design.

### The Wnt/β-catenin transcriptional programme is not induced

Of the 37 genes in the curated Wnt panel, 36 were quantified in at least two cohorts (30 in all three) and could be pooled. The largest pooled effect in either direction was +0.27 log2 for *Ccnd2* (1.20-fold; 95% CI 0.00 to 0.53; pooled from C2 and C3 only, as it is absent from C1). *Jun* (+0.26), *Birc5* (+0.24) and *Cd44* (+0.19) followed. No pooled point estimate exceeded 0.3 log2 in either direction.

How far this supports equivalence depends on the interval, and at k = 2–3 the choice matters. Thirty-one of the 36 pooled intervals lay wholly within ±0.5 log2 under the normal-quantile DerSimonian–Laird interval, but only 18 of 36 did so under the Hartung–Knapp–Sidik–Jonkman correction (Table 1, Fig. 1C). We treat 31/36 as an upper bound and report both: for 18 genes an induction of 1.41-fold or more is excluded, and the remaining 18 are indeterminate rather than null. Five genes — *Apc*, *Ccnd2*, *Birc5*, *Wnt2* and *Wnt3a* — fail the criterion under either interval.

The core components were flat. *Ctnnb1* (β-catenin) itself was +0.08 log2 (pooled; *p* = 0.16), *Axin2*, the most responsive canonical target, −0.02, *Lef1* −0.01, *Tcf4* −0.12 and *Gsk3b* −0.05 (Fig. 1C). The largest decrease was *Wnt2* at −0.17, whose interval (−0.75 to 0.41) is too wide to support equivalence in either direction.

At the set level the result was the same. For HALLMARK_WNT_BETA_CATENIN (29 to 33 of 38 symbols passing the expression filter) GSEA gave NES −1.29 (*p* = 0.100, *q* = 0.100) in C1, +1.02 (*p* = 0.304, *q* = 0.370) in C2 and +0.48 (*p* = 0.931, *q* = 0.931) in C3 (Fig. 2A). None reached significance and the sign differed between cohorts; because the null distributions themselves differ between cohorts (mean −0.18 in C1, +0.12 in C2, +0.30 in C3), these three values are descriptive rather than strictly comparable. In C1 the negative score was produced by mild depletion of *Axin2*, *Tcf7*, *Tcf7l1*, *Lef1*, *Pax2*, *Wif1* and *Cxxc4* among 13 leading-edge genes: the direction is depletion, not activation (Fig. 2B).

Because the Wnt null is the paper's main claim, its dependence on pre-processing was examined directly. Across four variants in C1 the NES ranged from −1.29 to −1.01 with *p* between 0.10 and 0.25; the variant most likely to move a *t*-based ranking, removing the expression filter altogether, moved NES only from −1.29 to −1.28 (*p* = 0.105) (Additional file 2). The null is not a product of the pre-processing choices that were varied.

### Excluding the cohort that failed its positive control does not change the Wnt result

C3 did not behave as a replication cohort. Five of the six gene sets were negative at *q* < 0.05 (INJURY_RESPONSE_POSCTRL NES −1.30, *q* = 0.018; MICROGLIA_ACTIVATION −2.42, *q* = 0.001; NFKB_CYTOKINE −1.63, *q* = 0.010; ASTROCYTE_REACTIVITY −1.39, *q* = 0.018; PYROPTOSIS_CORE −1.11, *q* = 0.040), its *Atf3* response was 2.04-fold and did not reach significance (*p* = 0.119), and its glial markers were at or below sham (*Aif1* 0.80-fold, *C1qa* 0.90-fold, *Gfap* 1.01-fold, *Nlrp3* 0.97-fold). The panels share members, and the gene-set-permutation null in C3 is itself displaced positive (mean +0.30, 89% of draws positive, so that a set with ES = 0 already lies at the 11th percentile). We therefore read this as a global property of the C3 ranking rather than as five independent biological findings, and do not interpret the direction of individual sets in C3. Two independent metadata defects support treating it as non-confirmatory on grounds that do not depend on its results: the GEO record states "Assembly: GRCh38.p13", a human assembly, for a rat dataset, and the post-injury interval is not recorded for C3 whereas it is day 10 for C1.

Re-pooling the Wnt panel over C1 and C2 left the conclusion unchanged in point estimate but not in precision: 31 genes were poolable, the largest induction was *Tcf3* at +0.30 log2 (1.23-fold), and the only gene whose point estimate exceeded 0.5 was *Wnt2*, which moved *down* (−0.50; 95% CI −0.86 to −0.15, *p* = 0.0056). Under the normal-quantile interval 20 of 31 genes met the equivalence criterion, but under the Hartung–Knapp correction only 2 of 31 did. The drop is a loss of precision from having one fewer cohort and a *k* of 2 throughout, not a larger effect: no point estimate crossed the bound. C3 is retained in the analysis because the design was fixed in advance and is reported in full; it is not counted as independent support.

### Pyroptosis transcripts are coordinately up-regulated in one cohort only

The pyroptosis panel behaved differently from Wnt. In C1 it was strongly enriched (NES +2.68, *p* ≤ 1 × 10^−4^, the permutation floor: no one of 10,000 random gene sets of size 19 reached the observed score, so the *p*-value is right-censored and must not be read as a measured magnitude). Thirteen of 19 quantified members were in the leading edge, all within the top 844 of 14,257 ranked genes; rank clustering was far beyond chance (13 of 19 in the top 5.9%, expected 1.14; hypergeometric *p* = 2.3 × 10^−12^, and 3.7 × 10^−11^ after Bonferroni correction over the 19 cut-points examined). Removing each leading-edge gene in turn reduced NES from 2.68 to no lower than 2.60, and all 13 remained at the permutation floor. Because the *p*-value cannot fall below that floor, this shows only that no single gene is necessary for an enrichment score of this magnitude. The enrichment is nevertheless carried by the leading edge as a block: removing all 13 at once leaves six genes with NES −0.48 (*p* = 0.69) (Additional file 2).

What the gene-level estimates show is narrower than any statement about pathway stage. The leading edge consisted of *Pycard*, *Nlrc4*, *Tlr4*, *Il6*, *Il1b*, *Gsdmd*, *Casp1*, *Txnip*, *Casp4*, *Il18*, *Rela*, *Aim2* and *Casp8*; *Nlrp3* was absent from it. The inducible components moved most (*Il6* +1.98, *p* = 0.021; *Il1b* +1.73, *p* = 0.063; *Pycard* +0.58, *p* = 0.003), while the sensor and one effector did not move at all (NLRP3 +0.01, *p* = 0.97; NEK7 −0.01, *p* = 0.95) and the effectors showed small changes (CASP1 +0.28, *p* = 0.11; GSDMD +0.22, *p* = 0.06; GSDME +0.06, *p* = 0.69). We attempted to formalise the distinction by splitting the panel into execution and priming subsets; it did not separate them (NES +2.34 for 9 execution genes versus +2.12 for 10 priming genes), and 6 of the 9 execution genes also fell within the top 6% of the ranking (expected 0.54; hypergeometric *p* = 3.3 × 10^−6^; 5 of 9, *p* = 7.9 × 10^−5^, without the expression filter). The three that did not were *Gsdme* (17th percentile), *Nlrp3* (25th) and *Nek7* (28th). The set-level statistic therefore carries no information about which stage of the pathway is engaged, and we draw no conclusion about that.

Two further points keep this result within its evidence. First, by the equivalence criterion used elsewhere in this paper only GSDMD ([−0.02, +0.45]) is equivalent to no change; NLRP3 ([−0.86, +0.89]), NEK7 ([−0.50, +0.48]) and CASP1 ([−0.09, +0.65]) have intervals too wide to support equivalence, and CASP1's does not exclude an induction of 1.41-fold. Second, inflammasome activation is rate-limited by post-translational events — assembly, caspase-1 autoprocessing and gasdermin D cleavage — that follow a second signal and are largely independent of transcript abundance [32]. Transcript levels index the priming step, not execution.

The pattern did not replicate. In C2 the panel gave NES +0.76 (*p* = 0.514) and in C3 −1.11 (*p* = 0.033, opposite in sign, within a globally displaced ranking). The coordinated up-regulation is a property of the CCI cohort, not a general feature of these models.

### The single-cell reference cannot localise the signal to astrocytes

GSE134003 is a mouse atlas (18,856 cells, lumbar L1–6, Drop-seq, ~1,200 transcripts per cell, 14 days after injury) [24–25]. It contains 735 microglia (385 SI, 350 SNI; 10 and 9 animals) and 40 astrocytes (27 SI, 13 SNI; 9 and 5 animals); one SNI sample was excluded by the original authors for dorsal-root-ganglion contamination, which is why the injured group has 9 animals. The astrocyte compartment averages 2.6 cells per injured animal. At that sampling depth — low-depth Drop-seq over whole-segment dissociations — the compartment is too shallow to support cell-type attribution, and we make none. This is a statement about the depth of this atlas. Deeper or astrocyte-enriched preparations do resolve spinal astrocytes after nerve injury [33], so the claim that Wnt and pyroptosis converge in astrocytes is untested here rather than untestable in principle; it is also a mouse atlas, whereas the bulk cohorts are rat. Glial heterogeneity and myelin pathology after spinal injury are well documented in single-cell work [34–35], but none of that material resolves the rat astrocyte compartment at the depth this question requires.

---

## Discussion

### Principal findings

In three independent rat spinal cord cohorts, a Wnt/β-catenin transcriptional programme was not induced after nerve injury: the largest pooled effect across 36 genes was 1.20-fold, set-level enrichment was non-significant in all three cohorts and under all four pre-processing variants, and no point estimate exceeded 0.3 log2. At the same time a gene that is unambiguously induced, *Atf3* at 16-fold, failed correction, and no gene in two of three cohorts survived false-discovery control under any of six pre-processing regimes. The second finding explains why the first is worth stating: at this sample size a pathway can be reported as activated on evidence that the same data cannot support.

### What a transcriptional null does and does not exclude

β-catenin activation is predominantly post-translational: destruction-complex disassembly, phosphorylation, nuclear translocation and co-factor recruitment all occur without changes in transcript abundance [2]. A transcriptomic null therefore cannot exclude Wnt signalling activity, and we do not claim that it does. What it excludes is the specific and commonly asserted proposition that nerve injury induces a Wnt/β-catenin transcriptional programme in spinal cord — the coordinated induction of canonical targets such as *Axin2*, *Ccnd1*, *Lef1* and *Tcf7*. That proposition is testable in bulk tissue and it fails here: *Axin2* moved −0.02 log2 and the β-catenin transcript +0.08.

The anatomical limitation is more serious than the post-translational one and we did not state it in the first instance. All three cohorts profiled whole spinal cord or a whole lumbar segment, whereas the studies that report Wnt/β-catenin engagement after nerve injury localise it to the dorsal horn [3–4]. A programme confined to the superficial dorsal horn would be diluted in whole-cord bulk RNA-seq by the volume and RNA mass of the remaining tissue, and its apparent fold change would fall in exactly the range we observe (pooled effects ≤ 0.27 log2). Our null therefore excludes a spinal-cord-wide induction of the magnitude specified, not a dorsal-horn-restricted one. Resolving this requires dorsal-horn-microdissected or spatially resolved material, which we do not have. Species is a second restriction: the motivating studies were conducted mainly in mouse [3–5,8–9] and our data are rat.

The distinction is not pedantic. Four studies motivated this analysis and all report pathway-dependence or ligand and target induction in tissue-level assays [3–6]. None of them reports the coordinated induction of canonical TCF/LEF target genes that a transcriptional programme would produce; their readouts are predominantly protein-level or pharmacological. Our result is therefore consistent with those reports at the level they were made and inconsistent only with the stronger reading that a transcriptional programme is engaged.

### Why we do not interpret the pyroptosis signal as inflammasome activation

The pyroptosis result initially appeared to contradict the Wnt result — an enriched inflammatory panel alongside a flat Wnt panel — and it is worth being explicit about how it was handled. The set-level enrichment is real and robust: it survives leave-one-out at every leading-edge gene, it holds under all four pre-processing variants, and it is confirmed by a permutation-free hypergeometric test. It is also, on inspection, not informative about pathway stage: the execution/priming split did not separate the two (NES +2.34 versus +2.12), and 6 of 9 execution genes sit in the top 6% of the ranking alongside the priming genes.

The biology points the same way. Inflammasome activation is rate-limited by post-translational events that follow a second signal and are largely independent of the transcript levels of the proteins involved [32]; transcript abundance is informative about priming, which is NF-κB-dependent, and not about execution. Our data therefore neither support nor refute inflammasome execution. They show that in the CCI cohort the NF-κB-dependent transcriptional footprint of priming is present, that the transcripts encoding the sensor and one pore-forming effector are not induced, and that — by our own equivalence criterion, applied honestly — the intervals for NLRP3, NEK7 and CASP1 are too wide to call those genes unchanged. We explicitly do not conclude that pyroptosis does not occur. The protein-level evidence in this literature, including NLRP3, caspase-1 and GSDMD-N measured by immunoblot and immunofluorescence in rat and mouse nerve-injury models [7–10], lies in the layer our design cannot see and is not contradicted by anything we report. Indeed the source study of C1 reports that levo-tetrahydropalmatine lowers NLRP3 and caspase-1 p20 at protein level [17]; our finding that the corresponding transcripts are not induced is complementary, not contradictory.

### Why the literature contains the opposite conclusion

The two problems documented above are mutually reinforcing. A study that tests 19,000 genes at n = 3–4 and reports nominal *p* < 0.05 will find roughly 1,000 "differentially expressed" genes by chance alone; among 20–40 genes of a candidate pathway several will fall below 0.05 and the pathway will appear activated. Our C1 data illustrate this exactly: 1,371 genes at nominal *p* < 0.05, none after correction, and a positive control at *q* = 0.47.

The consequence is not that individual studies are dishonest, but that the field's inferential standard is below its design's resolution. Two remedies are cheap and neither requires new experiments: report effect sizes with confidence intervals against a stated equivalence bound, and report the minimum detectable effect alongside non-significant results. The first converts "not significant" into "excludes effects above 1.41-fold"; the second makes the claim's resolution explicit. Published estimates of the replicate number required for RNA-seq [11] and of the consequences of underpowered neuroscience [12] have been available for a decade and have largely not been adopted here.

### Limitations

First, all data are transcriptomic. No protein, phosphorylation or localisation data were generated, and the post-translational events that define both Wnt and inflammasome activity are invisible to this design. The conclusions are restricted to transcript abundance.

Second, the cohorts are small (n = 3–4 per group), and this is the object of study as much as a limitation. Power is sufficient to detect a 16-fold change and insufficient, after correction, to detect a 1.4-fold one; genes whose intervals exceed the equivalence bound are indeterminate rather than null, and under the Hartung–Knapp correction that applies to 18 of 36 Wnt genes.

Third, the anatomical restriction is the most serious one (see above): a dorsal-horn-restricted programme would be diluted into the range we observe and is not excluded. None of the three cohorts microdissected dorsal horn.

Fourth, C3 is not confirmatory, for reasons that include two metadata defects independent of its results. It is retained because the design was fixed in advance, and it is not counted as support. The post-injury interval is recorded for C1 (day 10) but not for C2 or C3, so a between-cohort difference in interval cannot be excluded.

Fifth, most panels are curated rather than drawn from MSigDB, and the *q*-values apply to the six sets tested, not to a genome-wide collection; a full Hallmark scan would give larger *q*-values and was not performed. The gene-set-permutation null is anti-conservative as a test of phenotype association, and its location differs between cohorts, so the NES comparisons across cohorts in this paper are descriptive. The curated and Hallmark Wnt sets overlap at 14 genes only and are reported as separate lines of evidence, not as two views of one object.

Sixth, the equivalence bound of ±0.5 log2 is a convention, not a measurement. A reader who prefers ±0.3 would count fewer genes as equivalent under either interval; the direction would not change, because no pooled point estimate exceeds 0.3.

Seventh, we analysed one contrast per dataset; other contrasts, time points and spinal levels exist in the same series and were outside the pre-specified scope. This analysis was not formally pre-registered.

### Conclusions

A Wnt/β-catenin transcriptional programme is not detectable in bulk rat spinal cord after peripheral nerve injury, across three independent cohorts, four pre-processing variants and six normalisation regimes; a dorsal-horn-restricted programme is not excluded. The pyroptosis literature's central observation, coordinated induction of inflammasome-associated transcripts, is reproducible in the CCI cohort only, and it indexes NF-κB-dependent priming rather than inflammasome execution — which transcript abundance cannot address in either direction. More consequentially, the design on which these claims rest cannot support transcriptome-wide inference at the effect sizes being claimed: at n = 3–4 the smallest effect detectable after correction is 5.6-fold in the largest cohort. Progress requires larger samples or explicit equivalence and power reporting, and the second is available immediately.

---

## Declarations

**Ethics approval and consent to participate.** Not applicable. This study used only publicly available, de-identified rodent data deposited in the Gene Expression Omnibus; no animals were used and no human participants were involved. Responsibility for the original animal experiments, including ethical approval, rests with the investigators who generated and deposited each dataset [17,22,24–25].

**Consent for publication.** Not applicable.

**Availability of data and materials.** All datasets analysed are publicly available from the Gene Expression Omnibus [16]: GSE217932 [18], GSE256472 [20], GSE205906 [23] and GSE134003 [25]. Analysis code, gene panels, per-cohort result tables, GSEA outputs and a SHA-256 manifest of every file are available at [31] under tag v1.1.4. No third party generated or analysed data for this study.

**Competing interests.** The author declares no competing interests.

**Funding.** This work received no specific grant from any funding agency in the public, commercial or not-for-profit sectors.

**Authors' contributions.** YY conceived the study, designed the analysis, performed all analyses, interpreted the results, wrote the manuscript and read and approved the final manuscript.

**Acknowledgements.** Not applicable.

**Authors' information.** YY is a physician in the Department of Anesthesiology, The Second Affiliated Hospital of Fujian University of Traditional Chinese Medicine. ORCID 0009-0004-9698-6552.

**Use of AI tools.** No generative AI tool was used to produce the analyses, the figures or the text. Language-editing assistance was not used.

---

## List of abbreviations

BH, Benjamini–Hochberg; CCI, chronic constriction injury; CI, confidence interval; DL, DerSimonian–Laird; ES, enrichment score; FDR, false discovery rate; FPKM, fragments per kilobase of transcript per million mapped reads; GEO, Gene Expression Omnibus; GSEA, gene set enrichment analysis; HKSJ, Hartung–Knapp–Sidik–Jonkman; l-THP, levo-tetrahydropalmatine; MDE, minimum detectable effect; MSigDB, Molecular Signatures Database; NES, normalised enrichment score; PLM, paeoniflorin–liquiritin; RPKM, reads per kilobase of transcript per million mapped reads; SNI, spared nerve injury; TCF/LEF, T-cell factor / lymphoid enhancer-binding factor.

---

## Figure legends

**Fig. 1** Three independent rat spinal cord cohorts, one pipeline, and no induction of the Wnt/β-catenin transcriptional programme. (**A**) Pre-specified design: three bulk RNA-seq cohorts analysed with one pipeline; intervention arms (CCI + levo-tetrahydropalmatine, n = 4; paeoniflorin–liquiritin, n = 3) excluded before any comparison. (**B**) Positive control. *Atf3* rose 16.3-fold in the discovery cohort (nominal *p* = 0.0044), establishing that injury-induced transcription is detectable at this sample size; the same gene nevertheless failed Benjamini–Hochberg correction (*q* = 0.466), as it did in the two further cohorts (*q* = 0.700 and 0.761). (**C**) All 36 genes of the 37-gene curated Wnt/β-catenin panel quantified in at least two cohorts. Points show per-cohort log2 fold change (C1 red, C2 blue, C3 green); diamonds and bars show the DerSimonian–Laird pooled estimate and its 95% normal-quantile confidence interval. The grey band marks the pre-specified equivalence bound of ±0.5 log2 (1.41-fold); 31 of 36 pooled intervals lie wholly within it under this interval and 18 of 36 under the Hartung–Knapp–Sidik–Jonkman correction, the remainder being indeterminate. The largest pooled effect is 0.27 log2. (**D**) Genes reaching nominal *p* < 0.05 versus *q* < 0.05 in each cohort. The single gene passing correction (C3) is a near-zero-variance transcript (LOC688442, 1.09-fold). *Note*: re-pooling over C1 and C2 alone leaves the largest induction at 1.23-fold (*Tcf3*) and no point estimate beyond the bound, but reduces the equivalence count to 20 of 31 (normal-quantile) or 2 of 31 (Hartung–Knapp).

**Fig. 2** Preranked gene set enrichment in three cohorts. (**A**) Normalised enrichment scores for the six gene sets, per cohort. Ranking metric: Welch *t* (case versus control), winsorised at ±50, genes below the cohort's lower-quartile mean expression removed; 10,000 gene-set permutations under seed 20261029; *q* values are Benjamini–Hochberg across the six sets tested in each cohort, and nominal *p* is bounded below at 1 × 10^−4^. NES values are comparable within a cohort, sign and set size only. HALLMARK_WNT_BETA_CATENIN is non-significant and direction-inconsistent across cohorts. (**B**) Running enrichment score for the Wnt panel in C1 (NES −1.29, *p* = 0.100). The negative score is generated by mild depletion of *Axin2*, *Tcf7*, *Tcf7l1*, *Lef1*, *Pax2*, *Wif1* and *Cxxc4*, whose positions are marked; the direction is depletion, not activation.

---

## Additional files

**Additional file 1.** File format: XLSX. Title: Per-cohort gene-level statistics, pooled estimates with both interval methods, GSEA outputs and gene-set membership. Description: Table S1, per-cohort gene-level statistics (log2 fold change, Welch *t*, df, nominal *p*, Benjamini–Hochberg *q*, 95% CI, expression summaries) for every gene tested in each cohort. Table S2, DerSimonian–Laird random-effects and fixed-effect pooled estimates for the curated Wnt panel, the pyroptosis panel and an 18-gene injury-marker panel including *Atf3*, with normal-quantile and Hartung–Knapp–Sidik–Jonkman intervals. Table S3, GSEA output for all cohorts and sets with the four pre-processing variants and the execution/priming split in C1. Table S4, membership of all six gene sets. Table S5, sample-to-GSM mapping including excluded libraries.

**Additional file 2.** File format: PDF. Title: Sensitivity of the enrichment scores to pre-processing, leave-one-out analysis, and library-normalisation regimes. Description: Figure S1, dependence of the Wnt and pyroptosis enrichment scores on the four pre-processing settings in C1, and on removal of the entire leading edge. Figure S2, leave-one-out analysis of the 13 leading-edge pyroptosis genes. Figure S3, number of genes surviving Benjamini–Hochberg correction and the *Atf3* effect under raw, median-of-ratios and upper-quartile normalisation, with and without independent filtering.

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

16. Gene Expression Omnibus. National Center for Biotechnology Information. https://www.ncbi.nlm.nih.gov/geo/. Accessed 30 Oct 2026.

17. Wu D, Wang P, Zhao C, Su J, Zhang J, Ma W, et al. Levo-tetrahydropalmatine ameliorates neuropathic pain by inhibiting the activation of the Clec7a-MAPK/NF-κB-NLRP3 inflammasome axis. Phytomedicine. 2023;121:155075. doi:10.1016/j.phymed.2023.155075

18. Wu D, Wang P, Zhao C, Su J, Zhang J, Ma W, et al. [dataset]. Levo-tetrahydropalmatine Ameliorates Neuropathic Pain via Suppressing the Activation of the Clec7a-MAPK/NF-κB-NLRP3 Inflammasome Axis. Gene Expression Omnibus; 2023. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217932. Accessed 30 Oct 2026.

19. Bennett GJ, Xie YK. A peripheral mononeuropathy in rat that produces disorders of pain sensation like those seen in man. Pain. 1988;33(1):87-107. doi:10.1016/0304-3959(88)90209-6

20. [dataset]. Spinal cord RNA-seq of control and spared nerve injury rats. Gene Expression Omnibus; 2024. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE256472. Accessed 30 Oct 2026.

21. Decosterd I, Woolf CJ. Spared nerve injury: an animal model of persistent peripheral neuropathic pain. Pain. 2000;87(2):149-158. doi:10.1016/S0304-3959(00)00276-1

22. Chen YY, Feng LM, Xu DQ, Yue SJ, Fu RJ, Zhang MM, et al. Combination of paeoniflorin and liquiritin alleviates neuropathic pain by lipid metabolism and calcium signaling coordination. Front Pharmacol. 2022;13:944386. doi:10.3389/fphar.2022.944386

23. Chen YY, Feng LM, Xu DQ, Yue SJ, Fu RJ, Zhang MM, et al. [dataset]. Combination of paeoniflorin and liquiritin alleviates neuropathic pain by lipid metabolism and calcium signaling coordination. Gene Expression Omnibus; 2022. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205906. Accessed 30 Oct 2026.

24. Niehaus JK, Taylor-Blake B, Loo L, Simon JM, Zylka MJ. Spinal macrophages resolve nociceptive hypersensitivity after peripheral injury. Neuron. 2021;109(8):1274-1282.e6. doi:10.1016/j.neuron.2021.02.018

25. Niehaus JK, Taylor-Blake B, Loo L, Simon JM, Zylka MJ. [dataset]. Spinal macrophages resolve nociceptive hypersensitivity after peripheral injury. Gene Expression Omnibus; 2021. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134003. Accessed 30 Oct 2026.

26. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B Stat Methodol. 1995;57(1):289-300. doi:10.1111/j.2517-6161.1995.tb02031.x

27. DerSimonian R, Laird N. Meta-analysis in clinical trials. Control Clin Trials. 1986;7(3):177-88. doi:10.1016/0197-2456(86)90046-2

28. Röver C, Knapp G, Friede T. Hartung-Knapp-Sidik-Jonkman approach and its modification for random-effects meta-analysis with few studies. BMC Med Res Methodol. 2015;15:99. doi:10.1186/s12874-015-0091-1

29. Subramanian A, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, Gillette MA, et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proc Natl Acad Sci U S A. 2005;102(43):15545-50. doi:10.1073/pnas.0506580102

30. Liberzon A, Birger C, Thorvaldsdóttir H, Ghandi M, Mesirov JP, Tamayo P. The Molecular Signatures Database (MSigDB) hallmark gene set collection. Cell Syst. 2015;1(6):417-425. doi:10.1016/j.cels.2015.12.004

31. Yang Y. astrocytic-wnt-pyroptosis-np: three-cohort reanalysis of Wnt/β-catenin and pyroptosis transcription in neuropathic pain spinal cord (v1.1.4). GitHub; 2026. https://github.com/yyx-4113/astrocytic-wnt-pyroptosis-np. Accessed 30 Oct 2026.

32. Swanson KV, Deng M, Ting JP. The NLRP3 inflammasome: molecular activation and regulation to therapeutics. Nat Rev Immunol. 2019;19(8):477-489. doi:10.1038/s41577-019-0165-0

33. Xu J, Yan Z, Bang S, Velmeshev D, Dodge K, Lopez MA, et al. GPR37L1 identifies spinal cord astrocytes and protects neuropathic pain after nerve injury. Neuron. 2025;113(8):1206-1222.e6. doi:10.1016/j.neuron.2025.01.012

34. Li D, Yang K, Li J, Xu X, Gong L, Yue S, et al. Single-cell sequencing reveals glial cell involvement in development of neuropathic pain via myelin sheath lesion formation in the spinal cord. J Neuroinflammation. 2024;21(1):213. doi:10.1186/s12974-024-03207-3

35. Milich LM, Choi JS, Ryan C, Cerqueira SR, Benavides S, Yahn SL, et al. Single-cell analysis of the cellular heterogeneity and interactions in the injured mouse spinal cord. J Exp Med. 2021;218(8):e20210040. doi:10.1084/jem.20210040

---

## Table 1

**Table 1** Cohort characteristics, transcriptome-wide results, and the resolving power of each design

| | C1 GSE217932 | C2 GSE256472 | C3 GSE205906 |
|---|---|---|---|
| Model | CCI | SNI | SNI |
| Groups (control vs injury) | 4 vs 4 | 3 vs 3 | 3 vs 3 |
| Quantification | raw counts | raw counts (deposit also has FPKM; counts used) | processed RPKM |
| Excluded arms (GSM) | CCI + l-THP, n = 4 (6730109–6730112) | none | PLM, n = 3 (6235298–6235300) |
| Spinal segments | L4–L6 | not recorded in GEO | L3–L6 |
| Interval after injury | postoperative day 10 | not recorded in GEO | not recorded in GEO |
| Tissue dissected | whole segment (no dorsal horn / DRG) | same | same |
| Genes tested | 19,009 | 17,412 | 19,864 |
| Nominal *p* < 0.05 | 1,371 | 1,210 | 734 |
| *q* < 0.05 | 0 | 0 | 1^a^ |
| *q* < 0.05 over normalisation regimes^b^ | 0–6 | 0–0 | 0–1 |
| Median log2FC (depth diagnostic) | −0.073 | +0.029 | +0.080 |
| Intervals within ±0.5 log2 | 9,149 (48.1%) | 3,373 (19.4%) | 8,902 (44.8%) |
| *Atf3* log2FC (95% CI) | +4.03 (2.07 to 5.98) | +1.33 (0.75 to 1.91) | +1.03 (−0.43 to 2.48) |
| *Atf3* *p* / *q* | 0.0044 / 0.466 | 0.0053 / 0.700 | 0.119 / 0.761 |
| *Atf3* *q* over regimes^b^ | 0.335–0.466 | 0.550–0.750 | 0.761 |
| MDE, nominal *p* < 0.05^c^ | 1.38-fold | 1.79-fold | 1.38-fold |
| MDE, surviving BH *q* < 0.05^c^ | 5.6-fold | 428-fold | 31.7-fold |
| Wnt NES (*q*) | −1.29 (0.100) | +1.02 (0.370) | +0.48 (0.931) |
| Pyroptosis NES (*q*) | +2.68 (≤1 × 10^−4^) | +0.76 (0.514) | −1.11 (0.040)^d^ |

^a^ Near-zero-variance transcript (LOC688442), 1.09-fold change; not biologically interpretable.
^b^ Raw counts, median-of-ratios and upper-quartile normalisation, each with and without independent filtering (C1, C2); C3 contributes two regimes because it is already normalised abundance.
^c^ Minimum detectable effect for the median gene at 80% power, α = 0.05, from the median within-group SD of log2 expression (C1 0.197, C2 0.276, C3 0.152).
^d^ Five of six sets are negative at *q* < 0.05 in C3; see text — this is a global displacement of the ranking, not a set-specific result.
