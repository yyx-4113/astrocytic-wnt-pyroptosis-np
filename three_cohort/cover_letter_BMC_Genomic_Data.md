# Cover letter — BMC Genomic Data

Dear Editors,

Please find enclosed our manuscript, "No induction of a Wnt/β-catenin transcriptional programme in bulk rat spinal cord after peripheral nerve injury: a three-cohort reanalysis with equivalence testing", which I submit for consideration as a Research article in *BMC Genomic Data*.

**Why this journal.** *BMC Genomic Data* states in its Aims and Scope that it "reports new analyses of genomic data" and, as a BMC Series journal, "does not make editorial decisions based on the perceived interest or potential impact of a study," considering manuscripts when they are scientifically valid with appropriate methodology and analysis. This manuscript is precisely a new analysis of public genomic (bulk RNA-seq) data: three independent rat transcriptomic cohorts (GSE217932, GSE256472, GSE205906) reanalysed through one pre-specified pipeline. The analysis — not a new wet-lab experiment — is the contribution, and that is squarely within this journal's scope. We note that a sibling submission was declined by a clinical neuroscience journal on the grounds that it reported "solely bioinformatics … not accompanied by independent validation"; this journal's stated scope accommodates exactly the kind of deposited-data reanalysis we present, with all inputs public and reusable for validation by others.

**What the paper does.** The proposition that peripheral nerve injury induces a Wnt/β-catenin transcriptional programme in spinal cord — and that pyroptosis proceeds alongside it in astrocytes — rests on bulk transcriptomic comparisons of three or four animals per group interpreted at nominal *p* < 0.05. We tested that proposition as a hypothesis capable of failing, using three independent rat cohorts from different laboratories run through one pipeline: effect sizes with 95% confidence intervals against a pre-specified ±0.5 log2 equivalence bound, random-effects pooling, and preranked gene set enrichment calibrated on *Atf3*.

**Main findings.** *Atf3* rose 16.3-fold and still reached only *q* = 0.466, and no gene survived correction in two of three cohorts; three library-normalisation regimes left this unchanged, and the smallest detectable effect after correction is 5.6-fold — so the null is a property of the design, not of the pipeline. Across 36 Wnt/β-catenin genes the largest pooled effect was +0.27 log2; 31 of 36 intervals lie within the equivalence bound under the normal-quantile interval but only 18 of 36 under the Hartung–Knapp correction appropriate at k = 2–3, and both counts are reported. Pyroptosis transcripts were coordinately up-regulated in the chronic constriction cohort only and did not replicate in either spared-nerve-injury cohort.

**On originality.** To our knowledge this is the first test designed so that the proposition could fail: the studies supporting it [3–6 of the manuscript] report activation from nominal p-values, state no equivalence bound, report no interval estimates, and analyse no more than one cohort. This manuscript does not duplicate them — it shows that at the sample sizes used, the same design cannot support any transcriptome-wide claim in either direction.

**On pooled analyses.** Although the manuscript reports a random-effects synthesis, it is not a pooled analysis of published research within the meaning of your criteria: no published effect size enters the analysis, and all inputs were recomputed from primary deposited data through a single pipeline implemented by the author. Cohort selection was by model and data availability, not by reported outcome.

**On experimental validation.** We wish to be explicit about a point that ended a previous submission. This is a re-analysis of deposited transcriptomic data by a single author without laboratory capacity; protein-level or in vivo validation cannot be added at revision stage. Every claim is therefore confined to transcript abundance, and that boundary is stated in the Methods, the Limitations and the Conclusions. The underlying datasets are public and reusable precisely so that others may validate or extend the analysis. If the handling editor considers external validation a precondition for review, we would rather learn that now than after peer review.

**Declarations.** This manuscript is not under consideration elsewhere, has not been published previously, and all material is original. The sole author drafted, analysed and approved the final manuscript. There are no competing interests. No new animal experiments or human participants were involved: all data are publicly available, de-identified rodent datasets from GEO. Use of AI-assisted tools is disclosed in the manuscript's Declarations, in the section "Use of AI tools". Analysis code, gene panels, per-cohort results, GSEA outputs and a SHA-256 manifest of every file are publicly available at https://github.com/yyx-4113/astrocytic-wnt-pyroptosis-np (tag v1.1.4), and every number in the manuscript is verified by a 113-check audit script against those artefacts.

Suggested reviewers can be provided on request; we have none with a conflict of interest to declare.

Yours sincerely,

Yongxin Yang
Department of Anesthesiology, The Second Affiliated Hospital of Fujian University of Traditional Chinese Medicine, Fuzhou, Fujian 350003, China
960856791@qq.com; ORCID 0009-0004-9698-6552
