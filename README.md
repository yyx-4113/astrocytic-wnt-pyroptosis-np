# astrocytic-wnt-pyroptosis-np

**No induction of a Wnt/β-catenin transcriptional programme in bulk rat spinal cord
after peripheral nerve injury: a three-cohort reanalysis with equivalence testing.**

Single-author study (Yongxin Yang, ORCID 0009-0004-9698-6552). All analyses use
publicly available, de-identified rodent data; no new experimental data were
generated.

This repository is the reproducibility package for the manuscript. The current
release is **v1.1.4**, which replaces the earlier single-cohort analysis in
`archive_v1.0.3/`. v1.1.4 adds the Hartung–Knapp correction, the library-
normalisation sensitivity analysis, the observed power calculation and the
two-cohort (C1+C2) re-pooling requested during review.

## What the study claims, and what it does not

1. **The Wnt/β-catenin transcriptional programme is not induced** after nerve
   injury in rat spinal cord. Across a 37-gene curated panel the largest pooled
   effect over three cohorts was +0.27 log2 (*Ccnd2*, 1.20-fold). 31 of 36 pooled
   confidence intervals lay wholly within a pre-specified ±0.5 log2 equivalence
   bound under the normal-quantile interval, but only **18 of 36** under the
   **Hartung–Knapp–Sidik–Jonkman** correction appropriate at k = 2–3; the other
   18 are indeterminate, not equivalent. Both counts are reported. Set-level
   testing with HALLMARK_WNT_BETA_CATENIN was non-significant in every cohort and
   under four pre-processing variants.
2. **At n = 3–4 the design lacks the resolution to support transcriptome-wide
   inference** — this is a *power* statement, not a false-discovery statement.
   *Atf3*, an unambiguously induced injury gene, rose 16.3-fold and still reached
   only *q* = 0.466; no gene survived correction in two of three cohorts, and the
   single survivor in the third was a near-zero-variance artefact. Three library
   normalisation regimes (raw, DESeq2 median-of-ratios, upper-quartile), with and
   without independent filtering, left this unchanged (0–6 of ~19,000 genes), so
   the null is a property of the design and not of the pipeline: after correction
   the smallest detectable effect is **5.6-fold** even in the largest cohort.
3. **Pyroptosis transcripts were coordinately up-regulated in the CCI cohort
   only** (NES +2.68, *p* ≤ 1 × 10⁻⁴, the permutation floor). This did not
   replicate in either SNI cohort. *Nlrp3* itself was unchanged (+0.01 log2,
   *p* = 0.97) and the effector transcripts changed only 1.16–1.21-fold, with
   intervals that do not exclude induction. Transcript level indexes priming, not
   execution; no protein, cleavage or pore-formation data exist here.

The study does **not** claim that Wnt signalling is inactive: β-catenin
activation is largely post-translational and invisible to bulk transcriptomics.
Nor does it claim the absence of a dorsal-horn-restricted programme: all three
cohorts sampled whole spinal segments, so a lamina-restricted programme would be
diluted into exactly the range observed. See the manuscript Discussion for the
exact scope of the negative result.

## Important honesty disclosures (read before reuse)

- **One cohort failed its own positive control.** GSE205906 (C3) showed a
  *negative* injury-response enrichment score and glial markers at or below sham
  levels. It is reported, downgraded, and excluded from the primary meta-analysis
  as a source of replication. Do not count it as independent support.
- **Two Wnt gene sets are used and they are not interchangeable.** Gene-level and
  equivalence testing use a 37-gene curated panel; set-level testing uses
  HALLMARK_WNT_BETA_CATENIN as retrieved (38 symbols). The two overlap at only
  14 genes. This is stated in the manuscript Methods and in the Limitations.
- **GSE186237 is a byte-identical duplicate of GSE217932** (identical MD5) and is
  excluded. It was used in earlier drafts as an independent replication; that
  claim is withdrawn.
- **The Mendelian-randomisation component of earlier drafts was removed.** No MR
  results, instruments or causal estimates are reported anywhere in this release.
- **The WGCNA module-preservation file from earlier drafts was forged**
  (self-comparison yielding `Z_summary = 1.0`) and has been deleted. Module
  assignment is reported descriptively, as exploratory only.
- **All data are bulk transcriptomic.** No protein, phosphorylation or
  localisation data were generated. Statements about pathway *activity* are
  outside what this design can support.
- **The headline equivalence count depends on the interval method.** 31/36 is an
  *upper bound* obtained with the normal-quantile interval, which under-covers at
  k = 2–3 because it treats τ² as known. Under Hartung–Knapp the count is 18/36.
  `meta_focus.csv` carries both (`equiv_normal`, `equiv_HKSJ`) so the reader can
  reproduce either. Do not quote 31/36 without the qualifier.
- **The pyroptosis GSEA *p*-value is right-censored.** With 10,000 gene-set
  permutations the smallest attainable *p* is 1 × 10⁻⁴; the reported value sits
  on that floor. The leave-one-out result ("all 13 re-runs significant") is
  therefore weak evidence — removing all 13 leading-edge genes at once gives
  NES −0.48, *p* = 0.69, and that is reported too.
- **Earlier drafts said "FDR control is unattainable at n = 3–4".** That is
  wrong and was withdrawn: Schurch et al. (2016) show the major tools *do* hold
  FDR near 5% at any replicate number; what is lost is power. The manuscript now
  makes the power claim, with the observed minimum detectable effect.

## Repository structure

```
three_cohort/                    # v1.1.4 analysis (this release)
  MANUSCRIPT_BMC.md              # manuscript source, BMC format
  MANUSCRIPT_BMC.docx            # submission-ready Word file
  Additional_File_1_tables.xlsx  # Tables S1-S5 (per-cohort stats, pooled + HKSJ,
                                 #   GSEA, panels, per-sample inventory)
  Additional_File_2_figures.pdf  # Figures S1-S3 (pre-processing variants,
                                 #   leave-one-out, library-normalisation)
  gene_panels.json               # exact membership of all six gene sets
  scripts/
    run_three_cohort.py          # gene-level testing, equivalence, meta-analysis,
                                 #   Hartung-Knapp intervals, C1+C2 re-pooling
    run_gsea.py                  # preranked GSEA, 10,000 gene-set permutations
    run_sensitivity.py           # leave-one-out, priming/execution split, variants
    run_round2_addendum.py       # HKSJ per gene, execution rank counts, observed power
    run_round2_normalisation.py  # raw / median-of-ratios / upper-quartile x filtering
    make_figure.py               # forest plot (three-cohort validation)
    make_gsea_figure.py          # NES heatmap + C1 running enrichment score
    make_figure1.py              # Figure 1 (design, positive control, Wnt panel)
    make_additional_file1.py     # builds Additional file 1
    make_additional_file2.py     # builds Additional file 2
    audit_manuscript.py          # 100 numeric checks of the manuscript against the results
  results/
    cohort_C{1,2,3}_*_gene_level.csv   # per-cohort gene-level statistics
    meta_focus.csv                     # random-effects and fixed-effect pooling
    meta_focus_hksj.csv                # Hartung-Knapp intervals, gene by gene
    sensitivity_two_cohort.csv         # C1+C2 re-pooling (C3 excluded)
    three_cohort_focus.csv             # long-format focus-gene table
    gsea_three_cohort.csv              # all cohorts x all sets x NES/q
    gsea_results.json                  # full GSEA output
    sensitivity_results.json           # all sensitivity analyses
    round2_results.json                # HKSJ, execution ranks, power/MDE
    normalisation_sensitivity.json     # 3 regimes x 2 filtering states x 3 cohorts
    run_summary.txt                    # console log of the main run
  figures/                       # Figure 1, forest plot, GSEA figure (PNG + PDF)

data/                            # v1.0.3 processed tables (legacy, see archive note)
figures/                         # v1.0.3 figures (legacy)
scripts/                         # v1.0.3 rebuild scripts (legacy)
archive_v1.0.3/                  # superseded manuscript and supplements
MANIFEST_sha256.txt              # SHA-256 of every file in this repository
requirements.txt                 # pinned package versions used for this release
```

## Data sources

- Bulk RNA-seq, rat spinal cord: GSE217932 (CCI, 4 vs 4), GSE256472 (SNI, 3 vs 3),
  GSE205906 (SNI, 3 vs 3) — Gene Expression Omnibus
- Single-cell RNA-seq, used only to count available astrocytes and microglia:
  GSE134003 (GEO)

Intervention arms were excluded before any comparison (CCI + levo-tetrahydropalmatine,
n = 4 in GSE217932; paeoniflorin–liquiritin, n = 3 in GSE205906).

## How to reproduce

Python 3.13.14 with the packages pinned in `requirements.txt` (numpy 2.5.4,
scipy 1.18.1, openpyxl 3.1.5, matplotlib 3.11.1; nothing imports pandas).
Every script resolves its output directory without hard-coded absolute paths: it
writes next to itself when run from this repository, and honours two environment
variables — `NP_OUT` (where to write results/figures) and `NP_DATA` (where the
GEO series matrices live, default `D:/省自然/生信分析/data`).

```
python three_cohort/scripts/run_three_cohort.py    # -> results/cohort_C*_gene_level.csv, meta_focus.csv,
                                                   #    sensitivity_two_cohort.csv
python three_cohort/scripts/run_gsea.py            # -> results/gsea_three_cohort.csv, gsea_results.json
python three_cohort/scripts/run_sensitivity.py     # -> results/sensitivity_results.json
python three_cohort/scripts/run_round2_addendum.py     # -> results/round2_results.json, meta_focus_hksj.csv
python three_cohort/scripts/run_round2_normalisation.py # -> results/normalisation_sensitivity.json
python three_cohort/scripts/make_figure1.py        # -> figures/Figure1_design_and_Wnt.{png,pdf}
python three_cohort/scripts/make_figure.py         # -> figures/Figure_three_cohort_validation.{png,pdf}
python three_cohort/scripts/make_gsea_figure.py    # -> figures/Figure_gsea_three_cohort.{png,pdf}
python three_cohort/scripts/make_additional_file1.py  # -> Additional_File_1_tables.xlsx
python three_cohort/scripts/make_additional_file2.py  # -> Additional_File_2_figures.{pdf,png}
python three_cohort/scripts/audit_manuscript.py    # -> 100/100 numeric checks against the results
```

`run_round2_addendum.py` and `run_round2_normalisation.py` read the per-cohort
gene-level CSVs produced by `run_three_cohort.py`; run that one first.
`audit_manuscript.py` re-derives every number it compares from the result files
rather than comparing against literals, so tampering with a result changes the
expectation as well as the observation — which is exactly what makes it useless
as a rubber stamp. A tamper test is documented in the review record.

Key settings, reported because they change output: gene-level collapse keeps the
transcript with the highest mean expression per symbol (never the smallest
*p*-value); Welch t-test with Welch–Satterthwaite df; Benjamini–Hochberg across
all tested genes; DerSimonian–Laird random effects with τ² estimated by the
DerSimonian–Laird moment estimator and intervals additionally computed under
Hartung–Knapp–Sidik–Jonkman; GSEA ranking metric is the Welch *t* winsorised at
±50 after removing genes below the cohort's lower-quartile mean expression;
10,000 gene-set permutations with the ranking held fixed (seed 20261029);
*q* values are Benjamini–Hochberg across the six sets tested in that cohort, not
genome-wide; the permutation floor is therefore *p* ≥ 1 × 10⁻⁴.

## Archive note

`archive_v1.0.3/` holds the superseded single-cohort manuscript and its
supplements. They are retained for provenance only. Their positive claims about
convergent Wnt/pyroptosis activation, Wnt enrichment and a 4.86-fold *Nlrp3*
induction do not replicate in any of the three cohorts analysed here, and the
4.86 figure was traced to a single low-expression transcript. **Do not cite the
archived version.**

## License

MIT — see LICENSE.
