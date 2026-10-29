# astrocytic-wnt-pyroptosis-np

**The Wnt/β-catenin transcriptional programme is not induced in neuropathic pain
spinal cord: evidence from three independent rat cohorts.**

Single-author study (Yongxin Yang, ORCID 0009-0004-9698-6552). All analyses use
publicly available, de-identified rodent data; no new experimental data were
generated.

This repository is the reproducibility package for the manuscript. The current
release is **v1.1.2**, which replaces the earlier single-cohort analysis in
`archive_v1.0.3/`.

## What the study claims, and what it does not

1. **The Wnt/β-catenin transcriptional programme is not induced** after nerve
   injury in rat spinal cord. Across a 37-gene curated panel, the largest pooled
   effect over three cohorts was +0.27 log2 (*Ccnd2*, 1.20-fold); 31 of 36 pooled
   confidence intervals lay wholly within a pre-specified ±0.5 log2 equivalence
   bound. Set-level testing with HALLMARK_WNT_BETA_CATENIN was non-significant in
   every cohort and under four pre-processing variants.
2. **Transcriptome-wide false-discovery control is unattainable at n = 3–4.**
   *Atf3*, an unambiguously induced injury gene, rose 16.3-fold and still reached
   only *q* = 0.466; no gene survived correction in two of three cohorts, and the
   single survivor in the third was a zero-variance transcript.
3. **Pyroptosis transcripts are primed but the execution machinery is not
   engaged**, and this pattern was confined to the CCI cohort. It did not
   replicate in the SNI cohort.

The study does **not** claim that Wnt signalling is inactive: β-catenin
activation is largely post-translational and invisible to bulk transcriptomics.
See the manuscript Discussion for the exact scope of the negative result.

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

## Repository structure

```
three_cohort/                    # v1.1.2 analysis (this release)
  MANUSCRIPT_BMC.md              # manuscript source, BMC format
  MANUSCRIPT_BMC.docx            # submission-ready Word file
  Additional_File_1_tables.xlsx  # Tables S1-S4 (per-cohort stats, pooled, GSEA, panels)
  Additional_File_2_figures.pdf  # Figures S1-S2 (pre-processing variants, leave-one-out)
  gene_panels.json               # exact membership of all six gene sets
  scripts/
    run_three_cohort.py          # gene-level testing, equivalence, meta-analysis
    run_gsea.py                  # preranked GSEA, 10,000 gene-set permutations
    run_sensitivity.py           # leave-one-out, priming/execution split, variants
    make_figure.py               # forest plot (three-cohort validation)
    make_gsea_figure.py          # NES heatmap + C1 running enrichment score
    make_figure1.py              # Figure 1 (design, positive control, Wnt panel)
    make_additional_file1.py     # builds Additional file 1
    make_additional_file2.py     # builds Additional file 2
    audit_manuscript.py          # checks every number in the manuscript against the results
  results/
    cohort_C{1,2,3}_*_gene_level.csv   # per-cohort gene-level statistics
    meta_focus.csv                     # random-effects and fixed-effect pooling
    three_cohort_focus.csv             # long-format focus-gene table
    gsea_three_cohort.csv              # all cohorts x all sets x NES/q
    gsea_results.json                  # full GSEA output
    sensitivity_results.json           # all sensitivity analyses
    run_summary.txt                    # console log of the main run
  figures/                       # Figure 1, forest plot, GSEA figure (PNG + PDF)

data/                            # v1.0.3 processed tables (legacy, see archive note)
figures/                         # v1.0.3 figures (legacy)
scripts/                         # v1.0.3 rebuild scripts (legacy)
archive_v1.0.3/                  # superseded manuscript and supplements
MANIFEST_sha256.txt              # SHA-256 of every file in this repository
```

## Data sources

- Bulk RNA-seq, rat spinal cord: GSE217932 (CCI, 4 vs 4), GSE256472 (SNI, 3 vs 3),
  GSE205906 (SNI, 3 vs 3) — Gene Expression Omnibus
- Single-cell RNA-seq, used only to count available astrocytes and microglia:
  GSE134003 (GEO)

Intervention arms were excluded before any comparison (CCI + levo-tetrahydropalmatine,
n = 4 in GSE217932; paeoniflorin–liquiritin, n = 3 in GSE205906).

## How to reproduce

Python 3.11 with NumPy and SciPy. Place the three GEO series matrices in a `data/`
directory and adjust the paths at the top of each script.

```
python three_cohort/scripts/run_three_cohort.py    # -> results/cohort_C*_gene_level.csv, meta_focus.csv
python three_cohort/scripts/run_gsea.py            # -> results/gsea_three_cohort.csv, gsea_results.json
python three_cohort/scripts/run_sensitivity.py     # -> results/sensitivity_results.json
python three_cohort/scripts/make_figure1.py        # -> figures/Figure1_design_and_Wnt.{png,pdf}
python three_cohort/scripts/make_figure.py         # -> figures/Figure_three_cohort_validation.{png,pdf}
python three_cohort/scripts/make_gsea_figure.py    # -> figures/Figure_gsea_three_cohort.{png,pdf}
python three_cohort/scripts/make_additional_file1.py  # -> Additional_File_1_tables.xlsx
python three_cohort/scripts/make_additional_file2.py  # -> Additional_File_2_figures.{pdf,png}
python three_cohort/scripts/audit_manuscript.py    # -> 45/45 numeric checks against the results
```

Key settings, reported because they change output: gene-level collapse keeps the
transcript with the highest mean expression per symbol (never the smallest
*p*-value); Welch t-test with Welch–Satterthwaite df; Benjamini–Hochberg across
all tested genes; DerSimonian–Laird random effects; GSEA ranking metric is the
Welch *t* winsorised at ±50 after removing genes below the cohort's lower-quartile
mean expression; 10,000 gene-set permutations with the ranking held fixed;
*q* values are Benjamini–Hochberg across the six sets tested in that cohort, not
genome-wide.

## Archive note

`archive_v1.0.3/` holds the superseded single-cohort manuscript and its
supplements. They are retained for provenance only. Their positive claims about
convergent Wnt/pyroptosis activation, Wnt enrichment and a 4.86-fold *Nlrp3*
induction do not replicate in any of the three cohorts analysed here, and the
4.86 figure was traced to a single low-expression transcript. **Do not cite the
archived version.**

## License

MIT — see LICENSE.
