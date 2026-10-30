# -*- coding: utf-8 -*-
"""
Preranked GSEA across the three independent cohorts.

Why this exists
---------------
The submitted manuscript reports a GSEA result (Hallmark Wnt/beta-catenin,
NES approx -1.6, nominal p approx 0.03, "not significant after correction").
That number was never recomputed on a documented ranking. This script rebuilds
GSEA from the gene-level cohort tables produced by run_three_cohort.py, using
one fixed procedure for all three cohorts, so the enrichment claim either holds
or does not.

Method (documented, defensible)
-------------------------------
1. Ranking metric: Welch t statistic (case vs control), descending. t is used
   instead of raw log2FC so that genes with large but noisy fold changes do not
   dominate the top of the list.
2. Weighted Kolmogorov-Smirnov enrichment score (Subramanian et al. 2005, p=1).
3. Null distribution: gene-set permutation, 10,000 draws of a random gene set of
   identical size drawn from the same ranked list. The ranking is never
   permuted (permuting the ranking, as in classic phenotype permutation, is not
   valid for preranked data with unequal group sizes).
4. NES = ES / mean(|ES_null| on the matching sign side).
5. Nominal p = fraction of null draws at least as extreme, one-sided on the
   observed sign.
6. FDR q = Benjamini-Hochberg across the gene sets tested in that cohort. This
   is explicit: it corrects across the SETS TESTED HERE (7), not across all of
   MSigDB. A genome-wide Hallmark scan would give a different, larger q.

Gene sets
---------
HALLMARK_WNT is taken verbatim from the panel hardcoded in the project's own
check_fig2_gsea.py (annotated there as MSigDB v7 human symbols; rat orthologs
share symbols). The remaining sets are curated canonical marker panels, NOT
MSigDB sets, and are labelled as such in the output.
"""
import csv, os, math, json

def _default_out():
    """Resolve the output directory without hard-coding a machine-specific
    path.  If the script lives in <pkg>/scripts/, write to <pkg>/; otherwise
    write next to the script.  Set NP_OUT to override."""
    here = os.path.dirname(os.path.abspath(__file__))
    parent = os.path.dirname(here)
    if os.path.basename(here).lower() == 'scripts' and os.path.isdir(parent):
        return parent
    return here



import numpy as np

OUT = os.environ.get('NP_OUT') or _default_out()
# Repository layout: per-cohort CSVs and GSEA outputs live in <pkg>/results/.
RESULTS = os.path.join(OUT, 'results') if os.path.isdir(os.path.join(OUT, 'results')) else OUT


def _find(fn):
    """Locate an input file in the results dir, the package dir, or the CWD."""
    for d in (RESULTS, OUT, os.getcwd()):
        p = os.path.join(d, fn)
        if os.path.exists(p):
            return p
    return os.path.join(RESULTS, fn)  # not found; fail with a clear path
COHORTS = [
    ('C1', 'cohort_C1_GSE217932_CCI_gene_level.csv', 'GSE217932  CCI vs Sham (n=4v4)'),
    ('C2', 'cohort_C2_GSE256472_SNI_gene_level.csv', 'GSE256472  SNI vs Ctrl (n=3v3)'),
    ('C3', 'cohort_C3_GSE205906_SNI_gene_level.csv', 'GSE205906  SNI vs Sham (n=3v3)'),
]

# ---------------------------------------------------------------- gene sets
# MSigDB v7 HALLMARK_WNT_BETA_CATENIN_SIGNALING (37 symbols; as used by the
# project's own check_fig2_gsea.py, unchanged)
WNT = ["MYC","CCND1","CTNNB1","LEF1","TCF7","FOSL1","JAG2","NOTUM","PPP2R5B","WIF1",
       "AXIN2","NKD1","TBL1XR1","CD44","KLF4","PAX2","CXXC4","SOX17","TCF7L1","TLE1",
       "USP34","GLIS2","SFRP1","DKK1","DKK4","FZD7","CER1","FST","ASCL2","LGR4","LRP6",
       "OLFM1","PITX2","SPRY2","T","WNT11","WNT3","WNT5A"]

# curated (not MSigDB): core pyroptosis machinery
PYRO = ["NLRP3","NLRC4","AIM2","PYCARD","CASP1","CASP4","CASP8","CASP3","GSDMD","GSDME",
        "IL1B","IL18","IL1A","HMGB1","TNF","IL6","NEK7","TXNIP","TLR4","NFKB1","RELA","GZMB"]

# curated (not MSigDB): injury-response programme, positive control panel.
# These are the canonical transcriptional markers of axonal injury / reactive
# gliosis in rodent spinal cord and DRG.
INJURY = ["ATF3","SOX11","SPRR1A","GAP43","GAL","NPY","VIP","CACNA2D1","CCL2","CCL7",
          "TIMP1","SERPINE1","LCN2","GFAP","VIM","NES","CD44","BDNF","PTGS2","MMP9",
          "ICAM1","SPP1","CD68"]

# curated (not MSigDB): microglial activation / disease-associated microglia
MICRO = ["AIF1","CX3CR1","TMEM119","ITGAM","CD68","C1QA","C1QB","C1QC","TYROBP","TREM2",
         "APOE","LGALS3","GPNMB","SPP1","CTSS","FCGR1A","LYZ2","CYBB","CD9","ITGB2"]

# curated (not MSigDB): astrocyte reactivity
ASTRO = ["GFAP","VIM","NES","LCN2","SERPINA3","GJA1","AQP4","SLC1A3","SLC1A2","S100B",
         "ALDH1L1","SOX9","STAT3","CXCL10","CCL2","TIMP1","CD44","MT1A","MT2A"]

# curated (not MSigDB): NF-kB / cytokine signalling (a Wnt-independent
# inflammatory readout; if any inflammatory programme moves, this one should)
NFKB = ["NFKB1","RELA","RELB","NFKB2","NFKBIA","NFKBIB","IL6","TNF","IL1B","CCL2","CXCL1",
        "CXCL2","ICAM1","VCAM1","PTGS2","NOS2","MMP9","TLR2","TLR4","MYD88","TRAF6","IRAK1"]

SETS = [
    ('HALLMARK_WNT_BETA_CATENIN', WNT, 'MSigDB v7 Hallmark (37 genes)'),
    ('PYROPTOSIS_CORE', PYRO, 'curated core machinery (22)'),
    ('INJURY_RESPONSE_POSCTRL', INJURY, 'curated injury programme (23) - POSITIVE CONTROL'),
    ('MICROGLIA_ACTIVATION', MICRO, 'curated microglial activation (20)'),
    ('ASTROCYTE_REACTIVITY', ASTRO, 'curated astrocyte reactivity (19)'),
    ('NFKB_CYTOKINE_SIGNALLING', NFKB, 'curated NF-kB / cytokine (22)'),
]

NP_PERM = 10000
SEED = 20261029


def read_cohort(fn, expr_quantile=0.25, t_cap=50.0):
    """Load a cohort table, drop low-expression genes, rank by Welch t.

    Expression filtering is a mandatory GSEA pre-processing step, not a
    convenience: without it the ranking is occupied by near-zero-variance
    transcripts whose t statistics explode (observed here: t = 1.3e16 for a
    gene with mean log2 expression 0.03) and every gene set then appears
    'enriched' or 'depleted' purely as an artefact. We keep genes at or above
    the cohort's own lower-quartile mean expression.
    """
    raw = []
    with open(_find(fn), encoding='utf-8') as fh:
        for r in csv.DictReader(fh):
            try:
                t = float(r['t_stat'])
                p = float(r['pvalue'])
                lfc = float(r['log2FC'])
                expr = float(r['mean_expr'])
            except (ValueError, KeyError):
                continue
            if math.isnan(t) or math.isnan(lfc) or math.isnan(expr):
                continue
            raw.append((r['symbol'].upper(), t, lfc, p, expr))
    ex = sorted(x[4] for x in raw)
    thr = ex[int(len(ex) * expr_quantile)]
    rows = [(s, max(-t_cap, min(t_cap, t)), l, p) for s, t, l, p, e in raw if e >= thr]
    rows.sort(key=lambda x: -x[1])
    return rows, thr, len(raw)


def es_vectorized(metric_sorted, hit_mask):
    """Weighted KS enrichment score (Subramanian 2005, p=1)."""
    w = np.where(hit_mask, np.abs(metric_sorted), 0.0)
    sw = w.sum()
    nh = int(hit_mask.sum())
    nm = len(metric_sorted) - nh
    if sw == 0 or nh == 0 or nm == 0:
        return 0.0
    cum = np.cumsum(w) / sw - np.cumsum(~hit_mask) / nm
    i = int(np.argmax(np.abs(cum)))
    return float(cum[i])


def run_gsea(rows, setname, geneset, nperm=NP_PERM, seed=SEED):
    syms = [r[0] for r in rows]
    metric = np.array([r[1] for r in rows], dtype=float)
    setu = {g.upper() for g in geneset}
    hit = np.array([s in setu for s in syms])
    k = int(hit.sum())
    if k < 5:
        return None
    es_obs = es_vectorized(metric, hit)
    rng = np.random.default_rng(seed)
    n = len(syms)
    null = np.empty(nperm)
    for i in range(nperm):
        idx = rng.choice(n, size=k, replace=False)
        m = np.zeros(n, dtype=bool)
        m[idx] = True
        null[i] = es_vectorized(metric, m)
    if es_obs >= 0:
        side = null[null >= 0]
        nes = es_obs / side.mean() if len(side) and side.mean() != 0 else float('nan')
        pval = (np.sum(null >= es_obs) + 1) / (nperm + 1)
    else:
        side = null[null < 0]
        nes = es_obs / abs(side.mean()) if len(side) and side.mean() != 0 else float('nan')
        pval = (np.sum(null <= es_obs) + 1) / (nperm + 1)
    return dict(set=setname, k=k, N=n, es=es_obs, nes=float(nes), p=float(pval),
                frac_genes=100.0 * k / len(geneset))


def bh(pvals):
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    adj = [0.0] * m
    prev = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        val = min(prev, pvals[i] * m / rank)
        adj[i] = val
        prev = val
    return adj


def main():
    log = []

    def say(s):
        print(s, flush=True)
        log.append(s)

    say('=' * 92)
    say('PRERANKED GSEA — three cohorts, one procedure')
    say(f'ranking = Welch t; {NP_PERM} gene-set permutations; BH across the 6 sets tested')
    say('=' * 92)

    all_res = {}
    for cid, fn, label in COHORTS:
        rows, thr, n_raw = read_cohort(fn)
        say(f'\n--- {cid} {label} ---')
        say(f'  expression filter: mean log2 expr >= {thr:.2f}; '
            f'{n_raw} -> {len(rows)} genes retained')
        say(f'  top 5 after filtering: ' + ', '.join(f'{g}(t={t:.1f})' for g, t, _, _ in rows[:5]))
        res = []
        for sname, genes, note in SETS:
            r = run_gsea(rows, sname, genes)
            if r is None:
                say(f'  {sname:28s} SKIPPED (<5 genes present)')
                continue
            r['note'] = note
            res.append(r)
        qs = bh([r['p'] for r in res])
        for r, q in zip(res, qs):
            r['q'] = q
        all_res[cid] = res
        say(f'  {"gene set":30s} {"k":>4s} {"ES":>7s} {"NES":>7s} {"nom p":>8s} {"FDR q":>8s}')
        for r in sorted(res, key=lambda x: -abs(x['nes'])):
            say(f'  {r["set"]:30s} {r["k"]:4d} {r["es"]:+7.3f} {r["nes"]:+7.2f} '
                f'{r["p"]:8.4f} {r["q"]:8.3f}')

    # -------- cross-cohort consistency table
    say('\n' + '=' * 92)
    say('CROSS-COHORT CONSISTENCY (NES per cohort)')
    say('=' * 92)
    say(f'{"gene set":30s} {"C1 NES":>8s} {"C2 NES":>8s} {"C3 NES":>8s} {"all q<0.05?":>12s}')
    for sname, _, _ in SETS:
        line = f'{sname:30s}'
        nes_vals, sig = [], []
        for cid, _, _ in COHORTS:
            r = next((x for x in all_res[cid] if x['set'] == sname), None)
            if r is None:
                line += f'{"n/a":>8s}'
                continue
            line += f'{r["nes"]:+8.2f}'
            nes_vals.append(r['nes'])
            sig.append(r['q'] < 0.05)
        consistent = 'yes' if all(sig) and len(sig) == 3 else ('no' if sig else 'n/a')
        line += f'{consistent:>12s}'
        say(line)

    # -------- write csv
    path = os.path.join(RESULTS, 'gsea_three_cohort.csv')
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['cohort', 'gene_set', 'set_source', 'n_genes_in_set', 'n_ranked',
                    'ES', 'NES', 'nominal_p', 'FDR_q_across_sets_in_cohort', 'pct_of_set_found'])
        for cid, _, _ in COHORTS:
            for r in all_res[cid]:
                w.writerow([cid, r['set'], r['note'], r['k'], r['N'],
                            f"{r['es']:.4f}", f"{r['nes']:.3f}", f"{r['p']:.4e}",
                            f"{r['q']:.4f}", f"{r['frac_genes']:.1f}"])
    say(f'\nwrote {path}')

    # -------- headline
    say('\n' + '=' * 92)
    say('HEADLINE')
    say('=' * 92)
    for sname in ['INJURY_RESPONSE_POSCTRL', 'HALLMARK_WNT_BETA_CATENIN', 'PYROPTOSIS_CORE']:
        vals = []
        for cid, _, _ in COHORTS:
            r = next((x for x in all_res[cid] if x['set'] == sname), None)
            if r:
                vals.append(f'{cid}: NES={r["nes"]:+.2f} p={r["p"]:.3f} q={r["q"]:.2f}')
        say(f'  {sname}')
        for v in vals:
            say(f'    {v}')

    with open(os.path.join(RESULTS, 'gsea_summary.txt'), 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(log) + '\n')
    with open(os.path.join(RESULTS, 'gsea_results.json'), 'w', encoding='utf-8') as fh:
        json.dump({cid: all_res[cid] for cid in all_res}, fh, indent=2)
    say('\nwrote gsea_summary.txt, gsea_results.json')


if __name__ == '__main__':
    main()
