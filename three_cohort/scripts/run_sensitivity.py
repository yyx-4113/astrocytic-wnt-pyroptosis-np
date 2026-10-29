# -*- coding: utf-8 -*-
"""
Sensitivity analysis for the GSEA result reported in GSEA_ADDENDUM.md.

Three questions, each of which a reviewer would ask and none of which the
headline NES answers on its own:

Q1  Is the C1 pyroptosis enrichment carried by one or two genes?
    -> leave-one-out on the leading edge: drop each LE gene, recompute.
    -> also a "no leading edge" run that drops all of them at once.

Q2  Is the enrichment execution machinery or priming substrate?
    -> split the 22-gene curated set into disjoint EXECUTION (pore-forming /
       inflammasome core) and PRIMING (NF-kB-dependent inducible substrates)
       halves and run each separately. If only PRIMING moves, the word
       "pyroptosis" is being used for what is really signal-1 priming.

Q3  Does the conclusion depend on arbitrary analysis choices?
    -> rerun with the ranking metric changed (Welch t -> log2FC) and with the
       expression filter relaxed (lower quartile -> median). A conclusion that
       survives neither is not a conclusion.

Everything reuses the exact functions in run_gsea.py (same ES definition, same
10,000-draw gene-set permutation null, same seed) so the numbers are directly
comparable to gsea_three_cohort.csv.
"""
import os, sys, csv, json, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import run_gsea as G

OUT = G.OUT
SAY = []


def say(s=''):
    print(s, flush=True)
    SAY.append(s)


# ---------------------------------------------------------------- gene sets
# Disjoint split of the 22-gene curated PYROPTOSIS_CORE set.
# EXECUTION: the machinery that actually forms the pore / activates caspase-1.
# Everything here must move for "pyroptosis" to be the correct word.
EXEC = ["NLRP3", "NLRC4", "AIM2", "PYCARD", "NEK7", "CASP1", "CASP4", "GSDMD", "GSDME"]

# PRIMING: signal-1 / NF-kB-dependent inducible substrates and cytokines.
# These rise whenever ANY inflammatory stimulus is present and say nothing
# about whether the pore is being assembled.
PRIME = ["NFKB1", "RELA", "TLR4", "IL1B", "IL18", "IL1A", "TNF", "IL6",
         "HMGB1", "TXNIP", "CASP8", "CASP3", "GZMB"]

assert not (set(EXEC) & set(PRIME)), 'split must be disjoint'
assert set(EXEC) | set(PRIME) == set(G.PYRO), 'split must cover the curated set exactly'


def run_one(rows, name, genes, nperm=G.NP_PERM, seed=G.SEED):
    r = G.run_gsea(rows, name, genes, nperm=nperm, seed=seed)
    return r


def fmt(r):
    if r is None:
        return 'k<5, not testable'
    return (f'ES={r["es"]:+.4f}  NES={r["nes"]:+.3f}  p={r["p"]:.4f}  k={r["k"]}/{r["N"]}')


def leading_edge(rows, genes):
    """Return LE genes and their ranks, computed exactly as es_vectorized does."""
    syms = [r[0] for r in rows]
    metric = np.array([r[1] for r in rows], dtype=float)
    setu = {g.upper() for g in genes}
    hit = np.array([s in setu for s in syms])
    w = np.where(hit, np.abs(metric), 0.0)
    sw = w.sum()
    nh = int(hit.sum())
    nm = len(metric) - nh
    cum = np.cumsum(w) / sw - np.cumsum(~hit) / nm
    i = int(np.argmax(np.abs(cum)))
    es = float(cum[i])
    if es >= 0:
        idx = [j for j in range(i + 1) if hit[j]]
    else:
        idx = [j for j in range(i, len(syms)) if hit[j]]
    return es, i, [(syms[j], j + 1) for j in idx], cum


def main():
    c1_rows, c1_thr, _ = G.read_cohort('cohort_C1_GSE217932_CCI_gene_level.csv')
    n_all = len(c1_rows)

    say('=' * 96)
    say('SENSITIVITY ANALYSIS — C1 (GSE217932, CCI vs Sham, n=4v4)')
    say(f'ranked genes after expression filter: {n_all}  (mean log2 expr >= {c1_thr:.2f})')
    say('=' * 96)

    res = {}

    # ---------------------------------------------------------------- Q1 LOO
    say('\n[Q1] Leave-one-out on the pyroptosis leading edge')
    es_full, peak, le, cum = leading_edge(c1_rows, G.PYRO)
    base = run_one(c1_rows, 'PYROPTOSIS_CORE', G.PYRO)
    say(f'  full set: {fmt(base)}   peak rank={peak + 1}   LE={len(le)}/{base["k"]}')
    say(f'  LE genes (symbol, rank): ' + ', '.join(f'{g}({r})' for g, r in le))

    loo = []
    say(f'\n  {"gene dropped":12s} {"k":>3s} {"ES":>9s} {"NES":>8s} {"p":>8s} {"dNES":>7s}  verdict')
    for g, rk in le:
        sub = [x for x in G.PYRO if x != g]
        r = run_one(c1_rows, 'PYRO_no_' + g, sub)
        if r is None:
            say(f'  {g:12s}  -     not testable')
            continue
        d = r['nes'] - base['nes']
        verdict = 'still q<0.05' if r['p'] < 0.05 else 'DROPS to p>=0.05'
        say(f'  {g:12s} {r["k"]:3d} {r["es"]:+9.4f} {r["nes"]:+8.3f} {r["p"]:8.4f} {d:+7.3f}  {verdict}')
        loo.append(dict(dropped=g, rank=rk, k=r['k'], es=r['es'], nes=r['nes'],
                        p=r['p'], dNES=d))
    res['loo'] = loo

    # drop all LE at once — the honest "how much of this is the LE" question
    le_genes = {g for g, _ in le}
    sub = [x for x in G.PYRO if x not in le_genes]
    r_nole = run_one(c1_rows, 'PYRO_no_leadingEdge', sub)
    say(f'\n  drop all {len(le)} LE genes at once  ->  {fmt(r_nole)}')
    res['no_le'] = r_nole

    # ------------------------------------------------------- Q2 exec vs prime
    say('\n[Q2] Execution machinery vs priming substrate (disjoint halves)')
    r_exec = run_one(c1_rows, 'PYRO_EXECUTION_MACHINERY', EXEC)
    r_prime = run_one(c1_rows, 'PYRO_PRIMING_SUBSTRATE', PRIME)
    say(f'  EXECUTION ({len(EXEC)} genes: {", ".join(EXEC)})')
    say(f'      {fmt(r_exec)}')
    say(f'  PRIMING   ({len(PRIME)} genes: {", ".join(PRIME)})')
    say(f'      {fmt(r_prime)}')
    res['exec'] = r_exec
    res['prime'] = r_prime

    # per-gene position of every set member, so the split is auditable
    syms = [x[0] for x in c1_rows]
    rank = {s: i + 1 for i, s in enumerate(syms)}
    say('\n  per-gene position in the C1 ranking (percentile of 1 = top):')
    say(f'  {"gene":9s} {"rank":>7s} {"pctile":>7s} {"log2FC":>8s} {"t":>7s}   group')
    c1_by = {x[0]: x for x in c1_rows}
    for g in sorted(set(EXEC) | set(PRIME), key=lambda y: rank.get(y, 10 ** 9)):
        if g not in c1_by:
            say(f'  {g:9s} {"absent from filtered table":>40s}')
            continue
        _, t, lfc, _ = c1_by[g]
        pc = 100.0 * rank[g] / len(syms)
        grp = 'EXEC' if g in EXEC else 'PRIME'
        say(f'  {g:9s} {rank[g]:7d} {pc:6.1f}% {lfc:+8.2f} {t:+7.2f}   {grp}')

    # ------------------------------------------------ Q3 analysis-choice tests
    say('\n[Q3] Dependence on analysis choices')

    # (a) ranking metric: t -> log2FC (re-read with the same expression filter)
    raw = []
    with open(os.path.join(OUT, 'cohort_C1_GSE217932_CCI_gene_level.csv'), encoding='utf-8') as fh:
        for r in csv.DictReader(fh):
            try:
                lfc = float(r['log2FC']); expr = float(r['mean_expr'])
            except (ValueError, KeyError):
                continue
            if math.isnan(lfc) or math.isnan(expr):
                continue
            raw.append((r['symbol'].upper(), lfc, expr))
    alt = [(s, l, 0.0, 0.0) for s, l, e in raw if e >= c1_thr]
    alt.sort(key=lambda x: -x[1])
    r_lfc_pyro = run_one(alt, 'PYRO_lfc_ranked', G.PYRO)
    r_lfc_exec = run_one(alt, 'EXEC_lfc_ranked', EXEC)
    r_lfc_wnt = run_one(alt, 'WNT_lfc_ranked', G.WNT)
    say(f'  ranking by log2FC instead of Welch t ({len(alt)} genes):')
    say(f'    PYROPTOSIS_CORE  {fmt(r_lfc_pyro)}')
    say(f'    EXECUTION        {fmt(r_lfc_exec)}')
    say(f'    HALLMARK_WNT     {fmt(r_lfc_wnt)}')
    res['lfc'] = dict(pyro=r_lfc_pyro, exec=r_lfc_exec, wnt=r_lfc_wnt)

    # (b) expression filter relaxed to the median
    rows_med, thr_med, _ = G.read_cohort('cohort_C1_GSE217932_CCI_gene_level.csv',
                                         expr_quantile=0.50)
    r_med_pyro = run_one(rows_med, 'PYRO_median_filter', G.PYRO)
    r_med_exec = run_one(rows_med, 'EXEC_median_filter', EXEC)
    r_med_wnt = run_one(rows_med, 'WNT_median_filter', G.WNT)
    say(f'  expression filter at the median (>= {thr_med:.2f}, {len(rows_med)} genes):')
    say(f'    PYROPTOSIS_CORE  {fmt(r_med_pyro)}')
    say(f'    EXECUTION        {fmt(r_med_exec)}')
    say(f'    HALLMARK_WNT     {fmt(r_med_wnt)}')
    res['med'] = dict(pyro=r_med_pyro, exec=r_med_exec, wnt=r_med_wnt)

    # (c) no expression filter at all, to show what the artefact looks like
    rows_nof, _, _ = G.read_cohort('cohort_C1_GSE217932_CCI_gene_level.csv',
                                   expr_quantile=0.0)
    r_nf_pyro = run_one(rows_nof, 'PYRO_no_filter', G.PYRO)
    r_nf_exec = run_one(rows_nof, 'EXEC_no_filter', EXEC)
    say(f'  no expression filter ({len(rows_nof)} genes — the artefact condition):')
    say(f'    PYROPTOSIS_CORE  {fmt(r_nf_pyro)}')
    say(f'    EXECUTION        {fmt(r_nf_exec)}')
    res['nofilter'] = dict(pyro=r_nf_pyro, exec=r_nf_exec)

    # ------------------------------------------- same split in C2 and C3
    say('\n[Q2b] The same execution/priming split in the two other cohorts')
    for cid, fn in [('C2', 'cohort_C2_GSE256472_SNI_gene_level.csv'),
                    ('C3', 'cohort_C3_GSE205906_SNI_gene_level.csv')]:
        rows, thr, _ = G.read_cohort(fn)
        a = run_one(rows, 'PYRO_EXECUTION_MACHINERY', EXEC)
        b = run_one(rows, 'PYRO_PRIMING_SUBSTRATE', PRIME)
        say(f'  {cid} (filter {thr:.2f}, {len(rows)} genes)')
        say(f'    EXECUTION  {fmt(a)}')
        say(f'    PRIMING    {fmt(b)}')
        res[cid] = dict(exec=a, prime=b)

    # ---------------------------------------------------------------- verdict
    say('\n' + '=' * 96)
    say('VERDICT')
    say('=' * 96)
    worst = min(x['nes'] for x in loo)
    n_drop = sum(1 for x in loo if x['p'] >= 0.05)
    say(f'  LOO: {len(loo)} leading-edge genes dropped singly; NES range '
        f'{worst:+.2f} to {base["nes"]:+.2f}; {n_drop} of {len(loo)} single '
        f'drops push p above 0.05')
    say(f'  drop all LE: NES {r_nole["nes"]:+.2f}, p={r_nole["p"]:.3f}')
    say(f'  execution vs priming: NES {r_exec["nes"] if r_exec else float("nan"):+.2f} '
        f'vs {r_prime["nes"]:+0.2f}' if r_exec else '  execution: not testable')
    say('  -> see SENSITIVITY_ADDENDUM.md for the interpretation')

    with open(os.path.join(OUT, 'sensitivity_summary.txt'), 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(SAY) + '\n')
    with open(os.path.join(OUT, 'sensitivity_results.json'), 'w', encoding='utf-8') as fh:
        json.dump(res, fh, indent=2, default=str)
    print('\nwrote sensitivity_summary.txt, sensitivity_results.json')


if __name__ == '__main__':
    main()
