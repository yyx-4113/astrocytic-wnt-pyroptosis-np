#!/usr/bin/env python3
"""Round-2 review: is the "no gene survives correction" result a property of
the sample size or of our pipeline?

R3-M3 pointed out that the primary pipeline applies log2(x+1) straight to raw
counts, with no library-size normalisation and no independent filtering.  Both
are standard and both increase power, so the review required that this
alternative explanation be excluded or quantified rather than left implicit.

This script re-runs the C1 and C2 cohorts from their raw count matrices under
four pre-processing regimes and reports, for each, the number of genes
surviving Benjamini-Hochberg correction and the effect on Atf3 and on the Wnt
panel.  C3 is already normalised abundance (RPKM) and cannot be renormalised;
only independent filtering is applied to it.

Regimes
  raw            log2(count + 1)                      - the primary analysis
  median_ratio   DESeq2-style median-of-ratios size factors
  upper_quartile upper-quartile size factors
  each of the above optionally with independent filtering before BH

Outputs: normalisation_sensitivity.json
"""
import csv
import gzip
import json
import math
import os
import pickle
import sys

import numpy as np
import openpyxl
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_three_cohort as R  # noqa: E402  (WNT panel, collapse, Welch, BH)

DATA = os.environ.get('NP_DATA') or r'D:/省自然/生信分析/data'
OUT = R.RESULTS  # <pkg>/results in the repository layout; script dir when flat
WNT = {s.upper() for s in R.WNT}
ATF3 = 'ATF3'


# --------------------------------------------------------------- loaders
def load_c1():
    src = os.path.join(DATA, 'GSE217932_Raw_gene_counts_matrix.xlsx')
    wb = openpyxl.load_workbook(src, read_only=True)
    ws = wb.active
    it = ws.iter_rows(values_only=True)
    header = list(next(it))
    cols = header[1:]
    sham = [i for i, c in enumerate(cols) if 'Sham' in str(c)]
    cci = [i for i, c in enumerate(cols)
           if str(c).startswith('CCI-') and 'THP' not in str(c)]
    gmap = pickle.load(open(os.path.join(DATA, 'rat_gene_map.pkl'), 'rb'))
    gids, mat = [], []
    for row in it:
        gid = row[0]
        if gid is None:
            continue
        sym = gmap.get(str(gid).split('.')[0])
        if not sym:
            continue
        gids.append(sym)
        mat.append([0.0 if v is None else float(v) for v in row[1:]])
    return np.array(mat), gids, sham, cci


def load_c2():
    src = os.path.join(DATA, 'GSE256472_all.genes.expression.annot.xls.gz')
    with gzip.open(src, 'rt', encoding='utf-8', errors='replace') as fh:
        header = fh.readline().rstrip('\n').split('\t')
        cols = [c for c in header if c.endswith('_count')]
        idx = [header.index(c) for c in cols]
        sym_i = header.index('Symbol')
        ctrl = [i for i, c in enumerate(cols) if c.startswith('Control')]
        sni = [i for i, c in enumerate(cols) if c.startswith('SNI')]
        gids, mat = [], []
        for line in fh:
            f = line.rstrip('\n').split('\t')
            sym = f[sym_i].strip()
            if not sym:
                continue
            gids.append(sym)
            mat.append([float(f[i]) if f[i] not in ('', 'NA') else 0.0 for i in idx])
    return np.array(mat), gids, ctrl, sni


def load_c3():
    """Already-normalised abundance; returned on the log2 scale as deposited."""
    src = os.path.join(DATA, 'GSE205906_processed_data.xls.gz')
    with gzip.open(src, 'rt', encoding='utf-8', errors='replace') as fh:
        header = fh.readline().rstrip('\n').split('\t')
        ja = [i for i, c in enumerate(header) if c.startswith('JA')]
        sni = [i for i, c in enumerate(header) if c.startswith('SNI')]
        gids, mat = [], []
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 2:
                continue
            sym = f[1].strip() if len(f) > 1 else ''
            if not sym:
                continue
            gids.append(sym)
            try:
                mat.append([math.log2(float(f[i]) + 1) if float(f[i]) > 0 else 0.0
                            for i in ja + sni])
            except ValueError:
                continue
    return np.array(mat), gids, list(range(len(ja))), [len(ja) + i for i in range(len(sni))]


# ---------------------------------------------------- normalisation schemes
def size_factors(mat):
    """DESeq2 median-of-ratios and upper-quartile factors."""
    pos = np.all(mat > 0, axis=1)
    ref = np.exp(np.mean(np.log(mat[pos]), axis=1))
    ratios = mat[pos] / ref[:, None]
    mor = np.median(ratios, axis=0)
    nz = np.where(mat > 0, mat, np.nan)
    uq = np.nanpercentile(nz, 75, axis=0)
    uq = uq / np.exp(np.mean(np.log(uq)))
    return {'median_ratio': mor / np.exp(np.mean(np.log(mor))),
            'upper_quartile': uq}


def collapse(genes, logmat):
    """Highest mean expression per symbol - the primary pipeline's rule."""
    best = {}
    for i, g in enumerate(genes):
        m = float(np.mean(logmat[i]))
        if g not in best or m > best[g][1]:
            best[g] = (logmat[i], m)
    return {g: v[0] for g, v in best.items()}, {g: v[1] for g, v in best.items()}


def welch(a, b):
    na, nb = len(a), len(b)
    va, vb = np.var(a, ddof=1), np.var(b, ddof=1)
    se = math.sqrt(va / na + vb / nb)
    if not np.isfinite(se) or se <= 0:
        return None
    denom = (va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1)
    if denom <= 0:
        return None
    df = (va / na + vb / nb) ** 2 / denom
    if not np.isfinite(df) or df <= 0:
        return None
    t = float(np.mean(a) - np.mean(b)) / se
    p = float(2 * stats.t.sf(abs(t), df))
    lfc = float(np.mean(a) - np.mean(b))
    crit = float(stats.t.ppf(0.975, df))
    return dict(log2FC=lfc, t=t, df=float(df), p=p, se=se,
                ci=(lfc - crit * se, lfc + crit * se))


def bh(ps):
    ps = np.asarray(ps, dtype=float)
    n = len(ps)
    order = np.argsort(ps)
    q = np.empty(n)
    prev = 1.0
    for rank, i in enumerate(order[::-1]):
        r = n - rank
        prev = min(prev, ps[i] * n / r)
        q[i] = prev
    return q


def run_cohort(mat, genes, ctrl_idx, case_idx, label, normalised_deposited=False):
    """Return per-regime summaries for one cohort."""
    out = {}
    if normalised_deposited:
        regimes = {'as_deposited': mat}
    else:
        sf = size_factors(mat)
        regimes = {'raw': np.log2(mat + 1),
                   'median_ratio': np.log2(mat / sf['median_ratio'] + 1),
                   'upper_quartile': np.log2(mat / sf['upper_quartile'] + 1)}
        out['size_factors'] = {k: [round(float(x), 4) for x in v] for k, v in sf.items()}
        out['library_totals'] = [int(x) for x in mat.sum(axis=0)]

    for regime, logmat in regimes.items():
        genes_c, mean_expr = collapse(genes, logmat)
        res = {}
        for g, vec in genes_c.items():
            r = welch(vec[case_idx], vec[ctrl_idx])
            if r is not None:
                res[g] = r
        syms = list(res)
        ps = [res[g]['p'] for g in syms]
        qs = bh(ps)
        for g, qv in zip(syms, qs):
            res[g]['q'] = float(qv)
            res[g]['mean_expr'] = mean_expr[g]

        for filt in (False, True):
            keep = syms
            if filt:
                thr = float(np.median([mean_expr[g] for g in syms]))
                keep = [g for g in syms if mean_expr[g] >= thr]
                qf = bh([res[g]['p'] for g in keep])
                for g, qv in zip(keep, qf):
                    res[g]['q_filt'] = float(qv)
            qkey = 'q_filt' if filt else 'q'
            n_surv = sum(1 for g in keep if res[g][qkey] < 0.05)
            n_nom = sum(1 for g in keep if res[g]['p'] < 0.05)
            wnt = {g: round(res[g]['log2FC'], 4) for g in keep if g.upper() in WNT}
            med_lfc = float(np.median([res[g]['log2FC'] for g in keep]))
            atf3 = next((g for g in keep if g.upper() == ATF3), None)
            survivors = sorted((g for g in keep if res[g][qkey] < 0.05),
                               key=lambda g: res[g][qkey])
            entry = dict(n_genes=len(keep), nominal_p05=n_nom, q05_survivors=n_surv,
                         median_log2FC=med_lfc,
                         Atf3=(dict(log2FC=round(res[atf3]['log2FC'], 4),
                                    p=float(res[atf3]['p']), q=float(res[atf3][qkey]),
                                    fold=round(2 ** res[atf3]['log2FC'], 3))
                                if atf3 else None),
                         survivors=survivors[:10],
                         Wnt_n=len(wnt))
            out[f'{regime}{"_indepfiltered" if filt else ""}'] = entry
            out.setdefault('_wnt', {})[f'{regime}{"_indepfiltered" if filt else ""}'] = wnt
    out['_label'] = label
    return out


def main():
    report = {}
    mat, genes, ctrl, case = load_c1()
    print(f'C1 raw counts matrix {mat.shape}, sham={len(ctrl)} cci={len(case)}')
    report['C1'] = run_cohort(mat, genes, ctrl, case, 'C1 GSE217932 CCI')
    mat, genes, ctrl, case = load_c2()
    print(f'C2 raw counts matrix {mat.shape}, control={len(ctrl)} sni={len(case)}')
    report['C2'] = run_cohort(mat, genes, ctrl, case, 'C2 GSE256472 SNI')
    mat, genes, ctrl, case = load_c3()
    print(f'C3 deposited abundance {mat.shape}, ja={len(ctrl)} sni={len(case)}')
    report['C3'] = run_cohort(mat, genes, ctrl, case,
                              'C3 GSE205906 SNI', normalised_deposited=True)

    with open(os.path.join(OUT, 'normalisation_sensitivity.json'), 'w', encoding='utf-8') as fh:
        json.dump(report, fh, indent=2, ensure_ascii=False)

    print()
    print('=' * 88)
    print('Survivors of Benjamini-Hochberg correction by pre-processing regime')
    print('=' * 88)
    for c in ('C1', 'C2', 'C3'):
        print(f'\n{c}: {report[c]["_label"]}')
        if 'size_factors' in report[c]:
            print('  library totals      :', report[c]['library_totals'])
            print('  median-ratio factors:', report[c]['size_factors']['median_ratio'])
        for k in sorted(k for k in report[c]
                        if isinstance(report[c][k], dict) and 'n_genes' in report[c][k]):
            e = report[c][k]
            a = e['Atf3']
            at = (f"Atf3 {a['fold']:>7.2f}-fold p={a['p']:.4f} q={a['q']:.3f}"
                  if a else 'Atf3 n/a')
            print(f'  {k:32s} genes={e["n_genes"]:6d} nominal={e["nominal_p05"]:5d} '
                  f'q<0.05={e["q05_survivors"]:3d}  median log2FC={e["median_log2FC"]:+.3f}  {at}')


if __name__ == '__main__':
    main()
