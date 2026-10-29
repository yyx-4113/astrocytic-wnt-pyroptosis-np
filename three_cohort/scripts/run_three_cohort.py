# -*- coding: utf-8 -*-
"""
Three-cohort validation of Wnt / pyroptosis transcriptional responses in
neuropathic-pain spinal cord.

Cohorts (one contrast each, intervention arms excluded):
  C1  GSE217932  Rattus norvegicus  CCI  vs Sham   n=4 vs 4   raw counts
  C2  GSE256472  Rattus norvegicus  SNI  vs Control n=3 vs 3   raw counts
  C3  GSE205906  Rattus norvegicus  SNI  vs Sham(JA) n=3 vs 3  FPKM (processed)

Unified pipeline (identical for all three cohorts):
  1. log2(x + 1)
  2. gene-level collapse by symbol: among transcripts of one gene choose the
     transcript with the highest mean log2 expression (deterministic, unbiased
     w.r.t. the tested contrast)
  3. Welch two-sample t-test on log2 values (unequal variance)
  4. Benjamini-Hochberg FDR across all tested genes
  5. 95% CI of log2FC from the Welch standard error and Welch-Satterthwaite df

Outputs (RESCUE_2026-10-29/):
  cohort_C1_gene_level.csv, cohort_C2_gene_level.csv, cohort_C3_gene_level.csv
  three_cohort_focus.csv      (Wnt, pyroptosis, positive-control genes)
  meta_focus.csv              (DerSimonian-Laird random effects + fixed effect)
  run_summary.txt
"""
import os, sys, csv, gzip, math, pickle
import numpy as np
from scipy import stats

DATA = r'D:/省自然/生信分析/data'
OUT = r'D:/省自然/生信分析/submission/RESCUE_2026-10-29'
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------- gene sets
WNT = ['Ctnnb1','Tcf4','Tcf3','Tcf7','Tcf7l1','Tcf7l2','Lef1','Apc','Axin1','Axin2',
       'Gsk3a','Gsk3b','Ccnd1','Ccnd2','Ccnd3','Myc','Mycn','Jun','Fosl1','Cd44',
       'Mmp7','Birc5','Sox9','Wnt2','Wnt3','Wnt3a','Wnt5a','Wnt10b','Fzd1','Fzd4',
       'Fzd7','Lrp5','Lrp6','Dkk1','Dvl1','Dvl2','Dvl3']

PYRO = ['Nlrp3','Nlrc4','Aim2','Pycard','Casp1','Casp4','Casp8','Casp3','Gsdmd','Gsdme',
        'Il1b','Il18','Il1a','Hmgb1','Tnf','Il6','Nek7','Txnip','Tlr4','Nfkb1','Rela','Gzmb']

# positive controls: injury-response genes expected to move in any working
# neuropathic-pain model (Atf3 = regenerative/injury program; Cacna2d1 = alpha2delta-1,
# the gabapentinoid target and a canonical neuropathic-pain marker; Ccl7 = microglial
# chemokine; Limk1, Eln = largest observed effects in C1)
POSCTRL = ['Atf3','Cacna2d1','Ccl7','Limk1','Eln','Gfap','Aif1','Iba1','Cx3cr1','Tmem119',
           'Gja1','Sox2','Bdnf','Nefl','Snap25','Syp','Mbp','Plp1']

FOCUS = []
for g in WNT + PYRO + POSCTRL:
    if g not in FOCUS:
        FOCUS.append(g)

# ---------------------------------------------------------------- helpers
def bh(pvals):
    m = len(pvals)
    if m == 0:
        return []
    order = sorted(range(m), key=lambda i: pvals[i])
    adj = [0.0] * m
    prev = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        val = min(prev, pvals[i] * m / rank)
        adj[i] = val
        prev = val
    return adj


def welch_ci(a, b):
    """Welch t, p, se, df, 95% CI of (mean_b - mean_a)."""
    na, nb = len(a), len(b)
    ma, mb = float(np.mean(a)), float(np.mean(b))
    va, vb = float(np.var(a, ddof=1)), float(np.var(b, ddof=1))
    se = math.sqrt(va / na + vb / nb)
    lfc = mb - ma
    if se == 0 or math.isnan(se):
        return lfc, float('nan'), float('nan'), float('nan'), float('nan'), float('nan')
    t = lfc / se
    num = (va / na + vb / nb) ** 2
    den = (va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1)
    df = num / den if den > 0 else float(na + nb - 2)
    p = 2 * stats.t.sf(abs(t), df)
    tc = stats.t.ppf(0.975, df)
    return lfc, p, t, df, lfc - tc * se, lfc + tc * se


def collapse_and_test(records, group_a_idx, group_b_idx, sample_names, label):
    """records: list of (symbol, np.array(log2 values across all samples), mean_expr)

    Gene-level collapse: per symbol keep the transcript with max mean expression.
    Then Welch test between group indices.
    """
    best = {}
    for sym, vals, m in records:
        if sym is None:
            continue
        s = sym.strip()
        if not s or s == '-':
            continue
        k = s.upper()
        if k not in best or m > best[k][2]:
            best[k] = (s, vals, m)
    rows = []
    pv = []
    for k, (sym, vals, m) in best.items():
        a = [vals[i] for i in group_a_idx]
        b = [vals[i] for i in group_b_idx]
        if len(a) < 2 or len(b) < 2:
            continue
        lfc, p, t, df, lo, hi = welch_ci(a, b)
        if math.isnan(p):
            continue
        rows.append({'symbol': sym, 'key': k, 'lfc': lfc, 'p': p, 't': t, 'df': df,
                     'ci_lo': lo, 'ci_hi': hi, 'mean_expr': m,
                     'mean_a': float(np.mean(a)), 'mean_b': float(np.mean(b))})
        pv.append(p)
    adj = bh(pv)
    for r, q in zip(rows, adj):
        r['padj'] = q
    n_sig = sum(1 for r in rows if r['padj'] < 0.05)
    n_nom = sum(1 for r in rows if r['p'] < 0.05)
    n_equiv = sum(1 for r in rows
                  if r['ci_lo'] is not None and not math.isnan(r['ci_lo'])
                  and r['ci_lo'] > -0.5 and r['ci_hi'] < 0.5)
    print(f'  [{label}] genes={len(rows)} nominal_p<0.05={n_nom} '
          f'padj<0.05={n_sig} equivalence(|CI|<0.5)={n_equiv}')
    path = os.path.join(OUT, f'cohort_{label}_gene_level.csv')
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['symbol','log2FC','pvalue','t_stat','df','padj',
                    'ci_lo','ci_hi','mean_expr','mean_control','mean_case'])
        for r in sorted(rows, key=lambda x: x['p']):
            w.writerow([r['symbol'], f"{r['lfc']:.6f}", f"{r['p']:.6e}",
                        f"{r['t']:.6f}", f"{r['df']:.3f}", f"{r['padj']:.6f}",
                        f"{r['ci_lo']:.6f}", f"{r['ci_hi']:.6f}",
                        f"{r['mean_expr']:.6f}", f"{r['mean_a']:.6f}", f"{r['mean_b']:.6f}"])
    print(f'  [{label}] wrote {path}')
    return {r['key']: r for r in rows}, {'n': len(rows), 'nom': n_nom,
                                          'sig': n_sig, 'equiv': n_equiv}


# ---------------------------------------------------------------- cohort 1
def load_c1():
    import openpyxl
    src = os.path.join(DATA, 'GSE217932_Raw_gene_counts_matrix.xlsx')
    wb = openpyxl.load_workbook(src, read_only=True)
    ws = wb.active
    it = ws.iter_rows(values_only=True)
    header = list(next(it))
    cols = header[1:]
    sham = [i for i, c in enumerate(cols) if 'Sham' in str(c)]
    cci = [i for i, c in enumerate(cols) if str(c).startswith('CCI-') and 'THP' not in str(c)]
    print(f'  [C1] samples={len(cols)} sham={len(sham)} cci={len(cci)} '
          f'(CCI+THP excluded: {len(cols)-len(sham)-len(cci)} wells)')
    gmap = pickle.load(open(os.path.join(DATA, 'rat_gene_map.pkl'), 'rb'))
    recs = []
    n = 0
    for row in it:
        gid = row[0]
        if gid is None:
            continue
        n += 1
        sym = gmap.get(str(gid).split('.')[0])
        vals = np.array([0.0 if v is None else math.log2(float(v) + 1) for v in row[1:]],
                        dtype=float)
        recs.append((sym, vals, float(np.mean(vals))))
    print(f'  [C1] transcripts read={n}, with symbol={sum(1 for r in recs if r[0])}')
    return recs, sham, cci, cols


# ---------------------------------------------------------------- cohort 2
def load_c2():
    src = os.path.join(DATA, 'GSE256472_all.genes.expression.annot.xls.gz')
    print('  [C2] reading counts columns (fpkm columns ignored)')
    recs = []
    with gzip.open(src, 'rt', encoding='utf-8', errors='replace') as fh:
        header = fh.readline().rstrip('\n').split('\t')
        ctrl = [i for i, c in enumerate(header) if c.endswith('_count') and c.startswith('Control')]
        sni = [i for i, c in enumerate(header) if c.endswith('_count') and c.startswith('SNI')]
        si = header.index('Symbol')
        print(f'  [C2] control={len(ctrl)} sni={len(sni)} symbol_col={si}')
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) <= si:
                continue
            vals_all = [0.0] * len(header)
            for i in range(len(header)):
                if i >= len(p):
                    continue
                try:
                    vals_all[i] = math.log2(float(p[i]) + 1)
                except (ValueError, TypeError):
                    vals_all[i] = 0.0
            v = np.array(vals_all, dtype=float)
            recs.append((p[si], v, float(np.mean(v))))
    print(f'  [C2] transcripts read={len(recs)}')
    return recs, ctrl, sni, header


# ---------------------------------------------------------------- cohort 3
def load_c3():
    src = os.path.join(DATA, 'GSE205906_processed_data.xls.gz')
    print('  [C3] reading processed FPKM matrix (JA=sham, SNI=model, PLM=treated, excluded)')
    recs = []
    with gzip.open(src, 'rt', encoding='utf-8', errors='replace') as fh:
        header = fh.readline().rstrip('\n').split('\t')
        ja = [i for i, c in enumerate(header) if c.startswith('JA')]
        sni = [i for i, c in enumerate(header) if c.startswith('SNI')]
        plm = [i for i, c in enumerate(header) if c.startswith('PLM')]
        si = 1  # gene_name
        print(f'  [C3] JA(sham)={len(ja)} SNI={len(sni)} PLM(excluded)={len(plm)}')
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) <= max(max(ja), max(sni)):
                continue
            vals_all = [0.0] * len(header)
            for i in range(2, len(header)):
                if i >= len(p):
                    continue
                try:
                    vals_all[i] = math.log2(float(p[i]) + 1)
                except (ValueError, TypeError):
                    vals_all[i] = 0.0
            v = np.array(vals_all, dtype=float)
            recs.append((p[si], v, float(np.mean(v))))
    print(f'  [C3] transcripts read={len(recs)}')
    return recs, ja, sni, header


# ---------------------------------------------------------------- meta
def dersimonian_laird(effects, ses):
    """Random-effects (DL) + fixed-effect pooling. Returns dict."""
    k = len(effects)
    if k == 0:
        return None
    w = [1.0 / (s * s) for s in ses if s and s > 0]
    eff = [e for e, s in zip(effects, ses) if s and s > 0]
    if len(eff) < 1:
        return None
    k = len(eff)
    wsum = sum(w)
    fe = sum(wi * ei for wi, ei in zip(w, eff)) / wsum
    Q = sum(wi * (ei - fe) ** 2 for wi, ei in zip(w, eff))
    df = k - 1
    C = wsum - sum(wi * wi for wi in w) / wsum
    tau2 = max(0.0, (Q - df) / C) if C > 0 and df > 0 else 0.0
    I2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 and df > 0 else 0.0
    wr = [1.0 / (s * s + tau2) for s in ses if s and s > 0]
    wrs = sum(wr)
    re = sum(wi * ei for wi, ei in zip(wr, eff)) / wrs
    se_re = math.sqrt(1.0 / wrs)
    z = re / se_re if se_re > 0 else float('nan')
    p = 2 * stats.norm.sf(abs(z)) if not math.isnan(z) else float('nan')
    se_fe = math.sqrt(1.0 / wsum)
    return {'k': k, 'fe': fe, 'fe_se': se_fe, 're': re, 're_se': se_re,
            're_lo': re - 1.96 * se_re, 're_hi': re + 1.96 * se_re,
            'p': p, 'Q': Q, 'I2': I2, 'tau2': tau2}


# ---------------------------------------------------------------- main
def main():
    log = []

    def say(s):
        print(s, flush=True)
        log.append(s)

    say('=' * 78)
    say('THREE-COHORT VALIDATION — Wnt / pyroptosis in neuropathic-pain spinal cord')
    say('=' * 78)

    cohorts = {}

    say('\n--- C1 GSE217932 (CCI vs Sham, counts) ---')
    recs, ga, gb, names = load_c1()
    t1, s1 = collapse_and_test(recs, ga, gb, names, 'C1_GSE217932_CCI')
    cohorts['C1'] = (t1, s1, 'CCI vs Sham', ga, gb)
    del recs

    say('\n--- C2 GSE256472 (SNI vs Control, counts) ---')
    recs, ga, gb, names = load_c2()
    t2, s2 = collapse_and_test(recs, ga, gb, names, 'C2_GSE256472_SNI')
    cohorts['C2'] = (t2, s2, 'SNI vs Control', ga, gb)
    del recs

    say('\n--- C3 GSE205906 (SNI vs Sham, FPKM) ---')
    recs, ga, gb, names = load_c3()
    t3, s3 = collapse_and_test(recs, ga, gb, names, 'C3_GSE205906_SNI')
    cohorts['C3'] = (t3, s3, 'SNI vs Sham', ga, gb)
    del recs

    # ---------------- per-gene-set summary
    say('\n' + '=' * 78)
    say('GENE-SET LEVEL SUMMARY (gene-level, Welch + BH, equivalence bound +/-0.5)')
    say('=' * 78)
    header = f'{"set":14s} {"cohort":28s} {"tested":>7s} {"nom p<.05":>10s} {"padj<.05":>10s} {"equiv":>7s}'
    say(header)
    for setname, genes in [('Wnt', WNT), ('Pyroptosis', PYRO), ('PositiveCtrl', POSCTRL)]:
        keys = set(g.upper() for g in genes)
        for cid, (tab, st, contrast, _, _) in cohorts.items():
            sub = [tab[k] for k in tab if k in keys]
            nom = sum(1 for r in sub if r['p'] < 0.05)
            sig = sum(1 for r in sub if r['padj'] < 0.05)
            eq = sum(1 for r in sub if r['ci_lo'] > -0.5 and r['ci_hi'] < 0.5)
            say(f'{setname:14s} {cid+" "+contrast:28s} {len(sub):7d} {nom:10d} {sig:10d} {eq:7d}')

    # ---------------- meta per focus gene
    say('\n' + '=' * 78)
    say('META-ANALYSIS (DerSimonian-Laird random effects over 3 cohorts)')
    say('=' * 78)
    rows_out = []
    for g in FOCUS:
        k = g.upper()
        eff, ses, present = [], [], []
        for cid in ['C1', 'C2', 'C3']:
            tab = cohorts[cid][0]
            r = tab.get(k)
            if r is not None:
                t = r['t']
                if t is None or t == 0 or math.isnan(t):
                    continue
                se = abs(r['lfc'] / t)
                if se > 0 and not math.isnan(se):
                    eff.append(r['lfc']); ses.append(se); present.append(cid)
        if len(eff) < 2:
            rows_out.append({'symbol': g, 'k': len(eff), 'note': 'insufficient cohorts'})
            continue
        m = dersimonian_laird(eff, ses)
        if not m:
            continue
        rows_out.append({'symbol': g, 'k': m['k'], 'cohorts': '+'.join(present),
                         'lfc_C1': cohorts['C1'][0].get(k, {}).get('lfc', ''),
                         'lfc_C2': cohorts['C2'][0].get(k, {}).get('lfc', ''),
                         'lfc_C3': cohorts['C3'][0].get(k, {}).get('lfc', ''),
                         're': m['re'], 're_lo': m['re_lo'], 're_hi': m['re_hi'],
                         'p': m['p'], 'I2': m['I2'], 'tau2': m['tau2'], 'fe': m['fe']})
    with open(os.path.join(OUT, 'meta_focus.csv'), 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['symbol','k','cohorts','lfc_C1_CCI','lfc_C2_SNI','lfc_C3_SNI',
                    'pooled_RE_log2FC','RE_ci_lo','RE_ci_hi','p_RE','I2_pct','tau2',
                    'pooled_FE_log2FC'])
        for r in rows_out:
            if 're' not in r:
                w.writerow([r['symbol'], r['k'], r.get('note','')] + [''] * 10)
                continue
            fmt = lambda v: f'{v:.4f}' if isinstance(v, float) else (f'{v:.4f}' if isinstance(v, (int, np.floating)) else v)
            w.writerow([r['symbol'], r['k'], r['cohorts'],
                        fmt(r['lfc_C1']), fmt(r['lfc_C2']), fmt(r['lfc_C3']),
                        fmt(r['re']), fmt(r['re_lo']), fmt(r['re_hi']),
                        f"{r['p']:.4e}", fmt(r['I2']), fmt(r['tau2']), fmt(r['fe'])])
    say(f'  wrote meta_focus.csv ({sum(1 for r in rows_out if "re" in r)} genes pooled)')

    # ---------------- focus table (long format)
    with open(os.path.join(OUT, 'three_cohort_focus.csv'), 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['symbol','cohort','contrast','log2FC','ci_lo','ci_hi','pvalue','padj',
                    'equiv_within_0.5','direction'])
        for g in FOCUS:
            k = g.upper()
            for cid in ['C1', 'C2', 'C3']:
                tab, st, contrast = cohorts[cid][0], cohorts[cid][1], cohorts[cid][2]
                r = tab.get(k)
                if not r:
                    w.writerow([g, cid, contrast, 'NA', '', '', '', '', 'absent_in_cohort', ''])
                    continue
                eq = 'yes' if (r['ci_lo'] > -0.5 and r['ci_hi'] < 0.5) else 'no'
                d = 'up' if r['lfc'] > 0 else 'down'
                w.writerow([g, cid, contrast, f"{r['lfc']:.4f}", f"{r['ci_lo']:.4f}",
                            f"{r['ci_hi']:.4f}", f"{r['p']:.4e}", f"{r['padj']:.4f}",
                            eq, d])
    say('  wrote three_cohort_focus.csv')

    # ---------------- headline numbers
    say('\n' + '=' * 78)
    say('HEADLINE: positive controls (must move) vs Wnt (must stay flat)')
    say('=' * 78)
    for g in ['Atf3', 'Cacna2d1', 'Ccl7', 'Limk1', 'Eln']:
        k = g.upper()
        line = f'  {g:10s}'
        for cid in ['C1', 'C2', 'C3']:
            r = cohorts[cid][0].get(k)
            if r:
                line += f' | {cid} lfc={r["lfc"]:+6.2f} CI=[{r["ci_lo"]:+6.2f},{r["ci_hi"]:+6.2f}]'
            else:
                line += f' | {cid} absent'
        say(line)
    say('')
    for g in ['Ctnnb1', 'Axin2', 'Ccnd1', 'Tcf4', 'Lef1', 'Gsk3b', 'Apc', 'Nlrp3', 'Casp1', 'Gsdmd']:
        k = g.upper()
        line = f'  {g:10s}'
        for cid in ['C1', 'C2', 'C3']:
            r = cohorts[cid][0].get(k)
            if r:
                line += f' | {cid} lfc={r["lfc"]:+6.2f} CI=[{r["ci_lo"]:+6.2f},{r["ci_hi"]:+6.2f}]'
            else:
                line += f' | {cid} absent'
        say(line)

    with open(os.path.join(OUT, 'run_summary.txt'), 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(log) + '\n')
    say('\nwrote run_summary.txt')


if __name__ == '__main__':
    main()
