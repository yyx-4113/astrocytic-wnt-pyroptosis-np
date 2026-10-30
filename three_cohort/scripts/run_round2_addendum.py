#!/usr/bin/env python3
"""Round-2 review addendum.

Computes what the second review panel showed to be missing or wrong:

  A. Hartung-Knapp-Sidik-Jonkman intervals for every pooled Wnt gene
     (the headline 31/36 rests on a normal-quantile interval at k = 2-3).
  B. The pyroptosis execution-subset rank-clustering count, recomputed from
     the actual ranking (the manuscript says 8 of 9; the panel measured 6).
  C. A power / minimum-detectable-effect calculation read off these data.
  D. Reconciliation of three quantities whose provenance was challenged:
     Atf3 and Gfap in C3, and the size of the priming subset.

Everything is derived from the published result tables; no number here is
typed in by hand.  Outputs: round2_results.json, meta_focus_hksj.csv.
"""
import csv
import json
import math
import os
import sys

from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_three_cohort as R  # noqa: E402  (gene panels and DL implementation)

SEARCH = [".", "results", os.path.join(HERE, ".."),
          os.path.join(HERE, "..", "results"), os.path.join(HERE, "results")]


def P(fn):
    for d in SEARCH:
        c = os.path.join(d, fn)
        if os.path.exists(c):
            return c
    return fn


WNT = {s.upper() for s in R.WNT}
OUT = R.RESULTS  # <pkg>/results in the repository layout; script dir when flat
Z = 1.959963985


def load(fn):
    """symbol -> (log2FC, se).  SE is back-derived from the t statistic,
    which is the least rounded quantity in the published tables."""
    d = {}
    with open(P(fn), encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if not r["log2FC"] or not r["t_stat"]:
                continue
            lfc, t = float(r["log2FC"]), float(r["t_stat"])
            if abs(t) > 1e-9:
                d[r["symbol"].upper()] = (lfc, abs(lfc / t), float(r["pvalue"]),
                                          float(r["padj"]), float(r["mean_expr"]))
    return d


C = {"C1": load("cohort_C1_GSE217932_CCI_gene_level.csv"),
     "C2": load("cohort_C2_GSE256472_SNI_gene_level.csv"),
     "C3": load("cohort_C3_GSE205906_SNI_gene_level.csv")}
N_PER_GROUP = {"C1": 4, "C2": 3, "C3": 3}


def dl_hksj(y, se):
    """DerSimonian-Laird point estimate with both the normal-quantile and
    the Hartung-Knapp-Sidik-Jonkman interval."""
    k = len(y)
    v = [s * s for s in se]
    w = [1.0 / x for x in v]
    yf = sum(wi * yi for wi, yi in zip(w, y)) / sum(w)
    Q = sum(wi * (yi - yf) ** 2 for wi, yi in zip(w, y))
    Cw = sum(w) - sum(wi * wi for wi in w) / sum(w)
    tau2 = max(0.0, (Q - (k - 1)) / Cw) if Cw > 0 else 0.0
    wr = [1.0 / (vi + tau2) for vi in v]
    sw = sum(wr)
    ybar = sum(wi * yi for wi, yi in zip(wr, y)) / sw
    se_re = math.sqrt(1.0 / sw)
    # normal-quantile interval (what the manuscript currently reports)
    lo_z, hi_z = ybar - Z * se_re, ybar + Z * se_re
    # HKSJ: rescale by sqrt(q/(k-1)) and use the t_{k-1} quantile
    q = sum(wi * (yi - ybar) ** 2 for wi, yi in zip(wr, y)) / (k - 1) if k > 1 else 0.0
    se_h = math.sqrt(max(q, 0.0)) * se_re
    tcrit = stats.t.ppf(0.975, k - 1) if k > 1 else float("inf")
    lo_h, hi_h = ybar - tcrit * se_h, ybar + tcrit * se_h
    I2 = max(0.0, (Q - (k - 1)) / Q) * 100.0 if Q > 0 else 0.0
    return dict(k=k, re=ybar, se_re=se_re, tau2=tau2, I2=I2,
                lo_z=lo_z, hi_z=hi_z, lo_hksj=lo_h, hi_hksj=hi_h,
                se_hksj=se_h, tcrit=tcrit,
                p_z=2 * stats.norm.sf(abs(ybar / se_re)) if se_re > 0 else 1.0,
                p_hksj=(2 * stats.t.sf(abs(ybar / se_h), k - 1)
                        if se_h > 0 and k > 1 else 1.0))


def equivalence(lo, hi, bound=0.5):
    return bool(float(lo) > -bound and float(hi) < bound)


# ---------------------------------------------------------------- A. HKSJ
def pooled(panel_genes, cohorts, label):
    rows = []
    for g in sorted(panel_genes):
        present = [c for c in cohorts if g in C[c]]
        if len(present) < 2:
            continue
        y = [C[c][g][0] for c in present]
        se = [C[c][g][1] for c in present]
        m = dl_hksj(y, se)
        m["gene"] = g
        m["cohorts"] = "+".join(present)
        m["eq_z"] = equivalence(m["lo_z"], m["hi_z"])
        m["eq_hksj"] = equivalence(m["lo_hksj"], m["hi_hksj"])
        rows.append(m)
    nz = sum(1 for r in rows if r["eq_z"])
    nh = sum(1 for r in rows if r["eq_hksj"])
    return dict(label=label, n=len(rows), eq_z=nz, eq_hksj=nh, rows=rows)


three = pooled(WNT, ["C1", "C2", "C3"], "C1+C2+C3")
two = pooled(WNT, ["C1", "C2"], "C1+C2 (C3 excluded)")

with open(os.path.join(OUT, "meta_focus_hksj.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["analysis", "gene", "k", "cohorts", "pooled_RE_log2FC",
                "RE_lo_normal", "RE_hi_normal", "RE_lo_HKSJ", "RE_hi_HKSJ",
                "tau2", "I2_pct", "equivalent_normal", "equivalent_HKSJ"])
    for block in (three, two):
        for r in block["rows"]:
            w.writerow([block["label"], r["gene"], r["k"], r["cohorts"],
                        f"{r['re']:.4f}", f"{r['lo_z']:.4f}", f"{r['hi_z']:.4f}",
                        f"{r['lo_hksj']:.4f}", f"{r['hi_hksj']:.4f}",
                        f"{r['tau2']:.4f}", f"{r['I2']:.1f}",
                        int(r["eq_z"]), int(r["eq_hksj"])])

# ------------------------------------------------- B. execution rank counts
exec_genes = ["NLRP3", "NLRC4", "AIM2", "PYCARD", "NEK7", "CASP1",
              "CASP4", "GSDMD", "GSDME"]


def rank_by_t(cohort, expr_quantile=0.25):
    rows = []
    for g, (lfc, se, p, padj, mean_expr) in C[cohort].items():
        if se <= 0:
            continue
        rows.append((g, lfc / se, mean_expr, lfc))
    thr = sorted(r[2] for r in rows)[int(expr_quantile * len(rows))]
    kept = [(g, t, lfc) for g, t, me, lfc in rows if me >= thr]
    for i in range(len(kept)):
        g, t, lfc = kept[i]
        kept[i] = (g, max(-50.0, min(50.0, t)), lfc)
    kept.sort(key=lambda x: -x[1])
    return kept, thr


def cluster_count(ranked, genes, frac=0.06):
    N = len(ranked)
    K = int(round(frac * N))
    pos = {g: i + 1 for i, (g, _, _) in enumerate(ranked)}
    obs = sum(1 for g in genes if g in pos and pos[g] <= K)
    exp = K / N * len(genes)
    p = stats.hypergeom.sf(obs - 1, N, K, len(genes))
    return dict(N=N, K=K, observed=obs, n_genes=len(genes), expected=exp,
                p=p, ranks={g: pos.get(g) for g in genes},
                pct={g: (100.0 * pos[g] / N if g in pos else None) for g in genes})


ranked_c1, thr_c1 = rank_by_t("C1", 0.25)
ranked_c1_nf, thr_c1_nf = rank_by_t("C1", 0.0)
exec_primary = cluster_count(ranked_c1, exec_genes)
exec_nofilter = cluster_count(ranked_c1_nf, exec_genes)
# pre-specified alternatives, to show the cut-point is not doing the work
exec_top5 = cluster_count(ranked_c1, exec_genes, 0.05)
exec_top10 = cluster_count(ranked_c1, exec_genes, 0.10)

# ---------------------------------------------------------- C. power / MDE
def mde(se_median, n, alpha=0.05, power=0.80):
    df = 2 * n - 2
    t_alpha = stats.t.ppf(1 - alpha / 2, df)
    t_beta = stats.t.ppf(power, df)
    return (t_alpha + t_beta) * se_median


power_rows = []
for c in ("C1", "C2", "C3"):
    ses = [v[1] for v in C[c].values() if v[1] > 0]
    ses.sort()
    se_med = ses[len(ses) // 2]
    n = N_PER_GROUP[c]
    sd = se_med * math.sqrt(n / 2.0)          # se = SD*sqrt(2/n)
    ntot = len(C[c])
    # nominal alpha = 0.05
    m_nom = mde(se_med, n)
    # Benjamini-Hochberg: smallest p that can survive q<0.05 is 0.05/ntot
    alpha_bh = 0.05 / ntot
    df = 2 * n - 2
    t_bh = stats.t.ppf(1 - alpha_bh / 2, df)
    m_bh = (t_bh + stats.t.ppf(0.80, df)) * se_med
    power_rows.append(dict(cohort=c, n_per_group=n, genes=ntot,
                           median_within_group_SD=sd, median_SE=se_med,
                           MDE_log2_nominal=m_nom,
                           MDE_fold_nominal=2 ** m_nom,
                           MDE_log2_BHsurviving=m_bh,
                           MDE_fold_BHsurviving=2 ** m_bh))

# --------------------------------------------------- D. provenance checks
def gene_row(cohort, gene):
    v = C[cohort].get(gene)
    if v is None:
        return None
    return dict(log2FC=v[0], fold=2 ** v[0], p=v[2], q=v[3], mean_expr=v[4])


prov = {f"{g}_C3": gene_row("C3", g) for g in ("ATF3", "GFAP", "AIF1", "C1QA", "NLRP3")}
prov["Atf3_C1"] = gene_row("C1", "ATF3")
prov["Atf3_C2"] = gene_row("C2", "ATF3")

# priming subset: how many of the 13 defined genes are actually measurable?
import run_gsea as G  # noqa: E402
panels = {n: gs for n, gs, _src in G.SETS}
priming = None
for name, genes in panels.items():
    if "PRIM" in name.upper():
        priming = dict(name=name, defined=len(genes),
                       genes=list(genes),
                       measurable=sum(1 for g in genes if g.upper() in C["C1"]))
exec_panel = None
for name, genes in panels.items():
    if "EXEC" in name.upper() and exec_panel is None:
        exec_panel = dict(name=name, defined=len(genes), genes=list(genes),
                          measurable=sum(1 for g in genes if g.upper() in C["C1"]))

out = dict(
    hksj=dict(three_cohorts={k: three[k] for k in ("label", "n", "eq_z", "eq_hksj")},
              two_cohorts={k: two[k] for k in ("label", "n", "eq_z", "eq_hksj")},
              genes=[{k: r[k] for k in ("gene", "k", "cohorts", "re", "lo_z", "hi_z",
                                        "lo_hksj", "hi_hksj", "eq_z", "eq_hksj",
                                        "tau2", "I2", "p_z", "p_hksj")}
                     for r in three["rows"]],
              genes_two=[{k: r[k] for k in ("gene", "k", "re", "lo_z", "hi_z",
                                            "lo_hksj", "hi_hksj", "eq_z", "eq_hksj")}
                         for r in two["rows"]]),
    execution_cluster=dict(primary=exec_primary, no_filter=exec_nofilter,
                           top5=exec_top5, top10=exec_top10),
    power=power_rows,
    provenance=prov,
    subsets=dict(priming=priming, execution=exec_panel),
)

with open(os.path.join(OUT, "round2_results.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, indent=2, ensure_ascii=False)

print("=" * 72)
print("A. Pooled Wnt intervals: normal-quantile vs Hartung-Knapp-Sidik-Jonkman")
print("=" * 72)
for b in (three, two):
    print(f"  {b['label']:24s} poolable={b['n']:2d}  "
          f"equivalent: normal {b['eq_z']}/{b['n']}   HKSJ {b['eq_hksj']}/{b['n']}")
print()
print("B. Pyroptosis execution subset in the C1 ranking")
for lbl, d in [("primary (6%)", exec_primary), ("no filter (6%)", exec_nofilter),
               ("top 5%", exec_top5), ("top 10%", exec_top10)]:
    print(f"  {lbl:16s} {d['observed']}/{d['n_genes']} in top {d['K']:5d} of {d['N']:5d} "
          f"(expected {d['expected']:.2f})  p={d['p']:.3e}")
print("  ranks:", {g: exec_primary["ranks"][g] for g in exec_genes})
print("  pct  :", {g: (round(exec_primary["pct"][g], 1) if exec_primary["pct"][g] else None)
                   for g in exec_genes})
print()
print("C. Power / minimum detectable effect (median gene, 80% power)")
for r in power_rows:
    print(f"  {r['cohort']} n={r['n_per_group']}/group  median within-group SD={r['median_within_group_SD']:.3f} log2")
    print(f"      MDE nominal p<0.05      : {r['MDE_log2_nominal']:.2f} log2 = {r['MDE_fold_nominal']:.2f}-fold")
    print(f"      MDE surviving BH q<0.05 : {r['MDE_log2_BHsurviving']:.2f} log2 = {r['MDE_fold_BHsurviving']:.1f}-fold")
print()
print("D. Provenance reconciliation")
for k, v in prov.items():
    print(f"  {k:12s} {v}")
print()
