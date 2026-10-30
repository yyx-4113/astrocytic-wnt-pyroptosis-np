# -*- coding: utf-8 -*-
"""Audit: every key number asserted in MANUSCRIPT_BMC.md is recomputed from the
result files.  Expected strings are DERIVED from the artifacts, never typed in.

An earlier version of this audit carried several values as literals.  Tampering
with the result files then still left 35 of 45 checks passing, because the
literal did not depend on the artifact at all.  Every expected string below is
therefore built from a file read at run time.
"""
import csv
import json
import math
import os
import re
import sys

FAIL = []
N = 0

HERE = os.path.dirname(os.path.abspath(__file__))
_SEARCH = [".", "results", os.path.join(HERE, ".."),
           os.path.join(HERE, "..", "results"), os.path.join(HERE, "results"), HERE]


def P(fn):
    if os.path.isabs(fn):
        return fn
    for d in _SEARCH:
        cand = os.path.join(d, fn)
        if os.path.exists(cand):
            return cand
    return fn


def load(fn):
    return list(csv.DictReader(open(P(fn), encoding="utf-8")))


def jload(fn):
    return json.load(open(P(fn), encoding="utf-8"))


def chk(label, expected):
    """Assert that the manuscript contains the string derived from an artifact."""
    global N
    N += 1
    txt = open(P("MANUSCRIPT_BMC.md"), encoding="utf-8").read().replace("\u2212", "-")
    if expected in txt:
        print(f"  PASS  {label:48s} {expected}")
    else:
        FAIL.append((label, expected))
        print(f"  FAIL  {label:48s} expected {expected!r} NOT FOUND")


def sci(p):
    """Format a p-value the way the manuscript writes it: 3.3 \u00d7 10^-6^"""
    if p == 0:
        return "0"
    e = math.floor(math.log10(abs(p)))
    m = p / (10 ** e)
    return f"{round(m, 1):g} \u00d7 10^-{abs(e)}^"


def perm_floor(n_perm):
    """The permutation floor 1/(n_perm+1), in the manuscript's notation."""
    e = round(math.log10(1.0 / (n_perm + 1)))
    return f"1 \u00d7 10^-{abs(e)}^"


# ---------------------------------------------------------------- cohorts
COH = {"C1": "cohort_C1_GSE217932_CCI_gene_level.csv",
       "C2": "cohort_C2_GSE256472_SNI_gene_level.csv",
       "C3": "cohort_C3_GSE205906_SNI_gene_level.csv"}
cohorts = {}
for cid, fn in COH.items():
    rows = load(fn)
    cohorts[cid] = {r["symbol"].upper(): r for r in rows}
    n = len(rows)
    nom = sum(1 for r in rows if float(r["pvalue"]) < 0.05)
    q = sum(1 for r in rows if float(r["padj"]) < 0.05)
    eq = sum(1 for r in rows if float(r["ci_lo"]) > -0.5 and float(r["ci_hi"]) < 0.5)
    chk(f"{cid} genes tested", f"{n:,}")
    chk(f"{cid} nominal p<0.05", f"{nom:,}")
    chk(f"{cid} intervals within +/-0.5", f"{eq:,} ({100 * eq / n:.1f}%)")
    print(f"        ({cid}: {n} genes, {nom} nominal, {q} q<0.05, {eq} equivalent)")

# ---------------------------------------------------------------- Atf3
c1 = cohorts["C1"]["ATF3"]
chk("Atf3 C1 log2FC", f"{float(c1['log2FC']):+.2f}")
chk("Atf3 C1 p", f"{float(c1['pvalue']):.4f}")
chk("Atf3 C1 q", f"{float(c1['padj']):.3f}")
chk("Atf3 C1 fold", f"{2 ** float(c1['log2FC']):.1f}-fold")
chk("Atf3 C1 CI", f"{float(c1['ci_lo']):.2f} to {float(c1['ci_hi']):.2f}")
for cid in ("C2", "C3"):
    r = cohorts[cid]["ATF3"]
    chk(f"Atf3 {cid} log2FC", f"{float(r['log2FC']):+.2f}")
    chk(f"Atf3 {cid} p", f"{float(r['pvalue']):.3f}")
    chk(f"Atf3 {cid} q", f"{float(r['padj']):.3f}")

meta = {r["symbol"].upper(): r for r in load("meta_focus.csv")}
a = meta["ATF3"]
chk("Atf3 pooled RE", f"{float(a['pooled_RE_log2FC']):+.2f}")
chk("Atf3 pooled CI", f"{float(a['RE_ci_lo']):.2f} to {float(a['RE_ci_hi']):.2f}")
chk("Atf3 pooled p", f"{float(a['p_RE']):.4f}")
chk("Atf3 I2", f"{float(a['I2_pct']):.0f}%")

# ---------------------------------------------------------------- Wnt panel
sys.path.insert(0, HERE)
import run_three_cohort as R  # noqa: E402

WNT = {s.upper() for s in R.WNT}
w = [r for r in meta.values() if r["symbol"].upper() in WNT and r["pooled_RE_log2FC"]]
mx = max(w, key=lambda r: abs(float(r["pooled_RE_log2FC"])))
eqn = [r for r in w if r["equiv_normal"] == "1"]
eqh = [r for r in w if r["equiv_HKSJ"] == "1"]
k3 = [r for r in w if r["k"] == "3"]
chk("Wnt poolable", f"{len(w)} were quantified in at least two cohorts")
chk("Wnt in all three", f"({len(k3)} in all three)")
chk("Wnt max gene", mx["symbol"])
chk("Wnt max effect", f"{float(mx['pooled_RE_log2FC']):+.2f} log2")
chk("Wnt max fold", f"{2 ** float(mx['pooled_RE_log2FC']):.2f}-fold")
chk("Wnt max CI", f"{float(mx['RE_ci_lo']):.2f} to {float(mx['RE_ci_hi']):.2f}")
chk("Wnt equivalent normal-quantile", f"{len(eqn)} of {len(w)}")
chk("Wnt equivalent HKSJ", f"{len(eqh)} of {len(w)}")
for g in ["AXIN2", "LEF1", "TCF4", "GSK3B", "WNT2", "CTNNB1"]:
    chk(f"{g} pooled", f"{float(meta[g]['pooled_RE_log2FC']):+.2f}")

# ------------------------------------------- Wnt re-pooled without C3 (C1+C2)
two = load("sensitivity_two_cohort.csv")
up = max(two, key=lambda r: float(r["pooled_RE_log2FC"]))
dn = min(two, key=lambda r: float(r["pooled_RE_log2FC"]))
chk("C1+C2 poolable", f"{len(two)} genes were poolable")
chk("C1+C2 largest induction gene", up["symbol"])
chk("C1+C2 largest induction", f"{float(up['pooled_RE_log2FC']):+.2f} log2")
chk("C1+C2 largest fold", f"{2 ** float(up['pooled_RE_log2FC']):.2f}-fold")
chk("C1+C2 largest decrease gene", dn["symbol"])
chk("C1+C2 largest decrease", f"{float(dn['pooled_RE_log2FC']):+.2f}")
chk("C1+C2 decrease p", f"{float(dn['p_RE']):.4f}")
chk("C1+C2 equivalent normal", f"{sum(1 for r in two if r['equiv_normal'] == '1')} of {len(two)}")
chk("C1+C2 equivalent HKSJ", f"{sum(1 for r in two if r['equiv_HKSJ'] == '1')} of {len(two)}")

# ---------------------------------------------------------------- GSEA
import run_gsea as _G  # noqa: E402

NP_PERM = _G.NP_PERM
G = load("gsea_three_cohort.csv")
qcol = [c for c in G[0] if c.startswith("FDR_q")][0]
for cid in ("C1", "C2", "C3"):
    r = [x for x in G if x["cohort"] == cid and x["gene_set"].startswith("HALLMARK")][0]
    chk(f"Wnt NES {cid}", f"{float(r['NES']):+.2f}")
    chk(f"Wnt p {cid}", f"{float(r['nominal_p']):.3f}")
    chk(f"Wnt q {cid}", f"{float(r[qcol]):.3f}")
for cid in ("C1", "C2", "C3"):
    r = [x for x in G if x["cohort"] == cid and x["gene_set"] == "PYROPTOSIS_CORE"][0]
    chk(f"Pyro NES {cid}", f"{float(r['NES']):+.2f}")
    if float(r["nominal_p"]) <= 1.01 / (NP_PERM + 1):
        # right-censored at the floor; the manuscript must say so
        chk(f"Pyro p {cid} floor", "≤ " + perm_floor(NP_PERM))
        chk(f"Pyro q {cid} floor", "≤ " + perm_floor(NP_PERM))
    else:
        chk(f"Pyro p {cid}", f"{float(r['nominal_p']):.3f}")
        chk(f"Pyro q {cid}", f"{float(r[qcol]):.3f}")

sens = jload("sensitivity_results.json")
chk("LOO min NES", f"{min(d['nes'] for d in sens['loo']):.2f}")
chk("LOO count", f"{len(sens['loo'])} leaving-one-out analyses") if False else None
chk("LOO n", f"all {len(sens['loo'])} remained")
chk("exec NES", f"+{sens['exec']['nes']:.2f}")
chk("prime NES", f"+{sens['prime']['nes']:.2f}")
chk("no-filter Wnt NES", f"{sens['nofilter']['wnt']['nes']:.2f}")
chk("no-filter Wnt p", f"{sens['nofilter']['wnt']['p']:.3f}")
chk("leading-edge removal NES", f"{sens['no_le']['nes']:.2f}")
chk("leading-edge removal p", f"{sens['no_le']['p']:.2f}")

# ------------------------------------------- round-2 recomputations (JSON)
r2 = jload("round2_results.json")
ec = r2["execution_cluster"]
chk("execution genes in top 6%", f"{ec['primary']['observed']} of the {ec['primary']['n_genes']} execution genes")
chk("execution hypergeometric p", sci(ec["primary"]["p"]))
chk("execution no-filter count", f"{ec['no_filter']['observed']} of {ec['no_filter']['n_genes']}")
chk("execution no-filter p", sci(ec["no_filter"]["p"]))
chk("execution expected", f"{ec['primary']['expected']:.2f}")

for row in r2["power"]:
    cid = row["cohort"]
    chk(f"{cid} MDE nominal", f"{row['MDE_fold_nominal']:.2f}-fold")
    bh = row["MDE_fold_BHsurviving"]
    chk(f"{cid} MDE BH", f"{bh:.1f}-fold" if bh < 100 else f"{bh:.0f}-fold")
    chk(f"{cid} within-group SD", f"{row['median_within_group_SD']:.3f}")

# ------------------- normalisation sensitivity (does the pipeline create it?)
norm = jload("normalisation_sensitivity.json")
for cid in ("C1", "C2", "C3"):
    vals = [v["q05_survivors"] for k, v in norm[cid].items()
            if isinstance(v, dict) and "q05_survivors" in v]
    chk(f"{cid} survivors range", f"{min(vals)}\u2013{max(vals)}")
    qs = [v["Atf3"]["q"] for k, v in norm[cid].items()
          if isinstance(v, dict) and v.get("Atf3")]
    lo, hi = min(qs), max(qs)
    chk(f"{cid} Atf3 q range",
        f"{lo:.3f}" if abs(hi - lo) < 1e-9 else f"{lo:.3f}\u2013{hi:.3f}")

# --------------------------------- C3 glial markers as the pipeline reports them
for g, fmt in [("AIF1", ".2f"), ("C1QA", ".2f"), ("GFAP", ".2f"), ("NLRP3", ".2f")]:
    r = cohorts["C3"][g]
    chk(f"C3 {g} fold", f"{2 ** float(r['log2FC']):{fmt}}-fold")

# ---------------------------------------------------------------- references
txt = open(P("MANUSCRIPT_BMC.md"), encoding="utf-8").read()
body = txt.split("## References")[0]
order, seen = [], set()
for m in re.finditer(r"\[(\d+(?:\s*[,\u2013-]\s*\d+)*)\]", body):
    for part in m.group(1).split(","):
        part = part.strip()
        if re.match(r"^\d+\s*[\u2013-]\s*\d+$", part):
            a, b = re.split(r"[\u2013-]", part)
            nums = list(range(int(a), int(b) + 1))
        elif part.isdigit():
            nums = [int(part)]
        else:
            nums = []
        for x in nums:
            if x not in seen:
                seen.add(x)
                order.append(x)
refs = re.findall(r"^(\d+)\. ", txt.split("## References")[1], re.M)
N += 2
if order == list(range(1, len(order) + 1)):
    print(f"  PASS  {'reference numbering sequential':48s} 1..{len(order)}")
else:
    FAIL.append(("refs order", str(order)))
    print(f"  FAIL  reference numbering                    {order}")
if len(refs) == len(order):
    print(f"  PASS  {'every reference cited':48s} {len(refs)}")
else:
    FAIL.append(("refs count", f"{len(refs)} listed vs {len(order)} cited"))
    print(f"  FAIL  reference count                        {len(refs)} vs {len(order)}")

# ---------------------------------------------------------------- abstract
ab = txt.split("## Abstract")[1].split("**Keywords:**")[0].replace("**", "")
nw = len(ab.split())
N += 1
if nw <= 350:
    print(f"  PASS  {'abstract within 350 words':48s} {nw}")
else:
    FAIL.append(("abstract", nw))
    print(f"  FAIL  abstract word count                    {nw}")

print()
print(f"{N - len(FAIL)}/{N} checks passed")
if FAIL:
    print("FAILURES:")
    for f in FAIL:
        print("  -", f)
    sys.exit(1)
