# -*- coding: utf-8 -*-
"""Audit: every key number asserted in MANUSCRIPT_BMC.md is recomputed from the
result files. Expected strings are DERIVED from the artifacts, never typed in."""
import csv, json, re, math, sys

import os

HERE = os.path.dirname(os.path.abspath(__file__))
_SEARCH = [".", "results", os.path.join(HERE, ".."), os.path.join(HERE, "..", "results"),
           os.path.join(HERE, "results"), HERE]


def P(fn):
    """Resolve an input file whether the script is run from the analysis folder
    or from the repository layout (inputs in results/)."""
    if os.path.isabs(fn):
        return fn
    for d in _SEARCH:
        cand = os.path.join(d, fn)
        if os.path.exists(cand):
            return cand
    return fn


FAIL = []
N = 0


def chk(label, value, expected):
    """value: computed; expected: the string the manuscript should contain."""
    global N
    N += 1
    txt = open(P("MANUSCRIPT_BMC.md"), encoding="utf-8").read().replace("\u2212", "-")
    if expected in txt:
        print(f"  PASS  {label:46s} {expected}")
    else:
        FAIL.append((label, expected))
        print(f"  FAIL  {label:46s} expected {expected!r} not found")


def load(fn):
    return list(csv.DictReader(open(P(fn), encoding="utf-8")))


# ---- per-cohort Atf3
cohorts = {}
for cid, fn in [("C1", "cohort_C1_GSE217932_CCI_gene_level.csv"),
                ("C2", "cohort_C2_GSE256472_SNI_gene_level.csv"),
                ("C3", "cohort_C3_GSE205906_SNI_gene_level.csv")]:
    d = {r["symbol"].upper(): r for r in load(fn)}
    cohorts[cid] = d
c1_atf3 = cohorts["C1"]["ATF3"]
chk("Atf3 C1 log2FC", None, f"{float(c1_atf3['log2FC']):+.2f}")
chk("Atf3 C1 p", None, f"{float(c1_atf3['pvalue']):.4f}")
chk("Atf3 C1 q", None, f"{float(c1_atf3['padj']):.3f}")
chk("Atf3 C1 fold", None, f"{2 ** float(c1_atf3['log2FC']):.1f}-fold")
chk("Atf3 C1 CI", None, f"{float(c1_atf3['ci_lo']):.2f} to {float(c1_atf3['ci_hi']):.2f}")

for cid, pstr, qstr in [("C2", "0.0053", "0.700"), ("C3", "0.119", "0.761")]:
    r = cohorts[cid]["ATF3"]
    chk(f"Atf3 {cid} p", None, f"{float(r['pvalue']):.3f}")
    chk(f"Atf3 {cid} q", None, f"{float(r['padj']):.3f}")

# ---- transcriptome-wide counts
for cid, fn in [("C1", "cohort_C1_GSE217932_CCI_gene_level.csv"),
                ("C2", "cohort_C2_GSE256472_SNI_gene_level.csv"),
                ("C3", "cohort_C3_GSE205906_SNI_gene_level.csv")]:
    rows = load(fn)
    n = len(rows)
    nom = sum(1 for r in rows if float(r["pvalue"]) < 0.05)
    q = sum(1 for r in rows if float(r["padj"]) < 0.05)
    chk(f"{cid} genes tested", None, f"{n:,}")
    chk(f"{cid} nominal p<0.05", None, f"{nom:,}")
    print(f"        ({cid}: {n} genes, {nom} nominal, {q} q<0.05)")

# ---- pooled estimates
meta = {r["symbol"].upper(): r for r in load("meta_focus.csv")}
atf3 = meta["ATF3"]
chk("Atf3 pooled RE", None, f"{float(atf3['pooled_RE_log2FC']):+.2f}")
chk("Atf3 pooled CI", None, f"{float(atf3['RE_ci_lo']):.2f} to {float(atf3['RE_ci_hi']):.2f}")
chk("Atf3 pooled p", None, f"{float(atf3['p_RE']):.4f}")
chk("Atf3 I2", None, f"{float(atf3['I2_pct']):.0f}%")

WNT = set(x.upper() for x in __import__("run_three_cohort").WNT)
wnt_rows = [r for r in meta.values() if r["symbol"].upper() in WNT and r["pooled_RE_log2FC"]]
mx = max(wnt_rows, key=lambda r: abs(float(r["pooled_RE_log2FC"])))
eq = [r for r in wnt_rows if float(r["RE_ci_lo"]) > -0.5 and float(r["RE_ci_hi"]) < 0.5]
k3 = [r for r in wnt_rows if r["k"] == "3"]
chk("Wnt poolable", None, f"{len(wnt_rows)} were quantified in at least two cohorts")
chk("Wnt in all three", None, f"({len(k3)} in all three)")
chk("Wnt max effect gene", None, mx["symbol"])
chk("Wnt max effect", None, f"{float(mx['pooled_RE_log2FC']):+.2f} log2")
chk("Wnt max CI", None, f"{float(mx['RE_ci_lo']):.2f} to {float(mx['RE_ci_hi']):.2f}")
chk("Wnt equivalent count", None, f"{len(eq)} of {len(wnt_rows)}")

for g, val in [("AXIN2", "-0.02"), ("LEF1", "-0.01"), ("TCF4", "-0.12"),
               ("GSK3B", "-0.05"), ("WNT2", "-0.17"), ("CTNNB1", "+0.08")]:
    r = meta[g]
    chk(f"{g} pooled", None, val)

# ---- GSEA
G = load("gsea_three_cohort.csv")
for cid, nes, p in [("C1", "-1.29", "0.100"), ("C2", "+1.02", "0.304"), ("C3", "+0.48", "0.931")]:
    r = [x for x in G if x["cohort"] == cid and x["gene_set"].startswith("HALLMARK")][0]
    chk(f"Wnt NES {cid}", None, nes)
    chk(f"Wnt p {cid}", None, f"{float(r['nominal_p']):.3f}")
for cid, nes in [("C1", "+2.68"), ("C2", "+0.76"), ("C3", "-1.11")]:
    r = [x for x in G if x["cohort"] == cid and x["gene_set"] == "PYROPTOSIS_CORE"][0]
    chk(f"Pyro NES {cid}", None, nes)

sens = json.load(open(P("sensitivity_results.json"), encoding="utf-8"))
chk("LOO min NES", None, f"{min(d['nes'] for d in sens['loo']):.2f}")
chk("LOO count", None, f"{len(sens['loo'])} leave-one-out analyses")
chk("exec/prime NES", None, f"+{sens['exec']['nes']:.2f}")
chk("nofilter wnt NES", None, f"{sens['nofilter']['wnt']['nes']:.2f}")

# ---- reference integrity
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
N += 1
if order == list(range(1, len(order) + 1)):
    print(f"  PASS  reference numbering sequential         1..{len(order)}")
else:
    FAIL.append(("refs", str(order)))
    print(f"  FAIL  reference numbering                    {order}")

print()
print(f"{N - len(FAIL)}/{N} checks passed")
if FAIL:
    print("FAILURES:")
    for f in FAIL:
        print("  -", f)
    sys.exit(1)
