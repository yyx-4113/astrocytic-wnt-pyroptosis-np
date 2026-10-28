# -*- coding: utf-8 -*-
"""Rebuild 01 — faithfully re-run bulk DEG for GSE217932 (validate method
against the existing valid CSV) and GSE186237 (fresh, correct grouping),
then determine whether the two datasets genuinely diverge.

Context: pipeline_full.py stored GSE186237_DEGs.csv byte-identical to
GSE217932_DEGs.csv (MD5 a7b4df4f...). This script checks whether, using the
CURRENT distinct raw xlsx files, the two datasets produce different DEGs.
"""
import openpyxl
import numpy as np
from scipy import stats

DATA = r"D:\省自然\生信分析\data"
RESULTS = r"D:\省自然\生信分析\results"


def load_counts(path):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb.active
    rows = ws.iter_rows(values_only=True)
    header = list(next(rows))
    genes, mat = [], []
    for r in rows:
        genes.append(r[0])
        mat.append([0.0 if v is None else float(v) for v in r[1:]])
    return header, genes, np.array(mat, dtype=float)


def group_idx(header, sham_pred, cci_pred):
    sham = [i for i, c in enumerate(header[1:]) if sham_pred(c)]
    cci = [i for i, c in enumerate(header[1:]) if cci_pred(c)]
    return sham, cci


def normalize(mat, sham_idx, cci_idx):
    cols = sorted(set(sham_idx) | set(cci_idx))
    sub = mat[:, cols]
    keep = (sub > 5).sum(axis=1) >= 3
    sub = sub[keep]
    gm = np.exp(np.log(sub + 1).mean(axis=1))
    sf = np.median(sub / gm[:, None], axis=0)
    norm = sub / sf[None, :]
    return norm, keep, cols


def run_de(norm, sham_local, cci_local):
    log2 = np.log2(norm + 1)
    c, s = log2[:, cci_local], log2[:, sham_local]
    log2fc = c.mean(axis=1) - s.mean(axis=1)
    _, p = stats.ttest_ind(c, s, axis=1, equal_var=False)
    p = np.where(np.isnan(p), 1.0, p)
    order = np.argsort(p)
    m = len(p)
    bh = p[order] * m / np.arange(1, m + 1)
    bh = np.minimum.accumulate(bh[::-1])[::-1]
    padj = np.empty(m)
    padj[order] = np.clip(bh, 0, 1)
    return log2fc, p, padj


def deg_table(header, genes, mat, sham_pred, cci_pred):
    sham, cci = group_idx(header, sham_pred, cci_pred)
    norm, keep, cols = normalize(mat, sham, cci)
    log2fc, p, padj = run_de(norm,
                             [cols.index(i) for i in sham],
                             [cols.index(i) for i in cci])
    return (np.array(genes)[keep], log2fc, p, padj,
            len(sham), len(cci))


def load_csv_dict(path):
    import csv
    d = {}
    with open(path, newline="", encoding="utf-8") as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            d[row[0]] = float(row[1])  # gene_id -> log2FC
    return d


print("Loading raw xlsx...")
h1, g1, m1 = load_counts(f"{DATA}\\GSE217932_Raw_gene_counts_matrix.xlsx")
h2, g2, m2 = load_counts(f"{DATA}\\GSE186237_Raw_gene_counts_matrix.xlsx")

# GSE217932 valid (Sham vs CCI, excludes CCI+THP)
g217 = deg_table(h1, g1, m1,
                 lambda c: "Sham" in c,
                 lambda c: c.startswith("CCI-"))
# GSE186237 correct grouping (S-* = Sham, CCI-* = injury)
g186 = deg_table(h2, g2, m2,
                 lambda c: c.startswith("S-"),
                 lambda c: c.startswith("CCI-"))

# --- validate method on GSE217932 vs existing CSV ---
csv217 = load_csv_dict(f"{RESULTS}\\GSE217932_DEGs.csv")
common = set(g217[0]) & set(csv217)
vals_mine, vals_csv = [], []
for g in common:
    i = list(g217[0]).index(g)
    vals_mine.append(g217[1][i])
    vals_csv.append(csv217[g])
vals_mine, vals_csv = map(np.array, (vals_mine, vals_csv))
print(f"\n[Validate] GSE217932: {len(common)} common genes vs existing CSV")
print(f"  Pearson r(log2FC) = {np.corrcoef(vals_mine, vals_csv)[0,1]:.4f}")
print(f"  max|diff| = {np.max(np.abs(vals_mine-vals_csv)):.4f}")

# --- compare GSE217932 vs freshly computed GSE186237 ---
genes1, fc1, p1, padj1, n1s, n1c = g217
genes2, fc2, p2, padj2, n2s, n2c = g186
print(f"\n[GSE217932] n_genes={len(genes1)} Sham={n1s} CCI={n1c}")
print(f"[GSE186237] n_genes={len(genes2)} Sham={n2s} CCI={n2c}")

# sig sets (padj<0.05 & |log2FC|>0.5)
sig1 = set(np.array(genes1)[(padj1 < 0.05) & (np.abs(fc1) > 0.5)])
up1 = set(np.array(genes1)[(padj1 < 0.05) & (fc1 > 0.5)])
dn1 = set(np.array(genes1)[(padj1 < 0.05) & (fc1 < -0.5)])
sig2 = set(np.array(genes2)[(padj2 < 0.05) & (np.abs(fc2) > 0.5)])
up2 = set(np.array(genes2)[(padj2 < 0.05) & (fc2 > 0.5)])
dn2 = set(np.array(genes2)[(padj2 < 0.05) & (fc2 < -0.5)])
print(f"\n[Significant DEGs |log2FC|>0.5, padj<0.05]")
print(f"  GSE217932: total={len(sig1)} up={len(up1)} down={len(dn1)}")
print(f"  GSE186237: total={len(sig2)} up={len(up2)} down={len(dn2)}")

# common-gene log2FC correlation
gc = set(genes1) & set(genes2)
d1 = dict(zip(genes1, fc1))
d2 = dict(zip(genes2, fc2))
x = np.array([d1[g] for g in gc])
y = np.array([d2[g] for g in gc])
r = np.corrcoef(x, y)[0, 1]
print(f"\n[Cross-dataset log2FC correlation over {len(gc)} common genes] r = {r:.4f}")

# overlap of significant sets
inter = sig1 & sig2
union = sig1 | sig2
print(f"[Significant-gene overlap] intersect={len(inter)} union={len(union)} "
      f"Jaccard={len(inter)/max(1,len(union)):.3f}")

# --- show the STALE artifact ---
stale = load_csv_dict(f"{RESULTS}\\GSE186237_DEGs.csv")
cs = set(g217[0]) & set(stale)
a = np.array([d1[g] for g in cs])
b = np.array([stale[g] for g in cs])
print(f"\n[STALE artifact check] GSE186237_DEGs.csv vs GSE217932, "
      f"r = {np.corrcoef(a,b)[0,1]:.4f} (==1.0 means identical/duplicated)")
print("\nDONE.")
