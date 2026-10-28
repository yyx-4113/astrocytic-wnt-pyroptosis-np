# -*- coding: utf-8 -*-
"""Decisive test: are GSE217932 and GSE186237 raw count matrices the SAME
data under different column labels? If yes, the "two independent bulk
datasets" claim is untenable at the raw level."""
import openpyxl
import numpy as np

DATA = r"D:\省自然\生信分析\data"


def load(path):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb.active
    rows = ws.iter_rows(values_only=True)
    header = list(next(rows))
    genes, mat = [], []
    for r in rows:
        genes.append(r[0])
        mat.append([0.0 if v is None else float(v) for v in r[1:]])
    return header, genes, np.array(mat, dtype=float)


print("loading...", flush=True)
h1, g1, m1 = load(f"{DATA}\\GSE217932_Raw_gene_counts_matrix.xlsx")
h2, g2, m2 = load(f"{DATA}\\GSE186237_Raw_gene_counts_matrix.xlsx")

# column index maps (0-based within sample block)
s217_sham = [h1[1:].index(c) for c in h1[1:] if "Sham" in c]
s217_cci = [h1[1:].index(c) for c in h1[1:] if c.startswith("CCI-")]
s186_s = [h2[1:].index(c) for c in h2[1:] if c.startswith("S-")]
s186_cci = [h2[1:].index(c) for c in h2[1:] if c.startswith("CCI-")]

g1i = {g: i for i, g in enumerate(g1)}
common = [g for g in g2 if g in g1i]
print(f"common genes = {len(common)}", flush=True)

# build paired 8-vector per gene: [Sham1-4, CCI1-4] vs [S1-4, CCI1-4]
A, B = [], []
for g in common:
    i1, i2 = g1i[g], g2.index(g)
    a = np.concatenate([m1[i1, s217_sham], m1[i1, s217_cci]])
    b = np.concatenate([m2[i2, s186_s], m2[i2, s186_cci]])
    A.append(a); B.append(b)
A, B = map(np.array, (A, B))
print(f"A shape {A.shape}, B shape {B.shape}", flush=True)

identical = np.all(np.abs(A - B) < 1e-6, axis=1)
print(f"\n[RAW DUPLICATE TEST]")
print(f"  genes where all 8 raw counts identical (<1e-6): "
      f"{identical.sum()}/{len(common)} = {identical.mean()*100:.2f}%")
flat_corr = np.corrcoef(A.ravel(), B.ravel())[0, 1]
print(f"  flattened raw-count correlation = {flat_corr:.6f}")

# per-column check
for label, a_cols, b_cols in [("Sham vs S", s217_sham, s186_s),
                              ("CCI vs CCI", s217_cci, s186_cci)]:
    ca = m1[[g1i[g] for g in common]][:, a_cols].ravel()
    cb = m2[[g2.index(g) for g in common]][:, b_cols].ravel()
    print(f"  {label} raw-count correlation = {np.corrcoef(ca, cb)[0,1]:.6f}")
print("\nDONE.", flush=True)
