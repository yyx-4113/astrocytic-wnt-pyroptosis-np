#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Definitive comparison between the GENUINE GSE186237 count matrix (downloaded from GEO)
and GSE217932 count matrix. No trust of prior summaries.
Reads both xlsx with openpyxl, aligns by gene_id, and reports:
  - shapes, column names
  - gene_id format (first 5 IDs of each)
  - fraction of identical count values across common genes
  - overall count correlation
  - per-column (sample) identity
"""
import openpyxl, math

def load(path):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb.active
    rows = ws.iter_rows(values_only=True)
    header = list(next(rows))
    gene_ids = []
    mat = []   # list of lists (rows)
    for r in rows:
        r = list(r)
        if r[0] is None:
            continue
        gene_ids.append(str(r[0]).strip())
        mat.append([ (float(x) if x is not None else float('nan')) for x in r[1:]])
    wb.close()
    return header, gene_ids, mat

P1 = r'D:/省自然/生信分析/data/GSE217932_Raw_gene_counts_matrix.xlsx'
P2 = r'D:/省自然/生信分析/data/GSE186237_genuine_GEO.xlsx'

h1, g1, m1 = load(P1)
h2, g2, m2 = load(P2)

print('=== GSE217932 ===')
print('shape(rows,cols):', len(m1), len(h1)-1)
print('col0:', repr(h1[0]))
print('sample cols:', h1[1:])
print('first 5 gene_ids:', g1[:5])
print('=== GSE186237 (genuine GEO) ===')
print('shape(rows,cols):', len(m2), len(h2)-1)
print('col0:', repr(h2[0]))
print('sample cols:', h2[1:])
print('first 5 gene_ids:', g2[:5])

# Build dict for g1
d1 = {}
for gid, row in zip(g1, m1):
    d1[gid] = row
common = [g for g in g2 if g in d1]
print('\ncommon gene_ids:', len(common), '| g1-only:', len(set(g1)-set(g2)), '| g2-only:', len(set(g2)-set(g1)))

# Compare counts
n = len(common)
ident = 0
total = 0
import statistics
diffs = []
for g in common:
    r2 = dict(zip(h2[1:], m2[g2.index(g)])) if False else None
# simpler: index by position
idx2 = {g:i for i,g in enumerate(g2)}
for g in common:
    row1 = d1[g]
    row2 = m2[idx2[g]]
    for a,b in zip(row1, row2):
        total += 1
        if a==b:
            ident += 1
        else:
            diffs.append(abs(a-b))
print('compared value cells:', total)
print('identical cells: %d (%.4f)' % (ident, ident/total))
if diffs:
    print('max abs diff:', max(diffs), '| mean abs diff:', sum(diffs)/len(diffs))
