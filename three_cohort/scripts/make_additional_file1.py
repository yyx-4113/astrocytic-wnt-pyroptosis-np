# -*- coding: utf-8 -*-
"""Build Additional file 1 (XLSX): S1 per-cohort statistics, S2 pooled estimates,
S3 GSEA output incl. pre-processing variants, S4 gene-panel membership."""
import csv, json, os
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill
from openpyxl.utils import get_column_letter

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE) if os.path.basename(HERE).lower() == "scripts" else HERE
OUT = os.path.join(os.environ.get("NP_OUT", PKG), "Additional_File_1_tables.xlsx")
RESULTS = os.path.join(PKG, "results") if os.path.isdir(os.path.join(PKG, "results")) else PKG


def _find(fn):
    for d in (RESULTS, PKG, os.getcwd()):
        p = os.path.join(d, fn)
        if os.path.exists(p):
            return p
    return os.path.join(RESULTS, fn)  # not found; fail with a clear path
COH = [("C1", "cohort_C1_GSE217932_CCI_gene_level.csv", "GSE217932 (CCI, 4 vs 4)"),
       ("C2", "cohort_C2_GSE256472_SNI_gene_level.csv", "GSE256472 (SNI, 3 vs 3)"),
       ("C3", "cohort_C3_GSE205906_SNI_gene_level.csv", "GSE205906 (SNI, 3 vs 3)")]

HDR_FILL = PatternFill("solid", fgColor="1F3864")
HDR_FONT = Font(bold=True, color="FFFFFF", size=9, name="Arial")
BODY_FONT = Font(size=9, name="Arial")


def write_sheet(wb, title, header, rows, widths=None, first=False):
    ws = wb.active if first else wb.create_sheet()
    ws.title = title
    ws.append(header)
    for c in range(1, len(header) + 1):
        cell = ws.cell(row=1, column=c)
        cell.fill = HDR_FILL
        cell.font = HDR_FONT
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
    for r in rows:
        ws.append(r)
    for row in ws.iter_rows(min_row=2, max_row=ws.max_row):
        for cell in row:
            cell.font = BODY_FONT
    if widths:
        for i, w in enumerate(widths, start=1):
            ws.column_dimensions[get_column_letter(i)].width = w
    ws.freeze_panes = "A2"
    ws.auto_filter.ref = ws.dimensions
    return ws


def main():
    wb = Workbook()

    # ---------------- S1 per-cohort gene-level statistics
    hdr = ["cohort", "dataset", "symbol", "log2FC", "t_stat", "df", "pvalue",
           "padj", "ci_lo", "ci_hi", "mean_expr", "mean_control", "mean_case"]
    rows = []
    for cid, fn, label in COH:
        with open(_find(fn), encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                rows.append([cid, label, r["symbol"], r["log2FC"], r["t_stat"],
                             r["df"], r["pvalue"], r["padj"], r["ci_lo"],
                             r["ci_hi"], r["mean_expr"], r["mean_control"],
                             r["mean_case"]])
    write_sheet(wb, "S1 per-cohort", hdr, rows,
                widths=[7, 26, 14, 10, 12, 8, 11, 11, 10, 10, 11, 13, 11],
                first=True)

    # ---------------- S2 pooled estimates
    hdr2 = ["symbol", "k", "cohorts", "lfc_C1_CCI", "lfc_C2_SNI", "lfc_C3_SNI",
            "pooled_RE_log2FC", "RE_ci_lo", "RE_ci_hi", "p_RE", "I2_pct",
            "tau2", "pooled_FE_log2FC", "HKSJ_ci_lo", "HKSJ_ci_hi", "p_HKSJ",
            "equiv_within_0.5_normal", "equiv_within_0.5_HKSJ"]
    rows2 = []
    with open(_find("meta_focus.csv"), encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            rows2.append([r["symbol"], r["k"], r["cohorts"], r["lfc_C1_CCI"],
                          r["lfc_C2_SNI"], r["lfc_C3_SNI"], r["pooled_RE_log2FC"],
                          r["RE_ci_lo"], r["RE_ci_hi"], r["p_RE"], r["I2_pct"],
                          r["tau2"], r["pooled_FE_log2FC"],
                          r.get("HKSJ_ci_lo", ""), r.get("HKSJ_ci_hi", ""),
                          r.get("p_HKSJ", ""),
                          r.get("equiv_normal", ""), r.get("equiv_HKSJ", "")])
    write_sheet(wb, "S2 pooled", hdr2, rows2,
                widths=[12, 5, 14, 11, 11, 11, 16, 10, 10, 10, 8, 8, 15,
                        10, 10, 10, 12, 12])

    # ---------------- S3 GSEA output
    hdr3 = ["cohort", "gene_set", "set_source", "k", "genes_ranked", "ES", "NES",
            "nominal_p", "q_across_sets", "pct_of_set_found"]
    rows3 = []
    with open(_find("gsea_three_cohort.csv"), encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            rows3.append([r["cohort"], r["gene_set"], r["set_source"],
                          r["n_genes_in_set"], r["n_ranked"], r["ES"], r["NES"],
                          r["nominal_p"], r["FDR_q_across_sets_in_cohort"],
                          r["pct_of_set_found"]])
    # pre-processing variants and the execution/priming split (C1 primary settings)
    sens = json.load(open(_find("sensitivity_results.json"), encoding="utf-8"))
    # primary-setting values are taken from gsea_three_cohort.csv itself, never hard-coded
    primary_rows = [r for r in rows3 if r[0] == "C1"]
    primary = {r[1]: r for r in primary_rows}
    variants = [("primary", None), ("log2FC ranking", "lfc"),
                ("median expression filter", "med"),
                ("no expression filter", "nofilter")]
    for label, key in variants:
        if key is None:
            for setname in ("PYROPTOSIS_CORE", "HALLMARK_WNT_BETA_CATENIN"):
                r = primary[setname]
                rows3.append(["C1", setname, "primary setting"] + r[3:])
            continue
        for setname, disp in [("pyro", "PYROPTOSIS_CORE"), ("wnt", "HALLMARK_WNT_BETA_CATENIN"),
                              ("exec", "PYRO_EXECUTION_MACHINERY")]:
            v = sens.get(key, {}).get(setname)
            if v is None:
                rows3.append(["C1", disp, f"variant: {label}", "n/a", "n/a",
                              "n/a", "not testable (k<5)", "n/a", "n/a", "n/a"])
            else:
                rows3.append(["C1", disp, f"variant: {label}", v["k"], v["N"],
                              f"{v['es']:.4f}", f"{v['nes']:.4f}", f"{v['p']:.4g}",
                              "n/a (single set)", f"{v['frac_genes']:.1f}"])
    for key, disp in [("exec", "PYRO_EXECUTION_MACHINERY"), ("prime", "PYRO_PRIMING_SUBSTRATE")]:
        v = sens[key]
        rows3.append(["C1", disp, "execution/priming split", v["k"], v["N"],
                      f"{v['es']:.4f}", f"{v['nes']:.4f}", f"{v['p']:.4g}",
                      "n/a (single set)", f"{v['frac_genes']:.1f}"])
    for cid in ("C2", "C3"):
        for key, disp in [("exec", "PYRO_EXECUTION_MACHINERY"), ("prime", "PYRO_PRIMING_SUBSTRATE")]:
            v = sens[cid][key]
            rows3.append([cid, disp, "execution/priming split", v["k"], v["N"],
                          f"{v['es']:.4f}", f"{v['nes']:.4f}", f"{v['p']:.4g}",
                          "n/a (single set)", f"{v['frac_genes']:.1f}"])
    write_sheet(wb, "S3 GSEA", hdr3, rows3,
                widths=[8, 30, 26, 6, 13, 9, 9, 10, 12, 16])

    # ---------------- S4 panel membership
    panels = json.load(open(_find("gene_panels.json"), encoding="utf-8"))
    rows4 = []
    for name, d in panels.items():
        for g in d["genes"]:
            rows4.append([name, d["source"], d["n"], g])
    write_sheet(wb, "S4 panels", ["gene_set", "source", "n_genes", "symbol"], rows4,
                widths=[30, 34, 8, 14])

    # ---------------- S5 sample-to-GSM mapping, including excluded libraries
    # GSM identifiers are taken from the GEO esummary record for each series.
    rows5 = [
        ["C1", "GSE217932", "Sham", "GSM6730101", "included", "sciatic nerve exposed, not ligated"],
        ["C1", "GSE217932", "Sham", "GSM6730102", "included", "sciatic nerve exposed, not ligated"],
        ["C1", "GSE217932", "Sham", "GSM6730103", "included", "sciatic nerve exposed, not ligated"],
        ["C1", "GSE217932", "Sham", "GSM6730104", "included", "sciatic nerve exposed, not ligated"],
        ["C1", "GSE217932", "CCI", "GSM6730105", "included", "chronic constriction injury"],
        ["C1", "GSE217932", "CCI", "GSM6730106", "included", "chronic constriction injury"],
        ["C1", "GSE217932", "CCI", "GSM6730107", "included", "chronic constriction injury"],
        ["C1", "GSE217932", "CCI", "GSM6730108", "included", "chronic constriction injury"],
        ["C1", "GSE217932", "CCI + l-THP", "GSM6730109", "excluded", "intervention arm"],
        ["C1", "GSE217932", "CCI + l-THP", "GSM6730110", "excluded", "intervention arm"],
        ["C1", "GSE217932", "CCI + l-THP", "GSM6730111", "excluded", "intervention arm"],
        ["C1", "GSE217932", "CCI + l-THP", "GSM6730112", "excluded", "intervention arm"],
        ["C2", "GSE256472", "Control", "GSM8099011", "included", ""],
        ["C2", "GSE256472", "Control", "GSM8099012", "included", ""],
        ["C2", "GSE256472", "Control", "GSM8099013", "included", ""],
        ["C2", "GSE256472", "SNI", "GSM8099014", "included", "spared nerve injury"],
        ["C2", "GSE256472", "SNI", "GSM8099015", "included", "spared nerve injury"],
        ["C2", "GSE256472", "SNI", "GSM8099016", "included", "spared nerve injury"],
        ["C3", "GSE205906", "JA (sham)", "GSM6235292", "included", "sham-operated, oral normal saline"],
        ["C3", "GSE205906", "JA (sham)", "GSM6235293", "included", "sham-operated, oral normal saline"],
        ["C3", "GSE205906", "JA (sham)", "GSM6235294", "included", "sham-operated, oral normal saline"],
        ["C3", "GSE205906", "SNI", "GSM6235295", "included", "spared nerve injury, oral normal saline"],
        ["C3", "GSE205906", "SNI", "GSM6235296", "included", "spared nerve injury, oral normal saline"],
        ["C3", "GSE205906", "SNI", "GSM6235297", "included", "spared nerve injury, oral normal saline"],
        ["C3", "GSE205906", "PLM", "GSM6235298", "excluded", "paeoniflorin-liquiritin arm"],
        ["C3", "GSE205906", "PLM", "GSM6235299", "excluded", "paeoniflorin-liquiritin arm"],
        ["C3", "GSE205906", "PLM", "GSM6235300", "excluded", "paeoniflorin-liquiritin arm"],
    ]
    write_sheet(wb, "S5 samples",
                ["cohort", "series", "group", "GSM", "status", "note"], rows5,
                widths=[8, 12, 14, 14, 10, 40])

    wb.save(OUT)
    print("saved", OUT)
    print("S1 rows:", len(rows), "S2 rows:", len(rows2), "S3 rows:", len(rows3),
          "S4 rows:", len(rows4), "S5 rows:", len(rows5))


if __name__ == "__main__":
    main()
