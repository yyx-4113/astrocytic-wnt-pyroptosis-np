#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Batch audit of reference DOIs against Crossref.
For each DOI present in manuscript REFERENCES, query Crossref and report the
resolved title + container + year, plus a match verdict against the cited item.
"""
import urllib.request, json, ssl, time

# (ref_label, cited_author_year_topic, doi)
REFS = [
 ("1","van Hecke 2014 Pain epidemiological","10.1016/j.pain.2013.11.013"),
 ("2","Finnerup 2015 Lancet Neurol pharmacotherapy","10.1016/S1474-4422(14)70251-0"),
 ("3","Feng 2022 J Neuroinflammation Wnt beta-catenin NLRP3","10.1186/s12974-022-02457-1"),
 ("4","Zhang 2013 J Clin Invest WNT neuropathic pain","10.1172/JCI65364"),
 ("5","Shi 2019 Int J Mol Sci Wnt signaling neuropathic pain","10.3390/ijms20246161"),
 ("6","Liu 2017 J Neurosci Dexmedetomidine NLRP3 astrocytes","10.1523/JNEUROSCI.3751-16.2017"),
 ("7","Bortolotti 2018 Trends Mol Med Inflammasomes trauma","10.1016/j.molmed.2018.02.005"),
 ("8","Wolf 2018 SCANPY Genome Biol","10.1186/s13059-017-1382-0"),
 ("9","Eastman 2017 OpenMM7 PLoS Comput Biol","10.1371/journal.pcbi.1005659"),
 ("10","Kolde 2012 Robust Rank Aggregation","10.1093/bioinformatics/btr709"),
 ("11","Langfelder 2008 WGCNA BMC Bioinformatics","10.1186/1471-2105-9-559"),
 ("13","Jin 2021 CellChat Nat Commun","10.1038/s41467-021-21246-9"),
 ("14","Cao 2019 organogenesis Nature","10.1038/s41586-019-0969-x"),
 ("15","Newman 2019 CIBERSORTx Nat Biotechnol","10.1038/s41587-019-0114-2"),
 ("16","Szklarczyk 2023 STRING Nucleic Acids Res","10.1093/nar/gkac1000"),
 ("17","Castro-Mondragon 2022 JASPAR Nucleic Acids Res","10.1093/nar/gkab1113"),
 ("18","Wishart 2018 DrugBank Nucleic Acids Res","10.1093/nar/gkx1037"),
 ("19","Margolin 2006 ARACNE BMC Bioinformatics","10.1186/1471-2105-7-S1-S7"),
 ("20","GTEx 2020 Science","10.1126/science.aaz1776"),
 ("21","Vosa 2021 Nat Genet eQTL","10.1038/s41588-021-00913-z"),
 ("22","Eastman 2023 Spice-OpenMM J Chem Theory Comput","10.1021/acs.jctc.3c00587"),
 ("29","Ji 2019 Astrocytes chronic pain Nat Rev Neurosci","10.1038/s41583-019-0218-1"),
]

ctx = ssl.create_default_context()
def query(doi):
    url = "https://api.crossref.org/works/" + doi
    try:
        req = urllib.request.Request(url, headers={"User-Agent":"refaudit/1.0"})
        with urllib.request.urlopen(req, timeout=25, context=ctx) as r:
            d = json.load(r)
        m = d.get("message", {})
        title = (m.get("title") or ["(no title)"])[0]
        cont = (m.get("container-title") or ["(no container)"])[0] if m.get("container-title") else "(no container)"
        yr = ""
        if m.get("published"):
            yr = str(m["published"].get("date-parts",[[None]])[0][0] or "")
        return "OK", title, cont, yr
    except urllib.error.HTTPError as e:
        return "HTTP_%d"%e.code, "(not found)", "", ""
    except Exception as e:
        return "ERR", str(e)[:60], "", ""

print("%-4s %-44s %-10s %-4s %s" % ("REF","CITED","STATUS","YEAR","RESOLVED TITLE"))
for lab, cited, doi in REFS:
    st, title, cont, yr = query(doi)
    print("%-4s %-44s %-10s %-4s %s | %s" % (lab, cited[:44], st, yr, title[:70], cont[:30]))
    time.sleep(0.3)
