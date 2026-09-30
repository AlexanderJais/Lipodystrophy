"""Over-representation analysis of DEG subsets via the Enrichr API.

Background is Enrichr's default (all annotated genes); brain-expressed genes
are therefore over-represented in every list, so treat terms as leads.
"""
import json
import time
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
DATA, OUT = ROOT / "data", ROOT / "results" / "enrichment"
OUT.mkdir(parents=True, exist_ok=True)
URL = "https://maayanlab.cloud/Enrichr"
LIBS = ["GO_Biological_Process_2023", "PanglaoDB_Augmented_2021",
        "WikiPathways_2019_Mouse", "KEGG_2019_Mouse"]
C = ["LDsaline_WT", "LDleptin_WT", "LDsaline_LDleptin"]


def lists():
    d = {(s, c): pd.read_csv(DATA / f"DEG{s}_{c}.csv").drop_duplicates("SYMBOL")
         .set_index("SYMBOL").log2FoldChange for s in ("male", "female") for c in C}
    ms, ml = d[("male", "LDsaline_WT")], d[("male", "LDleptin_WT")]
    fs = d[("female", "LDsaline_WT")]
    shared = ms.index.intersection(ml.index)
    ml_only = ml[~ml.index.isin(ms.index)]
    return {
        "male_saline_up": ms[ms > 0].index,
        "male_saline_down": ms[ms < 0].index,
        "male_normalized_up": ms[(ms > 0) & ~ms.index.isin(ml.index)].index,
        "male_normalized_down": ms[(ms < 0) & ~ms.index.isin(ml.index)].index,
        "male_persisting_down": shared[ms[shared] < 0],
        "male_leptin_only_up": ml_only[ml_only > 0].index,
        "male_leptin_only_down": ml_only[ml_only < 0].index,
        "female_saline_up": fs[fs > 0].index,
        "female_saline_down": fs[fs < 0].index,
        "female_saline_vs_leptin": d[("female", "LDsaline_LDleptin")].index,
    }


def enrich(genes, desc):
    r = requests.post(f"{URL}/addList", files={
        "list": (None, "\n".join(genes)), "description": (None, desc)}, timeout=60)
    uid = r.json()["userListId"]
    rows = []
    for lib in LIBS:
        res = requests.get(f"{URL}/enrich", params={
            "userListId": uid, "backgroundType": lib}, timeout=60).json()[lib]
        for t in res:
            rows.append({"library": lib, "term": t[1], "p": t[2],
                         "padj": t[6], "overlap": ";".join(t[5]),
                         "n": len(t[5])})
        time.sleep(0.5)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    for name, genes in lists().items():
        genes = [g for g in genes if not g.startswith(("Gm", "mt-", "Mir"))]
        df = enrich(genes, name)
        df.to_csv(OUT / f"{name}.csv", index=False)
        top = df[(df.padj < 0.05) & (df.n >= 3)].sort_values("padj")
        print(f"\n## {name} ({len(genes)} genes)")
        for lib, g in top.groupby("library", sort=False):
            for _, t in g.head(6).iterrows():
                print(f"  [{lib.split('_')[0]}] {t.term[:70]} | q={t.padj:.1e} | {t.overlap[:90]}")
