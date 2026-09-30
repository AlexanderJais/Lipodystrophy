"""Partition male LD saline-vs-WT DEGs by their fate under metreleptin and test
cell-type and transcription-factor-target enrichment of each part (Enrichr).

normalized  = DE in LD saline vs WT, not DE in LD leptin vs WT
persistent  = DE in LD saline vs WT and in LD leptin vs WT (same direction)
"""
import time
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "leptin_partition"
OUT.mkdir(parents=True, exist_ok=True)
URL = "https://maayanlab.cloud/Enrichr"
LIBS = ["PanglaoDB_Augmented_2021", "ChEA_2022",
        "TRRUST_Transcription_Factors_2019",
        "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]
ARTEFACTS = ("Gm", "mt-", "Mir", "Lars2", "Tcaf2")


def lfc(name):
    df = pd.read_csv(ROOT / "data" / f"DEG{name}.csv").drop_duplicates("SYMBOL")
    return df.set_index("SYMBOL").log2FoldChange


def partition():
    s, l = lfc("male_LDsaline_WT"), lfc("male_LDleptin_WT")
    s = s[~s.index.str.startswith(ARTEFACTS)]
    pers = s.index.intersection(l.index)
    norm = s.index.difference(l.index)
    return {"Normalized up": norm[s[norm] > 0], "Normalized down": norm[s[norm] < 0],
            "Persistent up": pers[s[pers] > 0], "Persistent down": pers[s[pers] < 0]}


def enrich(genes, desc):
    uid = requests.post(f"{URL}/addList", files={
        "list": (None, "\n".join(genes)), "description": (None, desc)},
        timeout=60).json()["userListId"]
    rows = []
    for lib in LIBS:
        res = requests.get(f"{URL}/enrich", params={
            "userListId": uid, "backgroundType": lib}, timeout=60).json()[lib]
        rows += [{"library": lib, "term": t[1], "p": t[2], "padj": t[6],
                  "n": len(t[5]), "genes": ";".join(t[5])} for t in res]
        time.sleep(0.5)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    parts = partition()
    pd.Series({k: len(v) for k, v in parts.items()}).to_csv(OUT / "sizes.csv")
    for name, genes in parts.items():
        df = enrich(list(genes), name)
        df.to_csv(OUT / f"{name.replace(' ', '_')}.csv", index=False)
        print(f"\n## {name} (n={len(genes)})")
        for lib, g in df[(df.padj < 0.05) & (df.n >= 3)].groupby("library", sort=False):
            for _, t in g.sort_values("padj").head(6).iterrows():
                print(f"  [{lib.split('_')[0]}] {t.term[:60]} q={t.padj:.1e} n={t.n}")
