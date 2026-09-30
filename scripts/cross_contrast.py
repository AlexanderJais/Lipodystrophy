"""Cross-contrast summary of hypothalamic DEG lists (DESeq2 output, pre-filtered).

Contrast naming: DEG<sex>_<A>_<B>.csv -> log2FoldChange = A vs B.
"""
from pathlib import Path
import pandas as pd

DATA = Path(__file__).resolve().parents[1] / "data"
OUT = Path(__file__).resolve().parents[1] / "results"
OUT.mkdir(exist_ok=True)

CONTRASTS = [f"{s}_{c}" for s in ("male", "female")
             for c in ("LDsaline_WT", "LDleptin_WT", "LDsaline_LDleptin")]

PANELS = {
    "arcuate_energy": ["Agrp", "Npy", "Pomc", "Cartpt", "Lepr", "Socs3", "Vgf",
                       "Ghrh", "Sst", "Th", "Kiss1", "Nhlh2"],
    "OPC_OL_lineage": ["Olig1", "Olig2", "Sox10", "Pdgfra", "Cspg4", "Myrf",
                       "Nkx6-2", "Bcas1", "Enpp6", "Tcf7l2", "Gpr17"],
    "myelin_structural": ["Mbp", "Plp1", "Mog", "Mag", "Mobp", "Cnp", "Cldn11",
                          "Mal", "Opalin", "Ermn", "Ugt8a", "Fa2h", "Gjc2",
                          "Gjb1", "Gjc3", "Pllp", "Tmem125", "Aspa", "Trf",
                          "Hapln2", "Prr18", "Tmem278", "Qdpr", "S1pr5", "Anln"],
    "Ch25h_NLRP3": ["Ch25h", "Cyp7b1", "Nlrp3", "Pycard", "Casp1", "Il1b",
                    "Il18", "Gsdmd", "Nr1h3", "Abca1"],
    "inflammation_glia": ["Lcn2", "Icam1", "Tgm2", "Gfap", "Serpina3n", "C4b",
                          "Cd74", "H2-Aa", "H2-Ab1", "Cx3cr1", "Aif1", "Il33",
                          "Il1r1", "Socs3"],
    "dissection_markers": ["Kcnj13", "Ttr", "Lhx8", "Vip", "Gpr88", "Slc18a2",
                           "Ddc", "Aldh1a1", "Ret", "Pvalb"],
}


def load():
    frames = {}
    for c in CONTRASTS:
        df = pd.read_csv(DATA / f"DEG{c}.csv")
        frames[c] = df.drop_duplicates("SYMBOL").set_index("SYMBOL")
    return frames


def main():
    frames = load()
    lfc = pd.DataFrame({c: f["log2FoldChange"] for c, f in frames.items()})
    padj = pd.DataFrame({c: f["padj"] for c, f in frames.items()})
    lfc.round(3).to_csv(OUT / "lfc_matrix_all_DEGs.csv")

    summary = []
    for c, f in frames.items():
        summary.append({"contrast": c, "n_DEG": len(f),
                        "up": int((f.log2FoldChange > 0).sum()),
                        "down": int((f.log2FoldChange < 0).sum()),
                        "max_padj": f.padj.max(),
                        "min_abs_lfc": f.log2FoldChange.abs().min()})
    pd.DataFrame(summary).to_csv(OUT / "contrast_summary.csv", index=False)
    print(pd.DataFrame(summary).to_string(index=False), "\n")

    rows = []
    for panel, genes in PANELS.items():
        for g in genes:
            r = {"panel": panel, "gene": g}
            for c in CONTRASTS:
                r[c] = round(lfc.at[g, c], 2) if g in lfc.index and pd.notna(lfc.at[g, c]) else None
            rows.append(r)
    panel_df = pd.DataFrame(rows)
    panel_df.to_csv(OUT / "panel_lfc.csv", index=False)
    with pd.option_context("display.width", 200, "display.max_rows", 200):
        print(panel_df.fillna("").to_string(index=False), "\n")

    # Reversal test: genes changed LDsaline vs WT -> behaviour in the other two contrasts.
    for sex in ("male", "female"):
        s = frames[f"{sex}_LDsaline_WT"]
        lw = frames[f"{sex}_LDleptin_WT"]
        sl = frames[f"{sex}_LDsaline_LDleptin"]
        persists = [g for g in s.index if g in lw.index
                    and (s.at[g, "log2FoldChange"] > 0) == (lw.at[g, "log2FoldChange"] > 0)]
        rescued = [g for g in s.index if g in sl.index
                   and (s.at[g, "log2FoldChange"] > 0) == (sl.at[g, "log2FoldChange"] > 0)]
        print(f"{sex}: LDsaline-vs-WT DEGs={len(s)}; still DE same-direction in "
              f"LDleptin-vs-WT={len(persists)}; significantly changed by leptin "
              f"(saline vs leptin, concordant)={len(rescued)}")
        print("  rescued:", ", ".join(rescued))
        pd.Series(persists, name="gene").to_csv(OUT / f"{sex}_persisting_despite_leptin.csv", index=False)
        pd.Series(rescued, name="gene").to_csv(OUT / f"{sex}_leptin_rescued.csv", index=False)


if __name__ == "__main__":
    main()
