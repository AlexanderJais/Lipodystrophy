"""log2FC of oligodendrocyte-stage markers and related gene sets across contrasts.

Stage markers follow Marques et al., Science 2016 (OPC -> COP -> NFOL -> MFOL -> MOL).
Blank = not in the DEG list (padj >= 0.1 or |log2FC| <= 0.2).
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
C = ["LDsaline_WT", "LDleptin_WT", "LDsaline_LDleptin"]
SETS = {
    "OPC": ["Pdgfra", "Cspg4", "Ptprz1", "Vcan", "Tnr", "Lhfpl3", "Matn4", "Sox6"],
    "COP": ["Neu4", "Bmp4", "Gpr17", "Sox4", "Sox11"],
    "NFOL": ["Bcas1", "Enpp6", "Tmem2", "Fyn", "Frmd4a", "Tcf7l2", "Itpr2", "Tns3"],
    "MFOL": ["Mal", "Opalin", "Plekhh1", "Ctps", "Mog", "Serinc5", "Prr18", "Pllp"],
    "MOL": ["Klk6", "Apod", "Grm3", "Slc5a11", "Ptgds", "Il33", "Anln", "Hapln2",
            "Car2", "Trf", "Cryab", "Mobp", "Aspa"],
    "OL transcription factors": ["Olig1", "Olig2", "Sox10", "Sox8", "Myrf", "Zeb2",
                                 "Qki", "Nkx2-2", "Nkx6-2"],
    "Myelin lipid synthesis": ["Ugt8a", "Fa2h", "Gal3st1", "Elovl7", "Elovl1",
                               "Cers2", "Gltp", "Pex5l"],
    "Cholesterol synthesis": ["Hmgcr", "Hmgcs1", "Sqle", "Dhcr24", "Dhcr7",
                              "Fdft1", "Cyp51", "Insig1", "Srebf2", "Scap"],
    "Axon / node": ["Nefl", "Nefm", "Nefh", "Nfasc", "Cntn1", "Cntn2",
                    "Cntnap1", "Kcna1", "Ank3"],
}

if __name__ == "__main__":
    d = {(s, c): pd.read_csv(ROOT / "data" / f"DEG{s}_{c}.csv")
         .drop_duplicates("SYMBOL").set_index("SYMBOL").log2FoldChange
         for s in ("male", "female") for c in C}
    rows = [{"set": k, "gene": g, **{f"{s}_{c}": d[(s, c)].get(g, np.nan)
                                     for s, c in d}}
            for k, gs in SETS.items() for g in gs]
    out = pd.DataFrame(rows).round(2)
    out.to_csv(ROOT / "results" / "oligodendrocyte_stages.csv", index=False)
    m, f = d[("male", "LDsaline_WT")], d[("female", "LDsaline_WT")]
    sh = m.index.intersection(f.index)
    print(out.fillna("").to_string(index=False))
    print(f"\nFemale vs male effect (shared LD saline vs WT genes, n={len(sh)}): "
          f"slope {np.polyfit(m[sh], f[sh], 1)[0]:.2f}")
