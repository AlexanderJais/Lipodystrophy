"""Recover per-mouse values from the vector phenotype PDFs (ggplot output).

Each dot is a filled path whose colour encodes the group; its centre is mapped
to data units with a linear fit to the y-axis tick labels of the same panel.
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pymupdf

ROOT = Path(__file__).resolve().parents[1]
GROUP = {(0.42, 0.56, 0.14): "WT", (1.0, 0.65, 0.0): "LD saline",
         (1.0, 0.84, 0.0): "LD leptin"}
TITLE = re.compile(r"(Female|Male) mice \(week (\d+)\)")


def panels(page):
    """Panel titles with their top-left anchor."""
    out = []
    for b in page.get_text("dict")["blocks"]:
        for line in b.get("lines", []):
            text = "".join(s["text"] for s in line["spans"])
            m = TITLE.search(text)
            if m:
                out.append((line["bbox"][0], line["bbox"][1], m.group(1),
                            int(m.group(2))))
    return out


def owner(x, y, titles):
    cand = [t for t in titles if t[1] < y and t[0] - 60 < x]
    return max(cand, key=lambda t: (t[1], t[0])) if cand else None


def extract(pdf):
    rows = []
    for page in pymupdf.open(pdf):
        titles = panels(page)
        words = page.get_text("words")
        # y-axis title is the only rotated text on the page
        measure = next(
            "".join(s["text"] for s in line["spans"]).strip()
            for b in page.get_text("dict")["blocks"]
            for line in b.get("lines", []) if line["dir"][1] != 0)
        ticks = {}
        for w in words:
            try:
                v = float(w[4].replace("−", "-"))
            except ValueError:
                continue
            t = owner(w[2], (w[1] + w[3]) / 2, titles)
            if t and w[2] < t[0] + 40:  # y-axis tick labels sit left of the plot
                ticks.setdefault(t, []).append(((w[1] + w[3]) / 2, v))
        for d in page.get_drawings():
            fill = tuple(round(c, 2) for c in d.get("fill") or ())
            if fill not in GROUP or d["rect"].width > 8:
                continue
            cx, cy = (d["rect"].x0 + d["rect"].x1) / 2, (d["rect"].y0 + d["rect"].y1) / 2
            t = owner(cx, cy, titles)
            if t is None or len(ticks.get(t, [])) < 2:
                continue
            ys, vs = zip(*ticks[t])
            a, b = np.polyfit(ys, vs, 1)
            rows.append({"measure": measure, "sex": t[2].lower(), "week": t[3],
                         "group": GROUP[fill], "value": a * cy + b})
    return rows


if __name__ == "__main__":
    rows = []
    for pdf in sys.argv[1:]:
        rows += extract(pdf)
    df = pd.DataFrame(rows)
    df.to_csv(ROOT / "data" / "phenotype_extracted.csv", index=False)
    print(df.groupby(["measure", "sex", "week", "group"]).value
          .agg(["count", "median"]).round(1).to_string())
