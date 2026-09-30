"""Publication figures from the hypothalamic DEG lists.

Figure 1 - leptin-reversible energy-deficit signature (males)
Figure 2 - oligodendrocyte / myelin deficit that persists under metreleptin

Colour encodes log2 fold change; dot area encodes -log10 padj. Genes absent from
a DEG list (padj >= 0.1 or |log2FC| <= 0.2) are drawn as small open grey rings.
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
DATA, OUT = ROOT / "data", ROOT / "figures"
OUT.mkdir(exist_ok=True)

MM = 1 / 25.4
INK, RULE = "#1a1a1a", "#bdbdbd"
DOWN, MID, UP = "#1D4E89", "#F2F2F0", "#A4243B"
CMAP = LinearSegmentedColormap.from_list(
    "div", ["#0F2F57", DOWN, "#8FB3D9", MID, "#E3A1A9", UP, "#5E0F1E"])

mpl.rcParams.update({
    "font.family": "Nimbus Sans", "font.size": 6, "axes.labelsize": 6,
    "axes.titlesize": 6.5, "xtick.labelsize": 5.5, "ytick.labelsize": 5.5,
    "legend.fontsize": 5.5, "axes.linewidth": 0.5, "xtick.major.width": 0.5,
    "ytick.major.width": 0.5, "xtick.major.size": 2, "ytick.major.size": 2,
    "axes.edgecolor": INK, "text.color": INK, "axes.labelcolor": INK,
    "xtick.color": INK, "ytick.color": INK, "axes.spines.top": False,
    "axes.spines.right": False, "pdf.fonttype": 42, "svg.fonttype": "none",
    "savefig.dpi": 600, "mathtext.fontset": "custom",
    "mathtext.rm": "Nimbus Sans", "mathtext.it": "Nimbus Sans:italic",
    "mathtext.bf": "Nimbus Sans:bold",
})

LABEL = {"LDsaline_WT": "LD saline\nvs WT", "LDleptin_WT": "LD leptin\nvs WT",
         "LDsaline_LDleptin": "LD saline\nvs LD leptin"}
COMPS = list(LABEL)


def load():
    out = {}
    for sex in ("male", "female"):
        for c in COMPS:
            df = pd.read_csv(DATA / f"DEG{sex}_{c}.csv").drop_duplicates("SYMBOL")
            out[(sex, c)] = df.set_index("SYMBOL")
    return out


def panel_letter(ax, s, dx=-0.02, dy=1.0):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=8, fontweight="bold",
            va="bottom", ha="right")


def dot_heatmap(ax, deg, modules, columns, vlim=2.0):
    """modules: list of (name, [genes]); columns: list of (sex, contrast)."""
    rows, bounds = [], []
    for name, genes in modules:
        bounds.append((name, len(rows), len(rows) + len(genes) - 1))
        rows += genes
    norm = TwoSlopeNorm(0, -vlim, vlim)
    for j, col in enumerate(columns):
        d = deg[col]
        for i, g in enumerate(rows):
            if g in d.index:
                lfc, q = d.at[g, "log2FoldChange"], d.at[g, "padj"]
                size = 6 + 9 * min(-np.log10(q), 8)
                ax.scatter(j, i, s=size, c=[CMAP(norm(lfc))], lw=0.3,
                           edgecolors=INK, zorder=3)
            else:
                ax.scatter(j, i, s=4, facecolors="none", edgecolors=RULE,
                           lw=0.5, zorder=2)
    ax.set_xlim(-0.6, len(columns) - 0.4)
    ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows, fontstyle="italic")
    ax.set_xticks(range(len(columns)))
    ax.set_xticklabels([LABEL[c] for _, c in columns], rotation=0)
    ax.tick_params(axis="x", length=0, pad=3)
    ax.tick_params(axis="y", length=0, pad=2)
    for s in ax.spines.values():
        s.set_visible(False)
    for k, (name, a, b) in enumerate(bounds):
        if k:
            ax.axhline(a - 0.5, color=RULE, lw=0.4, zorder=1)
        ax.annotate(name, xy=(len(columns) - 0.35, (a + b) / 2),
                    xycoords="data", ha="left", va="center", fontsize=5.5,
                    color=INK, annotation_clip=False)
    return norm


def size_colour_legends(fig, norm, rect_cb, rect_sz, vlim=2.0):
    cax = fig.add_axes(rect_cb)
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=CMAP), cax=cax,
                      orientation="horizontal", ticks=[-vlim, -1, 0, 1, vlim])
    cb.outline.set_linewidth(0.4)
    cb.ax.tick_params(length=1.5, width=0.4, labelsize=5)
    cb.set_label("log$_2$ fold change", fontsize=5.5, labelpad=1.5)
    sax = fig.add_axes(rect_sz)
    sax.axis("off")
    handles = [Line2D([], [], ls="", marker="o", mfc=MID, mec=INK, mew=0.3,
                      ms=np.sqrt(6 + 9 * v)) for v in (1, 4, 8)]
    handles.append(Line2D([], [], ls="", marker="o", mfc="none", mec=RULE,
                          mew=0.5, ms=2))
    sax.legend(handles, ["1", "4", "≥8", "n.s."], ncol=4, frameon=False,
               title="−log$_{10}$ P$_{adj}$", title_fontsize=5.5,
               loc="center", handletextpad=0.2, columnspacing=0.8)


def deg_counts(ax, deg, sex):
    y = np.arange(len(COMPS))
    up = [(deg[(sex, c)].log2FoldChange > 0).sum() for c in COMPS]
    dn = [(deg[(sex, c)].log2FoldChange < 0).sum() for c in COMPS]
    ax.barh(y, up, color=UP, height=0.62, lw=0)
    ax.barh(y, [-v for v in dn], color=DOWN, height=0.62, lw=0)
    for yi, u, d in zip(y, up, dn):
        ax.text(u + 8, yi, str(u), va="center", ha="left", fontsize=5.5)
        ax.text(-d - 8, yi, str(d), va="center", ha="right", fontsize=5.5)
    ax.axvline(0, color=INK, lw=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels([LABEL[c].replace("\n", " ") for c in COMPS])
    ax.invert_yaxis()
    ax.set_xlim(-400, 200)
    ax.set_xticks([-300, -200, -100, 0, 100])
    ax.set_xticklabels(["300", "200", "100", "0", "100"])
    ax.set_xlabel("DEGs (down | up)")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_title(sex.capitalize(), loc="left")


ENERGY = [
    ("Arcuate\nneuropeptides", ["Agrp", "Npy"]),
    ("Leptin\nsignalling", ["Socs3", "Nhlh2"]),
    ("Glucocorticoid\nresponse", ["Fkbp5", "Zbtb16", "Cdkn1a", "Sult1a1",
                                  "Plin4", "Hif3a"]),
    ("Cold-inducible", ["Rbm3", "Cirbp"]),
    ("Endothelial\nactivation", ["Vwf", "Icam1", "Serpine1", "Pecam1", "Nos3",
                                 "Klf4", "Osmr", "Tgm2"]),
]

MYELIN = [
    ("OL lineage\nTFs", ["Olig1", "Olig2", "Sox10", "Myrf"]),
    ("Premyelinating\nOL", ["Bcas1", "Enpp6", "Gpr17"]),
    ("Myelin\nstructure", ["Mbp", "Plp1", "Mog", "Mag", "Mobp", "Cnp", "Cldn11",
                           "Mal", "Opalin", "Pllp", "Tmem125"]),
    ("Myelin lipid\nsynthesis", ["Ugt8a", "Fa2h", "Aspa", "Qdpr"]),
    ("Glial gap\njunctions", ["Gjc2", "Gjb1", "Gjc3"]),
    ("Astrocyte", ["Aqp4", "Il33", "Fabp7", "Car2", "Gpr37l1"]),
]
MYELIN_GENES = {g for _, gs in MYELIN[:5] for g in gs} | {
    "Trf", "Ermn", "S1pr5", "Anln", "Prr18", "Tmem278", "Hapln2", "Plekhh1"}


def figure1(deg):
    fig = plt.figure(figsize=(180 * MM, 95 * MM))
    ax_a1 = fig.add_axes([0.10, 0.60, 0.27, 0.26])
    ax_a2 = fig.add_axes([0.10, 0.17, 0.27, 0.26])
    deg_counts(ax_a1, deg, "male")
    deg_counts(ax_a2, deg, "female")
    ax_a1.set_xlabel("")
    panel_letter(ax_a1, "a", dx=-0.30, dy=1.08)

    ax_b = fig.add_axes([0.56, 0.12, 0.22, 0.78])
    cols = [("male", c) for c in COMPS]
    norm = dot_heatmap(ax_b, deg, ENERGY, cols, vlim=2.0)
    panel_letter(ax_b, "b", dx=-0.30, dy=1.02)
    size_colour_legends(fig, norm, [0.885, 0.62, 0.09, 0.022],
                        [0.855, 0.36, 0.14, 0.14])
    save(fig, "Fig1_energy_deficit")


def figure2(deg):
    fig = plt.figure(figsize=(180 * MM, 120 * MM))
    ax_a = fig.add_axes([0.10, 0.11, 0.36, 0.80])
    cols = [(s, c) for s in ("male", "female") for c in COMPS]
    norm = dot_heatmap(ax_a, deg, MYELIN, cols, vlim=1.0)
    ax_a.axvline(2.5, color=INK, lw=0.5)
    for x, s in ((1, "Male"), (4, "Female")):
        ax_a.text(x, -1.6, s, ha="center", va="bottom",
                  fontsize=6.5)
    ax_a.set_xticklabels([LABEL[c] for _, c in cols], fontsize=5)
    panel_letter(ax_a, "a", dx=-0.13, dy=1.04)

    ax_b = fig.add_axes([0.64, 0.46, 0.30, 0.45])
    s = deg[("male", "LDsaline_WT")].log2FoldChange
    l = deg[("male", "LDleptin_WT")].log2FoldChange
    shared = s.index.intersection(l.index)
    x, y = s[shared], l[shared]
    x, y = x[x.abs() < 3], y[x.abs() < 3]
    is_my = x.index.isin(MYELIN_GENES)
    lim = (-1.4, 1.4)
    ax_b.plot(lim, lim, color=RULE, lw=0.5, ls=(0, (3, 2)), zorder=1)
    ax_b.axhline(0, color=RULE, lw=0.4, zorder=1)
    ax_b.axvline(0, color=RULE, lw=0.4, zorder=1)
    ax_b.scatter(x[~is_my], y[~is_my], s=7, c="#9a9a9a", lw=0.3,
                 edgecolors="white", zorder=2)
    ax_b.scatter(x[is_my], y[is_my], s=11, c=DOWN, lw=0.3, edgecolors="white",
                 zorder=3, label="Oligodendrocyte / myelin")
    ax_b.legend(loc="lower right", frameon=False, handletextpad=0.1,
                borderaxespad=0.2)
    r = np.corrcoef(x, y)[0, 1]
    ax_b.text(0.04, 0.96, f"r = {r:.2f}",
              transform=ax_b.transAxes, va="top", fontsize=5.5)
    ax_b.set_xlim(lim)
    ax_b.set_ylim(lim)
    ax_b.set_aspect("equal")
    ax_b.set_xlabel("log$_2$ FC, LD saline vs WT")
    ax_b.set_ylabel("log$_2$ FC, LD leptin vs WT")
    panel_letter(ax_b, "b", dx=-0.17, dy=1.03)

    size_colour_legends(fig, norm, [0.68, 0.30, 0.22, 0.02],
                        [0.64, 0.14, 0.30, 0.10], vlim=1.0)
    save(fig, "Fig2_myelin_persistence")


CORE = [
    ("Oligodendrocyte", ["Olig1", "Bcas1", "Mog", "Mal", "Cldn11", "Gjc3", "Qdpr",
                         "Prr18", "Tmem278", "Plekhh1", "Adamts4"]),
    ("Astrocyte", ["Aqp4", "Il33", "Gpr37l1"]),
    ("Endothelial", ["Vwf"]),
    ("Cold-inducible", ["Rbm3"]),
]

LIPID = [
    ("Desaturation", ["Scd1", "Scd2"]),
    ("Elongation", ["Elovl5", "Elovl6", "Elovl7", "Hsd17b12"]),
    ("NADPH supply", ["Me1"]),
]


def figure3(deg):
    fig = plt.figure(figsize=(180 * MM, 90 * MM))
    ax_a = fig.add_axes([0.08, 0.14, 0.30, 0.76])
    cols_a = [(s, c) for s in ("male", "female") for c in COMPS]
    norm = dot_heatmap(ax_a, deg, CORE, cols_a, vlim=1.0)
    ax_a.axvline(2.5, color=INK, lw=0.5)
    ax_a.set_xticklabels([LABEL[c] for _, c in cols_a], fontsize=5)
    for x, sx in ((1, "Male"), (4, "Female")):
        ax_a.text(x, -1.3, sx, ha="center", va="bottom",
                  fontsize=6.5)
    panel_letter(ax_a, "a", dx=-0.2, dy=1.04)

    ax_b = fig.add_axes([0.57, 0.44, 0.31, 0.46])
    cols_b = [(s, c) for s in ("male", "female") for c in COMPS]
    dot_heatmap(ax_b, deg, LIPID, cols_b, vlim=1.0)
    ax_b.axvline(2.5, color=INK, lw=0.5)
    ax_b.set_xticklabels([LABEL[c] for _, c in cols_b], fontsize=4.6)
    for x, sx in ((1, "Male"), (4, "Female")):
        ax_b.text(x, -1.3, sx, ha="center", va="bottom",
                  fontsize=6.5)
    panel_letter(ax_b, "b", dx=-0.12, dy=1.08)

    size_colour_legends(fig, norm, [0.63, 0.20, 0.19, 0.022],
                        [0.58, 0.02, 0.30, 0.10], vlim=1.0)
    save(fig, "Fig3_sex_conserved_core_and_lipid")


def figure0():
    """Study design: three groups and the three comparisons."""
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
    fig = plt.figure(figsize=(120 * MM, 62 * MM))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 120)
    ax.set_ylim(0, 62)
    ax.axis("off")
    w, h = 30, 14
    groups = {
        "WT": (60, 46, "WT", "$\\it{Ldlr}^{-/-}$"),
        "sal": (22, 10, "LD saline", "$\\it{Ldlr}^{-/-}$; aP2-nSrebp1c-Tg"),
        "lep": (98, 10, "LD leptin", "$\\it{Ldlr}^{-/-}$; aP2-nSrebp1c-Tg"),
    }
    for x, y, name, geno in groups.values():
        ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                    boxstyle="round,pad=0,rounding_size=1.5",
                                    fc="white", ec=INK, lw=0.6))
        ax.text(x, y + 2.4, name, ha="center", va="center", fontsize=8)
        ax.text(x, y - 2.6, geno, ha="center", va="center", fontsize=6.5)

    def arrow(a, b, label, rot, off):
        (x1, y1), (x2, y2) = a, b
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="<->",
                                     mutation_scale=6, lw=0.6, color=INK,
                                     shrinkA=0, shrinkB=0))
        mx, my = (x1 + x2) / 2 + off[0], (y1 + y2) / 2 + off[1]
        ax.text(mx, my, label, ha="center", va="center", fontsize=7,
                rotation=rot, rotation_mode="anchor")

    arrow((22, 17.5), (45, 39), "Lipodystrophy effect\nLD saline vs WT",
          43, (-5.5, 4.5))
    arrow((98, 17.5), (75, 39), "Normalization\nLD leptin vs WT",
          -43, (5.5, 4.5))
    arrow((37.5, 10), (82.5, 10), "Leptin effect\nLD saline vs LD leptin",
          0, (0, -6.5))
    ax.text(60, 25, "Hypothalamus\nbulk RNA-seq\nmales and females",
            ha="center", va="center", fontsize=6.5, color=INK,
            linespacing=1.3)
    save(fig, "Fig0_study_design")


STAGES = [
    ("OPC", ["Pdgfra", "Cspg4", "Ptprz1", "Vcan", "Tnr", "Lhfpl3", "Matn4",
             "Sox6"]),
    ("COP", ["Neu4", "Bmp4", "Gpr17", "Sox4", "Sox11"]),
    ("NFOL", ["Bcas1", "Enpp6", "Tmem2", "Fyn", "Frmd4a", "Tcf7l2", "Itpr2",
              "Tns3"]),
    ("MFOL", ["Mal", "Opalin", "Plekhh1", "Ctps", "Mog", "Serinc5", "Prr18",
              "Pllp"]),
    ("MOL", ["Klk6", "Apod", "Grm3", "Slc5a11", "Ptgds", "Il33", "Anln",
             "Hapln2", "Car2", "Trf", "Cryab", "Mobp", "Aspa"]),
]
OL_TFS = [("OL transcription\nfactors", ["Olig1", "Olig2", "Sox10", "Sox8",
                                          "Myrf", "Zeb2", "Qki"])]


def stage_bars(ax, deg):
    x = np.arange(len(STAGES))
    for k, (c, col, lab) in enumerate((
            ("LDsaline_WT", DOWN, "LD saline vs WT"),
            ("LDleptin_WT", "#8FB3D9", "LD leptin vs WT"))):
        d = deg[("male", c)].log2FoldChange
        frac = [100 * sum(g in d.index and d[g] < 0 for g in gs) / len(gs)
                for _, gs in STAGES]
        ax.bar(x + (k - 0.5) * 0.36, frac, width=0.34, color=col, lw=0,
               label=lab)
    ax.set_xticks(x)
    ax.set_xticklabels([n for n, _ in STAGES])
    ax.set_ylim(0, 100)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylabel("Stage markers decreased (%)")
    ax.legend(frameon=False, loc="upper left", handlelength=1,
              borderaxespad=0)
    ax.annotate("", xy=(len(STAGES) - 0.6, -24), xytext=(-0.4, -24),
                xycoords="data", annotation_clip=False,
                arrowprops=dict(arrowstyle="->", lw=0.5, color=INK))
    ax.text((len(STAGES) - 1) / 2, -30, "Differentiation", ha="center",
            va="top", fontsize=5.5)


def figure4(deg):
    fig = plt.figure(figsize=(180 * MM, 130 * MM))
    ax_a = fig.add_axes([0.08, 0.66, 0.26, 0.26])
    stage_bars(ax_a, deg)
    ax_a.set_title("Male", loc="left", pad=4)
    panel_letter(ax_a, "a", dx=-0.22, dy=1.05)

    ax_c = fig.add_axes([0.09, 0.10, 0.29, 0.30])
    cols = [(s, c) for s in ("male", "female") for c in COMPS]
    norm = dot_heatmap(ax_c, deg, OL_TFS, cols, vlim=1.0)
    ax_c.axvline(2.5, color=INK, lw=0.5)
    ax_c.set_xticklabels([LABEL[c] for _, c in cols], fontsize=4.6)
    for x, sx in ((1, "Male"), (4, "Female")):
        ax_c.text(x, -1.2, sx, ha="center", va="bottom", fontsize=6.5)
    panel_letter(ax_c, "c", dx=-0.2, dy=1.1)

    ax_b = fig.add_axes([0.53, 0.10, 0.30, 0.84])
    dot_heatmap(ax_b, deg, STAGES, cols, vlim=1.0)
    ax_b.axvline(2.5, color=INK, lw=0.5)
    ax_b.set_xticklabels([LABEL[c] for _, c in cols], fontsize=4.6)
    for x, sx in ((1, "Male"), (4, "Female")):
        ax_b.text(x, -1.5, sx, ha="center", va="bottom", fontsize=6.5)
    panel_letter(ax_b, "b", dx=-0.16, dy=1.02)

    size_colour_legends(fig, norm, [0.08, 0.53, 0.14, 0.018],
                        [0.27, 0.49, 0.18, 0.08], vlim=1.0)
    save(fig, "Fig4_oligodendrocyte_maturation")


def save(fig, name):
    for ext in ("pdf", "svg", "png"):
        fig.savefig(OUT / f"{name}.{ext}", bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


if __name__ == "__main__":
    deg = load()
    figure0()
    figure1(deg)
    figure2(deg)
    figure3(deg)
    figure4(deg)
    print("written to", OUT)
