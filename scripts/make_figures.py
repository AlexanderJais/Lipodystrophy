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
    save(fig, "Fig3_energy_deficit")


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
    save(fig, "Fig4_myelin_persistence")


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
    save(fig, "Fig6_glial_core_lipid_synthesis")


def draw_design(ax):
    """Study design: three groups and the three comparisons (axes in mm-like units)."""
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
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
        ax.text(x, y + 2.4, name, ha="center", va="center", fontsize=7.5)
        ax.text(x, y - 2.6, geno, ha="center", va="center", fontsize=6)

    def arrow(a, b, label, rot, off):
        (x1, y1), (x2, y2) = a, b
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="<->",
                                     mutation_scale=6, lw=0.6, color=INK,
                                     shrinkA=0, shrinkB=0))
        mx, my = (x1 + x2) / 2 + off[0], (y1 + y2) / 2 + off[1]
        ax.text(mx, my, label, ha="center", va="center", fontsize=6.5,
                rotation=rot, rotation_mode="anchor")

    arrow((22, 17.5), (45, 39), "Lipodystrophy effect\nLD saline vs WT",
          43, (-5.5, 4.5))
    arrow((98, 17.5), (75, 39), "Normalization\nLD leptin vs WT",
          -43, (5.5, 4.5))
    arrow((37.5, 10), (82.5, 10), "Leptin effect\nLD saline vs LD leptin",
          0, (0, -6.5))
    ax.text(60, 25, "Hypothalamus\nbulk RNA-seq\nmales and females",
            ha="center", va="center", fontsize=6, linespacing=1.3)


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


def stage_strip(ax, deg):
    """log2FC of each differentially expressed stage marker."""
    series = (("LDsaline_WT", DOWN, "LD saline vs WT"),
              ("LDleptin_WT", "#8FB3D9", "LD leptin vs WT"),
              ("LDsaline_LDleptin", "#E3A1A9", "LD saline vs LD leptin"))
    rng = np.random.default_rng(1)
    for si, (_, gs) in enumerate(STAGES):
        jit = rng.uniform(-0.07, 0.07, len(gs))
        for k, (c, col, lab) in enumerate(series):
            d = deg[("male", c)].log2FoldChange
            xc = si + (k - 1) * 0.28
            for g, dx in zip(gs, jit):
                if g in d.index:
                    ax.scatter(xc + dx, d[g], s=9, c=col, lw=0.3,
                               edgecolors=INK, zorder=3)
    ax.axhline(0, color=INK, lw=0.5)
    ax.spines["bottom"].set_visible(False)
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", length=0, pad=3)
    ax.set_xticks(range(len(STAGES)))
    ax.set_xticklabels([n for n, _ in STAGES])
    ax.set_xlim(-0.55, len(STAGES) - 0.45)
    ax.set_ylim(-1.4, 0)
    ax.set_yticks([-1.2, -0.8, -0.4, 0])
    ax.set_yticklabels(["\u22121.2", "\u22120.8", "\u22120.4", "0"])
    ax.set_ylabel("log$_2$ fold change")
    handles = [Line2D([], [], ls="", marker="o", mfc=col, mec=INK, mew=0.3,
                      ms=3) for _, col, _ in series]
    ax.legend(handles, [lab for _, _, lab in series], frameon=False, loc="lower left", handletextpad=0.1,
              borderaxespad=0.3)
    ax.annotate("", xy=(len(STAGES) - 0.6, 0.30), xytext=(-0.4, 0.30),
                xycoords="data", annotation_clip=False,
                arrowprops=dict(arrowstyle="->", lw=0.5, color=INK))
    ax.text((len(STAGES) - 1) / 2, 0.34, "Differentiation", ha="center",
            va="bottom", fontsize=5.5)


def figure4(deg):
    fig = plt.figure(figsize=(180 * MM, 130 * MM))
    ax_a = fig.add_axes([0.08, 0.62, 0.26, 0.26])
    stage_strip(ax_a, deg)
    ax_a.text(0.0, 1.30, "Male", transform=ax_a.transAxes, fontsize=6.5,
              va="bottom")
    panel_letter(ax_a, "a", dx=-0.22, dy=1.30)

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
    save(fig, "Fig5_oligodendrocyte_maturation")


BODY = [("Weight", "Body weight (g)"), ("VAT", "VAT (mg)"), ("SAT", "SAT (mg)"),
        ("Liver", "Liver (mg)"), ("Triglycerides", "Triglycerides"),
        ("VLDL triglycerides", "VLDL triglycerides"), ("fFA", "Free fatty acids"),
        ("ALAT", "ALAT"), ("Glucose", "Glucose")]
LIPIDS = [("VLDL triglycerides", "VLDL triglycerides"),
          ("LDL triglycerides", "LDL triglycerides"),
          ("HDL cholesterol", "HDL cholesterol")]
GROUPS = [("WT", INK), ("LD saline", UP), ("LD leptin", DOWN)]


def stars(p):
    return ("****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 0.01
            else "*" if p < 0.05 else "")


def phenotype_panels(axes, measures, letters, week=14):
    """Per-mouse dot plots, females and males side by side, week of the RNA-seq cohort."""
    from scipy.stats import mannwhitneyu
    ph = pd.read_csv(DATA / "phenotype_extracted.csv")
    ph = ph[ph.week == week]
    rng = np.random.default_rng(3)
    for ax, (m, ylab), letter in zip(axes, measures, letters):
        panel_letter(ax, letter, dx=-0.25, dy=1.0)
        top = ph[ph.measure == m].value.max()
        for si, sex in enumerate(("female", "male")):
            vals = {}
            for gi, (g, col) in enumerate(GROUPS):
                v = ph[(ph.measure == m) & (ph.sex == sex) & (ph.group == g)].value
                vals[g] = v.values
                x = si * 3.6 + gi
                ax.scatter(x + rng.uniform(-0.18, 0.18, len(v)), v, s=7, c=col,
                           lw=0, alpha=0.85, zorder=3)
                ax.plot([x - 0.3, x + 0.3], [v.median()] * 2, color=INK, lw=0.8,
                        zorder=4)
            if not len(vals["WT"]):
                ax.text(si * 3.6 + 1, top * 0.05, "n.d.", ha="center",
                        va="bottom", fontsize=5.5)
                continue
            for gi, (g, _) in enumerate(GROUPS[1:], start=1):
                mk = stars(mannwhitneyu(vals[g], vals["WT"]).pvalue)
                if g == "LD leptin":
                    q = stars(mannwhitneyu(vals[g], vals["LD saline"]).pvalue)
                    mk += ("\n" if mk else "") + "#" * len(q)
                ax.text(si * 3.6 + gi, top * 1.08, mk, ha="center", va="bottom",
                        fontsize=5, linespacing=0.9)
        ax.set_ylim(0, top * 1.28)
        ax.set_xlim(-0.6, 6.2)
        ax.set_xticks([1, 4.6])
        ax.set_xticklabels(["Female", "Male"])
        ax.tick_params(axis="x", length=0)
        ax.set_ylabel(ylab)


def group_legend(fig, y):
    handles = [Line2D([], [], ls="", marker="o", mfc=c, mec="none", ms=3.5)
               for _, c in GROUPS]
    fig.legend(handles, [g for g, _ in GROUPS], ncol=3, frameon=False,
               loc="upper center", bbox_to_anchor=(0.5, y), handletextpad=0.1)


def fig1_design_adipose():
    fig = plt.figure(figsize=(180 * MM, 112 * MM))
    ax_d = fig.add_axes([0.2, 0.47, 0.6, 0.53])
    draw_design(ax_d)
    panel_letter(ax_d, "a", dx=-0.25, dy=0.95)
    axes = [fig.add_axes([0.08 + i * 0.33, 0.05, 0.24, 0.32]) for i in range(3)]
    phenotype_panels(axes, BODY[:3], "bcd")
    group_legend(fig, 0.47)
    save(fig, "Fig1_design_adipose")


def fig2_metabolic():
    fig, axes = plt.subplots(2, 3, figsize=(180 * MM, 110 * MM))
    fig.subplots_adjust(left=0.08, right=0.98, top=0.89, bottom=0.07,
                        wspace=0.45, hspace=0.32)
    measures = [BODY[3], BODY[4], BODY[7]] + LIPIDS
    phenotype_panels(axes.flat, measures, "abcdef")
    group_legend(fig, 1.0)
    save(fig, "Fig2_liver_lipoproteins")


def figS1_ffa_glucose():
    fig, axes = plt.subplots(1, 3, figsize=(180 * MM, 62 * MM))
    fig.subplots_adjust(left=0.08, right=0.98, top=0.80, bottom=0.12,
                        wspace=0.45)
    phenotype_panels(axes[:2], [BODY[6], BODY[8]], "ab")
    axes[2].axis("off")
    group_legend(fig, 1.0)
    save(fig, "FigS1_fatty_acids_glucose")


PART_COLS = [("Normalized_up", "Normalized\nup", UP),
             ("Normalized_down", "Normalized\ndown", DOWN),
             ("Persistent_down", "Persistent\ndown", DOWN)]
CELLTYPES = [("Oligodendrocytes", "Oligodendrocytes"),
             ("Oligodendrocyte Progenitor Cells", "OPCs"),
             ("Astrocytes", "Astrocytes"), ("Neurons", "Neurons"),
             ("Endothelial Cells", "Endothelial cells"),
             ("Pericytes", "Pericytes"), ("Microglia", "Microglia")]
TFS = [("NR3C1 34362910 ChIP-Seq WistarRat Hippocampus Stress", "NR3C1 (GR)"),
       ("OLIG2 23332759 ChIP-Seq OLIGODENDROCYTES Mouse", "OLIG2"),
       ("SMARCA4 23332759 ChIP-Seq OLIGODENDROCYTES Mouse", "SMARCA4 (BRG1)"),
       ("SOX10 human", "SOX10")]


def enrich_dots(ax, tables, rows):
    for j, (key, _, col) in enumerate(PART_COLS):
        t = tables[key]
        for i, (term, _) in enumerate(rows):
            r = t[t.term == term]
            q = r.padj.iloc[0] if len(r) else 1.0
            if q < 0.05:
                ax.scatter(j, i, s=6 + 5 * min(-np.log10(q), 20), c=col,
                           lw=0.3, edgecolors=INK, zorder=3)
            else:
                ax.scatter(j, i, s=4, facecolors="none", edgecolors=RULE,
                           lw=0.5, zorder=2)
    ax.set_xlim(-0.6, len(PART_COLS) - 0.4)
    ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([lab for _, lab in rows])
    ax.set_xticks([])
    ax.tick_params(length=0, pad=3)
    yb = len(rows) - 0.5 + 0.25
    ax.plot([-0.25, 1.25], [yb, yb], color=INK, lw=0.5, clip_on=False)
    ax.plot([1.75, 2.25], [yb, yb], color=INK, lw=0.5, clip_on=False)
    ax.text(0.5, yb + 0.2, "Normalized\nby leptin", ha="center", va="top",
            fontsize=5.5)
    ax.text(2, yb + 0.2, "No leptin\neffect", ha="center", va="top",
            fontsize=5.5)
    for sp in ax.spines.values():
        sp.set_visible(False)


def figure7():
    part = ROOT / "results" / "leptin_partition"
    tables = {k: pd.read_csv(part / f"{k}.csv") for k, _, _ in PART_COLS}
    sizes = pd.read_csv(part / "sizes.csv", index_col=0).iloc[:, 0]
    fig = plt.figure(figsize=(180 * MM, 70 * MM))

    ax_a = fig.add_axes([0.07, 0.18, 0.22, 0.66])
    groups = [("Normalized\nby leptin", sizes["Normalized up"], sizes["Normalized down"]),
              ("No leptin\neffect", sizes["Persistent up"], sizes["Persistent down"])]
    y = np.arange(len(groups))
    for yi, (_, up, dn) in zip(y, groups):
        ax_a.barh(yi, up, color=UP, height=0.55, lw=0)
        ax_a.barh(yi, dn, left=up, color=DOWN, height=0.55, lw=0)
        ax_a.text(up + dn + 6, yi, str(up + dn), va="center", fontsize=5.5)
    ax_a.set_yticks(y)
    ax_a.set_yticklabels([g for g, _, _ in groups])
    ax_a.invert_yaxis()
    ax_a.set_xlim(0, 330)
    ax_a.set_xlabel("Male LD saline vs WT DEGs")
    ax_a.spines["left"].set_visible(False)
    ax_a.tick_params(axis="y", length=0)
    panel_letter(ax_a, "a", dx=-0.42, dy=1.06)

    ax_b = fig.add_axes([0.45, 0.18, 0.21, 0.66])
    enrich_dots(ax_b, tables, CELLTYPES)
    ax_b.set_title("Cell-type markers", fontsize=6, pad=4)
    panel_letter(ax_b, "b", dx=-0.62, dy=1.06)

    ax_c = fig.add_axes([0.78, 0.40, 0.21, 0.44])
    enrich_dots(ax_c, tables, TFS)
    ax_c.set_title("Transcription-factor targets", fontsize=6, pad=4)
    panel_letter(ax_c, "c", dx=-0.62, dy=1.1)

    sax = fig.add_axes([0.66, 0.02, 0.32, 0.14])
    sax.axis("off")
    handles = [Line2D([], [], ls="", marker="o", mfc="white", mec=INK, mew=0.3,
                      ms=np.sqrt(6 + 5 * v)) for v in (2, 10, 20)]
    handles.append(Line2D([], [], ls="", marker="o", mfc="none", mec=RULE,
                          mew=0.5, ms=2))
    sax.legend(handles, ["2", "10", "\u226520", "n.s."], ncol=4, frameon=False,
               title="\u2212log$_{10}$ q", title_fontsize=5.5, loc="center",
               handletextpad=0.2, columnspacing=0.8)
    save(fig, "Fig7_leptin_fixes_vs_persists")


MODULE_NAMES = {"Sox10": "OL transcription &\nmyelin core",
                "Mog": "Myelin membrane", "Nefl": "Axonal cytoskeleton",
                "S1pr5": "GPCR signalling", "Gad1": "GABAergic",
                "Gjc2": "Gap junctions", "Bcas1": "Premyelinating OL",
                "Ngfr": "Neurotrophin", "P2ry12": "Microglial\nhomeostasis",
                "Elovl6": "Fatty-acid\nsynthesis", "Plcb1": "Phospholipase C",
                "Creb1": "CREB / PGC-1\u03b1", "Pax6": "Pax6", "Fgfr2": "FGF / VEGF"}


def network_modules(g, min_size=3):
    import networkx as nx
    mods = []
    for comp in nx.connected_components(g):
        if len(comp) < min_size:
            continue
        sub = g.subgraph(comp)
        if len(comp) > 8:
            mods += [set(c) for c in nx.community.greedy_modularity_communities(sub)]
        else:
            mods.append(set(comp))
    return sorted(mods, key=len, reverse=True)


def draw_network(ax, name, color_col, ncol, vlim=1.0, fs=7):
    import networkx as nx
    net = ROOT / "results" / "network"
    e = pd.read_csv(net / f"{name}_edges.csv")
    nodes = pd.read_csv(net / f"{name}_nodes.csv").set_index("gene")
    g = nx.from_pandas_edgelist(e, "preferredName_A", "preferredName_B")
    mods = network_modules(g)
    keep = set().union(*mods)
    pos = {}
    for k, m in enumerate(mods):
        cx, cy = (k % ncol) * 1.45, -(k // ncol) * 1.45
        sub = g.subgraph(m)
        if len(m) > 9:
            order = sorted(m, key=lambda n: -sub.degree(n))
            p = nx.circular_layout(sub.subgraph(order))
        elif len(m) > 2:
            p = nx.kamada_kawai_layout(sub)
        else:
            p = nx.circular_layout(sub)
        r = 0.52 if len(m) > 9 else 0.44 if len(m) > 4 else 0.30
        for n, (x, y) in p.items():
            pos[n] = (cx + x * r, cy + y * r)
        label = next((MODULE_NAMES[n] for n in m if n in MODULE_NAMES), "")
        ax.text(cx, cy + r + 0.16, label, ha="center", va="bottom",
                fontsize=fs + 1)
    norm = TwoSlopeNorm(0, -vlim, vlim)
    for a, b in g.subgraph(keep).edges:
        ax.plot(*zip(pos[a], pos[b]), color=RULE, lw=0.5, zorder=1)
    for n in keep:
        x, y = pos[n]
        v = nodes.at[n, color_col]
        deg = g.degree(n)
        ax.scatter(x, y, s=25 + 9 * deg, c=[CMAP(norm(v)) if pd.notna(v) else "white"],
                   edgecolors=INK, lw=0.4, zorder=3)
        ax.text(x, y - 0.075, n, ha="center", va="top", fontsize=fs,
                fontstyle="italic", zorder=4)
    ax.set_aspect("equal")
    ax.axis("off")
    return norm


def figure8():
    net = ROOT / "results" / "network"
    summ = pd.read_csv(net / "summary.csv").set_index("set")
    fig = plt.figure(figsize=(250 * MM, 150 * MM))
    ax_a = fig.add_axes([0.0, 0.03, 0.60, 0.88])
    norm = draw_network(ax_a, "persistent_down", "lfc_leptin_vs_wt", ncol=3)
    fig.text(0.30, 0.995, "No leptin effect", ha="center", va="top", fontsize=10)
    fig.text(0.005, 0.995, "a", ha="left", va="top", fontsize=12, fontweight="bold")

    ax_b = fig.add_axes([0.63, 0.42, 0.36, 0.50])
    draw_network(ax_b, "normalized_down", "lfc_saline_vs_wt", ncol=3)
    fig.text(0.81, 0.995, "Normalized by leptin", ha="center", va="top", fontsize=10)
    fig.text(0.625, 0.995, "b", ha="left", va="top", fontsize=12, fontweight="bold")

    ax_c = fig.add_axes([0.70, 0.08, 0.14, 0.24])
    order = [("normalized_down", "Normalized\nby leptin"),
             ("persistent_down", "No leptin\neffect")]
    x = np.arange(len(order))
    obs = [summ.at[k, "edges"] for k, _ in order]
    exp = [summ.at[k, "expected_edges"] for k, _ in order]
    ax_c.bar(x - 0.18, exp, width=0.34, color=RULE, lw=0, label="Expected")
    ax_c.bar(x + 0.18, obs, width=0.34, color=DOWN, lw=0, label="Observed")
    ax_c.set_xticks(x)
    ax_c.set_xticklabels([l for _, l in order], fontsize=8)
    ax_c.tick_params(axis="x", length=0)
    ax_c.tick_params(axis="y", labelsize=8)
    ax_c.set_ylabel("Interactions", fontsize=9)
    ax_c.legend(frameon=False, loc="upper left", handlelength=1, borderaxespad=0,
                fontsize=8)
    fig.text(0.625, 0.36, "c", ha="left", va="top", fontsize=12, fontweight="bold")

    cax = fig.add_axes([0.87, 0.14, 0.10, 0.022])
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=CMAP), cax=cax,
                      orientation="horizontal", ticks=[-1, 0, 1])
    cb.outline.set_linewidth(0.4)
    cb.ax.tick_params(length=2, width=0.4, labelsize=8)
    cb.set_label("log$_2$ fold change", fontsize=9, labelpad=2)
    save(fig, "Fig8_network")

def save(fig, name):
    for ext in ("pdf", "svg", "png"):
        fig.savefig(OUT / f"{name}.{ext}", bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


if __name__ == "__main__":
    deg = load()
    fig1_design_adipose()
    fig2_metabolic()
    figure1(deg)   # Fig. 3
    figure2(deg)   # Fig. 4
    figure4(deg)   # Fig. 5
    figure3(deg)   # Fig. 6
    figure7()      # Fig. 7
    figure8()      # Fig. 8
    figS1_ffa_glucose()
    print("written to", OUT)
