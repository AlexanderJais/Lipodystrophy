# Preliminary cross-contrast findings (hypothalamus RNA-seq)

Inputs: six DESeq2 DEG lists (`data/`), pre-filtered at padj < 0.1 and |log2FC| > 0.2.
Reproduce with `python3 scripts/cross_contrast.py` → tables in `results/`.

| Contrast | Male DEGs | Female DEGs |
|---|---|---|
| LDsaline vs WT | 456 (137 up / 319 down) | 62 |
| LDleptin vs WT | 297 | 4 |
| LDsaline vs LDleptin | 11 | 35 |

## What the data support
- **Arcuate energy-deficit signature (males only):** Agrp +1.88, Npy +1.06 in male LDsaline vs WT;
  both absent from male LDleptin vs WT (consistent with normalisation). Not detected in females.
- **Glucocorticoid-response signature (males):** Fkbp5, Zbtb16, Cdkn1a, Sult1a1, Hif3a, Plin4,
  Serpine1 all up ~2-fold in LDsaline, gone with leptin. Not mentioned in the outline; worth adding.
- **Oligodendrocyte/myelin down-regulation:** strong and coordinated in males (Mbp, Plp1, Mog, Mag,
  Mobp, Cnp, Myrf, Sox10, Olig1/2, Ugt8a, Fa2h, Gjc2, Enpp6, Opalin …, log2FC −0.4 to −1.0);
  weaker in females (Mog, Mal, Bcas1, Cldn11, Olig1, Hapln2, Qdpr, Gjc3; −0.2 to −0.36).

## Where the data contradict the draft
1. **Myelin loss is NOT reversed by metreleptin in males.** 170/456 male LDsaline-vs-WT DEGs remain
   DE in LDleptin vs WT, all 170 in the same direction (r = 0.98), including essentially the whole
   myelin program. No myelin gene is in male LDsaline vs LDleptin. The myelin phenotype tracks
   the genotype, not leptin status, in males.
2. **In females, "reversal" is indirect.** LDleptin vs WT has only 4 DEGs, but no myelin gene reaches
   significance in LDsaline vs LDleptin. The only directly leptin-rescued genes are Vgf, Podxl, Plxnd1,
   Nr2f2, Htra1, Plch2, Prkcg, Cacnb3, C4b, Dlk1, Fat4. Failing to find a difference from WT is not
   evidence of reversal.
3. **Ch25h / 25-HC / NLRP3 axis: no support.** Ch25h, Cyp7b1, Nlrp3, Pycard, Casp1, Il1b, Il18,
   Gsdmd absent from all six lists. Only Abca1 (+0.27, male saline) appears.
4. Leptin-responsive genes in male saline vs leptin are vascular/inflammatory
   (Lcn2, Icam1, Tgm2, Anxa3, Sgk3), not myelin.

## Technical flags
- **Dissection variability:** Kcnj13 (choroid plexus) ±2.2–2.4; Lhx8 −3.6 and Vip −2.9 (male
  leptin vs WT); midbrain dopaminergic markers Th, Ddc, Slc18a2, Aldh1a1, Ret down in female
  saline vs WT. These point to differing region boundaries between groups, which can also move
  oligodendrocyte/white-matter content. Myelin genes should be normalised against a white-matter
  contamination check (e.g. fraction of reads, deconvolution) before being interpreted.
- **Tcaf2** +3 log2FC in LD vs WT in both treatment groups: likely a transgene/strain-linked
  passenger gene, not biology.
- **Mir6236 (−12.8), Lars2 (−6.8):** known mapping artefacts; exclude.
- Socs3 +0.99 in male LDleptin vs WT: positive control for leptin exposure.

## Needed next
- **Full, unfiltered DESeq2 tables** (all genes, with baseMean/lfcSE) to test Agrp/Npy in females,
  Ch25h/Nlrp3 below threshold, and do GSEA on the myelin set.
- Normalised counts per sample + sample sheet (sex, group, batch, dissection date) for a
  sex × genotype × treatment model and interaction test for leptin rescue.

## Data-supported narrative (draft)
1. **Leptin-reversible hypothalamic "starvation state" (male-dominant).** In saline-treated LD mice:
   Agrp/Npy up; glucocorticoid-target genes up (Fkbp5, Zbtb16, Cdkn1a, Sult1a1, Plin4, Hif3a);
   cold-inducible Rbm3/Cirbp up (Rbm3 also in females), consistent with hypothermia; endothelial
   activation (Vwf, Icam1, Serpine1, Pecam1, Nos3, Klf4, Osmr). None persist in LDleptin vs WT;
   Icam1, Tgm2, Anxa3, Sgk3 are directly lower with leptin.
2. **Leptin-independent glial deficit.** Oligodendrocyte/myelin genes plus astrocyte genes
   (Il33, Aqp4, Fabp7, Car2, Gpr37l1) are reduced; in males this persists unchanged under
   metreleptin (r = 0.98). Females show the same direction, milder, and without direct rescue.
3. **Target engagement.** Socs3 and Nhlh2 (leptin-induced arcuate genes) up in leptin-treated males.
4. **Sex difference.** Males 456 DEGs vs females 62.
