# HackBio Stage 2 – R Implementation 🔬

## Overview 📌  
Exploratory data analysis and visualization of gene expression, breast cancer, and immune response datasets in **R**, for the HackBio Internship Stage Two tasks.

---

## Datasets 📁  
- **HBR/UHR Normalized Counts** – Normalized counts for top differentially expressed genes.  
- **HBR/UHR DEG Chr22** – Differential expression results for chromosome 22 with significance labels.  
- **Breast Cancer Wisconsin** – Diagnostic features for breast cancer diagnosis.  
- **HackBio Stage 2 Excel (`hb_stage_2.xlsx`)** – Multi-sheet kinetics data.  

---

## Tasks & Methods 🧮  

### Part 1 – Gene Expression Analysis  
- **1a – Heatmap:**  
  - `pheatmap` on HBR/UHR normalized counts  
  - Row/column clustering, custom row labels, blue palette  

- **1b – Volcano plot:**  
  - Base R scatter of `log2FoldChange` vs `log10(Padj)`  
  - Points colored by significance, dashed cut-off lines, legend  

---

### Part 2 – Breast Cancer Exploration 🩺  
- **2c – Scatter:** `radius_mean` vs `texture_mean` with `ggplot2`, colored by diagnosis and fixed axes.  
- **2d – Correlation heatmap:**  
  - Correlation of 6 mean features  
- **2e – Scatter:** `smoothness_mean` vs `compactness_mean` using `ggplot2`, colored by diagnosis.  
- **2f – Density:** `area_mean` distributions by diagnosis using `geom_density`.  

---

### Part 3 – Reproduced Visual images 🧫  

- **2a (3a) – Boxplots:** Cell-type ratio distributions by `cell_type` with a custom qualitative palette.  
- **2b (3b) – Log2 scatter:**  
  - `log2(half_life)` vs `log2(alpha)`  
  - Regimes defined with `dplyr::case_when`, custom colors, quadrant lines, labels for selected genes.  

- **2c (3c) – Heatmap (genes × cell type × time):**  
  - `pheatmap` with gene clustering only  
  - Column annotations for cell type and time  

- **2d (3d) – Pathway enrichment map:**  
  - `pivot_longer` to long format  
  - Tile heatmap with fixed pathway order and a diverging red–white–blue palette.  

- **2e (3e) – Bubble plot:**  
  - `half_life` vs `alpha`  
  - Point size = `count`, color = `stage`, classic theme.  

- **2f (3f) – Stacked bar:**  
  - Subset to B and Plasma cells  
  - Stacked proportions by stage with manual colors.  

- **2g – Directed network:**  
  - Adjacency matrix → directed, weighted graph in `igraph`  
  - Self-loops removed, zero-weight edges dropped  
  - Custom node/edge styling to highlight interaction structure.  

---

## Requirements ⚙️  

Install needed packages:

```r
install.packages(c(
  "tidyr", "dplR", "pheatmap", "reshape",
  "readxl", "ggplot2", "igraph", "dplyr", "png"
))
