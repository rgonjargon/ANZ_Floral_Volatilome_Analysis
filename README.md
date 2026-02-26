# Flore VOC

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18793708.svg)](https://doi.org/10.5281/zenodo.18793708)

Analysis of floral volatile organic compounds (FVOCs) from New Zealand plants—natives, crops, and ornamentals—for *The chemical diversity of floral volatilomes in native and exotic plants shaping Aotearoa-New Zealand's odourscapes, implications for plant-insect interactions*.

**Cite:** [https://doi.org/10.5281/zenodo.18793708](https://doi.org/10.5281/zenodo.18793708)

## Project structure

| Folder / file | Contents |
|---------------|----------|
| **Scripts/** | `1_analysis.qmd` — main Quarto report; renders to `2_voc_analysis.html` |
| **Analysis/Data/** | Input data (Excel/CSV) used by the report |
| **Plots/** | Figures (Figures 1–7, ESM 1–7) saved by the report |

## How to run

From the project root:

```bash
quarto render Scripts/1_analysis.qmd
```

Output: `Scripts/2_voc_analysis.html`.

## Requirements

- [R](https://www.r-project.org/) and [Quarto](https://quarto.org/)
- R packages: see the setup chunk in `1_analysis.qmd` (tidyverse, readxl, here, vegan, janitor, tidytext, ggvenn, reshape2, ape, ggtree, patchwork, RColorBrewer, ggvegan, V.PhyloMaker2, etc.)

## Authors

Tom Moore & Flore Mas
