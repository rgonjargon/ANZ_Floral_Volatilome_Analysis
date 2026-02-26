# Flore VOC

NZ Floral Volatilome Analysis — analysis of floral volatile organic compounds (VOCs) from New Zealand plants (Natives, Crops, Ornamentals).

## Repo

- **Scripts/** — Quarto report `1_analysis.qmd`; render produces `2_voc_analysis.html`.
- **Analysis/Data/** — Input data (Excel); paths used in the report.

## Run the report

From the project root in R or a terminal:

```bash
quarto render Scripts/1_analysis.qmd
```

Output: `Scripts/2_voc_analysis.html`.

## Requirements

R with packages listed in the setup chunk of `1_analysis.qmd` (tidyverse, readxl, here, vegan, janitor, tidytext, ggvenn, reshape2, ape, ggtree, patchwork, RColorBrewer, etc.). Quarto: [quarto.org](https://quarto.org).

## Authors

Tom Moore & Flore Mas
