# DNMT3A-HOX-AML-multi-Omics-
Multi-omics analysis of DNMT3A-mutant AML reveals dual HOX epigenetic reprogramming.
# DNMT3A-HOX-AML

## Overview

This project investigates the impact of DNMT3A mutations on HOX gene regulation in Acute Myeloid Leukemia (AML) using TCGA-LAML multi-omics data.

## Main Findings

- HOXA9 and HOXA10 are activated in DNMT3A-mutant AML.
- HOXC/HOXD clusters are repressed.
- Evidence supports a dual HOX epigenetic switch.
- Potential biomarkers: HOXA9, HOXA10, HOXB4, HOXB13.

## Data Source

- TCGA-LAML
- GDC Data Portal

## Pipeline

1. Mutation analysis
2. Differential expression analysis
3. HOX cluster analysis
4. DNA methylation analysis
5. miRNA integration
6. Survival analysis

## Author

AmirHossein Yari


# 📄 scripts/01_load_data.R

```r
library(data.table)

# Expression
expr <- fread("../data/TCGA-LAML.star_counts.tsv.gz")

# Mutation
mut <- fread("../data/TCGA-LAML.somaticmutation_wxs.tsv.gz")

# Clinical
clinical <- fread("../data/TCGA-LAML.clinical.tsv.gz")

# Survival
survival <- fread("../data/TCGA-LAML.survival.tsv.gz")

save(
  expr,
  mut,
  clinical,
  survival,
  file = "../data/raw_data.RData"
)
```

---

# 📄 scripts/02_DNMT3A_grouping.R

```r
load("../data/raw_data.RData")

dnmt3a <- mut[gene == "DNMT3A"]

mut.samples <- unique(
  dnmt3a$sample
)

expr.samples <- colnames(expr)[-1]

group <- ifelse(
  expr.samples %in% mut.samples,
  "Mutant",
  "WT"
)

group.df <- data.frame(
  Sample = expr.samples,
  Group = group
)

write.csv(
  group.df,
  "../results/DNMT3A_groups.csv",
  row.names = FALSE
)
```

---

# 📄 scripts/03_DEG_analysis.R

```r
library(DESeq2)

load("../data/raw_data.RData")

group.df <- read.csv(
  "../results/DNMT3A_groups.csv"
)

counts <- as.data.frame(expr)

rownames(counts) <- counts[,1]
counts <- counts[,-1]

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),
  colData = group.df,
  design = ~ Group
)

dds <- DESeq(dds)

res <- results(dds)

deg <- as.data.frame(res)

deg$Gene <- rownames(deg)

write.csv(
  deg,
  "../results/All_DEGs.csv",
  row.names = FALSE
)
```

---

# 📄 scripts/04_HOX_analysis.R

```r
deg <- read.csv(
  "../results/All_DEGs.csv"
)

HOX <- c(
"HOXA1","HOXA2","HOXA3","HOXA4","HOXA5",
"HOXA6","HOXA7","HOXA9","HOXA10","HOXA11",
"HOXB1","HOXB2","HOXB3","HOXB4","HOXB5",
"HOXB6","HOXB7","HOXB8","HOXB9","HOXB13",
"HOXC4","HOXC5","HOXC8","HOXC9","HOXC10",
"HOXC11","HOXC12","HOXC13",
"HOXD3","HOXD4","HOXD8","HOXD9",
"HOXD10","HOXD11","HOXD12","HOXD13"
)

hox <- subset(
  deg,
  Gene %in% HOX
)

write.csv(
  hox,
  "../results/HOX_DEGs.csv",
  row.names = FALSE
)
```

---

# 📄 scripts/05_volcano_plot.R

```r
library(EnhancedVolcano)

deg <- read.csv(
  "../results/All_DEGs.csv"
)

pdf(
  "../figures/Figure3_Volcano.pdf",
  width = 8,
  height = 7
)

EnhancedVolcano(
  deg,
  lab = deg$Gene,
  x = "log2FoldChange",
  y = "pvalue",
  pCutoff = 0.05,
  FCcutoff = 1
)

dev.off()
```

---

# 📄 scripts/06_HOX_heatmap.R

```r
library(pheatmap)

load("../data/raw_data.RData")

hox <- read.csv(
  "../results/HOX_DEGs.csv"
)

genes <- hox$Gene

mat <- counts[
  rownames(counts) %in% genes,
]

pdf(
  "../figures/Figure4_HOX_Heatmap.pdf",
  width = 10,
  height = 8
)

pheatmap(
  mat,
  scale = "row"
)

dev.off()
```

---

# 📄 scripts/07_survival_analysis.R

```r
library(survival)
library(survminer)

load("../data/raw_data.RData")

fit <- survfit(
  Surv(OS.time, OS) ~ Group,
  data = survival
)

pdf(
  "../figures/Figure8_Survival.pdf",
  width = 7,
  height = 6
)

ggsurvplot(
  fit,
  pval = TRUE
)

dev.off()
```

---

# 📄 scripts/08_export_supplementary_tables.R

```r
library(openxlsx)

deg <- read.csv(
  "../results/All_DEGs.csv"
)

hox <- read.csv(
  "../results/HOX_DEGs.csv"
)

wb <- createWorkbook()

addWorksheet(wb,"S1_All_DEGs")
writeData(wb,"S1_All_DEGs",deg)

addWorksheet(wb,"S2_HOX")
writeData(wb,"S2_HOX",hox)

saveWorkbook(
  wb,
  "../results/Supplementary_Tables.xlsx",
  overwrite = TRUE
)
```

---

# 📄 scripts/run_all.R


```r
source("01_load_data.R")

source("02_DNMT3A_grouping.R")

source("03_DEG_analysis.R")

source("04_HOX_analysis.R")

source("05_volcano_plot.R")

source("06_HOX_heatmap.R")

source("07_survival_analysis.R")

source("08_export_supplementary_tables.R")
```

---

