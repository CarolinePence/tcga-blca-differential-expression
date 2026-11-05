# TCGA-BLCA Differential Expression by Smoking History

Differential gene expression analysis of stage IV bladder urothelial carcinoma (BLCA), comparing male lifelong non-smokers to male patients who had quit smoking more than 15 years before diagnosis. Built in R using TCGA data from the NCI Genomic Data Commons.

**Question:** Do tumors from long-term former smokers show gene expression differences from tumors in patients who never smoked?

------------------------------------------------------------------------

## Data

Source: [TCGA-BLCA](https://portal.gdc.cancer.gov/projects/TCGA-BLCA) project, accessed through the GDC Data Portal.

Cases were filtered to male patients with stage IV transitional cell carcinoma originating in the bladder, with open-access gene expression data (`STAR - Counts`) available. They were then split by tobacco smoking history into two groups:

| Group         | Description                          | n   |
|---------------|--------------------------------------|-----|
| 1 (reference) | Lifelong non-smokers                 | 11  |
| 2             | Current reformed smokers, \>15 years | 23  |

A positive log2 fold change therefore indicates **higher** expression in the former-smoker group.

Sample sheets exported from the GDC are included (`BLCA_non_smokers.csv`, `BLCA_reformed_smokers.csv`). Count data is not stored in this repo — the pipeline downloads it directly from the GDC via `TCGAbiolinks`.

------------------------------------------------------------------------

## Approach

1.  **Retrieval** — Queried the GDC with `TCGAbiolinks`, filtered results to the target sample barcodes, and structured the counts into a `SummarizedExperiment` object.
2.  **Sample QC** — Ran principal component analysis on variance-filtered normalized counts to screen for outlier samples. PC1 and PC2 explained 9.54% and 7.91% of variance; no samples were extreme enough to warrant removal.
3.  **Locus filtering** — Dropped loci where 90% or more of the samples in either group had zero counts, since these produce unstable fold-change estimates. 38,979 loci remained.
4.  **Differential expression** — Fit a negative binomial model with `DESeq2` using a single grouping factor, with non-smokers as the reference level. Significance threshold: adjusted p \< 0.05 (FDR).
5.  **Visualization** — Built PCA and volcano plots with `ggplot2` and `ggrepel`.

------------------------------------------------------------------------

## Results

Of 38,979 loci analyzed, 27,643 had non-missing adjusted p-values. At an adjusted p-value cutoff of 0.05:

- **274 differentially expressed genes** total
- 192 upregulated in former smokers
- 82 downregulated in former smokers

Top 10 by adjusted p-value: `SFTA2`, `SFTPA2`, `PTHLH`, `PEX5L`, `STXBP6`, `LINC02432`, `HEPHL1`, `FER1L6`, `SIGLEC5`, `PRSS21`.

### Interpretation

The two strongest hits, **SFTA2** and **SFTPA2**, both encode pulmonary surfactant-associated proteins normally expressed in lung tissue, and both were elevated in the former-smoker tumors. A plausible reading is not that bladder tissue becomes "lung-like," but that long-term smoke exposure altered how genes normally silenced outside the lung are regulated in the urothelium.

**PTHLH**, which encodes parathyroid hormone-related protein and is generally associated with tumor progression, was *lower* in former smokers — the opposite of the naive expectation. Prior work has reported that some smoking-related expression changes persist irreversibly in airway epithelium after quitting, so one possibility is a durable downregulation that outlasted smoking cessation by 15+ years. Another is that bladder cancers in smokers and non-smokers progress through partly different pathways.

### Limitations

- PC1 and PC2 together explained under 20% of variance, so most sample-to-sample variation is driven by something other than smoking history — likely tumor subtype and patient biology.
- Group 1 is small (n = 11), which limits statistical power.
- Restricted to stage IV, so effects of late-stage disease cannot be separated from smoking effects.
- Differential expression detects association, not causation.

Natural extensions: add a current-smoker group to compare across never / recently quit / long-quit, include matched non-tumor tissue as a control, and run pathway-level enrichment across all 274 genes rather than interpreting top hits individually.

------------------------------------------------------------------------

## Files

| File | Description |
|----|----|
| `pipeline.R` | Main analysis script — GDC query, download, QC filtering, and DESeq2 model |
| `functions.R` | Plotting and summary functions: `nicePCA()`, `niceVolcano()`, `resSummary()` |
| `BLCA_non_smokers.csv` | GDC sample sheet, group 1 |
| `BLCA_reformed_smokers.csv` | GDC sample sheet, group 2 |
| `Samples.csv` | Clinical metadata for the 34 samples retained after filtering |
| `project_paper.pdf` | Full write-up with figures and references |

------------------------------------------------------------------------

## Running it

Requires R with `TCGAbiolinks`, `DESeq2`, `SummarizedExperiment`, `ggplot2`, and `ggrepel`.

``` r
source("functions.R")
source("pipeline.R")
```

The pipeline downloads roughly 34 count files from the GDC on first run, which takes a few minutes.

------------------------------------------------------------------------

## Notes

The analysis workflow was adapted from a template provided in INB 321G (Spring 2026) at UT Austin. The experimental design, sample selection, plotting functions in `functions.R`, and all interpretations are my own.

Results are based in part on data generated by the TCGA Research Network: <https://www.cancer.gov/tcga>
