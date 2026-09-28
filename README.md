# circular_ST

Code and processed data for **Circular Data Analysis for Spatial Omics**
(Shin, Yoo, Cho, et al.). The repository reproduces the two case studies in the
manuscript and accompanies the R examples in the Supplementary Note.

## Contents

```
circular_ST/
├── supplementary_core_analysis.R   # main analysis script
├── simulate_toy_data.R             # optional synthetic inputs with the same columns
└── data/
    ├── DESeq2_region_age_gene_phase_all_genes.csv
    ├── xenium_cell_metadata_with_pathology_and_official_10x_annotation.csv.gz
    └── gp_output/test_with_pred.csv.gz    # archived projected-GP predictions
```

## Data

| File | Dataset | Circular quantity | Main columns |
|---|---|---|---|
| `DESeq2_region_age_gene_phase_all_genes.csv` | AD mouse-brain spatial transcriptomics, WT and APP23, 7 and 14 months (Gelber et al., 2026) | Derived circadian peak phase per gene and context (brain region, genotype, age), from DESeq2 harmonic regression | `phase_angle` (radians, ZT0 = 0), `amplitude_log2`, `padj`, `genotype`, `age_months`, `broad_region` |
| `xenium_cell_metadata_with_pathology_and_official_10x_annotation.csv.gz` | 10x Genomics Xenium Prime FFPE human ovarian cancer | Latent cell-cycle position per cell, inferred by `tricycle` | `tricyclePosition` (radians), `x_centroid`, `y_centroid` (µm), `pathology_region`, `Phase` (Seurat) |

Both files are processed tables. Raw expression matrices are not included.

## Quick start

```r
# R >= 4.1, run from the repository root
install.packages(c("circular", "nnet", "ranger", "rpart", "BAMBI", "loo"))
source("supplementary_core_analysis.R")
```

Results are written to `results_core/AD/` and `results_core/ovarian/` as CSV
tables and fitted models (`.rds`), together with `sessionInfo.txt`.
Switches at the top of the script:

| Option | Default | Effect |
|---|---|---|
| `RUN_MIXTURE` | `TRUE` | von Mises mixtures with `BAMBI` (K = 1 to 4, WAIC). The slowest step |
| `RUN_PN` | `FALSE` | Full-data projected-normal regression with `bpnreg` |
| `RUN_GP_SAVED` | `FALSE` | Scores the archived projected-GP predictions in `data/gp_output/` |

## Script blocks

| Block | Analysis | Manuscript |
|---|---|---|
| `inputs` | Read the two tables and apply analysis filters | §5 |
| `adsummary` | Mean direction, MRL, and Rayleigh test for AD peak phase | §2.2, §2.3, §5.1 |
| `predictor` | Harmonic regression and random forest with phase as predictor | §3.1, §3.3, §5.1 |
| `nn` | Neural network predicting (sin θ, cos θ) | §3.3, §5.1 |
| `ovsummary` | Summaries by pathology region and sample-size effect | §2.2, §2.3, §5.2 |
| `cart` | Classification tree on (cos θ, sin θ) | §3.3, §5.2 |
| `spatial` | Empirical circular semivariogram, pooled and tumor only | §4.1, §5.2 |
| `mixture` | von Mises mixture clustering | §2.4, §5.1, §5.2 |
| `pn` | Projected-normal regression | §3.2, §5.1 |
| `gp` | Evaluation of archived projected-GP predictions | §4.2, §5.2 |

Stand-alone examples for `lm.circular`, `brms`, and `CircSpaceTime` are in the
Supplementary Note. Fitting the projected GP requires `CircSpaceTime` 0.9.0 from
the [CRAN archive](https://cran.r-project.org/src/contrib/Archive/CircSpaceTime/).

## Notes

- Random steps use fixed seeds. MCMC-based results can differ slightly across platforms.
- `simulate_toy_data.R` writes synthetic files with the **same names** into `data/`
  and overwrites the real inputs. Run it only in a separate copy of the repository.

## Citation

Shin J, Yoo J, Cho Y, et al. Circular Data Analysis for Spatial Omics. Manuscript, 2026.

Questions and bug reports: [GitHub issues](https://github.com/c16267/circular_ST/issues).
