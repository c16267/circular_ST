# circular_ST

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23022286.svg)](https://doi.org/10.5281/zenodo.23022286)

Code for **Circular Data Analysis for Spatial Omics** (Shin, Yoo, Cho, et al.).
The repository reproduces the two case studies in the manuscript and accompanies
the R examples in the Supplementary Note.

The R code and processed data are archived together at Zenodo.

- Version v1, used in the manuscript: [10.5281/zenodo.23022287](https://doi.org/10.5281/zenodo.23022287)
- All versions, resolving to the latest: [10.5281/zenodo.23022286](https://doi.org/10.5281/zenodo.23022286)

## Contents

```
circular_ST/
├── supplementary_core_analysis.R   # main analysis script
├── simulate_toy_data.R             # optional synthetic inputs with the same columns
├── CITATION.cff
├── LICENSE
└── data/                           # download from Zenodo (not tracked by git)
    ├── DESeq2_region_age_gene_phase_all_genes.csv.gz
    ├── anatomical_pseudobulk_metadata.csv
    ├── xenium_cell_metadata_with_pathology_and_official_10x_annotation.csv.gz
    └── gp_output/test_with_pred.csv.gz    # archived projected-GP predictions
```

## Data

| File | Content | Main columns |
|---|---|---|
| `DESeq2_region_age_gene_phase_all_genes.csv[.gz]` | AD mouse brain, WT and APP23, 7 and 14 months (Gelber et al., 2026). Derived circadian peak phase per gene and context from DESeq2 harmonic regression | `phase_angle` (radians, ZT0 = 0), `amplitude_log2`, `padj`, `genotype`, `age_months`, `broad_region` |
| `anatomical_pseudobulk_metadata.csv` | AD pseudobulk profiles linked to samples, GEO sample accessions, and brain regions | `pseudobulk_id`, `sample_id`, `gsm_accession`, `anatomical_label`, `broad_region` |
| `xenium_cell_metadata_with_pathology_and_official_10x_annotation.csv.gz` | 10x Genomics Xenium Prime FFPE human ovarian cancer. Latent cell-cycle position per cell from `tricycle` | `tricyclePosition` (radians), `x_centroid`, `y_centroid` (µm), `pathology_region`, `Phase` (Seurat) |

The metadata file documents provenance and is not read by the main script.
Download the files from the [Zenodo record](https://zenodo.org/records/23022287)
into `data/`. The script reads either `.csv` or `.csv.gz`. From R:

```r
base  <- "https://zenodo.org/records/23022287/files/"
files <- c("DESeq2_region_age_gene_phase_all_genes.csv.gz",
           "anatomical_pseudobulk_metadata.csv",
           "xenium_cell_metadata_with_pathology_and_official_10x_annotation.csv.gz")
dir.create("data/gp_output", recursive = TRUE, showWarnings = FALSE)
options(timeout = 3600)   # about 0.9 GB in total
for (f in files)
  download.file(paste0(base, f, "?download=1"), file.path("data", f), mode = "wb")
download.file(paste0(base, "test_with_pred.csv.gz?download=1"),
              "data/gp_output/test_with_pred.csv.gz", mode = "wb")   # optional
```

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

If you use this code or data, please cite the manuscript and the Zenodo record.

- Shin J, Yoo J, Cho Y, et al. Circular Data Analysis for Spatial Omics. Manuscript, 2026.
- Shin J, Yoo J, Cho Y, et al. R code and processed data for Circular Data Analysis for
  Spatial Omics (v1). Zenodo, 2026. https://doi.org/10.5281/zenodo.23022287

GitHub also provides a formatted citation under **Cite this repository** (from `CITATION.cff`).

Questions and bug reports: [GitHub issues](https://github.com/c16267/circular_ST/issues).

## License

Code is released under the MIT License (see `LICENSE`). The processed data follow
the terms stated on the Zenodo record.
