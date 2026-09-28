# ==============================================================================
# simulate_toy_data.R
#
# Generates small synthetic versions of the two input files so that
# circular_spatial_omics_analysis.R can be test-run end to end without access
# to the real data. Column names and value ranges match the real files.
# NOT used for any manuscript result.
#
# Usage: Rscript simulate_toy_data.R      (writes into ./data/)
# ==============================================================================

set.seed(1)
dir.create("data", showWarnings = FALSE)

rvm <- function(n, mu, kappa) {           # simple von Mises sampler (wrapped)
  suppressMessages(requireNamespace("circular"))
  as.numeric(circular::rvonmises(n, circular::circular(mu), kappa))
}

## ---- AD: DESeq2_region_age_gene_phase_all_genes.csv -------------------------

n_ad <- 30000
regions <- c("cortex", "basal_ganglia", "hippocampus", "piriform_amygdala",
             "fiber_tracts", "basal_forebrain_basal_ganglia", "amygdala",
             "meninges", "thalamus", "retrohippocampus", "ventricle_choroid")
region_p  <- c(0.35, 0.16, 0.12, 0.08, 0.06, 0.06, 0.05, 0.04, 0.04, 0.02, 0.02)
region_mu <- setNames(runif(length(regions), 0, 2 * pi), regions)

ad <- data.frame(
  gene         = paste0("gene", sample(1:4000, n_ad, replace = TRUE)),
  broad_region = sample(regions, n_ad, replace = TRUE, prob = region_p),
  genotype     = sample(c("WT", "APP23"), n_ad, replace = TRUE, prob = c(0.4, 0.6)),
  age_months   = sample(c(7, 14), n_ad, replace = TRUE)
)
shift <- ifelse(ad$genotype == "APP23", -0.4, 0) + ifelse(ad$age_months == 14, 0.3, 0)
ad$phase_angle    <- (rvm(n_ad, 0, 1.2) + region_mu[ad$broad_region] + shift) %% (2 * pi)
ad$phase_hours    <- ad$phase_angle * 24 / (2 * pi)
ad$amplitude_log2 <- rlnorm(n_ad, log(0.25), 0.6)
ad$padj           <- rbeta(n_ad, 0.4, 2)
write.csv(ad, "data/DESeq2_region_age_gene_phase_all_genes.csv", row.names = FALSE)

## ---- Ovarian: xenium_cell_metadata_..._annotation.csv -----------------------

n_ov <- 40000
paths <- c("Tumor", "Blood vessels", "Walthard rests (benign)",
           "Fallopian tube", "Immune cells", "Necrosis")
path_p  <- c(0.62, 0.14, 0.09, 0.08, 0.04, 0.03)
path_mu <- setNames(c(1.6, 0.0, 1.9, 0.3, 0.4, 2.9), paths)
path_k  <- setNames(c(0.4, 2.5, 0.4, 1.0, 1.5, 1.5), paths)

# Blocky tissue layout so pathology regions are spatially organized
ov <- data.frame(
  x_centroid = runif(n_ov, 0, 6000),
  y_centroid = runif(n_ov, 0, 4000)
)
blk <- (ov$x_centroid %/% 1500 + 2 * (ov$y_centroid %/% 1300)) %% length(paths) + 1
ov$pathology_region <- sample(paths, n_ov, replace = TRUE, prob = path_p)
mixidx <- runif(n_ov) < 0.6
ov$pathology_region[mixidx] <- paths[blk[mixidx]]

ov$tricyclePosition <- mapply(
  function(m, k) rvm(1, m, k),
  path_mu[ov$pathology_region], path_k[ov$pathology_region]) %% (2 * pi)
ov$tricyclePosition_0_1 <- ov$tricyclePosition / (2 * pi)

# Seurat phase labels correlated with the angle
cut_phase <- function(th)
  ifelse(th < 1.0 | th > 5.5, "G1", ifelse(th < 2.6, "S", "G2M"))
ov$Phase <- cut_phase((ov$tricyclePosition + rnorm(n_ov, 0, 0.5)) %% (2 * pi))
ov$Phase[sample(n_ov, 200)] <- "Undecided"

ov$official_10x_cell_type <- sample(
  c("Tumor Cells", "Proliferative Tumor Cells", "Tumor Associated Fibroblasts",
    "Macrophages", "T and NK Cells", "Pericytes", "Ciliated Epithelial Cells"),
  n_ov, replace = TRUE)

write.csv(ov,
  "data/xenium_cell_metadata_with_pathology_and_official_10x_annotation.csv",
  row.names = FALSE)

message("Toy data written to ./data/")
