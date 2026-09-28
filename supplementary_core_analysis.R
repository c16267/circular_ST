# Circular Data Analysis for Spatial Omics: essential downstream analyses.
# The numbered blocks are identical to the R listings in the supplement.
# Run from the repository root with source("supplementary_core_analysis.R").
# Inputs are the two processed CSV tables, not raw expression matrices.

# --- BLOCK setup ---
DATA_DIR <- "data"       # change to your local data folder
OUT_DIR <- "results_core"
RUN_MIXTURE <- TRUE      # BAMBI: can take substantial time
RUN_PN <- FALSE          # optional full-data projected-normal regression
RUN_GP_SAVED <- FALSE    # optional scoring of archived GP predictions
SEED <- 16267L
N_CLUSTER <- 5000L
MIX_ITER <- 5000L
MIX_CHAINS <- 3L
N_SPATIAL <- 1500L
PN_ITS <- 3000L
PN_BURN <- 500L
PN_LAG <- 10L

# Install once, in an interactive R session, if needed:
# install.packages(c("circular", "nnet", "ranger", "rpart"))
# install.packages(c("BAMBI", "loo"))  # if RUN_MIXTURE
# install.packages("bpnreg")            # if RUN_PN
pkgs <- c("circular", "nnet", "ranger", "rpart",
          if (RUN_MIXTURE) c("BAMBI", "loo"),
          if (RUN_PN) "bpnreg")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1),
                            quietly = TRUE)]
if (length(missing_pkgs)) stop("Install packages: ",
                              paste(missing_pkgs, collapse = ", "))
stopifnot(dir.exists(DATA_DIR), MIX_ITER > 0L,
          PN_ITS > PN_BURN, PN_LAG >= 1L)
options(contrasts = c("contr.treatment", "contr.poly"))
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(SEED)
for (p in c(OUT_DIR, file.path(OUT_DIR, c("AD", "ovarian"))))
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
capture.output(sessionInfo(), file = file.path(OUT_DIR, "sessionInfo_start.txt"))
write_out <- function(x, name, dataset) {
  utils::write.csv(x, file.path(OUT_DIR, dataset, name),
                   row.names = FALSE, na = "NA")
  invisible(x)
}

read_input <- function(stem) {
  paths <- file.path(DATA_DIR, c(paste0(stem, ".gz"), stem))
  path <- paths[file.exists(paths)][1L]
  if (is.na(path)) stop("Missing input: ", file.path(DATA_DIR, stem))
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  utils::read.csv(con, stringsAsFactors = FALSE, check.names = FALSE)
}
check_cols <- function(d, cols) {
  absent <- setdiff(cols, names(d))
  if (length(absent)) stop("Missing columns: ", paste(absent, collapse = ", "))
}
as_number <- function(x, name) {
  y <- suppressWarnings(as.numeric(as.character(x)))
  if (any(!is.na(x) & is.na(y))) stop("Non-numeric values in ", name)
  y
}
# --- END setup ---

# --- BLOCK helpers ---
wrap_2pi <- function(x) x %% (2 * pi)
wrap_pi <- function(x) (x + pi) %% (2 * pi) - pi
rad_to_hour <- function(x) x * 24 / (2 * pi) # preserves signed shifts
circ_mean <- function(x) {
  z <- mean(exp(1i * x))
  if (Mod(z) < 1e-12) NA_real_ else wrap_2pi(Arg(z))
}
circ_summary <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) stop("No finite angles to summarize")
  R <- min(1, Mod(mean(exp(1i * x))))
  th <- circular::circular(x, units = "radians", modulo = "2pi",
                           zero = 0, rotation = "counter")
  p <- if (length(x) >= 3L) circular::rayleigh.test(th)$p.value else NA_real_
  data.frame(n = length(x), mean_rad = circ_mean(x),
             mean_over_pi = circ_mean(x) / pi,
             MRL = R, rayleigh_Z = length(x) * R^2,
             rayleigh_p_iid = p)
}
by_summary <- function(d, groups, angle) {
  good <- complete.cases(d[, groups, drop = FALSE])
  d <- d[good, , drop = FALSE]
  key <- do.call(interaction, c(d[groups], list(drop = TRUE)))
  ans <- lapply(split(seq_len(nrow(d)), key), function(i)
    cbind(d[i[1L], groups, drop = FALSE], circ_summary(d[[angle]][i])))
  out <- do.call(rbind, ans)
  rownames(out) <- NULL
  out
}
# Exact-size proportional sampling, without replacement.
# Largest remainders allocate the rows left after rounding down.
prop_index <- function(group, n) {
  stopifnot(!anyNA(group), n >= 1L)
  sets <- split(seq_along(group), as.character(group))
  n <- min(as.integer(n), length(group))
  quota <- n * lengths(sets) / length(group)
  size <- floor(quota)
  left <- n - sum(size)
  if (left > 0L) {
    j <- order(quota - size, decreasing = TRUE)[seq_len(left)]
    size[j] <- size[j] + 1L
  }
  unlist(Map(function(i, m) i[sample.int(length(i), m)],
             sets, size), use.names = FALSE)
}
circ_mae_h <- function(pred, obs)
  rad_to_hour(mean(abs(wrap_pi(pred - obs))))
# --- END helpers ---

# --- BLOCK inputs ---
ad <- read_input("DESeq2_region_age_gene_phase_all_genes.csv")
check_cols(ad, c("phase_angle", "amplitude_log2", "padj",
                 "genotype", "age_months", "broad_region"))
ad$row_id <- seq_len(nrow(ad)) # source-table row, before filtering
for (nm in c("phase_angle", "amplitude_log2", "padj", "age_months"))
  ad[[nm]] <- as_number(ad[[nm]], nm)
ad <- ad[is.finite(ad$phase_angle), , drop = FALSE]
ad$phase_angle <- wrap_2pi(ad$phase_angle)
if (any(!is.na(ad$genotype) & !ad$genotype %in% c("WT", "APP23")))
  stop("AD genotype must be WT or APP23")
if (any(!is.na(ad$age_months) & !ad$age_months %in% c(7, 14)))
  stop("AD age_months must be 7 or 14")
ad$genotype <- factor(ad$genotype, levels = c("WT", "APP23"))
ad$age_months <- factor(ad$age_months, levels = c(7, 14))
ad$broad_region <- factor(trimws(ad$broad_region))

# Same analysis-specific filters as the supplied GitHub script.
ad_desc <- subset(ad, is.finite(padj) & padj <= 0.10 &
                    is.finite(amplitude_log2) & amplitude_log2 >= 0.15)
amp_cut <- median(ad$amplitude_log2[is.finite(ad$amplitude_log2)])
ad_reg <- subset(ad, is.finite(amplitude_log2) &
                   amplitude_log2 > amp_cut & amplitude_log2 > 0)
ad_reg <- droplevels(ad_reg[complete.cases(ad_reg[, c(
  "genotype", "age_months", "broad_region")]), , drop = FALSE])
ad_reg$log_amp <- log2(ad_reg$amplitude_log2)
ad_reg$amp_response <- log2(1 + ad_reg$amplitude_log2)
ad_reg$sin_phase <- sin(ad_reg$phase_angle)
ad_reg$cos_phase <- cos(ad_reg$phase_angle)
stopifnot(nrow(ad_desc) > 0L, nrow(ad_reg) > 10L)

oc <- read_input(paste0("xenium_cell_metadata_with_pathology_",
                        "and_official_10x_annotation.csv"))
check_cols(oc, c("tricyclePosition", "x_centroid", "y_centroid",
                 "pathology_region", "Phase"))
oc$row_id <- seq_len(nrow(oc))
oc$theta <- as_number(oc$tricyclePosition, "tricyclePosition")
oc$x <- as_number(oc$x_centroid, "x_centroid")
oc$y <- as_number(oc$y_centroid, "y_centroid")
oc$pathology <- trimws(oc$pathology_region)
oc$seurat_phase <- trimws(oc$Phase)
oc <- subset(oc, is.finite(theta) & is.finite(x) & is.finite(y) &
               !is.na(pathology) &
               !tolower(pathology) %in% c("", "na", "unannotated"))
oc$theta <- wrap_2pi(oc$theta)
stopifnot(nrow(oc) > 10L)
write_out(data.frame(set = c("AD_all", "AD_descriptive", "AD_regression",
                            "ovarian_annotated"),
                     n = c(nrow(ad), nrow(ad_desc), nrow(ad_reg), nrow(oc))),
          "analysis_counts.csv", "")
write_out(data.frame(row_id = ad$row_id,
  descriptive = ad$row_id %in% ad_desc$row_id,
  regression = ad$row_id %in% ad_reg$row_id), "AD_analysis_rows.csv", "AD")
# --- END inputs ---

# --- BLOCK adsummary ---
ad_overall <- circ_summary(ad_desc$phase_angle)
ad_overall$mean_ZT_h <- rad_to_hour(ad_overall$mean_rad)
write_out(ad_overall, "AD_summary_overall.csv", "AD")
for (g in list("broad_region", c("broad_region", "genotype", "age_months"))) {
  s <- by_summary(ad_desc, g, "phase_angle")
  s$mean_ZT_h <- rad_to_hour(s$mean_rad)
  write_out(s, paste0("AD_summary_by_", paste(g, collapse = "_"), ".csv"), "AD")
}
print(ad_overall)
# --- END adsummary ---

# --- BLOCK predictor ---
# Harmonic association: the response transform matches block A10.
f0 <- stats::lm(amp_response ~ genotype + age_months + broad_region,
                data = ad_reg)
f1 <- update(f0, . ~ . + sin_phase + cos_phase)
a <- as.data.frame(stats::anova(f0, f1))
a$model <- c("context_only", "context_plus_phase")
write_out(a, "AD_harmonic_joint_test.csv", "AD")
write_out(data.frame(term = names(coef(f1)), estimate = unname(coef(f1))),
          "AD_harmonic_coefficients.csv", "AD")

set.seed(123)
i_rf <- sample.int(nrow(ad_reg), floor(0.8 * nrow(ad_reg)))
rf <- ranger::ranger(
  amp_response ~ sin_phase + cos_phase + genotype + age_months + broad_region,
  data = ad_reg[i_rf, ], num.trees = 500,
  importance = "permutation", seed = 123, num.threads = 1)
p_rf <- predict(rf, data = ad_reg[-i_rf, ])$predictions
y_rf <- ad_reg$amp_response[-i_rf]
write_out(data.frame(n_test = length(y_rf), cor = cor(p_rf, y_rf),
                     RMSE = sqrt(mean((p_rf - y_rf)^2))),
          "AD_rf_metrics.csv", "AD")
write_out(data.frame(variable = names(rf$variable.importance),
                     importance = unname(rf$variable.importance)),
          "AD_rf_importance.csv", "AD")
write_out(data.frame(row_id = ad_reg$row_id,
                     split = ifelse(seq_len(nrow(ad_reg)) %in% i_rf,
                                    "train", "test")), "AD_rf_split.csv", "AD")
write_out(data.frame(row_id = ad_reg$row_id[-i_rf], observed = y_rf,
                     predicted = p_rf), "AD_rf_predictions.csv", "AD")
saveRDS(rf, file.path(OUT_DIR, "AD", "AD_rf.rds"))
# --- END predictor ---

# --- BLOCK nn ---
nn_dat <- subset(ad, is.finite(padj) & padj < 0.05 &
                    is.finite(amplitude_log2))
nn_dat <- droplevels(nn_dat[complete.cases(nn_dat[, c(
  "genotype", "age_months", "broad_region")]), , drop = FALSE])
stopifnot(nrow(nn_dat) >= 20L)
X <- model.matrix(~ amplitude_log2 + age_months + broad_region + genotype,
                   data = nn_dat)
Y <- cbind(sin(nn_dat$phase_angle), cos(nn_dat$phase_angle))
set.seed(42)
n <- nrow(nn_dat)
idx <- sample.int(n)
tr <- idx[seq_len(floor(0.6 * n))]
va <- idx[seq.int(floor(0.6 * n) + 1L, floor(0.8 * n))]
te <- idx[seq.int(floor(0.8 * n) + 1L, n)]
# Learn scaling from training rows only, then freeze it.
rng <- range(X[tr, "amplitude_log2"])
span <- diff(rng)
if (span == 0) span <- 1
X[, "amplitude_log2"] <- (X[, "amplitude_log2"] - rng[1]) / span
fit_nn <- function(i, decay, maxit) nnet::nnet(
  x = X[i, , drop = FALSE], y = Y[i, , drop = FALSE],
  size = 32, linout = TRUE, decay = decay, maxit = maxit,
  MaxNWts = max(5000, (ncol(X) + 1) * 32 + 66), trace = FALSE)
angle_pred <- function(f, i) {
  p <- predict(f, X[i, , drop = FALSE])
  if (any(rowSums(p^2) < 1e-20))
    stop("Predicted vector is zero; its angle is undefined")
  wrap_2pi(atan2(p[, 1], p[, 2]))
}
decay_grid <- c(0.001, 0.01, 0.05)
val_loss <- vapply(decay_grid, function(d) {
  f <- fit_nn(tr, d, 400)
  circ_mae_h(angle_pred(f, va), nn_dat$phase_angle[va])
}, numeric(1))
best_decay <- decay_grid[which.min(val_loss)]
nn <- fit_nn(c(tr, va), best_decay, 600)
p_nn <- angle_pred(nn, te)
obs <- nn_dat$phase_angle[te]
# Original comparison baseline: means estimated from training only.
group <- interaction(nn_dat$broad_region, nn_dat$genotype, drop = TRUE)
base_means <- tapply(nn_dat$phase_angle[tr], group[tr], circ_mean)
p_base <- as.numeric(base_means[as.character(group[te])])
p_base[!is.finite(p_base)] <- circ_mean(nn_dat$phase_angle[tr])
if (any(!is.finite(p_base))) stop("Undefined baseline circular mean")
loss <- mean(1 - cos(p_nn - obs))
mu_test <- circ_mean(obs)
# If the test resultant is zero, every constant direction has loss one.
denom <- if (is.na(mu_test)) 1 else mean(1 - cos(obs - mu_test))
metrics <- data.frame(n_train = length(tr), n_validation = length(va),
  n_test = length(te), inputs = ncol(X), best_decay = best_decay,
  MAE_nn_h = circ_mae_h(p_nn, obs),
  MAE_baseline_h = circ_mae_h(p_base, obs), chord_loss = loss,
  circular_R2 = if (denom > 1e-12) 1 - loss / denom else NA_real_,
  optimizer_convergence = nn$convergence)
print(metrics)
write_out(metrics, "AD_nn_metrics.csv", "AD")
write_out(data.frame(decay = decay_grid, validation_MAE_h = val_loss),
          "AD_nn_tuning.csv", "AD")
write_out(data.frame(row_id = nn_dat$row_id[te], observed = obs,
  predicted = p_nn, baseline = p_base), "AD_nn_predictions.csv", "AD")
split <- rep("test", n); split[tr] <- "train"; split[va] <- "validation"
write_out(data.frame(row_id = nn_dat$row_id, split = split),
          "AD_nn_split.csv", "AD")
saveRDS(list(model = nn, amplitude_range = rng, columns = colnames(X)),
        file.path(OUT_DIR, "AD", "AD_nn.rds"))
# --- END nn ---

# --- BLOCK ovsummary ---
write_out(circ_summary(oc$theta), "ovarian_summary_overall.csv", "ovarian")
write_out(by_summary(oc, "pathology", "theta"),
          "ovarian_summary_by_pathology.csv", "ovarian")
# Sample-size illustration: use the spatial sample and a 5,000-cell sample.
set.seed(1)
i_sp <- prop_index(oc$pathology, N_SPATIAL)
set.seed(100)
i_5k <- sample.int(nrow(oc), min(5000L, nrow(oc)))
subsets <- list(proportional_1500 = i_sp, random_5000 = i_5k,
                all = seq_len(nrow(oc)))
size_effect <- do.call(rbind, lapply(names(subsets), function(nm) {
  s <- by_summary(oc[subsets[[nm]], ], "pathology", "theta")
  cbind(sample = nm, sample_n = length(subsets[[nm]]), s)
}))
write_out(size_effect, "ovarian_sample_size_effect.csv", "ovarian")
# 24-bin proportions underlying the pathology-phase heatmap.
bins <- cut(oc$theta, seq(0, 2 * pi, length.out = 25),
             right = FALSE, include.lowest = TRUE, labels = FALSE)
h <- as.data.frame(table(pathology = oc$pathology,
                         phase_bin = factor(bins, levels = 1:24)))
h$proportion <- h$Freq / ave(h$Freq, h$pathology, FUN = sum)
write_out(h, "ovarian_pathology_phase_bins.csv", "ovarian")
# --- END ovsummary ---

# --- BLOCK cart ---
d <- subset(oc, seurat_phase %in% c("G1", "S", "G2M"))
d$phase <- factor(d$seurat_phase, levels = c("G1", "S", "G2M"))
d$cos_theta <- cos(d$theta); d$sin_theta <- sin(d$theta)
sets <- split(seq_len(nrow(d)), d$phase)
n_class <- min(12000L, min(lengths(sets)))
if (n_class < 4L) stop("Too few cells in at least one Seurat phase")
set.seed(20260820)
i_bal <- unlist(lapply(sets, function(i)
  i[sample.int(length(i), n_class)]), use.names = FALSE)
bal <- d[i_bal, ]
tr_tree <- unlist(lapply(split(seq_len(nrow(bal)), bal$phase), function(i)
  i[sample.int(length(i), floor(0.7 * length(i)))]), use.names = FALSE)
train <- bal[tr_tree, ]; test <- bal[-tr_tree, ]
tree <- rpart::rpart(phase ~ sin_theta + cos_theta, data = train,
  method = "class", parms = list(split = "information"),
  control = rpart::rpart.control(maxdepth = 3, minbucket = 750,
                                 cp = 0.005, xval = 10))
pred <- factor(predict(tree, test, type = "class"), levels = levels(d$phase))
cm <- table(observed = test$phase, predicted = pred)
recall <- diag(cm) / rowSums(cm)
write_out(data.frame(n_train = nrow(train), n_test = nrow(test),
  accuracy = sum(diag(cm)) / sum(cm), balanced_accuracy = mean(recall)),
  "ovarian_CART_metrics.csv", "ovarian")
write_out(data.frame(phase = names(recall), recall = unname(recall)),
          "ovarian_CART_recall.csv", "ovarian")
write_out(as.data.frame(cm), "ovarian_CART_confusion.csv", "ovarian")
write_out(data.frame(row_id = bal$row_id, phase = bal$phase,
  split = ifelse(seq_len(nrow(bal)) %in% tr_tree, "train", "test")),
  "ovarian_CART_split.csv", "ovarian")
grid <- data.frame(theta = seq(0, 2 * pi, length.out = 721))
grid$cos_theta <- cos(grid$theta); grid$sin_theta <- sin(grid$theta)
grid$predicted_phase <- predict(tree, grid, type = "class")
write_out(grid, "ovarian_CART_circle_predictions.csv", "ovarian")
saveRDS(tree, file.path(OUT_DIR, "ovarian", "ovarian_CART.rds"))
# --- END cart ---

# --- BLOCK spatial ---
circ_vario <- function(d, n_bins, label) {
  if (nrow(d) < 3L) stop("At least three cells required for ", label)
  D <- as.matrix(dist(d[, c("x", "y")]))
  A <- 1 - cos(outer(d$theta, d$theta, "-"))
  ut <- upper.tri(D)
  distance <- D[ut]; dissimilarity <- A[ut]
  br <- unique(as.numeric(quantile(distance,
                    seq(0, 1, length.out = n_bins + 1L))))
  if (length(br) < 2L) stop("No variation in pair distances")
  b <- cut(distance, br, include.lowest = TRUE)
  R <- min(1, Mod(mean(exp(1i * d$theta))))
  indep <- 1 - R^2
  ans <- do.call(rbind, lapply(split(seq_along(distance), b,
                                   drop = TRUE), function(i)
    data.frame(mean_distance_um = mean(distance[i]),
      mean_dissimilarity = mean(dissimilarity[i]), n_pairs = length(i))))
  ans$group <- label; ans$n_cells <- nrow(d); ans$independence <- indep
  ans$relative_dissimilarity <- if (indep > 1e-12)
    ans$mean_dissimilarity / indep else NA_real_
  rownames(ans) <- NULL
  ans
}
spatial_sample <- oc[i_sp, ]
tumor <- subset(spatial_sample, pathology == "Tumor")
v <- rbind(circ_vario(spatial_sample, 15L, "pooled"),
           circ_vario(tumor, 12L, "tumor"))
write_out(v, "ovarian_semivariogram.csv", "ovarian")
write_out(spatial_sample[, c("row_id", "theta", "x", "y", "pathology")],
          "ovarian_spatial_sample.csv", "ovarian")
# --- END spatial ---

# --- BLOCK mixture ---
if (RUN_MIXTURE) {
  run_mixture <- function(d, angle, groups, dataset, seed) {
    set.seed(seed)
    i <- sample.int(nrow(d), min(N_CLUSTER, nrow(d)))
    s <- d[i, , drop = FALSE]
    fit_paths <- file.path(OUT_DIR, dataset,
                           paste0("mixture_K", 1:4, ".rds"))
    w <- numeric(4L)
    message("Mixture fits: ", dataset, " (n = ", nrow(s), ")")
    for (k in 1:4) {
      message("  Fitting K = ", k)
      set.seed(seed + k)
      f <- BAMBI::fit_angmix(model = "vm", data = s[[angle]], ncomp = k,
        n.iter = MIX_ITER, n.chains = MIX_CHAINS,
        chains_parallel = FALSE, return_llik_contri = TRUE,
        show.progress = FALSE)
      w[k] <- loo::waic(f)$estimates["waic", "Estimate"]
      saveRDS(f, fit_paths[k])
      rm(f); invisible(gc())
    }
    if (any(!is.finite(w))) stop("Non-finite mixture WAIC for ", dataset)
    best_k <- which.min(w)
    best <- readRDS(fit_paths[best_k])
    if (best_k > 1L) best <- BAMBI::fix_label(best)
    allocation <- as.integer(BAMBI::latent_allocation(best))
    stopifnot(length(allocation) == nrow(s))
    s$cluster <- factor(allocation, levels = seq_len(best_k))
    write_out(data.frame(K = 1:4, WAIC = w, selected = 1:4 == best_k),
              "mixture_WAIC.csv", dataset)
    # Model-based component estimates and empirical allocated-group summaries
    # are separate: a hard-assignment mean is not a fitted component mean.
    capture.output(BAMBI::pointest(best), file = file.path(
      OUT_DIR, dataset, "mixture_parameter_estimates.txt"))
    write_out(by_summary(s, "cluster", angle),
              "mixture_allocated_summaries.csv", dataset)
    write_out(s[, c("row_id", angle, groups, "cluster")],
              "mixture_sample_assignments.csv", dataset)
    for (g in groups) {
      tab <- as.data.frame(table(group = s[[g]], cluster = s$cluster))
      tab$proportion <- tab$Freq / ave(tab$Freq, tab$group, FUN = sum)
      write_out(tab, paste0("mixture_by_", g, ".csv"), dataset)
    }
    saveRDS(best, file.path(OUT_DIR, dataset, "mixture_best.rds"))
    # Candidate fits have already been saved one at a time.
    print(data.frame(dataset = dataset, selected_K = best_k))
  }
  run_mixture(ad_desc, "phase_angle",
    c("broad_region", "genotype", "age_months"), "AD", SEED)
  run_mixture(oc, "theta", c("pathology", "seurat_phase"), "ovarian", 100L)
}
# --- END mixture ---

# --- BLOCK pn ---
if (RUN_PN) {
  # Full regression set; no implicit subsampling or standardization.
  pn_data <- ad_reg[, c("phase_angle", "log_amp", "genotype",
                        "age_months", "broad_region")]
  pn0 <- bpnreg::bpnr(phase_angle ~ 1, data = pn_data,
    its = PN_ITS, burn = PN_BURN, n.lag = PN_LAG, seed = 1)
  pn1 <- bpnreg::bpnr(
    phase_angle ~ log_amp + genotype + age_months + broad_region,
    data = pn_data, its = PN_ITS, burn = PN_BURN, n.lag = PN_LAG, seed = 1)
  saveRDS(list(null = pn0, main = pn1),
          file.path(OUT_DIR, "AD", "AD_projected_normal.rds"))
  fit0 <- bpnreg::fit(pn0); fit1 <- bpnreg::fit(pn1)
  write_out(data.frame(criterion = rownames(fit1),
    null = fit0$Statistic, main = fit1$Statistic,
    difference = fit1$Statistic - fit0$Statistic),
    "AD_PN_model_comparison.csv", "AD")
  capture.output(bpnreg::coef_circ(pn1, type = "continuous"),
    file = file.path(OUT_DIR, "AD", "AD_PN_continuous_effects.txt"))
  capture.output(bpnreg::coef_circ(pn1, type = "categorical"),
    file = file.path(OUT_DIR, "AD", "AD_PN_categorical_effects.txt"))
  # Only categorical angular differences are converted to hour shifts.
  sh <- as.data.frame(pn1$circ.coef.cat)
  check_cols(sh, c("mean", "LB", "UB"))
  sh$term <- rownames(sh)
  # bpnreg 2.0.3 defines these contrasts as reference minus level.
  # Keep that orientation; do not silently reverse a reported shift.
  sh$mean_h <- rad_to_hour(wrap_pi(sh$mean))
  sh$lower_h <- sh$mean_h - rad_to_hour(wrap_2pi(sh$mean - sh$LB))
  sh$upper_h <- sh$mean_h + rad_to_hour(wrap_2pi(sh$UB - sh$mean))
  write_out(sh, "AD_PN_categorical_shifts.csv", "AD")
  # Plain trace plots for diagnostic use, without publication styling.
  grDevices::pdf(file.path(OUT_DIR, "AD", "AD_PN_traces.pdf"))
  for (nm in c("beta1", "beta2"))
    matplot(pn1[[nm]], type = "l", lty = 1,
            xlab = "Saved draw", ylab = nm)
  grDevices::dev.off()
}
# --- END pn ---

# --- BLOCK gp ---
if (RUN_GP_SAVED) {
  # This evaluates archived predictions; it does not fit a GP.
  gp_path <- file.path(DATA_DIR, "gp_output", "test_with_pred.csv.gz")
  if (!file.exists(gp_path)) gp_path <- sub("\\.gz$", "", gp_path)
  if (!file.exists(gp_path)) stop("Missing archived GP predictions: ", gp_path)
  con <- if (grepl("\\.gz$", gp_path)) gzfile(gp_path, "rt") else file(gp_path, "rt")
  gp <- tryCatch(utils::read.csv(con), finally = close(con))
  check_cols(gp, c("tricyclePosition", "pred_mean"))
  if (!all(is.finite(gp$tricyclePosition) & is.finite(gp$pred_mean)))
    stop("Archived GP angles must be finite numeric values")
  e <- wrap_pi(gp$pred_mean - gp$tricyclePosition)
  write_out(data.frame(n = nrow(gp), MAE_rad = mean(abs(e)),
                       chord_loss = mean(1 - cos(e))),
            "ovarian_GP_saved_metrics.csv", "ovarian")
}
capture.output(sessionInfo(), file = file.path(OUT_DIR, "sessionInfo.txt"))
saveRDS(list(DATA_DIR = normalizePath(DATA_DIR), SEED = SEED,
  RUN_MIXTURE = RUN_MIXTURE, RUN_PN = RUN_PN, RUN_GP_SAVED = RUN_GP_SAVED,
  N_CLUSTER = N_CLUSTER, MIX_ITER = MIX_ITER, MIX_CHAINS = MIX_CHAINS,
  N_SPATIAL = N_SPATIAL, PN_ITS = PN_ITS, PN_BURN = PN_BURN, PN_LAG = PN_LAG,
  amplitude_median_cutoff = amp_cut), file.path(OUT_DIR, "run_settings.rds"))
message("Finished: ", normalizePath(OUT_DIR))
# --- END gp ---
