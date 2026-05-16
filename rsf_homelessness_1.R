# Description
# -----------
# Fits a Random Survival Forest to predict first entry into homelessness.
# The person-month panel is collapsed to one row per person (last observed
# covariate snapshot before event or censoring), then fit with
# Surv(survival_time, event_observed) using the ranger package.
#
# ranger is used instead of randomForestSRC for performance: it is 10-20x
# faster on large datasets while producing identical results.
#
# Out-of-sample prediction uses a landmark approach: at the end of the
# training period, we snapshot each at-risk person's most recent covariates
# and predict their homeless-entry risk over the next HORIZON_MONTHS.
# This mirrors Figure 7 in the paper.
#
# Outputs (written to OUTPUT_DIR)
# --------------------------------
#   rsf_model.rds                -- saved model object
#   rsf_variable_importance.csv  -- VIMP scores (if computed)
#   rsf_variable_importance.png  -- bar chart of VIMP (if computed)
#   rsf_survival_curves.png      -- predicted survival curves by risk tertile
#   rsf_oos_predictions.csv      -- person-level OOS predicted risk scores
#   rsf_calibration_plot.png     -- predicted vs. realised risk (Figure 7 style)
#   rsf_diagnostics.txt          -- key diagnostics saved to file
#
# To adapt to a different dataset, edit only the CONFIG section below.
# =============================================================================


# =============================================================================
# 0.  INSTALL PACKAGES (run once; comment out afterwards)
# =============================================================================
# install.packages(c("ranger", "readstata13", "dplyr", "ggplot2",
#                    "survival", "tibble", "tidyr", "scales"))


# =============================================================================
# 1.  LOAD LIBRARIES
# =============================================================================
suppressPackageStartupMessages({
  library(ranger)        # fast RSF engine (replaces randomForestSRC)
  library(readstata13)   # read .dta files
  library(dplyr)
  library(ggplot2)
  library(survival)      # concordance(), Surv()
  library(tibble)
  library(tidyr)
  library(scales)
})

set.seed(42)


# =============================================================================
# 2.  CONFIG  <-- edit this section to adapt to a different dataset
# =============================================================================

# -- File paths ----------------------------------------------------------------
TRAIN_FILE <- "/Users/pierreloup/Library/CloudStorage/Dropbox/research/Homeless/lasso_simulation/sim_homelessness_pm.dta"
OOS_FILE   <- NULL     # path to external OOS file, or NULL to skip
OUTPUT_DIR <- "/Users/pierreloup/Library/CloudStorage/Dropbox/research/Homeless/lasso_simulation/RSF_outputs"

# -- Panel identifiers ---------------------------------------------------------
ID_VAR    <- "id"
TIME_VAR  <- "month"
EVENT_VAR <- "event"   # 1 = first homelessness entry, 0 = at risk
YEAR_VAR  <- "year"

# -- Train / test split --------------------------------------------------------
CUTOFF_YEAR    <- 2023   # years >= this are held out as the test period
HORIZON_MONTHS <- 12     # prediction horizon for OOS exercise (months)

# -- Predictor variables -------------------------------------------------------
PREDICTORS <- c(
  "age", "female", "race_black", "race_hispanic", "race_other",
  "hh_size", "has_children", "married",
  "disability", "log_income", "employed", "rent_burden",
  "mental_health", "substance_use",
  "eviction_3mo", "job_loss_3mo", "hosp_3mo"
)

# -- RSF hyperparameters -------------------------------------------------------
# For a quick test run : N_TREES = 50,  NODE_SIZE = 50
# For final results    : N_TREES = 500, NODE_SIZE = 15
N_TREES      <- 500     # number of trees
NODE_SIZE    <- 15     # min observations per terminal node (larger = faster)
MTRY         <- NULL   # features tried per split; NULL = sqrt(p) (recommended)
N_THREADS    <- parallel::detectCores() - 1   # cores for parallelisation

# -- VIMP switch ---------------------------------------------------------------
# Set to TRUE only when you want variable importance scores.
# Computing VIMP roughly doubles runtime so keep FALSE during development.
COMPUTE_VIMP <- TRUE


# =============================================================================
# 3.  HELPER FUNCTIONS
# =============================================================================

section <- function(title) {
  cat("\n", strrep("=", 70), "\n  ", title, "\n", strrep("=", 70), "\n",
      sep = "")
}

#' Load a person-month panel from .dta, .csv, or .rds
load_panel <- function(filepath) {
  ext <- tolower(tools::file_ext(filepath))
  df  <- switch(ext,
    dta = readstata13::read.dta13(filepath, convert.factors = FALSE),
    csv = read.csv(filepath),
    rds = readRDS(filepath),
    stop("Unsupported file format: ", ext)
  )
  df <- as.data.frame(df)
  cat(sprintf("  Loaded %s person-months, %s individuals.\n",
              format(nrow(df), big.mark = ","),
              format(length(unique(df[[ID_VAR]])), big.mark = ",")))
  df
}

#' Keep only rows where the person is still at risk (up to and including
#' the event row).
trim_to_risk_set <- function(df) {
  df %>%
    arrange(.data[[ID_VAR]], .data[[TIME_VAR]]) %>%
    group_by(.data[[ID_VAR]]) %>%
    filter(cumsum(.data[[EVENT_VAR]]) <= 1) %>%
    ungroup()
}

#' Collapse a person-month panel to one row per person.
#'
#' survival_time  = month of first event, OR last observed month if censored.
#' event_observed = 1 if the person ever became homeless, 0 otherwise.
#' Covariates     = values from the last observed month before event/censoring.
collapse_to_person <- function(df, predictors) {
  df %>%
    trim_to_risk_set() %>%
    arrange(.data[[ID_VAR]], .data[[TIME_VAR]]) %>%
    group_by(.data[[ID_VAR]]) %>%
    slice_tail(n = 1) %>%        # last at-risk row = covariate snapshot
    ungroup() %>%
    transmute(
      id             = .data[[ID_VAR]],
      survival_time  = as.numeric(.data[[TIME_VAR]]),
      event_observed = as.integer(.data[[EVENT_VAR]]),
      across(all_of(predictors), as.numeric)
    ) %>%
    as.data.frame()
}

#' For each individual still at risk at landmark_month, return their
#' most recent covariate snapshot (one row per person).
landmark_snapshot <- function(df, landmark_month, predictors) {
  df %>%
    filter(.data[[TIME_VAR]] <= landmark_month,
           .data[[EVENT_VAR]] == 0) %>%
    arrange(.data[[ID_VAR]], .data[[TIME_VAR]]) %>%
    group_by(.data[[ID_VAR]]) %>%
    slice_tail(n = 1) %>%
    ungroup() %>%
    transmute(
      id             = .data[[ID_VAR]],
      survival_time  = as.numeric(landmark_month + HORIZON_MONTHS),
      event_observed = 0L,
      across(all_of(predictors), as.numeric)
    ) %>%
    as.data.frame()
}

#' Convert a ranger survival prediction to a scalar risk score per person.
#' Risk score = mean cumulative hazard = mean(-log(S(t))) across time points.
#' Higher score = higher predicted risk. The small constant avoids log(0).
survival_to_risk <- function(surv_matrix) {
  rowMeans(-log(surv_matrix + 1e-10))
}


# =============================================================================
# 4.  LOAD AND PREPARE DATA
# =============================================================================
section("Loading and preparing data")

df_panel <- load_panel(TRAIN_FILE)

df_train_panel <- df_panel %>% filter(.data[[YEAR_VAR]] <  CUTOFF_YEAR)
df_test_panel  <- df_panel %>% filter(.data[[YEAR_VAR]] >= CUTOFF_YEAR)

cat(sprintf("\n  Training period : years < %d  (%s person-months)\n",
            CUTOFF_YEAR, format(nrow(df_train_panel), big.mark = ",")))
cat(sprintf("  Test period     : years >= %d  (%s person-months)\n",
            CUTOFF_YEAR, format(nrow(df_test_panel),  big.mark = ",")))

# Collapse to one row per person: last covariate snapshot before event/censoring
cat("\n  Collapsing panel to one row per person...\n")
df_train <- collapse_to_person(df_train_panel, PREDICTORS)

n_train  <- nrow(df_train)
n_events <- sum(df_train$event_observed)

cat(sprintf("  Individuals in training : %s\n", format(n_train,  big.mark = ",")))
cat(sprintf("  Events (homelessness)   : %s (%.1f%%)\n",
            format(n_events, big.mark = ","), 100 * n_events / n_train))
cat(sprintf("  Predictors              : %d\n", length(PREDICTORS)))
cat("  Survival time range     :",
    min(df_train$survival_time), "to", max(df_train$survival_time), "months\n")

# Sanity checks
stopifnot(
  "survival_time must be positive"        = all(df_train$survival_time > 0),
  "event_observed must be 0 or 1"         = all(df_train$event_observed %in% 0:1),
  "no missing values allowed in training" = !any(is.na(df_train))
)
cat("  Sanity checks passed.\n")


# =============================================================================
# 5.  FIT THE RANDOM SURVIVAL FOREST
# =============================================================================
section("Fitting Random Survival Forest (ranger)")

rsf_formula <- as.formula(
  paste("Surv(survival_time, event_observed) ~",
        paste(PREDICTORS, collapse = " + "))
)

mtry_val <- if (is.null(MTRY)) floor(sqrt(length(PREDICTORS))) else MTRY

cat("\n  Formula:", deparse(rsf_formula), "\n")
cat(sprintf("\n  Parameters: num.trees=%d  min.node.size=%d  mtry=%d  threads=%d\n",
            N_TREES, NODE_SIZE, mtry_val, N_THREADS))
cat(sprintf("  VIMP: %s\n",
            ifelse(COMPUTE_VIMP, "permutation (slower)", "off (fast mode)")))
cat("\n  Training...\n")

t0 <- proc.time()

rsf_model <- ranger(
  formula       = rsf_formula,
  data          = df_train,
  num.trees     = N_TREES,
  min.node.size = NODE_SIZE,
  mtry          = mtry_val,
  importance    = ifelse(COMPUTE_VIMP, "permutation", "none"),
  num.threads   = N_THREADS,
  seed          = 42
)

t_elapsed <- (proc.time() - t0)["elapsed"]
cat(sprintf("  Done in %.1f seconds.\n", t_elapsed))


# =============================================================================
# 6.  IN-SAMPLE DIAGNOSTICS
# =============================================================================
section("In-sample diagnostics")

# ranger stores OOB prediction error = 1 - Harrell's C for survival forests
oob_error <- rsf_model$prediction.error
oob_cstat <- 1 - oob_error

cat(sprintf("\n  OOB prediction error : %.4f\n", oob_error))
cat(sprintf("  OOB C-statistic      : %.4f\n",  oob_cstat))
cat("    (0.5 = random, 1.0 = perfect discrimination)\n")
cat("\n  Model summary:\n")
print(rsf_model)


# =============================================================================
# 7.  VARIABLE IMPORTANCE (VIMP)
# =============================================================================
section("Variable importance")

if (!COMPUTE_VIMP) {

  cat("\n  VIMP is off (COMPUTE_VIMP = FALSE in CONFIG).\n")
  cat("  Set COMPUTE_VIMP <- TRUE and re-run to get importance scores.\n")
  vimp_df <- data.frame(variable   = PREDICTORS,
                        importance = NA_real_,
                        stringsAsFactors = FALSE)

} else {

  vimp_vec <- rsf_model$variable.importance   # named numeric vector

  vimp_df <- data.frame(
    variable         = names(vimp_vec),
    importance       = as.numeric(vimp_vec),
    stringsAsFactors = FALSE
  ) %>% arrange(desc(importance))

  cat("\n  VIMP (permutation-based, drop in OOB C-statistic):\n")
  cat(sprintf("  %-5s %-25s %10s\n", "Rank", "Variable", "VIMP"))
  cat("  ", strrep("-", 43), "\n", sep = "")
  for (i in seq_len(nrow(vimp_df))) {
    cat(sprintf("  %-5d %-25s %10.4f\n",
                i, vimp_df$variable[i], vimp_df$importance[i]))
  }

  # Bar chart
  p_vimp <- ggplot(vimp_df,
                   aes(x    = reorder(variable, importance),
                       y    = importance,
                       fill = importance > 0)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, linewidth = 0.5, colour = "black") +
    scale_fill_manual(
      values = c("TRUE" = "#2166ac", "FALSE" = "#d73027"),
      guide  = "none"
    ) +
    coord_flip() +
    labs(title    = "RSF Variable Importance (VIMP)",
         subtitle = "Permutation importance: drop in OOB C-statistic",
         x = NULL, y = "VIMP") +
    theme_bw(base_size = 12)

  vimp_plot_path <- file.path(OUTPUT_DIR, "rsf_variable_importance.png")
  ggsave(vimp_plot_path, p_vimp, width = 8, height = 5, dpi = 150)
  cat(sprintf("\n  Plot saved to : %s\n", vimp_plot_path))
}

vimp_path <- file.path(OUTPUT_DIR, "rsf_variable_importance.csv")
write.csv(vimp_df, vimp_path, row.names = FALSE)
cat(sprintf("  Table saved to: %s\n", vimp_path))


# =============================================================================
# 8.  PREDICTED SURVIVAL CURVES BY RISK TERTILE
# =============================================================================
section("Survival curves by predicted risk tertile")

# Predict survival probabilities on training data
# type = "response" returns the full survival matrix (persons x time points)
cat("  Computing training predictions for survival curves...\n")
pred_train   <- predict(rsf_model, data = df_train, type = "response")
surv_train   <- pred_train$survival          # matrix: n_persons x n_timepoints
time_points  <- pred_train$unique.death.times

# Derive scalar risk score per person: mean cumulative hazard
risk_train <- survival_to_risk(surv_train)

# Assign each person to a risk tertile
risk_tertile <- cut(
  risk_train,
  breaks         = quantile(risk_train, probs = c(0, 1/3, 2/3, 1)),
  labels         = c("Low risk", "Medium risk", "High risk"),
  include.lowest = TRUE
)

# Average survival curve within each tertile
avg_surv_df <- bind_rows(lapply(levels(risk_tertile), function(grp) {
  idx <- which(risk_tertile == grp)
  data.frame(
    time             = time_points,
    survival         = colMeans(surv_train[idx, , drop = FALSE]),
    group            = grp,
    stringsAsFactors = FALSE
  )
}))

p_surv <- ggplot(avg_surv_df, aes(x = time, y = survival, colour = group)) +
  geom_step(linewidth = 1) +
  scale_colour_manual(
    values = c("Low risk"    = "#2166ac",
               "Medium risk" = "#f4a582",
               "High risk"   = "#d73027")
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(title    = "Predicted Survival Curves by Risk Tertile",
       subtitle = "Average RSF-predicted survival function within each group",
       x = "Month", y = "Probability of remaining housed", colour = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

surv_plot_path <- file.path(OUTPUT_DIR, "rsf_survival_curves.png")
ggsave(surv_plot_path, p_surv, width = 8, height = 5, dpi = 150)
cat(sprintf("  Plot saved to: %s\n", surv_plot_path))


# =============================================================================
# 9.  OUT-OF-SAMPLE PREDICTION (landmark approach, mirrors Figure 7)
# =============================================================================
section("Out-of-sample prediction")

# Step 1: snapshot of at-risk individuals at end of training period
landmark_month <- max(df_train_panel[[TIME_VAR]])
cat(sprintf("\n  Landmark month : %d (last month of training period)\n",
            landmark_month))

snap_df <- landmark_snapshot(df_train_panel, landmark_month, PREDICTORS)
cat(sprintf("  At-risk individuals at landmark : %s\n",
            format(nrow(snap_df), big.mark = ",")))

# Step 2: predict survival probabilities for each at-risk individual
cat("  Predicting...")
pred_oos  <- predict(rsf_model, data = snap_df, type = "response")
risk_oos  <- survival_to_risk(pred_oos$survival)
cat(" done.\n")

# Step 3: build results table
results_df <- data.frame(
  id             = snap_df$id,
  predicted_risk = risk_oos,
  stringsAsFactors = FALSE
)

# Step 4: attach actual outcomes from the test period
actual_test <- df_test_panel %>%
  group_by(.data[[ID_VAR]]) %>%
  summarise(became_homeless = as.integer(any(.data[[EVENT_VAR]] == 1)),
            .groups = "drop") %>%
  rename(id = .data[[ID_VAR]])

results_df <- results_df %>%
  left_join(actual_test, by = "id") %>%
  mutate(became_homeless = ifelse(is.na(became_homeless), 0L, became_homeless))

pop_rate     <- mean(results_df$became_homeless)
n_oos_events <- sum(results_df$became_homeless)

cat(sprintf("\n  OOS individuals : %s\n",   format(nrow(results_df), big.mark = ",")))
cat(sprintf("  OOS events      : %d (%.2f%% population rate)\n",
            n_oos_events, 100 * pop_rate))

# Step 5: C-statistic on OOS set
oos_conc <- survival::concordance(
  Surv(rep(HORIZON_MONTHS, nrow(results_df)), results_df$became_homeless) ~
    results_df$predicted_risk
)
cat(sprintf("\n  OOS C-statistic : %.4f\n", oos_conc$concordance))

# Step 6: top-percentile table (mirrors Figure 7 Panel a)
cat("\n  Top-percentile realised event rates:\n")
cat(sprintf("  %-12s  %-18s  %s\n", "Top group", "Realised rate", "vs. population"))
cat("  ", strrep("-", 50), "\n", sep = "")
for (pct in c(99, 95, 90, 75, 50)) {
  thresh  <- quantile(results_df$predicted_risk, pct / 100)
  top_grp <- results_df %>% filter(predicted_risk >= thresh)
  rate    <- mean(top_grp$became_homeless)
  ratio   <- ifelse(pop_rate > 0, rate / pop_rate, NA)
  cat(sprintf("  Top %2d%%        %18s  %.1fx\n",
              100 - pct, sprintf("%.1f%%", 100 * rate), ratio))
}

# Step 7: calibration plot (Figure 7 style)
results_df <- results_df %>% mutate(risk_bin = ntile(predicted_risk, 20))

cal_df <- results_df %>%
  group_by(risk_bin) %>%
  summarise(
    mean_predicted_risk = mean(predicted_risk),
    realised_rate       = mean(became_homeless),
    n                   = n(),
    .groups = "drop"
  )

p_cal <- ggplot(cal_df, aes(x = mean_predicted_risk, y = realised_rate)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey50", linewidth = 0.8) +
  geom_point(aes(size = n), colour = "#2166ac", alpha = 0.8) +
  geom_smooth(method = "loess", se = TRUE,
              colour = "#d73027", fill = "#fddbc7", linewidth = 0.8) +
  scale_size_continuous(range = c(2, 10), guide = "none") +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title    = "Out-of-sample Calibration: Predicted vs. Realised Risk",
    subtitle = sprintf(
      "Landmark month %d  |  Horizon %d months  |  20 bins  |  bubble size proportional to n",
      landmark_month, HORIZON_MONTHS),
    x = "Mean predicted risk score (mean cumulative hazard)",
    y = "Realised homelessness entry rate"
  ) +
  theme_bw(base_size = 12)

cal_plot_path <- file.path(OUTPUT_DIR, "rsf_calibration_plot.png")
ggsave(cal_plot_path, p_cal, width = 8, height = 6, dpi = 150)
cat(sprintf("\n  Calibration plot saved to: %s\n", cal_plot_path))

oos_path <- file.path(OUTPUT_DIR, "rsf_oos_predictions.csv")
write.csv(results_df, oos_path, row.names = FALSE)
cat(sprintf("  OOS predictions saved to : %s\n", oos_path))


# =============================================================================
# 10.  OUT-OF-SAMPLE PREDICTION (external file, optional)
# =============================================================================
if (!is.null(OOS_FILE)) {
  section(paste("External OOS prediction:", OOS_FILE))

  df_ext   <- load_panel(OOS_FILE)
  snap_ext <- landmark_snapshot(df_ext, max(df_ext[[TIME_VAR]]), PREDICTORS)
  pred_ext <- predict(rsf_model, data = snap_ext, type = "response")

  ext_results <- data.frame(
    id             = snap_ext$id,
    predicted_risk = survival_to_risk(pred_ext$survival),
    stringsAsFactors = FALSE
  )
  ext_path <- file.path(OUTPUT_DIR, "rsf_ext_oos_predictions.csv")
  write.csv(ext_results, ext_path, row.names = FALSE)
  cat(sprintf("\n  External OOS predictions saved to: %s\n", ext_path))
}


# =============================================================================
# 11.  SAVE MODEL AND DIAGNOSTICS
# =============================================================================
section("Saving model and diagnostics")

model_path <- file.path(OUTPUT_DIR, "rsf_model.rds")
saveRDS(rsf_model, model_path)
cat(sprintf("\n  Model saved to : %s\n", model_path))
cat(sprintf("  Reload with   : rsf_model <- readRDS(\"%s\")\n", model_path))

# Write diagnostics to text file
# tryCatch guarantees sink() is always closed even if an error occurs,
# which would otherwise permanently redirect your console to the file.
diag_path <- file.path(OUTPUT_DIR, "rsf_diagnostics.txt")
tryCatch({
  sink(diag_path)
  cat("Random Survival Forest (ranger) -- Homelessness Prediction\n")
  cat(strrep("=", 60), "\n\n")
  cat("Training file   :", TRAIN_FILE, "\n")
  cat("Training period : years <", CUTOFF_YEAR, "\n")
  cat("Individuals     :", n_train,  "\n")
  cat("Events          :", n_events, "\n")
  cat("Predictors      :", paste(PREDICTORS, collapse = ", "), "\n\n")
  cat("ranger parameters\n")
  cat("  num.trees     :", N_TREES,    "\n")
  cat("  min.node.size :", NODE_SIZE,  "\n")
  cat("  mtry          :", mtry_val,   "\n\n")
  cat("OOB C-statistic :", round(oob_cstat,            4), "\n")
  cat("OOS C-statistic :", round(oos_conc$concordance,  4), "\n\n")
  if (COMPUTE_VIMP) {
    cat("Variable importance (VIMP):\n")
    print(vimp_df)
  } else {
    cat("Variable importance: not computed (set COMPUTE_VIMP <- TRUE)\n")
  }
  sink()
}, error = function(e) {
  if (sink.number() > 0) sink()
  stop(e)
})
cat(sprintf("  Diagnostics saved to: %s\n", diag_path))


# =============================================================================
# 12.  SUMMARY
# =============================================================================
section("Summary")
cat(sprintf("
  Model          : Random Survival Forest (ranger)
  N individuals  : %s
  Events         : %d (%.1f%%)
  Predictors     : %d
  Trees          : %d  |  min.node.size : %d  |  mtry : %d

  OOB C-statistic : %.4f
  OOS C-statistic : %.4f

  Output files
    Model       : %s
    VIMP table  : %s
    VIMP plot   : %s
    Surv curves : %s
    Calibration : %s
    OOS preds   : %s
    Diagnostics : %s

Done.
",
  format(n_train, big.mark = ","),
  n_events, 100 * n_events / n_train,
  length(PREDICTORS),
  N_TREES, NODE_SIZE, mtry_val,
  oob_cstat,
  oos_conc$concordance,
  model_path,
  vimp_path,
  ifelse(COMPUTE_VIMP, vimp_plot_path, "skipped (COMPUTE_VIMP = FALSE)"),
  surv_plot_path,
  cal_plot_path,
  oos_path,
  diag_path
))
