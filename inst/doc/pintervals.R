## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  eval = TRUE,
  collapse = TRUE,
  comment = "#>"
)

## ----data-setup---------------------------------------------------------------
library(pintervals)
library(dplyr)
library(tidyr)
data("county_turnout", package = "pintervals")

# Split data into calibration and test sets
set.seed(20260106)
nobs <- nrow(county_turnout)
calib_indices <- sample(nobs, size = 0.5 * nobs)
calib_data <- county_turnout[calib_indices, ]
test_data <- county_turnout[-calib_indices, ]

## ----conformal-intervals------------------------------------------------------
# Standard CP intervals
conformal_intervals <- pinterval_conformal(
  pred = test_data$predicted_turnout,
  calib = calib_data$predicted_turnout,
  calib_truth = calib_data$turnout,
  alpha = 0.1
)

# Evaluate coverage and width
conformal_coverage <- interval_coverage(
  truth = test_data$turnout,
  lower_bound = conformal_intervals$lower_bound,
  upper_bound = conformal_intervals$upper_bound
)
conformal_coverage

## ----parametric-bootstrap-intervals-------------------------------------------
# Parametric prediction intervals assuming normal distribution
norm_intervals <- pinterval_parametric(
  pred = test_data$predicted_turnout,
  calib = calib_data$predicted_turnout,
  calib_truth = calib_data$turnout,
  dist = "norm",
  alpha = 0.1
)

# Evaluate coverage
norm_coverage <- interval_coverage(
  truth = test_data$turnout,
  lower_bound = norm_intervals$lower_bound,
  upper_bound = norm_intervals$upper_bound
)

# Parametric prediction intervals assuming logistic distribution
logis_intervals <- pinterval_parametric(
  pred = test_data$predicted_turnout,
  calib = calib_data$predicted_turnout,
  calib_truth = calib_data$turnout,
  dist = "logis",
  alpha = 0.1
)

# Evaluate coverage
logis_coverage <- interval_coverage(
  truth = test_data$turnout,
  lower_bound = logis_intervals$lower_bound,
  upper_bound = logis_intervals$upper_bound
)

# Bootstrapped prediction intervals
bootstrap_intervals <- pinterval_bootstrap(
  pred = test_data$predicted_turnout,
  calib = calib_data$predicted_turnout,
  calib_truth = calib_data$turnout,
  alpha = 0.1,
  n_bootstrap = 1000
)

# Evaluate coverage
bootstrap_coverage <- interval_coverage(
  truth = test_data$turnout,
  lower_bound = bootstrap_intervals$lower_bound,
  upper_bound = bootstrap_intervals$upper_bound
)

c(norm_coverage, logis_coverage, bootstrap_coverage)

## ----mondrian-ccp-dwcp--------------------------------------------------------
# MCP intervals
mondrian_intervals <- pinterval_mondrian(
  pred = test_data$predicted_turnout,
  pred_class = test_data$region,
  calib = calib_data$predicted_turnout,
  calib_truth = calib_data$turnout,
  calib_class = calib_data$region,
  alpha = 0.1
)

# Evaluate coverage
mondrian_coverage <- interval_coverage(
  truth = test_data$turnout,
  lower_bound = mondrian_intervals$lower_bound,
  upper_bound = mondrian_intervals$upper_bound
)

# Clustered CP intervals
clustered_conformal_intervals <- pinterval_ccp(
  pred = test_data$predicted_turnout,
  pred_class = test_data$division,
  calib = calib_data$predicted_turnout,
  calib_truth = calib_data$turnout,
  calib_class = calib_data$division,
  alpha = 0.1,
  optimize_n_clusters = TRUE,
  max_n_clusters = 5
)

# Evaluate coverage
clustered_conformal_coverage <- interval_coverage(
  truth = test_data$turnout,
  lower_bound = clustered_conformal_intervals$lower_bound,
  upper_bound = clustered_conformal_intervals$upper_bound
)

# Distance-weighted CP intervals using geographic distances
dw_cp_intervals <- pinterval_conformal(
  pred = test_data$predicted_turnout,
  calib = calib_data$predicted_turnout,
  calib_truth = calib_data$turnout,
  alpha = 0.1,
  distance_weighted_cp = TRUE,
  distance_features_calib = calib_data %>%
    select(latitude, longitude),
  distance_features_pred = test_data %>%
    select(latitude, longitude),
  normalize_distance = "sd"
)

# Evaluate coverage
dw_cp_coverage <- interval_coverage(
  truth = test_data$turnout,
  lower_bound = dw_cp_intervals$lower_bound,
  upper_bound = dw_cp_intervals$upper_bound
)

c(mondrian_coverage, clustered_conformal_coverage, dw_cp_coverage)

## ----group-wise-coverage------------------------------------------------------
# Adding grouping variable and true turnout values
conformal_intervals <- conformal_intervals %>%
  mutate(region = test_data$region,
         method = "SCP",
         turnout = test_data$turnout)

mondrian_intervals <- mondrian_intervals %>%
  mutate(region = test_data$region,
         method = "MCP",
         turnout = test_data$turnout)

clustered_conformal_intervals <- clustered_conformal_intervals %>%
  mutate(region = test_data$region,
         method = "CCP",
         turnout = test_data$turnout)

dw_cp_intervals <- dw_cp_intervals %>%
  mutate(region = test_data$region,
         method = "DWCP",
         turnout = test_data$turnout)

bootstrap_intervals <- bootstrap_intervals %>%
  mutate(region = test_data$region,
         method = "Bootstrap",
         turnout = test_data$turnout)

norm_intervals <- norm_intervals %>%
  mutate(region = test_data$region,
         method = "Normal",
         turnout = test_data$turnout)

logis_intervals <- logis_intervals %>%
  mutate(region = test_data$region,
         method = "Logistic",
         turnout = test_data$turnout)

# Combine all intervals for group-wise coverage evaluation
all_conformal_intervals <- bind_rows(
  conformal_intervals,
  mondrian_intervals,
  clustered_conformal_intervals,
  dw_cp_intervals,
  bootstrap_intervals,
  norm_intervals,
  logis_intervals
)

# Evaluate group-wise coverage
group_wise_coverage <- all_conformal_intervals %>%
  group_by(method, region) %>%
  summarise(
    coverage = interval_coverage(
      truth = turnout,
      lower_bound = lower_bound,
      upper_bound = upper_bound
    ),
    .groups = "drop"
  )

group_wise_coverage %>%
  pivot_wider(
    names_from = region,
    values_from = coverage
  )

## ----mae-coverage-regions-----------------------------------------------------
group_wise_coverage %>%
  group_by(method) %>%
  summarize(mae_coverage = mean(abs(coverage - 0.9))) %>%
  pivot_wider(
    names_from = method,
    values_from = mae_coverage
  )

## ----bccp-intervals-----------------------------------------------------------
# Create bins of turnout rates in the calibration set
calib_data <- calib_data %>%
  mutate(turnout_bin = case_when(
    turnout < 0.5 ~ 1,
    turnout < 0.6 ~ 2,
    turnout < 0.65 ~ 3,
    TRUE ~ 4
  ))

bccp_contiguized_intervals <- pinterval_bccp(
  pred = test_data$predicted_turnout,
  calib = calib_data$predicted_turnout,
  calib_truth = calib_data$turnout,
  calib_bins = calib_data$turnout_bin,
  alpha = 0.1,
  contiguize = TRUE
)

# Evaluate coverage
bccp_contiguized_coverage <- interval_coverage(
  truth = test_data$turnout,
  lower_bound = bccp_contiguized_intervals$lower_bound,
  upper_bound = bccp_contiguized_intervals$upper_bound
)

bccp_discontigous_intervals <- pinterval_bccp(
  pred = test_data$predicted_turnout,
  calib = calib_data$predicted_turnout,
  calib_truth = calib_data$turnout,
  calib_bins = calib_data$turnout_bin,
  alpha = 0.1,
  contiguize = FALSE
)

# Evaluate coverage
bccp_discontigous_coverage <- interval_coverage(
  truth = test_data$turnout,
  intervals = bccp_discontigous_intervals$intervals
)

c(bccp_contiguized_coverage, bccp_discontigous_coverage)

## ----bin-wise-coverage--------------------------------------------------------
# Adding bin variable to the intervals
test_data <- test_data %>%
  mutate(turnout_bin = case_when(
    turnout < 0.5 ~ 1,
    turnout < 0.6 ~ 2,
    turnout < 0.65 ~ 3,
    TRUE ~ 4
  ))

bccp_contiguized_intervals <- bccp_contiguized_intervals %>%
  mutate(turnout_bin = test_data$turnout_bin,
         method = "BCCP(c)",
         turnout = test_data$turnout)

bccp_discontigous_intervals <- bccp_discontigous_intervals %>%
  mutate(turnout_bin = test_data$turnout_bin,
         method = "BCCP(d)",
         turnout = test_data$turnout)

conformal_intervals <- conformal_intervals %>%
  mutate(turnout_bin = test_data$turnout_bin)

bootstrap_intervals <- bootstrap_intervals %>%
  mutate(turnout_bin = test_data$turnout_bin)

norm_intervals <- norm_intervals %>%
  mutate(turnout_bin = test_data$turnout_bin)

logis_intervals <- logis_intervals %>%
  mutate(turnout_bin = test_data$turnout_bin)

# Combine all intervals for bin-wise coverage evaluation
all_intervals_bins <- bind_rows(
  bccp_contiguized_intervals,
  bccp_discontigous_intervals,
  conformal_intervals,
  bootstrap_intervals,
  norm_intervals,
  logis_intervals
)

# Evaluate bin-wise coverage
bin_wise_coverage <- all_intervals_bins %>%
  group_by(method, turnout_bin) %>%
  summarise(
    coverage = interval_coverage(
      truth = turnout,
      lower_bound = lower_bound,
      upper_bound = upper_bound,
      intervals = intervals
    ),
    .groups = "drop"
  )

bin_wise_coverage %>%
  pivot_wider(
    names_from = turnout_bin,
    values_from = coverage
  )

## ----mae-coverage-bins--------------------------------------------------------
bin_wise_coverage %>%
  group_by(method) %>%
  summarize(mae_coverage = mean(abs(coverage - 0.9))) %>%
  pivot_wider(
    names_from = method,
    values_from = mae_coverage
  )

