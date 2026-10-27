## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  eval = requireNamespace("ranger", quietly = TRUE) &&
    requireNamespace("tidyr", quietly = TRUE),
  collapse = TRUE,
  message = FALSE,
  warning = FALSE,
  comment = "#>"
)

## ----data-setup---------------------------------------------------------------
library(pintervals)
library(dplyr)
library(tidyr)
library(ranger)
data("county_turnout_panel", package = "pintervals")

set.seed(101010)

train <- county_turnout_panel %>%
  filter(election_year < 2020) %>%
  mutate(
    rel_dem_share = dem_share - dem_share_electorate
  ) %>%
  na.omit()

calib <- county_turnout_panel %>%
  filter(election_year == 2020) %>%
  mutate(
    rel_dem_share = dem_share - dem_share_electorate
  ) %>%
  na.omit()

test <- county_turnout_panel %>%
  filter(election_year == 2024) %>%
  mutate(
    rel_dem_share = dem_share - dem_share_electorate
  ) %>%
  na.omit()

m <- ranger(
  rel_dem_share ~ nonwhite_pct + foreignborn_pct + female_pct +
    age29andunder_pct + median_hh_inc + clf_unemploy_pct +
    lesscollege_pct + rural_pct,
  data = train
)

calib$predicted <- predict(m, data = calib)$predictions
test$predicted <- predict(m, data = test)$predictions

## ----conformal-intervals------------------------------------------------------
conformal_intervals <- pinterval_conformal(
  pred = test$predicted,
  calib = calib$predicted,
  calib_truth = calib$rel_dem_share,
  alpha = 0.1
)

conformal_coverage <- interval_coverage(
  truth = test$rel_dem_share,
  lower_bound = conformal_intervals$lower_bound,
  upper_bound = conformal_intervals$upper_bound
)
conformal_coverage

## ----parametric-bootstrap-intervals-------------------------------------------
norm_intervals <- pinterval_parametric(
  pred = test$predicted,
  calib = calib$predicted,
  calib_truth = calib$rel_dem_share,
  dist = "norm",
  alpha = 0.1
)

norm_coverage <- interval_coverage(
  truth = test$rel_dem_share,
  lower_bound = norm_intervals$lower_bound,
  upper_bound = norm_intervals$upper_bound
)

bootstrap_intervals <- pinterval_bootstrap(
  pred = test$predicted,
  calib = calib$predicted,
  calib_truth = calib$rel_dem_share,
  alpha = 0.1,
  n_bootstrap = 1000
)

bootstrap_coverage <- interval_coverage(
  truth = test$rel_dem_share,
  lower_bound = bootstrap_intervals$lower_bound,
  upper_bound = bootstrap_intervals$upper_bound
)

c(norm_coverage, bootstrap_coverage)

## ----mondrian-ccp-dwcp--------------------------------------------------------
mondrian_intervals <- pinterval_mondrian(
  pred = test$predicted,
  pred_class = test$region,
  calib = calib$predicted,
  calib_truth = calib$rel_dem_share,
  calib_class = calib$region,
  alpha = 0.1
)

mondrian_coverage <- interval_coverage(
  truth = test$rel_dem_share,
  lower_bound = mondrian_intervals$lower_bound,
  upper_bound = mondrian_intervals$upper_bound
)

clustered_conformal_intervals <- pinterval_ccp(
  pred = test$predicted,
  pred_class = test$division,
  calib = calib$predicted,
  calib_truth = calib$rel_dem_share,
  calib_class = calib$division,
  alpha = 0.1,
  optimize_n_clusters = TRUE,
  max_n_clusters = 5
)

clustered_conformal_coverage <- interval_coverage(
  truth = test$rel_dem_share,
  lower_bound = clustered_conformal_intervals$lower_bound,
  upper_bound = clustered_conformal_intervals$upper_bound
)

distance_weighted_cp_intervals <- pinterval_conformal(
  pred = test$predicted,
  calib = calib$predicted,
  calib_truth = calib$rel_dem_share,
  alpha = 0.1,
  distance_weighted_cp = TRUE,
  distance_features_calib = calib %>%
    select(latitude, longitude),
  distance_features_pred = test %>%
    select(latitude, longitude),
  distance_type = "euclidean",
  normalize_distance = "sd"
)

distance_weighted_cp_coverage <- interval_coverage(
  truth = test$rel_dem_share,
  lower_bound = distance_weighted_cp_intervals$lower_bound,
  upper_bound = distance_weighted_cp_intervals$upper_bound
)

c(mondrian_coverage, clustered_conformal_coverage,
  distance_weighted_cp_coverage)

## ----group-wise-coverage------------------------------------------------------
conformal_intervals <- conformal_intervals %>%
  mutate(
    region = test$region,
    method = "SCP",
    rel_dem_share = test$rel_dem_share
  )

mondrian_intervals <- mondrian_intervals %>%
  mutate(
    region = test$region,
    method = "MCP",
    rel_dem_share = test$rel_dem_share
  )

clustered_conformal_intervals <- clustered_conformal_intervals %>%
  mutate(
    region = test$region,
    method = "CCP",
    rel_dem_share = test$rel_dem_share
  )

distance_weighted_cp_intervals <- distance_weighted_cp_intervals %>%
  mutate(
    region = test$region,
    method = "DWCP",
    rel_dem_share = test$rel_dem_share
  )

bootstrap_intervals <- bootstrap_intervals %>%
  mutate(
    region = test$region,
    method = "Bootstrap",
    rel_dem_share = test$rel_dem_share
  )

norm_intervals <- norm_intervals %>%
  mutate(
    region = test$region,
    method = "Normal",
    rel_dem_share = test$rel_dem_share
  )

all_conformal_intervals <- bind_rows(
  conformal_intervals,
  mondrian_intervals,
  clustered_conformal_intervals,
  distance_weighted_cp_intervals,
  bootstrap_intervals,
  norm_intervals
)

group_wise_coverage <- all_conformal_intervals %>%
  group_by(method, region) %>%
  summarise(
    coverage = interval_coverage(
      truth = rel_dem_share,
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
mae_coverage_regions <- group_wise_coverage %>%
  group_by(method) %>%
  summarize(mae_coverage = mean(abs(coverage - 0.9))) %>%
  pivot_wider(
    names_from = method,
    values_from = mae_coverage
  )
mae_coverage_regions

## ----bccp-intervals-----------------------------------------------------------
calib <- calib %>%
  mutate(
    truth_bin = case_when(
      rel_dem_share < -0.3 ~ 1,
      rel_dem_share < -0.2 ~ 2,
      rel_dem_share < -0.1 ~ 3,
      TRUE ~ 4
    )
  )

bccp_contiguized_intervals <- pinterval_bccp(
  pred = test$predicted,
  calib = calib$predicted,
  calib_truth = calib$rel_dem_share,
  calib_bins = calib$truth_bin,
  alpha = 0.1,
  contiguize = TRUE
)

bccp_contiguized_coverage <- interval_coverage(
  truth = test$rel_dem_share,
  lower_bound = bccp_contiguized_intervals$lower_bound,
  upper_bound = bccp_contiguized_intervals$upper_bound
)

bccp_discontigous_intervals <- pinterval_bccp(
  pred = test$predicted,
  calib = calib$predicted,
  calib_truth = calib$rel_dem_share,
  calib_bins = calib$truth_bin,
  alpha = 0.1,
  contiguize = FALSE
)

bccp_discontigous_coverage <- interval_coverage(
  truth = test$rel_dem_share,
  intervals = bccp_discontigous_intervals$intervals
)

c(bccp_contiguized_coverage, bccp_discontigous_coverage)

## ----bin-wise-coverage--------------------------------------------------------
test <- test %>%
  mutate(
    truth_bin = case_when(
      rel_dem_share < -0.3 ~ 1,
      rel_dem_share < -0.2 ~ 2,
      rel_dem_share < -0.1 ~ 3,
      TRUE ~ 4
    )
  )

bccp_contiguized_intervals <- bccp_contiguized_intervals %>%
  mutate(
    truth_bin = test$truth_bin,
    method = "BCCP(c)",
    rel_dem_share = test$rel_dem_share
  )

bccp_discontigous_intervals <- bccp_discontigous_intervals %>%
  mutate(
    truth_bin = test$truth_bin,
    method = "BCCP(d)",
    rel_dem_share = test$rel_dem_share
  )

conformal_intervals <- conformal_intervals %>%
  mutate(
    truth_bin = test$truth_bin,
    method = "SCP",
    rel_dem_share = test$rel_dem_share
  )

bootstrap_intervals <- bootstrap_intervals %>%
  mutate(
    truth_bin = test$truth_bin,
    method = "Bootstrap",
    rel_dem_share = test$rel_dem_share
  )

norm_intervals <- norm_intervals %>%
  mutate(
    truth_bin = test$truth_bin,
    method = "Normal",
    rel_dem_share = test$rel_dem_share
  )

all_intervals_bins <- bind_rows(
  bccp_contiguized_intervals,
  bccp_discontigous_intervals,
  conformal_intervals,
  bootstrap_intervals,
  norm_intervals
)

bin_wise_coverage <- all_intervals_bins %>%
  group_by(method, truth_bin) %>%
  summarise(
    coverage = interval_coverage(
      truth = rel_dem_share,
      lower_bound = lower_bound,
      upper_bound = upper_bound,
      intervals = intervals
    ),
    .groups = "drop"
  )

bin_wise_coverage_wide <- bin_wise_coverage %>%
  pivot_wider(
    names_from = truth_bin,
    values_from = coverage
  )
bin_wise_coverage_wide

## ----mae-coverage-bins--------------------------------------------------------
mae_coverage_bins <- bin_wise_coverage %>%
  group_by(method) %>%
  summarize(mae_coverage = mean(abs(coverage - 0.9))) %>%
  pivot_wider(
    names_from = method,
    values_from = mae_coverage
  )
mae_coverage_bins

## ----score-scp-norm-boot------------------------------------------------------
conformal_score <- interval_score(
  truth = test$rel_dem_share,
  lower_bound = conformal_intervals$lower_bound,
  upper_bound = conformal_intervals$upper_bound,
  alpha = 0.1
)

norm_score <- interval_score(
  truth = test$rel_dem_share,
  lower_bound = norm_intervals$lower_bound,
  upper_bound = norm_intervals$upper_bound,
  alpha = 0.1
)

bootstrap_score <- interval_score(
  truth = test$rel_dem_share,
  lower_bound = bootstrap_intervals$lower_bound,
  upper_bound = bootstrap_intervals$upper_bound,
  alpha = 0.1
)

c(conformal_score, norm_score, bootstrap_score)

## ----width-scp-norm-boot------------------------------------------------------
conformal_width <- interval_width(
  lower_bound = conformal_intervals$lower_bound,
  upper_bound = conformal_intervals$upper_bound
)

norm_width <- interval_width(
  lower_bound = norm_intervals$lower_bound,
  upper_bound = norm_intervals$upper_bound
)

bootstrap_width <- interval_width(
  lower_bound = bootstrap_intervals$lower_bound,
  upper_bound = bootstrap_intervals$upper_bound
)

c(conformal_width, norm_width, bootstrap_width)

## ----score-mcp-ccp-dwcp-------------------------------------------------------
mondrian_score <- interval_score(
  truth = test$rel_dem_share,
  lower_bound = mondrian_intervals$lower_bound,
  upper_bound = mondrian_intervals$upper_bound,
  alpha = 0.1
)

clustered_conformal_score <- interval_score(
  truth = test$rel_dem_share,
  lower_bound = clustered_conformal_intervals$lower_bound,
  upper_bound = clustered_conformal_intervals$upper_bound,
  alpha = 0.1
)

distance_weighted_cp_score <- interval_score(
  truth = test$rel_dem_share,
  lower_bound = distance_weighted_cp_intervals$lower_bound,
  upper_bound = distance_weighted_cp_intervals$upper_bound,
  alpha = 0.1
)

c(mondrian_score, clustered_conformal_score,
  distance_weighted_cp_score)

## ----width-mcp-ccp-dwcp-------------------------------------------------------
mondrian_width <- interval_width(
  lower_bound = mondrian_intervals$lower_bound,
  upper_bound = mondrian_intervals$upper_bound
)

clustered_conformal_width <- interval_width(
  lower_bound = clustered_conformal_intervals$lower_bound,
  upper_bound = clustered_conformal_intervals$upper_bound
)

distance_weighted_cp_width <- interval_width(
  lower_bound = distance_weighted_cp_intervals$lower_bound,
  upper_bound = distance_weighted_cp_intervals$upper_bound
)

c(mondrian_width, clustered_conformal_width,
  distance_weighted_cp_width)

## ----group-wise-score---------------------------------------------------------
group_wise_score <- all_conformal_intervals %>%
  group_by(method, region) %>%
  summarise(
    score = interval_score(
      truth = rel_dem_share,
      lower_bound = lower_bound,
      upper_bound = upper_bound,
      alpha = 0.1
    ),
    width = interval_width(
      lower_bound = lower_bound,
      upper_bound = upper_bound
    ),
    .groups = "drop"
  )

group_wise_score %>%
  group_by(method) %>%
  mutate(
    mean_score = mean(score),
    mean_width = mean(width)
  ) %>%
  ungroup() %>%
  select(-width) %>%
  pivot_wider(
    names_from = region,
    values_from = score
  )

## ----score-bccp---------------------------------------------------------------
bccp_contiguized_score <- interval_score(
  truth = test$rel_dem_share,
  lower_bound = bccp_contiguized_intervals$lower_bound,
  upper_bound = bccp_contiguized_intervals$upper_bound,
  alpha = 0.1
)

bccp_discontigous_score <- interval_score(
  truth = test$rel_dem_share,
  intervals = bccp_discontigous_intervals$intervals,
  alpha = 0.1
)

c(bccp_contiguized_score, bccp_discontigous_score)

## ----width-bccp---------------------------------------------------------------
bccp_contiguized_width <- interval_width(
  lower_bound = bccp_contiguized_intervals$lower_bound,
  upper_bound = bccp_contiguized_intervals$upper_bound
)

bccp_discontigous_width <- interval_width(
  intervals = bccp_discontigous_intervals$intervals
)

c(bccp_contiguized_width, bccp_discontigous_width)

## ----bin-wise-score-----------------------------------------------------------
bin_wise_score <- all_intervals_bins %>%
  group_by(method, truth_bin) %>%
  summarise(
    score = interval_score(
      truth = rel_dem_share,
      lower_bound = lower_bound,
      upper_bound = upper_bound,
      intervals = intervals,
      alpha = 0.1
    ),
    width = interval_width(
      lower_bound = lower_bound,
      upper_bound = upper_bound,
      intervals = intervals
    ),
    .groups = "drop"
  )

bin_wise_score %>%
  group_by(method) %>%
  mutate(
    mean_score = mean(score),
    mean_width = mean(width)
  ) %>%
  ungroup() %>%
  select(-width) %>%
  pivot_wider(
    names_from = truth_bin,
    values_from = score
  )

