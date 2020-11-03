#!/usr/bin/env Rscript

"Usage:
  cgd_fit.R [options]
  
Options:
  --num-runs=<runs>  Number of runs per month [default: 1]
  --num-months=<months>  Number of months to fit with the CGD data [default: 12]
  --output=<run-data-file>  Where to store output [default: cgd_optim.rds]
  --append-stage-2=<run-data-file>
  --no-vcov-moments
  --vcov-moments-weight=<weight>  Weighting matrix weight for covariance moments [default: 1]
  --parallel=<cores>
" -> opt_desc

script_options <- if (interactive()) {
  docopt::docopt(opt_desc, '--output=test.rds --num-runs=1')
} else {
  docopt::docopt(opt_desc)
}

library(data.table)
library(magrittr)
library(tidyverse)
library(readxl)
library(gmm)
library(furrr)
library(R6)

library(vaccineEarlyInvest)

source("model_fit_util.R")

script_options %<>% 
  modify_at(c("num_runs", "num_months"), as.integer) %>% 
  modify_at(c("vcov_moments_weight"), as.numeric)

# Load the CGD data and convert it ----------------------------------------

zero_month_date <- lubridate::as_date("2020-10-01") 
start_month_offset <- 9

raw_cgd_master_input <- read_xlsx(file.path("data", "Master Input Data - COVID-19 Vaccine Predictions publication-2.xlsx")) %>% 
  pivot_longer(
    cols = `Phase 1 country`:`Phase 3 trial number`, 
    names_to = c("phase", ".value"), 
    names_pattern = "Phase (\\d) ((?:(?:start|end) date)|(?:country)|(?:trial number))",
    names_transform = list(
      phase = as.integer
    ),
    values_drop_na = TRUE
  ) %>% 
  rename(
    phase_country = "country",
    phase_start_date = "start date",
    phase_end_date = "end date",
    phase_trial_number = "trial number") %>% 
  mutate_at(vars(starts_with("phase_")), ~ na_if(.x, "NA")) %>% 
  nest(phase_data = phase:phase_trial_number) %>% 
  mutate(phase_data = map(phase_data, filter, across(starts_with("phase_"), ~ !is.na(.x)))) 

cgd_master_input <- raw_cgd_master_input %>% 
  filter(!str_detect(Institutes, "Kentucky Bioprocessing")) %>% # Exists in out data as pre-clinical but also Unknown/Unknown
  transmute(
    cgd_vaccine_id = Number,
    Platform = str_replace(Platform, ".+viral vector", "Viral vector"), # Not distinguishing between replicating and non-replicating
    Subcategory,
    phase = map_chr(phase_data, ~ if (nrow(.x) > 0) str_c("Phase ", max(.x$phase)) else "Pre-clinical")
  )

cgd_trials <- read_csv(
  file.path("data", "cgd_trials.csv"), 
  skip = 1,
  col_names = c("try_id", "vaccine_id", "phase_mon_1", "phase_mon_2", "phase_mon_3", "phase_mon_approval"), 
) %>% 
  mutate_all(as.integer) 

# Settings ----------------------------------------------------------------

group_vaccines_by <- NULL 
maxcand <- 50 
candidate_data <- load_candidate_data(file.path("data", "vaccinesSummaryCGD.csv"))
dordered <- calc_dordered(candidate_data, maxcand = maxcand)

cgd_id_dict <- if (FALSE) {
  get_candidate_draws(
    candidate_data, replications = 1, dordered = dordered,
      maxcand = maxcand,
      group_vaccines_by = group_vaccines_by
    ) %>% 
    select(candInd, Platform, Subcategory, phase) %>% 
    nest(ids = candInd) %>% 
    right_join(
      cgd_master_input %>% 
        filter(fct_match(phase, c("Phase 2", "Phase 3"))) %>% 
        nest(cgd_ids = cgd_vaccine_id),
      by = c("Platform", "Subcategory", "phase")
    ) %>% 
    mutate(
      ids = map2(ids, cgd_ids, ~ bind_cols(.x, sample_frac(.y))) # Randomly match IDs
    ) %>% 
    unnest(ids)
} else {
  read_rds(file.path("data", "cgd_id_dict.rds"))
}

cgd_draws_month <- lubridate::as_date("2020/10/01")
fit_month_offsets <- seq(script_options$num_months) + 8

initial_par <- c(poverall = 0.5, psubcat = 0.5, 
                 pvector = 0.5, psubunit = 0.5, prna = 0.5, pdna = 0.5, pinactivated = 0.5, pattenuated = 0.5,
                 pphase2 = 0.5, pphase3 = 0.5)

# CGD Fit -----------------------------------------------------------------

cgd_data <- tibble(
  month_offset = fit_month_offsets,
  month = cgd_draws_month + months(fit_month_offsets)
) %>% 
  rowwise() %>% 
  mutate(
    cgd_converted = list(convert_cgd_trials(month_offset, cgd_trials, cgd_id_dict)) 
  ) %>% 
  ungroup()

cgd_optim_data <- cgd_data %>%  
  mutate(
    cgd_summary = future_map(cgd_converted, summarize_draws)
  ) %>% 
  select(-cgd_converted)

num_success_rates_moments <- cgd_optim_data$cgd_summary[[1]]$success_rates %>% nrow()
num_vcov_moments <- cgd_optim_data$cgd_summary[[1]]$success_vcov %>% length()

weighting_matrix <- if (!script_options$no_vcov_moments) 
  diag(c(rep(1, num_success_rates_moments), rep(script_options$vcov_moments_weight, num_vcov_moments))) 

cgd_optim_data %<>% 
  modelr::data_grid(run_id = seq(script_options$num_runs), month_offset) %>% 
  left_join(cgd_optim_data, "month_offset") 

if (is_null(script_options$append_stage_2)) {
  cat("First stage...")
  
  cgd_optim_data %<>% 
    mutate(param_data = future_map(
      cgd_summary,
      function(summary) {
        optim_run(
          NULL, 
          summary, 
          initial_par = initial_par,
          use_vcov_moments = !script_options$no_vcov_moments,
          weighting_matrix = weighting_matrix,
          maxcand = maxcand)
      },
      .progress = TRUE
    )
  ) 
  
  write_rds(cgd_optim_data, file.path("data", script_options$output))
  
  cat("done.\n")
} else {
  cgd_optim_data <- read_rds(script_options$append_stage_2)
}

cat("Second stage...")

cgd_optim_data %<>% 
  mutate(
    param_data = future_map2(
      cgd_summary, param_data,
      function(summary, prev_run_data) {
        optim_run(
          NULL,
          summary,
          prev_run_data = prev_run_data,
          initial_par = initial_par,
          use_vcov_moments = !script_options$no_vcov_moments,
          weighting_matrix = weighting_matrix,
          maxcand = maxcand,
          replications = 100e3, 
          ndeps = rep_along(initial_par, 2e-3)
        )
      },
      .progress = TRUE),
  ) %>% 
  rowwise() %>% 
  mutate(run_summary = list(last(param_data$summary))) %>% 
  ungroup()

write_rds(cgd_optim_data, file.path("data", script_options$output))

cat("...done.\n")