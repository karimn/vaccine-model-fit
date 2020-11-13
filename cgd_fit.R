#!/usr/bin/env Rscript

"Usage:
  cgd_fit.R [options]
  
Options:
  --num-runs=<runs>  Number of runs per month [default: 1]
  --num-months=<months>  Number of months to fit with the CGD data [default: 12]
  --output=<run-data-file>  Where to store output [default: cgd_optim]
  --append-stage-2=<run-data-file>
  --no-vcov-moments
  --vcov-moments-weight=<weight>  Weighting matrix weight for covariance moments [default: 1]
  --parallel=<cores>  In Rscript this is automatically used by {future} to set up a multicore run.
  --pfizer-pref  Half of the CGD sample should have the Pfizer vaccine succeed before March 2021.
" -> opt_desc

script_options <- if (interactive()) {
  docopt::docopt(opt_desc, '--output=test --num-runs=1 --num-months=1')
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

# Load AHT summary data ---------------------------------------------------

long_aht_summary <- read_csv(file.path("data", "vaccinesSummaryCGD.csv")) %>% 
  pivot_longer(Phase1candidates:RepurposedCandidates, names_to = "phase", values_to = "n") %>% 
  mutate(n = as.integer(n))

# Load the CGD data ----------------------------------------

cgd_draws_date <- lubridate::as_date("2020/11/09")

# Only used for the Platform/Subcategory columns
cgd_master_input <- read_xlsx(file.path("data", "Master Input Data - COVID-19 Vaccine Predictions publication-2.xlsx")) %>% 
  transmute(
    cgd_vaccine_id = Number, #name, institutes, country, funding, firm_size,
    Platform = str_replace(Platform, ".+viral vector", "Viral vector"), # Not distinguishing between replicating and non-replicating
    Subcategory,
  )

cgd_params_json_data <- jsonlite::fromJSON(file.path("data", "cgd_params.json")) 

cgd_json_data <- jsonlite::fromJSON(file.path("data", "cgd_vaccines.json")) %>% 
  mutate(
    cgd_vaccine_id = as.integer(number),
  ) %>% 
  select(cgd_vaccine_id, name, institutes, country, platform_key, funding = funding_category, funding_key, firm_size, 
         dev_phase = phase, start_date = start_dates, end_date = end_dates) %>% 
  mutate(cgd_platform = cgd_params_json_data$platforms[platform_key + 1]) %>% 
  rowwise() %>% 
  mutate(
    phase_dates = list(
      tibble(phase = 0:4, start_date, end_date) %>% 
        mutate(across(c(start_date, end_date), ~ lubridate::as_datetime(na_if(.x, 0)))) %>% 
        filter(!is.na(start_date))
    ),
  ) %>% 
  ungroup() %>% 
  mutate(phase = map_chr(
    phase_dates, ~ {
      phase <- .x %>% 
        filter(start_date <= cgd_draws_date) %>% 
        pull(phase) %>% 
        max()
      
      if (max(.x$end_date, na.rm = TRUE) < cgd_draws_date) {
        phase %<>% add(1L) 
      }
      
      if (phase > 0) { 
        str_c("Phase ", phase) 
      } else {
        "Pre-clinical"
      }
    })
  ) %>% 
  select(!c(start_date, end_date)) %>% 
  left_join(select(cgd_master_input, cgd_vaccine_id, Platform, Subcategory), by = "cgd_vaccine_id") %>% 
  mutate( # These have to be manually set after matching them with the AHT vaccines
    Subcategory = if_else(cgd_vaccine_id == 18, "Inactivated", Subcategory),
    Subcategory = if_else(cgd_vaccine_id == 115, "S Protein", Subcategory),
    Subcategory = if_else(cgd_vaccine_id == 139, "Measles (replicating)", Subcategory),
  )

cgd2aht_summary <- cgd_json_data %>% 
  filter(fct_match(phase, c("Phase 2", "Phase 3"))) %>% 
  count(Platform, Subcategory, phase) %>% 
  mutate(phase = str_remove(phase, " ") %>% str_c("candidates")) %>% 
  full_join(long_aht_summary, by = c("Platform", "Subcategory", "phase"), suffix = c("_cgd", "_aht")) %>% 
  mutate(n = coalesce(n_cgd, n_aht)) %>%  
  select(!c(n_aht, n_cgd)) %>% 
  pivot_wider(Platform:PreviousState, names_from = phase, values_from = n)

aht_vaccines_file <- file.path("data", str_c(script_options$output, "_aht_vaccines.csv"))
write_csv(cgd2aht_summary, aht_vaccines_file)

raw_cgd_trials <- read_csv(
  file.path("data", "cgd_trials.csv"), 
  skip = 1,
  col_names = c("try_id", "vaccine_id", "phase_mon_1", "phase_mon_2", "phase_mon_3", "phase_mon_approval"), 
) %>% 
  mutate_all(as.integer) 

# Settings ----------------------------------------------------------------

group_vaccines_by <- NULL 
maxcand <- 75 
candidate_data <- #withr::with_tempfile(new = "aht_summary", tmpdir = "temp", {
  # load_candidate_data(file.path("data", "vaccinesSummaryCGD.csv"))
  # write_csv(cgd2aht_summary, aht_summary)
  load_candidate_data(aht_vaccines_file)
# })

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
      cgd_json_data %>% 
        filter(
          fct_match(phase, c("Phase 2", "Phase 3")),
          across(c(Platform, Subcategory), ~ ! is.na(.x)) # BUG Need to figure out these columns in the CGD data
        ) %>% 
        nest(cgd_ids = c(cgd_vaccine_id, name, institutes, country, platform_key, cgd_platform, funding, funding_key, firm_size, dev_phase, phase_dates)),
      by = c("Platform", "Subcategory", "phase")
    ) %>% 
    mutate(
      ids = map2(ids, cgd_ids, ~ bind_cols(.x, sample_frac(.y))) # Randomly match IDs
    ) %>% 
    select(-cgd_ids) %>% 
    unnest(ids) 
} else {
  read_rds(file.path("data", "cgd_id_dict.rds"))
}

cgd_trials <- raw_cgd_trials %>% 
  rename(r = try_id) %>% 
  pivot_longer(phase_mon_1:phase_mon_approval, names_to = "phase", values_to = "approval_month", names_prefix = "phase_mon_") %>% 
  filter(!is.na(approval_month)) %>%
  group_by(r, vaccine_id) %>% 
  filter(approval_month == max(approval_month)) %>% 
  filter(phase == max(phase)) %>% # Sometimes two phases would have the same max(approval_month) 
  ungroup() %>% 
  complete(r, vaccine_id = cgd_id_dict$cgd_vaccine_id) %>% 
  mutate(
    phase = coalesce(phase, "none"),
  )

if (script_options$pfizer_pref) {
  num_cgd_trials <- max(cgd_trials$r)
  
  cgd_trials %<>% 
    nest(vacc_phase = -r) %>%
    mutate(
      pfizer_before_mar = future_map_lgl(vacc_phase, 
                                         ~ filter(.x, vaccine_id == 150, phase == "approval", approval_month <= 3) %>% 
                                           nrow() %>% 
                                           is_greater_than(0))
    ) %>% 
    group_by(pfizer_before_mar) %>% 
    sample_n(num_cgd_trials / 2, replace = TRUE) %>% 
    ungroup() %>% 
    mutate(r = seq(n())) %>% 
    unnest(vacc_phase) %>% 
    select(-pfizer_before_mar)
}

cgd_draws_month <- cgd_draws_date %>% subtract(lubridate::day(.) - 1) 
fit_month_offsets <- seq(script_options$num_months) + 1

initial_par <- c(poverall = 0.5, psubcat = 0.5, 
                 pvector = 0.5, psubunit = 0.5, prna = 0.5, pdna = 0.5, pinactivated = 0.5, pattenuated = 0.5,
                 pphase2 = 0.5, pphase3 = 0.5)

# CGD Fit -----------------------------------------------------------------

cgd_data <- tibble(
  month_offset = fit_month_offsets,
  month = cgd_draws_month + months(fit_month_offsets)
) %>% 
  mutate(
    cgd_converted = future_map(month_offset, convert_cgd_trials, cgd_trials = cgd_trials, id_dict = cgd_id_dict, .progress = TRUE)
  ) 

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
          maxcand = maxcand,
          replications = 100e3,
          ndeps = rep_along(initial_par, 1e-2)
        )
      },
      .progress = TRUE
    )
  ) 
  
  write_rds(cgd_optim_data, file.path("data", str_c(script_options$output, ".rds")))
  
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
      .progress = TRUE
    )
  ) %>% 
  rowwise() %>% 
  mutate(run_summary = list(last(param_data$summary))) %>% 
  ungroup()

write_rds(cgd_optim_data, file.path("data", str_c(script_options$output, ".rds")))

cat("...done.\n")