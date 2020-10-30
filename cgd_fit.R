#!/usr/bin/env Rscript

library(data.table)
library(magrittr)
library(tidyverse)
library(readxl)
library(gmm)
library(furrr)
library(R6)

library(vaccineEarlyInvest)

plan(multiprocess)

source("model_fit_util.R")

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
    # Subcategory = if_else(str_detect(Institutes, "Hong Kong"), "Influenza (replicating)", Subcategory),
    phase = map_chr(phase_data, ~ if (nrow(.x) > 0) str_c("Phase ", max(.x$phase)) else "Pre-clinical")
  )

cgd_trials <- read_csv(
  file.path("data", "cgd_trials.csv"), 
  skip = 1,
  col_names = c("try_id", "vaccine_id", "phase_mon_1", "phase_mon_2", "phase_mon_3", "phase_mon_approval"), 
) %>% 
  mutate_all(as.integer) 

  # filter(phase_mon_approval <= start_month_offset)

# Settings ----------------------------------------------------------------

group_vaccines_by <- NULL 
maxcand <- 50 
candidate_data <- load_candidate_data(file.path("data", "vaccinesSummaryCGD.csv"))
dordered <- calc_dordered(candidate_data, maxcand = maxcand)

cgd_id_dict <- get_candidate_draws(
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
  transmute(
    ids = map2(ids, cgd_ids, ~ bind_cols(.x, sample_frac(.y))) # Randomly match IDs
  ) %>% 
  unnest(ids)

# CGD Fit -----------------------------------------------------------------

cgd_summaries <- map(seq(12) + 8, convert_cgd_trials, cgd_trials, cgd_id_dict) %>% 
  map(summarize_draws)

cgd_optim_data <- future_map2(
  seq(cgd_summaries),
  cgd_summaries,
  ~ optim_run(
    .x, 
    initial_par = c(poverall = 0.5, psubcat = 0.5, 
                    pvector = 0.5, psubunit = 0.5, prna = 0.5, pdna = 0.5, pinactivated = 0.5, pattenuated = 0.5,
                    pphase1 = 0.5, pphase2 = 0.5, pphase3 = 0.5),
    summaries = .y,
    maxcand = maxcand
  ))

write_rds(cgd_optim_data, file.path("data", "cgd_optim.rds"))
