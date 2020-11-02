#!/usr/bin/env Rscript

"Usage:
  cgd_fit.R [--num-runs=<runs> --append-stage-2=<run-data-file> --output=<run-data-file> --no-vcov-moments]
  
Options:
  --num-runs=<runs>  Number of runs per month [default: 1]
  --output=<run-data-file>  Where to store output [default: cgd_optim.rds]
" -> opt_desc

script_options <- if (interactive()) {
  docopt::docopt(opt_desc, '')
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

plan(multiprocess)

source("model_fit_util.R")

script_options %<>% 
  modify_at(c("num_runs"), as.integer)

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
  mutate(
    ids = map2(ids, cgd_ids, ~ bind_cols(.x, sample_frac(.y))) # Randomly match IDs
  ) %>% 
  unnest(ids)

cgd_draws_month <- lubridate::as_date("2020/10/01")
fit_month_offsets <- seq(12) + 8

initial_par <- c(poverall = 0.5, psubcat = 0.5, 
                 pvector = 0.5, psubunit = 0.5, prna = 0.5, pdna = 0.5, pinactivated = 0.5, pattenuated = 0.5,
                 pphase2 = 0.5, pphase3 = 0.5)

# CGD Fit -----------------------------------------------------------------

cgd_optim_data <- tibble(
  month_offset = fit_month_offsets,
  month = cgd_draws_month + months(fit_month_offsets)
) %>% 
  rowwise() %>% 
  mutate(
    cgd_summary = list(convert_cgd_trials(month_offset, cgd_trials, cgd_id_dict)) %>% 
      future_map(summarize_draws)
  ) %>% 
  ungroup()

cgd_optim_data %<>% 
  modelr::data_grid(run_id = seq(script_options$num_runs), month_offset) %>% 
  left_join(cgd_optim_data, "month_offset") 

if (is_null(script_options$append_stage_2)) {
  cat("First stage...")
  
  cgd_optim_data %<>% 
    mutate(param_data = future_map(
      cgd_summary,
      ~ optim_run(
        NULL, 
        .x, 
        initial_par = initial_par,
        use_vcov_moments = !script_options$no_vcov_moments,
        maxcand = maxcand),
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
      ~ optim_run(
        NULL,
        .x,
        prev_run_data = .y,
        initial_par = initial_par,
        maxcand = maxcand,
        replications = 100e3, 
        ndeps = rep_along(initial_par, 2e-3)
      ),
      .progress = TRUE),
    run_summary = map(param_data, filter, step == n()) %>%  
      future_map(~ {
        get_candidate_draws(
          candidate_data = candidate_data, replications = 3e5, dordered = dordered,
          param = .x,
          maxcand = maxcand,
          group_vaccines_by = group_vaccines_by
        ) %>%
          get_summary_success_rates()
      },
      .progress = TRUE)
  )

cat("...done.\n")

write_rds(cgd_optim_data, file.path("data", script_options$output))

if (!interactive()) {
  quit(save = "no")
}

# Plot --------------------------------------------------------------------

cgd_optim_data %>% { 
    cowplot::plot_grid(
    #   pivot_longer(., -c(step, run_id), names_to = "param_name", values_to = "param_val") %>% 
    #     ggplot() +
    #     geom_line(aes(step, param_val, color = factor(run_id), group = factor(run_id)), alpha = 0.75, show.legend = TRUE) +
    #     facet_wrap(vars(param_name), scales = "free") +
    #     labs(y = "") +
    #     theme_minimal(),
     
      select(., run_id, month, param_data) %>% 
        unnest(param_data) %>% 
        group_by(run_id, month) %>%
        filter(step == n()) %>%
        ungroup() %>%
        select(-obj) %>%  
        pivot_longer(., -c(step, run_id, month), names_to = "param_name", values_to = "param_val") %>% 
        ggplot() +
        geom_line(aes(month, param_val, group = run_id), alpha = 0.5) +
        scale_x_date("", breaks = "2 months", labels = scales::date_format("%m/%y")) +
        labs(title = "Parameter fit solutions", y = "") +
        facet_wrap(vars(param_name)) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        NULL,
    
        ggplot(.) +
        geom_line(aes(month, success_rate, group = run_id), alpha = 0.5,
                  data = . %>%
                    mutate(run_summary = map2(run_summary, cgd_summary, ~ semi_join(.x, .y$success_rates, by = "vacc_group_id"))) %>% 
                    select(run_id, month, run_summary) %>% 
                    unnest(run_summary)) +
        geom_line(aes(month, success_rate), color = "red", 
                  data = . %>% 
                    filter(run_id == 1) %>% 
                    select(month, cgd_summary) %>% 
                    mutate(cgd_summary = map(cgd_summary, pluck, "success_rates")) %>% 
                    unnest(cgd_summary)) + 
                    # left_join(cgd_id_dict, by = c("vacc_group_id" = "candInd"))) + 
        scale_x_date("", breaks = "2 months", labels = scales::date_format("%m/%y")) +
        labs(title = "Solution moments",
             subtitle = "Restricted to phase 2 and 3 vaccines.",
             caption = "Using 10 runs.\nThe red line is the moment calculated from the CGD model data.", y = "") +
        facet_wrap(vars(vacc_group_id), 
                   labeller = labeller(
                     vacc_group_id = function(ids) { 
                       cgd_id_dict %>% 
                         right_join(tibble(candInd = as.integer(ids)), by = c("candInd")) %$% 
                         str_glue("CGD ID: {cgd_vaccine_id}\n{Platform}\n{Subcategory}\n{phase}") %>% 
                         as.character()
                     })) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        NULL
    )
  }

