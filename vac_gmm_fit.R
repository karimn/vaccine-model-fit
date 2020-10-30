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

raw_cgd_trials <- read_csv(
  file.path("data", "cgd_trials.csv"), 
  skip = 1,
  col_names = c("try_id", "vaccine_id", "phase_mon_1", "phase_mon_2", "phase_mon_3", "phase_mon_approval"), 
) 

cgd_trials <- raw_cgd_trials %>%  
  mutate_all(as.integer) %>% 
  filter(phase_mon_approval <= start_month_offset)

# Settings ----------------------------------------------------------------

group_vaccines_by <- NULL 
maxcand <- 50 
candidate_data <- load_candidate_data(file.path("data", "vaccinesSummaryCGD.csv"))
dordered <- calc_dordered(candidate_data, maxcand = maxcand)

true_param <- c(poverall=0.9, psubcat=0.9, pvector=0.8, psubunit=0.8, prna=0.6, pdna=0.4, pattenuated=0.8, pinactivated=0.8, ppreclinical=0.14, pphase1=0.23, pphase2=0.32, pphase3=0.5)

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

cgd_summaries <- convert_cgd_trials(cgd_trials, cgd_id_dict) %>% 
  summarize_draws()

cgd_optim_data <- future_map(1:4, optim_run,
                      summaries = cgd_summaries,
                      maxcand = maxcand)

# Test Data --------------------------------------------------------------------

test_draws <- get_candidate_draws(
    candidate_data, replications = 3e5, dordered = dordered,
    param = true_param,
    maxcand = maxcand,
    group_vaccines_by = group_vaccines_by
  )

test_summaries <- summarize_draws(test_draws)

# Test optimization -------------------------------------------------------

test_optim_data <- map(1:4, 
# future_walk(test_optim_logs,
  optim_run, summaries = test_summaries, maxcand = maxcand
) #, .options = furrr_options(seed = TRUE))

test_optim_data %<>%
  group_nest(run_id, .key = "run_data") %>% 
  pull(run_data) %>% 
  map2_dfr(1:4, ., optim_run, replications = 100e3, ndeps = rep(2e-3, 10))

# Plot test optimization dynamics -----------------------------------------

cowplot::plot_grid(
  test_optim_data %>%
    pivot_longer(-c(step, run_id), names_to = "param_name", values_to = "param_val") %>% 
    ggplot() +
    geom_line(aes(step, param_val, color = factor(run_id)), alpha = 0.75, show.legend = FALSE) +
    geom_hline(aes(yintercept = value), 
               linetype = "dotted",
               data = . %>% semi_join(enframe(true_param, name = "param_name"), ., by = "param_name")) +
    facet_wrap(vars(param_name), scales = "free") +
    labs(y = "") +
    theme_minimal(),

  test_optim_data %>%
    group_by(run_id) %>% 
    filter(step == n()) %>% 
    ungroup() %>% 
    mutate(!!!true_param[setdiff(names(true_param), names(.))]) %>% 
    select(run_id, all_of(names(true_param))) %>% 
    group_by(run_id) %>% 
    group_modify(quick_get_summary, summary_type = "success_rates") %>% 
    ungroup() %>% 
    ggplot() +
    geom_col(aes(run_id, success_rate, fill = factor(run_id)), alpha = 0.75, show.legend = FALSE) +
    geom_hline(aes(yintercept = success_rate), linetype = "dashed", data = test_summaries$success_rates) +
    facet_wrap(vars(vacc_group_id)) +
    theme_minimal() +
    NULL
)

# Linear combination of optim solution and true param ---------------------

model_combine <- test_optim_data %>%
  group_by(run_id) %>% 
  filter(step == n()) %>% 
  ungroup() %>% 
  mutate(!!!true_param[setdiff(names(true_param), names(.))]) %>% 
  select(run_id, all_of(names(true_param))) %>% 
  group_by(run_id) %>% 
  group_modify(function(fit_param, true_param, ...) {
    tibble(
      alpha = seq(-0.5, 1.5, 0.1),
      param_mix = map(alpha, ~ (1 - .x) * fit_param + .x * true_param)
    )
  }, true_param = true_param) %>% 
  ungroup() %>% 
  mutate(
    draw_summaries = map(param_mix, quick_get_summary, summary_type = NULL),
    moments = map(draw_summaries, calculate_moments, test_summaries),
    objective = map_dbl(moments, calculate_objective),
  ) 

cowplot::plot_grid(
  model_combine %>% 
    select(run_id, alpha, draw_summaries) %>% 
    mutate(draw_summaries = map(draw_summaries, pluck, "success_rates")) %>% 
    unnest(draw_summaries) %>% 
    ggplot() +
    geom_line(aes(alpha, success_rate, color = factor(run_id)), show.legend = FALSE) +
    geom_vline(xintercept = c(0, 1), linetype = "dotted", color = "darkgrey") +
    facet_wrap(vars(vacc_group_id)) +
    labs(title = "Success Rates", y = "") +
    theme_minimal() +
    theme(axis.text.x = element_blank()),
  
  model_combine %>% 
    select(run_id, alpha, objective) %>% 
    ggplot() +
    geom_line(aes(alpha, objective, color = factor(run_id)), show.legend = FALSE) +
    geom_vline(xintercept = c(0, 1), linetype = "dotted", color = "darkgrey") +
    labs(title = "Objective", y = "") +
    theme_minimal()
)

# Linear combination of optim solution and wrong param --------------------

rand_param <- test_optim_data %>% 
  select(any_of(names(true_param))) %>% 
  slice(1:5) %>% 
  mutate_all(~ runif(1)) %>% 
  mutate(!!!true_param[setdiff(names(true_param), names(.))]) 

rand_model_combine <- test_optim_data %>%
  filter(run_id == 1) %>% 
  group_by(run_id) %>% 
  filter(step == n()) %>% 
  ungroup() %>% 
  mutate(!!!true_param[setdiff(names(true_param), names(.))]) %>% 
  select(run_id, all_of(names(true_param))) %>% 
  group_by(run_id) %>% 
  group_modify(~ bind_cols(.x, tibble(alpha = seq(-0.5, 1.5, 0.1)))) %>%
  ungroup() %>% 
  group_nest(run_id, alpha, .key = "param_mix") %>% 
  mutate(
    param_mix = map2(param_mix, alpha, function(param_mix, alpha, rand_param) {
      as_tibble(param_mix[rep(1, nrow(rand_param)), names(rand_param)] * (1 - alpha) + rand_param * alpha)
    }, rand_param = rand_param),
    
    draw_summaries = map(param_mix, 
                         ~ mutate(.x, rand_id = seq(n())) %>%
                           group_by(rand_id) %>%
                           group_modify(~ tibble(draw_summaries = list(quick_get_summary(.x, summary_type = NULL))))),
  ) %>% 
  unnest(c(param_mix, draw_summaries)) %>% 
  mutate(
    moments = map(draw_summaries, calculate_moments, test_summaries),
    objective = map_dbl(moments, calculate_objective),
  ) 

rand_model_combine %>% 
  select(run_id, alpha, objective) %>% 
  ggplot() +
  geom_line(aes(alpha, objective, color = factor(run_id)), show.legend = FALSE) +
  geom_vline(xintercept = c(0, 1), linetype = "dotted", color = "darkgrey") +
  labs(title = "Objective", y = "") +
  theme_minimal()

# Test GMM run ------------------------------------------------------------

test_log <- OptimLogger$new()

gmm_g <- build_gmm_g(candidate_data, 10e3, maxcand, use_vcov_moments = FALSE, group_vaccines_by = group_vaccines_by, sim_seed = 123,
                     fixed_model_probs = true_param, 
                     logger = test_log) 

test_results <- gmm(gmm_g, 
                    test_draws,  
                    wmatrix = "ident",
                    t0 = c(poverall = 0.5, psubcat = 0.5) %>% qlogis(), 
                    # method = "Brent", lower = -10, upper = 10, 
                    control = lst(reltol = 0.0001))
                    # t0 = c(poverall=0.9, psubcat=0.9) %>% qlogis()) 
#, pvector=0.8, psubunit=0.8, prna=0.6, pdna=0.4, pattenuated=0.8, pinactivated=0.8, ppreclinical=0.14, pphase1=0.23, pphase2=0.32, pphase3=0.5) %>% qlogis()) 
                    # weightsMatrix = solve(test_vcov[moments_to_use, moments_to_use]))
                    #t0 = rep(0, 12)) # rnorm(12)) 
                    # type = "iterative") 

# Check continuity of moments ---------------------------------------------

moment_landscape <- expand.grid(poverall_beta = seq(-5, 5, 1), psubcat_beta = seq(-5, 5, 1)) %>% 
  mutate(moments = future_pmap(lst(poverall_beta, psubcat_beta), ~ gmm_g(c(.x, .y, rep(0, 10)), test_summaries))) 

moment_landscape %>% 
  mutate(moments = map(moments, as_tibble) %>% 
           map(~ mutate(.x, moment_id = seq(n())))) %>% 
  unnest(moments) %>% 
  ggplot() +
  geom_raster(aes(poverall_beta, psubcat_beta, fill = value)) +
  facet_wrap(vars(moment_id))

test_log$data %>% 
  ggplot() +
  geom_line(aes(poverall, obj))
  # geom_tile(aes(poverall, psubcat, fill = obj))
  # geom_point(aes(poverall, psubcat, color = obj))


# Test candInd uniqueness ------------------------------------------------------------

test_draws1 <- get_candidate_draws(
  candidate_data, replications = 1,
  param = true_param,
  maxcand = maxcand,
  group_vaccines_by = NULL 
)

test_draws2 <- get_candidate_draws(
  candidate_data, replications = 1,
  param = c(poverall=0.9, psubcat=0.5, pvector=0.8, psubunit=0.8, prna=0.3, pdna=0.8, pattenuated=0.8, pinactivated=0.8, ppreclinical=0.14, pphase1=0.6, pphase2=0.3, pphase3=0.2),
  maxcand = maxcand,
  group_vaccines_by = NULL 
)

full_join(
  test_draws1 %>% 
    select(candInd, Platform, Subcategory, Target, phase),
  test_draws2 %>% 
    select(candInd, Platform, Subcategory, Target, phase),
  by = "candInd")
