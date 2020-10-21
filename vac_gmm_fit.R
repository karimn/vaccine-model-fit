library(data.table)
library(magrittr)
library(tidyverse)
library(gmm)
library(furrr)
library(R6)

library(vaccineEarlyInvest)

# plan(multisession)

# Load the CGD data and convert it ----------------------------------------

zero_month_date <- lubridate::as_date("2020-10-01") 
start_month_offset <- 9

cgd_trials <- read_csv(
  file.path("data", "cgd_trials.csv"), 
  skip = 1,
  col_names = c("try_id", "vaccine_id", "phase_mon_1", "phase_mon_2", "phase_mon_3", "phase_mon_approval"), 
) %>% 
  mutate_all(as.integer) %>% 
  filter(phase_mon_approval <= start_month_offset)


# Moments -----------------------------------------------------------------

load_candidate_data <- function(data_file) {
  candidate_data <- loadData(par = NULL, data_file)
  
  candidate_data$Target <- "Other"
  candidate_data$Target[1:5]<-"Spike"
  candidate_data$Target[10:15]<-"Recombinant"
  
  return(candidate_data)
}

get_candidate_draws <- function(candidate_data, replications, 
                                ...,
                                maxcand,
                                group_vaccines_by = vars(Platform, Subcategory),
                                seed = NULL) {
  param <- rlang::list2(...)
  par0 <- Parameters$new(maxcand = maxcand)  
  par <- exec(Parameters$new, replications = replications, !!!param, maxcand = maxcand)  
  
  draws <- candidate_data %>% 
    candidatesFung(par0) %>% 
    pluck("dordered") %>% 
    select(1:11) %>% 
    candidateDraws(par, seed = seed) 
    
  if (is_null(group_vaccines_by)) {
    draws %<>% 
      mutate(
        vacc_group_id = candInd,
        any_success = success, 
        success_rate = success,
      )
  } else {
    vacc_group_ids <- candidate_data %>%
      pivot_longer(Phase1candidates:PreClinicalCandidates, names_to = "phase") %>% # Add phase to ID columns
      filter(value > 0) %>% 
      mutate(phase = fct_relabel(phase, str_replace_all, c("Phase(\\d)candidates" = "Phase \\1", "PreClinicalCandidates" = "Pre-clinical"))) %>% 
      group_by_at(group_vaccines_by) %>% 
      summarize(vacc_group_id = cur_group_id(), .groups = "drop")

    draws %<>% 
      left_join(vacc_group_ids, by = setdiff(names(vacc_group_ids), "vacc_group_id")) %>% # Get the vaccine group IDs
      group_by(r, vacc_group_id, .drop = FALSE) %>% 
      summarize(any_success = sum(success) > 0, 
                success_rate = mean(success),
                .groups = "drop") %>%  # Any success in group
      complete(r, vacc_group_id = seq(max(vacc_group_ids$vacc_group_id)), 
               fill = lst(any_success = FALSE, success_rate = 0)) # Make sure all groups are present in the draws
  } 
  
  return(draws)
}

summarize_draws <- function(draws) {
    success_rates <- draws %>% 
      group_by(vacc_group_id) %>% 
      summarize(
        success_rate = mean(success_rate),
        
        .groups = "drop"
      ) %>% 
      ungroup()
    
    success_corr <- draws %>% 
      pivot_wider(id_cols = r, names_from = vacc_group_id, values_from = success_rate) %>% 
      select(-r) %>% 
      as.data.frame() %>% 
      pcaPP::cor.fk() %>% # Faster than base cor()
      magrittr::extract(lower.tri(.)) # Get lower triangle
    
    success_vcov <- draws %>% 
      pivot_wider(id_cols = r, names_from = vacc_group_id, values_from = success_rate) %>% 
      select(-r) %>% 
      as.data.frame() %>% 
      cov() %>% 
      magrittr::extract(lower.tri(.)) # Get lower triangle
    
    lst(success_rates, success_corr, success_vcov)
}

calculate_draws_vcov <- function(draws) {
  draws %>% 
    pivot_wider(id_cols = r, names_from = vacc_group_id, values_from = success_rate) %>% 
    select(-r) %>% 
    cov() 
}

OptimLogger <- R6Class("OptimLogger", 
                       list(
                         data = tibble(),
                         log = function(...) { self$data %<>% bind_rows(tibble(!!!rlang::list2(...))) }))

build_gmm_g <- function(candidate_data, x, replications, maxcand, use_vcov_moments = FALSE, group_vaccines_by = vars(Platform, Subcategory), sim_seed = NULL, moments_to_use = everything(), fixed_model_probs = NULL, logger = NULL, calculate_objective = FALSE) {
  function(pbeta, ...) {
    model_probs <- plogis(pbeta)
    all_model_probs <- list_modify(as.list(fixed_model_probs), !!!model_probs)
    
    model_probs %>% 
      imap(~ str_c(.y, " = ", .x)) %>% 
      str_c(collapse = ", ") %>% 
      cat("\n")
    
    draws <- get_candidate_draws(candidate_data, replications = replications, !!!all_model_probs, maxcand = maxcand, group_vaccines_by = group_vaccines_by, seed = sim_seed)  

    model_summaries <- summarize_draws(draws) 
    
    moments <- full_join(model_summaries$success_rates, x$success_rates, by = "vacc_group_id") %>%
      mutate_at(vars(starts_with("success_rate")), coalesce, 0) %>%
      transmute(vacc_group_id,
                success_moment = success_rate.x - success_rate.y) %>%
      pull(success_moment)
    
    if (use_vcov_moments) {
      moments %<>% 
        c(model_summaries$success_vcov - x$success_vcov) 
    }
    
    moments %<>% t() 

   # moments <- full_join(x, model_summaries$success_rates, by = "vacc_group_id") %>%
   #    mutate_at(vars(starts_with("success_rate")), coalesce, 0) %>%
   #    transmute(r, vacc_group_id, m = success_rate.x - success_rate.y) %>%
   #    pivot_wider(names_from = vacc_group_id, values_from = m) %>%
   #    select(-r) %>%
   #    mutate_all(coalesce, 0) %>% 
   #    select(moments_to_use)
   
   
   objective <- moments %>% colMeans() %>% crossprod() %>% c()
    
    if (!is_null(logger)) {
      logger$log(!!!model_probs, obj = objective) 
    }
   
    if(calculate_objective) {
      return(objective)
    } else {
      return(moments)
    }
  }
}

# Settings ----------------------------------------------------------------

group_vaccines_by <- NULL #vars(Platform, Subcategory, Target, phase)
candidate_data <- load_candidate_data(file.path("data", "vaccinesSummaryOct2.csv"))

true_param <- c(poverall=0.9, psubcat=0.9, pvector=0.8, psubunit=0.8, prna=0.6, pdna=0.4, pattenuated=0.8, pinactivated=0.8, ppreclinical=0.14, pphase1=0.23, pphase2=0.32, pphase3=0.5)

# Test Data --------------------------------------------------------------------

test_draws <- get_candidate_draws(
  candidate_data, replications = 3e5,
  !!!true_param,
  maxcand = 50,
  group_vaccines_by = group_vaccines_by
)

test_summaries <- summarize_draws(test_draws)

# test_vcov <- calculate_draws_vcov(test_draws)
# 
# moments_to_use <- test_vcov %>% 
#   diag() %>% 
#   equals(0) %>% 
#   not() %>% 
#   which()

test_log <- OptimLogger$new()

gmm_g <- build_gmm_g(candidate_data, 10e3, 50, use_vcov_moments = FALSE, group_vaccines_by = group_vaccines_by, sim_seed = 123,
                     fixed_model_probs = true_param, 
                     logger = test_log) 

# Test optimization -------------------------------------------------------

# test_optim_log <- OptimLogger$new()

test_optim <- optim(
  fn = build_gmm_g(
    candidate_data, test_summaries, 20e3, 50, use_vcov_moments = TRUE, group_vaccines_by = group_vaccines_by, # sim_seed = 123,
    fixed_model_probs = true_param,
    calculate_objective = TRUE,
    logger = test_optim_log
  ),
  par = c(poverall = 0.5, psubcat = 0.5, 
          pvector = 0.5, psubunit = 0.5, prna = 0.5, pdna = 0.5, pinactivated = 0.5,
          pphase1 = 0.5, pphase2 = 0.5, pphase3 = 0.5) %>% qlogis(),
  control = lst(reltol = 0.00001)
  # method = "BFGS"
)

# Plot test optimization dynamics -----------------------------------------

test_optim_log$data %>%
  group_by(run_id) %>% 
  mutate(step = seq(n())) %>% 
  ungroup() %>% 
  pivot_longer(-c(step, run_id), names_to = "param_name", values_to = "param_val") %>% 
  ggplot() +
  geom_line(aes(step, param_val, color = factor(run_id)), alpha = 0.75) +
  geom_hline(aes(yintercept = value), 
             linetype = "dotted",
             data = . %>% semi_join(enframe(true_param, name = "param_name"), ., by = "param_name")) +
  facet_wrap(vars(param_name), scales = "free") +
  labs(y = "") +
  theme_minimal()

# Test GMM run ------------------------------------------------------------

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


# Test candInd ------------------------------------------------------------

test_draws1 <- get_candidate_draws(
  candidate_data, replications = 1,
  !!!true_param,
  maxcand = 50,
  group_vaccines_by = NULL 
)

test_draws2 <- get_candidate_draws(
  candidate_data, replications = 1,
  poverall=0.9, psubcat=0.5, pvector=0.8, psubunit=0.8, prna=0.3, pdna=0.8, pattenuated=0.8, pinactivated=0.8, ppreclinical=0.14, pphase1=0.6, pphase2=0.3, pphase3=0.2,
  maxcand = 50,
  group_vaccines_by = NULL 
)

full_join(
  test_draws1 %>% 
    select(candInd, Platform, Subcategory, Target, phase),
  test_draws2 %>% 
    select(candInd, Platform, Subcategory, Target, phase),
  by = "candInd")
