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

calc_dordered <- function(candidate_data, maxcand) {
  par0 <- Parameters$new(maxcand)
  
  candidatesFung(candidate_data, par0)$dordered
}

get_candidate_draws <- function(candidate_data, replications, dordered, 
                                ...,
                                maxcand,
                                group_vaccines_by = vars(Platform, Subcategory),
                                seed = NULL) {
  param <- rlang::list2(...)
  par <- exec(Parameters$new, replications = replications, !!!param, maxcand = maxcand)  
  
  dordered[Platform == "DNA", pplat := as.numeric(par$pdna)]
  dordered[Platform == "RNA", pplat := par$prna]
  dordered[Platform == "Live attenuated virus", pplat := par$pattenuated]
  dordered[Platform == "Viral vector", pplat := par$pvector]
  dordered[Platform == "Protein subunit", pplat := par$psubunit]
  dordered[Platform == "Inactivated", pplat := par$pinactivated]
  dordered[Platform == "VLP", pplat := par$pvlp]
  dordered[Platform == "Dendritic cells", pplat := par$pdendritic]
  dordered[Platform == "Self-assembling vaccine", pplat := par$psav]
  dordered[Platform == "Unknown", pplat := par$punknown]
  dordered[Platform == "Artificial antigen presenting cells", pplat := par$paapc]
  dordered[Platform == "Live-attenuated bacteria", pplat := par$plivebac]
  
  dordered[phase == "Pre-clinical", pcand := par$ppreclinical]
  dordered[phase == "Phase 1", pcand := par$pphase1]
  dordered[phase == "Phase 2", pcand := par$pphase2]
  dordered[phase == "Phase 3", pcand := par$pphase3]
  dordered[phase == "Repurposed", pcand := par$prepurposed]
  
  dordered <- dordered[,1:11]
  dcandidate <- copy(dordered)

  draws <- candidateDraws(dcandidate, par, seed = seed)
    
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

calculate_moments <- function(model_summaries, x, use_vcov_moments = FALSE) {
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
}

calculate_objective <- function(moments, weighting_matrix = diag(length(moments))) {
  colMeans(moments) %>% 
     tcrossprod(. %*% weighting_matrix, .) %>% 
     c()
}

build_gmm_g <- function(candidate_data, dordered, x, replications, maxcand, 
                        use_vcov_moments = FALSE, 
                        group_vaccines_by = vars(Platform, Subcategory), 
                        sim_seed = NULL, 
                        moments_to_use = everything(), 
                        fixed_model_probs = NULL, 
                        logger = NULL, 
                        calculate_objective = FALSE, 
                        weighting_matrix = diag(nrow(dordered) + use_vcov_moments * sum(seq(nrow(ordered) - 1))), 
                        verbose = FALSE) {
  function(pbeta, ...) {
    model_probs <- plogis(pbeta)
    all_model_probs <- list_modify(as.list(fixed_model_probs), !!!model_probs)
    
    if (verbose) {
      model_probs %>% 
        imap(~ str_c(.y, " = ", .x)) %>% 
        str_c(collapse = ", ") %>% 
        cat("\n")
    }
    
    draws <- get_candidate_draws(candidate_data, replications = replications, dordered = dordered, !!!all_model_probs, maxcand = maxcand, group_vaccines_by = group_vaccines_by, seed = sim_seed)  

    model_summaries <- summarize_draws(draws) 
    
    moments <- calculate_moments(model_summaries, x) 

   # moments <- full_join(x, model_summaries$success_rates, by = "vacc_group_id") %>%
   #    mutate_at(vars(starts_with("success_rate")), coalesce, 0) %>%
   #    transmute(r, vacc_group_id, m = success_rate.x - success_rate.y) %>%
   #    pivot_wider(names_from = vacc_group_id, values_from = m) %>%
   #    select(-r) %>%
   #    mutate_all(coalesce, 0) %>% 
   #    select(moments_to_use)
   
   objective <- calculate_objective(moments, weighting_matrix) 
    
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
dordered <- calc_dordered(candidate_data, 50)

true_param <- c(poverall=0.9, psubcat=0.9, pvector=0.8, psubunit=0.8, prna=0.6, pdna=0.4, pattenuated=0.8, pinactivated=0.8, ppreclinical=0.14, pphase1=0.23, pphase2=0.32, pphase3=0.5)

quick_get_summary <- function(param, ..., summary_type = "success_rates") {
  summ <- get_candidate_draws(
    candidate_data, replications = 3e5, dordered = dordered,
    !!!param,
    maxcand = 50,
    group_vaccines_by = group_vaccines_by
  ) %>% 
    summarize_draws()
  
  if (!is_null(summary_type)) {
    return(pluck(summ, summary_type))
  } else {
    return(summ)
  }
}

# Test Data --------------------------------------------------------------------

test_draws <- get_candidate_draws(
    candidate_data, replications = 3e5, dordered = dordered,
    !!!param,
    maxcand = 50,
    group_vaccines_by = group_vaccines_by
  )

test_summaries <- summarize_draws(test_draws)

# Test optimization -------------------------------------------------------

test_optim_data <- map_dfr(1:4, ~ {
# future_walk(test_optim_logs, ~ {
  test_optim_log <- OptimLogger$new()
  
  test_optim <- optim(
    fn = build_gmm_g(
      candidate_data, dordered = dordered, test_summaries, 20e3, 50, use_vcov_moments = FALSE, group_vaccines_by = group_vaccines_by, # sim_seed = 123,
      fixed_model_probs = true_param,
      calculate_objective = TRUE,
      logger = test_optim_log,
      verbose = FALSE
    ),
    par = c(poverall = 0.5, psubcat = 0.5, 
            pvector = 0.5, psubunit = 0.5, prna = 0.5, pdna = 0.5, pinactivated = 0.5,
            pphase1 = 0.5, pphase2 = 0.5, pphase3 = 0.5) %>% qlogis(),
    control = lst(reltol = 0.00001)
    # method = "BFGS"
  )
  
  test_optim_log$data %>%
    mutate(run_id = .x, step = seq(n())) 
}) #, .options = furrr_options(seed = TRUE))

# Plot test optimization dynamics -----------------------------------------

test_optim_data %>%
  pivot_longer(-c(step, run_id), names_to = "param_name", values_to = "param_val") %>% 
  ggplot() +
  geom_line(aes(step, param_val, color = factor(run_id)), alpha = 0.75, show.legend = FALSE) +
  geom_hline(aes(yintercept = value), 
             linetype = "dotted",
             data = . %>% semi_join(enframe(true_param, name = "param_name"), ., by = "param_name")) +
  facet_wrap(vars(param_name), scales = "free") +
  labs(y = "") +
  theme_minimal()

# Plot Success Rates for Solutions ----------------------------------------

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

gmm_g <- build_gmm_g(candidate_data, 10e3, 50, use_vcov_moments = FALSE, group_vaccines_by = group_vaccines_by, sim_seed = 123,
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
