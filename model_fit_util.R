param_names <- c("poverall", "psubunit", "prna", "pdna", "pattenuated", "pinactivated", "ppreclinical", "pphase1", "pphase2", "pphase3")

load_candidate_data <- function(data_file) {
  candidate_data <- loadData(par = NULL, data_file)
  
  candidate_data$Target <- "Other"
  candidate_data$Target[1:5]<-"Spike"
  candidate_data$Target[10:15]<-"Recombinant"
  
  return(candidate_data)
}

calc_dordered <- function(candidate_data, maxcand) {
  par0 <- Parameters$new(maxcand = maxcand)
  
  candidatesFung(candidate_data, par0)$dordered
}

get_candidate_draws <- function(candidate_data, replications, dordered, 
                                param = c(poverall=0.9, psubcat=0.9, pvector=0.8, psubunit=0.8, prna=0.6, pdna=0.4, pattenuated=0.8, pinactivated=0.8, ppreclinical=0.14, pphase1=0.23, pphase2=0.32, pphase3=0.5),
                                maxcand,
                                group_vaccines_by = vars(Platform, Subcategory),
                                seed = NULL) {
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
  
  draws <- withr::with_seed(seed, candidateDraws(dcandidate, par, seed = seed))
    
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

get_summary_success_rates <- function(draws) {
  draws %>% 
    group_by(vacc_group_id) %>% 
    summarize(
      success_rate = mean(success_rate),
      .groups = "drop"
    ) %>% 
    ungroup()
}

summarize_draws <- function(draws) {
    success_rates <- get_summary_success_rates(draws) 
    
    success_vcov <- draws %>% 
      pivot_wider(id_cols = r, names_from = vacc_group_id, values_from = success_rate) %>% 
      select(-r) %>% 
      as.data.frame() %>% 
      cov() %>% 
      magrittr::extract(lower.tri(.)) # Get lower triangle
    
    lst(success_rates, success_vcov)
}

calculate_draws_vcov <- function(draws) {
  draws %>% 
    pivot_wider(id_cols = r, names_from = vacc_group_id, values_from = success_rate) %>% 
    select(-r) %>% 
    cov() 
}

OptimLogger <- R6Class(
  "OptimLogger", 
   list(
     data = tibble(),
     draw_summary = tibble(),
     log_param = function(...) { self$data %<>% bind_rows(tibble(!!!rlang::list2(...))) },
     log_summary = function(summary) { self$draw_summary <- summary }
   )
)

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
                        weighting_matrix = NULL, 
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
    
    draws <- get_candidate_draws(candidate_data, replications = replications, dordered = dordered, param = all_model_probs, maxcand = maxcand, group_vaccines_by = group_vaccines_by, seed = sim_seed) %>% 
      semi_join(x$success_rates, by = "vacc_group_id") # Only use the vaccines in the target summaries

    model_summaries <- summarize_draws(draws) 
    
    moments <- calculate_moments(model_summaries, x, use_vcov_moments = use_vcov_moments) 

   # moments <- full_join(x, model_summaries$success_rates, by = "vacc_group_id") %>%
   #    mutate_at(vars(starts_with("success_rate")), coalesce, 0) %>%
   #    transmute(r, vacc_group_id, m = success_rate.x - success_rate.y) %>%
   #    pivot_wider(names_from = vacc_group_id, values_from = m) %>%
   #    select(-r) %>%
   #    mutate_all(coalesce, 0) %>% 
   #    select(moments_to_use)
    
    if (is_null(weighting_matrix)) {
      weighting_matrix <- diag(length(moments))
    }
    
   objective <- calculate_objective(moments, weighting_matrix) 
    
    if (!is_null(logger)) {
      logger$log_param(!!!all_model_probs, obj = objective) 
      logger$log_summary(model_summaries)
    }
   
    if(calculate_objective) {
      return(objective)
    } else {
      return(moments)
    }
  }
}

optim_run <- function(run_id, summaries, maxcand, prev_run_data = NULL, replications = 20e3, use_vcov_moments = TRUE, 
                      initial_par = c(poverall = 0.5, psubcat = 0.5, 
                                      pvector = 0.5, psubunit = 0.5, prna = 0.5, pdna = 0.5, pinactivated = 0.5,
                                      pphase1 = 0.5, pphase2 = 0.5, pphase3 = 0.5),
                      fixed_model_probs = c(poverall=0.9, psubcat=0.9, 
                                            pvector=0.8, psubunit=0.8, prna=0.6, pdna=0.4, pattenuated=0.8, pinactivated=0.8, 
                                            ppreclinical=0.14, pphase1=0.23, pphase2=0.32, pphase3=0.5),
                      weighting_matrix = NULL,
                      ndeps = rep_along(initial_par, 1e-2)) {
  test_optim_log <- OptimLogger$new()
  run_seed <- as.integer(Sys.time()) %% 1e5
  
  if (!is_null(prev_run_data)) {
    initial_par <- prev_run_data %>% 
      filter(step == n()) %>% 
      select(all_of(names(initial_par))) %>% 
      unlist()
  }
  
  test_optim <- optim(
    fn = build_gmm_g(
      candidate_data, dordered = dordered, x = summaries, 
      replications = replications, 
      maxcand = maxcand, 
      use_vcov_moments = use_vcov_moments, 
      weighting_matrix = weighting_matrix,
      group_vaccines_by = group_vaccines_by, 
      sim_seed = run_seed,
      fixed_model_probs = fixed_model_probs,
      calculate_objective = TRUE,
      logger = test_optim_log,
      verbose = FALSE
    ),
    par = qlogis(initial_par), 
    control = lst(reltol = 0.0001, ndeps = ndeps, REPORT = 1),
    method = "BFGS"
  )
 
  run_data <- test_optim_log$data %>%
    mutate(
      run_id, 
      step = seq(n()),
      summary = map_if(seq(n()), step == n(), ~ test_optim_log$draw_summary, .else = ~ NULL) 
    )
  
  if (!is_null(prev_run_data)) {
    run_data %<>%
      mutate(step = step + nrow(prev_run_data)) %>% 
      bind_rows(prev_run_data, .)
  }
  
  return(run_data)
}

quick_get_summary <- function(param, ..., candidate_data, dordered, maxcand, group_vaccines_by = NULL, summary_type = "success_rates") {
  summ <- get_candidate_draws(
    candidate_data = candidate_data, replications = 3e5, dordered = dordered,
    param = param,
    maxcand = maxcand,
    group_vaccines_by = group_vaccines_by
  ) %>% 
    summarize_draws()
  
  if (!is_null(summary_type)) {
    return(pluck(summ, summary_type))
  } else {
    return(summ)
  }
}


convert_cgd_trials <- function(start_month_offset, cgd_trials, id_dict) {
  cgd_trials %>% 
    rename(r = try_id) %>% 
    mutate(success_rate = !is.na(phase_mon_approval) & phase_mon_approval <= start_month_offset) %>% 
    complete(r, vaccine_id = id_dict$cgd_vaccine_id) %>% 
    mutate(success_rate = coalesce(success_rate, 0L)) %>% 
    inner_join(select(id_dict, candInd, cgd_vaccine_id), by = c("vaccine_id" = "cgd_vaccine_id")) %>% # inner_ to exclude vaccines we are not considering (e.g. pre-clinical) 
    rename(vacc_group_id = candInd)
}

