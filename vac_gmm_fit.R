library(data.table)
library(magrittr)
library(tidyverse)
library(gmm)

library(vaccineEarlyInvest)

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
  
  # candidate_data %>% 
  #   group_nest(Platform, Subcategory, Target, .key = "candidates") %>% 
  #   mutate(vacc_group_id = seq(n()))
  
  return(candidate_data)
}

get_candidate_draws <- function(candidate_data, replications, 
                                poverall, psubcat, pvector, psubunit, prna, pdna, pattenuated, pinactivated, ppreclinical, pphase1, pphase2, pphase3, 
                                maxcand) {
    par <- Parameters$new(replications = replications, 
                          poverall = poverall, psubcat = psubcat, pvector = pvector, psubunit = psubunit, prna = prna, pdna = pdna, pattenuated = pattenuated, 
                          pinactivated = pinactivated, ppreclinical = ppreclinical, pphase1 = pphase1, pphase2 = pphase2, pphase3 = pphase3, maxcand = maxcand)  
    # dordered <- unnest(candidate_data, candidates) %>% 
      # as.data.table() %>% 
    dordered <- candidate_data %>% 
      candidatesFung(par) %>% 
      pluck("dordered")
    
    dplatforms <- unique(dordered[, .(Platform, pplat)])
    setkey(dplatforms, Platform, pplat)
  
    dordered <- dordered[,1:11]
    dcandidate <- copy(dordered)
    
    vacc_group_ids <- candidate_data %>% 
      group_by(Platform, Subcategory) %>% 
      summarize(vacc_group_id = cur_group_id(), .groups = "drop")
   
    draws <- candidateDraws(dcandidate, par, seed = NULL) %>% 
      # left_join(select(candidate_data, -candidates)) %>% # Get the vaccine group IDs
      left_join(vacc_group_ids, by = setdiff(names(vacc_group_ids), "vacc_group_id")) %>% # Get the vaccine group IDs
      group_by(r, vacc_group_id, .drop = FALSE) %>% 
      summarize(any_success = sum(success) > 0, 
                success_rate = mean(success),
                .groups = "drop") %>%  # Any success in group
      complete(r, vacc_group_id = seq(max(vacc_group_ids$vacc_group_id)), 
               fill = lst(any_success = FALSE, success_rate = 0)) # Make sure all groups are present in the draws
}

summarize_draws <- function(draws) {
    success_rates <- draws %>% 
      group_by(vacc_group_id) %>% 
      summarize(
        # success_rate = mean(any_success),
        success_rate = mean(success_rate),
        success_rate_pred = qlogis(success_rate),
        
        .groups = "drop"
      ) %>% 
      ungroup()
    
    success_corr <- draws %>% 
      # pivot_wider(id_cols = r, names_from = vacc_group_id, values_from = any_success) %>% 
      pivot_wider(id_cols = r, names_from = vacc_group_id, values_from = success_rate) %>% 
      select(-r) %>% 
      as.data.frame() %>% 
      pcaPP::cor.fk() %>% # Faster than base cor()
      magrittr::extract(lower.tri(.)) # Get lower triangle
    
    lst(success_rates, success_corr)
}

build_gmm_g <- function(candidate_data, replications, maxcand, use_corr_moments = FALSE) {
  function(pbeta, cgd_trials) {
    model_probs <- plogis(pbeta) %>% 
      set_names(c("poverall", "psubcat", "pvector", "psubunit", "prna", "pdna", "pattenuated", "pinactivated", 
                  "ppreclinical", "pphase1", "pphase2", "pphase3"))
    
    print(unname(model_probs))
    
    draws <- rlang::exec(get_candidate_draws, candidate_data, replications = replications, !!!model_probs, maxcand = maxcand)  
 
    model_summaries <- summarize_draws(draws) 
    
    # browser()
 
    moments <- full_join(model_summaries$success_rates, cgd_trials$success_rates, by = "vacc_group_id") %>%
      mutate_at(vars(starts_with("success_rate")), coalesce, 0) %>%
      transmute(vacc_group_id, 
                success_moment = success_rate.x - success_rate.y) %>% 
      # print() %>% 
      # complete(vacc_group_id = seq(max(candidate_data$vacc_group_id)), fill = lst(success_moment = 0)) %>%  # Make sure all groups are present in the draws
      pull(success_moment)
   
    if (use_corr_moments) { 
      known_corr <- !is.na(cgd_trials$success_corr)
      
      moments %<>% 
        c(model_summaries$success_corr[known_corr] - cgd_trials$success_corr[known_corr]) %>%
        coalesce(0.0)
    }
    
    # print(moments)
      
    return(moments)
   
    # full_join(cgd_trials, model_summaries$success_rates, by = "vacc_group_id") %>% 
    #   mutate_at(vars(any_success, success_rate), coalesce, 0) %>% 
    #   transmute(r, vacc_group_id, m = any_success - success_rate) %>% 
    #   pivot_wider(names_from = vacc_group_id, values_from = m) %>% 
    #   select(-r) %>% 
    #   mutate_all(coalesce, 0) 
  }
}

# Test --------------------------------------------------------------------

candidate_data <- load_candidate_data(file.path("data", "vaccinesSummaryOct2.csv"))

test_draws <- get_candidate_draws(
  candidate_data, replications = 10000,
  poverall=0.9, psubcat=0.9, pvector=0.8, psubunit=0.8, prna=0.6, pdna=0.4, pattenuated=0.8, pinactivated=0.8, ppreclinical=0.14, pphase1=0.23, pphase2=0.32, pphase3=0.5,
  maxcand = 50
)

test_summaries <- summarize_draws(test_draws)

test_results <- gmm(build_gmm_g(candidate_data, 10e3, 50, use_corr_moments = FALSE), 
                    test_summaries, 
                    t0 = rep(0, 12)) # rnorm(12)) 
                    # type = "iterative") 
