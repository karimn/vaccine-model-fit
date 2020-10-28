# Clearing workspace, loading packages, setting working directory
rm(list=ls())

library("ggplot2")
library("data.table")
library("tidyverse")

setwd("~/Documents/covidIncentives/ShinyApps/AHTIDdesign")

# Install and load vaccineEarlyInvest package
devtools::install("../../../vaccineEarlyInvest")

library("vaccineEarlyInvest")

#' Monte Carlo runs of success
#' 
#' Runs Monte Carlo runs of the vaccine candidate success model.
#'
#' @param candidateFile String with the name of the file with data on candidates
#' @param replications Number of Monte Carlo replications
#' @param poverall Probability that no problem at the overall level prevents vaccine feasibility
#' @param psubcat Probability that there's no problem at subcategory level
#' @param pvector Probability that there's no problem at the viral vector platform level
#' @param psubunit Probability that there's no problem at the protein subunit platform level
#' @param prna Probability that there's no problem at the RNA platform level
#' @param pdna Probability that there's no problem at the DNA platform level
#' @param pattenuated Probability that there's no problem at the live attenuated platform level
#' @param pinactivated Probability that there's no problem at the inactivated platform level
#' @param ppreclinical Probability that there's no problem at the candidate level when a vaccine is in preclincal trials
#' @param pphase1 Probability that there's no problem at the candidate level when a vaccine is in phase 1 trials
#' @param pphase2 Probability that there's no problem at the candidate level when a vaccine is in phase 2 trials
#' @param pphase3 Probability that there's no problem at the candidate level when a vaccine is in phase 3 trials
#' @param maxcand Number of vaccine candidates to consider
#' @param seed Seed for the random number generater
#' @param ... Additional parameters to pass on to `Parameters` constructor
#' 
#' @return `data.table` With successes for each candidate in every run
#' @export
#'
#' @examples
getDraws <- function(candidateFile, replications=10000, dordered, 
                     poverall, psubcat, pvector, psubunit, prna, pdna, pattenuated, pinactivated,
                     ppreclinical, pphase1, pphase2, pphase3, maxcand=50, seed=1, ...) {
  
  t0 <- proc.time()
  
  par <- Parameters$new(replications=replications, poverall=poverall, psubcat=psubcat, pvector=pvector,
                        psubunit, prna=prna, pdna=pdna, pattenuated=pattenuated, pinactivated=pinactivated,
                        ppreclinical=ppreclinical, pphase1=pphase1, pphase2=pphase2, pphase3=pphase3,
                        maxcand=maxcand, ...)
  
  t1 <- proc.time()
  
   # dordered is now generated based on par0
  
  t2 <- proc.time()

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
  
  dplatforms <- unique(dordered[, .(Platform, pplat)])
  setkey(dplatforms, Platform, pplat)
  
  dordered <- dordered[,1:11]
  dcandidate <- copy(dordered)
  
  t3 <- proc.time()
  
  dcanddraws <- candidateDraws(dcandidate, par, seed=seed)
  
  t4 <- proc.time()
  
  # print((t1 - t0)["elapsed"])
  # print((t2 - t1)["elapsed"])
  # print((t3 - t2)["elapsed"])
  # print((t4 - t3)["elapsed"])
  
  return(dcanddraws)
}

# Karim: Lines 97-103 only need to be run once, and dordered can be passed as an argument to getDraws
par0 <- Parameters$new(maxcand=50)
d <- loadData(par0, candidateFile="Data/vaccinesSummaryOct2.csv")
d$Target <- "Other"
d$Target[1:5]<-"Spike"
d$Target[10:15]<-"Recombinant"
dordered <- candidatesFung(d, par0)$dordered

draws <- getDraws(candidateFile="Data/vaccinesSummaryOct2.csv", replications=20000, dordered=dordered,
                  poverall=0.9, psubcat=0.9, pvector=0.8, psubunit=0.8, prna=0.6, pdna=0.4, pattenuated=0.8, 
                  pinactivated=0.8, ppreclinical=0.14, pphase1=0.23, pphase2=0.32, pphase3=0.5)


