#Author Peter Spoelstra
#Code for (Generalized) Truncated History for Synthetic Control Approaches
#This code is for the paper Table 1, 2 and 3


#---- Install and load packages ----

# install.packages("synthdid")
# devtools::install_github("synth-inference/synthdid")
# install.packages("devtools")
# install.packages("latex2exp")
# install_github("susanathey/MCPanel")

library(devtools) 
library(synthdid)
library(MCPanel)
library(dplyr)
rm(list = ls())

#We use the package Synthdid from Arkhangelsky et al. (2021)

#Load California data
data('california_prop99')
truncated_hist_test = list()


#Add DIFP estimator, the rest is already in synthdid (see https://synth-inference.github.io/synthdid/articles/paper-results.html for more details)
difp_estimate = function(Y, N0, T0) {
  synthdid_estimate(Y, N0, T0, weights=list(lambda=rep(1/T0, T0)), eta.omega=1e-6)
}

estimators = list(did=did_estimate,
                  sc=sc_estimate,
                  sdid=synthdid_estimate,
                  difp=difp_estimate)

#---- Table 1 Code ----

#Choose for how many years to apply the Truncated History Test
start_years = 1970:1974
for (start_year in start_years) {
  start_year_char = as.character(start_year)
  california_trun = california_prop99[california_prop99$Year >= start_year, ]
  setup = panel.matrices(california_trun)
  
  estimates = lapply(estimators, function(estimator) { estimator(setup$Y, setup$N0, setup$T0) } )
  standard.errors = mapply(function(estimate, name) {
    set.seed(12345)
    if(name == 'mc') { mc_placebo_se(setup$Y, setup$N0, setup$T0) }
    else {             sqrt(vcov(estimate, method='placebo'))     }
  }, estimates, names(estimators))
  
  california.table = rbind(unlist(estimates), unlist(standard.errors))
  rownames(california.table) = c('estimate', 'standard error')
  colnames(california.table) = toupper(names(estimators))
  
  california_df = as.data.frame(california.table)
  
  #Save results in list
  truncated_hist_test[[start_year_char]] = california_df
  
}
#The results are stored in truncated_hist_test and presented in Table 1 of the paper 

#---- Table 2 Code ----

#Load California data
data('california_prop99')
truncated_hist_test_table2 = list()

#Leave-one-time-period-out
pretreatment_years = 1970:1988  #1970:1988
for (exclude_year in pretreatment_years) {
  
  exclude_year_char = as.character(exclude_year)
  california_trun = california_prop99[california_prop99$Year != exclude_year, ]
  setup = panel.matrices(california_trun)
  
  estimates = lapply(estimators, function(estimator) { estimator(setup$Y, setup$N0, setup$T0) } )

  california.table = rbind(unlist(estimates))
  rownames(california.table) = c('estimate')
  colnames(california.table) = toupper(names(estimators))
  
  california_df = as.data.frame(california.table)
  
  #Save results in list
  truncated_hist_test_table2[[exclude_year_char]] = california_df
}

#The results in Table 2 of the leave-one-time-period out are stored in truncated_hist_test_table2 
#the code below extract them. We use the same setup for leave-two-time periods-out,leave-one-unit out and right and left truncation

# Combine all the results into a single data frame
combined_results <- do.call(rbind, truncated_hist_test_table2)

min(combined_results$DID)
max(combined_results$DID)
mean(combined_results$DID)

#Leave-two-time-periods-out
#Load California data
data('california_prop99')
truncated_hist_test_table2 = list()

pretreatment_years = 1970:1988  #1970:1988
year_combinations = combn(pretreatment_years, 2, simplify = FALSE)

# Loop over each pair of years to exclude
for (exclude_years in year_combinations) {
  
  combo_label = paste(exclude_years, collapse = "_")
  california_trun = california_prop99[!california_prop99$Year %in% exclude_years, ]
  setup = panel.matrices(california_trun)
  
  estimates = lapply(estimators, function(estimator) {
    estimator(setup$Y, setup$N0, setup$T0)
  })
  
  california.table = rbind(unlist(estimates))
  rownames(california.table) = c('estimate')
  colnames(california.table) = toupper(names(estimators))
  
  california_df = as.data.frame(california.table)
  
  # Save results in list
  truncated_hist_test_table2[[combo_label]] = california_df
}

# Combine all the results into a single data frame
combined_results <- do.call(rbind, truncated_hist_test_table2)

#Do for every model
min(combined_results$SC)
max(combined_results$SC)
mean(combined_results$SC)


#Right truncation part
truncated_hist_right = list()

#Choose for how many years to apply the SynthDiD
right_years = 1984:1987
for (right_year in right_years) {
  
  right_yearr_char = as.character(right_year)
  california_trun = california_prop99[california_prop99$Year <= right_year | california_prop99$Year > 1988 , ]
  setup = panel.matrices(california_trun)
  
  estimates = lapply(estimators, function(estimator) { estimator(setup$Y, setup$N0, setup$T0) } )

  california.table = rbind(unlist(estimates))
  rownames(california.table) = c('estimate')
  colnames(california.table) = toupper(names(estimators))
  
  california_df = as.data.frame(california.table)
  
  #Save results in list
  truncated_hist_right[[right_yearr_char]] = california_df
  
}

combined_results <- do.call(rbind, truncated_hist_right)

mean(combined_results$DIFP)
min(combined_results$DIFP)
max(combined_results$DIFP)


#Leave one unit out part of Table 2

#Load California data
data('california_prop99')
truncated_hist_test_table3 = list()

#Choose for how many years to apply the SynthDiD
states = unique(california_prop99$State)
states = setdiff(states, "California")
for (state in states) {
  
  exclude_state_char = as.character(state)
  california_trun = california_prop99[california_prop99$State != state, ]
  setup = panel.matrices(california_trun)
  
  estimates = lapply(estimators, function(estimator) { estimator(setup$Y, setup$N0, setup$T0) } )
  
  california.table = rbind(unlist(estimates))
  rownames(california.table) = c('estimate')
  colnames(california.table) = toupper(names(estimators))
  
  california_df = as.data.frame(california.table)
  
  #Save results in list
  truncated_hist_test_table3[[exclude_state_char]] = california_df
  
}

# Combine all the results into a single data frame
combined_results3 <- do.call(rbind, truncated_hist_test_table3)

mean(combined_results$DIFP)
min(combined_results$DIFP)
max(combined_results$DIFP)


#---- Table 3 Code ----

#The sdid possible combinations from table 3
sdid1 = function(Y, N0, T0) {
  synthdid_estimate(setup$Y, setup$N0, setup$T0, omega.intercept = TRUE)
}
sdid2 = function(Y, N0, T0) {
  synthdid_estimate(setup$Y, setup$N0, setup$T0, eta.omega=1e-6, omega.intercept = TRUE) #eta.omega = 1e-6 is zero
}
sdid3 = function(Y, N0, T0) {
  synthdid_estimate(setup$Y, setup$N0, setup$T0, omega.intercept = FALSE)
}
sdid4 = function(Y, N0, T0) {
  synthdid_estimate(setup$Y, setup$N0, setup$T0, eta.omega=1e-6, omega.intercept = FALSE)
}

estimators = list(n1=sdid1,
                  n2=sdid2,
                  n3=sdid3,
                  n4=sdid4)



#Choose for how many years to apply the Truncated History Robustness Check
truncated_hist_test_sdid = list()
start_years = 1970:1974

for (start_year in start_years) {
  
  start_year_char = as.character(start_year)
  california_trun = california_prop99[california_prop99$Year >= start_year, ]
  setup = panel.matrices(california_trun)
  
  estimates = lapply(estimators, function(estimator) { estimator(setup$Y, setup$N0, setup$T0) } )
  # standard.errors = mapply(function(estimate, name) {
  #   set.seed(12345)
  #   sqrt(vcov(estimate, method='placebo'))
  # }, estimates, names(estimators))
  
  positive.weights = sapply(estimates, function(estimate) {
    weights = attr(estimate, "weights")$omega
    sum(weights > 0)
  })
  
  positive.weights_1 = sapply(estimates, function(estimate) {
    weights = attr(estimate, "weights")$omega
    sum(weights > 0.01)
  })
  
  
  california.table = rbind(
    unlist(estimates),
    #unlist(standard.errors),
    unlist(positive.weights),
    unlist(positive.weights_1)
  )
  
  #rownames(california.table) = c('estimate', 'standard error', "pos_weight_0", "pos_weight_1")
  rownames(california.table) = c('estimate', "pos_weight_0", "pos_weight_1")
  colnames(california.table) = toupper(names(estimators))
  
  california_df = as.data.frame(california.table)
  
  #Save results in list
  truncated_hist_test_sdid[[start_year_char]] = california_df
  
}

#The information stored in truncated_hist_test_sdid is presented in Table 3

