#SynthDid for Truncated History Test

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
#Code for Truncated History Test for SynthDid Paper

mc_estimate = function(Y, N0, T0) {
  N1=nrow(Y)-N0
  T1=ncol(Y)-T0
  W <- outer(c(rep(0,N0),rep(1,N1)),c(rep(0,T0),rep(1,T1)))
  mc_pred <- mcnnm_cv(Y, 1-W, num_lam_L = 20)
  mc_fit  <- mc_pred$L + outer(mc_pred$u, mc_pred$v, '+')
  mc_est <- sum(W*(Y-mc_fit))/sum(W)
  mc_est
}
mc_placebo_se = function(Y, N0, T0, replications=200) {
  N1 = nrow(Y) - N0
  theta = function(ind) { mc_estimate(Y[ind,], length(ind)-N1, T0) }
  sqrt((replications-1)/replications) * sd(replicate(replications, theta(sample(1:N0))))
}

difp_estimate = function(Y, N0, T0) {
  synthdid_estimate(Y, N0, T0, weights=list(lambda=rep(1/T0, T0)), eta.omega=1e-6)
}

sc_estimate_reg = function(Y, N0, T0) {
  sc_estimate(Y, N0, T0, eta.omega=((nrow(Y)-N0)*(ncol(Y)-T0))^(1/4))
}
difp_estimate_reg = function(Y, N0, T0) {
  synthdid_estimate(Y, N0, T0, weights=list(lambda=rep(1/T0, T0)))
}


estimators = list(did=did_estimate,
                  sc=sc_estimate,
                  sdid=synthdid_estimate,
                  difp=difp_estimate,
                  sc_reg = sc_estimate_reg)


#Load California data
data('california_prop99')
truncated_hist_test = list()

#---- Table 1 Code ----

#Choose for how many years to apply the SynthDiD
start_years = 1971:1974
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

#---- Table 2 Code ----

#Load California data
data('california_prop99')
truncated_hist_test_table2 = list()

#Choose for how many years to apply the SynthDiD
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

# Combine all the results into a single data frame
combined_results <- do.call(rbind, truncated_hist_test_table2)

min(combined_results$DID)
max(combined_results$DID)

#Leave one out part of Table 2

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


#---- Experimentation Code ----  

#Load California data

setup = panel.matrices(california_prop99)

tau.sdid = synthdid_estimate(setup$Y, setup$N0, setup$T0, eta.omega=1e-6)
print(summary(tau.sdid))

tau.difp = difp_estimate(setup$Y, setup$N0, setup$T0)
print(summary(tau.difp))

tau.difp2 = difp_estimate2(setup$Y, setup$N0, setup$T0)
print(summary(tau.difp2))

tau.difpreg = difp_estimate_reg(setup$Y, setup$N0, setup$T0)
print(summary(tau.difp))

tau.sc = sc_estimate(setup$Y, setup$N0, setup$T0, eta.omega = ((nrow(setup$Y) - setup$N0) * (ncol(setup$Y) - setup$T0))^(1/4) )
print(summary(tau.sc))

sc_estimate_reg = function(Y, N0, T0) {
  sc_estimate(Y, N0, T0, eta.omega=((nrow(Y)-N0)*(ncol(Y)-T0))^(1/4))}

tau_screg = sc_estimate_reg(setup$Y, setup$N0, setup$T0)
print(summary(tau_screg))




california_trun = california_prop99[california_prop99$Year >= 1974, ]
setup = panel.matrices(california_trun)

tau.hat = synthdid_estimate(setup$Y, setup$N0, setup$T0, eta.omega=1e-6)
print(summary(tau.hat))

tau.sc = sc_estimate(setup$Y, setup$N0, setup$T0)
print(summary(tau.sc))

california_out = california_prop99[california_prop99$State != "Utah",]
setup = panel.matrices(california_out)

tau.hat = synthdid_estimate(setup$Y, setup$N0, setup$T0)
print(summary(tau.hat))

tau.sc = sc_estimate(setup$Y, setup$N0, setup$T0, eta.omega = ((nrow(setup$Y) - setup$N0) * (ncol(setup$Y) - setup$T0))^(1/4) )
print(summary(tau.sc))

synthdid_controls(tau.hat)

#Plot conneticut 

# Filter for California and Connecticut
df_subset <- subset(california_prop99, State %in% c("California", "Connecticut", "Utah", "Montana", "Nevada"))

# Create plot
plot(california_prop99$Year[california_prop99$State == "California"],
     california_prop99$PacksPerCapita[california_prop99$State == "California"],
     type = "l", col = "blue", lwd = 2,
     ylim = c(50,200),
     xlab = "Year", ylab = "Packs Per Capita",
     main = "Cigarette Consumption: California vs. Connecticut")

# Add Connecticut line
lines(california_prop99$Year[california_prop99$State == "Connecticut"],
      california_prop99$PacksPerCapita[california_prop99$State == "Connecticut"],
      col = "red", lwd = 2)
lines(california_prop99$Year[california_prop99$State == "Utah"],
      california_prop99$PacksPerCapita[california_prop99$State == "Utah"],
      col = "red", lwd = 2)
lines(california_prop99$Year[california_prop99$State == "Montana"],
      california_prop99$PacksPerCapita[california_prop99$State == "Montana"],
      col = "red", lwd = 2)
lines(california_prop99$Year[california_prop99$State == "Nevada"],
      california_prop99$PacksPerCapita[california_prop99$State == "Nevada"],
      col = "red", lwd = 2)
lines(california_prop99$Year[california_prop99$State == "New Mexico"],
      california_prop99$PacksPerCapita[california_prop99$State == "New Mexico"],
      col = "red", lwd = 2)
lines(california_prop99$Year[california_prop99$State == "Idaho"],
      california_prop99$PacksPerCapita[california_prop99$State == "Idaho"],
      col = "red", lwd = 2)
lines(california_prop99$Year[california_prop99$State == "North Carolina"],
      california_prop99$PacksPerCapita[california_prop99$State == "North Carolina"],
      col = "red", lwd = 2)





abline(v = 1989, lty = 2, col = "darkgray")


# Add legend
legend("topright", legend = c("California", "Connecticut"),
       col = c("blue", "red"), lty = 1, lwd = 2)




#if I want to do something with confidence intervals
setup = panel.matrices(california_prop99)
tau.hat = synthdid_estimate(setup$Y, setup$N0, setup$T0)
se = sqrt(vcov(tau.hat, method='placebo'))
sprintf('point estimate: %1.2f', tau.hat)
sprintf('95%% CI (%1.2f, %1.2f)', tau.hat - 1.96 * se, tau.hat + 1.96 * se)
plot(tau.hat)




