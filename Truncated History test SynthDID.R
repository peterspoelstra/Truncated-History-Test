#SynthDid for Truncated History Test

install.packages("synthdid")
devtools::install_github("synth-inference/synthdid")

install.packages("devtools")
install.packages("latex2exp")
library(devtools) 
install_github("susanathey/MCPanel")


library(synthdid)
library(MCPanel)
library(dplyr)

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
                  mc = mc_estimate,
                  sc_reg = sc_estimate_reg)


#Load California data
data('california_prop99')
truncated_hist_test = list()

#Choose for how many years to apply the SynthDiD
start_years = 1971:1977
for (start_year in start_years) {
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
  
  round(california.table, digits=1)
  truncated_hist_test[[start_year]] = california.table
  
}


#if I want to do something with confidence intervals
setup = panel.matrices(california_prop99)
tau.hat = synthdid_estimate(setup$Y, setup$N0, setup$T0)
se = sqrt(vcov(tau.hat, method='placebo'))
sprintf('point estimate: %1.2f', tau.hat)
sprintf('95%% CI (%1.2f, %1.2f)', tau.hat - 1.96 * se, tau.hat + 1.96 * se)
plot(tau.hat)






