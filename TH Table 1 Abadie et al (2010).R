#Test file for truncated history test
library(Synth)
library(data.table)
library(SCtools)
rm(list = ls());gc()


#---- Synth pretreatment function ---- 
scm_only_pretreatment <- function(dt, y_variable, unit, time_variable, pretreatment_period, post_treatment_period){
  unit_numeric = paste0(unit, "_numeric")
  years <- seq(min(pretreatment_period), max(pretreatment_period))
  outcome_lags_all <- lapply(years, function(year) {
    list(y_variable, year, "mean")
  })
  
  dataprep.out <- dataprep(
    foo = dt,
    special.predictors = outcome_lags_all,
    dependent = y_variable,
    unit.variable = unit_numeric,
    time.variable = time_variable,
    treatment.identifier = treated_unit,
    controls.identifier = control_unit,
    time.predictors.prior = pretreatment_period,
    time.optimize.ssr = pretreatment_period,
    unit.names.variable = unit,
    #time.plot = 1970:2000
    time.plot = min(pretreatment_period):max(post_treatment_period)
  )
  
  synth.out = synth(dataprep.out)
  synth.table = synth.tab(dataprep.res = dataprep.out, synth.res = synth.out)
  print(synth.table)
  
  path.plot(
    synth.res = synth.out,
    dataprep.res = dataprep.out,
    Ylab = "Y",
    Xlab = "year"
  )
  abline(v = max(pretreatment_period), col = "black", lwd = 2, lty = 2)
  
  gaps.plot(
    synth.res = synth.out,
    dataprep.res = dataprep.out,
    Ylab = "Y",
    Xlab = "year",
    Z.plot = FALSE,
  )
  abline(v = max(pretreatment_period), col = "black", lwd = 2, lty = 2)
  
  gaps = dataprep.out$Y1plot - dataprep.out$Y0plot %*% synth.out$solution.w
  names(gaps) <- seq(min(pretreatment_period), max(post_treatment_period))
  average_treatment = mean(gaps[as.character(post_treatment_period)])
  
  return(list(synth_out = synth.out, dataprep_out = dataprep.out, average_treatment = average_treatment))
  
}









#---- Load in the data ----
#Load data California
california_data_url = "https://raw.githubusercontent.com/causalify-code/synthetic-control-replications/refs/heads/main/california-tobacco-control-program/data-abadie-diamond-hainmueller-california-proposition-99-tobacco-control-program-synthetic-control.csv"
california_dataframe = read.csv(california_data_url)

#---- Truncated History Paper Table 1 Abadie part----
#Code for Truncated History Test
california_dataframe = as.data.table(california_dataframe)
california_dataframe[, state_numeric := as.numeric(state)]
california_dataframe[, state := as.character(state)]

california_dataframe[, treatment := ifelse(state_numeric == 3, 1,0)] #state_numeric 3 is california

treated_unit <- california_dataframe[treatment == 1, unique(state_numeric)]
control_unit <- california_dataframe[treatment == 0, unique(state_numeric)]

result_TruncHist = data.frame(Unit = numeric(), mspe = numeric(),ate = numeric(), significance = numeric())
years = 1970:1974
for (year in years) {
  scm_tobacco = scm_only_pretreatment(california_dataframe, y_variable = "cigsale", unit = "state",
                                      time_variable = "year", pretreatment_period = year:1988, post_treatment_period = 1989:2000)
  ate_tobacco = scm_tobacco$average_treatment
  mspe = scm_tobacco[["synth_out"]][["loss.v"]]
  placebo_res = generate.placebos(scm_tobacco$dataprep_out, scm_tobacco$synth_out, Sigf.ipop = 3, strategy = "multicore")
  p_statistic = mspe.test(placebo_res, discard.extreme = FALSE)
  p_value = p_statistic$p.val
  result_TruncHist = rbind(result_TruncHist, data.frame(Unit = year, mspe = mspe, ate = ate_tobacco,  significance = p_value ))
}

#The results are stored in the dataframe: results_TruncHist and they are presented in Table 1




