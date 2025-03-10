#Test file for truncated history test
library(Synth)
library(data.table)
library(SCtools)
rm(list = ls());gc()

#Three Abadie paper cases:

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



#---- Abadie 2003 Basque ----
data("basque")
basque = as.data.table(basque)
basque[, regionno_numeric := as.numeric(regionno)]
basque[, regionno := as.character(regionno)]

basque[, treatment := ifelse(regionno_numeric == 17, 1,0)]

treated_unit <- basque[treatment == 1, unique(regionno_numeric)]
control_unit <- basque[treatment == 0, unique(regionno_numeric)]


scm_basque = scm_only_pretreatment(basque, y_variable = "gdpcap", unit = "regionno", time_variable = "year",
                                   pretreatment_period = 1959:1969, post_treatment_period = 1970:1997)

scm_basque$average_treatment

#In-Time Placebo test:
placebo_res = generate.placebos(scm_basque$dataprep_out, scm_basque$synth_out, Sigf.ipop = 3)
plot_placebos(placebo_res, discard.extreme = FALSE)
p_value = mspe.test(placebo_res, discard.extreme = FALSE)
mspe.plot(placebo_res, discard.extreme = FALSE)
mspe.plot(placebo_res, discard.extreme = TRUE, mspe.limit = 20)




#---- Abadie 2010 Tobocco ----

#Load data California
california_data_url = "https://raw.githubusercontent.com/causalify-code/synthetic-control-replications/refs/heads/main/california-tobacco-control-program/data-abadie-diamond-hainmueller-california-proposition-99-tobacco-control-program-synthetic-control.csv"
california_dataframe = read.csv(california_data_url)

#Replicate original paper results:

#Reproducing the paper
dataprep_out <- dataprep(
  foo = california_dataframe,
  predictors = c("lnincome", "age15to24", "beer", "retprice"),
  special.predictors = list(
    list("cigsale", 1975, c("mean")),
    list("cigsale", 1980, c("mean")),
    list("cigsale", 1988, c("mean"))
  ),
  dependent = "cigsale",
  unit.variable = "state",
  time.variable = "year",
  treatment.identifier = "California",
  controls.identifier = unique(california_dataframe$state_name[-which(california_dataframe$state_name == "California")]),
  time.predictors.prior = 1970:1988,
  time.optimize.ssr = 1970:1988,
  time.plot = 1970:2000,
  unit.names.variable = "state_name"
) 


synth_out <- synth(
  data.prep.obj = dataprep_out
)

# FIGURE 2: REAL VS SYNTH

synth_california <- dataprep_out$Y0 %*% synth_out$solution.w
plot(
  1970:2000,
  # Synthetic California: (control units matrix) X (weights column)
  synth_california,
  type="l",
  lty="dashed", ylim=c(0,140))

lines(1970:2000,
      dataprep_out$Y1
)
abline(v=1988,lty="dotted")


# GAPS PLOT
gaps <- dataprep_out$Y1 - synth_california

plot(
  1970:2000,
  gaps,
  type="l"
  
)
abline(h=0, lty="dotted")

#Get placebo test
placebo_res = generate.placebos(dataprep_out, synth_out, Sigf.ipop = 3)
plot_placebos(placebo_res, discard.extreme = FALSE)
p_value = mspe.test(placebo_res, discard.extreme = FALSE)
mspe.plot(placebo_res, discard.extreme = FALSE)
mspe.plot(placebo_res, discard.extreme = TRUE, mspe.limit = 20)

#Code for Truncated History Test
california_dataframe = as.data.table(california_dataframe)
california_dataframe[, state_numeric := as.numeric(state)]
california_dataframe[, state := as.character(state)]

california_dataframe[, treatment := ifelse(state_numeric == 3, 1,0)]

treated_unit <- california_dataframe[treatment == 1, unique(state_numeric)]
control_unit <- california_dataframe[treatment == 0, unique(state_numeric)]

result_TruncHist = data.frame(Unit = numeric(), ate = numeric(), significance = numeric())
years = 1970:1980
for (year in years) {
  scm_tobacco = scm_only_pretreatment(california_dataframe, y_variable = "cigsale", unit = "state",
                                      time_variable = "year", pretreatment_period = year:1988, post_treatment_period = 1989:2000)
  ate_tobacco = scm_tobacco$average_treatment
  placebo_res = generate.placebos(scm_tobacco$dataprep_out, scm_tobacco$synth_out, Sigf.ipop = 3)
  p_statistic = mspe.test(placebo_res, discard.extreme = FALSE)
  p_value = p_statistic$p.val
  result_TruncHist = rbind(result_TruncHist, data.frame(Unit = year, ate = ate_tobacco,  significance = p_value ))
}

scm_tobacco = scm_only_pretreatment(california_dataframe, y_variable = "cigsale", unit = "state",
                                    time_variable = "year", pretreatment_period = 1980:1988, post_treatment_period = 1989:2000)
scm_tobacco$average_treatment



#In-Space Placebo test:
placebo_res = generate.placebos(scm_tobacco$dataprep_out, scm_tobacco$synth_out, Sigf.ipop = 3)
plot_placebos(placebo_res, discard.extreme = FALSE)
p_value = mspe.test(placebo_res, discard.extreme = FALSE)

mspe.plot(placebo_res, discard.extreme = FALSE)
mspe.plot(placebo_res, discard.extreme = TRUE, mspe.limit = 20)



#---- Abadie 2015 German Reunification ----
german_reunification_url = "https://raw.githubusercontent.com/OscarEngelbrektson/SyntheticControlMethods/refs/heads/master/examples/datasets/german_reunification.csv"
german_dataframe = read.csv(german_reunification_url) 

german_dataframe = as.data.table(german_dataframe)
german_dataframe[, code_numeric := as.numeric(code)]
german_dataframe[, code := as.character(code)]

german_dataframe[, treatment := ifelse(code_numeric == 7, 1,0)]

treated_unit <- german_dataframe[treatment == 1, unique(code_numeric)]
control_unit <- german_dataframe[treatment == 0, unique(code_numeric)]


scm_german = scm_only_pretreatment(german_dataframe, y_variable = "gdp", unit = "code",
                                   time_variable = "year", pretreatment_period = 1961:1989, post_treatment_period = 1990:2003)
scm_german$average_treatment

#In-Time Placebo test:
placebo_res = generate.placebos(scm_german$dataprep_out, scm_german$synth_out, Sigf.ipop = 3)
plot_placebos(placebo_res, discard.extreme = FALSE)
p_value = mspe.test(placebo_res, discard.extreme = FALSE)
mspe.plot(placebo_res, discard.extreme = FALSE)
mspe.plot(placebo_res, discard.extreme = TRUE, mspe.limit = 20)


