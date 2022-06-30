
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

logit = function(x) {
  log(x / (1 - x))
}
inv_logit = function(x) {
  exp(x) / (1 + exp(x))
}

log_unif_prec = function(x, a, b) {
  lower = log(1 / b^2)
  upper = log(1 / a^2)
  log_dens = c(rep(- 1000, sum(x < lower)), - log(2) - log(b - a) - 1.5 * x[x >= lower & x <= upper], rep(- 1000, sum(x > upper)))
  log_jacobian = x
  return(log_dens + log_jacobian)
}

#gamma: half-t scale parameter
#ni: half-t dof
#x: threshold value

EPC = function(gamma, ni, x) {
  #loading actuar package
  #require(actuar)
  #defining the function
  - log(1 - actuar:::betaint((x^2 / (x^2 + gamma^2 * ni)), 1 / 2, ni / 2) / beta(1 / 2, ni / 2)) / x
}


#m: number of cohorts
#N: max number of patients enrolled per arm (same for each arm)
#p_true: response rate for generating the data
#ia1_fraction: define 1st IA as fraction of N, as numeric

# Ale:
# Ho aggiornato la funzione trial;
# ora p_true va passato come vettore e quindi possono essere
# gestiti tutti i differenti scenari.
# Lo scenario è quello che viene impostato in sim_basket
# e alcuni valori di input di operating_char
# (i.e., m, N, p_true, p_null, p_target) dovranno essere gli
# stessi usati per generare i dati
# La scelta spetta all'utente in base all'ordine in cui
# inserisce le probabilità di successo in p_true.

trial = function(m, N, p_true, ia1_fraction) {
  #first interim analysis
  ia1 = ceiling(N * ia1_fraction)

  #successive interim analyses
  step = ceiling(ia1 * m / 2)

  #simulating arms' responses
  yy = rbinom(rep(1, m), N, p_true)
  ee = NULL
  for (ii in 1:m) {
    ee = c(ee, c(rep(1, yy[ii]), rep(0, (N - yy[ii]))))
  }

  #randomizing the patients
  pp = sample(1:(N*m), (N*m))

  #arms indicator
  arms = NULL
  for (ii in 1:m) {
    arms = c(arms, rep(ii, N))
  }

  #patients assignments
  aa = arms[pp]
  delta = matrix(0, nrow = m, ncol = (N*m))
  for (ii in 1:(N*m)) {
    delta[aa[ii], ii] = 1
  }

  #patients responses
  responses = ee[pp]

  #finding the first interim analysis times
  first = max(sapply(1:m, function(y) {
    which(cumsum(aa == y) == ia1)[1]
  }))

  #first interim analysis accruals
  ia1_accruals = apply(delta[, 1:first], 1, sum)

  #returning the object
  return(list(assignments = delta, responses = responses, ia_accruals = as.data.frame(ia1_accruals)))
}

# OLD
# trial = function(m, N, p_true, ia1_fraction) {
#   #first interim analysis
#   ia1 = ceiling(N * ia1_fraction)
#
#   #successive interim analyses
#   step = ceiling(ia1 * m / 2)
#
#   #simulating arms' responses
#   yy = sapply(1:m, function(x) rbinom(1, N, p_true))
#   ee = NULL
#   for (ii in 1:m) {
#     ee = c(ee, c(rep(1, yy[ii]), rep(0, (N - yy[ii]))))
#   }
#
#   #randomizing the patients
#   pp = sample(1:(N*m), (N*m))
#
#   #arms indicator
#   arms = NULL
#   for (ii in 1:m) {
#     arms = c(arms, rep(ii, N))
#   }
#
#   #patients assignments
#   aa = arms[pp]
#   delta = matrix(0, nrow = m, ncol = (N*m))
#   for (ii in 1:(N*m)) {
#     delta[aa[ii], ii] = 1
#   }
#
#   #patients responses
#   responses = ee[pp]
#
#   #finding the first interim analysis times
#   first = max(sapply(1:m, function(y) {
#     which(cumsum(aa == y) == ia1)[1]
#   }))
#
#   #first interim analysis accruals
#   ia1_accruals = apply(delta[, 1:first], 1, sum)
#
#   #returning the object
#   return(list(assignments = delta, responses = responses, ia_accruals = as.data.frame(ia1_accruals)))
# }


#sim.data: a two-column dataframe with the number of positive responses in the first column and number of patients enrolled in each cohort in the second column; each row represents a cohort
#p_null: the null hypothesis uninteresting threshold
#p_target: the alternative hypothesis target threshold
#prior: prior distribution for the random term parameter handling the borrowing of information:
#   - inverse-Gamma (a, b) on the variance
#   - half-t(gamma, ni) on the standard deviation; when the number of groups is < 5
#   - uniform (a, b) on the standard deviation; when the number of subgroups is >= 5
#   - PC(sd_x) prior; for the EPC please use the EPC_sd function to compute the desired sd_x value

BHM_basket = function(sim.data, p_null, p_target, prior, parameters) {
  #loading INLA package
  # require(INLA)

  #loading logit and inverse logit functions
  # source("logit.R")
  # source("inv_logit.R")
  # source("log_uniform_precision.R")

  #checks on the dataframe
  if (class(sim.data) != "data.frame") {
    stop("Only data.frames allowed")
  }

  #checks on the dataframe dimensions
  if (ncol(sim.data) != 2) {
    stop("Wrong number of variables")
  }

  #checks on the dataframe values
  if (sum(sim.data[, 1] > sim.data[, 2]) != 0) {
    stop("The number of positive responses cannot be grater than the number of patients enrolled")
  }

  #storing prior's parameters
  prior_param = as.list(parameters)

  #defining the offset
  #DO WE NEED THIS???? in our paper this  value is equal to 0
  mu_par = logit(p_target) - logit(p_null)

  #check on the response rate null and target values
  if (p_null >= p_target) {
    stop("The null response rate cannot be greater or eaual than the target")
  }

  #PC prior
  if (prior == "PC") {
    #checks on the number of parameter
    if (length(prior_param) != 1) {
      stop("Wrong number of parameters for the PC prior")
    }

    #checks on the sign of the parameters
    if (prior_param <= 0) {
      stop("Negative values for parameters of the PC prior cannot be accepted")
    }

    #INLA code
    inla_output = INLA:::inla(y ~ 1 + f(x,
                                        model = "iid",
                                        hyper = list(prec =
                                                       list(initial = 0,
                                                            prior = "pc.prec",
                                                            param = c(prior_param[[1]] / 0.31, 0.01)))),
                              data = list(y = sim.data[, 1], x = 1:nrow(sim.data)),
                              family = "binomial",
                              Ntrials = sim.data[, 2], #histologies' sample size
                              offset = logit(p_target), #offset to be included (see BHMDesign.R)
                              control.fixed = list(mean.intercept = mu_par, prec.intercept = 1 / 100), #prior for the intercept
                              ### THINK about it! do we need correct=TRUE; it was giving an error when using it...
                              #control.inla = list(correct = TRUE),
                              control.compute=list(return.marginals.predictor=TRUE),
                              control.predictor = list(compute = TRUE),
                              verbose = FALSE)
  }

  #half-t prior
  else if (prior == "half-t") {
    #checks on the number of parameter
    if (length(prior_param) != 2) {
      stop("Wrong number of parameters for the Half-t prior")
    }

    #checks on the sign of the parameters
    if (prior_param[[1]] <= 0 | prior_param[[2]] <= 0) {
      stop("Negative values for parameters of the Half-t prior cannot be accepted")
    }

    #defining the half-t prior (from the Half-Cauchy in http://www.maths.bath.ac.uk/~jjf23/brinla/prior.html)
    log_half_t = paste("expression:
                   lambda = ", prior_param[[1]], " ;
                   ni = ", prior_param[[2]], " ;
                   precision = exp(log_precision);
                   logdens = lgamma((ni + 1) / 2) - lgamma(ni / 2) - 0.5 * log(ni * pi * lambda^2) - (ni + 1) / 2 * log(1 + 1 / (ni * precision * lambda^2)) - 1.5 * log_precision;
                   log_jacobian = log_precision;
                   return(logdens + log_jacobian);")
    half_t_prior = list(prec = list(prior = log_half_t))

    #INLA code
    inla_output = INLA:::inla(y ~ 1 + f(x,
                                        model = "iid",
                                        hyper = half_t_prior),
                              data = list(y = sim.data[, 1], x = 1:nrow(sim.data)),
                              family = "binomial",
                              Ntrials = sim.data[, 2], #histologies' sample size
                              offset = logit(p_target), #offset to be included (see BHMDesign.R)
                              control.fixed = list(mean.intercept = mu_par, prec.intercept = 1 / 100), #prior for the intercept
                              ### THINK about it! do we need correct=TRUE; it was giving an error when using it...
                              #control.inla = list(correct = TRUE),
                              control.compute=list(return.marginals.predictor=TRUE),
                              control.predictor = list(compute = TRUE),
                              verbose = FALSE)
  }

  #uniform prior
  else if (prior == "uniform") {
    #checks on the number of parameter
    if (length(prior_param) != 2) {
      stop("Wrong number of parameters for the Uniform prior")
    }

    #checks on the sign of the parameters
    if (prior_param[[1]] < 0 | prior_param[[2]] < 0) {
      stop("Negative values for parameters of the Uniform prior cannot be accepted")
    }

    #checks on the relation between the parameter
    if (prior_param[[1]] > prior_param[[2]]) {
      stop("Paramter a of the Uniform prior cannot be smaller than parameter b")
    }

    #tabling log-percision and log-uni pdf
    unif_prior_table = paste0("table: ",
                              paste(c(seq(- 10, 10, 0.01),
                                      log_unif_prec(seq(- 10, 10, 0.01),
                                                    prior_param[[1]], prior_param[[2]]))),
                              collapse = " ")

    #INLA code
    inla_output = INLA:::inla(y ~ 1 + f(x,
                                        model = "iid",
                                        hyper = list(prec = list(prior = unif_prior_table))),
                              data = list(y = sim.data[, 1], x = 1:nrow(sim.data)),
                              family = "binomial",
                              Ntrials = sim.data[, 2], #histologies' sample size
                              offset = logit(p_target), #offset to be included (see BHMDesign.R)
                              control.fixed = list(mean.intercept = mu_par, prec.intercept = 1 / 100), #prior for the intercept
                              ### THINK about it! do we need correct=TRUE; it was giving an error when using it...
                              #control.inla = list(correct = TRUE),
                              control.compute=list(return.marginals.predictor=TRUE),
                              control.predictor = list(compute = TRUE),
                              verbose = FALSE)
  }

  #inverse-gamma prior
  else if (prior == "inverse-gamma") {
    #checks on the number of parameter
    if (length(prior_param) != 2) {
      stop("Wrong number of parameters for the Inverse-Gamma prior")
    }

    #checks on the sign of the parameters
    if (prior_param[[1]] <= 0 | prior_param[[2]] <= 0) {
      stop("Negative values for parameters of the Inverse-Gamma prior cannot be accepted")
    }

    #INLA code
    inla_output = INLA:::inla(y ~ 1 + f(x,
                                        model = "iid",
                                        hyper = list(prec =
                                                       list(initial = 0,
                                                            prior = "loggamma",
                                                            param = c(prior_param[[1]],
                                                                      prior_param[[2]])))),
                              data = list(y = sim.data[, 1], x = 1:nrow(sim.data)),
                              family = "binomial",
                              Ntrials = sim.data[, 2], #histologies' sample size
                              offset = logit(p_target), #offset to be included (see BHMDesign.R)
                              control.fixed = list(mean.intercept = mu_par, prec.intercept = 1 / 100), #prior for the intercept
                              ### THINK about it! do we need correct=TRUE; it was giving an error when using it...
                              #control.inla = list(correct = TRUE),
                              control.compute=list(return.marginals.predictor=TRUE),
                              control.predictor = list(compute = TRUE),
                              verbose = FALSE)
  }

  #wrong specification of the prior
  else {print("Wrong prior specification")}

  #returning the results
  return(inla_output)
}
