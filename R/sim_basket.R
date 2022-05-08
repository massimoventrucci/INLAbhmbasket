#sim: number of simulated trials
#m: number of cohorts
#N: max number of patients enrolled per arm  (same for each arm)
#p_true: response rate for generating the data
#p_null:  the null hypothesis uninteresting threshold
#  PROBLEMA NEL CODICE!!!
#   bisogna fare in modo di permettere all'utente di scegliere  gli scenari,
#   occorre un input con il numero di coorti in cui il trattamento è efficace:
#    0, 1, 2, ... , m

#p_target: the alternative hypothesis target threshold
#ia1_fraction: define 1st IA as fraction of N, as numeric
#step: define subsequent IAs' timing as an increment of the patients enrolled at the 1st IA
#utility_threshold = the futility thereshold
#efficacy_threshold = the futility thereshold
#prior: prior distribution for the random term parameter handling the borrowing of information:
#   - inverse-Gamma (a, b) on the variance
#   - half-t(gamma, ni) on the standard deviation; when the number of groups is < 5
#   - uniform (a, b) on the standard deviation; when the number of subgroups is >= 5
#   - PC(sd_x) prior; for the EPC please use the EPC_sd function to compute the desired sd_x value
#parameters: prior distribution's parameters
### da documentare; TO DO

sim_basket = function(sim, m, N, p_true,
                      p_null, p_target,
                      ia1_fraction = 0.4, step = 0.5,
                      futility_threshold = 0.05,
                      efficacy_threshold = 0.90,
                      prior, parameters) {
  #require(pbapply)

  #options for the progress bar
  op = pbapply:::pboptions(type = "timer")

  tmp = pbapply:::pblapply(1:sim, function(xx) {
    #require(INLA)
    # source('BHM_basket.R')
    # source('EPC.R')
    # source('inv_logit.R')
    # source('log_uniform_precision.R')
    # source('logit.R')
    # source('posterior_summary.R')
    # source('trial.R')

    #simulating arrivals
    patients = trial(m, N, p_true, ia1_fraction)

    #accruals up the i-th interim analysis
    NN = patients$ia_accruals

    #defining the data for each interim analysis
    yy = sapply(1:m, function(x) {
      sum(patients$responses[1:sum(NN)] * patients$assignments[x, 1:sum(NN)])
    })

    #BHM
    inla.mod = BHM_basket(cbind.data.frame(yy, NN),
                          p_null, p_target, prior, parameters)

    #posterior parameters
    post = inla.mod$marginals.fitted.values

    #standard deviation posterior
    st.dev = INLA:::inla.tmarginal(function(x) {1 / sqrt(exp(x))}, inla.mod$internal.marginals.hyperpar[[1]])

    #posterior probability
    post.prob = 1 - posterior_probability(data = inla.mod, x = mean(p_null, p_target))

    #first interim analyses
    stop = as.matrix(sapply(post.prob, function(x) {
      ifelse(x < futility_threshold, 'fut', ifelse(x > efficacy_threshold, 'eff', 'no'))
    }))

    #adjusting for the NAs in stop
    for (kk in 1:m) {
      if(is.na(stop[kk, 1])) {
        if(INLA:::inla.smarginal(post[[kk]])$x < 1 & INLA:::inla.smarginal(post[[kk]])$y > 10^20) {
          stop[kk, 1] = 'fut'
        }
      }
    }

    #naming columns
    colnames(stop) = 'ia1'

    #modifying the assignments
    for (x in 1:m) {
      if (stop[x, 1] != 'no') {
        patients$assignments[x, ] = c(patients$assignments[x, 1:sum(patients$ia_accruals)], rep(0, sum(m*N) - sum(patients$ia_accruals)))
      }
    }

    #storing responses
    successes = as.matrix(yy)

    #initializing the check on the maximum sample size
    accr = rep(TRUE, m)

    #initializing counter
    jj = 2

    #successive interim analyses
    ias_step = ceiling(sum(patients$ia_accruals) * step)

    while (sum(stop[, (jj - 1)] == 'no') >= 1 & sum(accr) != 0) {
      #updating accruals
      index = which(cumsum(apply(patients$assignments, 2, sum)) == sum(patients$ia_accruals[, jj - 1]) + ias_step)[1]
      index = ifelse(is.na(index), sum(m*N), index)
      NN = apply(patients$assignments[, 1:index], 1, sum)
      patients$ia_accruals = cbind(patients$ia_accruals, NN)

      #updating responses
      yy = sapply(1:m, function(x) {
        sum(patients$responses[1:index] * patients$assignments[x, 1:index])
      })

      #BHM
      inla.mod = BHM_basket(cbind.data.frame(yy, NN), p_null, p_target, prior, parameters)

      #posterior parameters
      post = inla.mod$marginals.fitted.values

      #standard deviation posterior
      st.dev = INLA:::inla.tmarginal(function(x) {1 / sqrt(exp(x))}, inla.mod$internal.marginals.hyperpar[[1]])

      #interim analyses
      stop = cbind(stop, as.matrix(sapply(post.prob, function(x) {
        ifelse(x < futility_threshold, 'fut', ifelse(x > efficacy_threshold, 'eff', 'no'))
      })))

      #adjusting for the NAs in stop
      for (kk in 1:m) {
        if(is.na(stop[kk, 1])) {
          if(INLA:::inla.smarginal(post[[kk]])$x < 1 & INLA:::inla.smarginal(post[[kk]])$y > 10^20) {
            stop[kk, 1] = 'fut'
          }
        }
      }

      #naming stop columns
      colnames(stop) = paste0('ia', 1:jj)
      rownames(stop) = paste("arm", LETTERS[1:m])

      #modifying the assignments
      for (x in 1:m) {
        if (stop[x, 1] != 'no') {
          patients$assignments[x, ] = c(patients$assignments[x, 1:index], rep(0, sum(m*N) - index))
        }
      }

      #check on the maximum sample size
      accr = sapply(1:m, function(x) {
        stop[x, jj] == 'no' & patients$ia_accruals[x, jj] < N
      })

      #naming patients accruals columns
      colnames(patients$ia_accruals) = c(paste0('ia', 1:(jj - 1)), 'fa')
      rownames(patients$ia_accruals) = paste("arm", LETTERS[1:m])

      #updating the counter
      jj = jj + 1
    }

    stop = as.matrix(stop)
    colnames(stop) = paste0('ia', 1:ncol(stop))
    rownames(stop) = paste("arm", LETTERS[1:m])

    patients$ia_accruals = as.matrix(patients$ia_accruals)
    colnames(patients$ia_accruals) = paste0('ia', 1:ncol(patients$ia_accruals))
    rownames(patients$ia_accruals) = paste("arm", LETTERS[1:m])

    #returing data
    return(list(stop = as.matrix(stop[, 1:(jj - 2)]), post.par = post, post.hyper = st.dev, accruals = as.matrix(patients$ia_accruals[, 1:(jj - 1)]), successes = successes))
  })

  #"close" the options
  pbapply:::pboptions(op)

  #returning the simulated data
  return(tmp)
}
