### da documentare; TO DO


#sim.trials: simulated trials
#m: number of cohorts
#N: max number of patients enrolled per arm  (same for each arm)
#p_true: response rate for generating the data;
#   NOTE: I do not think p_true is needed for computing operational char!!!!!
#p_null:  the null hypothesis uninteresting threshold
#p_target: the alternative hypothesis target threshold
#   NOTE: I do not think p_target is needed for computing operational char!!!!!
#prob_cutoff: the probability cutoff for controlling the type-I error rate

operating_char = function(sim.trials,
                          m,
                          #N,
                          #p_true,
                          p_null,
                          #p_target,
                          prob_cutoff) {

  #probability cutoff
  #z = prob_cutoff #list(rep(prob_cutoff, m))

  #reorganizing the results
  results = lapply(seq_along(sim.trials), function(ii) {
    rr = sapply(1:m, function(jj) {
      #anti-logit transform
      post = sim.trials[[ii]]$post.par[[jj]]

      aa <- ifelse('try-error' %in% class(try(INLA:::inla.pmarginal(p_null, post), silent = TRUE)),
                   #checking for a Dirac delta in o or 1
                   ifelse(post[1, 2] > tail(post[, 2], n = 1), 0, 1),
                   #computing the posterior probability
                   ifelse((1 - INLA:::inla.pmarginal(p_null, post)) >= prob_cutoff, 1, 0))

      # old; by Ale
      # aa = sapply(1:length(z), function(ll) {
      #   #controlling for possible Dirac delta posterior
      #   ifelse('try-error' %in% class(try(INLA:::inla.pmarginal(p_null[jj], post), silent = TRUE)),
      #          #checking for a Dirac delta in o or 1
      #          ifelse(post[1, 2] > tail(post[, 2], n = 1), 0, 1),
      #          #computing the posterior probability
      #          ifelse((1 - INLA:::inla.pmarginal(p_null[jj], post)) >= z[[ll]][jj], 1, 0))
      # })

      return(aa)
    })

    tmp = list(
      pm = sapply(sim.trials[[ii]]$post.par, function(x) INLA:::inla.emarginal(function(y) y, x)),
      #pm = sapply(sim.trials[[1]]$post.par, function(x) INLA:::inla.emarginal(function(y) y, x)),
      prob.resp = rr,
      sample.size = sum(sim.trials[[ii]]$accruals[, ncol(sim.trials[[ii]]$accruals)]),
      stop = sim.trials[[ii]]$stop,
      st.dev = c(INLA:::inla.emarginal(function(x) x, sim.trials[[ii]]$post.hyper), sqrt(INLA:::inla.emarginal(function(x) x^2, sim.trials[[ii]]$post.hyper) - INLA:::inla.emarginal(function(x) x, sim.trials[[ii]]$post.hyper)^2), INLA:::inla.qmarginal(p = c(0.025, 0.975), sim.trials[[ii]]$post.hyper)))

    #returning results
    return(tmp)
  })

  #removing "problematic" trials
  issues = sapply(results, function(x) {
    !is.na(sum(x$prob.resp))
  })

  #posterior means
  post.par = round(apply(sapply(results, function(x) {
    x$pm
  })[, issues], 1, mean), 3)
  names(post.par) = paste("arm", LETTERS[1:m])

  #probability of response
  # prob.resp = round(apply(sapply(results, function(x) {
  #   x$prob.resp
  # })[, issues], 1, mean), 3)
  prob.resp = apply(sapply(results, function(x) {
    x$prob.resp
  })[, issues], 1, mean)
  names(prob.resp) = paste("arm", LETTERS[1:m])

  #expected sample size
  ess = as.numeric(round(mean(sapply(results, function(x) {
    x$sample.size
  })[issues]), 3))

  #futility stops
  fut.stops = round(apply(sapply(results, function(x) {
    t(sapply(1:m, function(y) {
      as.matrix(x$stop)[y, ncol(as.matrix(x$stop))] == 'fut'
    }))
  })[, issues], 1, sum, na.rm = TRUE) / length(sim.trials), 3)
  names(fut.stops) = paste("arm", LETTERS[1:m])

  #efficacy stops
  eff.stops = round(apply(sapply(results, function(x) {
    t(sapply(1:m, function(y) {
      as.matrix(x$stop)[y, ncol(as.matrix(x$stop))] == 'eff'
    }))
  })[, issues], 1, sum, na.rm = TRUE) / length(sim.trials), 3)
  names(eff.stops) = paste("arm", LETTERS[1:m])

  #standard deviation
  post.hyper = round(apply(t(sapply(results, function(x) {
    x$st.dev
  })), 2, mean, na.rm = TRUE), 3)
  names(post.hyper) = c('mean', 'st.dev', 'CrI.LB', 'CrI.UB')

  #returning results
  return(list(`posterior means` = post.par,
              `rejection probabilities` = prob.resp,
              `expected sample size` = ess,
              `futility stops` = fut.stops,
              `efficacy stops` = eff.stops,
              `posterior hyperparameter` = post.hyper,
              `issues` = issues))
}


### function to find the cut off;
# maybe not optimal but works fairly well
# for future; do another version where we search more accurately and make sure
# that all arms reach the desired alpha.target
find_cutoff <- function(alpha_target,
                        cutoff.range.ini = c(0.7, 0.999),
                        sim.trials,
                        m,
                        p_null, ...){
  cutoff <- optimize(f=INLAbhmbasket:::objective,
                     interval = cutoff.range.ini,
                     maximum = FALSE,
                     alpha_target = alpha_target,
                     sim.trials = sim.trials,
                     m=m,
                     p_null=p_null, ...)
  return(cutoff$minimum)
}
objective <- function(x, alpha_target, sim.trials, m, p_null){
  op <- operating_char(sim.trials = sim.trials,
                       m=m,
                       p_null=p_null,
                       prob_cutoff= x)
  dev <- (op$"rejection probabilities" - alpha_target)^2
  return(sum(dev))
}




