### da documentare; TO DO


#data: simulated trials
   # No, data dovrebbe essere l'output di sim_basket
#m: number of cohorts
#N: max number of patients enrolled per arm  (same for each arm)
#p_true: response rate for generating the data
#p_null:  the null hypothesis uninteresting threshold
#p_target: the alternative hypothesis target threshold
#prob_cutoff: the probability cutoff for controlling the type-I error rate

operating_char = function(data, m, N, p_true,
                          p_null, p_target, prob_cutoff) {

  #probability cutoff
  z = list(rep(prob_cutoff, m))

  #reorganizing the results
  results = lapply(seq_along(data), function(ii) {
    rr = sapply(1:m, function(jj) {
      #anti-logit transform
      post = data[[ii]]$post.par[[jj]]

      aa = sapply(1:length(z), function(ll) {
        #controlling for possible Dirac delta posterior
        ifelse('try-error' %in% class(try(INLA:::inla.pmarginal(p.null[jj], post), silent = TRUE)),
               #checking for a Dirac delta in o or 1
               ifelse(post[1, 2] > tail(post[, 2], n = 1), 0, 1),
               #computing the posterior probability
               ifelse((1 - INLA:::inla.pmarginal(p.null[jj], post)) >= z[[ll]][jj], 1, 0))
      })

      return(aa)
    })

    tmp = list(pm = sapply(data[[1]]$post.par, function(x) INLA:::inla.emarginal(function(y) y, x)),
               prob.resp = rr,
               sample.size = sum(data[[ii]]$accruals[, ncol(data[[ii]]$accruals)]),
               stop = data[[ii]]$stop,
               st.dev = c(INLA:::inla.emarginal(function(x) x, data[[ii]]$post.hyper), sqrt(INLA:::inla.emarginal(function(x) x^2, data[[ii]]$post.hyper) - INLA:::inla.emarginal(function(x) x, data[[ii]]$post.hyper)^2), INLA:::inla.qmarginal(p = c(0.025, 0.975), data[[ii]]$post.hyper)))

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
  prob.resp = round(apply(sapply(results, function(x) {
    x$prob.resp
  })[, issues], 1, mean), 3)
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
  })[, issues], 1, sum, na.rm = TRUE) / length(data), 3)
  names(fut.stops) = paste("arm", LETTERS[1:m])

  #efficacy stops
  eff.stops = round(apply(sapply(results, function(x) {
    t(sapply(1:m, function(y) {
      as.matrix(x$stop)[y, ncol(as.matrix(x$stop))] == 'eff'
    }))
  })[, issues], 1, sum, na.rm = TRUE) / length(data), 3)
  names(eff.stops) = paste("arm", LETTERS[1:m])

  #standard deviation
  post.hyper = round(apply(t(sapply(results, function(x) {
    x$st.dev
  })), 2, mean, na.rm = TRUE), 3)
  names(post.hyper) = c('mean', 'st.dev', 'CrI.LB', 'CrI.UB')

  #returning results
  return(list(`posterior means` = post.par, `rejection probabilities` = prob.resp, `expected sample size` = ess, `futility stops` = fut.stops, `efficacy stops` = eff.stops, `posterior hyperparameter` = post.hyper))
}

# # problema; funziona solo con m=4
find_cutoff <- function(res,
                        p_null,
                        prob_cutoff = 0.9,
                        desired_alpha = 0.1,
                        tol=0.05){
  n <- length(res)
  rej.prob1 <- rej.prob2 <- rej.prob3 <- rej.prob4  <- rep(1, n)
  while(sum(abs(c(mean(rej.prob1), mean(rej.prob2), mean(rej.prob3), mean(rej.prob4)) - desired_alpha) <= tol) < 4){
    prob_cutoff <- prob_cutoff-0.001
    for(i in 1:n) {
      rej.prob1[i] <- (1-inla.pmarginal(p_null, res[[i]]$post.par[[1]])) >= prob_cutoff
      rej.prob2[i] <- (1-inla.pmarginal(p_null, res[[i]]$post.par[[2]])) >= prob_cutoff
      rej.prob3[i] <- (1-inla.pmarginal(p_null, res[[i]]$post.par[[3]])) >= prob_cutoff
      rej.prob4[i] <- (1-inla.pmarginal(p_null, res[[i]]$post.par[[4]])) >= prob_cutoff
    }
    print(c(prob_cutoff, mean(rej.prob1), mean(rej.prob2), mean(rej.prob3), mean(rej.prob4)))
  }

}

