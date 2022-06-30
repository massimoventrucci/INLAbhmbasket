### NON mi è chiaro se va documentata anche questa,
### o se è una funzione interna

#summarizing the posterior distribution

#posterior summary for the arms
arms_posterior_summary = function(sim.trials,
                                  quantiles = FALSE,
                                  median = FALSE,
                                  mode = FALSE,
                                  HPD.CrI = FALSE,
                                  credibility.level = 0.95) {
  #basic summaries
  tmp = sim.trials$summary.fitted.values[, 1:2]

  #adding quantiles
  if (quantiles) {
    tmp = cbind.data.frame(tmp, sim.trials$summary.fitted.values[, c(3, 5)])
  }

  #adding median
  if (median) {
    tmp = cbind.data.frame(tmp, median = sim.trials$summary.fitted.values[, 4])
  }

  #adding mode
  if (mode) {
    tmp = cbind.data.frame(tmp, mode = sim.trials$summary.fitted.values[, 6])
  }

  #adding HPD CrI
  if (HPD.CrI) {
    tmp = cbind.data.frame(tmp, t(sapply(sim.trials$marginals.fitted.values, function(x) INLA:::inla.hpdmarginal(credibility.level, x))))

    #modifying the names
    names(tmp)[c(length(tmp) - 1, length(tmp))] = c("HPD.CrI.LB", "HPD.CrI.UB")
  }

  #modifing rownames
  rownames(tmp) = paste("arm", LETTERS[1:nrow(tmp)])

  #returning results
  return(tmp)
}

#posterior summary mixing parameter
mixing_posterior_summary = function(sim.trials,
                                    scale = "precision",
                                    quantiles = FALSE,
                                    median = FALSE,
                                    mode = FALSE,
                                    HPD.CrI = FALSE,
                                    credibility.level = 0.95) {
  if (scale == "precision") {
    #basic summaries
    tmp = sim.trials$summary.hyperpar[, 1:2]

    #adding quantiles
    if (quantiles) {
      tmp = cbind.data.frame(tmp, sim.trials$summary.hyperpar[, c(3, 5)])
    }

    #adding median
    if (median) {
      tmp = cbind.data.frame(tmp, median = sim.trials$summary.hyperpar[, 4])
    }

    #adding mode
    if (mode) {
      tmp = cbind.data.frame(tmp, mode = sim.trials$summary.hyperpar[, 6])
    }

    #adding HPD CrI
    if (HPD.CrI) {
      tmp = cbind.data.frame(tmp, t(sapply(sim.trials$marginals.hyperpar, function(x) INLA:::inla.hpdmarginal(0.95, x))))

      #modifying the names
      names(tmp)[c(length(tmp) - 1, length(tmp))] = c("HPD.CrI.LB", "HPD.CrI.UB")
    }

    #modigying the names
    rownames(tmp) = "precision"
  } else if (scale == "sd") {
    #basic summaries
    tmp = cbind.data.frame(mean = INLA:::inla.emarginal(function(x) {1 / sqrt(x)}, sim.trials$marginals.hyperpar[[1]]),
                           sd = sqrt(INLA:::inla.emarginal(function(x) 1 / x,
                                                           sim.trials$marginals.hyperpar[[1]]) -
                                       INLA:::inla.emarginal(function(x) 1 / sqrt(x),
                                                             sim.trials$marginals.hyperpar[[1]])^2))

    #adding quantiles
    if (quantiles) {
      tmp = cbind.data.frame(tmp, `0.025quant` = INLA:::inla.qmarginal(0.025, INLA:::inla.tmarginal(function(x) {1 / sqrt(x)},
                                                                                                    sim.trials$marginals.hyperpar[[1]])),
                             `0.975quant` = INLA:::inla.qmarginal(0.975, INLA:::inla.tmarginal(function(x) {1 / sqrt(x)},
                                                                                 sim.trials$marginals.hyperpar[[1]])))
    }

    #adding median
    if (median) {
      tmp = cbind.data.frame(tmp,
                             median = INLA:::inla.qmarginal(0.5,
                                                            INLA:::inla.tmarginal(function(x) {1 / sqrt(x)},
                                                                                  sim.trials$marginals.hyperpar[[1]])))
    }

    #adding mode
    if (mode) {
      tmp = cbind.data.frame(tmp, mode = INLA:::inla.mmarginal(INLA:::inla.tmarginal(function(x) {1 / sqrt(x)}, sim.trials$marginals.hyperpar[[1]])))
    }

    #adding HPD CrI
    if (HPD.CrI) {
      tmp = cbind.data.frame(tmp, INLA:::inla.hpdmarginal(credibility.level, INLA:::inla.tmarginal(function(x) {1 / sqrt(x)}, sim.trials$marginals.hyperpar[[1]])))

      #modifying the names
      names(tmp)[c(length(tmp) - 1, length(tmp))] = c("HPD.CrI.LB", "HPD.CrI.UB")
    }

    #modigying the names
    rownames(tmp) = "sd"
  }  else if (scale == "variance") {
    #basic summaries
    tmp = cbind.data.frame(mean = INLA:::inla.emarginal(function(x) {1 / x}, sim.trials$marginals.hyperpar[[1]]),
                           sd = sqrt(INLA:::inla.emarginal(function(x) 1 / x^2, sim.trials$marginals.hyperpar[[1]]) - INLA:::inla.emarginal(function(x) 1 / x, sim.trials$marginals.hyperpar[[1]])^2))

    #adding quantiles
    if (quantiles) {
      tmp = cbind.data.frame(tmp, `0.025quant` = INLA:::inla.qmarginal(0.025, INLA:::inla.tmarginal(function(x) {1 / x}, sim.trials$marginals.hyperpar[[1]])), `0.975quant` = INLA:::inla.qmarginal(0.975, INLA:::inla.tmarginal(function(x) {1 / x}, sim.trials$marginals.hyperpar[[1]])))
    }

    #adding median
    if (median) {
      tmp = cbind.data.frame(tmp, median = INLA:::inla.qmarginal(0.5, INLA:::inla.tmarginal(function(x) {1 / x}, sim.trials$marginals.hyperpar[[1]])))
    }

    #adding mode
    if (mode) {
      tmp = cbind.data.frame(tmp, mode = INLA:::inla.mmarginal(INLA:::inla.tmarginal(function(x) {1 / x}, sim.trials$marginals.hyperpar[[1]])))
    }

    #adding HPD CrI
    if (HPD.CrI) {
      tmp = cbind.data.frame(tmp, INLA:::inla.hpdmarginal(credibility.level, INLA:::inla.tmarginal(function(x) {1 / x}, sim.trials$marginals.hyperpar[[1]])))

      #modifying the names
      names(tmp)[c(length(tmp) - 1, length(tmp))] = c("HPD.CrI.LB", "HPD.CrI.UB")
    }

    #modigying the names
    rownames(tmp) = "variance"
  }

  #returning results
  return(tmp)
}

#posterior probability (CDF)
posterior_probability = function(inla.out, x) {
  sapply(inla.out$marginals.fitted.values, function(y) {
    #controlling for possible Dirac delta posterior
    if ('try-error' %in% class(try(INLA:::inla.pmarginal(x, y), silent = TRUE))) {
      #checking for a Dirac delta in 0 or 1
      prob = ifelse(y[1, 2] > tail(y[, 2], n = 1), 1, 0)
    } else {
      #computing the posterior probability
      prob = INLA:::inla.pmarginal(x, y)
    }
  })
}





