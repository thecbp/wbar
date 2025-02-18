###############################################################################

# Helper functions for running simulations

###############################################################################

precompileBrmsModel = function(NUM_TRTS) {
  # Precompile Stan model for simulation repetitions
  
  # Build design matrix for
  X = diag(NUM_TRTS) 
  X[, 1] = 1
  X = X[rep(seq_len(nrow(X)), each = 1),]
  X = X[rep(seq_len(nrow(X)), times = 1),]
  Y = rnorm(nrow(X))
  DATA = cbind(X, Y = Y) |> tibble::as_tibble()
  colnames(DATA) = c(paste0("X", 1:NUM_TRTS), "Y", "period", "trt")
  
  brms::brm(formula = paste0("Y~", paste0("X", 2:NUM_TRTS, collapse = "+")),
            data = DATA,
            family = gaussian(link = "identity"),
            prior = c(
              brms::set_prior("normal(0,100)", class = "Intercept"),
              brms::set_prior("normal(0,100)", class = "b")
            ),
            refresh = 0)
  
}

precompileAR1BrmsModel = function(NUM_TRTS) {
  # Precompile Stan model for simulation repetitions
  
  # Build design matrix for
  X = diag(NUM_TRTS) 
  X[, 1] = 1
  X = X[rep(seq_len(nrow(X)), each = 1),]
  X = X[rep(seq_len(nrow(X)), times = 1),]
  Y = rnorm(nrow(X))
  DATA = cbind(X, Y = Y) |> 
    tibble::as_tibble() |> 
    mutate(
      obs = row_number()
    )
  colnames(DATA) = c(paste0("X", 1:NUM_TRTS), "Y", "obs")
  
  brms::brm(formula = paste0("Y ~ ",                             # outcome
                             paste0("X", 2:3, collapse = " + "), # treatments
                             " + ar(time = obs, p = 1)"),        # AR1 term
            data = DATA,
            family = gaussian(link = "identity"),
            prior = c(
              brms::set_prior("normal(0,100)", class = "Intercept"),
              brms::set_prior("normal(0,100)", class = "b")
            ),
            refresh = 0)
  
}

generateOutcomes = function(X, EFFECT_SIZE, OBVS_FREQ, AR_PARAM, CONFIG) {
  
  # Initialize regression coefficients
  BETAS = rep(0, CONFIG[["NUM_TRTS"]])       # Initialize vector of null effects
  BETAS[1] = CONFIG[["INTERCEPT"]]
  BETAS[2] = EFFECT_SIZE * CONFIG[["SIGMA"]] # Set second treatment as effective treatment
  
  if (AR_PARAM == 0.0) {
    
    Y = ((X %*% BETAS) |> c()) + stats::rnorm(nrow(X), 0, sd = CONFIG[["SIGMA"]])
    
  } else if (AR_PARAM > 0.0) {
    
    Y_ar = arima.sim(n = OBVS_FREQ,
                     model = list(ar = c(AR_PARAM)), 
                     sd = CONFIG[["SIGMA"]])
    Y = Y_ar + (X %*% BETAS)
    
  }
  
  Y
  
}

selectNextTreatment = function(FIT, DATA, ALG, CONFIG, ACCFUNC) {
  
  if (ALG == "BAR") {
    
    # My original probability matching implementation of Thompson Sampling
    alloc_probs = BAR(FIT, DATA, CONFIG[["OBJECTIVE"]], STAB)
    sample(1:CONFIG[["NUM_TRTS"]], size = 1, prob = alloc_probs)
  
  } else if (ALG == "CB") {
    
    # Myopic current best algorithm
    FIT = update(CONFIG[["COMPILED_FIT"]], newdata = DATA)
    
    posterior_means = as_draws_df(FIT) |> 
      select(contains("b_Intercept"), contains("b_X")) |> 
      rename(b_X1 = b_Intercept) |> 
      mutate(across(-b_X1, ~ . + b_X1)) |> 
      pivot_longer(everything(), names_to = "arm", values_to = "sample") |> 
      group_by(arm) |> 
      summarize(
        postmean = mean(sample)
      ) |> 
      pull(postmean)
    
    if (CONFIG[["OBJECTIVE"]] == "Maximize") { objfunc = which.max }
    else { objfunc = which.min }
    
    next_trt = objfunc(posterior_means)
    
  } else if (ALG == "WBAR") {
    
    # Proposed algorithm 1: maximizing weighted average of TS and CB
    alloc_probs = WBAR(FIT, DATA, CONFIG, ACCFUNC)
    sample(1:CONFIG[["NUM_TRTS"]], size = 1, prob = alloc_probs)
    
    
  } 
  
}

getPosteriorSamples = function(FIT, DATA) {
  
  newfit = update(FIT, newdata = DATA)
  samples = posterior::as_draws_df(newfit)
  
  samples
}
