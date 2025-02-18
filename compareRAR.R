compareRAR = function(REP_ID, 
                      OBVS_FREQ, 
                      EFFECT_SIZE, 
                      AR_PARAM, 
                      CONFIG, 
                      ALG, 
                      NAME,
                      ACCFUNC = NULL) {
  
  # Start timer for simulation
  tic = lubridate::now()
  
  # Initialize output row with given simulation parameters
  OUT = tibble(
    REP_ID = REP_ID,
    OBVS_FREQ = OBVS_FREQ, 
    EFFECT_SIZE = EFFECT_SIZE, 
    AR_PARAM = AR_PARAM,
    CONFIG = CONFIG[["NAME"]],
    ALG = ALG,
    NAME = NAME
  )
  
  ##############################################################################
  # Initialization phase
  ##############################################################################
  
  # Design matrix of single person getting all of the treatments for n_obvs each
  X = diag(CONFIG[["NUM_TRTS"]])
  X[, 1] = 1
  X = X[rep(seq_len(nrow(X)), each = OBVS_FREQ),]
  
  # Generate outcome based on given AR model
  # For initialization, all 3 treatments are observed, so we must multiply to account for this
  Y = generateOutcomes(X, EFFECT_SIZE, CONFIG[["NUM_TRTS"]] * OBVS_FREQ, AR_PARAM, CONFIG)
  
  DATA = cbind(X, Y = Y) |>
    tibble::as_tibble() |>
    dplyr::mutate(
      period = rep(1:CONFIG[["NUM_TRTS"]], each = OBVS_FREQ),
      trt = rep(1:CONFIG[["NUM_TRTS"]], each = OBVS_FREQ)
    )
  
  colnames(DATA) = c(paste0("X", 1:CONFIG[["NUM_TRTS"]]), "Y", "period", "trt")
  
  # Calculate allocation probabilities via Thompson Sampling (or other algorithm)
  # stab_parameter = CONFIG[["NUM_TRTS"]] / (2 * CONFIG[["MAX_DURATION"]])
  
  FIT = update(CONFIG[["COMPILED_FIT"]], newdata = DATA)
  
  # Initialize tibbles to hold metrics for simulation
  PP_EFFECTIVE_TABLE = tibble()
  ESTIMATES_TABLE = tibble() # bias and credible interval of parameter
  
  ##############################################################################
  # Adaptive phase
  ##############################################################################
  
  adaptive_periods = (CONFIG[["NUM_TRTS"]] + 1):CONFIG[["MAX_DURATION"]]
  
  for (period in adaptive_periods) {
    
    ############################################################################
    # Generate data for an adaptive treatment assignment
    ############################################################################
    
    next_trt = selectNextTreatment(FIT, DATA, ALG, CONFIG, ACCFUNC)
    
    # Construct design matrix and "observe" the outcome under this treatment
    X = matrix(0, nrow = OBVS_FREQ, ncol = CONFIG[["NUM_TRTS"]])
    X[, 1] = 1        # intercept
    X[, next_trt] = 1 # setting next treatment
    Y = generateOutcomes(X, EFFECT_SIZE, OBVS_FREQ, AR_PARAM, CONFIG)
    
    period_cols = matrix(rep(period, OBVS_FREQ), nrow = OBVS_FREQ)
    assignment_cols = matrix(rep(next_trt, OBVS_FREQ), nrow = OBVS_FREQ)
    
    new_data = cbind(X, Y, period_cols, assignment_cols) |> as_tibble()
    colnames(new_data) = c(paste0("X", 1:CONFIG[["NUM_TRTS"]]), "Y", "period", "trt")
    
    # Append new data to current set
    DATA = rbind(DATA, new_data)
    
    ############################################################################
    # Calculate metrics of interest from the posterior distribution
    ############################################################################
    
    posterior_samples = getPosteriorSamples(FIT, DATA)
    
    # Calculate the posterior probability that a parameter > DELTA (since "Maximize")
    PP_EFFECTIVE_ROW = posterior_samples |> 
      select(contains("b_Intercept"), contains("b_X")) |> 
      rename(b_X1 = b_Intercept) |> 
      mutate(across(-b_X1, ~ . + b_X1)) |> 
      pivot_longer(everything(), names_to = "treatment", values_to = "sample") |> 
      group_by(treatment) |> 
      summarize(
        # prob_eff = mean(sample > CONFIG[["DELTA"]]), # For maximizing
        prob_eff = mean(sample < CONFIG[["DELTA"]])    # For minimizing
      ) |> 
      pivot_wider(values_from = "prob_eff", 
                  names_from = treatment) |> 
      mutate( period = period )
    
    # Calculate the posterior mean and the 95% credible interval for the regression parameters
    ESTIMATES_ROW = posterior_samples |> 
      select(contains("b_Intercept"), contains("b_X")) |> 
      rename(b_X1 = b_Intercept) |> 
      pivot_longer(everything(), names_to = "treatment", values_to = "sample") |> 
      group_by(treatment) |> 
      summarize(
        post_mean = mean(sample),
        q025 = quantile(sample, 0.025),
        q975 = quantile(sample, 0.975)
      ) |> 
      pivot_wider(values_from = c("post_mean", "q025", "q975"), 
                  names_from = treatment) |> 
      mutate( period = period )
    
    # Append results of simulation to the ongoing set of results
    PP_EFFECTIVE_TABLE = bind_rows(PP_EFFECTIVE_TABLE, PP_EFFECTIVE_ROW)
    ESTIMATES_TABLE = bind_rows(ESTIMATES_TABLE, ESTIMATES_ROW)
    
  } # end of adaptive phase
  
  # Convert the actual treatment assignment into wide
  ASSIGNMENTS_WIDE = DATA |> 
    transmute(
      col = glue::glue("assigned_{period}"),
      val = trt
    ) |> 
    distinct() |> 
    pivot_wider(names_from = col, values_from = val) 
  
  PP_EFFECTIVE_TABLE_WIDE = PP_EFFECTIVE_TABLE |> 
    pivot_wider(names_from = period, 
                values_from = contains("b_"))
  
  ESTIMATES_TABLE_WIDE = ESTIMATES_TABLE |> 
    pivot_wider(names_from = period, values_from = contains("b_"))
  
  # Combine everything into single output row
  OUT = bind_cols(OUT, ASSIGNMENTS_WIDE)
  OUT = bind_cols(OUT, PP_EFFECTIVE_TABLE_WIDE)
  OUT = bind_cols(OUT, ESTIMATES_TABLE_WIDE)
  
  # Finally, log the time it took to finish the 
  toc = lubridate::now()
  OUT[["DURATION"]] = toc - tic
  
  OUT
  
}