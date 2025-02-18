WBAR = function(fit, data, config, accfunc = NULL) {
  
  # Set the correct objective function to optimize for
  if (config[["OBJECTIVE"]] == "Maximize") { objfunc = which.max }
  else { objfunc = which.min }
  
  num_trts = max(data["trt"])
  
  updated_fit = update(fit, newdata = data)
  
  posterior_draws = as_draws_df(updated_fit) |> 
    select(contains("b_Intercept"), contains("b_X")) |> 
    rename(b_X1 = b_Intercept) |> 
    mutate(across(-b_X1, ~ . + b_X1))
  
  posterior_means = posterior_draws |> 
    pivot_longer(everything(), names_to = "arm", values_to = "sample") |> 
    group_by(arm) |> 
    summarize(
      postmean = mean(sample)
    ) |> 
    pull(postmean)
  
  # Calculate weights for combining posterior samples with current best
  # First weight is for
  if (!is.null(accfunc)) {
    w = accfunc(data, config)
  } else {
    w = c(0.5, 0.5)
  }
  
  weighted_posterior_outcomes = posterior_draws |> 
    mutate(across(everything(), ~ (w[1] * .) + (w[2] * posterior_means[which(colnames(posterior_draws) == cur_column())])))
  
  # Calculate probability of being optimal from the posterior sample
  optimal_arms = apply(weighted_posterior_outcomes, 1, objfunc)
  probs = c(table(factor(optimal_arms, levels = 1:num_trts )) / nrow(weighted_posterior_outcomes))
  
  probs
  
}

