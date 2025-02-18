BAR = function(fit, data, objective, stab = NULL) {
  
  # Set the correct objective function to optimize for
  if (objective == "Maximize") { objfunc = which.max }
  else { objfunc = which.min }
  
  num_trts = max(data["trt"])
  
  updated_fit = update(fit, newdata = data)
  
  posterior_outcomes = as_draws_df(updated_fit) |> 
    select(contains("b_Intercept"), contains("b_X")) |> 
    rename(b_X1 = b_Intercept) |> 
    mutate(across(-b_X1, ~ . + b_X1)) 
  
  # Calculate probability of being optimal from the posterior sample
  optimal_arms = apply(posterior_outcomes, 1, objfunc)
  probs = c(table(factor(optimal_arms, levels = 1:num_trts )) / nrow(posterior_outcomes))
  
  if (!is.null(stab)) {
    probs = probs^stab / sum(probs^stab)
  }
  
  probs
  
}

