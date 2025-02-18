################################################################################

# Functions for defining different weighitng functions for WBAR

################################################################################

tuneExponential = function(r) {
  
  weightingExponential = function(data, config, rho = r) {
    
    i = max(data |> pull(period))
    n = config[["MAX_DURATION"]]
    
    # Exponential acceleration function
    w2 = (i / n)^rho
    w1 = 1 - w2
    
    c(w1, w2)
    
  }
  
  return(weightingExponential)
  
}

tuneLogistic = function(b, mid) {
  
  weightingLogistic = function(data, config, B = b, x0 = mid) {
    
    i = max(data["period"])
    n = config[["MAX_DURATION"]]
    
    # Exponential acceleration function
    w2 = (1 + exp(-B * ((i / n) - x0)))^(-1)
    w1 = 1 - w2
    
    c(w1, w2)
    
  }
  
  return(weightingLogistic)
  
}