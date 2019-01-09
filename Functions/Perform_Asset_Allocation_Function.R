#'----
#' author: Gabriel Cabrera 
#' title: Asset Allocation Function
#' date: 28/12/2018 (last updated: 28/12/2018)
#'---

Perform_asset_allocation <- function(actual, risk_free, forecast, 
                                     volatility_forecast, gamma_MV){
  
  
  T <- matrix(0, nrow = NROW(actual), ncol = 1)
  weight_risky <<- matrix(0, nrow = NROW(T), ncol = 1)
  
  for(iter in seq_along(T)){
    
    weight_risky[iter] <- (1/gamma_MV)*(forecast[iter]/volatility_forecast[iter])

    ifelse(weight_risky[iter,] < 0,
           weight_risky[iter,] <- 0,
           ifelse(weight_risky[iter,] > 1.5,
                  weight_risky[iter,] <- 1.5,
                  weight_risky[iter,]))
  
  }
  
  return_portfolio <- risk_free + weight_risky*actual
  avg_utility <<- mean(return_portfolio) - 0.5*gamma_MV*(sd(return_portfolio))^2
  
}

