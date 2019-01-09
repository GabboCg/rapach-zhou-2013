#'----
#' author: Gabriel Cabrera 
#' title:
#' date: 17/12/2018 (last updated: 17/12/2018)
#'---

if(!require("pacman")) iinstall.packages("pacman")
p_load("tidyverse", "microbenchmark", "janitor", "lubridate", "zoo", "sandwich",
       "pracma")

handbook_data <- readxl:: read_xls("data/Returns_handbook_data.xls", sheet = 1)

data <- handbook_data  %>% 
        clean_names() %>% 
        filter(date_yyyymm >= "192611" & date_yyyymm <= "201012")

data <- as.data.frame(apply(data, 2, as.numeric)) %>% 
        rename(market_return = crsp_s_p_500_value_weighted_return_with_dividends, 
               D12 = x12_month_moving_sum_of_s_p_500_dividends,
               SP500 = s_p_500_index,
               E12 = x12_month_moving_sum_of_s_p_500_earnings,
               SVAR = monthly_sum_of_squared_daily_returns_on_s_p_500_index,
               BM = djia_book_to_market_value_ratio,
               NTIS = net_equity_expansion,
               TBL = x3_month_treasury_bill_yield_secondary_market,
               LTY = long_term_government_bond_yield,
               LTR = long_term_government_bond_return,
               AAA = moodys_aaa_rated_corporate_bond_yield,
               BAA = moodys_baa_rated_corporate_bond_yield,
               CORPR = long_term_corporate_bond_return,
               INFL = cpi_all_urban_consumers_inflation_rate) %>% 
        mutate(equity_premium = log(1 + market_return) - log(1 + dplyr::lag(risk_free_rate)),
               DP = log(D12) - log(SP500),
               EP = log(E12) - log(SP500),
               TMS = LTY - TBL,
               DFY = BAA -AAA,
               DFR = CORPR - LTR,
               INFL_lag = dplyr::lag(INFL),
               DY = log(D12) - log(dplyr::lag(SP500)),
               DE = log(D12) - log(E12)
        ) %>% 
        select(-date_yyyymm, -crsp_s_p_500_value_weighted_return_excluding_dividends,
               -nber_recession_dummies,-nber_recession_dummies_with_peak_included,
               -SP500, -INFL, -D12, -E12, -risk_free_rate, -BAA, -AAA) %>% 
        select(equity_premium, DP, DY, EP, DE, SVAR, BM, NTIS, TBL, LTY, LTR, TMS, DFY, DFR, INFL_lag) %>% 
        na.omit()


ECON_Sink <- handbook_data  %>% 
             clean_names() %>% 
             filter(date_yyyymm >= "192611" & date_yyyymm <= "201012")
   
ECON_Sink <- as.data.frame(apply(ECON_Sink, 2, as.numeric)) %>% 
             rename(market_return = crsp_s_p_500_value_weighted_return_with_dividends, 
                    D12 = x12_month_moving_sum_of_s_p_500_dividends,
                    SP500 = s_p_500_index,
                    E12 = x12_month_moving_sum_of_s_p_500_earnings,
                    SVAR = monthly_sum_of_squared_daily_returns_on_s_p_500_index,
                    BM = djia_book_to_market_value_ratio,
                    NTIS = net_equity_expansion,
                    TBL = x3_month_treasury_bill_yield_secondary_market,
                    LTY = long_term_government_bond_yield,
                    LTR = long_term_government_bond_return,
                    AAA = moodys_aaa_rated_corporate_bond_yield,
                    BAA = moodys_baa_rated_corporate_bond_yield,
                    CORPR = long_term_corporate_bond_return,
                    INFL = cpi_all_urban_consumers_inflation_rate) %>% 
              mutate(equity_premium = market_return - dplyr::lag(risk_free_rate),
                     DP = log(D12) - log(SP500),
                     EP = log(E12) - log(SP500),
                     TMS = LTY - TBL,
                     DFY = BAA -AAA,
                     DFR = CORPR - LTR,
                     INFL_lag = dplyr::lag(INFL),
                     DY = log(D12) - log(dplyr::lag(SP500)),
                     DE = log(D12) - log(E12)
              ) %>% 
              select(-date_yyyymm, -crsp_s_p_500_value_weighted_return_excluding_dividends,
                     -nber_recession_dummies,-nber_recession_dummies_with_peak_included,
                     -SP500, -INFL, -D12, -E12, -risk_free_rate, -BAA, -AAA) %>% 
             select(equity_premium, DP, DY, EP, SVAR, BM, NTIS, TBL, LTY, LTR, DFY, DFR, INFL_lag) %>% 
             na.omit()


# sum-of-the-parts variables ---------------------------------------------------
SopData <- handbook_data  %>% 
           clean_names() %>% 
           filter(date_yyyymm >= "192611" & date_yyyymm <= "201012") 

SopData <- as.data.frame(apply(SopData, 2, as.numeric)) %>% 
           rename(SP500 = s_p_500_index,
                  D12 = x12_month_moving_sum_of_s_p_500_dividends,
                  E12 = x12_month_moving_sum_of_s_p_500_earnings) %>% 
           mutate(E12_lag = lag(E12,1),
                  E_growth = (E12 - E12_lag)/E12_lag,
                  DP_SOP = (1/12)*D12/SP500,
                  r_f = risk_free_rate
           ) %>% 
           select(E12_lag, E_growth, DP_SOP, r_f) %>% 
           na.omit()

T <- NROW(data)
LHS <- data$equity_premium[2:T]

ols <- function(Y,X){
  
  betas <<- solve(t(X)%*%X)%*%(t(X)%*%Y) 
  
  e <<- Y - X%*%betas
  variance <- t(e)%*%e
  vcv <- 1/(NROW(X)-NCOL(X))*(as.numeric(variance)*solve(t(X)%*%X))
  betastd <<- sqrt(diag(vcv))
  
  regr <<- cbind(betas, as.matrix(betastd))

}

Y <- data$equity_premium

data <- data %>% select(-equity_premium)
ECON_Sink <- ECON_Sink %>% select(-equity_premium)

beta_full <- matrix(0, nrow = NCOL(data), ncol = 1)

for(i in 1:NCOL(data)){
  
  RHS_i <- cbind(as.matrix(data[1:(T-1),i]), matrix(1,nrow = (T-1), ncol = 1))
  results_i <- ols(LHS, RHS_i)
  beta_full[i] = results_i[1, 1]
  print(i) 
  
}

beta_full
beta_full[5] <- 1

T <- NROW(Y) 
N <- NCOL(data)
R <- (1946 - 1926)*12 + 1 
P_0 <-  (1956 - 1946)*12 
P <- T - (R + P_0)
theta <- 0.75
r <- 1
MA_SOP <- 20*12

FC_HA <- matrix(0, nrow = (P_0 + P), ncol = 1)
FC_ECON <- matrix(0, nrow = (P_0 + P), ncol = N)
beta_ECON <- array(0, c((P_0 + P),N,2))
FC_ECON_CT <- matrix(0, nrow = (P_0 + P), ncol = N)
omega_DMSFE <- matrix(0, nrow = (P_0 + P), ncol = N)
FC_OTHER <- matrix(0, nrow = (P_0 + P), ncol = (5 + NROW(theta)))
FC_OTHER_CT <- matrix(0, nrow = (P_0 + P), ncol = (5 + NROW(theta)))

for(i in 1:(P_0+P)){
  
  FC_HA[i] <-  mean(Y[1:(R+(i-1))])
  
  X_t <- data[1:(R+(i-1)-1),]
  Y_t <- Y[2:(R+(i-1))]
  
  for(j in 1:N){
    
    results_i_j <- ols(Y_t, cbind(X_t[,j], matrix(1, nrow = (R+(i-1)-1), ncol = 1)))
    FC_ECON[i,j] <- cbind(data[R+(i-1),j], 1)%*%betas # betas de results_i_j
    
    beta_ECON[i,j,1] <- betas[1,1]
    beta_ECON[i,j,2] <- betastd[1]
    
    ifelse(beta_full[j] > 0,
           ifelse(results_i_j[1] > 0,
                  FC_ECON_CT[i,j] <- FC_ECON[i,j],
                  FC_ECON_CT[i,j] <- FC_HA[i]),
           ifelse(beta_full[j] < 0,
                  ifelse(betas[1,1] < 0,
                         FC_ECON_CT[i,j] <- FC_ECON[i,j],
                         FC_ECON_CT[i,j] <- FC_HA[i]),
                  0))
    
    if(FC_ECON_CT[i,j]<0){
      FC_ECON_CT[i,j] <-  0
    }
  }
  
  if(i > P_0){
    
   # Kitchen sink forecast
   X_t_sink <- ECON_Sink[1:(R+(i-1)-1),]  
   results_t_sink <- ols(Y_t, cbind(as.matrix(X_t_sink), matrix(1, nrow = NROW(ECON_Sink[1:(R+(i-1)-1),]), ncol = 1)))
   FC_OTHER[i, 1] <- as.matrix(cbind(ECON_Sink[(R+(i-1)),],1))%*%as.matrix(results_t_sink[,1])
   
   ifelse(FC_OTHER[i,1] < 0,
          FC_OTHER_CT[i,1] <- 0,
          FC_OTHER_CT[i,1] <- FC_OTHER[i,1])
  
   # j_max <- 3
   # SIC_t <- matrix(ncol = 2)
   # 
   # for(j in 1:j_max){
   #   
   #   select_j <- t(combn(1:NCOL(X_t_sink), j))
   #   
   #   for(k in 1:NROW(select_j)){
   #     
   #     if(j == 1){
   #       
   #       X_t_j_k <- cbind(X_t_sink[1:(R+(i-1)-1), select_j[k]], matrix(1, nrow = (R+(i-1)-1), ncol = 1))
   #       results_t_j_k <- ols(Y_t, as.matrix(X_t_j_k))
   #       resid <-  Y_t - crossprod(t(X_t_j_k),results_t_j_k[,1])
   #       SIC_t_j_k <- log(t(resid)%*%resid/NROW(Y_t)) + log(NROW(Y_t))*NCOL(X_t_j_k)/NROW(Y_t)
   #       
   #       FC_t_j_k <- as.matrix(cbind(ECON_Sink[(R+(i-1)), select_j[k]], 1))%*%as.matrix(results_t_j_k[,1])
   #       print(cbind(SIC_t_j_k, FC_t_j_k))
   #       
   #     } else {
   #       
   #       X_t_j_k <- cbind(X_t_sink[1:(R+(i-1)-1), select_j[k,]], matrix(1, nrow = (R+(i-1)-1), ncol = 1))
   #       results_t_j_k <- ols(Y_t, as.matrix(X_t_j_k))
   #       resid <-  Y_t - crossprod(t(X_t_j_k),results_t_j_k[,1])
   #       SIC_t_j_k <- log(t(resid)%*%resid/NROW(Y_t)) + log(NROW(Y_t))*NCOL(X_t_j_k)/NROW(Y_t)
   #       
   #       FC_t_j_k <- as.matrix(cbind(ECON_Sink[(R+(i-1)), select_j[k,]], 1))%*%as.matrix(results_t_j_k[,1])
   #       print(cbind(SIC_t_j_k, FC_t_j_k))
   #       
   #     }
   #    
   #   }
   # }
   
   # Pooled forecast: simple average
   FC_OTHER[i,2] <- mean(FC_ECON[i,])
   
   ifelse(FC_OTHER[i,2] < 0,
          FC_OTHER_CT[i,2] <- 0,
          FC_OTHER_CT[i,2] <- FC_OTHER[i,2])
   
   # Pooled forecast: DMSFE
   # powers_t <- sort(seq(0, i-2, 1), decreasing = TRUE)
   # m <- t(sum((kronecker(matrix(1, nrow = 1, ncol = N), (theta*matrix(1, nrow = i-1, ncol = 1))^powers_t))*((kronecker(matrix(1, nrow = 1, ncol = N), Y[(R+1):(R+(i-1))]) - FC_ECON[1:(i-1),])^2)))
   # omega <- (m^(-1))/(sum(m^(-1)))
   # FC_OTHER[i,3] <- FC_ECON[i,]*omega
   # 
   # ifelse(FC_OTHER[i,3] < 0,
   #        FC_OTHER_CT[i,3] <- 0,
   #        FC_OTHER_CT[i,3] <- FC_OTHER[i,3])
   
   # Sum-of-the-parts forecast
   FC_OTHER[i,5] <- mean(SopData$E_growth[(R+(i-1)-MA_SOP + 1):(R+(i-1))]) + SopData$DP_SOP[(R+(i-1))] - SopData$r_f[(R+(i-1))]
   
   ifelse(FC_OTHER[i,5] < 0,
          FC_OTHER_CT[i,5] <- 0,
          FC_OTHER_CT[i,5] <- FC_OTHER[i,5])
   
  }
}

beta_ECON <- beta_ECON[(P_0+1):NROW(beta_ECON),,]
actual <- Y[(R+P_0+1):(NROW(Y))]

FC_HA <- FC_HA[(P_0+1):NROW(FC_HA)]
FC_ECON <- FC_ECON[(P_0+1):NROW(FC_ECON),]
FC_ECON_CT <- FC_ECON_CT[(P_0+1):NROW(FC_ECON_CT),]

FC_OTHER <- FC_OTHER[(P_0+1):NROW(FC_OTHER),]
FC_OTHER_CT <- FC_OTHER_CT[(P_0+1):NROW(FC_OTHER_CT),]

# Expansion and recession ------------------------------------------------------
REC <- handbook_data %>% 
       clean_names() %>% 
       select(date_yyyymm, nber_recession_dummies_with_peak_included) %>% 
       filter(nber_recession_dummies_with_peak_included != "NaN", date_yyyymm >=195701) %>% 
       rename(recession = nber_recession_dummies_with_peak_included) %>% 
       mutate(recession = as.numeric(recession))

EXP <- -1*(REC$recession - matrix(1, nrow = nrow(REC), ncol = 1))
index_EXP <- which(EXP == 1)
index_REC <- which(REC$recession == 1)

# ECON -------------------------------------------------------------------------
e_HA <- actual - FC_HA
e_ECON <- kronecker(matrix(1, nrow = 1, ncol = NCOL(FC_ECON)), actual) - FC_ECON
e_ECON_CT <- kronecker(matrix(1, nrow = 1, ncol = NCOL(FC_ECON_CT)), actual) - FC_ECON_CT

CSFE_HA <- cumsum(e_HA^2)
CSFE_ECON <- cumsum(e_ECON^2)
CSFE_ECON_CT <- cumsum(e_ECON_CT^2)

DCSFE_ECON <- kronecker(matrix(1, nrow = 1, ncol = NCOL(FC_ECON)), CSFE_HA) - CSFE_ECON
DCSFE_ECON_CT <- kronecker(matrix(1, nrow = 1, ncol = NCOL(FC_ECON)), CSFE_HA) - CSFE_ECON_CT

R2OS_ECON <- matrix(0, nrow = NCOL(FC_ECON), 6)
R2OS_ECON_CT <- matrix(0, nrow = NCOL(FC_ECON_CT), 6)

for(i in 1:NROW(R2OS_ECON)){
  
  # Overall
  R2OS_ECON[i,1] <- 100*(1 - (sum(e_ECON[,i]^2/sum(e_HA^2))))
  f_i <- e_HA^2 - (e_ECON[,i]^2 - (FC_HA-FC_ECON[,i])^2)
  results_i <- lm(f_i ~ matrix(1, nrow = NROW(f_i), ncol = 1))
  tstat <- coef(summary(results_i))[, "t value"]
  R2OS_ECON[i,2] <- 1 - pnorm(tstat, 0, 1)
  
  R2OS_ECON_CT[i,1] <- 100*(1 - (sum(e_ECON_CT[,i]^2/sum(e_HA^2))))
  f_i <- e_HA^2 - (e_ECON_CT[,i]^2 - (FC_HA-FC_ECON_CT[,i])^2)
  results_i <- lm(f_i ~ matrix(1, nrow = NROW(f_i), ncol = 1))
  tstat <- coef(summary(results_i))[, "t value"]
  R2OS_ECON_CT[i,2] <- 1 - pnorm(tstat, 0, 1)
  
  # Expansion
  R2OS_ECON[i,3] <- 100*(1 - (sum(e_ECON[index_EXP,i]^2/sum(e_HA[index_EXP]^2))))
  f_i <- e_HA[index_EXP]^2 - (e_ECON[index_EXP,i]^2 - (FC_HA[index_EXP]-FC_ECON[index_EXP,i])^2)
  results_i <- lm(f_i ~ matrix(1, nrow = NROW(f_i), ncol = 1))
  tstat <- coef(summary(results_i))[, "t value"]
  R2OS_ECON[i,4] <- 1 - pnorm(tstat, 0, 1)
  
  R2OS_ECON_CT[i,3] <- 100*(1 - (sum(e_ECON_CT[index_EXP,i]^2/sum(e_HA[index_EXP]^2))))
  f_i <- e_HA[index_EXP]^2 - (e_ECON_CT[index_EXP,i]^2 - (FC_HA[index_EXP]-FC_ECON_CT[index_EXP,i])^2)
  results_i <- lm(f_i ~ matrix(1, nrow = NROW(f_i), ncol = 1))
  tstat <- coef(summary(results_i))[, "t value"]
  R2OS_ECON_CT[i,4] <- 1 - pnorm(tstat, 0, 1)
  
  # Recession
  R2OS_ECON[i,5] <- 100*(1 - (sum(e_ECON[index_REC,i]^2/sum(e_HA[index_REC]^2))))
  f_i <- e_HA[index_REC]^2 - (e_ECON[index_REC,i]^2 - (FC_HA[index_REC]-FC_ECON[index_REC,i])^2)
  results_i <- lm(f_i ~ matrix(1, nrow = NROW(f_i), ncol = 1))
  tstat <- coef(summary(results_i))[, "t value"]
  R2OS_ECON[i,6] <- 1 - pnorm(tstat, 0, 1)
  
  R2OS_ECON_CT[i,5] <- 100*(1 - (sum(e_ECON_CT[index_REC,i]^2/sum(e_HA[index_REC]^2))))
  f_i <- e_HA[index_REC]^2 - (e_ECON_CT[index_REC,i]^2 - (FC_HA[index_REC]-FC_ECON_CT[index_REC,i])^2)
  results_i <- lm(f_i ~ matrix(1, nrow = NROW(f_i), ncol = 1))
  tstat <- coef(summary(results_i))[, "t value"]
  R2OS_ECON_CT[i,6] <- 1 - pnorm(tstat, 0, 1)
}

R2OS_ECON
R2OS_ECON_CT

# others -----------------------------------------------------------------------
e_OTHER <- kronecker(matrix(1, nrow = 1, ncol = NCOL(FC_OTHER)), actual) - FC_OTHER
e_OTHER_CT <- kronecker(matrix(1, nrow = 1, ncol = NCOL(FC_OTHER_CT)), actual) - FC_OTHER_CT

CSFE_OTHER <- cumsum(e_OTHER^2)
CSFE_OTHER_CT <- cumsum(e_OTHER_CT^2)

DCSFE_OTHER <- kronecker(matrix(1, nrow = 1, ncol = NCOL(FC_OTHER)), CSFE_HA) - CSFE_OTHER
DCSFE_OTHER_CT <- kronecker(matrix(1, nrow = 1, ncol = NCOL(FC_OTHER)), CSFE_HA) - CSFE_OTHER_CT

R2OS_OTHER <- matrix(0, nrow = NCOL(FC_OTHER) - 1, 6)
R2OS_OTHER_CT <- matrix(0, nrow = NCOL(FC_OTHER_CT) - 1, 6)

for(i in 1:NROW(R2OS_OTHER)){
  
  # Overall
  R2OS_OTHER[i,1] <- 100*(1 - (sum(e_OTHER[,i]^2/sum(e_HA^2))))
  f_i <- e_HA^2 - (e_OTHER[,i]^2 - (FC_HA-FC_OTHER[,i])^2)
  results_i <- lm(f_i ~ matrix(1, nrow = NROW(f_i), ncol = 1))
  tstat <- coef(summary(results_i))[, "t value"]
  R2OS_OTHER[i,2] <- 1 - pnorm(tstat, 0, 1)
  
  R2OS_OTHER_CT[i,1] <- 100*(1 - (sum(e_OTHER_CT[,i]^2/sum(e_HA^2))))
  f_i <- e_HA^2 - (e_OTHER_CT[,i]^2 - (FC_HA-FC_OTHER_CT[,i])^2)
  results_i <- lm(f_i ~ matrix(1, nrow = NROW(f_i), ncol = 1))
  tstat <- coef(summary(results_i))[, "t value"]
  R2OS_OTHER_CT[i,2] <- 1 - pnorm(tstat, 0, 1)
  
  # Expansion
  R2OS_OTHER[i,3] <- 100*(1 - (sum(e_OTHER[index_EXP,i]^2/sum(e_HA[index_EXP]^2))))
  f_i <- e_HA[index_EXP]^2 - (e_OTHER[index_EXP,i]^2 - (FC_HA[index_EXP]-FC_OTHER[index_EXP,i])^2)
  results_i <- lm(f_i ~ matrix(1, nrow = NROW(f_i), ncol = 1))
  tstat <- coef(summary(results_i))[, "t value"]
  R2OS_OTHER[i,4] <- 1 - pnorm(tstat, 0, 1)
  
  R2OS_OTHER_CT[i,3] <- 100*(1 - (sum(e_OTHER_CT[index_EXP,i]^2/sum(e_HA[index_EXP]^2))))
  f_i <- e_HA[index_EXP]^2 - (e_OTHER_CT[index_EXP,i]^2 - (FC_HA[index_EXP]-FC_OTHER_CT[index_EXP,i])^2)
  results_i <- lm(f_i ~ matrix(1, nrow = NROW(f_i), ncol = 1))
  tstat <- coef(summary(results_i))[, "t value"]
  R2OS_OTHER_CT[i,4] <- 1 - pnorm(tstat, 0, 1)
  
  # Recession
  R2OS_OTHER[i,5] <- 100*(1 - (sum(e_OTHER[index_REC,i]^2/sum(e_HA[index_REC]^2))))
  f_i <- e_HA[index_REC]^2 - (e_OTHER[index_REC,i]^2 - (FC_HA[index_REC]-FC_OTHER[index_REC,i])^2)
  results_i <- lm(f_i ~ matrix(1, nrow = NROW(f_i), ncol = 1))
  tstat <- coef(summary(results_i))[, "t value"]
  R2OS_OTHER[i,6] <- 1 - pnorm(tstat, 0, 1)
  
  R2OS_OTHER_CT[i,5] <- 100*(1 - (sum(e_OTHER_CT[index_REC,i]^2/sum(e_HA[index_REC]^2))))
  f_i <- e_HA[index_REC]^2 - (e_OTHER_CT[index_REC,i]^2 - (FC_HA[index_REC]-FC_OTHER_CT[index_REC,i])^2)
  results_i <- lm(f_i ~ matrix(1, nrow = NROW(f_i), ncol = 1))
  tstat <- coef(summary(results_i))[, "t value"]
  R2OS_OTHER_CT[i,6] <- 1 - pnorm(tstat, 0, 1)
}

R2OS_OTHER
R2OS_OTHER_CT
