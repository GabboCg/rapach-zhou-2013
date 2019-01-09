#'----
#' author: Gabriel Cabrera 
#' title:
#' date: 17/12/2018 (last updated: 17/12/2018)
#'---

if(!require("pacman")) iinstall.packages("pacman")
p_load("tidyverse", "microbenchmark", "janitor", "lubridate", "zoo", "sandwich")

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

data <- data %>% 
        select(-equity_premium)

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
  
}

beta_ECON <- beta_ECON[(P_0+1):NROW(beta_ECON),,]
actual <- Y[(R+P_0+1):(NROW(Y))]

FC_HA <- FC_HA[(P_0+1):NROW(FC_HA)]
FC_ECON <- FC_ECON[(P_0+1):NROW(FC_ECON),]
FC_ECON_CT <- FC_ECON_CT[(P_0+1):NROW(FC_ECON_CT),]

REC <- handbook_data %>% 
       clean_names() %>% 
       select(date_yyyymm, nber_recession_dummies_with_peak_included) %>% 
       filter(nber_recession_dummies_with_peak_included != "NaN", date_yyyymm >=195701) %>% 
       rename(recession = nber_recession_dummies_with_peak_included) %>% 
       mutate(recession = as.numeric(recession))

EXP <- -1*(REC$recession - matrix(1, nrow = nrow(REC), ncol = 1))
index_EXP <- which(EXP == 1)
index_REC <- which(REC$recession == 1)

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


