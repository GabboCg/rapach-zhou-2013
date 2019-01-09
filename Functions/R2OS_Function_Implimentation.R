#'---
#' author: Gabriel Cabrera
#' title: R2OS Welch & Goyal implementation 
#' date: 18/12/2018 (last updated: 19/12/2018)
#'---

if(!require("pacman")) install.packages("pacman")
p_load("tidyverse", "microbenchmark", "janitor", "lubridate", "gt")

source("R2OS_WG_Function.R")

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
        select(equity_premium, DP, DY, EP, DE, SVAR, BM, NTIS, TBL, LTY, LTR, 
               TMS, DFY, DFR, INFL_lag) %>% 
        na.omit()

Y <- data %>% 
     select(equity_premium) %>% 
     as.matrix()

matrix_value <- matrix(0, nrow = NCOL(data)-1, ncol = 4)

for(i in 2:NCOL(data)){
  
    R2OS(Y,data[i], 1, 1926, 1946, 10)
    matrix_value[i-1,] <- as.matrix(R2OS_results)

}

label_row <- c("log(DP)","log(DY)","log(EP)","log(DE)","SVAR","BM",
               "NTIS","TBL","LTY","LTR","TMS","DFY","DFR","INFL")

label_col <- c("without_R2OS", "without_R2OS_p_value",
               "with_R2OS", "with_R2OS_p_value")

data.frame(predictor = label_row, matrix_value) %>% 
  gt() %>% 
  tab_header(title = md("**Forecasting Stock Returns**"),
    subtitle = md("*Replication Table 1*")) %>% 
  tab_spanner(
    label = html("With <br> Restriction"),
    columns = vars(without_R2OS, without_R2OS_p_value)
  ) %>%
  tab_spanner(
    label = html("Without <br> Restriction"),
    columns = vars(with_R2OS, with_R2OS_p_value)
  ) %>% 
  cols_label(
    predictor = md("*Predictors*"),
    without_R2OS = html("R<sup>2</sup><sub>OS</sub>"),
    without_R2OS_p_value = "P-Value",
    with_R2OS = html("R<sup>2</sup><sub>OS</sub>"),
    with_R2OS_p_value = "P-Value"
  ) %>% 
  fmt_number(
    columns = 2:5,
    decimals = 2,
    use_seps = TRUE
  ) %>%   
  tab_source_note(
    source_note = html("<i> Reference </i>: 
                        Rapach, D., & Zhou, G. (2013). Forecasting stock returns. <br>
                        In Handbook of economic forecasting (Vol. 2, pp. 328-383). 
                        Elsevier.")
  )

