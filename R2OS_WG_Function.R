#'---
#' author: Gabriel Cabrera
#' title: R2OS Welch & Goyal
#' date: 18/12/2018 (last updated: 19/12/2018)
#'---

R2OS <- function(Y, X, H, init, fin, space){
  
  ols <- function(Y,X){
    
    betas <<- solve(t(X)%*%X)%*%(t(X)%*%Y) 
    
    e <<- Y - X%*%betas
    variance <- t(e)%*%e
    vcv <- 1/(NROW(X)-NCOL(X))*(as.numeric(variance)*solve(t(X)%*%X))
    betastd <<- sqrt(diag(vcv))
    
    regr <<- cbind(betas, as.matrix(betastd))
    
  }
  
  Y <- Y
  
  T <- NROW(Y)
  LHS <- Y[(1+H):T]
  
  beta_full <- matrix(0, nrow = NCOL(Y), ncol = 1)
  
  for(i in 1:NCOL(X)){
    
    RHS_i <- cbind(as.matrix(X[1:(T-H),i]), matrix(1,nrow = (T-H), ncol = 1))
    results_i <- ols(LHS, RHS_i)
    beta_full[i] = results_i[1, 1]
    # print(i) 
    
  }
  
  T <- NROW(Y) 
  N <- NCOL(X)
  R <- (fin - init)*12 + 1 
  P_0 <-  ((fin + space) - fin)*12 
  P <- T - (R + P_0)
  theta <- 0.75
  r <- 1
  MA_SOP <- (fin - init)*12
  
  FC_HA <- matrix(0, nrow = (P_0 + P), ncol = 1)
  FC_ECON <- matrix(0, nrow = (P_0 + P), ncol = N)
  beta_ECON <- array(0, c((P_0 + P),N,2))
  FC_ECON_CT <- matrix(0, nrow = (P_0 + P), ncol = N)
  
  FC_OTHER <- matrix(0, nrow = (P_0 + P), ncol = (5 + NROW(theta)))
  FC_OTHER_CT <- matrix(0, nrow = (P_0 + P), ncol = (5 + NROW(theta)))
  
  for(i in 1:(P_0+P)){
    
    FC_HA[i] <-  mean(Y[1:(R+(i-H))])
    
    X_t <- X[1:(R+(i-H)-H),]
    Y_t <- Y[(1+H):(R+(i-H))]
    
    for(j in 1:N){
      
      results_i_j <- ols(Y_t, cbind(X_t, matrix(1, nrow = (R+(i-H)-H), ncol = 1))) # X_t[,j]
      FC_ECON[i,j] <- cbind(X[R+(i-H),j], 1)%*%betas # betas de results_i_j
      
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
  
  beta_ECON_AUX <- beta_ECON[(P_0+1):NROW(beta_ECON),,]
  
  beta_ECON <- array(0,c(NROW(beta_ECON_AUX),1,NCOL(beta_ECON_AUX)))
  
  beta_ECON[,,1] <- beta_ECON_AUX[,1]
  beta_ECON[,,2] <- beta_ECON_AUX[,2]
  
  actual <- Y[(R+P_0+1):(NROW(Y))]
  
  FC_HA <- FC_HA[(P_0+1):NROW(FC_HA)]
  FC_ECON <- FC_ECON[(P_0+1):NROW(FC_ECON),]
  FC_ECON_CT <- FC_ECON_CT[(P_0+1):NROW(FC_ECON_CT),]
  
  e_HA <- actual - FC_HA
  e_ECON <- actual - FC_ECON
  e_ECON_CT <- actual - FC_ECON_CT
  
  CSFE_HA <- cumsum(e_HA^2)
  CSFE_ECON <- cumsum(e_ECON^2)
  CSFE_ECON_CT <- cumsum(e_ECON_CT^2)
  
  DCSFE_ECON <- CSFE_HA - CSFE_ECON
  DCSFE_ECON_CT <-  CSFE_HA - CSFE_ECON_CT
  
  R2OS_ECON <- matrix(0, nrow = NCOL(FC_ECON), 2)
  R2OS_ECON_CT <- matrix(0, nrow = NCOL(FC_ECON_CT), 2)
  
  for(i in 1:NROW(R2OS_ECON)){
    
    R2OS_ECON[i,1] <- 100*(1 - (sum(e_ECON^2/sum(e_HA^2))))
    f_i <- e_HA^2 - (e_ECON^2 - (FC_HA-FC_ECON)^2) # FC_ECON[,i]
    results_i <- lm(f_i ~ matrix(1, nrow = NROW(f_i), ncol = 1))
    tstat <- coef(summary(results_i))[, "t value"]
    R2OS_ECON[i,2] <- 1 - pnorm(tstat, 0, 1)
    
    R2OS_ECON_CT[i,1] <- 100*(1 - (sum(e_ECON_CT^2/sum(e_HA^2))))
    f_i <- e_HA^2 - (e_ECON_CT^2 - (FC_HA-FC_ECON_CT)^2) # FC_ECON_CT[,i]
    results_i <- lm(f_i ~ matrix(1, nrow = NROW(f_i), ncol = 1))
    tstat <- coef(summary(results_i))[, "t value"]
    R2OS_ECON_CT[i,2] <- 1 - pnorm(tstat, 0, 1)
    
  }
  
  R2OS_results <<- cbind.data.frame(without = R2OS_ECON, with = R2OS_ECON_CT) %>% 
                   rename(without_R2OS = "without.1",
                          without_p_value = "without.2",
                          with_R2OS = "with.1",
                          with_p_value = "with.2")
  
}



