library(MASS)
library(sandwich)

### data <- mocc_v2
### given the MOCC individual-level data, create aggregated data (AD)
### define variables for ITS analysis
### D = dummy variable indicating pre-MOCC (=0) and post-MOCC (=1)
### P = time passed since MOCC initiation (P is 0 for pre-MOCC)
### week2 = week**2, quadratic term
### P2 = P**2, quadratic term
### Analysis: negative binomial regression with Newey-West standard error estimation
### Return values:
### fit_nb_qd_coefficients: coefficients and s.e. estimates under quadratic model
### fit_nb_coefficients: coefficients and s.e. estimates with linear terms of week
### W: Wald test statistic for the quadratic terms of week and P
### pvalue: p-value of the above Wald test
### Two figures: predicted values of the IFTs with the MOCC initiation (blue)
###   and predicted values of the potential IFTs had the MOCC not been initiated (red)

### include interaction effects between surge and long term MOCC effects (P or P2)
### include surge variable and season in the model
mocc_routines_v4 <- function(data, figurename1, 
                             figurename2, model_selected)
{
  ### exclude observations with missing surge data
  data <- data[!is.na(data$surge),]
  tabweek <- table(data$week)
  rownames(tabweek)
  week <- as.integer(rownames(tabweek))
  ### number of transfers from week -13 to 70
  ntrans <- tabweek
  
  
  ### file name to save the results 
  flename_quad <- paste(substring(figurename1, 1, 
                                  nchar(figurename1)-5),
                        ".txt",sep="")
  
  ### distribution of season in each week from -13 to 70
  season <- table(data$week, data$season)
  ##colnames(season)
  
  ### create AD: week, number of transfers
  data_ad <- data.frame(cbind(week,ntrans))
  
  ## season: winter as reference
  data_ad$Spring <- as.vector(season[,2]/ntrans)
  data_ad$Summer <- as.vector(season[,3]/ntrans)
  data_ad$Fall <- as.vector(season[,1]/ntrans)
  
  ### D = dummy variable indicating pre-MOCC (=0) and post-MOCC (=1)
  data_ad$D <- 1*(data_ad$week >= 1)
  
  ### P = time passed since MOCC initiation (P is 0 for pre-MOCC)
  data_ad$P <- data_ad$week * (data_ad$week >= 1)
  
  ### week2 = week**2, quadratic term
  data_ad$week2 <- data_ad$week**2
  
  ### P2 = P**2, quadratic term
  data_ad$P2 <- data_ad$P**2
  
  ### surge is the average of the surge values in each relative week 
  data_ad$surge <- rep(0, nrow(data_ad))
  for (i in 1:nrow(data_ad))
  {
    data_ad$surge[i] <- mean(data$surge[data$week == data_ad$week[i]])
  }
  
  ### create design matrix
  X <- cbind(rep(1, nrow(data_ad)), data_ad$week, 
             data_ad$week2, data_ad$D, 
             data_ad$P, data_ad$P2, data_ad$surge, data_ad$Spring,
             data_ad$Summer, data_ad$Fall)
  ### Model 1: include interaction effects between surge and (D, P, P2)
  ###           quadratic term of P
  fit_nb1 <- glm.nb(ntrans ~ week + week2 + D + P + P2 
                   + surge + D*surge + P*surge + P2*surge 
                   + Spring + Summer + Fall, #+ P*surge + P2*surge, 
                   data = data_ad)
  
  ### Model 2: include interaction effects between surge and D
  ###           quadratic term of P
  fit_nb2 <- glm.nb(ntrans ~ week + week2 + D + P + P2 
                   + surge + D*surge  
                   + Spring + Summer + Fall, #+ P*surge + P2*surge, 
                   data = data_ad)
  
  ### Model 3: include interaction effects between surge and (D, P)
  ###           linear term P
  fit_nb3 <- glm.nb(ntrans ~ week + week2 + D + P  
                    + surge + D*surge + P*surge  
                    + Spring + Summer + Fall, #+ P*surge + P2*surge, 
                    data = data_ad)
  ### Model 4: include interaction effects between surge and D
  ###           linear term P
  fit_nb4 <- glm.nb(ntrans ~ week + week2 + D + P  
                    + surge + D*surge   
                    + Spring + Summer + Fall, #+ P*surge + P2*surge, 
                    data = data_ad)
  
  ### use AIC to choose the model
  minaic <- min(fit_nb1$aic, fit_nb2$aic, fit_nb3$aic, fit_nb4$aic)
  
  model1_par <- c("D", "P", "P2", "surge", "D:surge", "P:surge", "P2:surge")
  model2_par <- c("D", "P", "P2", "surge", "D:surge")
  model3_par <- c("D", "P", "surge", "D:surge", "P:surge")
  model4_par <- c("D", "P", "surge", "D:surge")
  
  ### select the model with the smallest AIC
  if (is.null(model_selected)) {
  model_selected <- 0
  if (minaic == fit_nb1$aic) {
    model_selected <- 1
    fit_nb <- fit_nb1
    model_par <- model1_par
  } else if (minaic == fit_nb2$aic) {
    model_selected <- 2
    fit_nb <- fit_nb2
    model_par <- model2_par
  } else if (minaic == fit_nb3$aic) {
    model_selected <- 3
    fit_nb <- fit_nb3
    model_par <- model3_par
  } else  {
    model_selected <- 4
    fit_nb <- fit_nb4
    model_par <- model4_par
  }
  } else {
    if (model_selected == 1) {
      fit_nb <- fit_nb1
      model_par <- model1_par
    } else if (model_selected == 2) {
      fit_nb <- fit_nb2
      model_par <- model2_par
    } else if (model_selected == 3) {
      fit_nb <- fit_nb3
      model_par <- model3_par
    } else if (model_selected == 4) {
      fit_nb <- fit_nb4
      model_par <- model4_par
    } else {
      print("Error: must choose a model number between 1 and 4.")
    }
  }
  ### fit model 3 (for all subgroup analysis)
  # fit_nb <- fit_nb3
  # model_par <- model3_par
  
  ### use the Newey West method to estimate standard errors
  covmat <- NeweyWest(fit_nb)
  sqrt(diag(covmat))
  sumry_fit_nb <- summary(fit_nb)
  sumry_fit_nb$coefficients
  sumry_fit_nb$coefficients[,2] <- sqrt(diag(covmat))
  sumry_fit_nb$coefficients[,3] <- sumry_fit_nb$coefficients[,1]/
    sumry_fit_nb$coefficients[,2]
  sumry_fit_nb$coefficients[,4] <- 2*pnorm(-abs(sumry_fit_nb$coefficients[,3]))
  
  ### coefficients estimates of the selected model
  fit_nb_coefficients <- sumry_fit_nb$coefficients
  
  ### Rate ratio estimates, lower and upper limits of the 95% CI
  fit_nb_exp_coefficients <- exp(fit_nb_coefficients[,1:3])
  fit_nb_exp_coefficients[,2] <- exp(fit_nb_coefficients[,1]
                                        -1.96*fit_nb_coefficients[,2])
  fit_nb_exp_coefficients[,3] <- exp(fit_nb_coefficients[,1]
                                        +1.96*fit_nb_coefficients[,2])
  colnames(fit_nb_exp_coefficients)[2:3] <- c("95% CI (lower)",
                                                 "95% CI (upper)")
  
  ### other parameter of interest
  ### 1. surge effect before MOCC
  ### 2. surge effect after MOCC at week 0
  ### 3. surge effect after MOCC at week 10
  ### 4. surge effect after MOCC at week 20
  ### 5. surge effect after MOCC at week 30
  ### 6. surge effect after MOCC at week 40
  ### 7. surge effect after MOCC at week 50
  ### 8. surge effect after MOCC at week 60
  ### 9. surge effect after MOCC at week 70
  ### 10: MOCC immediate effect with surge=1 (1 - mean(surge))
  ### 11: MOCC immediate effect with surge=2 (2 - mean(surge))
  ###...
  ### 19: MOCC immediate effect with surge=10 (10 - mean(surge))
  ### 20: MOCC long term effect with surge=1 (1 - mean(surge))
  ### ...
  ### 29: MOCC long term effect with surge=10 (10 - mean(surge))
  
  model_par
  betaest <- fit_nb_coefficients[model_par,1]
  beta_covest <- covmat[model_par,
                        model_par]
  ### Model 1: 
  ### > model1_par
  ### "D"        "P"        "P2"       "surge"    "D:surge"  "P:surge"  "P2:surge"
  ### > model2_par
  ###  "D"       "P"       "P2"      "surge"   "D:surge"
  ### > model3_par
  ### "D"       "P"       "surge"   "D:surge" "P:surge"
  ### model4_par
  ### "D"       "P"       "surge"   "D:surge"
  C <- matrix(nrow=29, ncol=length(model_par), byrow=T)
  colnames(C) <- model_par
  C[1,1:3] <- c(0,0,1) # "D" "P" "surge" 
  C[2:9,1:3] <- cbind(rep(1,8), (0:7)*10, rep(1,8))
  C[10:19,1] <- rep(1, 10)
  C[10:19,2] <- rep(0, 10)
  C[10:19,3] <- c(1:10) - meansurge
  C[20:29,1] <- rep(0, 10)
  C[20:29,2] <- rep(1, 10)
  C[20:29,3] <- c(1:10) - meansurge
  
  if (model_selected == 1) {
   C[,4] <- C[,3]
   C[,3] <- C[,2]*C[,2]
   C[,5] <- C[,1]*C[,4]
   C[,6] <- C[,2]*C[,4]
   C[,7] <- C[,3]*C[,4]
   C[2:9,1:3] <- 0
   C[10:19,4] <- 0
   C[20:29,c(3,4,7)] <- 0
  } else if (model_selected == 2) {
    C[,4] <- C[,3]
    C[,3] <- C[,2]*C[,2]
    C[,5] <- C[,1]*C[,4]
    C[2:9,1:3] <- 0
    C[10:19,4] <- 0
    C[20:29,c(3,4)] <- 0
  } else if (model_selected == 3) {
    C[,4] <- C[,1]*C[,3]
    C[,5] <- C[,2]*C[,3]
    C[2:9,1:2] <- 0
    C[10:19,3] <- 0
    C[20:29,c(3)] <- 0
  } else if (model_selected == 4) {
    C[,4] <- C[,1]*C[,3]
    C[2:9,1:2] <- 0
    C[10:19,3] <- 0
    C[20:29,c(3)] <- 0
  }
  
  C_betaest <- C %*% betaest
  C_covest <- C %*% beta_covest %*% t(C)
  
  parest <- matrix(nrow=29, ncol=4, byrow = T)
  colnames(parest) <- colnames(fit_nb_coefficients)
  rownames(parest) <- c("surge (before MOCC)",
                           "surge (after MOCC at week 0)",
                        "surge (after MOCC at week 10)",
                        "surge (after MOCC at week 20)",
                        "surge (after MOCC at week 30)",
                        "surge (after MOCC at week 40)",
                        "surge (after MOCC at week 50)",
                        "surge (after MOCC at week 60)",
                        "surge (after MOCC at week 70)",
                           "Immed. MOCC with surge=1",
                           "Immed. MOCC with surge=2",
                           "Immed. MOCC with surge=3",
                           "Immed. MOCC with surge=4",
                           "Immed. MOCC with surge=5",
                           "Immed. MOCC with surge=6",
                           "Immed. MOCC with surge=7",
                           "Immed. MOCC with surge=8",
                           "Immed. MOCC with surge=9",
                           "Immed. MOCC with surge=10",
                        "Long-term MOCC with surge=1",
                        "Long-term MOCC with surge=2",
                        "Long-term MOCC with surge=3",
                        "Long-term MOCC with surge=4",
                        "Long-term MOCC with surge=5",
                        "Long-term MOCC with surge=6",
                        "Long-term MOCC with surge=7",
                        "Long-term MOCC with surge=8",
                        "Long-term MOCC with surge=9",
                        "Long-term MOCC with surge=10")
  parest[,1] <- C_betaest
  parest[,2] <- sqrt(diag(C_covest))
  parest[,3] <- parest[,1]/parest[,2]
  parest[,4] <- 2*pnorm(-abs(parest[,3]))
  
  ### rate ratio                        
  parest_exp <- exp(parest[,1:3])
  parest_exp[,2] <- exp(parest[,1]-1.96*parest[,2])
  parest_exp[,3] <- exp(parest[,1]+1.96*parest[,2])
  colnames(parest_exp)[2:3] <- c("95% CI (lower)",
                                    "95% CI (upper)")
  
  
  ### plot observed and predicted values
  predicted <- exp(predict.glm(fit_nb))
  
  # fit_nb$coefficients[c(4,5,10,11)] <- 0
  # predicted2 <- exp(X[,-6] %*% fit_nb$coefficients[1:9])
  # 
  data_ad$log10_ntrans <- log(data_ad$ntrans)/log(10)
  write.csv(cbind(data_ad$week, data_ad$surge + meansurge, 
                  data_ad$ntrans, predicted),
            file=flename_quad, quote=F, 
            row.names = FALSE)
  
  # write.csv(cbind(data_ad$week, data_ad$surge + meansurge, 
  #                 data_ad$ntrans, predicted, predicted2),
  #           file=flename_quad, quote=F, 
  #           row.names = FALSE)
  
  ggplot(data_ad,aes(week, ntrans)) + 
    xlab("Relative time from MOCC initiation (week)") +
    ylab("Number of IFTs") +
    geom_point() +
    geom_point(aes(week, predicted), color="blue") 
  ggsave(figurename1)
  
  ggplot(data_ad,aes(week, log10_ntrans)) + 
    xlab("Relative time from MOCC initiation (week)") +
    ylab(expression(paste(log[10], "(Number of IFTs)",sep=""))) +
    geom_point() +
    geom_point(aes(week, log(predicted)/log(10)), color="blue") 
  ggsave(figurename1)
  
  return(list(fit_nb=fit_nb, covmat=covmat, data_ad=data_ad,
              fit_nb_coefficients=fit_nb_coefficients,
              fit_nb_exp_coefficients=fit_nb_exp_coefficients,
              parest=parest,
              parest_exp=parest_exp,
              aic=c(fit_nb1$aic, fit_nb2$aic, fit_nb3$aic, fit_nb4$aic)))
}

### calculate the predicted IFT rates and 95% CIs
mocc_pred <- function(fit_nb, covmat, data_ad, surge, meansurge, season)
{
  newdata <- cbind(rep(1, nrow(data_ad)),data_ad$week, data_ad$week2, 
                   data_ad$D, data_ad$P, 
                   data_ad$P2, rep(surge- meansurge, nrow(data_ad)),
                   rep(0, nrow(data_ad)), 
                   rep(0, nrow(data_ad)), rep(0, nrow(data_ad)))
  if (season == "Spring")
  { 
    newdata[,8] <- 1
  } else if (season == "Summer")
  {
    newdata[,9] <- 1
  } else if (season == "Fall")
  {
    newdata[,10] <- 1
  }
  
  newdata <- cbind(newdata, newdata[,4]*newdata[,7],
                   newdata[,5]*newdata[,7],
                   newdata[,6]*newdata[,7])
  colnames(newdata) <- c("(Intercept)", "week", "week2", "D", "P", "P2", 
                         "surge", "Spring", "Summer", "Fall", "D:surge",
                         "P:surge", "P2:surge")
  
  preds_lp <- newdata[,names(fit_nb$coefficients)] %*% fit_nb$coefficients
  preds_lp_cov <- newdata[,names(fit_nb$coefficients)] %*% covmat %*% 
    t(newdata[,names(fit_nb$coefficients)])
  sqrt(diag(preds_lp_cov))
  
  preds <- cbind(data_ad$week, exp(preds_lp), exp(preds_lp - 1.96*sqrt(diag(preds_lp_cov))),
                 exp(preds_lp + 1.96*sqrt(diag(preds_lp_cov))),
                 exp(preds_lp)*sqrt(diag(preds_lp_cov)))
  
  colnames(preds) <- c("week", "predicted", "ll (95%)", "ul (95%)", "s.e.")
  
  return(preds)
}


