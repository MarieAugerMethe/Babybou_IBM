##
# Date created: September 25th, 2013
# Authors: Marie Auger-Methe, Ulrike Schlaegel, Craig DeMars
# Please cite our paper if you use our script:
#   DeMars, C., M. Auger-Methe, U. Schlaegel, S. Boutin, (Published online) 
#   Inferring Parturition and Neonate Survival from Movement Patterns of Female Ungulates.
#   Ecology and Evolution. DOI: 10.1002/ece3.785
# For an in-depth explanation of the code see:
#   Appendix S2. Likelihood Functions Used in the Individual-based Method. 
#   Appendix S3. R Code for the Individual-based Method.
# from the supporting information:
#   ece3785-sup-0001-AppendixS1-S4.docx


# The main difference between this code and the one presented in IBM.R and the appendix S3
# is this function calculated the full likelihood instead of the approximation.
# Because, in this case, we are numerically estimating all of the parameters,
# this code will take longer to run.
# fullnll1 now replaces nllk and nll1 now optimizess fullnll1 instead of nllk.
# This mean that now both lambda and k are estimated simultaneously numerically.
# Now we have fullnll2 and changed nll2 so it optimizes simultaneously lambda and k.
# Also, nll1 and nll2 returns now both NLL value and the MLE values.
# In the function mnll3M we now have:
#  columns for BIC values in resCA
#	 small changes in sections "Calf survived" and "Calf lost"
#	    (MLE values return from function nll1 and nll2 used
# 		 instead of optimizing again for best BP)
#	calculate BIC values and get model with min BIC





#################################################################
# Functions that calculate the negative log likelihood (nll) of 2 of the 3 models
# and the nll of k, one of the parameter used to calculate the nll of these models.
# Note that all models assume 
# that the step length are distributed with an exponential distribution.
# They differ in how the main parameter of the exponential distribution
# varies through the time series.



# nll1 calculates the nll of a female with a calf that survived.
# This model assumes that the mean step length is constant prior to the birth of a calf.
# Once the calf is born the mean step length decrease to 1/(l_a * k), 
# which should be really close to 0.
# The mean step length then linearly increases 
# to reach the same mean step length as prior to calf birth.
# The speed of the increase depends on the time it takes to recover the normal movement,
# which is represented by k.
# Once the mean step length reaches the value from before the birth of the calf,
# it remains at this value.
# The parameters to estimate are:
# 1 BP, l_a, and k

# fullnll1 evaluates the negative log likelihood at parameter values (lambda, k) under model 1

fullnll1 <- function(param, SL_a, SL_b, ti_b){
  
  # param is the vector of parameters of interest
  l_a <- param[1]
  k <- param[2]
  
  # likelihood part t=1,...,BP
  LL_a <- dexp(SL_a, l_a, log=T)
  
  # likelihood part t=BP+1,...,T
  z <- 1/l_a
  z_b <- (ti_b-ti_b[1]+1)*z/k
  z_b[z_b > z ] <- z
  LL_k <- dexp(SL_b, 1/z_b, log=T)
  
  # sum all negative log likelihoods from all parts
  return(-sum(LL_a) - sum(LL_k))
}


# nll1 minimizes fullnll1 over possible values for (lambda, k)
# This results in the maximum likelihood estimates for (lambda, K)
nll1 <- function(BP, SL, ti, kc){
  
  n <- length(SL)
  ###
  # Divides the time series into two sections:
  # a: before the birth of the calf
  # b: after the birth of the calf
  SL_a <- SL[1:BP]
  s_b <- (BP+1):n
  SL_b <- SL[s_b] 
  ti_b <- ti[s_b]
  
  # Numerically estimate the MLE of lambda and k
  # and obtain the corresponding negative log-likelihood values
  # start values for optim are:
  # for lambda: the mean step length for t=1,...,BP
  # for k: 4 weeks (equivalent thereof in steps)
  startvalues <- c(length(SL_a)/sum(SL_a), 4*24/(dint*subs)*7)
  mnll_b <- optim(startvalues, fullnll1, method="L-BFGS-B", SL_a=SL_a, SL_b=SL_b, ti_b=ti_b, lower=c(0.0005, kc[1]), upper=c(1, kc[2]))
  
  ###
  
  # Getting the total nll for the whole time series
  nLL <- mnll_b$value
  l_a <- mnll_b$par[1]
  k <- mnll_b$par[2]
  return(list("nLL"=nLL, "l_a"=l_a, "k"=k))
}






# nll2 calculates the nll of a female with a calf that died.
# This model assumes that the mean step length is constant prior to the birth of a calf.
# Once the calf is born the mean step length decrease to 1/(l_a * k), 
# which should be really close to 0.
# The mean step length then linearly increases 
# to reach the same mean step length as prior to calf birth.
# The speed of the increase depends on the time it takes to recover the normal movement,
# which is represented by k.
# Once the mean step length reaches the value from before the birth of the calf,
# it remains at this value.
# The k represent the time it takes to recover the normal movement if the calf survives.
# The previous assumptions are the same as for nll1. 
# This model further assumes
# that the female immediately recovers her normal movement when the calf dies. 
# This abrupt switch is represented by a second BP.
# The parameters to estimate are:
# 2 BPs, l_a, and k


# fullnll2 evaluates the negative log likelihood at parameter values (lambda, k) under model 2
fullnll2 <- function(param, SL_a, SL_b, SL_c, ti_b){
  # param is the vector of parameters of interest
  l_a <- param[1]
  k <- param[2]
  
  # likelihood part t=1,...,BP1
  LL_a <- dexp(SL_a, l_a, log=T)
  
  # likelihood part t=BP1+1,...,BP2
  z <- 1/l_a
  z_b <- (ti_b-ti_b[1]+1)*z/k
  z_b[z_b > z ] <- z
  LL_k <- dexp(SL_b, 1/z_b, log=T)
  
  # likelihood part t=BP2+1,...,T
  LL_c <- dexp(SL_c, l_a, log=T)
  
  # sum all negative log likelihoods from all parts
  return(-sum(LL_a) - sum(LL_k) - sum(LL_c))
}


# nll2 minimizes fullnll2 over possible values for (lambda, k)
# This results in the maximum likelihood estimates for (lambda, K)
nll2 <- function(BP, SL, ti, kc){
  
  n <- length(SL)
  
  # Divides the time series into three sections:
  # a: before the birth of the calf
  # b: after the birth of the calf but before it dies
  # c: after the death of the calf
  SL_a <- SL[1:BP[1]]
  s_b <-  (BP[1]+1):BP[2]
  SL_b <- SL[s_b]
  ti_b <- ti[s_b]
  SL_c <- SL[(BP[2]+1):n]
  
  # Numerically estimate the MLE of lambda and k
  # and obtain the corresponding negative log-likelihood values
  # start values for optim are:
  # for lambda: the mean step length for t=1,...,BP
  # for k: 4 weeks (equivalent thereof in steps)
  startvalues <- c(length(SL_a)/sum(SL_a), 4*24/(dint*subs)*7)
  
  mnll_b <- optim(startvalues, fullnll2, method="L-BFGS-B", SL_a=SL_a, 
                  SL_b=SL_b, SL_c=SL_c, ti_b=ti_b, 
                  lower=c(0.0005, kc[1]), upper=c(1, kc[2]))
  
  # Getting the total nll for the whole time series

  nLL <- mnll_b$value
  l_a <- mnll_b$par[1]
  k <- mnll_b$par[2]
  return(list("nLL"=nLL, "l_a"=l_a, "k"=k))
}






# Function minimise the likelihood of the 3 Models

mnll3M <- function(movF, int, kcons){
  # movF: list with at least 3 elements
  #   1. SL: numeric vector that contains the step lengths
  #           measured at regular time intervals.
  #           NAs are not allowed.
  #   2. ti: integer vector that identifies
  #           the time of each step length of SL.
  #           It is the index of the times found in tp.
  #           ti and SL should be of exactly the same length.
  #           NAs are not allowed.
  #   3. tp: POSIXct vector that identifies the real date and time of the SL.
  #           The missing steps should be represented with NAs.
  #           For example if you have locations
  #           (0,0) (0,1) (0,3) (1,3) (7,3)
  #           taken at:
  #           01:00, 02:00, 03:00, 08:00, 09:00
  #           SL = (1,2,6)
  #           ti = (1,2,4)
  #           tp = (01:00, 02:00, NA, 08:00, NA) 
  #           (Although excluded for clarity, tp should include the date)
  #           We recommend that the time series only included the time period
  #           relevant to the birth and death of the calf.
  #           We only included movement from xx to xx.
  
  # int: integer value indicating the minimum number of steps needed between
  #     the beginning of the time series and the first BP (birth of calf),
  #     the last BP (death of calf) and the end of the time series, 
  #     and between the two BPs.
  #     This is required in part because you need some steps
  #     in most sections of the time series to estimate parameters.
  
  # kcons: numeric vector with two values that contains
  #     the minimum and maximum values c(kmin,kmax) for parameter k.
  #     k represents the time it takes for the female with a calf that survives
  #     to recover her normal movement. 
  #     We contrained our values to be the equivalent in steps of 3 and 6 weeks.
  #     kcons is required as it is needed by the function optimize()
  #     but to place no informed restriction, 
  #     you can make kcons <- c(0,vlv), where vlv is a very large value
  
  # Results for model comparison
  resCA <- matrix(NA, 1, ncol=12)
  colnames(resCA) <- c("n", "mnll_0","mnll_1","mnll_2",
                       "AIC0", "AIC1", "AIC2", "BM",
                       "BIC0", "BIC1", "BIC2", "BM_BIC")
  
  # BPs index and actual date and time
  BPs <- data.frame(matrix(NA, 1, ncol=6))
  colnames(BPs) <- c("BP1c", "BP2c", "BP2l", "iBP1c", "iBP2c", "iBP2l")
  
  # Parameters (other than the BP) of each model
  mpar <- matrix(NA, 1, 5)
  colnames(mpar) <- c('b0','b1','b2','k1', 'k2')
  
  # Sample size
  resCA[1] <- length(movF$SL)
  
  ##
  # M0: No calf
  # This model assumes that the movement pattern is constant for the whole time series
  # and thus the mean step length (b0) is constant.
  # The only parameter estimated is b0, which is the inverse of the rate.
  # It has a analytical solution and 
  # thus we can easily get the minimum negative log likelihood (mnll)
  mpar[1] <- sum(movF$SL)/resCA[1]  #b0
  resCA[2] <- -sum(dexp(movF$SL,1/mpar[1],log=T)) #mnll0
  
  ##
  # Calf survived
  BP1ser <- int:(resCA[1]-int) # Serie of all possible BPs
  temp_1 <- lapply(BP1ser, nll1, SL=movF$SL, ti=movF$ti, kc=kcons)
  NLL_1 <- as.vector(unlist(temp_1)[seq(1,length(unlist(temp_1)),3)]) 
  MNLL_1_i <- which.min(NLL_1) 
  resCA[3] <- NLL_1[MNLL_1_i] # mnll_1
  BPs[4] <- BP1ser[MNLL_1_i]
  BPs[1] <- as.character(movF$tp[movF$ti[BPs[[4]]]]) #BP1c
  mpar[2] <- 1/temp_1[[MNLL_1_i]]$l_a #b1
  mpar[4] <- temp_1[[MNLL_1_i]]$k
  
  ##
  # Calf lost
  # Getting all possible combination of BPs
  BP2ser <- combn(int:(resCA[1]-int),2)
  BP2ser <- BP2ser[,diff(BP2ser) >= int]
  BP2ser <- split(t(BP2ser),1:ncol(BP2ser))
  
  # Applying the model 2 on all possible BPs
  temp_2 <- lapply(BP2ser,nll2,SL=movF$SL, ti=movF$ti, kc=kcons)
  NLL_2 <- as.vector(unlist(temp_2)[seq(1,length(unlist(temp_2)),3)])
  MNLL_2_i <- which.min(NLL_2)
  resCA[4] <- NLL_2[MNLL_2_i] # mnll_2
  BPs[5:6] <- BP2ser[[MNLL_2_i]]
  BPs[2] <- as.character(movF$tp[movF$ti[BPs[[5]]]]) #BP2c
  BPs[3] <- as.character(movF$tp[movF$ti[BPs[[6]]]]) #BP2l
  mpar[3] <-  1/temp_2[[MNLL_2_i]]$l_a  #b2
  mpar[5] <- temp_2[[MNLL_1_i]]$k
  
  # Calculate AIC and compare models
  resCA[5] <- 2*(resCA[2] + 1) # Only lambda to estimate
  resCA[6] <- 2*(resCA[3] + 3) # lambda, k and BP
  resCA[7] <- 2*(resCA[4] + 4) # lambda, k and 2*BPs
  resCA[8] <- which.min(resCA[,5:7])-1
  resCA[9] <- 2*resCA[2] + 1*log(resCA[1])
  resCA[10] <- 2*resCA[3] + 3*log(resCA[1])
  resCA[11] <- 2*resCA[4] + 4*log(resCA[1])
  resCA[12] <- which.min(resCA[,9:11])-1
  
  return(list(resCA=resCA,BPs=BPs,mpar=mpar))
}