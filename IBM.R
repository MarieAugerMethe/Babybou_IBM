##
# Date created: September 25th, 2013
# Authors: Marie Auger-Methe, Ulrike Schlaegel, Craig DeMars
# Please cite our paper if you use our script:
#   DeMars, C., M. Auger-Méthé, U. Schlägel, S. Boutin, (Published online) 
#   Inferring Parturition and Neonate Survival from Movement Patterns of Female Ungulates.
#   Ecology and Evolution. DOI: 10.1002/ece3.785
# For an in-depth explanation of the code see:
#   Appendix S2. Likelihood Functions Used in the Individual-based Method. 
#   Appendix S3. R Code for the Individual-based Method.
# from the supporting information:
#   ece3785-sup-0001-AppendixS1-S4.docx

#################################################################
# Functions that calculate the negative log likelihood (nll) of 2 of the 3 models
# and the nll of k, one of the parameter used to calculate the nll of these models.
# Note that all models assume 
# that the step length are distributed with an exponential distribution.
# They differ in how the main parameter of the exponential distribution
# varies through the time series.

# nllk calculates the nll of the section of the time series
# when the female is with a calf
# in both the models where the calf survived and the one where the calf died.
# It is used in part to estimate k.
# k represents the time it takes, in number of steps,
# for the female to recover her normal movement.
# nllk is minimised within both nll1 and nll2.
# The function assumes a linear increase in movement after the birth of a calf.
# The slope of this linear increase is ba/k,
# where ba represents the mean step length of the prior birth movement of a female.

nllk <- function(k, SLb, tib, ba, tiBP1){
  bb <- (tib - tiBP1) * ba/k
  bb[ bb > ba ] <- ba
  nllb <- -sum(dexp(SLb, 1/bb, log=TRUE))
  return(nllb)
}

# nll1 calculates the nll of a female with a calf that survived.
# This model assumes that the mean step length is constant prior to the birth of a calf.
# Once the calf is born the mean step length decrease to ba/k, 
# which should be really close to 0.
# The mean step length then linearly increases 
# to reach the same mean step length as prior to calf birth.
# The speed of the increase depends on the time it takes to recover the normal movement,
# which is represented by k.
# Once the mean step length reaches the value from before the birth of the calf,
# it remains at this value.
# The parameters to estimate are:
# 1 BP, ba, and k

nll1 <- function(BP, SL, ti, kc){
  
  ###
  # Divides the time series into two sections:
  # a: before the birth of the calf
  # b: after the birth of the calf
  SLa <- SL[1:BP]
  n <- length(SL)
  SLb <- SL[(BP+1):n] 
  tib <- ti[(BP+1):n]
  
  # maximum likelihood estimate (MLE) of prior birth mean step length
  ba <- mean(SLa)
  # Probability of the SL before the birth of the calf when ba is used
  mnlla <- -sum(dexp(SLa, 1/ba, log=TRUE))
  
  # Numericaly estimating the MLE of k and the probability of SL using this k
  mnllb <- optimize(nllk, kc, SLb=SLb, tib=tib, ba=ba, tiBP1=ti[BP])$objective
  
  # Getting the total nll for the whole time series
  nll <- mnlla + mnllb
  return(nll)
}

# nll2 calculates the nll of a female with a calf that dies.
# This model assumes that the mean step length is constant prior to the birth of a calf.
# Once the calf is born the mean step length decrease to ba/k, 
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
# 2 BPs, ba, and k

nll2 <- function(BP, SL, ti, kc){
  
  # Divides the time series into three sections:
  # a: before the birth of the calf
  # b: after the birth of the calf but before it dies
  # c: after the death of the calf
  SLa <- SL[1:BP[1]]
  SLb <- SL[(BP[1]+1):BP[2]]
  tib <- ti[(BP[1]+1):BP[2]]
  SLc <- SL[(BP[2]+1):length(SL)]
  
  # maximum likelihood estimate (MLE) of prior birth mean step length
  ba <- mean(SLa)
  
  # Probability of the SL before the birth and after the death of the calf
  # when l_a is used
  mnlla <- -sum(dexp(SLa, 1/ba, log=TRUE))
  mnllc <- -sum(dexp(SLc, 1/ba, log=TRUE))
  
  # Numericaly estimating the MLE of k and the probability of SL using this k
  mnllb <- optimize(nllk, kc, SLb=SLb, tib=tib, ba=ba, tiBP1=ti[BP[1]])$objective
  
  # Getting the total nll for the whole time series
  nll <- sum(mnlla, mnllc, mnllb)
  return(nll)
}

# Function minimise the likelihood of the 3 Models

mnll3M <- function(SL, ti, tp, int, kcons){
  # SL: numeric vector that contains the step lengths
  #     measured at regular time intervals.
  #     NAs are not allowed.
  # ti: integer vector that identifies
  #     the time of each step length of SL.
  #     It is the index of the times found in tp.
  #     ti and SL should be of exactly the same length.
  #     NAs are not allowed.
  # tp: POSIXct vector that identifies the real date and time of the SL.
  #     The missing steps should be represented with NAs.
  #     For example if you have locations
  #     (0,0) (0,1) (0,3) (1,3) (7,3)
  #     taken at:
  #     01:00, 02:00, 03:00, 08:00, 09:00
  #     SL = (1,2,6)
  #     ti = (1,2,8)
  #     tp = (01:00, 02:00, NA, NA, NA, NA, NA, 08:00, NA)
  #     (Although excluded for clarity, tp should include the date)
  #     We recommend that the time series only included the time period
  #     relevant to the birth and death of the calf.
  #     We only included movement from April 15th to June 30th.
  
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
  resCA <- matrix(NA, 1, ncol=8)
  colnames(resCA) <- c("n", "mnll_0","mnll_1","mnll_2",
                       "AIC0", "AIC1", "AIC2", "BM")
  
  # BPs index and actual date and time
  BPs <- data.frame(matrix(NA, 1, ncol=6))
  colnames(BPs) <- c("BP1c", "BP2c", "BP2l", "iBP1c", "iBP2c", "iBP2l")
  
  # Parameters (other than the BP) of each model
  mpar <- matrix(NA, 1, 5)
  colnames(mpar) <- c("b0", "b1", "b2", "k1", "k2")
  
  # Sample size
  resCA[1] <- length(SL)
  
  ##
  # M0: No calf
  # This model assumes that the movement pattern is constant for the whole time series
  # and thus the mean step length (b0) is constant.
  # The only parameter estimated is b0, which is the inverse of the rate.
  # It has a analytical solution and 
  # thus we can easily get the minimum negative log likelihood (mnll)
  mpar[1] <- mean(SL)  #b0
  resCA[2] <- -sum(dexp(SL, 1/mpar[1], log=TRUE)) #mnll0
  
  ##
  # Calf survived
  BP1ser <- int:(resCA[1]-int) # All possible BP1c
  NLL1 <- lapply(BP1ser, nll1, SL=SL, ti=ti, kc=kcons)
  MNLL1i <- which.min(NLL1) 
  resCA[3] <- NLL1[[MNLL1i]] # mnll1
  BPs[4] <- BP1ser[MNLL1i] # mle of BP1c in terms of index of SL
  BPs[1] <- as.character(tp[ti[BPs[[4]]]]) #mle of BP1c in real date and time
  mpar[2] <- mean(SL[1:BPs[[4]]]) #b1
  mpar[4] <- optimize(nllk, kcons, SLb=SL[(BPs[[4]]+1):resCA[1]], 
                      tib=ti[(BPs[[4]]+1):resCA[1]],
                      ba=mpar[2], tiBP1=ti[BPs[[4]]])$minimum #k1
  
  ##
  # Calf lost
  # Getting all possible combination of BPs
  # Note that BP are constrained to be int number of non-missing steps apart.
  # To make the code run faster, the BPs are also limited to be less than
  # maximum number of steps it takes for the female to recover her movement apart
  BP2ser <- combn(int:(resCA[1]-int), 2)
  BP2ser <- BP2ser[,diff(BP2ser) >= int]
  BP2ser <- BP2ser[,diff(BP2ser) <= kcons[2]]
  BP2ser <- split(t(BP2ser),1:ncol(BP2ser))
  
  # Applying the model 2 on all possible BPs
  NLL2 <- lapply(BP2ser, nll2,SL=SL, ti=ti, kc=kcons)
  MNLL2i <- which.min(NLL2)
  resCA[4] <- NLL2[[MNLL2i]] # mnll2
  BPs[5:6] <- BP2ser[[MNLL2i]] #mle of iBP2c and iBP2l in terms of index of SL
  BPs[2] <- as.character(tp[ti[BPs[[5]]]]) #mle BP2c in real date and time
  BPs[3] <- as.character(tp[ti[BPs[[6]]]]) #mle BP2l in real date and time
  mpar[3] <- mean(SL[1:BPs[[5]]]) #b2
  mpar[5] <- optimize(nllk, kcons,
                      SLb=SL[(BPs[[5]]+1):BPs[[6]]], 
                      tib=ti[(BPs[[5]]+1):BPs[[6]]],
                      ba=mpar[3], tiBP1=ti[BPs[[5]]])$minimum #k2
  
  # Calculate AIC and compare models
  resCA[5] <- 2*(resCA[2] + 1) # Only b0 to estimate
  resCA[6] <- 2*(resCA[3] + 3) # b1, k1 and BP1c
  resCA[7] <- 2*(resCA[4] + 4) # b2, k2, BP2c, and BP2l
  resCA[8] <- which.min(resCA[,5:7])-1
  
  return(list(resCA=resCA,BPs=BPs,mpar=mpar))
}