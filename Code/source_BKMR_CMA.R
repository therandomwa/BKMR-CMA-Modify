##########################################################
####      Posterior/Bootstrap Summary Function        ####
##########################################################

##### function to calculate mean/median/CI/sd of a vector of posterior/bootstrap samples
#' Calculate mean/median/CI/sd of a vector of posterior/bootstrap samples
#' 
#' @param posteriorsamp sample prediction
#' @param alpha 1-confidence interval
#' @return Mean. median , CI and sd of the sample prediction
postresults <- function(posteriorsamp, alpha){
  toreturn <- vector()
  toreturn["mean"]   <- mean(posteriorsamp, na.rm=TRUE)
  toreturn["sd"]     <- sd(posteriorsamp, na.rm=TRUE)
  toreturn["lower"]  <- quantile(posteriorsamp, probs=alpha/2, na.rm=TRUE)
  toreturn["median"]  <- quantile(posteriorsamp, probs=0.5, na.rm=TRUE)
  toreturn["upper"] <- quantile(posteriorsamp, probs=1-alpha/2, na.rm=TRUE)
  return(toreturn)
}


##########################################################
####             Estimate TE for BKMR                 ####
##########################################################

#' Get Ya and Yastar used in effects caculation
#' 
#' @param a exposure variables at current level
#' @param astar exposure variables at counterfactual level
#' @param e.y effect modifier for the outcome variable
#' @param fit.y.TE total effect model fit regressing outcome on exposures, effect modifiers and confounders on outcome
#' @param X.predict.Y counfounders for outcome
#' @param sel a vector selecting which iterations of the fit should be retained or inference
#' @param seed the random seed to use to evaluate the code
#' @return A list containing the sample prediction for Ya and Yastar
YaYastar.SamplePred <- function(a, astar, e.y, fit.y.TE, X.predict.Y, sel, seed){
  set.seed(seed)
  z.y = c(a, e.y)
  zstar.y = c(astar, e.y)
  newz <- rbind(z.y, zstar.y)
  
  # give prediction Y for both a and astar
  TE.mat <- SamplePred(fit.y.TE, Znew = newz, Xnew = X.predict.Y, sel = sel)
  Ya <- TE.mat[,"znew1"]
  Yastar <- TE.mat[,"znew2"]
  
  return(list(Ya = Ya, Yastar = Yastar))
}

# a <- apply(A, 2, quantile, probs=0.75)
# astar <- apply(A, 2, quantile, probs=0.25)
# e.y <- quantile(E.Y, probs=0.1)
# TE.age10 <- TE.bkmr(a=a, astar=astar, e.y=e.y, fit.y.TE=fit.y.TE, X.predict.Y=X.predict, alpha = 0.01, sel=sel, seed=122)

#' Estimate total effect for BKMR
#' 
#' @param a exposure variables at current level
#' @param astar exposure variables at counterfactual level
#' @param e.y effect modifier for the outcome variable
#' @param fit.y.TE total effect model fit regressing outcome on exposures, effect modifiers and confounders on outcome
#' @param X.predict.Y counfounders for outcome
#' @param alpha 1-confidence interval
#' @param sel a vector selecting which iterations of the fit should be retained or inference
#' @param seed the random seed to use to evaluate the code
#' @return Totak effect for BKMR
TE.bkmr <- function(a, astar, e.y, fit.y.TE, X.predict.Y, alpha=0.05, sel, seed){

  toreturn <- list()
  
  YaYastar <- YaYastar.SamplePred(a, astar, e.y, fit.y.TE, X.predict.Y, sel, seed)
  Ya     <- YaYastar$Ya
  Yastar <- YaYastar$Yastar
  
  toreturn$TE.samp     <- as.vector(Ya - Yastar)
  toreturn$Ya.samp     <- as.vector(Ya) 
  toreturn$Yastar.samp <- as.vector(Yastar) 
  
  TE.sum <- postresults(toreturn$TE, alpha=alpha)
  toreturn$est  <- TE.sum[c("mean","median","lower","upper")]
  return(toreturn)
}



##########################################################
####           Estimate CDE for BKMR                  ####
##########################################################
# CDE.age10 <- CDE.bkmr(a=a, astar=astar, e.y=e.y, m.quant=c(0.1,0.5,0.75), fit.y=fit.y, sel=sel, seed=777)


#' Estimate controlled direct effect for BKMR
#' 
#' @param a exposure variables at current level
#' @param astar exposure variables at counterfactual level
#' @param e.y effect modifier for the outcome variable
#' @param m.quant values of the quantile that the mediator is set to 
#' @param fit.y model fit regressing outcome on exposures, effect modifiers, mediator and confounders on outcome
#' @param alpha 1-confidence interval
#' @param sel a vector selecting which iterations of the fit should be retained or inference
#' @param seed the random seed to use to evaluate the code
#' @return Controlled direct effect for BKMR
CDE.bkmr <- function(a, astar, e.y, m.quant, fit.y, alpha=0.05, sel, seed){
  
  toreturn <- list()
  
  m <- fit.y$Z[,ncol(fit.y$Z)]  ### okay as long as m is the LAST variable in Zm birthlength
  Z <- fit.y$Z[,-ncol(fit.y$Z)] # exposure + effect modifier
  X.predict<-rep(0,ncol(fit.y$X)) # in the calculation for CDE this value doesn't matter since it will get cancelled out
  
  toreturn$est <- matrix(NA, nrow=length(m.quant), ncol=4, dimnames=list(paste0("CDE",m.quant*100), c("mean","median","lower","upper")))
  for(i in seq_along(m.quant)){
    print(paste(i, "out of", length(m.quant)))
    mnew <- quantile(m, probs=m.quant[i]) # mediator set to certain quantile 
    
    set.seed(seed)
    z.y = c(a, e.y)
    zstar.y = c(astar, e.y)
    newZm  <- rbind(c(z.y,mnew),c(zstar.y,mnew)) # exposure(3), 10% age, 10% mediator 
    CDE.mat <- SamplePred(fit.y, Znew = newZm, Xnew = X.predict, sel=sel) # X.predict is mean of confounders
    
    Yam     <- CDE.mat[,"znew1"]
    Yastarm <- CDE.mat[,"znew2"]
    
    CDE <- as.vector(Yam - Yastarm)
    toreturn[[paste0("CDE",m.quant[i]*100,".samp")]] <- CDE
    toreturn$est[paste0("CDE",m.quant[i]*100),] <- postresults(CDE, alpha=alpha)[c("mean","median","lower","upper")]
  }
  
  return(toreturn)
}



  
##########################################################
####           YaMastar counterfactual                ####
##########################################################

#' Get YaMastar used in NDE and NIE caculation
#' 
#' @param a exposure variables and effect modifier for outcome at current level
#' @param astar.m exposure variables and effect modifier for mediator at counterfactual level 
#' @param fit.m model fit regressing mediator on exposures and confounders on mediator
#' @param fit.y model fit regressing outcome on exposures, effect modifiers, mediator and confounders on outcome
#' @param X.predict.M counfounders for mediator
#' @param X.predict.Y counfounders for outcome
#' @param sel a vector selecting which iterations of the fit should be retained or inference
#' @param seed the random seed to use to evaluate the code
#' @param K number of samples to generate for each MCMC iteration
#' @return A vector containing the sample prediction for YaMastar
YaMastar.SamplePred <- function(a, astar, e.y, fit.m, fit.y, X.predict.M, X.predict.Y, sel, seed, K){
  start.time <- proc.time()
 
  set.seed(seed)
  EM.samp <- SamplePred(fit.m, Znew = astar.m, Xnew = X.predict.M, sel=sel) 
  Mastar     <- as.vector(EM.samp)
  
  sigma.samp  <- sqrt(fit.m$sigsq.eps[sel])
  random.samp <- matrix(rnorm(length(sel)*K),nrow=length(sel),ncol=K)
  
  Mastar.samp  <- Mastar + sigma.samp*random.samp
  
  YaMastar.samp.mat     <- matrix(NA,nrow=length(sel),ncol=K)
  z.y = c(a, e.y)
  for(j in 1:length(sel)){
    Mastar.j <-  Mastar.samp[j,]
    aMastar.j <- cbind(matrix(z.y, nrow=K, ncol=length(z.y), byrow=TRUE), Mastar.j)
    YaMastar.j <- SamplePred(fit.y, Znew = aMastar.j, Xnew = X.predict.Y, sel=sel[j])
    YaMastar.samp.mat[j,] <- as.vector(YaMastar.j)
    
    end.time.temp <- proc.time() 
    if(j%%50==0) print(paste("iter", j, "time: ", round((end.time.temp - start.time)["elapsed"]/60,2),"min")) 
  }
  toreturn <- apply(YaMastar.samp.mat,1,mean)
  
  return(toreturn)
}

YaMastar.SamplePred(a=a, astar=astar.M, e.y = e.y, fit.m=fit.m, fit.y=fit.y,
                    X.predict.M=X.predict.M, X.predict.Y=X.predict.Y, sel=sel, seed=seed, K=K)

##########################################################
####      Estimate NDE/NIE for BKMR(plus TE)          ####
##########################################################

# a <- apply(A, 2, quantile, probs=0.75)
# astar <- apply(A, 2, quantile, probs=0.25)
# e.y <- quantile(E.Y, probs=0.1)
# TE.age10 <- TE.bkmr(a=a, astar=astar, e.y=e.y, fit.y.TE=fit.y.TE, X.predict.Y=X.predict, alpha = 0.01, sel=sel, seed=122)

#' Estimate controlled direct effect for BKMR
#' 
#' @param a.Y exposure variables and effect modifier for outcome at current level
#' @param astar.Y exposure variables and effect modifier for outcome at counterfactual level
#' @param astar.M exposure variables and effect modifier for mediator at counterfactual level 
#' @param m.quant values of the quantile that the mediator is set to 
#' @param fit.m model fit regressing mediator on exposures and confounders on mediator
#' @param fit.y model fit regressing outcome on exposures, effect modifiers, mediator and confounders on outcome
#' @param fit.y.TE total effect model fit regressing outcome on exposures, effect modifiers and confounders on outcome
#' @param X.predict.M counfounders for mediator
#' @param X.predict.Y counfounders for outcome
#' @param alpha 1-confidence interval
#' @param sel a vector selecting which iterations of the fit should be retained or inference
#' @param seed the random seed to use to evaluate the code
#' @param K number of samples to generate for each MCMC iteration in YaMastar calculation
#' @return A list contaning the sample prediction for TE, NDE, NIE and their summary statistics
mediation.bkmr <- function(a.Y, astar.Y, e.y, astar.M, fit.m, fit.y, fit.y.TE, X.predict.M, X.predict.Y, alpha = 0.05, sel, seed, K){

  toreturn <- list()

  # there should be a simpler way of running
  TE <- TE.bkmr(a=a.Y, astar=astar.Y, e.y=e.y, fit.y.TE=fit.y.TE, X.predict=X.predict.Y, alpha=alpha, sel=sel, seed=(seed+100))

  Ya     <- TE$Ya.samp
  Yastar <- TE$Yastar.samp


  YaMastar <- YaMastar.SamplePred(a=a.Y, astar=astar.M, e.y = e.y, fit.m=fit.m, fit.y=fit.y,
                                        X.predict.M=X.predict.M, X.predict.Y=X.predict.Y, sel=sel, seed=seed, K=K)
  
  NDE <- YaMastar - Yastar
  NIE <- Ya - YaMastar


  toreturn$TE.samp <- TE$TE.samp
  toreturn$NDE.samp <- NDE
  toreturn$NIE.samp <- NIE

  toreturn$est <- matrix(NA, nrow=3, ncol=4, dimnames=list(c("TE","NDE","NIE"), c("mean","median","lower","upper")))
  toreturn$est[c("TE","NDE","NIE"),] <- rbind(postresults(TE$TE.samp, alpha=alpha) [c("mean","median","lower","upper")],
                                              postresults(NDE, alpha=alpha)[c("mean","median","lower","upper")],
                                              postresults(NIE, alpha=alpha)[c("mean","median","lower","upper")])

  return(toreturn)
}

