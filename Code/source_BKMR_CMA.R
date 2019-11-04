##########################################################
####      Posterior/Bootstrap Summary Function        ####
##########################################################


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


#' Estimate controlled direct effect for BKMR
#' 
#' @param a exposure variables at current level
#' @param astar exposure variables at counterfactual level
#' @param e.y effect modifier for the outcome variable
#' @param m.value values that the mediator is set to
#' @param m.quant values of the quantile that the mediator is set to 
#' @param fit.y model fit regressing outcome on exposures, effect modifiers, mediator and confounders on outcome
#' @param alpha 1-confidence interval
#' @param sel a vector selecting which iterations of the fit should be retained or inference
#' @param seed the random seed to use to evaluate the code
#' @return Controlled direct effect for BKMR
CDE.bkmr <- function(a, astar, e.y, m.value=NULL, m.quant=NULL, fit.y, alpha=0.05, sel, seed){
  
  toreturn <- list()
  m <- fit.y$Z[,ncol(fit.y$Z)]  ### okay as long as m is the LAST variable in Zm birthlength
  Z <- fit.y$Z[,-ncol(fit.y$Z)] # exposure + effect modifier
  X.predict<-rep(0,ncol(fit.y$X)) # in the calculation for CDE this value doesn't matter since it will get cancelled out
  if (is.null(m.value)){
    toreturn$est <- matrix(NA, nrow=length(m.quant), ncol=4, dimnames=list(paste0("CDE",m.quant*100, "%"), c("mean","median","lower","upper")))
    print(paste("Running", length(m.quant), "mediator values for CDE:"))
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
      toreturn[[paste0("CDE",m.quant[i]*100,"%.samp")]] <- CDE
      toreturn$est[paste0("CDE",m.quant[i]*100,"%"),] <- postresults(CDE, alpha=alpha)[c("mean","median","lower","upper")]
    }
  }
  else if (is.null(m.quant)){
    toreturn$est <- matrix(NA, nrow=length(m.value), ncol=4, dimnames=list(paste0("CDE",m.value), c("mean","median","lower","upper")))
    for(i in seq_along(m.value)){
      print(paste(i, "out of", length(m.value)))
      mnew = m.value[i]
      
      set.seed(seed)
      z.y = c(a, e.y)
      zstar.y = c(astar, e.y)
      newZm  <- rbind(c(z.y,mnew),c(zstar.y,mnew)) # exposure(3), 10% age, 10% mediator 
      CDE.mat <- SamplePred(fit.y, Znew = newZm, Xnew = X.predict, sel=sel) # X.predict is mean of confounders
      
      Yam     <- CDE.mat[,"znew1"]
      Yastarm <- CDE.mat[,"znew2"]
      
      CDE <- as.vector(Yam - Yastarm)
      toreturn[[paste0("CDE",mnew,".samp")]] <- CDE
      toreturn$est[paste0("CDE",mnew),] <- postresults(CDE, alpha=alpha)[c("mean","median","lower","upper")]
    }
  }
  
  return(toreturn)
}



  
##########################################################
####           YaMastar counterfactual                ####
##########################################################


#' Get YaMastar used in NDE and NIE caculation
#' 
#' @param a exposure variables at current level
#' @param astar exposure variables at counterfactual level
#' @param e.y effect modifier for the outcome variable
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
  EM.samp <- SamplePred(fit.m, Znew = astar, Xnew = X.predict.M, sel=sel) 
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


##########################################################
####      Estimate NDE/NIE for BKMR(plus TE)          ####
##########################################################


#' Estimate controlled direct effect for BKMR
#' 
#' @param a exposure variables at current level
#' @param astar exposure variables at counterfactual level
#' @param e.y effect modifier for the outcome variable
#' @param fit.m model fit regressing mediator on exposures and confounders on mediator
#' @param fit.y model fit regressing outcome on exposures, effect modifiers, mediator and confounders on outcome
#' @param fit.y.TE total effect model fit regressing outcome on exposures, effect modifiers and confounders on outcome
#' @param X.predict.M counfounders for mediator
#' @param X.predict.Y counfounders for outcome
#' @param effects type(s) of effects that users want to output
#' @param m.value values that the mediator is set to
#' @param m.quant values of the quantile that the mediator is set to 
#' @param alpha 1-confidence interval
#' @param sel a vector selecting which iterations of the fit should be retained or inference
#' @param seed the random seed to use to evaluate the code
#' @param K number of samples to generate for each MCMC iteration in YaMastar calculation
#' @return A list contaning the sample prediction for TE, NDE, NIE and their summary statistics
mediation.bkmr <- function(a, astar, e.y, fit.m=NULL, fit.y=NULL, fit.y.TE=NULL, 
                           X.predict.M=NULL, X.predict.Y=NULL, 
                           effects = "all",  # c("all", "TE", "CDE", "NDE", "NIE")
                           m.quant=NULL, # c(0.1,0.5,0.75), 
                           m.value=NULL,
                           alpha = 0.05, sel, seed, K){
  
  if (sum(!effects %in% c("all", "TE", "CDE", "NDE", "NIE"))) {
    stop("effects must be in c('all', 'TE', 'CDE', 'NDE', 'NIE')")
  }
  else{
    if ("all" %in% effects){
      effects = c("TE", "CDE", "NDE", "NIE")
    }
    if ("TE" %in% effects){
      if (is.null(fit.y.TE)){
        stop("Must specify 'fit.y.TE'")
      }
      if (is.null(X.predict.Y)){
        stop("Must specify 'X.predict.Y'")
      }
    }
    if (sum(c("NIE", "NDE") %in% effects)){
      if (is.null(fit.y.TE) | is.null(fit.m) | is.null(fit.y)){
        stop("Must specify all three fits: 'fit.y.TE', 'fit.y', 'fit.m'")
      }
      if (is.null(X.predict.Y)){
        stop("Must specify 'X.predict.Y'")
      }
      if (is.null(X.predict.M)){
        stop("Must specify 'X.predict.M'")
      }
    }
    if ("CDE" %in% effects){
      if (is.null(fit.y)){
        stop("Must specify 'fit.y'")
      }
      if (is.null(m.value) & is.null(m.quant)){
        stop("Must specify either 'm.value' or 'm.quant'")
      }
      if (!is.null(m.value) & !is.null(m.quant)){
        stop("Must only specify one of 'm.value' and 'm.quant'")
      }
    }
    
    # start code
    TE = NULL; CDE = NULL; NDE = NULL; NIE = NULL;
    toreturn <- list()
    effects.temp = effects
    if ("CDE" %in% effects){
      effects.temp = effects[effects!="CDE"]
    }
    toreturn$est <- matrix(NA, nrow=length(effects.temp), ncol=4, dimnames=list(effects.temp, c("mean","median","lower","upper")))
  
    if ("TE" %in% effects){
      TE <- TE.bkmr(a=a, astar=astar, e.y=e.y, fit.y.TE=fit.y.TE, X.predict=X.predict.Y, alpha=alpha, sel=sel, seed=(seed+100))
      toreturn$TE.samp <- TE$TE.samp
      toreturn$est["TE",] <- postresults(TE$TE.samp, alpha=alpha) [c("mean","median","lower","upper")]
    }
    if (sum(c("NIE", "NDE") %in% effects)){
      if (is.null(TE)){
        TE <- TE.bkmr(a=a, astar=astar, e.y=e.y, fit.y.TE=fit.y.TE, X.predict=X.predict.Y, alpha=alpha, sel=sel, seed=(seed+100))
      }
      
      Ya     <- TE$Ya.samp
      Yastar <- TE$Yastar.samp
      
      YaMastar <- YaMastar.SamplePred(a=a, astar=astar, e.y = e.y, fit.m=fit.m, fit.y=fit.y,
                                      X.predict.M=X.predict.M, X.predict.Y=X.predict.Y, sel=sel, seed=seed, K=K)
      if ("NDE" %in% effects){
        NDE <- YaMastar - Yastar
        toreturn$NDE.samp <- NDE
        toreturn$est["NDE",] <- postresults(NDE, alpha=alpha) [c("mean","median","lower","upper")]
      }
      if ("NIE" %in% effects){
        NIE <- Ya - YaMastar
        toreturn$NIE.samp <- NIE
        toreturn$est["NIE",] <- postresults(NIE, alpha=alpha) [c("mean","median","lower","upper")]
      }
    }
    if ("CDE" %in% effects){
      CDE <- CDE.bkmr(a, astar, e.y, m.value=m.value, m.quant=m.quant, fit.y, alpha, sel, seed)
      toreturn$est = rbind(toreturn$est, CDE$est)
      toreturn = append(toreturn, CDE[2:length(CDE)])
    }
    return (toreturn)
  }
}

