source("source_BKMR_CMA.R")

TERiskSummaries.CMA <- function(fit.TE,
                                e.y=NULL, e.y.names=NULL,
                                qs = seq(0.25, 0.75, by = 0.05), q.fixed = 0.5, 
                                alpha = 0.05, sel, seed = 122) {
  
  toreturn <- data.frame(quantile=qs, 
                         est=rep(NA,times=length(qs)), 
                         sd=rep(NA,times=length(qs)))
  fit <- fit.TE
  Z <- fit$Z
  if (!is.null(e.y.names)){
    Z = Z[,-which(e.y.names == colnames(Z))]
  }
  X <- fit$X
  X.predict <- matrix(colMeans(X),nrow=1)
  for(i in 1:length(qs)){
    quant <- qs[i]
    astar <- apply(Z, 2, quantile, q.fixed)
    a <- apply(Z, 2, quantile, quant)
    preds = TE.bkmr(a=a, astar=astar, e.y=e.y, fit.y.TE=fit.TE, 
                    X.predict.Y=X.predict, alpha = alpha, sel=sel, seed=seed)
    toreturn[i, c(2,3)] = c(preds$est["mean"], preds$est["sd"])
  }
  return(toreturn)
}



CDERiskSummaries.CMA <- function(fit.y,
                                 e.y=NULL, e.y.names=NULL,
                                 m.value = NULL, m.quant = c(0.1, 0.5, 0.75), m.name,
                                 qs = seq(0.25, 0.75, by = 0.05), q.fixed = 0.5, 
                                 alpha = 0.05, sel, seed = 122) {
  if (!is.null(m.value) & !is.null(m.quant)){
    m.quant = NULL # if both m.value and m.quant are specified, default set to m.value
  }
  df <- dplyr::tibble()
  fit <- fit.y
  Z <- fit$Z
  m = Z[,m.name]
  if (!is.null(e.y.names)){
    Z = Z[,-which(e.y.names == colnames(Z))]
  }
  Z = Z[,-which(m.name == colnames(Z))]
  X <- fit$X
  X.predict <- matrix(colMeans(X),nrow=1)
  for(i in seq_along(qs)) {
    print(paste("Running", i, "out of", length(qs), "quantile values:"))
    quant <- qs[i]
    astar <- apply(Z, 2, quantile, q.fixed)
    a <- apply(Z, 2, quantile, quant)
    preds = CDE.bkmr(a, astar, e.y, m.value, m.quant, fit.y, alpha, sel, seed)
    if (is.null(m.quant)){newm = m.value} else(newm = m.quant)
    df0 <-dplyr::data_frame(quantile = quant,
                            m = newm,
                            est = preds$est[,"mean"],
                            sd = preds$est[,"sd"])
    df <- dplyr::bind_rows(df, df0)
  }
  df$m = as.factor(df$m)
  df$quantile = as.factor(df$quantile)
  return(df)
}



### ***** a combination function of VarRiskSummary and riskSummary.approx for MI BKMR fits ******
#Compare estimated \code{h} function when a single variable (or a set of variables) is at the 75th versus 25th percentile, when all of the other variables are fixed at a particular percentile
VarRiskSummary.CMA <-function(whichz = 1,
                             BKMRfits,
                             e.y = NULL, e.y.names = NULL, 
                             qs.diff = c(0.25, 0.75), q.fixed = 0.5,
                             alpha, sel, seed) {
  
  fit <- BKMRfits
  y <- fit$y
  Z <- fit$Z
  X <- fit$X
  if (!is.null(e.y.names)){
    Z = Z[,-which(e.y.names == colnames(Z))]
  }
  X.predict <- matrix(colMeans(X),nrow=1)
  
  a <- astar <- apply(Z, 2, quantile, q.fixed)
  astar[whichz] <- apply(Z[, whichz, drop = FALSE], 2, quantile, qs.diff[1]) # fix z-th variable to the quantile
  a[whichz] <- apply(Z[, whichz, drop = FALSE], 2, quantile, qs.diff[2])
  
  preds = TE.bkmr(a=a, astar=astar, e.y=e.y, fit.y.TE=fit, 
                  X.predict.Y=X.predict, alpha = alpha, sel=sel, seed=seed)
  toreturn = c(preds$est["mean"], preds$est["sd"])
  names(toreturn) = c("est", "sd")
  return(toreturn)
}


# Single Variable Risk Summaries
# 
# Compute summaries of the risks associated with a change in a single variable in \code{Z} from a single level (quantile) to a second level (quantile), for the other variables in \code{Z} fixed to a specific level (quantile)
# **** for MI BKMR fits ****

SingVarRiskSummaries.CMA <-function(BKMRfits,
                                    e.y = NULL, e.y.names = NULL,
                                    which.z = 1:ncol(BKMRfits$Z),
                                    z.names = colnames(BKMRfits$Z),
                                    qs.diff = c(0.25, 0.75),
                                    q.fixed = c(0.25, 0.50, 0.75),
                                    alpha = 0.05, sel, seed = 122) {
  Z = BKMRfits$Z
  if(is.null(z.names)) z.names <- paste0("z", 1:ncol(Z))
  if (!is.null(e.y.names)){
    Z = Z[,-which(e.y.names == colnames(Z))]
    which.z = 1:ncol(Z)
    z.names = colnames(Z)
  }
  df <- tibble::tibble()
  for(i in seq_along(q.fixed)) {
    for(j in seq_along(which.z)) {
      risk = VarRiskSummary.CMA(whichz = which.z[j],
                                BKMRfits = BKMRfits,
                                e.y = e.y, e.y.names = e.y.names,
                                qs.diff = qs.diff, q.fixed = q.fixed[i],
                                alpha = alpha, sel = sel, seed = seed)
      df0 <- tibble::tibble(q.fixed = q.fixed[i],
                            variable = z.names[which.z[j]],
                            est = risk["est"],
                            sd = risk["sd"])
      df <- dplyr::bind_rows(df, df0)
    }
  }
  df$variable <- factor(df$variable, levels = z.names[which.z])
  df$q.fixed = as.factor(df$q.fixed)
  attr(df, "qs.diff") <- qs.diff
  return(df)
}


###### Single variable plot for cde #######
### ***** a combination function of VarRiskSummary and riskSummary.approx for MI BKMR fits ******
#Compare estimated \code{h} function when a single variable (or a set of variables) is at the 75th versus 25th percentile, when all of the other variables are fixed at a particular percentile
CDEVarRiskSummary.CMA <-function(whichz = 1,
                              BKMRfits,
                              e.y = NULL, e.y.names = NULL, 
                              m.value = NULL, m.quant = c(0.1, 0.5, 0.75), m.name,
                              qs.diff = c(0.25, 0.75), q.fixed = 0.5,
                              alpha, sel, seed) {
  
  # if (!is.null(m.value) & !is.null(m.quant)){
  #   m.quant = NULL
  # }
  fit <- BKMRfits
  y <- fit$y
  Z <- fit$Z
  # X <- fit$X
  m = Z[,m.name]
  if (!is.null(e.y.names)){
    Z = Z[,-which(e.y.names == colnames(Z))]
  }
  Z = Z[,-which(m.name == colnames(Z))]
  # X.predict <- matrix(colMeans(X),nrow=1)
  
  a <- astar <- apply(Z, 2, quantile, q.fixed)
  astar[whichz] <- apply(Z[, whichz, drop = FALSE], 2, quantile, qs.diff[1]) # fix z-th variable to the quantile
  a[whichz] <- apply(Z[, whichz, drop = FALSE], 2, quantile, qs.diff[2])
  
  preds = CDE.bkmr(a=a, astar=astar, e.y=e.y, m.value = m.value, m.quant = m.quant, 
                   fit.y=fit, alpha = alpha, sel=sel, seed=seed)
  toreturn = preds$est[,c("mean","sd")]
  colnames(toreturn) = c("est", "sd")
  return(toreturn)
}


# Single Variable Risk Summaries
# 
# Compute summaries of the risks associated with a change in a single variable in \code{Z} from a single level (quantile) to a second level (quantile), for the other variables in \code{Z} fixed to a specific level (quantile)
# **** for MI BKMR fits ****
CDESingVarRiskSummaries.CMA <-function(BKMRfits,
                                       e.y = NULL, e.y.names = NULL,
                                       which.z = 1:ncol(BKMRfits$Z),
                                       z.names = colnames(BKMRfits$Z),
                                       m.value = NULL, m.quant = c(0.1, 0.5, 0.75), m.name,
                                       qs.diff = c(0.25, 0.75),
                                       q.fixed = c(0.25, 0.50, 0.75),
                                       alpha = 0.05, sel, seed = 122) {
  if (!is.null(m.value) & !is.null(m.quant)){
    m.quant = NULL
  }
  Z = BKMRfits$Z
  if(is.null(z.names)) z.names <- paste0("z", 1:ncol(Z))
  if (!is.null(e.y.names)){
    Z = Z[,-which(e.y.names == colnames(Z))]
  }
  Z = Z[,-which(m.name == colnames(Z))]
  which.z = 1:ncol(Z)
  z.names = colnames(Z)
  df <- tibble::tibble()
  for(i in seq_along(q.fixed)) {
    for(j in seq_along(which.z)) {
      risk = CDEVarRiskSummary.CMA(whichz = which.z[j],
                                BKMRfits = BKMRfits,
                                e.y = e.y, e.y.names = e.y.names,
                                m.value = m.value, m.quant = m.quant, m.name = m.name,
                                qs.diff = qs.diff, q.fixed = q.fixed[i],
                                alpha = alpha, sel = sel, seed = seed)
      df0 <- tibble::tibble(q.fixed = q.fixed[i],
                            m.fixed = rownames(risk),
                            variable = z.names[which.z[j]],
                            est = risk[,"est"],
                            sd = risk[,"sd"])
      df <- dplyr::bind_rows(df, df0)
    }
  }
  df$m.fixed = gsub("CDE","", df$m.fixed)
  df$m.fixed = as.factor(df$m.fixed)
  df$variable = as.factor(df$variable)
  # df$variable <- factor(df$variable, levels = z.names[which.z])
  df$q.fixed = as.factor(df$q.fixed)
  attr(df, "qs.diff") <- qs.diff
  return(df)
}

###############################
# PredictorResponseBivarPair.MI <-
#   function(fit,
#            y,Z,X,
#            whichz1 = 1,whichz2 = 2,whichz3 = NULL,
#            prob = 0.5,
#            q.fixed = 0.5,
#            sel = NULL,
#            ngrid = 50,
#            min.plot.dist = 0.5,
#            center = TRUE,
#            Z.MI) {
#     
#   if(ncol(Z) < 3) stop("requires there to be at least 3 Z variables")
#   
#   if(is.null(colnames(Z))) colnames(Z) <- paste0("z", 1:ncol(Z))
#   
#   if(is.null(whichz3)) {
#     ord <- c(whichz1, whichz2, setdiff(1:ncol(Z), c(whichz1, whichz2)))
#   } else {
#     ord <- c(whichz1, whichz2, whichz3, setdiff(1:ncol(Z), c(whichz1, whichz2, whichz3)))
#   }
#   z1 <- seq(min(Z[,ord[1]]), max(Z[,ord[1]]), length=ngrid) # As
#   z2 <- seq(min(Z[,ord[2]]), max(Z[,ord[2]]), length=ngrid) # Mn 
#   z3 <- quantile(Z[, ord[3]], probs = prob) # Pb
#   z.all <- c(list(z1), list(z2), list(z3))
#   if(ncol(Z) > 3) {
#     z.others <- lapply(4:ncol(Z), function(x) quantile(Z[,ord[x]], q.fixed))
#     z.all <- c(z.all, z.others) # + age
#   }
#   newz.grid <- expand.grid(z.all) # Fix Pb(prob) and Age(q.fixed) to 50%
#   z1save <- newz.grid[, 1]
#   z2save <- newz.grid[, 2]
#   colnames(newz.grid) <- colnames(Z)[ord]
#   newz.grid <- newz.grid[,colnames(Z)]
#   
#   if(!is.null(min.plot.dist)) {
#     mindists <- rep(NA, nrow(newz.grid))
#     for(k in seq_along(mindists)) {
#       pt <- as.numeric(newz.grid[k,c(colnames(Z)[ord[1]],colnames(Z)[ord[2]])]) # get each row of As and Mn
#       dists <- fields::rdist(matrix(pt, nrow = 1), Z[, c(colnames(Z)[ord[1]],colnames(Z)[ord[2]])])
#       mindists[k] <- min(dists)
#     }
#   }
#   
#   
#     preds <- ComputePostmeanHnew(fit = fit, y = y, Z = Z, X = X, Znew = newz.grid, sel = sel, method = method)
#     preds.plot <- preds$postmean
#     se.plot <- sqrt(diag(preds$postvar))
#   
#   if(center) preds.plot <- preds.plot - mean(preds.plot)
#   if(!is.null(min.plot.dist)) {
#     preds.plot[mindists > min.plot.dist] <- NA
#     se.plot[mindists > min.plot.dist] <- NA
#   }
#   #     hgrid <- matrix(preds.plot, ngrid, ngrid, dimnames=list(z1=round(z1,2), z2=round(z2,2)))
#   #     se.grid <- matrix(se.plot, ngrid, ngrid, dimnames=list(z1=round(z1,2), z2=round(z2,2)))
#   
#   res <- dplyr::data_frame(z1 = z1save, z2 = z2save, est = preds.plot, se = se.plot)
# }

#' Predict the exposure-response function at a new grid of points
#'
#' Predict the exposure-response function at a new grid of points
#'

# PredictorResponseBivar.singfit.MI <- function(fit, y = NULL, Z = NULL, X = NULL, z.pairs = NULL, method = "approx", ngrid = 50, q.fixed = 0.5, sel = NULL, min.plot.dist = 0.5, center = TRUE, z.names = colnames(Z), verbose = TRUE, Z.MI,k=1,K=1, ...) {
#   
#   if (inherits(fit, "bkmrfit")) {
#     if (is.null(y)) y <- fit$y
#     if (is.null(Z)) Z <- fit$Z
#     if (is.null(X)) X <- fit$X
#   }
#   
#   if (is.null(z.names)) {
#     z.names <- colnames(Z.MI)
#     if (is.null(z.names)) {
#       z.names <- paste0("z", 1:ncol(Z))
#     }
#   }
#   
#   if (is.null(z.pairs)) {
#     z.pairs <- expand.grid(z1 = 1:ncol(Z), z2 = 1:ncol(Z))
#     z.pairs <- z.pairs[z.pairs$z1 < z.pairs$z2, ]
#   }
#   
#   df <- dplyr::data_frame()
#   for(i in 1:nrow(z.pairs)) {
#     compute <- TRUE
#     whichz1 <- z.pairs[i, 1] %>% unlist %>% unname
#     whichz2 <- z.pairs[i, 2] %>% unlist %>% unname
#     if(whichz1 == whichz2) compute <- FALSE
#     z.name1 <- z.names[whichz1]
#     z.name2 <- z.names[whichz2]
#     names.pair <- c(z.name1, z.name2)
#     if(nrow(df) > 0) { ## determine whether the current pair of variables has already been done
#       completed.pairs <- df %>%
#         dplyr::select_('variable1', 'variable2') %>%
#         dplyr::distinct() %>%
#         dplyr::transmute(z.pair = paste('variable1', 'variable2', sep = ":")) %>%
#         unlist %>% unname
#       if(paste(names.pair, collapse = ":") %in% completed.pairs | paste(rev(names.pair), collapse = ":") %in% completed.pairs) compute <- FALSE
#     }
#     if(compute) {
#       if(verbose) message("MI fit ", k," out of ",K, ":  Pair ", i, " out of ", nrow(z.pairs))
#       res <- PredictorResponseBivarPair.MI(fit = fit, y = y, Z = Z, X = X, whichz1 = whichz1, whichz2 = whichz2, method = method, ngrid = ngrid, q.fixed = q.fixed, sel = sel, min.plot.dist = min.plot.dist, center = center, z.names = z.names, Z.MI=Z.MI, ...)
#       df0 <- res
#       df0$variable1 <- z.name1
#       df0$variable2 <- z.name2
#       df0 %<>%
#         dplyr::select_(~variable1, ~variable2, ~z1, ~z2, ~est, ~se)
#       df <- dplyr::bind_rows(df, df0)
#     }
#   }
#   df$variable1 <- factor(df$variable1, levels = z.names)
#   df$variable2 <- factor(df$variable2, levels = z.names)
#   df
# }


#' Plot cross-sections of the bivariate predictor-response function
#' 
#' Function to plot the \code{h} function of a particular variable at different levels (quantiles) of a second variable
#' 


# PredictorResponseBivar.MI <- function(BKMRfits, z.pairs = NULL, method = "approx", ngrid = 50, q.fixed = 0.5, sel = NULL, min.plot.dist = 0.5, center = TRUE, z.names = colnames(BKMRfits[[1]]$Z), verbose = TRUE, ...) {
#   
#   start.time <- proc.time()["elapsed"]
#   Z.MI <- Z.complete.MI(BKMRfits)
#   l <- ncol(Z.MI)
#   z.names <- colnames(Z.MI)
#   
#   K <- length(BKMRfits)  
#   
#   npairs <- length(z.pairs)
#   if(is.null(z.pairs)) npairs <- factorial(l)/factorial(l-2)/factorial(2)
#   est.matrix <- matrix(NA,nrow=ngrid*ngrid*npairs, ncol=K)
#   for(k in 1:K){
#     fit <- BKMRfits[[k]]
#     y <- fit$y
#     Z <- fit$Z
#     X <- fit$X	
#     
#     temp <- PredictorResponseBivar.singfit.MI(fit=fit, z.pairs = z.pairs, method = method, ngrid = ngrid, q.fixed = q.fixed, sel = sel, min.plot.dist = min.plot.dist, center = center, z.names = z.names, verbose = TRUE, Z.MI=Z.MI,k=k,K=K, ...)
#     
#     est.matrix[,k] <- temp %>% select(est) %>% unlist(use.names=FALSE)
#     
#     end.time <- proc.time()["elapsed"]
#     message(paste("MI fit", k, "out of", K, "complete: ", round((end.time - start.time)/60, digit=2), "min run time" ))	
#   }
#   
#   est.MI <- apply(est.matrix,1,mean)
#   
#   
#   data.toreturn <- temp ### okay since the grid numbers and variable order are the same
#   data.toreturn[,"est"] <- est.MI 
#   
#   data.toreturn
# }
# 
# 
