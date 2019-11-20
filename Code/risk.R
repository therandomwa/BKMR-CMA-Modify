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
  a[whichz] <- apply(Z[, whichz, drop = FALSE], 2, quantile, qs.diff[1]) # fix z-th variable to the quantile
  astar[whichz] <- apply(Z[, whichz, drop = FALSE], 2, quantile, qs.diff[2])
  
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
  
  if(is.null(z.names)) z.names <- paste0("z", 1:ncol(BKMRfits$Z))
  if (!is.null(e.y.names)){
    Z = BKMRfits$Z
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
