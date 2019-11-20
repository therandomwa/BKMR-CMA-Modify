source("source_BKMR_CMA.R")

TERiskSummaries.CMA <- function(fit.TE,
                                e.y=NULL,
                                e.y.name=NULL,
                                qs = seq(0.25, 0.75, by = 0.05), 
                                q.fixed = 0.5, 
                                alpha = 0.05, 
                                sel,
                                seed = 122) {
  
  toreturn <- data.frame(quantile=qs, 
                         est=rep(NA,times=length(qs)), 
                         sd=rep(NA,times=length(qs)))
  fit <- fit.TE
  Z <- fit$Z
  if (!is.null(e.y.name)){
    Z = Z[,-which(e.y.name == colnames(Z))]
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
  toreturn  
}

CDERiskSummaries.CMA <- function(fit.y,
                                e.y=NULL,
                                e.y.name=NULL,
                                m.value = NULL,
                                m.quant = c(0.1, 0.5, 0.75), 
                                m.name,
                                qs = seq(0.25, 0.75, by = 0.05), 
                                q.fixed = 0.5, 
                                alpha = 0.05, 
                                sel,
                                seed = 122) {
  if (!is.null(m.value) & !is.null(m.quant)){
    m.quant = NULL # if both m.value and m.quant are specified, default set to m.value
  }
  df <- dplyr::tibble()
  fit <- fit.y
  Z <- fit$Z
  m = Z[,m.name]
  if (!is.null(e.y.name)){
    Z = Z[,-which(e.y.name == colnames(Z))]
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
  df
}
