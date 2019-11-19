source("source_BKMR_CMA.R")

TERiskSummaries.CMA <- function(fit.TE, 
                                    e.y, 
                                    e.y.name,
                                    qs = seq(0.25, 0.75, by = 0.05), 
                                    q.fixed = 0.5, 
                                    sel) {
  
  toreturn <- data.frame(quantile=qs, est=rep(NA,times=length(qs)), sd=rep(NA,times=length(qs)))
  
    for(i in 1:length(qs)){
      quant <- qs[i]
      fit <- fit.TE
      Z <- fit$Z
      Z = Z[,-which(e.y.name == colnames(Z))]
      X <- fit$X
      X.predict <- matrix(colMeans(X),nrow=1)
      
      astar <- apply(Z, 2, quantile, q.fixed)
      a <- apply(Z, 2, quantile, quant)
      
      preds = TE.bkmr(a=a, astar=astar, e.y=e.y, fit.y.TE=BKMRfits, X.predict.Y=X.predict, alpha = 0.01, sel=sel, seed=122)
      toreturn[i, c(2,3)] = c(preds$est["mean"], preds$est["sd"])
      }
  toreturn  
}

# CDERiskSummaries <- function(fit.y, e.y, e.y.name, 
#                              qs = seq(0.25, 0.75, by = 0.05), q.fixed = 0.5, sel)
