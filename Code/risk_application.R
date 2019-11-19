load("../Data/test_data.RData")
source("source_BKMR_CMA.R")
source("risk.R")

sel = seq(5001,20000,by=1000)
riskSummary10 = TERiskSummaries.CMA(fit.TE = fit.y.TE, e.y=e.y10, e.y.name = "age", sel=sel)
ggplot(riskSummary10, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) + 
  geom_pointrange()

riskSummary90 = TERiskSummaries.CMA(fit.TE = fit.y.TE, e.y=e.y90, e.y.name = "age", sel=sel)
ggplot(riskSummary90, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) + 
  geom_pointrange()
