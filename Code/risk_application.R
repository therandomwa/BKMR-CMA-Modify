library(bkmr)

load("../Data/test_data.RData")
source("source_BKMR_CMA.R")
source("risk.R")

sel = seq(5001,20000,by=1000)

# overall effects
riskSummary10 = TERiskSummaries.CMA(fit.TE = fit.y.TE, e.y=e.y10, e.y.name = "age", sel=sel)
ggplot(riskSummary10,
       aes(quantile,
           est,
           ymin = est - 1.96 * sd,
           ymax = est + 1.96 * sd)) +
  geom_pointrange()


riskSummary90 = TERiskSummaries.CMA(fit.TE = fit.y.TE, e.y=e.y90, e.y.name = "age", sel=sel)
ggplot(riskSummary90,
       aes(quantile,
           est,
           ymin = est - 1.96 * sd,
           ymax = est + 1.96 * sd)) +
  geom_pointrange()



# CDE 
CDEriskSummary10 = CDERiskSummaries.CMA(fit.y = fit.y, e.y = e.y10, e.y.name = "age", m.name = "BL", sel = sel)
ggplot(CDEriskSummary10, aes(quantile, est, ymin = est - 1.96*sd, 
                             ymax = est + 1.96*sd, col = m)) + 
  geom_pointrange(position = position_dodge(width = 0.75))


CDEriskSummary90 = CDERiskSummaries.CMA(fit.y = fit.y, e.y = e.y90, e.y.name = "age", m.name = "BL", sel = sel)
ggplot(CDEriskSummary90, aes(quantile, est, ymin = est - 1.96*sd, 
                             ymax = est + 1.96*sd, col = m)) + 
  geom_pointrange(position = position_dodge(width = 0.75))



# single variable total effects
risks.singvar10 = SingVarRiskSummaries.CMA(BKMRfits = fit.y.TE,
                                           e.y=e.y10, e.y.names="age",
                                           sel=sel)
ggplot(risks.singvar10, aes(variable, est, ymin = est - 1.96*sd, 
                            ymax = est + 1.96*sd, col = q.fixed)) + 
  geom_pointrange(position = position_dodge(width = 0.75)) + 
  coord_flip()


risks.singvar90 = SingVarRiskSummaries.CMA(BKMRfits = fit.y.TE,
                                           e.y=e.y90, e.y.names="age",
                                           sel=sel)
ggplot(risks.singvar90, aes(variable, est, ymin = est - 1.96*sd, 
                            ymax = est + 1.96*sd, col = q.fixed)) + 
  geom_pointrange(position = position_dodge(width = 0.75)) + 
  coord_flip()
