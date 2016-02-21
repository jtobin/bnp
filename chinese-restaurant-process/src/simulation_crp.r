require(dplyr)
require(ggplot2)
require(reshape2)

source('crp.r')

design = expand.grid(epochs = 100, n = 1000, a = c(1, 10, 100))

simulate = function(epochs, n, a) replicate(epochs, list(crp(n, a)))

experiment = apply(
    design
  , MARGIN = 1
  , function(row) { simulate(row[1], row[2], row[3]) }
  )

results          = melt(experiment, id.vars = 'table')
by_concentration = group_by(results, table = table, settings = factor(L1))

log_log_plot =
  ggplot(by_concentration, aes(table, value, fill = settings, colour = settings)) +
  geom_jitter(width = 0.5, height = 0.45, alpha = 0.2) +
  scale_x_log10() + scale_y_log10()

summarised = summarise(by_concentration, customers = mean(value))

mean_log_log_plot =
  ggplot(summarised, aes(table, customers, fill = settings, colour = settings)) +
  geom_point(alpha = 0.8) +
  scale_x_log10() + scale_y_log10()

