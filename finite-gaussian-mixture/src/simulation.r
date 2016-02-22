set.seed(42)

require(ggplot2)
require(reshape2)

source('fmm_generative.r')

config = list(
    k = 4
  , n = 10000
  )

d      = model(config$k, config$n)
melted = melt(d)

p = ggplot(melted, aes(value, colour = factor(L1), fill = factor(L1))) +
      geom_density(alpha = 0.5)

