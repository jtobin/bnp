set.seed(42)

require(ggplot2)
require(reshape2)

source('fmm_multivariate_generative.r')

config = list(
    m = 2
  , k = 4
  , n = 10000
  )

d      = model(config$m, config$k, config$n)
melted = melt(d, id.vars = c('x', 'y'))

p = ggplot(melted, aes(x, y, colour = factor(L1))) + geom_point(alpha = 0.2)

