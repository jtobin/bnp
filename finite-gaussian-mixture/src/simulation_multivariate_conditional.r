set.seed(222)

require(ggplot2)
require(gridExtra)
require(reshape2)

source('fmm_multivariate_conditional.r')

dimension = 2

config = list(
    k = 3
  , m = dimension
  , a = 1
  , l = rep(0, dimension)
  , r = diag(0.05, dimension)
  , b = 2
  , w = diag(1, dimension)
  , n = 1000
  )

origin = list(
    p = mixing_model(config$k, config$a)
  , m = location_model(config$k, config$l, config$r)
  , s = precision_model(config$k, config$b, config$w)
  )

d = melt(model(config$m, config$k, config$n), id.vars = c('x', 'y'))

set.seed(990909)

params = inverse_model(
    config$n, config$k, d[, c('x', 'y')]
  , config$a
  , config$l, config$r
  , config$b, config$w
  )


m_ts_plot = function(j) {
  melted = as.data.frame(j, id.vars = c('V1', 'V2'))
  ggplot(
      melted
    , aes(x = V1, y = V2, alpha = seq_along(V1), colour = seq_along(V1))) +
    geom_line() + xlim(-2, 2) + ylim(-10, 10)
}


