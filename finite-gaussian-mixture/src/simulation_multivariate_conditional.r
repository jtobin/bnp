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

set.seed(222)

d = list(
    t(replicate(250, rnorm(2, c(5, 5))))
  , t(replicate(250, rnorm(2, c(-5, -5))))
  , t(replicate(500, rnorm(2))))
dn = lapply(d, function(j) { data.frame(x = j[,1], y = j[,2]) })
m  = melt(dn, id.vars = c('x', 'y'))

set.seed(990909)

params = inverse_model(
    config$n, config$k, m[, c('x', 'y')]
  , config$a
  , config$l, config$r
  , config$b, config$w
  )

#m_ts_plot = function(j) {
#  melted = as.data.frame(j, id.vars = c('V1', 'V2'))
#  ggplot(
#      melted
#    , aes(x = V1, y = V2, alpha = seq_along(V1), colour = seq_along(V1))) +
#    geom_line() + xlim(-2, 2) + ylim(-10, 10)
#}


