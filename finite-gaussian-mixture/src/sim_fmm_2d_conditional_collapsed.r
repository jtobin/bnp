require(ggplot2)
require(gridExtra)
require(reshape2)

source('fmm_multivariate_conditional_collapsed.r')

dimension = 2

config = list(
    k = 3
  , m = dimension
  , a = 1
  , l = rep(0, dimension)
  , r = diag(0.05, dimension)
  , b = 2
  , w = diag(1, dimension)
  , n = 50
  )

set.seed(222)

d = list(
    t(replicate(100, rnorm(2, c(5, 5))))
  , t(replicate(100, rnorm(2, c(-5, -5))))
  , t(replicate(200, rnorm(2))))
dn = lapply(d, function(j) { data.frame(x = j[,1], y = j[,2]) })
m  = melt(dn, id.vars = c('x', 'y'))

set.seed(990909)

params = inverse_model(
    config$n, config$k, m[, c('x', 'y')]
  , config$a
  , config$l, config$r
  , config$b, config$w
  )

early = data.frame(x = m$x, y = m$y, variable = params$z[1,])
mid   = data.frame(x = m$x, y = m$y, variable = params$z[round(config$n * 1 / 2),])
late  = data.frame(x = m$x, y = m$y, variable = params$z[config$n - 1,])

p_early =
  ggplot(early, aes(x, y, colour = factor(variable), fill = factor(variable))) +
    geom_point(alpha = 0.5)

p_mid =
  ggplot(mid, aes(x, y, colour = factor(variable), fill = factor(variable))) +
    geom_point(alpha = 0.5)

p_late =
  ggplot(late, aes(x, y, value, colour = factor(variable), fill = factor(variable))) +
    geom_point(alpha = 0.5)

inferred_plots = grid.arrange(p_early, p_mid, p_late, ncol = 3)

