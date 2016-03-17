require(ggplot2)
require(gridExtra)
require(reshape2)

source('fmm_multivariate_conditional_collapsed.r')

dimension = 1

config = list(
    k = 4
  , m = dimension
  , a = 1
  , l = rep(0, dimension)
  , b = dimension
  , w = diag(1, dimension)
  , n = 50
  )

set.seed(222)

d = list(
    as.matrix(replicate(100, rnorm(config$m, 5)))
  , as.matrix(replicate(100, rnorm(config$m, -5)))
  , as.matrix(replicate(100, rnorm(config$m, 10)))
  , as.matrix(replicate(200, rnorm(config$m))))

dn = lapply(d, function(j) { data.frame(x = j) })

m  = melt(dn, id.vars = c('x'))

set.seed(990909)

params = inverse_model(
    config$n, config$k, as.matrix(m[, c('x')])
  , config$a
  , config$l
  , config$b, config$w
  )

early = data.frame(x = m$x, variable = params$z[1,])
mid   = data.frame(x = m$x, variable = params$z[round(config$n * 1 / 2),])
late  = data.frame(x = m$x, variable = params$z[config$n - 1,])

p_early =
  ggplot(early, aes(x, colour = factor(variable), fill = factor(variable))) +
    geom_histogram(alpha = 0.5)

p_mid =
  ggplot(mid, aes(x, colour = factor(variable), fill = factor(variable))) +
    geom_histogram(alpha = 0.5)

p_late =
  ggplot(late, aes(x, colour = factor(variable), fill = factor(variable))) +
    geom_histogram(alpha = 0.5)

inferred_plots = grid.arrange(p_early, p_mid, p_late, ncol = 3)

