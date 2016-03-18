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
  , w = diag(0.05, dimension)
  , n = 1000
  )

set.seed(222)

d = list(
    t(replicate(250, rnorm(dimension, c(5, 5))))
  , t(replicate(250, rnorm(dimension, c(-5, -5))))
  , t(replicate(500, rnorm(dimension))))
dn = lapply(d, function(j) { data.frame(x = j[,1], y = j[,2]) })
m  = melt(dn, id.vars = c('x', 'y'))

set.seed(222)

params = inverse_model(
    config$n, config$k, m[, c('x', 'y')]
  , config$a
  , config$l, config$r
  , config$b, config$w
  )

dp = melt(data.frame(params$p))
dm = melt(lapply(params$m, data.frame), id.vars = c('x', 'y'))
dl = melt(as.data.frame(params$l))

py = ggplot(m, aes(x, y)) + geom_point()

pp = ggplot(dp, aes(seq_along(value), value, colour = variable)) +
       geom_line() + facet_grid(. ~ variable)

pm = ggplot(dm, aes(x, y, colour = factor(L1), fill = factor(L1))) +
       geom_point(alpha = 0.5)

pl = ggplot(dl, aes(x = seq_along(value), y = value)) +
       geom_line(colour = 'darkblue')

early = data.frame(x = m$x, y = m$y, variable = params$z[1,])
mid   = data.frame(x = m$x, y = m$y, variable = params$z[round(config$n / 2),])
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

mean_convergence_plots =
  ggplot(dm, aes(x, y, colour = factor(L1), fill = factor(L1))) +
    geom_point(alpha = 0.2) + facet_grid(. ~ L1)

chain_plots = grid.arrange(pp, mean_convergence_plots, nrow = 2)

inferred_plots = grid.arrange(py, p_early, p_mid, p_late, nrow = 2, ncol = 2)

