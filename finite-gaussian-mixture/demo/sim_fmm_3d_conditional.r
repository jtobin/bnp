require(ggplot2)
require(gridExtra)
require(reshape2)

source('../src/fmm_multivariate_conditional.r')

dimension = 3

config = list(
    k = 3
  , m = dimension
  , a = 1
  , l = rep(0, dimension)
  , r = diag(0.05, dimension)
  , b = dimension
  , w = diag(1, dimension)
  , n = 5000
  )

set.seed(222)

d = list(
    t(replicate(250, rnorm(config$m, c(5, 5))))
  , t(replicate(250, rnorm(config$m, c(-5, -5))))
  , t(replicate(500, rnorm(config$m))))
dn = lapply(d, function(j) { data.frame(x = j[,1], y = j[,2], z = j[,3]) })
m  = melt(dn, id.vars = c('x', 'y', 'z'))

set.seed(990909)

params = inverse_model(
    config$n, config$k, m[, c('x', 'y', 'z')]
  , config$a
  , config$l, config$r
  , config$b, config$w
  )

dp = melt(data.frame(params$p))

dm = melt(lapply(params$m, data.frame), id.vars = c('x', 'y', 'z'))

py = ggplot(m, aes(x, y)) + geom_point()

pp = ggplot(dp, aes(seq_along(value), value, colour = variable)) +
       geom_line() + facet_grid(. ~ variable)

pm = ggplot(dm, aes(x, y, colour = factor(L1), fill = factor(L1))) +
       geom_point(alpha = 0.5)

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

