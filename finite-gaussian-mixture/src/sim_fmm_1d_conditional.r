require(ggplot2)
require(gridExtra)
require(reshape2)

source('fmm_conditional.r')

config = list(
    k = 3
  , a = 1
  , l = 0
  , r = 0.01
  , b = 1
  , w = 1
  , n = 1000
  )

set.seed(990909)

d = data.frame(
  value = c(rnorm(250, -3, 0.25), rnorm(500, 0, 0.25), rnorm(250, 3, 0.25)))

set.seed(990909)

params = inverse_model(
    config$n, config$k, d$value
  , config$a, config$l, config$r
  , config$b, config$w
  )

dp = melt(as.data.frame(params$p))
dm = melt(as.data.frame(params$m))
ds = melt(as.data.frame(params$s))
dl = melt(as.data.frame(params$l))

py = ggplot(d, aes(value)) + geom_histogram(alpha = 0.5, fill = 'darkblue')

pp = ggplot(dp, aes(x = seq_along(value), y = value, colour = variable)) +
       geom_line()

pm = ggplot(dm, aes(x = seq_along(value), y = value, colour = variable)) +
       geom_line()

ps = ggplot(ds, aes(x = seq_along(value), y = log(value), colour = variable)) +
       geom_line()

pl = ggplot(dl, aes(x = seq_along(value), y = value)) +
       geom_line(colour = 'darkblue')

early = data.frame(value = d$value, variable = params$z[1,])
mid   = data.frame(value = d$value, variable = params$z[round(config$n / 2),])
late  = data.frame(value = d$value, variable = params$z[config$n - 1,])

p_early =
  ggplot(early, aes(value, colour = factor(variable), fill = factor(variable))) +
    geom_histogram(alpha = 0.5)

p_mid =
  ggplot(mid, aes(value, colour = factor(variable), fill = factor(variable))) +
    geom_histogram(alpha = 0.5)

p_late =
  ggplot(late, aes(value, colour = factor(variable), fill = factor(variable))) +
    geom_histogram(alpha = 0.5)

chain_plots    = grid.arrange(py, pp, pm, ps, nrow = 2, ncol = 2)

inferred_plots = grid.arrange(py, p_early, p_mid, p_late, nrow = 2, ncol = 2)

