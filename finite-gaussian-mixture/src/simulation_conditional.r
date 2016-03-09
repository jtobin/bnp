set.seed(42)

require(ggplot2)
require(reshape2)

source('fmm_conditional.r')

config = list(
    k = 3
  , a = 1
  , l = 0
  , r = 0.1
  , b = 1
  , w = 1
  , n = 1000
  )

origin = list(
    p = mixing_model(config$k, config$a)
  , m = location_model(config$k, config$l, config$r)
  , s = precision_model(config$k, config$b, config$w)
  )

d = melt(model(config$k, config$n))

params = inverse_model(
    config$n, config$k, d$value
  , config$a, config$l, config$r
  , config$b, config$w
  )

dp = melt(as.data.frame(params$p))
dm = melt(as.data.frame(params$m))
ds = melt(as.data.frame(params$s))
dl = melt(as.data.frame(params$l))

py = ggplot(d, aes(value, colour = factor(L1), fill = factor(L1))) +
       geom_density(alpha = 0.5)

pp = ggplot(dp, aes(x = seq_along(value), y = value, colour = variable)) +
       geom_line()

pm = ggplot(dm, aes(x = seq_along(value), y = value, colour = variable)) +
       geom_line()

ps = ggplot(ds, aes(x = seq_along(value), y = value, colour = variable)) +
       geom_line()

pl = ggplot(dl, aes(x = seq_along(value), y = value, colour = variable)) +
       geom_line()

early = data.frame(value = d$value, variable = params$z[1,])
mid   = data.frame(value = d$value, variable = params$z[round(config$n / 2),])
late  = data.frame(value = d$value, variable = params$z[config$n - 1,])

p_early =
  ggplot(early, aes(value, colour = factor(variable), fill = factor(variable))) +
    geom_density(alpha = 0.5)

p_mid =
  ggplot(mid, aes(value, colour = factor(variable), fill = factor(variable))) +
    geom_density(alpha = 0.5)

p_late =
  ggplot(late, aes(value, colour = factor(variable), fill = factor(variable))) +
    geom_density(alpha = 0.5)

