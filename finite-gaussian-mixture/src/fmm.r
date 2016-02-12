set.seed(42)

require(dplyr)
require(gtools)

#  basic finite gaussian mixture model
#
#  p               ~ symmetric-dirichlet(a)
#  c | p           ~ multinomial(p)
#  mu | c          ~ gaussian(l, 1 / r)
#  s  | c          ~ gamma(b, 1 / w)
#  y | p, c, mu, s ~ <p, normal(mu, 1 / s)>

# mixing probabilities

rp_model = function(k, a) drop(rdirichlet(1, (rep(a, k))))

rp_conditional = function(n, a0) {
    k  = length(n)
    a1 = sapply(n, function(x) { x + a0 / k })
    drop(rdirichlet(1, a1))
  }

# cluster labels

rc_model = function(n, p) drop(rmultinom(1, size = n, prob = p))

# FIXME too easy to generate NaN probabilities when precisions are too large
#       may not be correct at all.  perplexing.
rc_conditional = function(y, p0, mu, s) {
  k        = length(p0)
  reducer  = function(p, m, prec) { p * dnorm(y, m, 1 / prec) }
  elements = mapply(reducer, p, m, s)
  p1       = elements / sum(elements)

  sample(1:k, size = 1, prob = p1)
  }

# cluster locations

rmu_model = function(k, l, r) rnorm(k, l, 1 / r)

rmu_conditional = function(y, s, l, r) {
  k    = length(s)
  c    = sapply(y, length)
  ybar = sapply(y, safe_mean)
  m    = (ybar * c * s + l * r) / (c * s + r)
  v    = 1 / (c * s + r)
  rnorm(k, m, v)
  }

# cluster precisions

rs_model = function(k, b, w) rgamma(k, b, 1 / w)

# FIXME seems to be generating enormous precisions
rs_conditional = function(y, mu, b, w) {
  k = length(y)
  c = sapply(y, length)

  centered = mapply("-", y, mu)
  squared  = lapply(centered, function(x) { x ^ 2 })
  ss       = unlist(lapply(squared, sum))

  a   = b + c
  bet = a / (w * b + ss)

  reducer = function(a, b) { rgamma(1, a, b) }
  mapply(reducer, a, bet)
  }

# parameter model (prior)

rconfig_model = function(k, n) {
  p  = rp_model(k, 1)
  c  = rc_model(n, p)
  mu = rmu_model(k, 0, 0.1)
  s  = rs_model(k, 1, 1)
  data.frame(c = c, mu = mu, s = s)
  }

# data model (likelihood)

ry_model = function(config) {
  sampler = function(row) rnorm(row[1], row[2], 1 / row[3])
  apply(config, MARGIN = 1, sampler)
  }

# model

rmodel = function(k, n) {
  config = rconfig_model(k, n)
  ry_model(config)
  }

rmodel_conditional = function(n, y) {
  k = 3
  a = 1
  l = 0
  r = 0.1
  b = 1
  w = 1

  raw = unlist(y)

  kernel = function(p0, mu0, s0) {
    z = sapply(raw, function(x) rc_conditional(x, p0, mu0, s0)) # FIXME slow
    labelled = data.frame(label = z, y = raw)
    id_query = function(c) filter(labelled, label == c)$y

    clustered = sapply(1:k, id_query)
    counts    = sapply(clustered, length)

    p1  = rp_conditional(counts, a)
    mu1 = rmu_conditional(clustered, s0, l, r)
    s1  = rs_conditional(clustered, mu1, b, w)
    list(p = p1, mu = mu1, s = s1)
    }

  gibbs = function(epochs, acc, p0, mu0, s0) {
    if (epochs <= 0) {
        return(acc)
      } else {
          result = kernel(p0, mu0, s0)
          acc$p  = rbind(acc$p, result$p)
          acc$mu = rbind(acc$mu, result$mu)
          acc$s  = rbind(acc$s, result$s)
          gibbs(epochs - 1, acc, result$p, result$mu, result$s)
      }
    }

  p0   = rp_model(k, a)
  mu0  = rmu_model(k, l, r)
  s0   = rs_model(k, b, w)
  init = list(p = p0, mu = mu0, s = s0)

  gibbs(n, init, p0, mu0, s0)
  }

# utilities

safe_mean = function(x) {
  if (is.null(x) || (length(x) == 0)) {
    return(0)
    } else {
    mean(x)
    }
  }

# debug

test_data = list(
    rnorm(101, 10.5, 1)
  , rnorm(38, 0.3, 1)
  , rnorm(90, -8.2, 0.5)
  )

# y = list(
#     rnorm(801, 3.5, 1)
#   , rnorm(300, 0.3, 0.8)
#   , rnorm(722, -4.2, 0.5)
#   )
# k = length(y)
# l = 0
# r = 0.1
# mu = rmu_model(k, l, r)
# w = 1
# b = 1
# s = rs_model(k, b, w)
# p = rp_model(k, 1)
# c = drop(rmultinom(1, 1, prob = p))
# p0 = p
# mu0 = mu
# s0 = s

