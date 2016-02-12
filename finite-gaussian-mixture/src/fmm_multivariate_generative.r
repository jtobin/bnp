set.seed(42)

require(gtools)
require(mvtnorm)

mixing_model    = function(k, a) drop(rdirichlet(1, (rep(a, k))))
label_model     = function(n, p) drop(rmultinom(1, size = n, prob = p))
location_model  = function(k, l, r) rmvnorm(k, l, solve(r))
precision_model = function(k, b, w) rinvwishart(k, b, solve(w))

parameter_model = function(k, n) {
  p  = mixing_model(k, 1)
  c  = lapply(label_model(n, p), list)
  mu = apply(location_model(k, rep(0, k), diag(10, k)), MARGIN = 1, list)
  s  = precision_model(k, 10, diag(1, k))
  list(c, mu, s)
  }

# FIXME mapply not working here
data_model = function(config) {
  sampler = function(c, m, s) rmvnorm(c, m, solve(s))
  mapply(sampler, config[[1]], config[[2]], config[[3]])
  }

model = function(k, n) {
  config = parameter_model(k, n)
  data_model(config)
  }

rinvwishart = function(n, v, S) {
  wishes = rWishart(n, v, solve(S))
  invs   = apply(wishes, MARGIN = 3, function(x) list(solve(x)))
  delabel(invs)
  }

