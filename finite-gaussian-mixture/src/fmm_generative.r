set.seed(42)

require(gtools)

mixing_model    = function(k, a) drop(rdirichlet(1, (rep(a, k))))
label_model     = function(n, p) drop(rmultinom(1, size = n, prob = p))
location_model  = function(k, l, r) rnorm(k, l, 1 / r)
precision_model = function(k, b, w) rgamma(k, b, 1 / w)

parameter_model = function(k, n) {
  p  = mixing_model(k, 1)
  c  = label_model(n, p)
  mu = location_model(k, 0, 0.1)
  s  = precision_model(k, 1, 1)
  list(c, mu, s)
  }

data_model = function(config) {
  sampler = function(y, m, s) rnorm(y, m, 1 / s) # FIXME this may not do what i expect
  mapply(sampler, config[[1]], config[[2]], config[[3]])
  }

model = function(k, n) {
  config = parameter_model(k, n)
  data_model(config)
  }

