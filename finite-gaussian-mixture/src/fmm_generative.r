require(gtools)
require(magrittr)

mixing_model    = function(k, a) drop(rdirichlet(1, (rep(a, k))))
label_model     = function(n, p) drop(rmultinom(1, size = n, prob = p))
location_model  = function(k, l, r) rnorm(k, l, sqrt(1 / r))
precision_model = function(k, b, w) rgamma(k, b, 1 / w)

parameter_model = function(k, n) {
  p  = mixing_model(k, 1)
  c  = label_model(n, p)
  mu = location_model(k, 0, 0.1)
  s  = precision_model(k, 1, 1)
  list(c, mu, s)
}

data_model = function(config) {
  sampler = function(y, m, s) rnorm(y, m, 1 / s)
  mapply(sampler, config[[1]], config[[2]], config[[3]])
}

model = function(k, n) parameter_model(k, n) %>% data_model

lmodel = function(y, p, m, s) {
  score      = function(pr, mu, prec) { pr * dnorm(y, mu, sqrt(1 / prec)) }
  by_cluster = mapply(score, p, m, s)
  totalled   = apply(by_cluster, MARGIN = 1, sum)

  # NOTE (jtobin): adjusted for numerical stability
  small    = 1.379783e-316
  adjusted = totalled
  adjusted[which(adjusted == 0)] = small
  sum(log(adjusted))
}

