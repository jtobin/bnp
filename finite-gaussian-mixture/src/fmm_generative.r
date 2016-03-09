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

lmodel = function(y, z, p, m, t) {
  clustered = data.frame(value = y, L1 = z)
  cluster   = clustered$L1
  score     = log(p[cluster]) +
    dnorm(clustered$value, m[cluster], sqrt(1 / p[cluster]), log = T)
  sum(score)
}

