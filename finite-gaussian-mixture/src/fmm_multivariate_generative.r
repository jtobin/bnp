require(gtools)
require(magrittr)
require(mvtnorm)

mixing_model = function(k, a) drop(rdirichlet(1, (rep(a, k))))

label_model = function(n, p) {
  vals = drop(rmultinom(1, size = n, prob = p))
  delabel(lapply(vals, list))
}

location_model  = function(k, l, r) {
  vals = rmvnorm(k, l, solve(r))
  delabel(apply(vals, MARGIN = 1, list))
}

precision_model = function(k, b, w) rinvwishart(k, b, solve(w))

parameter_model = function(m, k, b, n) {
  p  = mixing_model(k, 1)
  c  = label_model(n, p)
  mu = location_model(k, rep(0, m), diag(0.05, m))
  s  = precision_model(k, b, diag(1, m))
  list(n = c, m = mu, s = s)
}

data_model = function(config) {
  mapply(safe_rmvnorm, config$n, config$m, config$s)
}

model = function(m, k, b, n) parameter_model(m, k, b, n) %>% data_model

# FIXME (jtobin): checkme, not correct
lmodel = function(y, z, p, m, s) {

  clustered = cbind(y, L1 = z)
  cluster   = clustered$L1

  score     = log(p[cluster]) +
    dmvnorm(clustered$value, m[cluster], solve(s[cluster]), log = T)

  sum(score)
}

# utilities

rinvwishart = function(n, v, S) {
  wishes = rWishart(n, v, solve(S))
  delabel(apply(wishes, MARGIN = 3, function(x) list(solve(x))))
  }

delabel = function(x) lapply(x, "[[", 1)

safe_rmvnorm = function(c, m, s) {
  if (c <= 0) return(numeric(0))
  else rmvnorm(c, m, solve(s))
  }
