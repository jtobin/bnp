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

parameter_model = function(m, k, n) {
  p  = mixing_model(k, 1)
  c  = label_model(n, p)
  mu = location_model(k, rep(0, m), diag(0.05, m))
  s  = precision_model(k, 2, diag(1, m))
  list(c, mu, s)
  }

data_model = function(config) {
  raw    = mapply(safe_rmvnorm, config[[1]], config[[2]], config[[3]])
  frame  = function(m) data.frame(x = m[,1], y = m[,2])
  lapply(raw, frame)
  }

model = function(m, k, n) parameter_model(m, k, n) %>% data_model

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
