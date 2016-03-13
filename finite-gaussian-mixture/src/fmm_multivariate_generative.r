require(gtools)
require(magrittr)
require(mvtnorm)
require(plyr)

# Cluster probabilities via a symmetric Dirichlet distribution.
#
# k : integer, > 0
# a : numeric, > 0
#
# Returns a vector of probabilities of size k.
mixing_model = function(k, a) drop(rdirichlet(1, (rep(a, k))))

# Number of observations per cluster, by indicator.
#
# n : integer, > 0
# p : numeric, probability
#
# Returns a list of integer sizes corresponding to the given cluster.  The
# number of clusters is determined by the length of 'p'.
label_model = function(n, p) {
  vals = drop(rmultinom(1, size = n, prob = p))
  as.list(vals)
}

# Location, by cluster.
#
# k : integer, > 0
# l : numeric
# r : numeric, positive definite
#
# Returns a list of 'k' locations, each having same dimension as 'l'.
location_model  = function(k, l, r) {
  vals = rmvnorm(k, l, solve(r))
  alply(vals, 1)
}

# Precision matrix, by cluster.
#
# k : integer, > 0
# b : numeric, >= 1
# w : numeric, positive definite
#
# Returns a list of 'k' precision matrices with same dimension as 'w'.
precision_model = function(k, b, w) {
  vals = rWishart(k, b, w)
  alply(vals, 3)
}

# Parameter model for the finite Gaussian mixture model.
#
# k : integer, > 0
# l : numeric
# r : numeric, positive definite
# b : numeric, >= 1
# w : numeric, positive definite
# n : integer, > 0
#
# Returns a list of three components:
#   * n : list of length 'k' containing the size of the kth cluster
#   * m : list of length 'k' containing the location of the kth cluster
#   * s : list of length 'k' containing the precision of the kth cluster
parameter_model = function(k, l, r, b, w, n) {
  p  = mixing_model(k, 1)
  c  = label_model(n, p)
  mu = location_model(k, l, r)
  s  = precision_model(k, b, w)
  list(n = c, m = mu, s = s)
}

# Data model for the finite Gaussian mixture model.
#
# params : output type of 'paramter_model'
#
# Returns observations by cluster as a list.
data_model = function(params) {
  safe_rmvnorm = function(c, m, s) {
    if (c <= 0) { numeric(0) } else { rmvnorm(c, m, solve(s)) }
  }
  mapply(safe_rmvnorm, params$n, params$m, params$s)
}

# The finite Gaussian mixture model.
model = function(k, l, r, b, w, n) {
  parameter_model(k, l, r, b, w, n) %>% data_model
}

# Log-likelihood for the finite Gaussian mixture model.
#
# y : numeric
# p : numeric, probability
# m : numeric
# s : numeric, positive definite
#
# 'y' is a matrix of observations. 'p', 'm', and 's' are a probability vector,
# list of location vectors, and list of precision matrices of the appropriate
# dimensions.
lmodel = function(y, p, m, s) {
  score      = function(pr, mu, prec) { pr * dmvnorm(y, mu, solve(prec)) }
  by_cluster = mapply(score, p, m, s)
  totalled   = apply(by_cluster, MARGIN = 1, sum)

  # NOTE (jtobin): adjusted for numerical stability
  small    = 1.379783e-316
  adjusted = totalled
  adjusted[which(adjusted == 0)] = small
  sum(log(adjusted))
}

