BNP_DIR = "/Users/jtobin/projects/bnp"
SBP_SRC = paste(BNP_DIR, "stick-breaking-process/src/sbp.r", sep = "/")

require(gtools)
require(magrittr)
source(SBP_SRC)

mixing_model = function(n, a) {
  if (n <= 1) {
    stop("need > 1 observation.")
    } else {
    sbp(n - 1, a)
    }
  }

label_model     = function(n, p) drop(rmultinom(1, size = n, prob = p))
location_model  = function(k, l, r) rnorm(k, l, 1 / r)
precision_model = function(k, b, w) rgamma(k, b, 1 / w)

parameter_model = function(n, a) {
  p  = mixing_model(n, a)
  k  = length(p)
  c  = label_model(n, p)
  mu = location_model(k, 0, 0.1)
  s  = precision_model(k, 1, 1)
  list(c, mu, s)
  }

data_model = function(config) {
  sampler = function(y, m, s) rnorm(y, m, 1 / s)
  mapply(sampler, config[[1]], config[[2]], config[[3]])
  }

model = function(n, a) {
  clusters = parameter_model(n, a) %>% data_model
  nonempty = unlist(lapply(clusters, function(x) length(x) != 0))
  clusters[nonempty]
  }

