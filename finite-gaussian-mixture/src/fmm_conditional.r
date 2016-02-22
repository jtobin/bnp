set.seed(42)

require(dplyr)
require(gtools)
require(reshape2)

source('fmm_generative.r')

# FIXME move to simulation module
config = list(
    k = 3
  , a = 1
  , l = 0
  , r = 1
  , b = 1
  , w = 1
  , n = 1000
  )

# FIXME move to simulation module
origin = list(
    p = mixing_model(config$k, config$a)
  , m = location_model(config$k, config$l, config$r)
  , s = precision_model(config$k, config$b, config$w)
  )

# FIXME move to simulation module
data = melt(model(config$k, config$n))

conditional_label_model = function(y, p, m, s) {
  scorer     = function(mix, mu, prec) { mix * dnorm(y, mu, 1 / prec) }
  unweighted = mapply(scorer, p, m, s)
  weights    = 1 / apply(unweighted, MARGIN = 1, sum)
  probs      = weights * unweighted
  apply(probs
    , MARGIN = 1
    , function(row) { sample(seq_along(m), size = 1, prob = row) }
    )
  }

conditional_mixing_model = function(y, k, z, a) {
  labelled = data.frame(value = y, L1 = z)
  counts   = summarise(group_by(labelled, L1), count = length(value))

  concentration = sapply(seq(k),
    function(cluster) {
      idx = which(counts$L1 == cluster)
      if (length(idx) != 0) {
        counts$count[idx] + a / k
        } else {
        0
        }
      })

  rdirichlet(1, concentration)
  }





# mixing_model    = function(k, a) drop(rdirichlet(1, (rep(a, k))))
# label_model     = function(n, p) drop(rmultinom(1, size = n, prob = p))
# location_model  = function(k, l, r) rnorm(k, l, 1 / r)
# precision_model = function(k, b, w) rgamma(k, b, 1 / w)
#
# parameter_model = function(k, n) {
#   p  = mixing_model(k, 1)
#   c  = label_model(n, p)
#   mu = location_model(k, 0, 0.1)
#   s  = precision_model(k, 1, 1)
#   list(c, mu, s)
#   }
#
# data_model = function(config) {
#   sampler = function(y, m, s) rnorm(y, m, 1 / s)
#   mapply(sampler, config[[1]], config[[2]], config[[3]])
#   }
#
# model = function(k, n) parameter_model(k, n) %>% data_model
#
