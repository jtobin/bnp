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

conditional_location_model = function(y, z, s, l, r) {
  clustered = group_by(data.frame(value = y, L1 = z), L1)
  lengths   = summarise(clustered, value = n())
  sums      = summarise(clustered, value = sum(value))

  n = sapply(seq_along(s),
    function(cluster) {
      idx = which(lengths$L1 == cluster)
      if (length(idx) != 0) {
        lengths$value[idx]
        } else {
        0
        }
      })

  yt = sapply(seq_along(s),
    function(cluster) {
      idx = which(sums$L1 == cluster)
      if (length(idx) != 0) {
        sums$value[idx]
        } else {
        0
        }
      })

  m  = (yt * s + l * r) / (n * s + r)
  v  = 1 / (n * s + r)

  mapply(rnorm, 1, m, v)
  }

conditional_precision_model = function(y, z, m, b, w) {
  labelled  = data.frame(value = y, L1 = z)
  clustered = group_by(labelled, L1)

  acc = list()
  for (j in seq_along(m)) {
    acc[[j]] = labelled[which(labelled$L1 == j), 'value']
    }

  centered = mapply("-", acc, m)
  squared  = lapply(centered, function(x) x ^ 2)
  ss       = unlist(lapply(squared, sum))

  n = sapply(seq_along(s),
    function(cluster) {
      lengths = summarise(clustered, value = n())
      idx     = which(lengths$L1 == cluster)
      if (length(idx) != 0) {
        lengths$value[idx]
        } else {
        0
        }
      })

  a   = b + n
  bet = (w * b + ss) / a

  mapply(function(a, b) rgamma(1, a, b), a, bet)
  }


