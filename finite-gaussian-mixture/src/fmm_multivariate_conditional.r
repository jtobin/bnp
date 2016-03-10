require(dplyr)
require(gtools)
require(mvtnorm)
require(reshape2) # FIXME move to sim module

source('fmm_multivariate_generative.r')

# FIXME (jtobin): move to simulation module
set.seed(42)

# FIXME (jtobin): move to simulation module
dimension = 2

# FIXME (jtobin): move to simulation module
config = list(
    k = 3
  , m = dimension
  , a = 1
  , l = rep(0, dimension)
  , r = diag(0.05, dimension)
  , b = 2
  , w = diag(1, dimension)
  , n = 1000
  )

# FIXME (jtobin): move to simulation module
origin = list(
    p = mixing_model(config$k, config$a)
  , m = location_model(config$k, config$l, config$r)
  , s = precision_model(config$k, config$b, config$w)
  )

# FIXME (jtobin): move to simulation module
d = melt(model(config$m, config$k, config$n), id.vars = c('x', 'y'))

# y is a nxm matrix, z is a nx1 vector
conditional_mixing_model = function(y, k, z, a) {
  labelled      = cbind(y, L1 = z)
  counts        = summarise(group_by(labelled, L1), count = n())
  concentration = sapply(seq(k),
    function(cluster) {
      idx = which(counts$L1 == cluster)
      if (length(idx) != 0) {
        counts$count[idx] + a / k
      } else {
        a / k
      }
    })
  drop(rdirichlet(1, concentration))
}

# FIXME (jtobin): seems ok but need to check
conditional_label_model = function(y, p, m, s) {
  scorer     = function(mix, mu, prec) {
    exp(log(mix) + dmvnorm(y, mu, solve(prec), log = T))
  }

  unweighted = mapply(scorer, p, m, s)
  weights    = 1 / apply(unweighted, MARGIN = 1, sum)
  weighted   = weights * unweighted

  probabilize = function(row) {
    rs = sum(row)
    if (rs == 0 || is.na(rs) || is.nan(rs)) {
      drop(rdirichlet(1, rep(1, length(m))))
    } else {
      row
    }
  }

  probs    = t(apply(weighted, MARGIN = 1, probabilize))
  clusters = apply(
      probs
    , MARGIN = 1
    , function(row) { sample(seq_along(m), size = 1, prob = row) }
    )
  unname(clusters)
  }

# FIXME (jtobin): this will change quite a bit
conditional_location_model = function(y, z, s, l, r) {
  labelled = cbind(y, L1 = z)

  cluster = function(d, j) {
    vals = d[which(d$L1 == j), !(names(d) %in% 'L1')]
  }

  clustered = lapply(seq_along(s), function(j) { cluster(labelled, j) })
  lengths   = lapply(clustered, nrow)
  sums      = lapply(clustered, function(foo) { apply(foo, MARGIN = 2, sum) })

  n  = lengths
  yt = sums

  # FIXME (jtobin): these must be multivariate quantities
  m  = (yt * s + l * r) / (n * s + r)
  v  = 1 / (n * s + r)

  # FIXME (jtobin): needs to be rmvnorm
  mapply(rnorm, 1, m, v)
  }

# FIXME (jtobin): this will change quite a bit
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

  n = sapply(seq_along(m),
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

inverse_model = function(n, k, y, a, l, r, b, w) {
  gibbs = function(p0, m0, s0) {
    z  = conditional_label_model(y, p0, m0, s0)
    p1 = conditional_mixing_model(y, k, z, a)
    m1 = conditional_location_model(y, z, s0, l, r)
    s1 = conditional_precision_model(y, z, m1, b, w)
    l  = lmodel(y, z, p1, m1, s1)
    list(p = p1, m = m1, s = s1, z = z, l = l)
    }

  p0     = mixing_model(k, a)
  m0     = location_model(k, l, r)
  s0     = precision_model(k, b, w)
  params = list(p = p0, m = m0, s = s0)

  acc = params
  for (j in seq(n - 1)) {
      params = gibbs(params$p, params$m, params$s)
      acc$p  = rbind(acc$p, params$p)
      acc$m  = rbind(acc$m, params$m)
      acc$s  = rbind(acc$s, params$s)
      acc$z  = rbind(acc$z, params$z)
      acc$l  = c(acc$l, params$l)
    }
  acc
  }

