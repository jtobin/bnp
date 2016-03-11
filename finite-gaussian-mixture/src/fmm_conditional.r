require(dplyr)
require(gtools)

source('fmm_generative.r')

conditional_mixing_model = function(y, k, z, a) {
  labelled = data.frame(value = y, L1 = z)
  counts   = summarise(group_by(labelled, L1), count = n())

  concentration = sapply(
      seq(k)
    , function(cluster) {
        idx = which(counts$L1 == cluster)
        if (length(idx) != 0) {
          counts$count[idx] + a / k
        } else {
          a / k
        }
      })

  drop(rdirichlet(1, concentration))
  }

conditional_label_model = function(y, p, m, s) {
  scorer     = function(mix, mu, prec) {
    exp(log(mix) + dnorm(y, mu, sqrt(1 / prec), log = T))
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

  probs = t(apply(weighted, MARGIN = 1, probabilize))
  apply(
      probs
    , MARGIN = 1
    , function(row) { sample(seq_along(m), size = 1, prob = row) }
    )
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

  mapply(rnorm, 1, m, sqrt(v))
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
    l  = lmodel(y, p1, m1, s1)
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

