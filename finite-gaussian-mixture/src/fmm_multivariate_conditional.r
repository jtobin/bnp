require(dplyr)
require(gtools)
require(mvtnorm)

source('fmm_multivariate_generative.r')

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

conditional_location_model = function(y, z, s, l, r) {
  labelled = cbind(y, L1 = z)
  cluster  = function(d, j) {
    vals = d[which(d$L1 == j), !(names(d) %in% 'L1')]
  }

  clustered = lapply(seq_along(s), function(j) { cluster(labelled, j) })
  n    = lapply(clustered, nrow)
  yt   = lapply(clustered, function(j) { apply(j, MARGIN = 2, sum) })
  num0 = mapply('%*%', yt, s, SIMPLIFY = F)
  num  = lapply(num0, function(z) { z + (l %*% r) })
  den0 = mapply('*', n, s, SIMPLIFY = F)
  den  = lapply(den0, function(z) z + r)

  v = lapply(den, solve)
  m = mapply('%*%', num, v, SIMPLIFY = F)
  mapply(rmvnorm, 1, m, v, SIMPLIFY = F)
}

conditional_precision_model = function(y, z, m, b, w) {
  labelled = cbind(y, L1 = z)
  cluster  = function(d, j) {
    vals = d[which(d$L1 == j), !(names(d) %in% 'L1')]
  }

  clustered = lapply(seq_along(m), function(j) { cluster(labelled, j) })
  yt = lapply(clustered, function(foo) { apply(foo, MARGIN = 2, sum) })

  center = function(i, j) {
    if (nrow(i) == 0) { as.matrix(i) } else { as.matrix(i - j) }
  }

  centered = mapply(center, clustered, m, SIMPLIFY = F)
  ss       = lapply(centered, function(x) t(x) %*% x)
  n        = lapply(clustered, nrow)
  a        = lapply(n, function(j) j + b)
  bet0     = lapply(ss, function(j) { (j + w * b) })
  bet      = mapply('/', bet0, a, SIMPLIFY = F)

  mapply(function(i, j) drop(rWishart(1, i, j)), a, bet, SIMPLIFY = F)
}

inverse_model = function(n, k, y, a, l, r, b, w) {
  gibbs = function(p0, m0, s0) {
    z  = conditional_label_model(y, p0, m0, s0)
    p1 = conditional_mixing_model(y, k, z, a)
    m1 = conditional_location_model(y, z, s0, l, r)
    s1 = conditional_precision_model(y, z, m1, b, w)
    # FIXME (jtobin): log scores
    # l  = lmodel(y, z, p1, m1, s1)
    list(p = p1, m = m1, s = s1, z = z) # l = l)
    }

  p0     = mixing_model(k, a)
  m0     = location_model(k, l, r)
  s0     = precision_model(k, b, w)
  params = list(
      p = p0
    , m = lapply(m0, function(j) { matrix(j, ncol = length(j)) })
    , s = s0
    )

  acc = params
  for (j in seq(n - 1)) {
      params = gibbs(params$p, params$m, params$s)

      acc$p  = rbind(acc$p, params$p)
      acc$m  = mapply(rbind, acc$m, params$m, SIMPLIFY = F)
      # FIXME (jtobin): not logging intermediate covariances
      #                 might be desirable to log some reduced ellipse dims
      acc$s  = params$s
      acc$z  = rbind(acc$z, params$z)
      # FIXME (jtobin): log scores
      # acc$l  = c(acc$l, params$l)
    }
  acc
  }

