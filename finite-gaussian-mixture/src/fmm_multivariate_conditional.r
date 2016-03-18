require(mvtnorm)

source('fmm_multivariate_generative.r')
source('fmm_utils.r')

# NOTE (jtobin): must load dplyr after plyr
require(dplyr)

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
  scorer = function(mix, mu, prec) {
    exp(log(mix) + dmvnorm(y, mu, solve(prec), log = T))
  }

  unweighted = mapply(scorer, p, m, s)
  weights    = 1 / apply(unweighted, MARGIN = 1, sum)
  probs      = weights * unweighted

  clusters = apply(
      probs
    , MARGIN = 1
    , function(row) { sample(seq_along(m), size = 1, prob = row) }
    )
  unname(clusters)
}

conditional_cluster_parameters_model = function(y, k, z, l, r, b, w) {
  labelled  = data.frame(y, L1 = z)
  clustered = lapply(seq(k),
    function(j) {
      vals = labelled[which(labelled$L1 == j), !(names(labelled) %in% 'L1')]
      as.matrix(vals)
    })

  ybar     = lapply(clustered, colMeans)
  n        = lapply(clustered, nrow)
  pl       = function(lj, nj, ybarj) {
    if (nj == 0) {
      lj
    } else {
      (lj + nj * ybarj) / (1 + nj)
    }
  }
  ln       = mapply(pl, list(l), n, ybar, SIMPLIFY = F)
  centered = mapply('-', clustered, ybar, SIMPLIFY = F)
  ss       = lapply(centered, function(x) t(x) %*% x)

  # NOTE (jtobin): the extra 'solve' calls here helped; came from
  # http://thaines.com/content/misc/gaussian_conjugate_prior_cheat_sheet.pdf
  # murphy's famous reference at
  # http://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf is incorrect.
  pt = function(wj, ssj, nj, ybarj) {
    if (nj == 0) { wj } else {
      solve(solve(wj) + ssj + nj / (1 + nj) * ((l - ybarj) %*% t(l - ybarj)))
    }
  }

  tn   = mapply(pt, list(w), ss, n, ybar, SIMPLIFY = F)
  bn   = lapply(n, function(x) x + b)
  prec = mapply(function(i, j) drop(rWishart(1, i, j)), bn, tn, SIMPLIFY = F)
  cov  = mapply(function(i, j) solve((i + 1) * j), n, tn, SIMPLIFY = F)
  loc  = mapply(rmvnorm, 1, ln, cov, SIMPLIFY = F)
  if (any(is.nan(unlist(loc)))) { browser() }
  list(m = loc, s = prec)
}

inverse_model = function(n, k, y, a, l, r, b, w) {
  gibbs = function(p0, m0, s0) {
    z  = conditional_label_model(y, p0, m0, s0)
    p1 = conditional_mixing_model(y, k, z, a)
    ps = conditional_cluster_parameters_model(y, k, z, l, r, b, w)
    m1 = ps$m
    s1 = ps$s
    ll = lmodel(y, p1, m1, s1)

    list(p = p1, m = m1, s = s1, z = z, l = ll)
  }

  params = list(
      p = mixing_model(k, a)
    , m = lapply(
              location_model(k, l, r)
            , function(j) { matrix(j, ncol = length(j)) })
    , s = precision_model(k, b, w)
    )

  acc   = params
  acc$s = list(acc$s)
  for (j in seq(n - 1)) {
      params = gibbs(params$p, params$m, params$s)

      acc$p  = rbind(acc$p, params$p)
      acc$m  = mapply(rbind, acc$m, params$m, SIMPLIFY = F)
      acc$s  = c(acc$s, list(params$s))
      acc$z  = rbind(acc$z, params$z)
      acc$l  = c(acc$l, params$l)
    }
  acc
}

