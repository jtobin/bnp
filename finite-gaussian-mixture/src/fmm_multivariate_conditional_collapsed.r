require(mvtnorm)

source('fmm_multivariate_generative.r')

cluster_statistics = function(cluster, l, b, w) {
  mclust   =
    # R, seriously?
    if (is.null(dim(cluster))) {
      matrix(cluster, ncol = ncol(w))
    } else {
      as.matrix(cluster)
    }
  m        = ncol(mclust)
  n        = nrow(mclust)
  ybar     = colMeans(mclust)
  centered = mclust - ybar
  ss       = t(centered) %*% centered
  ln       = (l + n * ybar) / (1 + n)
  tn       =
    if (n == 0) {
      w
    } else {
      solve(solve(w) + ss + n / (1 + n) * ((l - ybar) %*% t(l - ybar)))
    }
  bn       = b + n
  df       = bn - m + 1
  mu       = ln
  coef     = n / ((n + 1) * df)
  v        = coef * solve(tn)
  list(
      n = n, ln = ln, tn = tn, bn = bn
    , df = df, mu = mu, v = v
    )
}

conditional_label_model = function(y, k, z, a, l, b, w) {
  cluster_labels = seq(k)
  rows           = seq(nrow(y))
  m              = ncol(y)

  initial_clusters = sapply(
      cluster_labels
    , function(j) { y[which(z == j),] }
    , simplify = F)

  sufficient_statistics = lapply(
      initial_clusters
    , function(c) { cluster_statistics(c, l, b, w) })

  relabel = function(i) {
    old_label  = z[i]
    val        = y[i,]
    y_censored = as.matrix(y[-i,])
    z_censored = z[-i]
    n_censored = sapply(
        cluster_labels
      , function(j) { length(which(z_censored == j)) })

    score_by_cluster = function(j) {
      sufficient_stats = if (j == old_label) {
          cluster = y_censored[which(z_censored == j), ]
          sufficient_statistics[[j]] <<- cluster_statistics(cluster, l, b, w)
          sufficient_statistics[[j]]
        } else {
          sufficient_statistics[[j]]
        }
      dmvt(
          val
        , df    = sufficient_stats$df
        , sigma = sufficient_stats$v
        , delta = sufficient_stats$mu
        , log   = T
        )
    }

    scores    = exp(sapply(cluster_labels, score_by_cluster))
    weight    = n_censored + a / k
    probs     = scores * weight / sum(scores * weight)
    new_label = sample(cluster_labels, size = 1, prob = probs)

    z[i] <<- new_label
    new_stats = cluster_statistics(y[which(z == new_label),], l, b, w)
    sufficient_statistics[[new_label]] <<- new_stats

    new_label
  }
  sapply(rows, relabel)
}

inverse_model = function(n, k, y, a, l, b, w) {
  gibbs = function(z0) {
    z = conditional_label_model(y, k, z0, a, l, b, w)
    clustered = lapply(seq(k),
      function(j) {
        vals = y[which(z == j),]
        as.matrix(vals)
      })

    ps    = lapply(clustered, function(j) { nrow(j) / nrow(y) })
    mus   = lapply(clustered, colMeans)
    precs = lapply(clustered, function(j) (solve(cov(j))))
    ll    = lmodel(y, ps, mus, precs)
    list(z = z, ll = ll)
  }

  params = list(z = sample(seq(k), size = nrow(y), replace = T))
  acc    = params
  for (j in seq(n - 1)) {
    params = gibbs(params$z)
    acc$z  = rbind(acc$z, params$z)
    acc$ll = c(acc$ll, params$ll)
  }
  acc
}

