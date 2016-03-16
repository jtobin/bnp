require(mvtnorm)

cluster_statistics = function(cluster, l, b, w) {
  n        = nrow(cluster)
  ybar     = colMeans(cluster)
  centered = as.matrix(cluster) - ybar
  ss       = t(centered) %*% centered
  ln       = (l + n * ybar) / (1 + n)
  tn       = solve(solve(w) + ss + n / (1 + n) * ((l - ybar) %*% t(l - ybar)))
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

# FIXME (jtobin): more efficient to cache sufficient statistics in gibbs loop
conditional_label_model = function(y, k, z, a, l, r, b, w) {
  m = ncol(y)
  cluster_labels = seq(k)
  rows           = sample(seq(nrow(y)))

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
    y_censored = y[-i,]
    z_censored = z[-i]
    n_censored = sapply(
        cluster_labels
      , function(j) { length(which(z_censored == j)) })

    score_by_cluster = function(j) {
      sufficient_stats = if (j == old_label) {
          cluster = y_censored[which(z_censored == j), ]
          cluster_statistics(cluster, l, b, w)
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

    # MUTATION
    z[i] <- new_label
    new_stats = cluster_statistics(y[which(z == new_label),], l, b, w)
    sufficient_statistics[[new_label]] <- new_stats

    new_label
  }
  sapply(rows, relabel)
}

inverse_model = function(n, y, k, a, l, r, b, w) {
  # FIXME (jtobin): add likelihood calculation
  gibbs = function(z0) {
    list(z = conditional_label_model(y, k, z0, a, l, r, b, w))
  }
  params = list(z = sample(seq(k), size = nrow(y), replace = T))
  acc    = params
  # FIXME (jtobin): can use replicate
  for (j in seq(n - 1)) {
    params = gibbs(params$z)
    acc$z  = rbind(acc$z, params$z)
  }
  acc
}



# development

require(reshape2) # FIXME move to sim
require(ggplot2)
require(gridExtra)

d = list(
    t(replicate(250, rnorm(2, c(5, 5))))
  , t(replicate(250, rnorm(2, c(-5, -5))))
  , t(replicate(500, rnorm(2))))
dn = lapply(d, function(j) { data.frame(x = j[,1], y = j[,2]) })
m  = melt(dn, id.vars = c('x', 'y'))

dimension = 2

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

foo = inverse_model(100, y, 3, a, l, r, b, w)

early = data.frame(x = y$x, y = y$y, variable = foo$z[1,])
mid   = data.frame(x = y$x, y = y$y, variable = foo$z[round(80),])
late  = data.frame(x = y$x, y = y$y, variable = foo$z[100 - 1,])

p_early =
  ggplot(early, aes(x, y, colour = factor(variable), fill = factor(variable))) +
    geom_point(alpha = 0.5)

p_mid =
  ggplot(mid, aes(x, y, colour = factor(variable), fill = factor(variable))) +
    geom_point(alpha = 0.5)

p_late =
  ggplot(late, aes(x, y, value, colour = factor(variable), fill = factor(variable))) +
    geom_point(alpha = 0.5)

inferred_plots = grid.arrange(p_early, p_mid, p_late, ncol = 3)

