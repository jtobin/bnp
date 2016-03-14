require(ggplot2)
require(reshape2)
require(scatterplot3d)

source('fmm_multivariate_generative.r')

# 2d example

config = list(
    k = 4
  , l = rep(0, 2)
  , r = diag(0.05, 2)
  , b = 2
  , w = diag(1, 2)
  , n = 10000
  )

set.seed(42)

d = model(
    config$k, config$l, config$r
  , config$b, config$w, config$n
  )

framed = lapply(d, function(mat) { data.frame(x = mat[,1], y = mat[,2]) })
melted = melt(framed, id.vars = c('x', 'y'))
p = ggplot(melted, aes(x, y, colour = factor(L1))) + geom_point(alpha = 0.2)

# 3d example

config_3d = list(
    k = 4
  , l = rep(0, 3)
  , r = diag(0.05, 3)
  , b = 3
  , w = diag(1, 3)
  , n = 10000
  )

set.seed(42)

d_3d = model(
    config_3d$k, config_3d$l, config_3d$r
  , config_3d$b, config_3d$w, config_3d$n
  )

framed_3d = lapply(d_3d,
  function(mat) { data.frame(x = mat[,1], y = mat[,2], z = mat[,3]) })
melted_3d = do.call(rbind, framed_3d)
scatterplot3d(melted_3d, highlight.3d = T, pch = 19)

