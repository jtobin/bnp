require(reshape2)
require(scatterplot3d)

source('../src/fmm_multivariate_generative.r')

dimension = 3

config_3d = list(
    k = 4
  , l = rep(0, dimension)
  , r = diag(0.05, dimension)
  , b = dimension
  , w = diag(1, dimension)
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

