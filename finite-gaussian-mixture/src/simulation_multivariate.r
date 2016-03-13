require(ggplot2)
require(reshape2)
require(scatterplot3d)

source('fmm_multivariate_generative.r')

# 2d

config = list(
    m = 2
  , v = 2
  , k = 4
  , n = 10000
  )

set.seed(42)

d      = model(config$m, config$k, config$v, config$n)
framed = lapply(d, function(mat) { data.frame(x = mat[,1], y = mat[,2]) })
melted = melt(framed, id.vars = c('x', 'y'))
p = ggplot(melted, aes(x, y, colour = factor(L1))) + geom_point(alpha = 0.2)

# 3d

set.seed(42)

config_3d = list(
    m = 3
  , v = 3
  , k = 4
  , n = 10000
  )

d_3d      = model(config_3d$m, config_3d$k, config_3d$v, config_3d$n)
framed_3d = lapply(d_3d,
  function(mat) { data.frame(x = mat[,1], y = mat[,2], z = mat[,3]) })
melted_3d = do.call(rbind, framed_3d)

scatterplot3d(melted_3d, highlight.3d = T, pch = 19)

