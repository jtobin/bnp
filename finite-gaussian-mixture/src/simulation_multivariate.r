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

alt_config = list(
    m = 3
  , v = 3
  , k = 4
  , n = 10000
  )

alt_d      = model(alt_config$m, alt_config$k, alt_config$v, alt_config$n)
alt_framed = lapply(alt_d,
  function(mat) { data.frame(x = mat[,1], y = mat[,2], z = mat[,3]) })
alt_melted = do.call(rbind, alt_framed)
scatterplot3d(alt_melted, highlight.3d = T, pch = 19)
