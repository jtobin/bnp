require(ggplot2)
require(reshape2)

source('../src/fmm_multivariate_generative.r')

dimension = 2

config = list(
    k = 4
  , l = rep(0, dimension)
  , r = diag(0.05, dimension)
  , b = dimension
  , w = diag(1, dimension)
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

