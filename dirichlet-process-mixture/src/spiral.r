require(MASS)
require(scatterplot3d)

set.seed(42)

n = 800
t = sort(runif(n) * 4 * pi)

x = (13 - 0.5 * t) * cos(t)
y = (13 - 0.5 * t) * sin(t)
Z = mvrnorm(n, mu = rep(0, 3), Sigma = 0.5 * diag(3))

X = matrix(data = c(x, y, t), nrow = n, ncol = 3) + Z

# visualization

# quartz()
# scatterplot3d(X, highlight.3d = T, pch = 19)
