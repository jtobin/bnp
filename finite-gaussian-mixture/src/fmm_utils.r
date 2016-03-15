rdirichlet = function(n, a) {
  l = length(a)
  x = matrix(rgamma(n * l, a), ncol = l, byrow = T)
  s = x %*% rep(1, l)
  x / as.vector(s)
}
