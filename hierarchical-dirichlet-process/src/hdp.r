BNP_DIR = "/Users/jtobin/projects/bnp"
DP_SRC  = paste(BNP_DIR, "dirichlet-process/src/dp.r", sep = "/")

source(DP_SRC)

hdp = function(n, a, h, n1, n2, a1, a2) {
  g0 = dp(n, a, h)
  h0 = function() { sample(g0[[1]], size = 1) }
  g1 = dp(n1, a1, h0)
  g2 = dp(n1, a1, h0)
  h1 = function() { sample(g1[[1]], size = 1) }
  h2 = function() { sample(g2[[1]], size = 1) }
  list(dp(n2, a2, h1), dp(n2, a2, h2))
  }

