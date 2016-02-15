
pitman_yor = function(n, a, b) {
  bundle = rbeta(1, 1 - a, b + a)
  for (j in seq(n)) {
    bundle = snap(bundle, a, b, j + 1)
  }
  bundle
  }

snap = function(acc, a, b, k) {
  v = rbeta(1, 1 - a, b + k * a)
  p = v * (1 - sum(acc))
  c(acc, p)
  }


