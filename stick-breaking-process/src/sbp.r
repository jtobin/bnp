sbp = function(n, a) {
  bundle = list(0, numeric(0))
  for (j in seq(n)) {
    bundle = snap(bundle[[1]], bundle[[2]], a)
    }
  bundle[[2]]
  }

snap = function(acc, bun, a) {
  b      = rbeta(1, 1, a)
  stick  = exp(log(b) + acc)

  nacc = log (1 - b) + acc
  nbun = c(bun, stick)
  list(nacc, nbun)
  }

