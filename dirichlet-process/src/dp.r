BNP_DIR = "/Users/jtobin/projects/bnp"
SBP_SRC = paste(BNP_DIR, "stick-breaking-process/src/sbp.r", sep = "/")

source(SBP_SRC)

# ex: gaussian base measure
#
# > dp(10, 1, function() { rnorm(1) })
dp = function(n, a, h) {
  p = sbp(n - 1, a)
  g = replicate(length(p), h())
  list(p, g)
  }


