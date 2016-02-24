BNP_DIR = "/Users/jtobin/projects/bnp"
CRP_SRC = paste(BNP_DIR, "chinese-restaurant-process/src/crp.r", sep = "/")

require(dplyr)
source(CRP_SRC)

frequencies = function(restaurant) {
  summarised =
    restaurant %>%
    group_by(table) %>%
    summarise(n = sum(customers)) %>%
    mutate(freq = n / sum(n))
  summarised$freq
  }

# FIXME (jtobin): this is not correct
crf = function(n, a, a1) {
  g0 = frequencies(crp(n, a))
  h0 = function() { sample(g0, size = 1) }
  g1 = frequencies(crp(n, a1))
  g2 = frequencies(crp(n, a1))
  list(g1 * g0, g2 * g0)
  }

