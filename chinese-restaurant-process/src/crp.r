crp = function(n, a) {
  restaurant = data.frame(table = 1, customers = 1)
  for (j in seq(n - 1)) {
    restaurant = arrival(restaurant, a)
    }
  restaurant
  }

arrival = function(r, a) {
  p = 1 - a / (sum(r$customers) + a)
  if (rbinom(1, 1, p)) {
    join_table(r, a)
    } else {
    start_table(r)
    }
  }

join_table = function(r, a) {
  probs = r$customers / sum(r$customers)
  table = sample(1:nrow(r), size = 1, prob = probs)
  r[table, 'customers'] = r[table, 'customers'] + 1
  r
  }

start_table = function(r) {
  rbind(r, c(nrow(r) + 1, 1))
  }

