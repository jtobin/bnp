
ibp = function(n, a) {

  dishes = max(1, rpois(1, a))
  diners = data.frame(dish = seq(dishes), diners = rep(1, dishes))
  buffet = list(buffet = diners, choices = list(seq(dishes)))

  for (j in seq(n - 1)) {
      buffet = arrival(buffet, a)
    }
  buffet
  }

arrival = function(b, a) {
  config  = b$buffet
  existing_choices = b$choices
  j = length(existing_choices) + 1
  n = nrow(config)

  num_new_dishes = rpois(1, a / j)
  new_dishes     =
      if (num_new_dishes > 0) {
        (n + 1):(n + num_new_dishes)
      } else {
        numeric(0)
      }

  probs     = sapply(config$diners, function(n) { n / j })
  selection = rbinom(n, 1, prob = probs) == 1
  choices   = list(c(as.integer(row.names(config[selection,])), new_dishes))
  diners    = with(config, replace(diners, selection, diners[selection] + 1))

  buffet = data.frame(
      dish   = seq(n + num_new_dishes)
    , diners = c(diners, rep(1, num_new_dishes))
    )

  list(buffet = buffet, choices = c(existing_choices, choices))
  }

ibp_matrix = function(buffet) {
  dishes = do.call(max, buffet$choices)
  index  = function(j) { as.numeric(seq(dishes) %in% j) }
  rows   = lapply(buffet$choices, index)

  do.call(rbind, rows)
  }

