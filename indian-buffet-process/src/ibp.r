# FIXME handle sampled zero values

ibp = function(n, a) {
  dishes = rpois(1, a)
  diners = data.frame(dish = seq(dishes), diners = rep(1, dishes))
  buffet = list(buffet = diners, choices = list(seq(dishes)))

  for (j in seq(n - 1)) {
      buffet = arrival(buffet, a)
    }
  buffet
  }

arrival = function(b, a) {

  config  = b$buffet
  choices = b$choices

  j          = length(choices) + 1
  n          = nrow(config)
  probs      = sapply(config$diners, function(n) { n / j })
  selection  = rbinom(n, 1, prob = probs) == 1

  existing_diners = config[selection, 'diners']
  new_diners = with(config, replace(diners, selection, diners[selection] + 1))

  existing_dishes = data.frame(dish = config$dish, diners = new_diners)
  num_new_dishes  = rpois(1, a / j )
  new_dishes = data.frame(
      dish   = (n + 1):(n + num_new_dishes)
    , diners = rep(1, num_new_dishes)
    )

  list(
      buffet  = rbind(existing_dishes, new_dishes)
    , choices = list(choices, rep(TRUE, num_new_dishes)) # FIXME
    )
  }
