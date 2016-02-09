#  finite gaussian mixture model
#
#  a               ~ inverse-gamma(1, 1)
#  p | a           ~ symmetric-dirichlet(a)
#  c | p           ~ multinomial(p)
#  l               ~ gaussian(mu_y, var_y)
#  r               ~ gamma(1, prec_y)
#  mu | c, l, r    ~ gaussian(l, r^-1)
#  s  | c, b, w    ~ gamma(b, w^-1)
#  b               ~ inverse-gamma(1, 1)
#  w               ~ gamma(1, var_y)
#  y | p, c, mu, s ~ <p, normal(mu, s^-1)>

# conditional posterior densities

dmeans = function(mu, y, s, l, r) {
    cluster_count = unlist(lapply(y, length))
    cluster_mean  = unlist(lapply(y, mean))

    m = (cluster_mean * cluster_count * s + l * r) / (cluster_count * s + r)
    v = 1 / (cluster_count * s + r)

    dnorm(x, m, v)
  }

dl = function(l, y, mu, r) {
    mu_y   = mean(unlist(y))
    var_y  = var(unlist(y))
    prec_y = 1 / var_y
    k      = length(y)

    m = (mu_y * prec_y + r * sum(mu)) / (prec_y + k * r)
    v = 1 / (prec_y + k * r)

    dnorm(l, m, v)
  }

dr = function(r, y, mu, l) {
    var_y = var(unlist(y))
    k     = length(y)

    a = k + 1
    b = a / (var_y + sum((mu - l)^2))

    dgamma(r, a, b)
  }

ds = function(s, y, mu, b, w) {
    cluster_count = unlist(lapply(y, length))

    squares = unlist(mapply("-", y, as.list(mu))) ^ 2

    a   = b + cluster_count
    bet = a / (w * b + sum(squares))

    dgamma(s, a, bet)
  }

dw = function(w, y, s, b) {
    k      = length(y)
    var_y  = var(unlist(y))
    prec_y = 1 / var_y

    a   = k * b + 1
    bet = a / (prec_y + b * sum(s))

    dgamma(w, a, bet)
  }

db = function(b, y, s, w) {
    k = length(y)

    t0 = gamma(b / 2) ^ (- k)
    t1 = exp(-1 / (2 * b))
    t2 = (b / 2) ^ ((k * b - 3) / 2)
    t3 = prod((s * w) ^ (b / 2) * exp(- b * s * w / 2))

    t0 * t1 * t2 * t3
  }

dn = function(n, a) {
    nt = sum(n)
    k  = length(n)
    a ^ k * gamma(a) / gamma(nt + a)
  }

da = function(a, n, k) {
    a ^ (k - 3 / 2) * exp(-1 / (2 * a)) * gamma(a) / gamma(n + a)
  }

# sampling functions


