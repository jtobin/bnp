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

# conditional posterior densities / samplers

dr_mu = function(mu, y, s, l, r) {
    cluster_count = unlist(lapply(y, length))
    cluster_mean  = unlist(lapply(y, mean))

    m = (cluster_mean * cluster_count * s + l * r) / (cluster_count * s + r)
    v = 1 / (cluster_count * s + r)

    list(dnorm(x, m, v), rnorm(1, m, v))
  }

dmu = function(mu, y, s, l, r) { dr_mu(mu, y, s, l, r)[[1]] }
rmu = function(mu, y, s, l, r) { dr_mu(mu, y, s, l, r)[[2]] }

dr_l = function(l, y, mu, r) {
    mu_y   = mean(unlist(y))
    var_y  = var(unlist(y))
    prec_y = 1 / var_y
    k      = length(y)

    m = (mu_y * prec_y + r * sum(mu)) / (prec_y + k * r)
    v = 1 / (prec_y + k * r)

    list(dnorm(l, m, v), rnorm(1, m, v))
  }

dl = function(l, y, mu, r) { dr_l(l, y, mu, r)[[1]] }
rl = function(l, y, mu, r) { dr_l(l, y, mu, r)[[2]] }

dr_r = function(r, y, mu, l) {
    var_y = var(unlist(y))
    k     = length(y)

    a = k + 1
    b = a / (var_y + sum((mu - l)^2))

    list(dgamma(r, a, b), rgamma(1, a, b))
  }

dr = function(r, y, mu, l) { dr_r(r, y, mu, l)[[1]] }
rr = function(r, y, mu, l) { dr_r(r, y, mu, l)[[2]] }

dr_s = function(s, y, mu, b, w) {
    cluster_count = unlist(lapply(y, length))

    squares = unlist(mapply("-", y, as.list(mu))) ^ 2

    a   = b + cluster_count
    bet = a / (w * b + sum(squares))

    list(dgamma(s, a, bet), rgamma(1, a, bet))
  }

ds = function(s, y, mu, b, w) { ds_r(s, y, mu, b, w)[[1]] }
rs = function(s, y, mu, b, w) { ds_r(s, y, mu, b, w)[[2]] }

dr_w = function(w, y, s, b) {
    k      = length(y)
    var_y  = var(unlist(y))
    prec_y = 1 / var_y

    a   = k * b + 1
    bet = a / (prec_y + b * sum(s))

    lis(dgamma(w, a, bet), rgamma(1, a, bet))
  }

dw = function(w, y, s, b) { dr_w(w, y, s, b)[[1]] }
rw = function(w, y, s, b) { dr_w(w, y, s, b)[[2]] }

dr_b = function(b, y, s, w) {
    k = length(y)

    t0 = gamma(b / 2) ^ (- k)
    t1 = exp(-1 / (2 * b))
    t2 = (b / 2) ^ ((k * b - 3) / 2)
    t3 = prod((s * w) ^ (b / 2) * exp(- b * s * w / 2))

    list(t0 * t1 * t2 * t3, NULL) # FIXME sampling function
  }

db = function(b, y, s, w) { dr_b(b, y, s, w)[[1]] }

dr_n = function(n, a) {
    nt = sum(n)
    k  = length(n)
    list(a ^ k * gamma(a) / gamma(nt + a), NULL) # FIXME sampling function
  }

dn = function(n, a) { dr_n(n, a)[[1]] }

dr_a = function(a, n, k) {
    list(a ^ (k - 3 / 2) * exp(-1 / (2 * a)) * gamma(a) / gamma(n + a), NULL)
    # FIXME sampling function
  }

da = function(a, n, k) { dr_a(a, n, k)[[1]] }

