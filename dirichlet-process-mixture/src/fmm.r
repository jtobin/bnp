require(gtools)

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

test_data = list(
    c(1.2, 1.1, 1.23, 0.93, 0.89, 1.01)
  , c(0.3, 0.34, 0.32, 0.29, 0.31)
  , c(-0.5, -0.1, -0.21, 0.02, -0.81, 0.08, -0.33)
  )

# observations

y_prior = function(mu, s, p) {
  density = function(y) { sum(p * dnorm(y, mu, 1 / s)) }
  sampler = function() {
      c = t(rmultinom(1, 1, prob = p))
      v = as.vector(c * rnorm(1, mu, 1 / s))
      v[which(v != 0)]
    }
  list(density, sampler)
  }

dy_prior = function(y, mu, s, p) { y_prior(mu, s, p)[[1]](y) }
ry_prior = function(mu, s, p) { y_prior(mu, s, p)[[2]] }

# mu (component means) and parameter models

l_obj_prior = function(y) {
  mu_y   = mean(unlist(y))
  var_y  = var(unlist(y))

  density = function(l) { dnorm(l, mu_y, var_y) }
  sampler = rnorm(1, mu_y, var_y)
  list(density, sampler)
  }

dl_obj_prior = function(l, y) { l_obj_prior(y)[[1]](l) }
rl_obj_prior = function(y) { l_obj_prior(y)[[2]] }

l_posterior = function(y, mu, r) {
  mu_y   = mean(unlist(y))
  var_y  = var(unlist(y))
  prec_y = 1 / var_y
  k      = length(y)
  m      = (mu_y * prec_y + r * sum(mu)) / (prec_y + k * r)
  v      = 1 / (prec_y + k * r)

  density = function(l) { dnorm(l, m, v) }
  sampler = rnorm(1, m, v)
  list(density, sampler)
  }

dl_posterior = function(l, y, mu, r) { l_posterior(y, mu, r)[[1]](l) }
rl_posterior = function(y, mu, r) { l_posterior(y, mu, r)[[2]] }

r_obj_prior = function(y) {
  var_y  = var(unlist(y))
  prec_y = 1 / var_y

  density = function(r) { dgamma(r, 1, prec_y) }
  sampler = rgamma(1, 1, prec_y)
  list(density, sampler)
  }

dr_obj_prior = function(r, y) { r_obj_prior(y)[[1]](r) }
rr_obj_prior = function(y) { r_obj_prior(y)[[2]] }

r_posterior = function(y, mu, l) {
  var_y = var(unlist(y))
  k     = length(y)
  a     = k + 1
  b     = a / (var_y + sum((mu - l)^2))

  density = function(r) { dgamma(r, a, b) }
  sampler = rgamma(1, a, b)
  list(density, sampler)
  }

dr_posterior = function(r, y, mu, l) { r_posterior(y, mu, l)[[1]](r) }
rr_posterior = function(y, mu, l) { r_posterior(y, mu, l)[[2]] }

mu_prior = function(l, r) {
  density = function(mu) { dnorm(mu, l, 1 / r) }
  sampler = function(k) { rnorm(k, l, 1 / r) }
  list(density, sampler)
  }

dmu_prior = function(mu, l, r) { mu_prior(l, r)[[1]](mu) }
rmu_prior = function(k, l, r) { mu_prior(l, r)[[2]](k) }

mu_posterior = function(y, s, l, r) {
  n    = unlist(lapply(y, length))
  ybar = unlist(lapply(y, mean))
  m    = (ybar * n * s + l * r) / (n * s + r)
  v    = 1 / (n * s + r)

  density = function(x) { dnorm(x, m, v) }
  sampler = rnorm(length(y), m, v)
  list(density, sampler)
  }

dmu_posterior = function(mu, y, s, l, r) { mu_posterior(y, s, l, r)[[1]](mu) }
rmu_posterior = function(y, s, l, r) { mu_posterior(y, s, l, r)[[2]] }

# s (component precisions) and parameter models

s_prior = function(b, w) {
  density = function(s) { dgamma(s, b, 1 / w) }
  sampler = function(k) { rgamma(k, b, 1 / w) }
  list(density, sampler)
  }

ds_prior = function(s, b, w) { s_prior(b, w)[[1]](s) }
rs_prior = function(k, b, w) { s_prior(b, w)[[2]](k) }

s_posterior = function(y, mu, b, w) {
  n       = unlist(lapply(y, length))
  squares = unlist(mapply("-", y, as.list(mu))) ^ 2
  a       = b + n
  bet     = a / (w * b + sum(squares))

  density = function(s) { dgamma(s, a, bet) }
  sampler = rgamma(1, a, bet)
  list(density, sampler)
  }

ds_posterior = function(s, y, mu, b, w) { s_posterior(y, mu, b, w)[[1]](s) }
rs_posterior = function(y, mu, b, w) { s_posterior(y, mu, b, w)[[2]] }

w_obj_prior = function(y) {
  var_y = var(unlist(y))
  density = function(w) { dgamma(w, 1, var_y) }
  sampler = rgamma(1, 1, var_y)
  list(density, sampler)
  }

dw_obj_prior = function(w, y) { w_obj_prior(y)[[1]](w) }
rw_obj_prior = function(y) { w_obj_prior(y)[[2]] }

w_posterior = function(y, s, b) {
  k      = length(y)
  var_y  = var(unlist(y))
  prec_y = 1 / var_y
  a      = k * b + 1
  bet    = a / (prec_y + b * sum(s))

  density = function(w) { dgamma(w, a, bet) }
  sampler = dgamma(1, a, bet)
  list(density, sampler)
  }

dw_posterior = function(w, y, s, b) { w_posterior(y, s, b)[[1]](w) }
rw_posterior = function(y, s, b) { w_posterior(y, s, b)[[2]] }

b_prior = function() {
  density = function(b) { dgamma(1 / b, 1, 1) }
  sampler = 1 / rgamma(1, 1, 1)
  list(density, sampler)
  }

db_prior = function(b) { b_prior()[[1]](b) }
rb_prior = function() { b_prior()[[2]] }

b_posterior = function(y, s, w) {
  k  = length(y)

  density = function(b) {
    t0 = gamma(b / 2) ^ (- k)
    t1 = exp(-1 / (2 * b))
    t2 = (b / 2) ^ ((k * b - 3) / 2)
    t3 = prod((s * w) ^ (b / 2) * exp(- b * s * w / 2))

    t0 + t1 + t2 + t3
    }

  sampler = NULL
  list(density, sampler)
  }

db_posterior = function(b, y, s, w) { b_posterior(y, s, w)[[1]](b) }
rb_posterior = function(y, s, w) { b_posterior(y, s, w)[[2]] }

# p (mixing probabilities), n (occupation numbers) and parameter models
# note that the distribution of cluster labels can be written in terms of the
#  dirichlet concentration parameter

p_prior = function(a) {
  density = function(p) { ddirichlet(p, a) }
  sampler = rdirichlet(1, a)
  list(density, sampler)
  }

dp_prior = function(p, a) { p_prior(a)[[1]](p) }
rp_prior = function(a) { p_prior(a)[[2]] }

n_prior = function(p) {
  density = function(n) { dmultinom(n, prob = p) }
  sampler = rmultinom(1, size = 1, prob = p)
  list(density, sampler)
  }

dn_prior = function(n, p) { n_prior(p)[[1]](n) }
rn_prior = function(p) { n_prior(p)[[2]] }

a_prior = function() {
  density = function(a) { dgamma(1 / a, 1, 1) }
  sampler = 1 / rgamma(1, 1, 1)
  list(density, sampler)
  }

da_prior = function(a) { a_prior()[[1]](a) }
ra_prior = function() { a_prior()[[2]] }

a_posterior = function(k, n) {
  density = function(a) {
    a ^ (k - 3 / 2) * exp(-1 / (2 * a)) * gamma(a) / gamma(n + a)
    }

  sampler = NULL
  list(density, sampler)
  }

da_posterior = function(a, k, n) { a_posterior(k, n)[[1]](a) }

# model

model_prior = function(y) {
  k  = length(y)
  l  = rl_obj_prior(y)
  r  = rr_obj_prior(y)
  mu = rmu_prior(k, l, r)

  b = rb_prior()
  w = rw_obj_prior(y)
  s = rs_prior(k, b, w)

  a = ra_prior()
  p = rp_prior(rep(a, k))

  ry_prior(mu, s, p)
  }
