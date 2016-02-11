set.seed(42)

require(actuar)
require(ars)
require(gtools)

#  finite gaussian mixture model
#
#  a               ~ inverse-gamma(1, 1)
#  p | a           ~ symmetric-dirichlet(a)
#  c | p           ~ multinomial(p)
#  l               ~ gaussian(mu_y, var_y)
#  r               ~ gamma(1, prec_y)
#  mu | c, l, r    ~ gaussian(l, 1 / r)
#  s  | c, b, w    ~ gamma(b, 1 / w)
#  b               ~ inverse-gamma(1, 1)
#  w               ~ gamma(1, var_y)
#  y | p, c, mu, s ~ <p, normal(mu, 1 / s)>

test_data = list(
    rnorm(801, 3.5, 1)
  , rnorm(300, 0.3, 0.8)
  , rnorm(722, -4.2, 0.5)
  )

# observations

dy_prior = function(y, mu, s, p) sum(p * dnorm(y, mu, 1 / s))

ry_prior = function(mu, s, p) {
   c   = drop(rmultinom(1, 1, prob = p))
   muc = drop(c %*% mu)
   sc  = drop(c %*% s)
   rnorm(1, muc, 1 / sc)
  }

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

dmu_prior = function(mu, l, r) dnorm(mu, l, 1 / r)
rmu_prior = function(k, l, r) rnorm(k, l, 1 / r)

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

ds_prior = function(s, b, w) dgamma(s, b, 1 / w)
rs_prior = function(k, b, w) rgamma(k, b, 1 / w)

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

db_prior = function(b) dinvgamma(b, 1, 1)
rb_prior = function() rinvgamma(1, 1, 1)

b_posterior = function(y, s, w) {
  k = length(y)

  log_density = function(b) {
    -k * lgamma(b / 2) -
    (1 / (2 * b)) +
    ((k * b - 3) / 2) * log(b / 2) +
    b / 2 * sum((log(s * w)) - w * s)
    }

  density = function(b) { exp(log_density(b)) }

  sampler = function() {

    # CHECK do i really need to transform this?  is it log-concave already?
    density_xformed     = function(b) { b * density(b) }
    log_density_xformed = function(b) { log(density_xformed(b)) }

    grad_log_density_xformed = function(b) {
        1 / b - k * digamma(b / 2) +
        1 / (2 * b ^ 2) +
        k / 2 * log(b / 2) + (k * b - 3) / b +
        sum((log(s * w)) - w * s) / 2
      }

    # error due to non-log-concavity detection
    val = ars(
        n = 1
      , f = log_density_xformed
      , fprima = grad_log_density_xformed
      , x = c(0.01, 0.5, 1)
      , lb = T
      , xlb = 0.01
      , ub = T
      , xub = 10)

      exp(sample(val, 1))
  }

  list(density, sampler)
  }

db_posterior = function(b, y, s, w) { b_posterior(y, s, w)[[1]](b) }
rb_posterior = function(y, s, w) { b_posterior(y, s, w)[[2]]() }

# p (mixing probabilities), n (occupation numbers) and parameter models

dp_prior = function(p, a) ddirichlet(p, a)
rp_prior = function(a) drop(rdirichlet(1, a))

dn_prior = function(n, p) dmultinom(n, prob = p)
rn_prior = function(p) drop(rmultinom(1, size = 1, prob = p))

da_prior = function(a) dinvgamma(a, 1, 1)
ra_prior = rinvgamma(1, 1, 1)

a_posterior = function(k, n) {
  log_density = function(a) {
    (k - 3 / 2) * log(a) - 1 / (2 * a) + lgamma(a) - lgamma(n + a)
    }

  density = function(a) { exp(log_density(a)) }

  # log(a) ~ log-concave and can thus be sampled via ARS
  #
  # note that f(log(a)) = a * f(a), or log(f(log(a))) = log a + log f(a)
  #
  # ARS flakes out when k / n is sufficiently small (i.e. when 'a' would be
  # extremely large).  probably rare for this to happen in practice, so 'a' is
  # capped at 20 here
  sampler = function() {
    log_xformed = expression(
      (k - 1 / 2) * log(a) - 1 / (2 * a) + lgamma(a) - lgamma(n + a))
    dens_log_xformed = function(a) eval(log_xformed)

    grad_log_xformed = D(log_xformed, 'a')
    grad_dens_log_xformed = function(a) eval(grad_log_xformed)

    val = ars(
        n = 1
      , f = dens_log_xformed
      , fprima = grad_dens_log_xformed
      , x = c(0.51, 1, 1.49)
      , lb = T
      , xlb = 0.001
      , ub = T
      , xub = 20.0)

    exp(val)
    }
  list(density, sampler)
  }

da_posterior = function(a, k, n) { a_posterior(k, n)[[1]](a) }
ra_posterior = function(k, n) { a_posterior(k, n)[[2]]() }

# model

model_prior = function(y) {
  k  = length(y)
  l  = rl_obj_prior(y)
  r  = rr_obj_prior(y)
  mu = rmu_prior(k, l, r)
  b  = rb_prior()
  w  = rw_obj_prior(y)
  s  = rs_prior(k, b, w)
  a  = ra_prior()
  p  = rp_prior(rep(a, k))

  sampler = ry_prior(mu, s, p)
  conditional_sampler = function(n) { replicate(n, ry_prior(mu, s, p)()) }
  list(sampler, conditional_sampler)
  }

rmodel_prior = function(y) { model_prior(y)[[1]] }

# conditional model prior; sample a cluster configuration and then sample many
# observations from that configuration

rmodel_conditional_prior = function(n, y) { model_prior(y)[[2]](n) }

# debug

# y = test_data
# k = length(y)
# l = rl_obj_prior(y)
# r = rr_obj_prior(y)
# mu = rmu_prior(k, l, r)
# w = rw_obj_prior(y)
# b = rb_prior()
# s = rs_prior(k, b, w)


