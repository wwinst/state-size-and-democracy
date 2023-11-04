## TSLS with vcovPL

mediate_tsls_PL=function (model.m, model.y, treat = "treat.name", conf.level = 0.95, 
          robustSE = FALSE, cluster = NULL, boot = FALSE, sims = 1000, 
          est_se = TRUE, ...) 
{
  if (!inherits(model.m, "lm") | !inherits(model.y, "lm")) 
    stop("both mediator and outcome models must be of class `lm'.")
  m_var <- all.vars(formula(model.m)[[2]])
  y_var <- all.vars(formula(model.y)[[2]])
  t_var <- treat
  if (length(y_var) > 1L || length(m_var) > 1L) 
    stop("Left-hand side of model must only have one variable.")
  n_y <- nobs(model.y)
  n_m <- nobs(model.m)
  if (n_y != n_m) 
    stop("number of observations in both models must be identical.")
  if (!is.null(cluster)) {
    if (NROW(cluster) != n_y) 
      stop("length of `cluster' must be equal to number of observations in models.")
  }
  else {
    cluster <- seq(n_y)
  }
  .dat <- eval(getCall(model.y)$data)
  .dat <- .dat[names(model.m$fitted.values), ]
  .dat[[m_var]] <- predict(model.m)
  mod.y <- my_update(model.y, data = .dat)
  d <- coef(mod.y)[m_var] * coef(model.m)[t_var]
  z <- coef(mod.y)[t_var]
  tau.coef <- d + z
  nu <- d/tau.coef
  if (!est_se) {
    se_d <- se_z <- se_tau <- se_n <- NA
    d.ci <- z.ci <- tau.ci <- n.ci <- NA
    d.p <- z.p <- tau.p <- n.p <- NA
  }
  else {
    if (!boot) {
      sims <- NA
      if (!is.null(cluster)) {
        vcv_y <- sandwich::vcovCL(mod.y, cluster = cluster, 
                                  ...)
        vcv_m <- sandwich::vcovCL(model.m, cluster = cluster, 
                                  ...)
      }
      else if (robustSE) {
        vcv_y <- sandwich::vcovPL(mod.y, aggregate=FALSE,cluster=~country,order.by=~year)
        vcv_m <- sandwich::vcovPL(model.m, aggregate=FALSE,cluster=~country,order.by=~year)
      }
      else {
        vcv_y <- vcov(mod.y)
        vcv_m <- vcov(model.m)
      }
      se_d <- sqrt(coef(mod.y)[m_var]^2 * vcv_m[t_var, 
                                                t_var] + coef(model.m)[t_var]^2 * vcv_y[m_var, 
                                                                                        m_var] + vcv_m[t_var, t_var] * vcv_y[m_var, 
                                                                                                                             m_var])
      se_z <- sqrt(vcv_y[t_var, t_var])
      se_tau <- sqrt(vcv_y[t_var, t_var] + (se_d)^2 + 
                       2 * vcv_y[t_var, m_var] * coef(model.m)[t_var])
      delta <- function(f, B, Sigma) {
        ff <- deriv(f, names(B), func = TRUE)
        x <- do.call(ff, as.list(B))
        grad <- as.matrix(attr(x, "gradient"), nr = 1)
        sqrt(grad %*% Sigma %*% t(grad))
      }
      Coefs <- c(coef(model.m)[t_var], coef(mod.y)[t_var], 
                 coef(mod.y)[m_var])
      Coefs <- setNames(Coefs, c("b2", "b3", "gamma"))
      Sigma <- diag(c(vcv_m[t_var, t_var], diag(vcv_y)[c(t_var, 
                                                         m_var)]))
      Sigma[3, 2] <- Sigma[2, 3] <- vcv_y[t_var, m_var]
      f <- ~b2 * gamma/(b2 * gamma + b3)
      se_n <- as.vector(delta(f, Coefs, Sigma))
      qq <- (1 - conf.level)/2
      qq <- setNames(c(qq, 1 - qq), c("low", "high"))
      d.ci <- d + qnorm(qq) * se_d
      z.ci <- z + qnorm(qq) * se_z
      tau.ci <- tau.coef + qnorm(qq) * se_tau
      n.ci <- nu + qnorm(qq) * se_n
      d.p <- pnorm(-abs(d), sd = se_d)
      z.p <- pnorm(-abs(z), sd = se_z)
      tau.p <- pnorm(-abs(tau.coef), sd = se_tau)
      n.p <- pnorm(-abs(nu), sd = se_n)
    }
    else {
      cl <- split(seq_along(cluster), cluster)
      cf <- matrix(rep.int(0, 4 * sims), ncol = 4, dimnames = list(NULL, 
                                                                   c("delta", "zeta", "tau", "nu")))
      for (i in 1:sims) {
        .subset <- unlist(cl[sample(names(cl), length(cl), 
                                    replace = TRUE)])
        .dat_y <- eval(getCall(model.y)$data)[.subset, 
        ]
        .dat_m <- eval(getCall(model.m)$data)[.subset, 
        ]
        out <- tryCatch({
          up_y <- my_update(model.y, data = .dat_y)
          up_m <- my_update(model.m, data = .dat_m)
          mediate_tsls(up_m, up_y, treat = treat, cluster = NULL, 
                       est_se = FALSE)[c("d1", "z0", "tau.coef", 
                                         "n0")]
        }, error = function(e) {
          setNames(rep(list(NA), 4), c("d1", "z0", "tau.coef", 
                                       "n0"))
        })
        cf[i, ] <- unlist(out)
      }
      se_d <- sd(cf[, "delta"], na.rm = TRUE)
      se_z <- sd(cf[, "zeta"], na.rm = TRUE)
      se_tau <- sd(cf[, "tau"], na.rm = TRUE)
      se_n <- sd(cf[, "nu"], na.rm = TRUE)
      qq <- (1 - conf.level)/2
      qq <- setNames(c(qq, 1 - qq), c("low", "high"))
      d.ci <- quantile(cf[, "delta"], qq, na.rm = TRUE)
      z.ci <- quantile(cf[, "zeta"], qq, na.rm = TRUE)
      tau.ci <- quantile(cf[, "tau"], qq, na.rm = TRUE)
      n.ci <- quantile(cf[, "nu"], qq, na.rm = TRUE)
      d.p <- pval(cf[, "delta"], d)
      z.p <- pval(cf[, "zeta"], z)
      tau.p <- pval(cf[, "tau"], tau.coef)
      n.p <- pval(cf[, "nu"], nu)
    }
  }
  out <- list(d1 = unname(d), d1.se = se_d, d1.p = d.p, d1.ci = d.ci, 
              d0 = unname(d), d0.se = se_d, d0.p = d.p, d0.ci = d.ci, 
              z1 = unname(z), z1.se = se_z, z1.p = z.p, z1.ci = z.ci, 
              z0 = unname(z), z0.se = se_z, z0.p = z.p, z0.ci = z.ci, 
              tau.coef = unname(tau.coef), tau.se = se_tau, tau.ci = tau.ci, 
              tau.p = tau.p, n0 = unname(nu), n0.se = se_n, n0.ci = n.ci, 
              n0.p = n.p, boot = boot, boot.ci.type = "perc", treat = treat, 
              mediator = m_var, nobs = nobs(model.y), sims = sims, 
              INT = FALSE, conf.level = conf.level, model.y = model.y, 
              model.m = model.m)
  class(out) <- c("mediate", "mediate.tsls")
  return(out)
}

mediate_tsls_SCC=function (model.m, model.y, treat = "treat.name", conf.level = 0.95, 
                          robustSE = FALSE, cluster = NULL, boot = FALSE, sims = 1000, 
                          est_se = TRUE, ...) 
{
  if (!inherits(model.m, "lm") | !inherits(model.y, "lm")) 
    stop("both mediator and outcome models must be of class `lm'.")
  m_var <- all.vars(formula(model.m)[[2]])
  y_var <- all.vars(formula(model.y)[[2]])
  t_var <- treat
  if (length(y_var) > 1L || length(m_var) > 1L) 
    stop("Left-hand side of model must only have one variable.")
  n_y <- nobs(model.y)
  n_m <- nobs(model.m)
  if (n_y != n_m) 
    stop("number of observations in both models must be identical.")
  if (!is.null(cluster)) {
    if (NROW(cluster) != n_y) 
      stop("length of `cluster' must be equal to number of observations in models.")
  }
  else {
    cluster <- seq(n_y)
  }
  .dat <- eval(getCall(model.y)$data)
  .dat <- .dat[names(model.m$fitted.values), ]
  .dat[[m_var]] <- predict(model.m)
  mod.y <- my_update(model.y, data = .dat)
  d <- coef(mod.y)[m_var] * coef(model.m)[t_var]
  z <- coef(mod.y)[t_var]
  tau.coef <- d + z
  nu <- d/tau.coef
  if (!est_se) {
    se_d <- se_z <- se_tau <- se_n <- NA
    d.ci <- z.ci <- tau.ci <- n.ci <- NA
    d.p <- z.p <- tau.p <- n.p <- NA
  }
  else {
    if (!boot) {
      sims <- NA
      if (!is.null(cluster)) {
        vcv_y <- sandwich::vcovCL(mod.y, cluster = cluster, 
                                  ...)
        vcv_m <- sandwich::vcovCL(model.m, cluster = cluster, 
                                  ...)
      }
      else if (robustSE) {
        vcv_y <- sandwich::vcovPL(mod.y, aggregate=TRUE,cluster=~country,order.by=~year)
        vcv_m <- sandwich::vcovPL(model.m, aggregate=TRUE,cluster=~country,order.by=~year)
      }
      else {
        vcv_y <- vcov(mod.y)
        vcv_m <- vcov(model.m)
      }
      se_d <- sqrt(coef(mod.y)[m_var]^2 * vcv_m[t_var, 
                                                t_var] + coef(model.m)[t_var]^2 * vcv_y[m_var, 
                                                                                        m_var] + vcv_m[t_var, t_var] * vcv_y[m_var, 
                                                                                                                             m_var])
      se_z <- sqrt(vcv_y[t_var, t_var])
      se_tau <- sqrt(vcv_y[t_var, t_var] + (se_d)^2 + 
                       2 * vcv_y[t_var, m_var] * coef(model.m)[t_var])
      delta <- function(f, B, Sigma) {
        ff <- deriv(f, names(B), func = TRUE)
        x <- do.call(ff, as.list(B))
        grad <- as.matrix(attr(x, "gradient"), nr = 1)
        sqrt(grad %*% Sigma %*% t(grad))
      }
      Coefs <- c(coef(model.m)[t_var], coef(mod.y)[t_var], 
                 coef(mod.y)[m_var])
      Coefs <- setNames(Coefs, c("b2", "b3", "gamma"))
      Sigma <- diag(c(vcv_m[t_var, t_var], diag(vcv_y)[c(t_var, 
                                                         m_var)]))
      Sigma[3, 2] <- Sigma[2, 3] <- vcv_y[t_var, m_var]
      f <- ~b2 * gamma/(b2 * gamma + b3)
      se_n <- as.vector(delta(f, Coefs, Sigma))
      qq <- (1 - conf.level)/2
      qq <- setNames(c(qq, 1 - qq), c("low", "high"))
      d.ci <- d + qnorm(qq) * se_d
      z.ci <- z + qnorm(qq) * se_z
      tau.ci <- tau.coef + qnorm(qq) * se_tau
      n.ci <- nu + qnorm(qq) * se_n
      d.p <- pnorm(-abs(d), sd = se_d)
      z.p <- pnorm(-abs(z), sd = se_z)
      tau.p <- pnorm(-abs(tau.coef), sd = se_tau)
      n.p <- pnorm(-abs(nu), sd = se_n)
    }
    else {
      cl <- split(seq_along(cluster), cluster)
      cf <- matrix(rep.int(0, 4 * sims), ncol = 4, dimnames = list(NULL, 
                                                                   c("delta", "zeta", "tau", "nu")))
      for (i in 1:sims) {
        .subset <- unlist(cl[sample(names(cl), length(cl), 
                                    replace = TRUE)])
        .dat_y <- eval(getCall(model.y)$data)[.subset, 
        ]
        .dat_m <- eval(getCall(model.m)$data)[.subset, 
        ]
        out <- tryCatch({
          up_y <- my_update(model.y, data = .dat_y)
          up_m <- my_update(model.m, data = .dat_m)
          mediate_tsls(up_m, up_y, treat = treat, cluster = NULL, 
                       est_se = FALSE)[c("d1", "z0", "tau.coef", 
                                         "n0")]
        }, error = function(e) {
          setNames(rep(list(NA), 4), c("d1", "z0", "tau.coef", 
                                       "n0"))
        })
        cf[i, ] <- unlist(out)
      }
      se_d <- sd(cf[, "delta"], na.rm = TRUE)
      se_z <- sd(cf[, "zeta"], na.rm = TRUE)
      se_tau <- sd(cf[, "tau"], na.rm = TRUE)
      se_n <- sd(cf[, "nu"], na.rm = TRUE)
      qq <- (1 - conf.level)/2
      qq <- setNames(c(qq, 1 - qq), c("low", "high"))
      d.ci <- quantile(cf[, "delta"], qq, na.rm = TRUE)
      z.ci <- quantile(cf[, "zeta"], qq, na.rm = TRUE)
      tau.ci <- quantile(cf[, "tau"], qq, na.rm = TRUE)
      n.ci <- quantile(cf[, "nu"], qq, na.rm = TRUE)
      d.p <- pval(cf[, "delta"], d)
      z.p <- pval(cf[, "zeta"], z)
      tau.p <- pval(cf[, "tau"], tau.coef)
      n.p <- pval(cf[, "nu"], nu)
    }
  }
  out <- list(d1 = unname(d), d1.se = se_d, d1.p = d.p, d1.ci = d.ci, 
              d0 = unname(d), d0.se = se_d, d0.p = d.p, d0.ci = d.ci, 
              z1 = unname(z), z1.se = se_z, z1.p = z.p, z1.ci = z.ci, 
              z0 = unname(z), z0.se = se_z, z0.p = z.p, z0.ci = z.ci, 
              tau.coef = unname(tau.coef), tau.se = se_tau, tau.ci = tau.ci, 
              tau.p = tau.p, n0 = unname(nu), n0.se = se_n, n0.ci = n.ci, 
              n0.p = n.p, boot = boot, boot.ci.type = "perc", treat = treat, 
              mediator = m_var, nobs = nobs(model.y), sims = sims, 
              INT = FALSE, conf.level = conf.level, model.y = model.y, 
              model.m = model.m)
  class(out) <- c("mediate", "mediate.tsls")
  return(out)
}

mediate_tsls_PCSE=function (model.m, model.y, treat = "treat.name", conf.level = 0.95, 
                          robustSE = FALSE, cluster = NULL, boot = FALSE, sims = 1000, 
                          est_se = TRUE, ...) 
{
  if (!inherits(model.m, "lm") | !inherits(model.y, "lm")) 
    stop("both mediator and outcome models must be of class `lm'.")
  m_var <- all.vars(formula(model.m)[[2]])
  y_var <- all.vars(formula(model.y)[[2]])
  t_var <- treat
  if (length(y_var) > 1L || length(m_var) > 1L) 
    stop("Left-hand side of model must only have one variable.")
  n_y <- nobs(model.y)
  n_m <- nobs(model.m)
  if (n_y != n_m) 
    stop("number of observations in both models must be identical.")
  if (!is.null(cluster)) {
    if (NROW(cluster) != n_y) 
      stop("length of `cluster' must be equal to number of observations in models.")
  }
  else {
    cluster <- seq(n_y)
  }
  .dat <- eval(getCall(model.y)$data)
  .dat <- .dat[names(model.m$fitted.values), ]
  .dat[[m_var]] <- predict(model.m)
  mod.y <- my_update(model.y, data = .dat)
  d <- coef(mod.y)[m_var] * coef(model.m)[t_var]
  z <- coef(mod.y)[t_var]
  tau.coef <- d + z
  nu <- d/tau.coef
  if (!est_se) {
    se_d <- se_z <- se_tau <- se_n <- NA
    d.ci <- z.ci <- tau.ci <- n.ci <- NA
    d.p <- z.p <- tau.p <- n.p <- NA
  }
  else {
    if (!boot) {
      sims <- NA
      if (!is.null(cluster)) {
        vcv_y <- sandwich::vcovCL(mod.y, cluster = cluster, 
                                  ...)
        vcv_m <- sandwich::vcovCL(model.m, cluster = cluster, 
                                  ...)
      }
      else if (robustSE) {
        vcv_y <- sandwich::vcovPC(mod.y, pairwise=TRUE,cluster=~country,order.by=~year)
        vcv_m <- sandwich::vcovPC(model.m, pairwise=TRUE,cluster=~country,order.by=~year)
      }
      else {
        vcv_y <- vcov(mod.y)
        vcv_m <- vcov(model.m)
      }
      se_d <- sqrt(coef(mod.y)[m_var]^2 * vcv_m[t_var, 
                                                t_var] + coef(model.m)[t_var]^2 * vcv_y[m_var, 
                                                                                        m_var] + vcv_m[t_var, t_var] * vcv_y[m_var, 
                                                                                                                             m_var])
      se_z <- sqrt(vcv_y[t_var, t_var])
      se_tau <- sqrt(vcv_y[t_var, t_var] + (se_d)^2 + 
                       2 * vcv_y[t_var, m_var] * coef(model.m)[t_var])
      delta <- function(f, B, Sigma) {
        ff <- deriv(f, names(B), func = TRUE)
        x <- do.call(ff, as.list(B))
        grad <- as.matrix(attr(x, "gradient"), nr = 1)
        sqrt(grad %*% Sigma %*% t(grad))
      }
      Coefs <- c(coef(model.m)[t_var], coef(mod.y)[t_var], 
                 coef(mod.y)[m_var])
      Coefs <- setNames(Coefs, c("b2", "b3", "gamma"))
      Sigma <- diag(c(vcv_m[t_var, t_var], diag(vcv_y)[c(t_var, 
                                                         m_var)]))
      Sigma[3, 2] <- Sigma[2, 3] <- vcv_y[t_var, m_var]
      f <- ~b2 * gamma/(b2 * gamma + b3)
      se_n <- as.vector(delta(f, Coefs, Sigma))
      qq <- (1 - conf.level)/2
      qq <- setNames(c(qq, 1 - qq), c("low", "high"))
      d.ci <- d + qnorm(qq) * se_d
      z.ci <- z + qnorm(qq) * se_z
      tau.ci <- tau.coef + qnorm(qq) * se_tau
      n.ci <- nu + qnorm(qq) * se_n
      d.p <- pnorm(-abs(d), sd = se_d)
      z.p <- pnorm(-abs(z), sd = se_z)
      tau.p <- pnorm(-abs(tau.coef), sd = se_tau)
      n.p <- pnorm(-abs(nu), sd = se_n)
    }
    else {
      cl <- split(seq_along(cluster), cluster)
      cf <- matrix(rep.int(0, 4 * sims), ncol = 4, dimnames = list(NULL, 
                                                                   c("delta", "zeta", "tau", "nu")))
      for (i in 1:sims) {
        .subset <- unlist(cl[sample(names(cl), length(cl), 
                                    replace = TRUE)])
        .dat_y <- eval(getCall(model.y)$data)[.subset, 
        ]
        .dat_m <- eval(getCall(model.m)$data)[.subset, 
        ]
        out <- tryCatch({
          up_y <- my_update(model.y, data = .dat_y)
          up_m <- my_update(model.m, data = .dat_m)
          mediate_tsls(up_m, up_y, treat = treat, cluster = NULL, 
                       est_se = FALSE)[c("d1", "z0", "tau.coef", 
                                         "n0")]
        }, error = function(e) {
          setNames(rep(list(NA), 4), c("d1", "z0", "tau.coef", 
                                       "n0"))
        })
        cf[i, ] <- unlist(out)
      }
      se_d <- sd(cf[, "delta"], na.rm = TRUE)
      se_z <- sd(cf[, "zeta"], na.rm = TRUE)
      se_tau <- sd(cf[, "tau"], na.rm = TRUE)
      se_n <- sd(cf[, "nu"], na.rm = TRUE)
      qq <- (1 - conf.level)/2
      qq <- setNames(c(qq, 1 - qq), c("low", "high"))
      d.ci <- quantile(cf[, "delta"], qq, na.rm = TRUE)
      z.ci <- quantile(cf[, "zeta"], qq, na.rm = TRUE)
      tau.ci <- quantile(cf[, "tau"], qq, na.rm = TRUE)
      n.ci <- quantile(cf[, "nu"], qq, na.rm = TRUE)
      d.p <- pval(cf[, "delta"], d)
      z.p <- pval(cf[, "zeta"], z)
      tau.p <- pval(cf[, "tau"], tau.coef)
      n.p <- pval(cf[, "nu"], nu)
    }
  }
  out <- list(d1 = unname(d), d1.se = se_d, d1.p = d.p, d1.ci = d.ci, 
              d0 = unname(d), d0.se = se_d, d0.p = d.p, d0.ci = d.ci, 
              z1 = unname(z), z1.se = se_z, z1.p = z.p, z1.ci = z.ci, 
              z0 = unname(z), z0.se = se_z, z0.p = z.p, z0.ci = z.ci, 
              tau.coef = unname(tau.coef), tau.se = se_tau, tau.ci = tau.ci, 
              tau.p = tau.p, n0 = unname(nu), n0.se = se_n, n0.ci = n.ci, 
              n0.p = n.p, boot = boot, boot.ci.type = "perc", treat = treat, 
              mediator = m_var, nobs = nobs(model.y), sims = sims, 
              INT = FALSE, conf.level = conf.level, model.y = model.y, 
              model.m = model.m)
  class(out) <- c("mediate", "mediate.tsls")
  return(out)
}
