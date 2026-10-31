# This script runs prior-sensitivity and non-ALD error experiments.
# Results are saved as CSV files in the output directory.

suppressPackageStartupMessages({
  library(quantreg)
  library(GIGrvg)
  library(MASS)
  library(coda)
})

set.seed(2026)

options(digits = 4, width = 180)

# Output is written to a folder under the current working directory.
OUTPUT_DIR <- file.path(getwd(), "BFLQR_Supplementary_Results")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(OUTPUT_DIR)) stop("Could not create output directory: ", OUTPUT_DIR)

save_csv <- function(x, filename) {
  if (!is.null(x) && nrow(x) > 0L) {
    utils::write.csv(x, file.path(OUTPUT_DIR, filename), row.names = FALSE, na = "")
  }
  invisible(NULL)
}

# 1. Simulation settings
N    <- c(50, 100, 200)
TAU  <- c(0.25, 0.50, 0.75)
R    <- 1               
ITER <- 12000
BURN <- 6000
NTEST <- 5000               
p    <- 10
rho  <- 0.5

RUN_PRIOR_SENSITIVITY <- TRUE
RUN_NON_ALD_EXPERIMENT <- TRUE
R_SUPP <- R
if (length(R_SUPP) != 1L || !is.finite(R_SUPP) || R_SUPP < 1 ||
    R_SUPP != as.integer(R_SUPP)) {
  stop("R_SUPP must be a single positive integer.")
}
# Prior sensitivity settings: baseline Gamma(1,1) hyperpriors versus Gamma(0.5,0.5)
# for the global shrinkage precisions. This is a sensitivity scenario
# that either prior is universally preferable.
PRIOR_ALT <- list(a1 = 0.5, b1 = 0.5, a2 = 0.5, b2 = 0.5,
                  a_sigma = 1, b_sigma = 1)

b0 <- 1.5
bt <- c(2, 2, 1.5, 0, 0, -1.5, -1.5, 0, 1, 1)
true_active <- bt != 0
true_zero   <- bt == 0
S <- outer(seq_len(p), seq_len(p), function(i, j) rho^abs(i - j))

# Hyperparameters and INITIAL log-scale random-walk proposal SDs
# (SDs are adapted during burn-in, then frozen)
a1 <- 1; b1 <- 1
a2 <- 1; b2 <- 1
a_sigma <- 1; b_sigma <- 1
MH_SD_TAU   <- 0.8
MH_SD_OMEGA <- 0.8
MH_TARGET_ACC <- 0.35 

# 2. Data generation and contamination
rALD <- function(n, tau, sigma = 1) {
  phi1 <- (1 - 2 * tau) / (tau * (1 - tau))
  phi2 <- 2 / (tau * (1 - tau))
  v <- rexp(n, rate = 1 / sigma)
  z <- rnorm(n)
  phi1 * v + sqrt(phi2 * sigma * v) * z
}

r_error <- function(n, tau, error_dist = "ALD") {
  if (identical(error_dist, "ALD")) return(rALD(n, tau))
  if (identical(error_dist, "t3_shifted")) {
    # Heavy-tailed alternative. Shift t_3 errors so their tau-quantile is zero,
    # preserving the target conditional quantile b0 + X beta.
    return(rt(n, df = 3) - qt(tau, df = 3))
  }
  stop("Unknown error_dist: ", error_dist)
}

gendata <- function(n, tau, error_dist = "ALD") {
  X <- matrix(rnorm(n * p), nrow = n, ncol = p) %*% chol(S)
  y <- as.numeric(b0 + X %*% bt + r_error(n, tau, error_dist))
  list(X = X, y = y)
}

contaminate <- function(y) {
  m <- max(1L, round(0.10 * length(y)))
  id <- sample(seq_along(y), m, replace = FALSE)
  y[id] <- y[id] + 8 * sd(y) * sample(c(-1, 1), m, replace = TRUE)
  y
}

QL <- function(y, pred, tau) {
  u <- as.numeric(y - pred)
  mean(ifelse(u >= 0, tau * u, (tau - 1) * u))
}

# First-difference matrix (p-1) x p; used by FLQR and BFLQR
make_D <- function(pp) {
  D <- matrix(0, nrow = pp - 1L, ncol = pp)
  for (j in seq_len(pp - 1L)) {
    D[j, j] <- -1
    D[j, j + 1L] <- 1
  }
  D
}

# 3. Frequentist comparator models
# NOTE: both rq.fit.lasso and the FLQR objective use the SUMMED check loss,
# so penalty grids are scaled by n. In CV, lambda is rescaled by n_train / n
# so that the per-observation penalty strength matches the final fit.

# Lasso quantile regression: select lambda using held-out 5-fold pinball loss.
fitLQR <- function(X, y, tau, nfolds = 5L, lambda_grid = NULL) {
  n <- nrow(X)
  if (is.null(lambda_grid))
    lambda_grid <- n * exp(seq(log(1e-3), log(0.5), length.out = 30))
  fold_id <- sample(rep(seq_len(min(nfolds, n)), length.out = n))
  Z <- cbind(1, X)
  cv_loss <- matrix(NA_real_, nrow = length(lambda_grid),
                    ncol = max(fold_id))

  for (k in seq_len(max(fold_id))) {
    test_id <- which(fold_id == k)
    train_id <- which(fold_id != k)
    scale_k <- length(train_id) / n
    for (j in seq_along(lambda_grid)) {
      penalty <- c(0, rep(lambda_grid[j] * scale_k, ncol(X))) # intercept unpenalized
      fit <- rq.fit.lasso(Z[train_id, , drop = FALSE],
                          y[train_id], tau = tau,
                          lambda = penalty)
      pred <- as.numeric(Z[test_id, , drop = FALSE] %*% fit$coefficients)
      cv_loss[j, k] <- QL(y[test_id], pred, tau)
    }
  }

  mean_cv <- rowMeans(cv_loss)
  best_lambda <- lambda_grid[which.min(mean_cv)]
  penalty <- c(0, rep(best_lambda, ncol(X))) # intercept unpenalized
  fit <- rq.fit.lasso(Z, y, tau = tau, lambda = penalty)
  ans <- as.numeric(fit$coefficients)
  attr(ans, "lambda") <- best_lambda
  ans
}

# FLQR: min sum rho_tau(y - b0 - X b) + l1*sum|b| + l2*sum|D b|, intercept unpenalized.
# Solved exactly as an LP. A pseudo-row (x = a, y = 0) paired with its negative
# (x = -a, y = 0) contributes tau*|a'b| + (1-tau)*|a'b| = |a'b| for ANY tau.
fitFLQR <- function(X, y, tau, nfolds = 5L,
                    l1_grid = NULL, l2_grid = NULL) {
  pp <- ncol(X)
  n <- nrow(X)
  if (is.null(l1_grid)) l1_grid <- n * c(0.002, 0.01, 0.03, 0.1, 0.3)
  if (is.null(l2_grid)) l2_grid <- n * c(0.002, 0.01, 0.03, 0.1, 0.3)
  D <- make_D(pp)

  fit_one <- function(xx, yy, l1, l2) {
    A <- rbind(l1 * diag(pp), -l1 * diag(pp), l2 * D, -l2 * D)
    Zaug <- rbind(cbind(1, xx), cbind(0, A))
    yaug <- c(yy, rep(0, nrow(A)))
    as.numeric(rq.fit.br(Zaug, yaug, tau = tau)$coefficients)
  }

  K <- min(as.integer(nfolds), n)
  fold_id <- sample(rep(seq_len(K), length.out = n))
  grid <- expand.grid(l1 = l1_grid, l2 = l2_grid)
  cv <- rep(NA_real_, nrow(grid))
  for (g in seq_len(nrow(grid))) {
    losses <- numeric(K)
    for (k in seq_len(K)) {
      itest <- which(fold_id == k)
      itrain <- which(fold_id != k)
      sc <- length(itrain) / n
      cf <- fit_one(X[itrain, , drop = FALSE], y[itrain],
                    grid$l1[g] * sc, grid$l2[g] * sc)
      pred <- as.numeric(cf[1] + X[itest, , drop = FALSE] %*% cf[-1])
      losses[k] <- QL(y[itest], pred, tau)
    }
    cv[g] <- mean(losses)
  }
  best <- which.min(cv)
  ans <- fit_one(X, y, grid$l1[best], grid$l2[best])
  attr(ans, "lambda1") <- grid$l1[best]
  attr(ans, "lambda2") <- grid$l2[best]
  ans
}

# 4. Matrix and determinant utilities for the Bayesian sampler

make_Qbeta <- function(tau2, omega2, D) {
  diag(1 / tau2, nrow = length(tau2)) +
    crossprod(D, sweep(D, 1L, omega2, FUN = function(x, w) x / w))
}

# Strict Cholesky factorization. The sampler's rank-one updates require the
# factor, inverse and determinant to refer to exactly the same matrix. Do not
# silently add jitter here: if an SPD matrix fails, stop and report the cause.
strict_chol <- function(A, label = "matrix") {
  A <- as.matrix(A)
  A <- (A + t(A)) / 2
  if (any(!is.finite(A))) {
    stop(label, " contains non-finite values.")
  }
  tryCatch(
    chol(A),
    error = function(e) {
      stop("Cholesky factorization failed for ", label,
           ". Check local scales, conditioning, and model inputs. Original error: ",
           conditionMessage(e), call. = FALSE)
    }
  )
}

logdet_spd <- function(A) {
  Rc <- strict_chol(A, "SPD matrix")
  2 * sum(log(diag(Rc)))
}

# Log targets for one log-local-scale proposal, including:
# (i) the Gaussian-prior determinant 0.5*log|Q| (supplied as `ld`, the log-determinant
#     of Q_beta WITH the proposed value, computed by a rank-1 update),
# (ii) quadratic term, (iii) exponential prior, (iv) Jacobian +s under u=exp(s).
log_tau2_target <- function(s, j, beta, sigma, lambda1sq, ld) {
  u <- exp(s)
  0.5 * ld -
    beta[j]^2 / (2 * sigma * u) -
    0.5 * lambda1sq * u + s
}

log_omega2_target <- function(s, j, beta, sigma, lambda2sq, ld, D) {
  u <- exp(s)
  delta_j <- as.numeric(D[j, , drop = FALSE] %*% beta)
  0.5 * ld -
    delta_j^2 / (2 * sigma * u) -
    0.5 * lambda2sq * u + s
}

# Exact univariate slice sampler on the log-local-scale. This avoids random-walk
# MH getting stuck while a local scale explores a highly skewed/heavy-tailed target.
# The supplied log-density must include the log-scale Jacobian.
slice_sample_1d <- function(x0, logf, width = 1, max_steps = 100L, max_shrink = 300L) {
  f0 <- logf(x0)
  if (!is.finite(f0)) stop("Non-finite slice density at current local scale.")
  logy <- f0 + log(runif(1))
  L <- x0 - runif(1, 0, width)
  R <- L + width
  j <- sample.int(max_steps + 1L, 1L) - 1L
  k <- max_steps - j
  while (j > 0L && logf(L) > logy) { L <- L - width; j <- j - 1L }
  while (k > 0L && logf(R) > logy) { R <- R + width; k <- k - 1L }
  for (iter_shrink in seq_len(max_shrink)) {
    x1 <- runif(1, L, R)
    f1 <- logf(x1)
    if (is.finite(f1) && f1 >= logy) return(x1)
    if (x1 < x0) L <- x1 else R <- x1
    if ((R - L) < 1e-12) break
  }
  # A failed slice update is an explicit error rather than a silent stale draw.
  stop("Slice sampler failed to find an acceptable point; inspect local-scale target.")
}

# 5. BFLQR Metropolis-within-Gibbs sampler
fitBFLQR <- function(X, y, tau, iter = ITER, burn = BURN,
                     mh_sd_tau = MH_SD_TAU,
                     mh_sd_omega = MH_SD_OMEGA,
                     a1 = 1, b1 = 1, a2 = 1, b2 = 1,
                     a_sigma = 1, b_sigma = 1,
                     keep_diagnostics = TRUE,
                     keep_trace = FALSE,
                     init_jitter = FALSE) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- nrow(X)
  pp <- ncol(X)

  if (length(y) != n) stop("nrow(X) must equal length(y).")
  if (pp < 2L) stop("BFLQR requires at least two predictors.")
  if (tau <= 0 || tau >= 1) stop("tau must lie strictly between 0 and 1.")
  if (length(iter) != 1L || !is.finite(iter) || iter < 2 || iter != as.integer(iter))
    stop("iter must be a single integer >= 2.")
  if (length(burn) != 1L || !is.finite(burn) || burn < 0 ||
      burn >= iter || burn != as.integer(burn))
    stop("Require burn to be an integer satisfying 0 <= burn < iter.")
  if (length(mh_sd_tau) != 1L || !is.finite(mh_sd_tau) || mh_sd_tau <= 0 ||
      length(mh_sd_omega) != 1L || !is.finite(mh_sd_omega) || mh_sd_omega <= 0)
    stop("MH proposal standard deviations must be finite and positive.")
  if (any(!is.finite(X)) || any(!is.finite(y))) stop("X and y must be finite.")

  phi1 <- (1 - 2 * tau) / (tau * (1 - tau))
  phi2 <- 2 / (tau * (1 - tau))
  Z <- cbind(1, X)
  D <- make_D(pp)

  # Start near a data-informed quantile-regression fit rather than beta=0.
  # This is only an initialization; all parameters remain sampled afterward.
  rq_init <- tryCatch(quantreg::rq.fit.br(Z, y, tau = tau)$coefficients,
                      error = function(e) NULL)
  if (is.null(rq_init) || length(rq_init) != pp + 1L || any(!is.finite(rq_init))) {
    beta0 <- as.numeric(quantile(y, probs = tau, names = FALSE))
    beta <- rep(0, pp)
  } else {
    beta0 <- as.numeric(rq_init[1L])
    beta <- as.numeric(rq_init[-1L])
  }
  if (init_jitter) {
    beta0 <- beta0 + rnorm(1, 0, 0.5 * sd(y))
    beta <- beta + rnorm(pp, 0, 0.5 * max(sd(y), 1e-3))
  }
  sigma <- max(var(y), 1e-3)
  v <- rep(max(sigma, 1e-3), n)
  tau2 <- pmax(abs(beta) * sqrt(sigma) / max(mean(abs(beta)), 0.1), 0.1)
  omega2 <- pmax(abs(diff(beta)) * sqrt(sigma) / max(mean(abs(diff(beta))), 0.1), 0.1)
  if (init_jitter) {
    tau2 <- tau2 * exp(rnorm(pp, 0, 0.8))
    omega2 <- omega2 * exp(rnorm(pp - 1L, 0, 0.8))
  }
  lambda1sq <- 1
  lambda2sq <- 1

  # Slice sampling replaces random-walk MH for local scales. The legacy mh_sd
  # arguments remain accepted for call compatibility but are not used.
  slice_width_tau <- 1
  slice_width_omega <- 1

  K <- iter - burn
  S0 <- numeric(K)
  S1 <- matrix(NA_real_, nrow = K, ncol = pp)
  Ssigma <- numeric(K)
  Slam1 <- numeric(K)
  Slam2 <- numeric(K)
  Stau2 <- if (keep_trace) matrix(NA_real_, nrow = K, ncol = pp) else NULL
  Somega2 <- if (keep_trace) matrix(NA_real_, nrow = K, ncol = pp - 1L) else NULL
  move_tau <- integer(pp)      # retained updates that changed the scale
  move_omega <- integer(pp - 1L)

  for (m in seq_len(iter)) {
    if (m %% 1000L == 0L) {
      cat(sprintf("    BFLQR MCMC iteration %d/%d (burn-in %d)\n", m, iter, burn))
      flush.console()
    }
    # Step 1: current local-scale matrices and precision Q_beta
    Qbeta <- make_Qbeta(tau2, omega2, D)
    Qhat <- matrix(0, nrow = pp + 1L, ncol = pp + 1L)
    Qhat[2:(pp + 1L), 2:(pp + 1L)] <- Qbeta

    # Step 2: multivariate-normal update for intercept and beta
    # A = Z' V^-1 Z / phi2 + Qhat; covariance = sigma * A^-1.
    rV <- 1 / pmax(v, 1e-12)
    A <- crossprod(Z, Z * rV) / phi2 + Qhat
    A <- (A + t(A)) / 2
    rhs <- as.numeric(crossprod(Z, (y - phi1 * v) * rV)) / phi2
    Rchol <- strict_chol(A, "coefficient precision matrix A")
    mu_xi <- backsolve(Rchol, forwardsolve(t(Rchol), rhs))
    # Draw N(mu, sigma*A^{-1}) using the Cholesky precision factor.
    xi <- as.numeric(mu_xi + sqrt(sigma) *
                     backsolve(Rchol, rnorm(pp + 1L)))
    beta0 <- xi[1]
    beta <- xi[-1]

    # Step 3: ALD latent variables v_i | rest ~ GIG(1/2, chi_i, psi)
    resid <- as.numeric(y - beta0 - X %*% beta)
    chi <- pmax(resid^2 / (phi2 * sigma), 1e-12)
    psi <- (phi1^2 / phi2 + 2) / sigma
    v <- vapply(seq_len(n), function(i) {
      GIGrvg::rgig(1, lambda = 0.5, chi = chi[i], psi = psi)
    }, numeric(1))

    # Cache Q_beta^{-1} and log|Q_beta| once per iteration (exact, refreshed each
    # iteration so rank-1 update drift cannot accumulate).
    Rq <- strict_chol(Qbeta, "local-scale precision matrix Q_beta")
    Qinv <- chol2inv(Rq)
    logdetQ <- 2 * sum(log(diag(Rq)))

    # Step 4: update coefficient local scales. The determinant lemma gives
    # the exact proposed log-determinant without refactorizing Q_beta at every
    # slice evaluation. After each draw, update its inverse and log-determinant.
    for (j in seq_len(pp)) {
      old_u <- tau2[j]
      s_current <- log(old_u)
      logf_tau <- function(s) {
        if (!is.finite(s) || s > 700 || s < -700) return(-Inf)
        new_u <- exp(s)
        change <- 1 / new_u - 1 / old_u
        denom <- 1 + change * Qinv[j, j]
        if (!is.finite(denom) || denom <= 0) return(-Inf)
        ld <- logdetQ + log(denom)
        log_tau2_target(s, j, beta, sigma, lambda1sq, ld)
      }
      s_new <- slice_sample_1d(s_current, logf_tau, width = slice_width_tau)
      new_u <- exp(s_new)
      change <- 1 / new_u - 1 / old_u
      denom <- 1 + change * Qinv[j, j]
      qcol <- Qinv[, j]
      Qinv <- Qinv - (change / denom) * tcrossprod(qcol)
      Qbeta[j, j] <- Qbeta[j, j] + change
      logdetQ <- logdetQ + log(denom)
      tau2[j] <- new_u
      if (m > burn && abs(s_new - s_current) > 1e-10) move_tau[j] <- move_tau[j] + 1L
    }

    # Step 5: update fusion local scales using the same exact rank-one identity.
    for (j in seq_len(pp - 1L)) {
      old_u <- omega2[j]
      s_current <- log(old_u)
      d <- as.numeric(D[j, ])
      qd <- as.numeric(Qinv %*% d)
      dqd <- sum(d * qd)
      logf_omega <- function(s) {
        if (!is.finite(s) || s > 700 || s < -700) return(-Inf)
        new_u <- exp(s)
        change <- 1 / new_u - 1 / old_u
        denom <- 1 + change * dqd
        if (!is.finite(denom) || denom <= 0) return(-Inf)
        ld <- logdetQ + log(denom)
        log_omega2_target(s, j, beta, sigma, lambda2sq, ld, D)
      }
      s_new <- slice_sample_1d(s_current, logf_omega, width = slice_width_omega)
      new_u <- exp(s_new)
      change <- 1 / new_u - 1 / old_u
      denom <- 1 + change * dqd
      Qinv <- Qinv - (change / denom) * tcrossprod(qd)
      Qbeta <- Qbeta + change * tcrossprod(d)
      logdetQ <- logdetQ + log(denom)
      omega2[j] <- new_u
      if (m > burn && abs(s_new - s_current) > 1e-10) move_omega[j] <- move_omega[j] + 1L
    }

    # Step 6: lambda_1^2 | rest ~ Gamma(a1+p, b1+sum(tau2)/2)
    lambda1sq <- rgamma(1, shape = a1 + pp,
                        rate = b1 + 0.5 * sum(tau2))

    # Step 7: lambda_2^2 | rest ~ Gamma(a2+p-1, b2+sum(omega2)/2)
    lambda2sq <- rgamma(1, shape = a2 + pp - 1L,
                        rate = b2 + 0.5 * sum(omega2))

    # Step 8: update sigma | rest using Q_beta from the UPDATED local scales
    Qbeta <- make_Qbeta(tau2, omega2, D)
    vf <- pmax(v, 1e-12)
    delta <- as.numeric(y - beta0 - X %*% beta - phi1 * vf)
    a_star <- a_sigma + (3 * n + pp) / 2
    b_star <- b_sigma +
      sum(delta^2 / (2 * phi2 * vf)) +
      sum(vf) +
      0.5 * as.numeric(crossprod(beta, Qbeta %*% beta))
    sigma <- 1 / rgamma(1, shape = a_star, rate = b_star)

    if (m > burn) {
      k <- m - burn
      S0[k] <- beta0
      S1[k, ] <- beta
      Ssigma[k] <- sigma
      Slam1[k] <- lambda1sq
      Slam2[k] <- lambda2sq
      if (keep_trace) {
        Stau2[k, ] <- tau2
        Somega2[k, ] <- omega2
      }
    }
  }

  ans <- list(
    coef = c(mean(S0), colMeans(S1)),
    intercept_draws = S0,
    draws = S1,
    sigma_draws = Ssigma,
    lambda1sq_draws = Slam1,
    lambda2sq_draws = Slam2,
    tau2_draws = Stau2,
    omega2_draws = Somega2,
    move_rate_tau = move_tau / K,
    move_rate_omega = move_omega / K,
    slice_width_tau = slice_width_tau,
    slice_width_omega = slice_width_omega
  )

  if (keep_diagnostics) {
    ans$diagnostics <- list(
      move_rate_tau = ans$move_rate_tau,
      move_rate_omega = ans$move_rate_omega,
      ess_beta = as.numeric(coda::effectiveSize(coda::mcmc(S1))),
      ess_sigma = as.numeric(coda::effectiveSize(coda::mcmc(Ssigma))),
      ess_lambda1sq = as.numeric(coda::effectiveSize(coda::mcmc(Slam1))),
      ess_lambda2sq = as.numeric(coda::effectiveSize(coda::mcmc(Slam2)))
    )
  }
  ans
}

# 5b. Bayesian Lasso Quantile Regression (BLQR) comparator
# beta_j | tau2_j,sigma ~ N(0,sigma*tau2_j); tau2_j ~ Exp(lambda1sq/2).
# The GIG local-scale conditional includes the normalized Gaussian-prior factor.
fitBLQR <- function(X, y, tau, iter = ITER, burn = BURN,
                    a1 = 1, b1 = 1, a_sigma = 1, b_sigma = 1,
                    keep_diagnostics = TRUE) {
  X <- as.matrix(X); y <- as.numeric(y)
  n <- nrow(X); pp <- ncol(X)
  if (length(y) != n || any(!is.finite(X)) || any(!is.finite(y)))
    stop("BLQR requires finite X and y with matching dimensions.")
  if (tau <= 0 || tau >= 1 || burn < 0 || iter <= burn)
    stop("Invalid tau, iteration, or burn-in setting.")
  phi1 <- (1 - 2 * tau) / (tau * (1 - tau))
  phi2 <- 2 / (tau * (1 - tau))
  Z <- cbind(1, X)
  init <- tryCatch(quantreg::rq.fit.br(Z, y, tau=tau)$coefficients,
                   error=function(e) NULL)
  if (is.null(init) || length(init) != pp+1L || any(!is.finite(init))) {
    beta0 <- as.numeric(quantile(y, tau, names=FALSE)); beta <- rep(0, pp)
  } else { beta0 <- as.numeric(init[1]); beta <- as.numeric(init[-1]) }
  sigma <- max(var(y), 1e-3); v <- rep(sigma, n)
  tau2 <- pmax(abs(beta)*sqrt(sigma)/max(mean(abs(beta)), 0.1), 0.1)
  lambda1sq <- 1; K <- iter-burn
  S0 <- numeric(K); S1 <- matrix(NA_real_, K, pp)
  Ssigma <- numeric(K); Slam1 <- numeric(K)
  for (m in seq_len(iter)) {
    Qbeta <- diag(1/tau2, nrow=pp)
    Qhat <- matrix(0, pp+1L, pp+1L); Qhat[2:(pp+1L),2:(pp+1L)] <- Qbeta
    rV <- 1/pmax(v, 1e-12)
    A <- crossprod(Z, Z*rV)/phi2 + Qhat; A <- (A+t(A))/2
    rhs <- as.numeric(crossprod(Z,(y-phi1*v)*rV))/phi2
    Rc <- strict_chol(A, 'BLQR coefficient precision matrix')
    mu <- backsolve(Rc, forwardsolve(t(Rc), rhs))
    xi <- as.numeric(mu + sqrt(sigma)*backsolve(Rc, rnorm(pp+1L)))
    beta0 <- xi[1]; beta <- xi[-1]
    resid <- as.numeric(y-beta0-X%*%beta)
    chi <- pmax(resid^2/(phi2*sigma), 1e-12)
    psi <- (phi1^2/phi2+2)/sigma
    v <- vapply(seq_len(n), function(i)
      GIGrvg::rgig(1, lambda=0.5, chi=chi[i], psi=psi), numeric(1))
    # tau2_j | rest ~ GIG(1/2, beta_j^2/sigma, lambda1sq).
    for (j in seq_len(pp)) tau2[j] <- GIGrvg::rgig(
      1, lambda=0.5, chi=max(beta[j]^2/sigma, 1e-12), psi=lambda1sq)
    lambda1sq <- rgamma(1, shape=a1+pp, rate=b1+0.5*sum(tau2))
    Qbeta <- diag(1/tau2, nrow=pp); vf <- pmax(v, 1e-12)
    delta <- as.numeric(y-beta0-X%*%beta-phi1*vf)
    a_star <- a_sigma+(3*n+pp)/2
    b_star <- b_sigma+sum(delta^2/(2*phi2*vf))+sum(vf)+
      0.5*as.numeric(crossprod(beta,Qbeta%*%beta))
    sigma <- 1/rgamma(1, shape=a_star, rate=b_star)
    if (m>burn) { k<-m-burn; S0[k]<-beta0; S1[k,]<-beta
      Ssigma[k]<-sigma; Slam1[k]<-lambda1sq }
  }
  ans <- list(coef=c(mean(S0),colMeans(S1)), intercept_draws=S0,
              draws=S1, sigma_draws=Ssigma, lambda1sq_draws=Slam1)
  if (keep_diagnostics) ans$diagnostics <- list(
    ess_beta=as.numeric(coda::effectiveSize(coda::mcmc(S1))),
    ess_sigma=as.numeric(coda::effectiveSize(coda::mcmc(Ssigma))),
    ess_lambda1sq=as.numeric(coda::effectiveSize(coda::mcmc(Slam1))))
  ans
}

# 6. Variable-selection metrics
# Frequentist primary threshold = 0.10, with 0.05 and 0.20 sensitivity rules.
selectionFromCoef <- function(coef, threshold = 0.10) {
  selected <- abs(as.numeric(coef[-1])) > threshold
  TP <- sum(selected & true_active)
  TN <- sum(!selected & true_zero)
  FP <- sum(selected & true_zero)
  FN <- sum(!selected & true_active)
  c(Selection_Accuracy = 100 * (TP + TN) / p,
    True_Zero_Recovery = 100 * TN / sum(true_zero),
    True_Nonzero_Recovery = 100 * TP / sum(true_active),
    False_Positive_Rate = 100 * FP / sum(true_zero),
    False_Negative_Rate = 100 * FN / sum(true_active),
    Selected_Model_Size = sum(selected))
}

# Bayesian selection uses a 95% posterior credible interval excluding zero.
selectionBFLQR <- function(draws) {
  lower <- apply(draws, 2, quantile, probs = 0.025, names = FALSE)
  upper <- apply(draws, 2, quantile, probs = 0.975, names = FALSE)
  selected <- lower > 0 | upper < 0
  TP <- sum(selected & true_active)
  TN <- sum(!selected & true_zero)
  FP <- sum(selected & true_zero)
  FN <- sum(!selected & true_active)
  c(Selection_Accuracy = 100 * (TP + TN) / p,
    True_Zero_Recovery = 100 * TN / sum(true_zero),
    True_Nonzero_Recovery = 100 * TP / sum(true_active),
    False_Positive_Rate = 100 * FP / sum(true_zero),
    False_Negative_Rate = 100 * FN / sum(true_active),
    Selected_Model_Size = sum(selected))
}

# Bayesian selection by 95% posterior credible intervals.
selectionBayesCI <- function(draws) {
  lower <- apply(draws, 2, quantile, probs=0.025, names=FALSE)
  upper <- apply(draws, 2, quantile, probs=0.975, names=FALSE)
  selected <- lower > 0 | upper < 0
  TP <- sum(selected & true_active); TN <- sum(!selected & true_zero)
  FP <- sum(selected & true_zero); FN <- sum(!selected & true_active)
  c(Selection_Accuracy=100*(TP+TN)/p, True_Zero_Recovery=100*TN/sum(true_zero),
    True_Nonzero_Recovery=100*TP/sum(true_active),
    False_Positive_Rate=100*FP/sum(true_zero),
    False_Negative_Rate=100*FN/sum(true_active), Selected_Model_Size=sum(selected))
}

# 7. One Monte Carlo replication: LQR, BLQR, FLQR and BFLQR
oneRun <- function(n, tau, r, error_dist = "ALD", prior_settings = NULL) {
  tr <- gendata(n, tau, error_dist = error_dist)
  te <- gendata(NTEST, tau, error_dist = error_dist) # large common test set
  tro <- tr
  tro$y <- contaminate(tr$y)

  # Fit on clean training data
  cat("  Fitting clean-data LQR (5-fold CV)...\n"); flush.console()
  l <- fitLQR(tr$X, tr$y, tau)
  lqr_lambda <- attr(l, "lambda")
  cat("  Fitting clean-data BLQR...\n"); flush.console()
  if (is.null(prior_settings)) bl <- fitBLQR(tr$X, tr$y, tau)
  else bl <- fitBLQR(tr$X, tr$y, tau,
                     a1=prior_settings$a1, b1=prior_settings$b1,
                     a_sigma=prior_settings$a_sigma, b_sigma=prior_settings$b_sigma)
  cat("  Fitting clean-data FLQR (penalty CV)...\n"); flush.console()
  f <- fitFLQR(tr$X, tr$y, tau)
  flqr_l1 <- attr(f, "lambda1")
  flqr_l2 <- attr(f, "lambda2")
  # Retain local-scale traces only for one representative fit to keep memory manageable.
  keep_trace_this_fit <- (n == min(N) && tau == min(TAU) && r == 1L)
  cat("  Fitting clean-data BFLQR...\n"); flush.console()
  if (is.null(prior_settings)) {
    b <- fitBFLQR(tr$X, tr$y, tau, keep_trace = keep_trace_this_fit)
  } else {
    b <- fitBFLQR(tr$X, tr$y, tau, keep_trace = keep_trace_this_fit,
                  a1=prior_settings$a1, b1=prior_settings$b1,
                  a2=prior_settings$a2, b2=prior_settings$b2,
                  a_sigma=prior_settings$a_sigma, b_sigma=prior_settings$b_sigma)
  }

  # Fit on contaminated training data for the robustness comparison
  cat("  Fitting contaminated-data LQR (5-fold CV)...\n"); flush.console()
  lo <- fitLQR(tro$X, tro$y, tau)
  cat("  Fitting contaminated-data BLQR...\n"); flush.console()
  if (is.null(prior_settings)) blo <- fitBLQR(tro$X, tro$y, tau)
  else blo <- fitBLQR(tro$X, tro$y, tau,
                      a1=prior_settings$a1, b1=prior_settings$b1,
                      a_sigma=prior_settings$a_sigma, b_sigma=prior_settings$b_sigma)
  cat("  Fitting contaminated-data FLQR (penalty CV)...\n"); flush.console()
  fo <- fitFLQR(tro$X, tro$y, tau)
  cat("  Fitting contaminated-data BFLQR...\n"); flush.console()
  if (is.null(prior_settings)) bo <- fitBFLQR(tro$X, tro$y, tau)
  else bo <- fitBFLQR(tro$X, tro$y, tau,
                      a1=prior_settings$a1, b1=prior_settings$b1,
                      a2=prior_settings$a2, b2=prior_settings$b2,
                      a_sigma=prior_settings$a_sigma, b_sigma=prior_settings$b_sigma)

  Zte <- cbind(1, te$X)
  pred_clean <- list(
    LQR=as.numeric(Zte %*% as.numeric(l)), BLQR=as.numeric(Zte %*% bl$coef),
    FLQR=as.numeric(Zte %*% as.numeric(f)), BFLQR=as.numeric(Zte %*% b$coef))
  pred_cont <- list(
    LQR=as.numeric(Zte %*% as.numeric(lo)), BLQR=as.numeric(Zte %*% blo$coef),
    FLQR=as.numeric(Zte %*% as.numeric(fo)), BFLQR=as.numeric(Zte %*% bo$coef))
  cleanQL <- c(LQR=QL(te$y,pred_clean$LQR,tau), BLQR=QL(te$y,pred_clean$BLQR,tau),
               FLQR=QL(te$y,pred_clean$FLQR,tau), BFLQR=QL(te$y,pred_clean$BFLQR,tau))
  outQL <- c(LQR=QL(te$y,pred_cont$LQR,tau), BLQR=QL(te$y,pred_cont$BLQR,tau),
             FLQR=QL(te$y,pred_cont$FLQR,tau), BFLQR=QL(te$y,pred_cont$BFLQR,tau))
  deterioration <- 100 * (outQL - cleanQL) /
    pmax(cleanQL, .Machine$double.eps)

  selections <- rbind(
    LQR=selectionFromCoef(l, threshold=0.10),
    BLQR=selectionBayesCI(bl$draws),
    FLQR=selectionFromCoef(f, threshold=0.10),
    BFLQR=selectionBayesCI(b$draws))
  mse <- c(LQR=mean((as.numeric(l)[-1]-bt)^2),
           BLQR=mean((bl$coef[-1]-bt)^2),
           FLQR=mean((as.numeric(f)[-1]-bt)^2),
           BFLQR=mean((b$coef[-1]-bt)^2))

  # Empirical coverage indicators for 95% posterior intervals. These are
  # calculated for each Bayesian fit and summarized across Monte Carlo reps.
  bl_lower <- apply(bl$draws, 2, quantile, probs=0.025, names=FALSE)
  bl_upper <- apply(bl$draws, 2, quantile, probs=0.975, names=FALSE)
  b_lower <- apply(b$draws, 2, quantile, probs=0.025, names=FALSE)
  b_upper <- apply(b$draws, 2, quantile, probs=0.975, names=FALSE)
  bl_coverage <- mean(bt >= bl_lower & bt <= bl_upper)
  b_coverage <- mean(bt >= b_lower & bt <= b_upper)
  bl_active_coverage <- mean(bt[true_active] >= bl_lower[true_active] &
                             bt[true_active] <= bl_upper[true_active])
  b_active_coverage <- mean(bt[true_active] >= b_lower[true_active] &
                            bt[true_active] <= b_upper[true_active])
  bl_zero_coverage <- mean(bt[true_zero] >= bl_lower[true_zero] &
                           bt[true_zero] <= bl_upper[true_zero])
  b_zero_coverage <- mean(bt[true_zero] >= b_lower[true_zero] &
                          bt[true_zero] <= b_upper[true_zero])

  main <- data.frame(
    n=n, tau=tau, Replication=r,
    Model=c("LQR","BLQR","FLQR","BFLQR"),
    LQR_lambda=c(lqr_lambda,NA_real_,NA_real_,NA_real_),
    BLQR_lambda1sq_mean=c(NA_real_,mean(bl$lambda1sq_draws),NA_real_,NA_real_),
    FLQR_lambda1=c(NA_real_,NA_real_,flqr_l1,NA_real_),
    FLQR_lambda2=c(NA_real_,NA_real_,flqr_l2,NA_real_),
    Coefficient_Recovery_MSE = as.numeric(mse),
    BLQR_CI_Coverage = c(NA_real_, 100*bl_coverage, NA_real_, NA_real_),
    BFLQR_CI_Coverage = c(NA_real_, NA_real_, NA_real_, 100*b_coverage),
    BLQR_Active_CI_Coverage = c(NA_real_, 100*bl_active_coverage, NA_real_, NA_real_),
    BFLQR_Active_CI_Coverage = c(NA_real_, NA_real_, NA_real_, 100*b_active_coverage),
    BLQR_Zero_CI_Coverage = c(NA_real_, 100*bl_zero_coverage, NA_real_, NA_real_),
    BFLQR_Zero_CI_Coverage = c(NA_real_, NA_real_, NA_real_, 100*b_zero_coverage),
    Selection_Accuracy = selections[, "Selection_Accuracy"],
    True_Zero_Recovery = selections[, "True_Zero_Recovery"],
    True_Nonzero_Recovery = selections[, "True_Nonzero_Recovery"],
    False_Positive_Rate = selections[, "False_Positive_Rate"],
    False_Negative_Rate = selections[, "False_Negative_Rate"],
    Selected_Model_Size = selections[, "Selected_Model_Size"],
    Clean_Test_QL = as.numeric(cleanQL),
    Contaminated_Training_Test_QL = as.numeric(outQL),
    Outlier_Deterioration_Percent = as.numeric(deterioration)
  )

  # Selection-threshold sensitivity results.
  # Bayesian median-threshold sensitivity is reported separately from credible-interval selection.
  bmed <- c(NA_real_, apply(b$draws, 2, median))
  sens <- rbind(
    data.frame(n=n, tau=tau, Replication=r, Model="LQR",
               Threshold=0.05, t(selectionFromCoef(l, 0.05))),
    data.frame(n=n, tau=tau, Replication=r, Model="LQR",
               Threshold=0.20, t(selectionFromCoef(l, 0.20))),
    data.frame(n=n, tau=tau, Replication=r, Model="FLQR",
               Threshold=0.05, t(selectionFromCoef(f, 0.05))),
    data.frame(n=n, tau=tau, Replication=r, Model="FLQR",
               Threshold=0.20, t(selectionFromCoef(f, 0.20))),
    data.frame(n=n, tau=tau, Replication=r, Model="BFLQR(median)",
               Threshold=0.05, t(selectionFromCoef(bmed, 0.05))),
    data.frame(n=n, tau=tau, Replication=r, Model="BFLQR(median)",
               Threshold=0.10, t(selectionFromCoef(bmed, 0.10))),
    data.frame(n=n, tau=tau, Replication=r, Model="BFLQR(median)",
               Threshold=0.20, t(selectionFromCoef(bmed, 0.20))),
    data.frame(n=n, tau=tau, Replication=r, Model="BLQR(median)",
               Threshold=0.05, t(selectionFromCoef(c(NA_real_, apply(bl$draws,2,median)), 0.05))),
    data.frame(n=n, tau=tau, Replication=r, Model="BLQR(median)",
               Threshold=0.10, t(selectionFromCoef(c(NA_real_, apply(bl$draws,2,median)), 0.10))),
    data.frame(n=n, tau=tau, Replication=r, Model="BLQR(median)",
               Threshold=0.20, t(selectionFromCoef(c(NA_real_, apply(bl$draws,2,median)), 0.20)))
  )

  required_diag <- list(b$diagnostics, bo$diagnostics, bl$diagnostics, blo$diagnostics)
  if (any(vapply(required_diag, is.null, logical(1)))) {
    stop("Sampler diagnostics were not returned for one or more Bayesian fits.")
  }
  diag <- data.frame(
    n=n, tau=tau, Replication=r,
    Fit=c("Clean_BFLQR", "Contaminated_BFLQR", "Clean_BLQR", "Contaminated_BLQR"),
    Mean_move_rate_tau=c(mean(b$move_rate_tau), mean(bo$move_rate_tau), NA_real_, NA_real_),
    Mean_move_rate_omega=c(mean(b$move_rate_omega), mean(bo$move_rate_omega), NA_real_, NA_real_),
    Min_ESS_beta=c(min(b$diagnostics$ess_beta), min(bo$diagnostics$ess_beta),
                   min(bl$diagnostics$ess_beta), min(blo$diagnostics$ess_beta)),
    ESS_sigma=c(b$diagnostics$ess_sigma, bo$diagnostics$ess_sigma,
                bl$diagnostics$ess_sigma, blo$diagnostics$ess_sigma),
    ESS_lambda1sq=c(b$diagnostics$ess_lambda1sq, bo$diagnostics$ess_lambda1sq,
                    bl$diagnostics$ess_lambda1sq, blo$diagnostics$ess_lambda1sq),
    ESS_lambda2sq=c(b$diagnostics$ess_lambda2sq, bo$diagnostics$ess_lambda2sq,
                    NA_real_, NA_real_)
  )

  list(performance = main, sensitivity = sens, diagnostics = diag,
       representative_fit = if (keep_trace_this_fit) b else NULL)
}

# 8. Optional convergence checks (disabled by default for supplementary-only runs)
# The simulation will not proceed unless all monitored parameters pass diagnostics.
# Numerical rank-one determinant/inverse identity check.
check_rank_one_updates <- function(tol = 1e-8) {
  set.seed(771)
  pp0 <- 6L; D0 <- make_D(pp0)
  t0 <- exp(rnorm(pp0)); o0 <- exp(rnorm(pp0 - 1L))
  Q0 <- make_Qbeta(t0, o0, D0)
  R0 <- strict_chol(Q0, "rank-one self-check Q")
  i0 <- chol2inv(R0); ld0 <- 2 * sum(log(diag(R0)))
  for (kind in c("tau", "omega")) {
    j <- 2L; u <- if (kind == "tau") t0[j] * 1.7 else o0[j] * 0.6
    if (kind == "tau") {
      cc <- 1/u - 1/t0[j]; q <- i0[, j]
      Q1 <- Q0; Q1[j, j] <- Q1[j, j] + cc
      den <- 1 + cc*q[j]
      i1 <- i0 - (cc/den) * tcrossprod(q)
      ld1 <- ld0 + log(den)
    } else {
      d <- D0[j, ]; q <- as.numeric(i0 %*% d); cc <- 1/u - 1/o0[j]
      Q1 <- Q0 + cc * tcrossprod(d)
      den <- 1 + cc*sum(d*q)
      i1 <- i0 - (cc/den) * tcrossprod(q)
      ld1 <- ld0 + log(den)
    }
    i_direct <- solve(Q1)
    ld_direct <- as.numeric(determinant(Q1, logarithm = TRUE)$modulus)
    if (!is.finite(den) || den <= 0 || max(abs(i1-i_direct)) > tol || abs(ld1-ld_direct) > tol)
      stop("Rank-one inverse/determinant self-check failed for ", kind, " update.")
  }
  message("Rank-one inverse/determinant self-check passed.")
  invisible(TRUE)
}
check_rank_one_updates()

# Rank-normalized split-Rhat.
# rank normalization must be performed on the POOLED draws across
# all split chains, not separately within each chain.
rank_normalize_matrix <- function(xmat) {
  vals <- as.numeric(xmat)
  n_all <- length(vals)
  r <- rank(vals, ties.method = "average")
  z <- qnorm((r - 3/8) / (n_all + 1/4))
  matrix(z, nrow = nrow(xmat), ncol = ncol(xmat))
}

classic_rhat <- function(chains) {
  # chains: iterations x chains matrix, already split and transformed
  chains <- as.matrix(chains)
  n <- nrow(chains)
  m <- ncol(chains)
  if (n < 2L || m < 2L) stop("At least 2 draws and 2 chains are required.")
  chain_means <- colMeans(chains)
  W <- mean(apply(chains, 2, var))
  if (!is.finite(W) || W <= 0) return(1)
  B <- n * stats::var(chain_means)
  var_plus <- ((n - 1) / n) * W + B / n
  sqrt(var_plus / W)
}

rank_split_rhat <- function(chain_list) {
  # Split each original chain into first and last halves.
  half <- min(vapply(chain_list,
                     function(z) as.integer(floor(length(z) / 2)),
                     integer(1)))
  if (half < 4L) stop("Chains are too short for split-Rhat.")

  split_chains <- list()
  for (z in chain_list) {
    z <- as.numeric(z)
    split_chains[[length(split_chains) + 1L]] <- z[seq_len(half)]
    split_chains[[length(split_chains) + 1L]] <-
      z[(length(z) - half + 1L):length(z)]
  }
  mat <- do.call(cbind, split_chains)

  # Rank-normalize all draws together, preserving between-chain differences.
  rz <- rank_normalize_matrix(mat)

  # Fold around the pooled median, then again rank-normalize pooled draws.
  pooled_median <- stats::median(as.numeric(mat))
  folded <- abs(mat - pooled_median)
  frz <- rank_normalize_matrix(folded)

  max(classic_rhat(rz), classic_rhat(frz))
}

RUN_CHAIN_CHECK <- FALSE  # Set TRUE to rerun four-chain checks before supplementary experiments
RHAT_LIMIT <- 1.01
ESS_MIN <- 400
if (RUN_CHAIN_CHECK) {
  cat("\n=== BFLQR FOUR-CHAIN CONVERGENCE CHECK ===\n")
  set.seed(83001)
  diagnostic_data <- gendata(min(N), min(TAU))
  chain_mats <- vector("list", 4L)
  chain_acceptance <- vector("list", 4L)
  for (ch in seq_len(4L)) {
    set.seed(83001 + ch)
    fit_ch <- fitBFLQR(diagnostic_data$X, diagnostic_data$y, min(TAU),
                       iter = ITER, burn = BURN, keep_trace = TRUE, init_jitter = TRUE)
    z <- cbind(intercept = fit_ch$intercept_draws,
               fit_ch$draws, sigma = fit_ch$sigma_draws,
               lambda1sq = fit_ch$lambda1sq_draws,
               lambda2sq = fit_ch$lambda2sq_draws,
               fit_ch$tau2_draws, fit_ch$omega2_draws)
    colnames(z) <- c("intercept", paste0("beta", seq_len(p)), "sigma",
                     "lambda1sq", "lambda2sq", paste0("tau2_", seq_len(p)),
                     paste0("omega2_", seq_len(p-1L)))
    chain_mats[[ch]] <- z
    chain_acceptance[[ch]] <- data.frame(
      chain = ch,
      mean_move_rate_tau = mean(fit_ch$move_rate_tau),
      min_move_rate_tau = min(fit_ch$move_rate_tau),
      max_move_rate_tau = max(fit_ch$move_rate_tau),
      mean_move_rate_omega = mean(fit_ch$move_rate_omega),
      min_move_rate_omega = min(fit_ch$move_rate_omega),
      max_move_rate_omega = max(fit_ch$move_rate_omega),
      slice_width_tau = fit_ch$slice_width_tau,
      slice_width_omega = fit_ch$slice_width_omega
    )
    cat("Completed BFLQR chain", ch, "of 4\n"); flush.console()
  }
  param_names <- colnames(chain_mats[[1]])
  rhat_vals <- vapply(param_names, function(nm)
    rank_split_rhat(lapply(chain_mats, function(z) z[, nm])), numeric(1))
  ess_vals <- vapply(param_names, function(nm) {
    mc_list <- do.call(coda::mcmc.list,
                       lapply(chain_mats, function(z) coda::mcmc(z[, nm])))
    as.numeric(coda::effectiveSize(mc_list))
  }, numeric(1))
  convergence_check <- data.frame(parameter = param_names, rank_split_Rhat = rhat_vals,
                          ESS = ess_vals, row.names = NULL)
  cat("\n=== BFLQR R-HAT AND ESS FOR ALL PARAMETERS ===\n")
  print(convergence_check, row.names = FALSE, digits = 6)
  cat("Maximum BFLQR R-hat:", format(max(convergence_check$rank_split_Rhat, na.rm = TRUE), digits = 6),
      " (parameter:", convergence_check$parameter[which.max(convergence_check$rank_split_Rhat)], ")\n")
  cat("Minimum BFLQR ESS:", format(min(convergence_check$ESS, na.rm = TRUE), digits = 6),
      " (parameter:", convergence_check$parameter[which.min(convergence_check$ESS)], ")\n")

  failed <- convergence_check$rank_split_Rhat >= RHAT_LIMIT | convergence_check$ESS < ESS_MIN |
            !is.finite(convergence_check$rank_split_Rhat) | !is.finite(convergence_check$ESS)
  if (any(failed)) {
    failed_table <- convergence_check[failed, , drop = FALSE]
    cat("\nCONVERGENCE CHECK FAILED. Main simulation was NOT started.\n")
    cat("Parameters failing Rhat < ", RHAT_LIMIT,
        " and/or ESS >= ", ESS_MIN, ":\n", sep = "")
    print(failed_table, row.names = FALSE)
    stop("BFLQR convergence gate failed: review the printed convergence diagnostics before changing the sampler.", call. = FALSE)
  }
  cat("Convergence check passed: all monitored BFLQR parameters satisfy R-hat < ", RHAT_LIMIT,
      " and ESS >= ", ESS_MIN, ".\n", sep = "")
}

# BLQR four-chain convergence gate. The BLQR model has no fusion local scales;
# monitor intercept, coefficients, sigma, and its global lasso penalty.
if (RUN_CHAIN_CHECK) {
  cat("\n=== BLQR FOUR-CHAIN CONVERGENCE CHECK ===\n")
  blqr_chain_mats <- vector("list", 4L)
  for (ch in seq_len(4L)) {
    set.seed(84001 + ch)
    fit_blqr_ch <- fitBLQR(diagnostic_data$X, diagnostic_data$y, min(TAU),
                           iter=ITER, burn=BURN, keep_diagnostics=FALSE)
    z <- cbind(intercept=fit_blqr_ch$intercept_draws, fit_blqr_ch$draws,
               sigma=fit_blqr_ch$sigma_draws,
               lambda1sq=fit_blqr_ch$lambda1sq_draws)
    colnames(z) <- c("intercept", paste0("beta", seq_len(p)), "sigma", "lambda1sq")
    blqr_chain_mats[[ch]] <- z
    cat("Completed BLQR chain", ch, "of 4\n"); flush.console()
  }
  blqr_names <- colnames(blqr_chain_mats[[1]])
  blqr_rhat <- vapply(blqr_names, function(nm)
    rank_split_rhat(lapply(blqr_chain_mats, function(z) z[, nm])), numeric(1))
  blqr_ess <- vapply(blqr_names, function(nm) {
    mc <- do.call(coda::mcmc.list,
                  lapply(blqr_chain_mats, function(z) coda::mcmc(z[, nm])))
    as.numeric(coda::effectiveSize(mc))
  }, numeric(1))
  blqr_convergence_check <- data.frame(parameter=blqr_names,
                               rank_split_Rhat=blqr_rhat, ESS=blqr_ess,
                               row.names=NULL)
  cat("\n=== BLQR R-HAT AND ESS FOR ALL PARAMETERS ===\n")
  print(blqr_convergence_check, row.names=FALSE, digits=6)
  cat("Maximum BLQR R-hat:", format(max(blqr_convergence_check$rank_split_Rhat, na.rm=TRUE), digits=6),
      " (parameter:", blqr_convergence_check$parameter[which.max(blqr_convergence_check$rank_split_Rhat)], ")\n")
  cat("Minimum BLQR ESS:", format(min(blqr_convergence_check$ESS, na.rm=TRUE), digits=6),
      " (parameter:", blqr_convergence_check$parameter[which.min(blqr_convergence_check$ESS)], ")\n")
  blqr_failed <- blqr_convergence_check$rank_split_Rhat >= RHAT_LIMIT |
    blqr_convergence_check$ESS < ESS_MIN | !is.finite(blqr_convergence_check$rank_split_Rhat) |
    !is.finite(blqr_convergence_check$ESS)
  if (any(blqr_failed)) {
    stop("BLQR convergence gate failed: review the printed convergence diagnostics before simulation.",
         call.=FALSE)
  }
  cat("BLQR convergence check passed: all monitored parameters satisfy R-hat < ", RHAT_LIMIT,
      " and ESS >= ", ESS_MIN, ".\n", sep="")
}

# Main simulation omitted. The supplementary experiments generate their own data.

# Metrics are displayed as compact tables by topic to avoid very wide console output.
selection_metrics <- c(
  "Coefficient_Recovery_MSE", "Selection_Accuracy", "True_Zero_Recovery",
  "True_Nonzero_Recovery", "False_Positive_Rate", "False_Negative_Rate",
  "Selected_Model_Size"
)
prediction_metrics <- c(
  "Clean_Test_QL", "Contaminated_Training_Test_QL",
  "Outlier_Deterioration_Percent"
)
coverage_metrics <- c(
  "BLQR_CI_Coverage", "BFLQR_CI_Coverage",
  "BLQR_Active_CI_Coverage", "BFLQR_Active_CI_Coverage",
  "BLQR_Zero_CI_Coverage", "BFLQR_Zero_CI_Coverage"
)

summarise_by_group <- function(dat, group_cols, metric_cols) {
  metric_cols <- metric_cols[metric_cols %in% names(dat)]
  if (nrow(dat) == 0L) stop("Cannot create tables from an empty results dataset.")
  if (length(metric_cols) == 0L) return(dat[FALSE, group_cols, drop = FALSE])
  aggregate(
    dat[, metric_cols, drop = FALSE],
    by = dat[, group_cols, drop = FALSE],
    FUN = function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  )
}

print_result_tables <- function(dat, group_cols, label) {
  groups <- c(group_cols)
  tab_select <- summarise_by_group(dat, groups, selection_metrics)
  tab_predict <- summarise_by_group(dat, groups, prediction_metrics)
  tab_cover <- summarise_by_group(dat, groups, coverage_metrics)

  cat("\n", label, " — Coefficient recovery and variable selection\n", sep = "")
  if (ncol(tab_select) > length(groups)) print(tab_select, row.names = FALSE, digits = 4) else cat("No selection metrics available.\n")

  cat("\n", label, " — Predictive performance\n", sep = "")
  if (ncol(tab_predict) > length(groups)) print(tab_predict, row.names = FALSE, digits = 4) else cat("No predictive metrics available.\n")

  cat("\n", label, " — Credible-interval coverage\n", sep = "")
  if (ncol(tab_cover) > length(groups)) print(tab_cover, row.names = FALSE, digits = 4) else cat("No coverage metrics available.\n")
}

# Shared failure log for both supplementary experiments.
failure_rows <- list()
idx_failure <- 1L

# Prior sensitivity: compare baseline and alternative priors on matched datasets.
if (RUN_PRIOR_SENSITIVITY) {
  cat("\n=== PAIRED PRIOR SENSITIVITY EXPERIMENT ===\n")
  prior_rows <- list(); idx_prior <- 1L
  prior_sensitivity_rows <- list(); idx_prior_sens <- 1L
  prior_diagnostic_rows <- list(); idx_prior_diag <- 1L
  for (nn in N) for (tt in TAU) for (rr in seq_len(R_SUPP)) {
    cat(sprintf("Prior sensitivity: n=%d tau=%.2f rep=%d/%d\n", nn, tt, rr, R_SUPP))
    paired_seed <- 920000 + nn*10000 + round(tt*100)*100 + rr

    # Baseline prior: Gamma(1,1) global shrinkage and Gamma(1,1) sigma prior.
    set.seed(paired_seed)
    base <- tryCatch(
      oneRun(nn, tt, rr, error_dist="ALD", prior_settings=NULL),
      error=function(e) {
        message("Baseline-prior sensitivity fit failed: ", conditionMessage(e))
        failure_rows[[idx_failure]] <<- data.frame(Experiment="Prior sensitivity",
          n=nn, tau=tt, Replication=rr, Scenario="Baseline prior",
          Error=conditionMessage(e), stringsAsFactors=FALSE)
        idx_failure <<- idx_failure + 1L
        NULL
      }
    )
    if (!is.null(base)) {
      z <- base$performance
      z$PriorScenario <- "Baseline prior"
      prior_rows[[idx_prior]] <- z; idx_prior <- idx_prior + 1L
      if (!is.null(base$sensitivity)) {
        q <- base$sensitivity; q$PriorScenario <- "Baseline prior"
        prior_sensitivity_rows[[idx_prior_sens]] <- q; idx_prior_sens <- idx_prior_sens + 1L
      }
      if (!is.null(base$diagnostics)) {
        q <- base$diagnostics; q$PriorScenario <- "Baseline prior"
        prior_diagnostic_rows[[idx_prior_diag]] <- q; idx_prior_diag <- idx_prior_diag + 1L
      }
    }

    # Reset seed to reuse the same generated train/test data and contamination.
    set.seed(paired_seed)
    alt <- tryCatch(
      oneRun(nn, tt, rr, error_dist="ALD", prior_settings=PRIOR_ALT),
      error=function(e) {
        message("Alternative-prior sensitivity fit failed: ", conditionMessage(e))
        failure_rows[[idx_failure]] <<- data.frame(Experiment="Prior sensitivity",
          n=nn, tau=tt, Replication=rr, Scenario="Alternative prior",
          Error=conditionMessage(e), stringsAsFactors=FALSE)
        idx_failure <<- idx_failure + 1L
        NULL
      }
    )
    if (!is.null(alt)) {
      z <- alt$performance
      z$PriorScenario <- "Alternative prior"
      prior_rows[[idx_prior]] <- z; idx_prior <- idx_prior + 1L
      if (!is.null(alt$sensitivity)) {
        q <- alt$sensitivity; q$PriorScenario <- "Alternative prior"
        prior_sensitivity_rows[[idx_prior_sens]] <- q; idx_prior_sens <- idx_prior_sens + 1L
      }
      if (!is.null(alt$diagnostics)) {
        q <- alt$diagnostics; q$PriorScenario <- "Alternative prior"
        prior_diagnostic_rows[[idx_prior_diag]] <- q; idx_prior_diag <- idx_prior_diag + 1L
      }
    }
  }
  if (length(prior_rows)) {
    PRIOR_SENS_SIM <- do.call(rbind, prior_rows)
    print_result_tables(
      PRIOR_SENS_SIM,
      c("n", "tau", "Model", "PriorScenario"),
      "PRIOR SENSITIVITY RESULTS (means across replications)"
    )
    save_csv(PRIOR_SENS_SIM, "Prior_Sensitivity_Raw.csv")
    save_csv(summarise_by_group(PRIOR_SENS_SIM, c("n", "tau", "Model", "PriorScenario"), selection_metrics),
             "Prior_Selection.csv")
    save_csv(summarise_by_group(PRIOR_SENS_SIM, c("n", "tau", "Model", "PriorScenario"), prediction_metrics),
             "Prior_Prediction.csv")
    save_csv(summarise_by_group(PRIOR_SENS_SIM, c("n", "tau", "Model", "PriorScenario"), coverage_metrics),
             "Prior_Coverage.csv")
    if (length(prior_sensitivity_rows)) {
      PRIOR_THRESHOLD_SIM <- do.call(rbind, prior_sensitivity_rows)
      save_csv(PRIOR_THRESHOLD_SIM, "Prior_Threshold_Sensitivity.csv")
    }
    if (length(prior_diagnostic_rows)) {
      PRIOR_DIAGNOSTICS <- do.call(rbind, prior_diagnostic_rows)
      save_csv(PRIOR_DIAGNOSTICS, "Prior_Sampler_Diagnostics.csv")
    }
  } else {
    warning("No prior-sensitivity replications completed; no result tables were produced.")
  }
}

# Robustness experiment with shifted Student-t errors (3 degrees of freedom).
if (RUN_NON_ALD_EXPERIMENT) {
  cat("\n=== NON-ALD ERROR EXPERIMENT (SHIFTED t WITH 3 DF) ===\n")
  nonald_rows <- list(); idx_nonald <- 1L
  nonald_sensitivity_rows <- list(); idx_nonald_sens <- 1L
  nonald_diagnostic_rows <- list(); idx_nonald_diag <- 1L
  for (nn in N) for (tt in TAU) for (rr in seq_len(R_SUPP)) {
    cat(sprintf("Non-ALD: n=%d tau=%.2f rep=%d/%d\n", nn, tt, rr, R_SUPP))
    set.seed(930000 + nn*10000 + round(tt*100)*100 + rr)
    res <- tryCatch(oneRun(nn, tt, rr, error_dist="t3_shifted"),
                    error=function(e) {
                      message("Non-ALD replication failed: ", conditionMessage(e))
                      failure_rows[[idx_failure]] <<- data.frame(Experiment="Non-ALD robustness",
                        n=nn, tau=tt, Replication=rr, Scenario="Shifted t(3)",
                        Error=conditionMessage(e), stringsAsFactors=FALSE)
                      idx_failure <<- idx_failure + 1L
                      NULL
                    })
    if (!is.null(res)) {
      z <- res$performance
      z$ErrorDistribution <- "Shifted t(3)"
      nonald_rows[[idx_nonald]] <- z; idx_nonald <- idx_nonald+1L
      if (!is.null(res$sensitivity)) {
        q <- res$sensitivity; q$ErrorDistribution <- "Shifted t(3)"
        nonald_sensitivity_rows[[idx_nonald_sens]] <- q; idx_nonald_sens <- idx_nonald_sens + 1L
      }
      if (!is.null(res$diagnostics)) {
        q <- res$diagnostics; q$ErrorDistribution <- "Shifted t(3)"
        nonald_diagnostic_rows[[idx_nonald_diag]] <- q; idx_nonald_diag <- idx_nonald_diag + 1L
      }
    }
  }
  if (length(nonald_rows)) {
    NONALD_SIM <- do.call(rbind, nonald_rows)
    print_result_tables(
      NONALD_SIM,
      c("n", "tau", "Model", "ErrorDistribution"),
      "NON-ALD ERROR RESULTS (means across replications)"
    )
    save_csv(NONALD_SIM, "Non_ALD_Raw.csv")
    save_csv(summarise_by_group(NONALD_SIM, c("n", "tau", "Model", "ErrorDistribution"), selection_metrics),
             "Non_ALD_Selection.csv")
    save_csv(summarise_by_group(NONALD_SIM, c("n", "tau", "Model", "ErrorDistribution"), prediction_metrics),
             "Non_ALD_Prediction.csv")
    save_csv(summarise_by_group(NONALD_SIM, c("n", "tau", "Model", "ErrorDistribution"), coverage_metrics),
             "Non_ALD_Coverage.csv")
    if (length(nonald_sensitivity_rows)) {
      NONALD_THRESHOLD_SIM <- do.call(rbind, nonald_sensitivity_rows)
      save_csv(NONALD_THRESHOLD_SIM, "Non_ALD_Threshold_Sensitivity.csv")
    }
    if (length(nonald_diagnostic_rows)) {
      NONALD_DIAGNOSTICS <- do.call(rbind, nonald_diagnostic_rows)
      save_csv(NONALD_DIAGNOSTICS, "Non_ALD_Sampler_Diagnostics.csv")
    }
  }
}

# Save a failure log so that incomplete Monte Carlo cells are visible and auditable.
if (exists("failure_rows") && length(failure_rows)) {
  SUPPLEMENTARY_FAILURES <- do.call(rbind, failure_rows)
  save_csv(SUPPLEMENTARY_FAILURES, "Supplementary_Failures.csv")
} else {
  SUPPLEMENTARY_FAILURES <- data.frame(Experiment=character(), n=integer(), tau=numeric(),
    Replication=integer(), Scenario=character(), Error=character())
  save_csv(SUPPLEMENTARY_FAILURES, "Supplementary_Failures.csv")
}

# Record software and package versions for reproducibility.
writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "sessionInfo.txt"))
writeLines(c(
  "BFLQR supplementary simulation analyses",
  paste0("N = ", paste(N, collapse = ", ")),
  paste0("TAU = ", paste(TAU, collapse = ", ")),
  paste0("Replications per cell = ", R_SUPP),
  paste0("MCMC iterations = ", ITER),
  paste0("Burn-in = ", BURN),
  paste0("Test sample size = ", NTEST),
  paste0("Predictors = ", p),
  paste0("AR(1) correlation = ", rho),
  paste0("Prior sensitivity enabled = ", RUN_PRIOR_SENSITIVITY),
  paste0("Non-ALD experiment enabled = ", RUN_NON_ALD_EXPERIMENT),
  paste0("Completed at = ", Sys.time())
), file.path(OUTPUT_DIR, "Supplementary_Run_Info.txt"))
cat("\nSupplementary analyses finished. Results saved in:\n",
    normalizePath(OUTPUT_DIR, winslash = "/", mustWork = TRUE), "\n", sep = "")
