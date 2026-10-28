# Real-data analysis: LQR, BLQR, FLQR, and BFLQR
# Predictors X1-X10 follow the order and definitions given in the manuscript.

suppressPackageStartupMessages({
  library(quantreg)
  library(MASS)
  library(GIGrvg)
  library(coda)
})

set.seed(2026)
TAU <- c(0.25, 0.50, 0.75)
ITER <- 24000
BURN <- 12000
MCMC_MAX_ITER <- 36000
MCMC_MAX_BURN <- 18000
MCMC_RHAT_LIMIT <- 1.01
MCMC_ESS_MIN <- 400
BOOT_B <- 500

# Read the supplied, already-standardized annual dataset.
DATA_FILE <- "thesisdata.csv"
D0 <- read.csv(DATA_FILE, stringsAsFactors = FALSE)
req <- c("Y", paste0("X", 1:10))
D <- D0[, req, drop = FALSE]
y <- D$Y
Xraw <- as.matrix(D[, paste0("X", 1:10), drop = FALSE])
storage.mode(Xraw) <- "double"
# Trend and lag-1 serial-dependence summaries for the annual series.
trend_diag <- do.call(rbind,lapply(names(D),function(v){z<-D[[v]]
ct<-suppressWarnings(cor.test(seq_along(z),z,method='spearman',exact=FALSE))
data.frame(Variable=v,Spearman_Time_Correlation=unname(ct$estimate),
    Trend_P_Value=ct$p.value,Lag1_ACF=as.numeric(acf(z,plot=FALSE,lag.max=1)$acf[2L]))}))
cat('\nTREND AND SERIAL-DEPENDENCE DIAGNOSTICS\n')
print(trend_diag, row.names=FALSE)
cat('Annual observations:',nrow(D),'\nTraining observations:',nrow(D)-max(1L,ceiling(.20*nrow(Xraw))),'\nChronological holdout observations:',max(1L,ceiling(.20*nrow(Xraw))),'\n')
VarNames <- c('GDP per capita','GDP growth','Industry','Foreign Direct Investment',
    'Natural Resources Rent','Energy Use per Capita','Electric Power Consumption',
    'Access to Electricity','Renewable Energy Consumption','Urban Population Growth')
PN <- c('Intercept',VarNames)
# The supplied Y and X1-X10 are already standardized. Use them as supplied;
# do not standardize a second time or back-transform only the predictors.
nt <- max(1L,ceiling(.20*nrow(Xraw)))
tr <- seq_len(nrow(Xraw)-nt)
te <- (max(tr)+1L):nrow(Xraw)
Xt <- Xraw[tr,,drop=FALSE]
yt <- y[tr]
Xv <- Xraw[te,,drop=FALSE]
yv <- y[te]
QL <- function(y,yp,tau){u<-y-as.numeric(yp)
mean(ifelse(u>=0,tau*u,(tau-1)*u))}
MAE <- function(y,yp) mean(abs(y-as.numeric(yp)))
pred <- function(X,b) as.numeric(cbind(1,X)%*%b)

# Expanding-window, one-step-ahead validation splits
splits <- function(n){m<-max(10L,floor(.60*n))
if(m>=n)stop('Insufficient data for rolling-origin tuning.')
lapply(seq.int(m,n-1L),function(t)list(a=seq_len(t),v=t+1L))}

# LQR: tune the lasso penalty by expanding-window, one-step-ahead validation.
# The intercept is unpenalized. Penalties are defined on the summed check loss.
lqr_fixed <- function(X, y, tau, lambda) {
  Z <- cbind(1, X)
  penalty <- c(0, rep(lambda, ncol(X)))
  fit <- quantreg::rq.fit.lasso(Z, y, tau = tau, lambda = penalty)
  as.numeric(fit$coefficients)
}

fitLQR <- function(X, y, tau) {
  validation_splits <- splits(length(y))
  lambda_grid <- length(y) * exp(seq(log(1e-3), log(0.5), length.out = 30))

  cv_loss <- vapply(lambda_grid, function(lambda) {
    losses <- vapply(validation_splits, function(s) {
      coefficient <- tryCatch(
        lqr_fixed(X[s$a, , drop = FALSE], y[s$a], tau,
                  lambda * length(s$a) / length(y)),
        error = function(e) rep(NA_real_, ncol(X) + 1L)
      )
      if (any(!is.finite(coefficient))) return(NA_real_)
      QL(y[s$v], pred(X[s$v, , drop = FALSE], coefficient), tau)
    }, numeric(1))
    if (all(!is.finite(losses))) Inf else mean(losses, na.rm = TRUE)
  }, numeric(1))

  if (!any(is.finite(cv_loss))) stop("All LQR cross-validation fits failed.")
  best_lambda <- lambda_grid[which.min(cv_loss)]
  coefficient <- lqr_fixed(X, y, tau, best_lambda)
  list(coef = coefficient, tuning = best_lambda, cv_loss = min(cv_loss))
}

# FLQR: fused-lasso quantile regression, estimated through linear programming.
# Paired pseudo-rows encode l1*sum(abs(beta)) + l2*sum(abs(D beta)) exactly
# under the summed check-loss objective. The intercept is unpenalized.
flqr_fixed <- function(X, y, tau, lambda1, lambda2) {
  p <- ncol(X)
  D <- make_D(p)
  penalty_matrix <- rbind(
    lambda1 * diag(p),
    -lambda1 * diag(p),
    lambda2 * D,
    -lambda2 * D
  )
  augmented_X <- rbind(
    cbind(1, X),
    cbind(0, penalty_matrix)
  )
  augmented_y <- c(y, rep(0, nrow(penalty_matrix)))
  fit <- quantreg::rq.fit.br(augmented_X, augmented_y, tau = tau)
  as.numeric(fit$coefficients)
}

fitFLQR <- function(X, y, tau) {
  n <- nrow(X)
  validation_splits <- splits(n)
  penalty_grid <- expand.grid(
    lambda1 = n * c(0.002, 0.01, 0.03, 0.1, 0.3),
    lambda2 = n * c(0.002, 0.01, 0.03, 0.1, 0.3)
  )

  cv_loss <- rep(Inf, nrow(penalty_grid))
  for (g in seq_len(nrow(penalty_grid))) {
    fold_loss <- vapply(validation_splits, function(s) {
      scale_factor <- length(s$a) / n
      coefficient <- tryCatch(
        flqr_fixed(
          X[s$a, , drop = FALSE], y[s$a], tau,
          penalty_grid$lambda1[g] * scale_factor,
          penalty_grid$lambda2[g] * scale_factor
        ),
        error = function(e) rep(NA_real_, ncol(X) + 1L)
      )
      if (any(!is.finite(coefficient))) return(NA_real_)
      QL(y[s$v], pred(X[s$v, , drop = FALSE], coefficient), tau)
    }, numeric(1))
    if (any(is.finite(fold_loss))) cv_loss[g] <- mean(fold_loss, na.rm = TRUE)
  }

  if (!any(is.finite(cv_loss))) stop("All FLQR cross-validation fits failed.")
  best <- which.min(cv_loss)
  lambda1 <- penalty_grid$lambda1[best]
  lambda2 <- penalty_grid$lambda2[best]
  coefficient <- flqr_fixed(X, y, tau, lambda1, lambda2)
  list(
    coef = coefficient,
    tuning = c(lambda1 = lambda1, lambda2 = lambda2),
    cv_loss = cv_loss[best]
  )
}

# BFLQR sampler matched to the latest simulation code:
# ALD latent-variable update -> coefficient update -> determinant-adjusted
# log-scale slice updates for tau2 and omega2 -> Gamma global-scale updates
# -> inverse-Gamma sigma update. The rank-one determinant/inverse identities
# are the same as in the simulation sampler.
make_D <- function(pp) {
  D <- matrix(0, nrow = pp - 1L, ncol = pp)
  for (j in seq_len(pp - 1L)) { D[j,j] <- -1
D[j,j+1L] <- 1 }
  D
}
make_Qbeta <- function(tau2, omega2, D) {
  diag(1/tau2, nrow=length(tau2)) +
    crossprod(D, sweep(D, 1L, omega2, FUN=function(x,w) x/w))
}
strict_chol <- function(A, label='matrix') {
  A <- as.matrix(A)
A <- (A+t(A))/2
  if (any(!is.finite(A))) stop(label, ' contains non-finite values.')
  tryCatch(chol(A), error=function(e) stop('Cholesky factorization failed for ', label,
    '. Original error: ', conditionMessage(e), call.=FALSE))
}
log_tau2_target <- function(s,j,beta,sigma,lambda1sq,ld) {
  if (!is.finite(s) || s < -700 || s > 700) return(-Inf)
  u <- exp(s)
  0.5*ld - beta[j]^2/(2*sigma*u) - 0.5*lambda1sq*u + s
}
log_omega2_target <- function(s,j,beta,sigma,lambda2sq,ld,D) {
  if (!is.finite(s) || s < -700 || s > 700) return(-Inf)
  u <- exp(s)
delta_j <- as.numeric(D[j,,drop=FALSE] %*% beta)
  0.5*ld - delta_j^2/(2*sigma*u) - 0.5*lambda2sq*u + s
}
slice_sample_1d <- function(x0, logf, width=1, max_steps=100L, max_shrink=300L) {
  f0 <- logf(x0)
  if (!is.finite(f0)) stop('Non-finite slice density at current local scale.')
  logy <- f0 + log(runif(1))
L <- x0-runif(1,0,width)
R <- L+width
  j <- sample.int(max_steps+1L,1L)-1L
k <- max_steps-j
  while (j>0L && logf(L)>logy) { L <- L-width
j <- j-1L }
  while (k>0L && logf(R)>logy) { R <- R+width
k <- k-1L }
  for (ii in seq_len(max_shrink)) {
    x1 <- runif(1,L,R)
f1 <- logf(x1)
    if (is.finite(f1) && f1>=logy) return(x1)
    if (x1<x0) L <- x1 else R <- x1
    if ((R-L)<1e-12) break
  }
  stop('Slice sampler failed to find an acceptable point; inspect local-scale target.')
}
# BLQR sampler: Bayesian lasso quantile regression without the fusion penalty.
# It uses the same ALD normal-exponential mixture as BFLQR, but the coefficient
# prior precision is diagonal; therefore no adjacent-difference scales enter.
fitBLQR <- function(X,y,tau,iter=ITER,burn=BURN,
                    a1=1,b1=1,a_sigma=1,b_sigma=1) {
  X <- as.matrix(X)
y <- as.numeric(y)
n <- nrow(X)
pp <- ncol(X)
  if (length(y)!=n) stop('nrow(X) must equal length(y).')
  if (tau<=0 || tau>=1) stop('tau must lie strictly between 0 and 1.')
  if (iter<2 || burn<0 || burn>=iter) stop('Require 0 <= burn < iter.')
  if (any(!is.finite(X)) || any(!is.finite(y))) stop('X and y must be finite.')
  phi1 <- (1-2*tau)/(tau*(1-tau))
phi2 <- 2/(tau*(1-tau))
  Z <- cbind(1,X)
  rq_init <- tryCatch(quantreg::rq.fit.br(Z,y,tau=tau)$coefficients,error=function(e) NULL)
  if (is.null(rq_init) || length(rq_init)!=pp+1L || any(!is.finite(rq_init))) {
    beta0 <- as.numeric(quantile(y,probs=tau,names=FALSE))
beta <- rep(0,pp)
  } else { beta0 <- as.numeric(rq_init[1L])
beta <- as.numeric(rq_init[-1L]) }
  # Disperse chain starting values so R-hat is not based on identical starts.
  y_scale <- max(sd(y), 1e-3)
  beta0 <- beta0 + rnorm(1L, 0, 0.10*y_scale)
  beta <- beta + rnorm(pp, 0, 0.10*y_scale/sqrt(pp))
  sigma <- max(var(y),1e-3)
v <- rep(max(sigma,1e-3),n)
  tau2 <- pmax(abs(beta)*sqrt(sigma)/max(mean(abs(beta)),0.1),0.1)
  lambda1sq <- exp(rnorm(1L, mean=0, sd=0.35))
K <- iter-burn
  S0 <- numeric(K)
S1 <- matrix(NA_real_,K,pp)
Ssigma <- numeric(K)
  Slam1 <- numeric(K)
  for (m in seq_len(iter)) {
    Qbeta <- diag(1/tau2, nrow=pp)
    Qhat <- matrix(0,pp+1L,pp+1L)
Qhat[2:(pp+1L),2:(pp+1L)] <- Qbeta
    rV <- 1/pmax(v,1e-12)
    A <- crossprod(Z,Z*rV)/phi2 + Qhat
A <- (A+t(A))/2
    rhs <- as.numeric(crossprod(Z,(y-phi1*v)*rV))/phi2
    Rchol <- strict_chol(A,'BLQR coefficient precision matrix')
    mu_xi <- backsolve(Rchol,forwardsolve(t(Rchol),rhs))
    xi <- as.numeric(mu_xi + sqrt(sigma)*backsolve(Rchol,rnorm(pp+1L)))
    beta0 <- xi[1L]
beta <- xi[-1L]
    resid <- as.numeric(y-beta0-X%*%beta)
    chi <- pmax(resid^2/(phi2*sigma),1e-12)
    psi <- (phi1^2/phi2+2)/sigma
    v <- vapply(seq_len(n),function(i) GIGrvg::rgig(1,lambda=0.5,chi=chi[i],psi=psi),numeric(1))
    # Update local coefficient scales on the log scale, including the
    # Gaussian-prior determinant and transformation Jacobian.
    Qbeta <- diag(1/tau2, nrow=pp)
    Rq <- strict_chol(Qbeta,'BLQR local-scale precision matrix')
    Qinv <- chol2inv(Rq)
logdetQ <- 2*sum(log(diag(Rq)))
    for (j in seq_len(pp)) {
      old_u <- tau2[j]
s_current <- log(old_u)
      logf <- function(s) {
        if (!is.finite(s)||s>700||s< -700) return(-Inf)
        new_u <- exp(s)
change <- 1/new_u-1/old_u
        denom <- 1+change*Qinv[j,j]
        if (!is.finite(denom)||denom<=0) return(-Inf)
        0.5*(logdetQ+log(denom))-beta[j]^2/(2*sigma*new_u)-
          0.5*lambda1sq*new_u+s
      }
      s_new <- slice_sample_1d(s_current,logf,width=1)
new_u <- exp(s_new)
      change <- 1/new_u-1/old_u
denom <- 1+change*Qinv[j,j]
      qcol <- Qinv[,j]
Qinv <- Qinv-(change/denom)*tcrossprod(qcol)
      Qbeta[j,j] <- Qbeta[j,j]+change
logdetQ <- logdetQ+log(denom)
tau2[j] <- new_u
    }
    lambda1sq <- rgamma(1,shape=a1+pp,rate=b1+0.5*sum(tau2))
    Qbeta <- diag(1/tau2,nrow=pp)
vf <- pmax(v,1e-12)
    delta <- as.numeric(y-beta0-X%*%beta-phi1*vf)
    a_star <- a_sigma+(3*n+pp)/2
    b_star <- b_sigma+sum(delta^2/(2*phi2*vf))+sum(vf)+
      0.5*as.numeric(crossprod(beta,Qbeta%*%beta))
    sigma <- 1/rgamma(1,shape=a_star,rate=b_star)
    if (m %% 1000L==0L) message('BLQR iteration ',m,'/',iter,' at tau=',tau)
    if (m>burn) { k <- m-burn
S0[k] <- beta0
S1[k,] <- beta
Ssigma[k] <- sigma
Slam1[k] <- lambda1sq }
  }
  draws <- cbind(Intercept=S0,S1)
colnames(draws) <- c('Intercept',paste0('X',seq_len(pp)))
  list(coef=colMeans(draws),lower=apply(draws,2,quantile,probs=.025),
       upper=apply(draws,2,quantile,probs=.975),draws=draws,
       sigma_draws=Ssigma,lambda1sq_draws=Slam1)
}

fitBFLQR <- function(X,y,tau,iter=ITER,burn=BURN,
                     a1=1,b1=1,a2=1,b2=1,a_sigma=1,b_sigma=1,
                     keep_trace=FALSE) {
  X <- as.matrix(X)
y <- as.numeric(y)
n <- nrow(X)
pp <- ncol(X)
  if (length(y)!=n) stop('nrow(X) must equal length(y).')
  if (pp<2L) stop('BFLQR requires at least two predictors.')
  if (tau<=0 || tau>=1) stop('tau must lie strictly between 0 and 1.')
  if (iter<2 || burn<0 || burn>=iter) stop('Require 0 <= burn < iter.')
  if (any(!is.finite(X)) || any(!is.finite(y))) stop('X and y must be finite.')
  phi1 <- (1-2*tau)/(tau*(1-tau))
phi2 <- 2/(tau*(1-tau))
  Z <- cbind(1,X)
D <- make_D(pp)
  rq_init <- tryCatch(quantreg::rq.fit.br(Z,y,tau=tau)$coefficients,error=function(e) NULL)
  if (is.null(rq_init) || length(rq_init)!=pp+1L || any(!is.finite(rq_init))) {
    beta0 <- as.numeric(quantile(y,probs=tau,names=FALSE))
beta <- rep(0,pp)
  } else { beta0 <- as.numeric(rq_init[1L])
beta <- as.numeric(rq_init[-1L]) }
  # Disperse chain starting values so R-hat is not based on identical starts.
  y_scale <- max(sd(y), 1e-3)
  beta0 <- beta0 + rnorm(1L, 0, 0.10*y_scale)
  beta <- beta + rnorm(pp, 0, 0.10*y_scale/sqrt(pp))
  sigma <- max(var(y),1e-3)
v <- rep(max(sigma,1e-3),n)
  tau2 <- pmax(abs(beta)*sqrt(sigma)/max(mean(abs(beta)),0.1),0.1)
  omega2 <- pmax(abs(diff(beta))*sqrt(sigma)/max(mean(abs(diff(beta))),0.1),0.1)
  lambda1sq <- exp(rnorm(1L, mean=0, sd=0.35))
lambda2sq <- exp(rnorm(1L, mean=0, sd=0.35))
K <- iter-burn
  S0 <- numeric(K)
S1 <- matrix(NA_real_,K,pp)
Ssigma <- numeric(K)
  Slam1 <- numeric(K)
Slam2 <- numeric(K)
  Stau2 <- if (keep_trace) matrix(NA_real_,K,pp) else NULL
  Somega2 <- if (keep_trace) matrix(NA_real_,K,pp-1L) else NULL
  for (m in seq_len(iter)) {
    # Step 1: coefficient update using current local scales and latent v.
    Qbeta <- make_Qbeta(tau2,omega2,D)
    Qhat <- matrix(0,pp+1L,pp+1L)
Qhat[2:(pp+1L),2:(pp+1L)] <- Qbeta
    rV <- 1/pmax(v,1e-12)
    A <- crossprod(Z,Z*rV)/phi2 + Qhat
A <- (A+t(A))/2
    rhs <- as.numeric(crossprod(Z,(y-phi1*v)*rV))/phi2
    Rchol <- strict_chol(A,'coefficient precision matrix A')
    mu_xi <- backsolve(Rchol,forwardsolve(t(Rchol),rhs))
    xi <- as.numeric(mu_xi + sqrt(sigma)*backsolve(Rchol,rnorm(pp+1L)))
    beta0 <- xi[1L]
beta <- xi[-1L]
    # Step 2: ALD latent variables, exactly as in the simulation code.
    resid <- as.numeric(y-beta0-X%*%beta)
    chi <- pmax(resid^2/(phi2*sigma),1e-12)
    psi <- (phi1^2/phi2+2)/sigma
    v <- vapply(seq_len(n),function(i) GIGrvg::rgig(1,lambda=0.5,chi=chi[i],psi=psi),numeric(1))
    # Cache inverse and log determinant, refreshed each iteration.
    Rq <- strict_chol(Qbeta,'local-scale precision matrix Q_beta')
    Qinv <- chol2inv(Rq)
logdetQ <- 2*sum(log(diag(Rq)))
    # Step 3: coefficient local scales tau2 via exact log-scale slice updates.
    for (j in seq_len(pp)) {
      old_u <- tau2[j]
s_current <- log(old_u)
      logf <- function(s) {
        if (!is.finite(s)||s>700||s< -700) return(-Inf)
        new_u <- exp(s)
change <- 1/new_u-1/old_u
        denom <- 1+change*Qinv[j,j]
        if (!is.finite(denom)||denom<=0) return(-Inf)
        log_tau2_target(s,j,beta,sigma,lambda1sq,logdetQ+log(denom))
      }
      s_new <- slice_sample_1d(s_current,logf,width=1)
new_u <- exp(s_new)
      change <- 1/new_u-1/old_u
denom <- 1+change*Qinv[j,j]
      qcol <- Qinv[,j]
Qinv <- Qinv-(change/denom)*tcrossprod(qcol)
      Qbeta[j,j] <- Qbeta[j,j]+change
logdetQ <- logdetQ+log(denom)
tau2[j] <- new_u
    }
    # Step 4: fusion local scales omega2 via the same rank-one identities.
    for (j in seq_len(pp-1L)) {
      old_u <- omega2[j]
s_current <- log(old_u)
d <- as.numeric(D[j,])
      qd <- as.numeric(Qinv%*%d)
dqd <- sum(d*qd)
      logf <- function(s) {
        if (!is.finite(s)||s>700||s< -700) return(-Inf)
        new_u <- exp(s)
change <- 1/new_u-1/old_u
denom <- 1+change*dqd
        if (!is.finite(denom)||denom<=0) return(-Inf)
        log_omega2_target(s,j,beta,sigma,lambda2sq,logdetQ+log(denom),D)
      }
      s_new <- slice_sample_1d(s_current,logf,width=1)
new_u <- exp(s_new)
      change <- 1/new_u-1/old_u
denom <- 1+change*dqd
      Qinv <- Qinv-(change/denom)*tcrossprod(qd)
      Qbeta <- Qbeta+change*tcrossprod(d)
logdetQ <- logdetQ+log(denom)
omega2[j] <- new_u
    }
    # Step 5: global shrinkage parameters.
    lambda1sq <- rgamma(1,shape=a1+pp,rate=b1+0.5*sum(tau2))
    lambda2sq <- rgamma(1,shape=a2+pp-1L,rate=b2+0.5*sum(omega2))
    # Step 6: sigma update using Q_beta built from UPDATED local scales.
    Qbeta <- make_Qbeta(tau2,omega2,D)
vf <- pmax(v,1e-12)
    delta <- as.numeric(y-beta0-X%*%beta-phi1*vf)
    a_star <- a_sigma+(3*n+pp)/2
    b_star <- b_sigma+sum(delta^2/(2*phi2*vf))+sum(vf)+
      0.5*as.numeric(crossprod(beta,Qbeta%*%beta))
    sigma <- 1/rgamma(1,shape=a_star,rate=b_star)
    if (m %% 1000L==0L) message('BFLQR iteration ',m,'/',iter,' at tau=',tau)
    if (m>burn) {
      k <- m-burn
S0[k] <- beta0
S1[k,] <- beta
Ssigma[k] <- sigma
      Slam1[k] <- lambda1sq
Slam2[k] <- lambda2sq
      if (keep_trace) { Stau2[k,] <- tau2
Somega2[k,] <- omega2 }
    }
  }
  draws <- cbind(Intercept=S0,S1)
colnames(draws) <- c('Intercept',paste0('X',seq_len(pp)))
  list(coef=colMeans(draws),lower=apply(draws,2,quantile,probs=.025),
       upper=apply(draws,2,quantile,probs=.975),draws=draws,
       sigma_draws=Ssigma,lambda1sq_draws=Slam1,lambda2sq_draws=Slam2,
       tau2_draws=Stau2,omega2_draws=Somega2)
}

# Moving-block bootstrap for LQR/FLQR intervals; paired rows remain together.
blockIds <- function(n,L){st<-sample.int(n,ceiling(n/L),replace=TRUE)
id<-unlist(lapply(st,function(s)((s-1L+0:(L-1L))%%n)+1L))
id[seq_len(n)]}
bootCI <- function(X,y,tau,model,tuning,B=BOOT_B){p<-ncol(X)
A<-matrix(NA_real_,B,p+1)
L<-max(2L,floor(length(y)^(1/3)))
for(i in seq_len(B)){id<-blockIds(length(y),L)
A[i,]<-tryCatch(if (model == 'LQR') lqr_fixed(X[id, , drop = FALSE],
     y[id], tau, tuning) else flqr_fixed(X[id, , drop = FALSE], y[id],
     tau, tuning['lambda1'], tuning['lambda2']),error=function(e)rep(NA_real_,
    p+1))}
A<-A[complete.cases(A),,drop=FALSE]
if(nrow(A)<max(20,.7*B))warning(model,': many bootstrap fits failed')
if(nrow(A)<20)stop(model,': too few successful bootstrap fits')
list(lower=apply(A,2,quantile,.025),upper=apply(A,2,quantile,.975),draws=A,successful=nrow(A),block_length=L)}


# Run independent chains and extend them only when convergence criteria are not met.
run_bayesian_chains <- function(fit_function, X, y, tau, model, nchains=4L) {
  if (nchains < 4L) warning('At least four chains are recommended for robust convergence assessment.')
  stage_settings <- list(
    c(iter=ITER, burn=BURN),
    c(iter=MCMC_MAX_ITER, burn=MCMC_MAX_BURN)
  )
  final <- NULL
  for (stage in seq_along(stage_settings)) {
    settings <- stage_settings[[stage]]
    chains <- vector('list', nchains)
    message(model, ': running ', nchains, ' independent chains at tau=', tau,
            ' (', settings[['iter']], ' iterations; ', settings[['burn']], ' burn-in).')
    for (ch in seq_len(nchains)) {
      # Distinct reproducible seeds for model, quantile, stage, and chain.
      model_offset <- if (model == 'BFLQR') 70000L else 30000L
      set.seed(20261021L + model_offset + 10000L*stage + 1000L*ch + round(100*tau))
      chains[[ch]] <- fit_function(X, y, tau,
                                   iter=as.integer(settings[['iter']]),
                                   burn=as.integer(settings[['burn']]))
    }
    chain_mats <- lapply(chains, function(z) {
      out <- cbind(z$draws, sigma=z$sigma_draws,
                   lambda1sq=z$lambda1sq_draws)
      if (model == 'BFLQR') out <- cbind(out, lambda2sq=z$lambda2sq_draws)
      out
    })
    mlist <- do.call(coda::mcmc.list, lapply(chain_mats, coda::mcmc))
    rh <- tryCatch(coda::gelman.diag(mlist, autoburnin=FALSE,
                                     multivariate=FALSE)$psrf[, 1L],
                   error=function(e) rep(NA_real_, ncol(chain_mats[[1L]])))
    names(rh) <- colnames(chain_mats[[1L]])
    ess <- tryCatch(coda::effectiveSize(mlist),
                    error=function(e) rep(NA_real_, ncol(chain_mats[[1L]])))
    names(ess) <- colnames(chain_mats[[1L]])
    diagnostics <- data.frame(
      Model=model, Tau=tau, Parameter=colnames(chain_mats[[1L]]),
      Rhat=as.numeric(rh[colnames(chain_mats[[1L]])]),
      ESS=as.numeric(ess[colnames(chain_mats[[1L]])]),
      Rhat_Pass=as.numeric(rh[colnames(chain_mats[[1L]])]) <= MCMC_RHAT_LIMIT,
      ESS_Pass=as.numeric(ess[colnames(chain_mats[[1L]])]) >= MCMC_ESS_MIN,
      Iterations=as.integer(settings[['iter']]),
      Burn_in=as.integer(settings[['burn']]),
      row.names=NULL
    )
    final <- list(chains=chains, diagnostics=diagnostics,
                  chain_mats=chain_mats, settings=settings)
    failed <- anyNA(diagnostics$Rhat) || anyNA(diagnostics$ESS) ||
      any(!diagnostics$Rhat_Pass | !diagnostics$ESS_Pass, na.rm=TRUE)
    if (!failed) {
      message(model, ': all monitored parameters passed R-hat <= ', MCMC_RHAT_LIMIT,
              ' and ESS >= ', MCMC_ESS_MIN, ' at tau=', tau, '.')
      break
    }
    if (stage < length(stage_settings)) {
      message('WARNING: ', model, ' diagnostics did not pass at tau=', tau,
              '. Re-running with longer chains; failed parameters are not being hidden.')
    } else {
      warning(model, ' at tau=', tau,
              ' still has failed convergence diagnostics after the maximum chain length. Inspect trace plots and do not report these posterior summaries as fully converged.')
    }
  }
  chains <- final$chains
  diagnostics <- final$diagnostics
  all_draws <- do.call(rbind, lapply(chains, `[[`, 'draws'))
  combined <- chains[[1L]]
  combined$draws <- all_draws
  combined$coef <- colMeans(all_draws)
  combined$lower <- apply(all_draws, 2L, quantile, probs=0.025)
  combined$upper <- apply(all_draws, 2L, quantile, probs=0.975)
  combined$sigma_draws <- unlist(lapply(chains, `[[`, 'sigma_draws'))
  combined$lambda1sq_draws <- unlist(lapply(chains, `[[`, 'lambda1sq_draws'))
  if (model == 'BFLQR')
    combined$lambda2sq_draws <- unlist(lapply(chains, `[[`, 'lambda2sq_draws'))
  list(fit=combined, diagnostics=diagnostics, chains=chains,
       final_iterations=final$settings[['iter']], final_burn=final$settings[['burn']])
}

# Fit the models at each quantile
coef_out<-list()
pred_out<-list()
tune_out<-list()
mcmc_out<-list()
convergence_out<-list()
bootstrap_out<-list()
for(tau in TAU){message('Fitting tau=',tau)
l<-fitLQR(Xt,yt,tau)
bl_chains <- run_bayesian_chains(fitBLQR, Xt, yt, tau, "BLQR", nchains=4L)
bl <- bl_chains$fit
bfl_chains <- run_bayesian_chains(fitBFLQR, Xt, yt, tau, "BFLQR", nchains=4L)
b <- bfl_chains$fit
convergence_out[[as.character(tau)]] <- rbind(bl_chains$diagnostics,
                                               bfl_chains$diagnostics)
f<-fitFLQR(Xt,yt,tau)
lc<-bootCI(Xt,yt,tau,'LQR',l$tuning)
fc<-bootCI(Xt,yt,tau,'FLQR',f$tuning)
bootstrap_out[[as.character(tau)]]<-data.frame(
  Tau=tau, Model=c('LQR','FLQR'), Requested_Replications=c(BOOT_B,BOOT_B),
  Successful_Replications=c(lc$successful,fc$successful),
  Block_Length=c(lc$block_length,fc$block_length),
  row.names=NULL
)
 # Coefficients and intervals remain on the supplied standardized scale.
 llo<-lc$lower
lhi<-lc$upper
flo<-fc$lower
fhi<-fc$upper
 rowfun<-function(model,est,lo,hi,kind){sel<-rep(NA,length(est))
sel[-1]<-lo[-1]>0|hi[-1]<0
data.frame(Tau=tau,Model=model,Parameter=PN,Estimate=est,Lower=lo,Upper=hi,Interval=kind,Selected=sel)}
 coef_out[[as.character(tau)]]<-rbind(rowfun('LQR',l$coef,llo,lhi,'Moving-block bootstrap CI'),
    rowfun('BLQR',bl$coef,bl$lower,bl$upper,'Posterior credible interval'),
    rowfun('FLQR',f$coef,flo,fhi,'Moving-block bootstrap CI'),rowfun('BFLQR',
    b$coef,b$lower,b$upper,'Posterior credible interval'))
 pl<-pred(Xv,l$coef)
pbl<-pred(Xv,bl$coef)
pf<-pred(Xv,f$coef)
pb<-pred(Xv,b$coef)
pred_out[[as.character(tau)]]<-data.frame(Tau=tau,Model=c('LQR','BLQR',
    'FLQR','BFLQR'),Test_QL=c(QL(yv,pl,tau),QL(yv,pbl,tau),QL(yv,pf,
    tau),QL(yv,pb,tau)),Test_MAE=c(MAE(yv,pl),MAE(yv,pbl),MAE(yv,pf),
    MAE(yv,pb)))
 tune_out[[as.character(tau)]]<-data.frame(
    Tau=tau, LQR_CV_Penalty=l$tuning,
    FLQR_CV_Lambda1=unname(f$tuning['lambda1']),
    FLQR_CV_Lambda2=unname(f$tuning['lambda2'])
  )
mcmc_out[[as.character(tau)]]<-rbind(data.frame(Tau=tau,Model='BLQR',
    Parameter=c('sigma','lambda1sq'),Posterior_Mean=c(mean(bl$sigma_draws),
    mean(bl$lambda1sq_draws)),Posterior_SD=c(sd(bl$sigma_draws),sd(bl$lambda1sq_draws))),
    data.frame(Tau=tau,Model='BFLQR',Parameter=c('sigma','lambda1sq',
    'lambda2sq'),Posterior_Mean=c(mean(b$sigma_draws),mean(b$lambda1sq_draws),
    mean(b$lambda2sq_draws)),Posterior_SD=c(sd(b$sigma_draws),sd(b$lambda1sq_draws),
    sd(b$lambda2sq_draws))))}
T1<-do.call(rbind,coef_out)
T1$Result<-sprintf('%.4f (%.4f, %.4f)',T1$Estimate,T1$Lower,T1$Upper)
T1$Selection<-ifelse(is.na(T1$Selected),NA,ifelse(T1$Selected,'Yes','No'))
T1<-T1[,c('Tau','Model','Parameter','Result','Interval','Selection')]
T2<-do.call(rbind,pred_out)
T2$Test_QL<-round(T2$Test_QL,4)
T2$Test_MAE<-round(T2$Test_MAE,4)
TUNE<-do.call(rbind,tune_out)
MCMC_SUMMARY<-do.call(rbind,mcmc_out)
MCMC_DIAGNOSTICS<-do.call(rbind,convergence_out)
BOOTSTRAP_SUMMARY<-do.call(rbind,bootstrap_out)
# Summarize interval-based selection for explanatory variables only; the
# intercept is intentionally excluded. This is not an exact lasso-zero rule.
SELECTION_SUMMARY<-aggregate(Selection ~ Tau + Model,
  data=T1[T1$Parameter!='Intercept', c('Tau','Model','Parameter','Selection')],
  FUN=function(z) sum(z == 'Yes', na.rm=TRUE))
names(SELECTION_SUMMARY)[names(SELECTION_SUMMARY)=='Selection'] <- 'Variables_Whose_95CI_Excludes_Zero'
SELECTION_SUMMARY$Total_Predictors <- length(VarNames)
SELECTION_SUMMARY <- SELECTION_SUMMARY[order(SELECTION_SUMMARY$Tau, SELECTION_SUMMARY$Model),]
cat('\nPARAMETER ESTIMATES AND INTERVAL-BASED SELECTION\n')
print(T1, row.names=FALSE)
cat('\nSELECTION SUMMARY (predictors only; interval excludes zero)\n')
print(SELECTION_SUMMARY, row.names=FALSE)
cat('\nHOLDOUT PREDICTIVE PERFORMANCE\n')
print(T2, row.names=FALSE)
cat('\nCROSS-VALIDATION-SELECTED PENALTIES (LQR/FLQR)\n')
print(TUNE, row.names=FALSE)
cat('\nBAYESIAN GLOBAL-SHRINKAGE SUMMARIES ARE POSTERIOR ESTIMATES, NOT CV-TUNED PENALTIES\n')
cat('\nMOVING-BLOCK BOOTSTRAP SUCCESS COUNTS\n')
print(BOOTSTRAP_SUMMARY, row.names=FALSE)
cat('\nBAYESIAN POSTERIOR SUMMARIES\n')
print(MCMC_SUMMARY, row.names=FALSE)
cat('\nMCMC CONVERGENCE DIAGNOSTICS (4 independent chains; adaptive extension if needed)\n')
print(MCMC_DIAGNOSTICS, row.names=FALSE)
cat('Convergence thresholds: R-hat <= ', MCMC_RHAT_LIMIT, '; ESS >= ', MCMC_ESS_MIN, '\n', sep='')
if (anyNA(MCMC_DIAGNOSTICS$Rhat) || anyNA(MCMC_DIAGNOSTICS$ESS) ||
    any(!MCMC_DIAGNOSTICS$Rhat_Pass | !MCMC_DIAGNOSTICS$ESS_Pass, na.rm=TRUE)) {
  cat('WARNING: Some R-hat/ESS values are missing or failed the chosen thresholds even after adaptive extension; inspect chains before reporting results.\n')
} else {
  cat('All reported R-hat and ESS checks passed the chosen thresholds.\n')
}
cat('\nAnalysis completed. All results are displayed in the console.\n')
