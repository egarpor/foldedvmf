
# Reproduces the empirical results of "The folded von Mises--Fisher
# distribution" (García-Portugués and Guardiola): Figures 1-3, Tables 1-2,
# and the numbers quoted in Sections 2, 3, and 5.

## Setup
{

# Single-threaded vecLib BLAS on macOS to avoid segfaults in parallel::mclapply
Sys.setenv(VECLIB_MAXIMUM_THREADS = "1")

# Cores for the parallel parts: all but one
cores <- max(1, parallel::detectCores() - 1)

# Needed packages
library(DirStats)
library(rotasym)
library(Matrix)
library(movMF)
library(slam)
library(parallel)
library(tools)
library(plot3D)
library(viridisLite)
library(magick)
stopifnot(packageVersion("DirStats") >= "1.0.0")
stopifnot(packageVersion("rotasym") >= "1.4.0")

# Packages with the Section 5 data
data_pkgs <- c("ade4", "analogue", "archdata", "coda4microbiome",
               "compositions", "microbiome", "mvoutlier", "pdfCluster",
               "phyloseq", "robCompositions", "zCompositions")
stopifnot(sapply(data_pkgs, function(pkg) nzchar(system.file(package = pkg))))

# Log with every number quoted in the paper
out <- "folded-vmf.txt"
invisible(file.create(out))
say <- function(...) {

  cat(sprintf(...), "\n", sep = "", file = out, append = TRUE)

}
say_df <- function(x) {

  say("%s", paste(capture.output(print(x, row.names = FALSE)),
                  collapse = "\n"))

}

}

## Functions
{

## Utilities
{

# Unit vector along |a|
unit_pos <- function(a) abs(a) / sqrt(sum(a^2))

# Row-wise log(sum(exp(L))) in a stable way
log_sum_exp <- function(L) {

  L <- rbind(L)
  m <- apply(L, 1, max)
  m + log(rowSums(exp(L - m)))

}

# log cosh(t) in a stable way
log_cosh <- function(t) {

  t <- abs(t)
  t + log1p(exp(-2 * t)) - log(2)

}

# Sizes of batches of at most B adding up to n
batches <- function(n, B) c(rep(B, n %/% B), if (n %% B > 0) n %% B)

}

## vMF density and fit
{

# Log-density of vMF(mu, kappa) in a stable way
log_f_vmf <- function(X, mu, kappa) {

  DirStats::log_c_vmf(q = length(mu) - 1, kappa = kappa) +
    kappa * as.numeric(X %*% mu)

}

# A_p(kappa) = I_{p/2}(kappa) / I_{p/2-1}(kappa) in a stable way
A_p <- function(kappa, p) {

  exp(DirStats::log_besselI_scaled(nu = p / 2, x = kappa) -
        DirStats::log_besselI_scaled(nu = p / 2 - 1, x = kappa))

}

# vMF fit by weighted maximum likelihood
fit_vmf <- function(X, w = rep(1, nrow(X))) {

  # Normalized weighted mean and kappa solving A_p(kappa) = R_weighted
  m <- colSums(w * X)
  R <- sqrt(sum(m^2))
  list(mu = m / R, kappa = DirStats::kappa_ml(data = X, R = R / sum(w),
                                              tol = 1e-12))

}

}

## vMF+ density, fit, mode, and sampling
{

# Stable log-density of vMF+(mu, kappa), product-cosh form, X dense or sparse
log_f_vmfp <- function(X, mu, kappa) {

  Y <- X %*% Matrix::Diagonal(x = kappa * mu)
  Y@x <- log_cosh(t = Y@x)
  DirStats::log_c_vmf(q = length(mu) - 1, kappa = kappa) +
    length(mu) * log(2) + Matrix::rowSums(Y)

}

# vMF+ fit by weighted maximum likelihood in xi = kappa * mu parametrization,
# started at xi0 (by default, the vMF fit), for dense or sparse X.
fit_vmfp <- function(X, w = rep(1, nrow(X)), xi0 = NULL) {

  # Obtain vMF fit
  p <- ncol(X)
  w <- w / sum(w)
  if (is.null(xi0)) {

    vmf_fit <- fit_vmf(X = X, w = w)
    xi0 <- vmf_fit$kappa * vmf_fit$mu

  }

  # Minus the weighted log-likelihood in xi and minus its score (Proposition 2),
  # on the nonzero entries of X (column j) only, as log cosh(0) = tanh(0) = 0.
  X <- as(X, "CsparseMatrix")
  j <- rep.int(seq_len(p), diff(X@p))
  wi <- w[X@i + 1]
  Y <- X
  fn <- function(xi) {

    kappa <- sqrt(sum(xi^2))
    -(DirStats::log_c_vmf(q = p - 1, kappa = kappa) + p * log(2) +
        sum(wi * log_cosh(t = xi[j] * X@x)))

  }
  gr <- function(xi) {

    kappa <- sqrt(sum(xi^2))
    Y@x <- wi * tanh(xi[j] * X@x) * X@x
    -(Matrix::colSums(Y) - A_p(kappa = kappa, p = p) * xi / kappa)

  }

  # Tight tolerance (factr = 1e3), at which L-BFGS-B may stop in a line search
  # at machine precision, so convergence is judged by the score vanishing.
  opt <- optim(par = xi0, fn = fn, gr = gr, method = "L-BFGS-B",
               control = list(maxit = 1000, factr = 1e3))
  list(mu = unit_pos(a = opt$par), kappa = sqrt(sum(opt$par^2)),
       converged = max(abs(gr(opt$par))) <= 1e-5)

}

# Mode m of vMF+(mu, kappa) by Proposition 6
mode_vmfp <- function(mu, kappa, tol = 1e-12) {

  # xi and maximum
  xi <- kappa * mu
  xi_max <- max(xi)

  # r_j(lambda): 0 if xi_j^2 <= lambda, else the root in (0, 1] of
  # xi_j * tanh(xi_j * r) / r = lambda, whose left-hand side is xi_j^2 at r = 0.
  r <- function(lambda) {

    sapply(xi, function(xi_j) {

      if (xi_j^2 <= lambda) return(0)
      uniroot(function(r) ifelse(r == 0, xi_j^2, xi_j * tanh(xi_j * r) / r) -
                lambda, interval = c(0, 1), tol = tol)$root

    })

  }

  # Solve sum_j r_j(lambda)^2 = 1 on [xi_max * tanh(xi_max), xi_max^2)
  lambda_m <- uniroot(function(lambda) sum(r(lambda = lambda)^2) - 1,
                      interval = c(xi_max * tanh(xi_max), xi_max^2),
                      tol = tol)$root

  # Mode
  unit_pos(a = r(lambda = lambda_m))

}

# Sampling from vMF+(mu, kappa) by folding a vMF sample
r_vmfp <- function(n, mu, kappa) {

  abs(rotasym::r_vMF(n = n, mu = mu, kappa = kappa))

}

}

## vMF orthant mass by importance sampling
{

# Find a > 0 solving a^2 E||V||^2 = kappa^2, for V ~ N_p(a * mu, I_p)
# truncated to R^p_+, where
# E||V||^2 = p + a^2 + a * sum_j mu_j * phi(a * mu_j) / Phi(a * mu_j).
a_tn <- function(kappa, mu) {

  # f(a) = a^2 E||V||^2 - kappa^2, increasing in a, so its root is unique
  p <- length(mu)
  f <- function(a) {

    a^2 * (p + a^2 + a * sum(mu * exp(dnorm(a * mu, log = TRUE) -
                                        pnorm(a * mu, log.p = TRUE)))) -
      kappa^2

  }

  # Upper bound a0, the root without the truncation (E||V||^2 = p + a^2; the
  # truncation only adds to it), written without cancellation for small kappa.
  a0 <- kappa * sqrt(2 / (sqrt(p^2 + 4 * kappa^2) + p))
  uniroot(f, lower = 0, upper = a0, extendInt = "upX", tol = 1e-10 * a0)$root

}

# log m_p(s) in the projected normal exactly computed by the recurrence
# m_{j + 1} = s * m_j + (j - 1) * m_{j - 1}, with
# m_1(s) = sqrt(2 * pi) * exp(s^2 / 2) * Phi(s) and m_2(s) = 1 + s * m_1(s).
log_m <- function(s, p) {

  # The m_j overflow, so the recurrence runs on
  # rho_j = m_{j + 1} / m_j = s + (j - 1) * m_{j - 1} / m_j
  #       = s + (j - 1) / rho_{j - 1},
  # accumulating log m_{j + 1} = log m_j + log rho_j.

  # l = log m_1, rho_1 = m_2 / m_1 = s + 1 / m_1, and l = log m_2
  l <- log(2 * pi) / 2 + s^2 / 2 + pnorm(s, log.p = TRUE)
  rho <- s + exp(-l)
  l <- l + log(rho)

  # Recurrence
  for (j in seq_len(p - 2) + 1) {

    rho <- s + (j - 1) / rho
    l <- l + log(rho)

  }
  l

}

# Log-density on S^{p-1}_+ of X = V / ||V||, for V ~ N_p(a * mu, I_p)
# truncated to R^p_+ and a = a_tn(kappa, mu),
# (2 * pi)^(-p / 2) * e^(-a^2 / 2) * m_p(a * mu'x) / prod_j Phi(a * mu_j).
log_d_tnorm_proj <- function(x, mu, kappa) {

  p <- length(mu)
  a <- a_tn(kappa = kappa, mu = mu)
  -p / 2 * log(2 * pi) - a^2 / 2 - sum(pnorm(a * mu, log.p = TRUE)) +
    log_m(s = a * as.numeric(x %*% mu), p = p)

}

# Simulates X = V / ||V||, of independent V_j ~ N(a * mu_j, 1) truncated to
# (0, Inf), with a = a_tn(kappa, mu). V_j is drawn by inverting its survival
# function, Phi(a * mu_j - v) / Phi(a * mu_j) = U.
r_tnorm_proj <- function(n, mu, kappa) {

  p <- length(mu)
  a <- a_tn(kappa = kappa, mu = mu)
  V <- t(a * mu - qnorm(matrix(runif(n * p), nrow = p) * pnorm(a * mu)))
  V / sqrt(rowSums(V^2))

}

# log P(Z in S^{p-1}_+) for Z ~ vMF(mu, kappa), mu >= 0, by importance sampling
# (Section 5.2): P = E{f_vmf(X) / d_tnorm_proj(X)}, X ~ d_tnorm_proj, with
# a = a_tn() making the weights nearly constant. M draws in batches.
log_P_orthant <- function(mu, kappa, M = 5000, seed = 987202226) {

  # Draws in batches of B, so that each B x p matrix of draws has at most 1e7
  # entries, bounding the memory in high dimension.
  set.seed(seed)
  B <- max(1, floor(1e7 / length(mu)))
  lw <- unlist(lapply(batches(n = M, B = B), function(b) {
    X <- r_tnorm_proj(n = b, mu = mu, kappa = kappa)
    log_f_vmf(X = X, mu = mu, kappa = kappa) -
      log_d_tnorm_proj(x = X, mu = mu, kappa = kappa)
  }))
  as.numeric(log_sum_exp(L = lw)) - log(M)

}

}

## Mixtures fitting
{

# Expectation--maximization for an MvMF+, started at init = list(pi, pars) (an
# MvMF fit). Both tolerances are relative to the log-likelihood to avoid
# rounding errors for large datasets. Returns converged = FALSE if the
# iteration cap was reached.
em_vmfp <- function(X, init, iters = 2000, tol = 1e-10) {

  # Initial checks
  if (is.null(init)) stop("EM needs an initial MvMF fit.")
  pis <- init$pi
  pars <- init$pars
  old <- -Inf
  converged <- FALSE

  # EM algorithm
  for (it in seq_len(iters)) {

    # E-step: the weighted log-likelihoods
    # wloglik_terms[i, k] = log pi_k + log f_k(X_i), the log-likelihood
    # loglik_i of each X_i, and the posterior probabilities
    # post[i, k] = pi_k * f_k(X_i) / f(X_i) of the components.
    wloglik_terms <- sapply(seq_along(pars), function(k) {
      log(pis[k]) + log_f_vmfp(X = X, mu = pars[[k]]$mu,
                               kappa = pars[[k]]$kappa)
    })
    loglik_i <- log_sum_exp(L = wloglik_terms)
    loglik <- sum(loglik_i)
    post <- exp(wloglik_terms - loglik_i)

    # Stop at a non-finite or decreased log-likelihood, as the M-step cannot
    # decrease it, and break when it increased by less than tol * |loglik|.
    if (!is.finite(loglik)) stop("EM log-likelihood not finite.")
    if (loglik < old - tol * abs(loglik)) stop("EM log-likelihood decreased.")
    converged <- loglik - old < tol * abs(loglik)
    if (converged) break
    old <- loglik

    # M-step: weights and weighted fits, warm-started
    pis <- colMeans(post)
    pars <- lapply(seq_along(pars), function(k) {
      fit_vmfp(X = X, w = post[, k], xi0 = pars[[k]]$kappa * pars[[k]]$mu)
    })

  }
  list(loglik = loglik, pi = pis, pars = pars, converged = converged)

}

# The MvMF fit at K for dense or sparse X out of n_start movMF runs from random
# prototypes. Runs that are degenerate (kappa > 1e10) or end in an error are
# skipped.
vmf_K <- function(X, K, cores, n_start = 25, seed = 20261016) {

  # movMF takes the data as a slam triplet matrix, built from a general matrix,
  # as Matrix may store a square symmetric X by only one of its triangles.
  Xs <- as(as(X, "generalMatrix"), "TsparseMatrix")
  Xs <- slam::simple_triplet_matrix(i = Xs@i + 1L, j = Xs@j + 1L, v = Xs@x,
                                    nrow = nrow(X), ncol = ncol(X))

  # n_start single-run calls, as movMF keeps its best run even if degenerate
  runs <- parallel::mclapply(seq_len(n_start), function(j) {
    set.seed(seed + 100 * K + j)
    tryCatch(movMF::movMF(Xs, k = K, control = list(nruns = 1,
                                                    maxiter = 2000)),
             error = conditionMessage)
  }, mc.cores = cores, mc.preschedule = FALSE)

  # A run that crashed (e.g., out of memory) returns NULL or a try-error
  if (any(sapply(runs, function(r) is.null(r) || inherits(r, "try-error")))) {

    stop("A movMF run crashed at K = ", K, ".")

  }

  # Largest kappa of each run, NA if it ended in an error; runs with an error
  # or with kappa > 1e10 are excluded, and counted.
  kap_max <- sapply(runs, function(r) {
    if (is.character(r)) NA else max(sqrt(rowSums(r$theta^2)))
  })
  ok <- !is.na(kap_max) & kap_max <= 1e10
  health <- c(runs = n_start, failed_vmf = sum(is.na(kap_max)),
              degenerate_vmf = sum(kap_max > 1e10, na.rm = TRUE),
              capped_vmf = FALSE)
  if (!any(ok)) return(list(fit = NULL, health = health))

  # Best of the remaining runs, as weights pi and components (mu, kappa). Runs
  # that reach the iteration cap without converging remain in the pool of
  # candidates; the selected fit is flagged if it did.
  fm <- runs[ok][[which.max(sapply(runs[ok], logLik))]]
  health[["capped_vmf"]] <- fm$details$iter[["iter"]] >= 2000
  kap <- unname(sqrt(rowSums(fm$theta^2)))
  vmf_fit <- list(pi = unname(fm$alpha), pars = lapply(seq_len(K), function(k) {
    list(mu = unname(fm$theta[k, ]) / kap[k], kappa = kap[k])
  }))

  # Log-likelihood with the density of log_f_vmf, as that of movMF is with
  # respect to another reference measure (a constant per dataset).
  vmf_fit$loglik <- sum(log_f_mix(X = X, fit = vmf_fit, log_dens = log_f_vmf))

  # Return fit and diagnostics
  list(fit = vmf_fit, health = health)

}

# MvMF fits at K = 1, 2, 3, 4, 6, 8, 12, 16, ... (powers of 2 and their
# midpoints) until the MvMF BIC increases, a K has no fit, or K would exceed
# n / 5. Returns the MvMF fit at K_hat, K_hat, and diagnostics of every K.
scan_K <- function(X, cores) {

  # Grid of K: powers of 2 and their midpoints, up to n / 5
  n <- nrow(X)
  p <- ncol(X)
  K_all <- sort(c(2^(0:12), 1.5 * 2^(1:12)))
  K_all <- K_all[K_all <= n / 5]
  bic <- Inf
  fit <- health <- NULL
  for (K in K_all) {

    # Fit at K and its safeguard counts; stop if no run is left
    r <- vmf_K(X = X, K = K, cores = cores)
    health <- cbind(health, r$health)
    if (is.null(r$fit)) break

    # BIC with K * p + K - 1 parameters; stop if it increases, else keep K
    b <- -2 * r$fit$loglik + (K * p + K - 1) * log(n)
    if (b > bic) break
    bic <- b
    fit <- r$fit
    K_hat <- K

  }

  # At least K = 1 has a fit
  stopifnot(!is.null(fit))

  # Return fit, K_hat, and diagnostics
  list(fit = fit, K_hat = K_hat, health = health)

}

# The three mixtures at the K of an MvMF fit vmf_fit: that fit, the MrvMF
# (that fit divided once by its orthant mass), and the MvMF+ (em_vmfp from
# it). Returns the fits, their log-likelihoods, whether the bound below forces
# the gain of the MvMF+ over the MrvMF, and the MvMF+ safeguard counts.
trio_K <- function(X, vmf_fit) {

  # MvMF+ fit, giving an error if degenerate
  vmfp_fit <- em_vmfp(X = X, init = vmf_fit)
  if (max(sapply(vmfp_fit$pars, `[[`, "kappa")) > 1e10) {

    stop("Degenerate MvMF+ fit.")

  }

  # The MvMF+ log-likelihood is the MvMF one plus at least log(2) * k_zero
  # (k_zero the number of zeros), and EM from the MvMF fit does not decrease it.
  n <- nrow(X)
  k_zero <- n * ncol(X) - sum(X@x > 0)
  stopifnot(vmfp_fit$loglik + 1e-6 >= vmf_fit$loglik + log(2) * k_zero)

  # MrvMF "fit", updating the MvMF log-likelihood minus n * log P_K
  log_p <- sapply(vmf_fit$pars, function(q) {
    log_P_orthant(mu = q$mu, kappa = q$kappa)
  })
  log_P <- as.numeric(log_sum_exp(L = log(vmf_fit$pi) + log_p))

  # The gain of the MvMF+ over the MrvMF is forced to be nonnegative if
  # log(2) * k_zero + n * log P_K >= 0, up to a rounding tolerance.
  forced_gain <- log(2) * k_zero + n * log_P + 1e-6 >= 0

  # Return trio
  list(fit_vmf = vmf_fit, fit_vmfp = vmfp_fit,
       loglik = c(vmf = vmf_fit$loglik, rvmf = vmf_fit$loglik - n * log_P,
                  vmfp = vmfp_fit$loglik),
       forced_gain = forced_gain,
       health = c(capped_vmfp = !vmfp_fit$converged))

}

}

## Mixture regions
{

# Log-density at the rows of X, dense or sparse, of a mixture fit with
# component log-density function log_dens.
log_f_mix <- function(X, fit, log_dens) {

  log_sum_exp(L = matrix(sapply(seq_along(fit$pars), function(k) {
    log(fit$pi[k]) + log_dens(X, mu = fit$pars[[k]]$mu,
                              kappa = fit$pars[[k]]$kappa)
  }), nrow = nrow(X)))

}

# Log-volume of A = {log f >= cut}, f = sum_k pi_k f_k a mixture with log(f_k)
# given by the function log_dens, by importance sampling (Section 5.3):
# vol(A) = E{1_A(Y) / q(Y)}, Y ~ q. q mixes equally the q_k ~ r_q, with
# concentration eta*_k = beta_k * kappa_k and N / K draws. beta_k is found by
# bisection so that ~50% of draws of q_k fall in {log pi_k + log f_k >= cut}.
log_vol <- function(cut, fit, log_dens, r_q, log_dens_q, N, B) {

  # Proposal with equal weights and components q_k with direction mu_k and
  # concentration eta*_k = beta_k * kappa_k, with beta_k searched for.
  K <- length(fit$pars)
  proposal <- list(pi = rep(1 / K, K), pars = lapply(seq_len(K), function(k) {

    # Search for beta_k
    par_k <- fit$pars[[k]]
    half <- function(log2_beta) {

      Y <- r_q(n = 400, mu = par_k$mu, kappa = 2^log2_beta * par_k$kappa)
      log_dens_k <- log(fit$pi[k]) +
        log_dens(X = Y, mu = par_k$mu, kappa = par_k$kappa)
      mean(log_dens_k >= cut) - 0.5

    }
    bracket <- c(-30, 10)
    for (i in 1:14) {

      mid <- mean(bracket)
      if (half(log2_beta = mid) < 0) bracket[1] <- mid else bracket[2] <- mid

    }
    beta_k <- 2^mean(bracket)

    # Component parameters
    list(mu = par_k$mu, kappa = beta_k * par_k$kappa)

  }))

  # Compute log f(Y) and log q(Y) with Y ~ q
  log_fq <- do.call(rbind, lapply(proposal$pars, function(q_k) {
    do.call(rbind, lapply(batches(n = N %/% K, B = B), function(b) {
      Y <- r_q(n = b, mu = q_k$mu, kappa = q_k$kappa)
      cbind(log_f = log_f_mix(X = Y, fit = fit, log_dens = log_dens),
            log_q = log_f_mix(X = Y, fit = proposal, log_dens = log_dens_q))
    }))
  }))

  # Log of the mean of 1_A(Y) / q(Y) over the draws
  inside <- log_fq[, "log_f"] >= cut
  if (!any(inside)) stop("Empty region.")
  as.numeric(log_sum_exp(L = -log_fq[inside, "log_q"])) - log(nrow(log_fq))

}

# Regions at level 1 - alpha of the trio_K fits on X_tr: R (MvMF+, on
# S^{p-1}_+), C (naive MvMF, on S^{p-1}), and C+ (MrvMF, on S^{p-1}_+).
# Thresholds are: (1) plug-in (Hyndman, 1996), computed from X_tr and
# evaluated on (X_ca, X_te); (2) conformal (Lei et al., 2015), computed from
# X_ca and evaluated on X_te.
regions_split <- function(X_tr, X_ca, X_te, fit, N = 20000, alpha = 0.10,
                          cores = 1) {

  # Fitted log-densities of R, C, and C+ at X_tr, X_ca, and X_te
  vmfp_fit <- fit$fit_vmfp
  vmf_fit <- fit$fit_vmf
  log_dens_mix <- function(X) {

    # As the MrvMF log-density is the MvMF one minus log P_K, a threshold from
    # data on the orthant makes C+ = C intersected with S^{p-1}_+, so C+ uses
    # the MvMF log-density, with the threshold and coverage of C.
    log_dens_C <- log_f_mix(X = X, fit = vmf_fit, log_dens = log_f_vmf)
    list(R = log_f_mix(X = X, fit = vmfp_fit, log_dens = log_f_vmfp),
         C = log_dens_C, Cp = log_dens_C)

  }
  log_dens_tr <- log_dens_mix(X = X_tr)
  log_dens_ca <- log_dens_mix(X = X_ca)
  log_dens_te <- log_dens_mix(X = X_te)

  # Thresholds
  kth <- function(l, k) {

    stopifnot(k >= 1)
    sort(l)[k]

  }
  cuts <- list(plugin = sapply(log_dens_tr, function(l) {
    kth(l = l, k = floor(length(l) * alpha))
  }), conformal = sapply(log_dens_ca, function(l) {
    kth(l = l, k = floor((length(l) + 1) * alpha))
  }))

  # Volume draws in batches of B, so that each B x p matrix of draws has at
  # most 5e6 entries, bounding the memory in high dimension.
  B <- max(1, floor(5e6 / ncol(X_tr)))

  # Volumes by importance sampling, parallelized by region and threshold.
  # Proposals: vMF+ for R, vMF for C, and r_tnorm_proj for C+.
  tasks <- expand.grid(r = c("R", "C", "Cp"), j = names(cuts),
                       stringsAsFactors = FALSE)
  res <- parallel::mclapply(seq_len(nrow(tasks)), function(i) {
    set.seed(20261016 + i)
    cut <- cuts[[tasks$j[i]]][[tasks$r[i]]]
    switch(tasks$r[i],
           R = log_vol(cut = cut, fit = vmfp_fit, log_dens = log_f_vmfp,
                       r_q = r_vmfp, log_dens_q = log_f_vmfp, N = N, B = B),
           C = log_vol(cut = cut, fit = vmf_fit, log_dens = log_f_vmf,
                       r_q = rotasym::r_vMF, log_dens_q = log_f_vmf,
                       N = N, B = B),
           Cp = log_vol(cut = cut, fit = vmf_fit, log_dens = log_f_vmf,
                        r_q = r_tnorm_proj, log_dens_q = log_d_tnorm_proj,
                        N = N, B = B))
  }, mc.cores = min(cores, nrow(tasks)), mc.preschedule = FALSE)
  bad <- sapply(res, inherits, "try-error")
  if (any(bad)) stop("A volume failed: ", paste(unlist(res[bad]),
                                                collapse = "; "))

  # Coverage on held-out data: (X_ca, X_te) for plug-in, and X_te for conformal
  log_dens_cover <- list(plugin = Map(c, log_dens_ca, log_dens_te),
                         conformal = log_dens_te)

  # Coverage and log-volumes by threshold
  out <- lapply(names(cuts), function(j) {
    lv <- setNames(unlist(res[tasks$j == j]), c("R", "C", "Cp"))
    list(cover = mapply(function(l, cut) mean(l >= cut), log_dens_cover[[j]],
                        cuts[[j]]), lv = lv)
  })
  names(out) <- names(cuts)
  list(plugin = out$plugin, conformal = out$conformal)

}

}

}

## Figure 1: estimation bias
{

# Mean ratios kappa_hat / kappa0 of the vMF and vMF+ fits
bias_cell <- function(mu0, kappa0, n = 1e3, M = 1e3) {

  k <- replicate(M, {
    X <- r_vmfp(n = n, mu = mu0, kappa = kappa0)
    c(kratio_vmf = fit_vmf(X = X)$kappa, kratio_vmfp = fit_vmfp(X = X)$kappa)
  })
  rowMeans(k) / kappa0

}

# Simulation design
ps <- c(3, 5, 10)
kappas <- exp(seq(log(1), log(500), length.out = 15))
mus <- lapply(ps, function(p) {
  z <- floor(p / 2)
  list(interior = unit_pos(a = 1:p),
       boundary = unit_pos(a = rep(0:1, c(z, p - z))))
})
names(mus) <- ps

# Quoted in Section 2.1: the limit (p - 1) / (p - 1 - 2 * k / pi) of the vMF
# ratio, for the boundary mu0 (k = floor(p / 2)) and for mu0 = e_j (k = p - 1).
say("## Figure 1: limits of the vMF ratio kappa* / kappa0")
plateau <- (ps - 1) / (ps - 1 - 2 * floor(ps / 2) / pi)
say("Boundary mu0: %s; mu0 = e_j: %.2f",
    paste(sprintf("%.2f (p = %g)", plateau, ps), collapse = ", "),
    1 / (1 - 2 / pi))

# Monte Carlo results; delete the results to recompute
results <- "folded-vmf-bias.rds"
if (file.exists(results)) {

  res <- readRDS(results)

} else {

  # Cells in parallel
  cells <- expand.grid(kappa0 = kappas, mu = c("interior", "boundary"),
                       p = ps, stringsAsFactors = FALSE)
  res <- parallel::mclapply(seq_len(nrow(cells)), function(i) {
    set.seed(20260530 + i)
    bias_cell(mu0 = mus[[as.character(cells$p[i])]][[cells$mu[i]]],
              kappa0 = cells$kappa0[i])
  }, mc.cores = cores, mc.preschedule = FALSE)
  res <- cbind(cells, do.call(rbind, res))
  res <- res[order(res$p, res$mu, res$kappa0), ]
  saveRDS(res, file = results)

}

# Square panels, equal narrow left margins, y-axis on the first only
# with its labels in the left outer margin. Height 3.35 in nearly fills the
# width, (8.6 - 7.4 * lh) / 3 + 5.6 * lh = 3.36 in for a margin line of
# lh = 0.158 in (cex = 0.66 at pointsize 14.4). Reference lines at one and at
# the boundary-mu0 limit (p - 1) / (p - 1 - 2 * k / pi) of the vMF ratio,
# k = floor(p / 2) zeros.
yl <- range(c(res$kratio_vmf, res$kratio_vmfp, 0.97))
xt <- c(1, 2, 5, 10, 50, 200, 500)
pchs <- c(17, 2, 16, 1)
cols <- rep(c("steelblue", "firebrick"), each = 2)
ltys <- c(1, 3, 1, 3)
pdf("../img/app-kappa-bias.pdf", width = 8.6, height = 3.35, pointsize = 14.4)
par(mfrow = c(1, 3), oma = c(0, 3.5, 0, 0), mar = c(3.6, 0.5, 2.0, 0.8),
    mgp = c(2.2, 0.7, 0), pty = "s")
for (p in ps) {

  r <- res[res$p == p, ]
  b <- r$mu == "boundary"
  matplot(r$kappa0[b], with(r, cbind(kratio_vmf[!b], kratio_vmf[b],
                                     kratio_vmfp[!b], kratio_vmfp[b])),
          type = "o", pch = pchs, col = cols, lty = ltys, lwd = 2, log = "xy",
          ylim = yl, xaxt = "n", yaxt = if (p == ps[1]) "s" else "n",
          xlab = expression(kappa[0]), ylab = "", main = bquote(p == .(p)),
          cex.main = 1.1)
  axis(1, at = xt, gap.axis = -1)
  abline(h = c(1, (p - 1) / (p - 1 - 2 * floor(p / 2) / pi)), col = "grey60",
         lty = 3)
  if (p == ps[1]) {

    legend("topright", legend = c("vMF MLE, interior", "vMF MLE, boundary",
                                  "vMF+ MLE, interior", "vMF+ MLE, boundary"),
           pch = pchs, col = cols, lty = ltys, lwd = 2, bty = "n", cex = 0.85)

  }

}
mtext(expression("Mean " * hat(kappa) / kappa[0]), side = 2, outer = TRUE,
      line = 2.1, cex = par("cex.lab") * par("cex"))
invisible(dev.off())

}

## Figure 2: modes
{

say("\n## Figure 2: samples, direction parameters, and modes on S^2_+")

# Two direction parameters and kappa = 10, 20, 50; 5000 draws per panel,
# coloured by within-panel density rank.
set.seed(20260705)
n <- 5000
kappas <- c(10, 20, 50)
mus <- list(unit_pos(a = c(0.20, 0.45, 0.87)),
            unit_pos(a = c(0.32, 0.74, 0.59)))
png("../img/modes.png", width = 9.4, height = 6.4, units = "in", res = 600,
    bg = "transparent")
par(mfrow = c(length(mus), length(kappas)), mar = c(1.2, 1.2, 2.0, 0.6),
    oma = c(0, 0, 0, 0))
pal <- viridisLite::viridis(256)
for (mu in mus) {

  for (k in kappas) {

    X <- r_vmfp(n = n, mu = mu, kappa = k)
    cv <- rank(log_f_vmfp(X = X, mu = mu, kappa = k), ties.method = "average")
    m <- mode_vmfp(mu = mu, kappa = k)
    pmat <- plot3D::scatter3D(X[, 1], X[, 2], X[, 3], colvar = cv, col = pal,
                              pch = 19, cex = 0.28, colkey = FALSE, bty = "b2",
                              theta = 135, phi = 35, ticktype = "detailed",
                              xlim = c(0, 1), ylim = c(0, 1), zlim = c(0, 1),
                              xlab = "", ylab = "", zlab = "")
    mu_xy <- plot3D::trans3D(mu[1], mu[2], mu[3], pmat = pmat)
    m_xy <- plot3D::trans3D(m[1], m[2], m[3], pmat = pmat)
    points(mu_xy$x, mu_xy$y, pch = 21, cex = 2.1, bg = "red", col = "black",
           lwd = 1.1, xpd = NA)
    points(m_xy$x, m_xy$y, pch = 21, cex = 2.1, bg = "blue", col = "black",
           lwd = 1.1, xpd = NA)
    title(main = bquote(kappa == .(k)), line = 0.3, cex.main = 1.25)

    # Activation condition of Proposition 6(iii), quoted in Section 3.3:
    # m_j > 0 if and only if xi_j^2 > lambda_m = ||grad g(m)||.
    xi <- k * mu
    say(paste0("kappa = %2d, mu = (%.3f, %.3f, %.3f): xi^2 = (%.1f, %.1f, ",
               "%.1f), lambda_m = %.1f"), k, mu[1], mu[2], mu[3], xi[1]^2,
        xi[2]^2, xi[3]^2, sqrt(sum((xi * tanh(xi * m))^2)))

  }

}
invisible(dev.off())

# Compressed copy: JPEG at 300 dpi and quality 85, on white as JPEG has no
# transparency (about 0.8 MB against 4.3 MB).
img <- magick::image_read(path = "../img/modes.png")
img <- magick::image_background(image = img, color = "white")
magick::image_write(image = magick::image_scale(image = img, geometry = "50%"),
                    path = "../img/modes.jpg", format = "jpeg", quality = 85)

}

## Figure 3: region shapes
{

say("\n## Figure 3: vMF+ region against the naive and renormalized vMF caps")

# Uniform draws on S^2_+ for the areas
set.seed(981587298)
U <- abs(matrix(rnorm(n = 3 * 5e6), ncol = 3))
U <- U / sqrt(rowSums(U^2))

# Three cases at level 1 - alpha = 0.5 on S^2_+, shown by the central
# projection x -> x / (x_1 + x_2 + x_3).
alpha <- 0.50
cases <- list(list(mu = unit_pos(a = c(0.99, 0.10, 0.07)), kappa = 5),
              list(mu = unit_pos(a = c(0.20, 0.45, 0.87)), kappa = 5),
              list(mu = unit_pos(a = c(1, 1, 1)), kappa = 2))

# Ternary grid over the triangle P with vertices e_1, e_2, e_3, lifted to the
# sphere by x = b / ||b|| for barycentric coordinates b.
P <- rbind(c(0, 0), c(1, 0), c(0.5, sqrt(3) / 2))
ng <- 900
uu <- seq(-0.02, 1.02, length.out = ng)
vv <- seq(-0.02, sqrt(3) / 2 + 0.02, length.out = ng)
B <- cbind(1, as.matrix(expand.grid(uu, vv))) %*% solve(cbind(1, P))
inside <- rowSums(B >= 0) == 3
S_x <- B[inside, ] / sqrt(rowSums(B[inside, ]^2))
col_R <- "#1b7837"
col_Cp <- "#762a83"
col_C <- "#b35806"

pdf("../img/app-region-shape.pdf", width = 12, height = 4.3)
par(mfrow = c(1, 3), mar = c(0.4, 0.4, 0.6, 0.4))
for (i in seq_along(cases)) {

  # Get parameters
  mu <- cases[[i]]$mu
  kappa <- cases[[i]]$kappa
  m <- mode_vmfp(mu = mu, kappa = kappa)

  # Thresholds, alpha-quantiles of the vMF+ log-density and of mu'x on the whole
  # sphere (naive cap) and on the vMF draws in the orthant (renormalized cap).
  Y <- rotasym::r_vMF(n = 1e6, mu = mu, kappa = kappa)
  X <- abs(Y)
  s_a <- quantile(log_f_vmfp(X = X, mu = mu, kappa = kappa), alpha)
  t_a <- quantile(Y %*% mu, alpha)
  t_ap <- quantile(Y[rowSums(Y >= 0) == 3, ] %*% mu, alpha)

  # Both fields on the grid, NA outside the triangle
  gf <- cf <- matrix(NA, ng, ng)
  gf[inside] <- log_f_vmfp(X = S_x, mu = mu, kappa = kappa)
  cf[inside] <- S_x %*% mu

  # Naive cap at the bottom, renormalized cap, vMF+ region on top
  plot.new()
  plot.window(xlim = c(-0.03, 1.03), ylim = c(-0.05, sqrt(3) / 2 + 0.05),
              asp = 1)
  polygon(P, col = "grey97", border = NA)
  image(uu, vv, ifelse(cf >= t_a, 1, NA), add = TRUE,
        col = adjustcolor(col_C, 0.20), useRaster = TRUE)
  image(uu, vv, ifelse(cf >= t_ap, 1, NA), add = TRUE,
        col = adjustcolor(col_Cp, 0.30), useRaster = TRUE)
  image(uu, vv, ifelse(gf >= s_a, 1, NA), add = TRUE,
        col = adjustcolor(col_R, 0.42), useRaster = TRUE)

  # Boundaries after all the fills, so that no fill covers them
  contour(uu, vv, cf, levels = c(t_a, t_ap), drawlabels = FALSE, add = TRUE,
          col = c(col_C, col_Cp), lwd = 2.8, lty = c(3, 2))
  contour(uu, vv, gf, levels = s_a, drawlabels = FALSE, add = TRUE,
          col = col_R, lwd = 2.8)
  polygon(P, border = "grey30", lwd = 1.6)

  # Direction mu (red), mode m (blue), and vertices e_1, e_2, e_3
  points(rbind(mu / sum(mu), m / sum(m)) %*% P, pch = 21,
         bg = c("red", "blue"), cex = 2.1, lwd = 1.2)
  text(P[1:2, ], labels = expression(e[1], e[2]), pos = 1, cex = 1.4)
  text(P[3, 1], P[3, 2], expression(e[3]), pos = 3, offset = 0.4, cex = 1.4)

  # Legend, in the first panel only
  if (i == 1) {

    legend("topleft", bty = "n", cex = 1.25,
           legend = c(expression(vMF["+"] ~ "region" ~ R[alpha]),
                      expression("Renormalized cap" ~ C[alpha]^"+"),
                      expression("Na\u00efve cap" ~ C[alpha]),
                      expression("Direction" ~ mu), expression("Mode" ~ m)),
           col = c(col_R, col_Cp, col_C, "black", "black"),
           lty = c(1, 2, 3, NA, NA), lwd = c(2.6, 2.6, 2.6, NA, NA),
           pch = c(NA, NA, NA, 21, 21), pt.bg = c(NA, NA, NA, "red", "blue"),
           pt.cex = 1.7)

  }

  # Quoted in Section 3.4: the vMF+ probabilities of the naive and
  # renormalized caps, and the area of the latter relative to that of R_alpha.
  a_R <- mean(log_f_vmfp(X = U, mu = mu, kappa = kappa) >= s_a)
  say(paste0("mu = (%.3f, %.3f, %.3f), kappa = %g: mode = (%.3f, %.3f, ",
             "%.3f); vMF+ probability of the naive cap = %.3f, of the ",
             "renormalized cap = %.3f; renormalized cap area = %.2f x ",
             "R_alpha"), mu[1], mu[2], mu[3], kappa, m[1], m[2], m[3],
      mean(X %*% mu >= t_a), mean(X %*% mu >= t_ap),
      mean(U %*% mu >= t_ap) / a_R)

}
invisible(dev.off())

}

## Application to text corpora and compositional data
{

say("\n## Text corpora and compositional data")

## Data
{

# Builds the 30 datasets of Section 5, 10 per group, in the order of Table 1.
# Each is a list with X, a dgCMatrix, and the group.
build_pool <- function() {

  ## Helper functions

  # Data object obj of a package, by default named as the item
  get_data <- function(pkg, item, obj = item) {

    data(list = item, package = pkg, envir = environment())
    get(obj)

  }

  # Term frequencies of documents, sparse: terms absent from every document
  # dropped, and documents normalized to unit length.
  docs <- function(X) {

    X <- X[, Matrix::colSums(X) > 0, drop = FALSE]
    list(X = Matrix::Diagonal(x = 1 / sqrt(Matrix::rowSums(X^2))) %*% X,
         group = "text")

  }

  # Compositions in one unit, normalized to unit sum, and square-rooted
  comp <- function(Y) {

    Y <- as.matrix(Y)
    Y <- Y[, colSums(Y) > 0, drop = FALSE]
    list(X = sqrt(Y / rowSums(Y)),
         group = if (any(Y == 0)) "compzero" else "compfree")

  }

  ## Text corpora

  # CLUTO term frequencies
  cluto <- "cluto.tar.gz"
  if (!file.exists(cluto)) {

    download.file(paste0("https://karypis.github.io/glaros/files/sw/cluto/",
                         "datasets.tar.gz"), destfile = cluto, mode = "wb")

  }
  stopifnot(unname(tools::md5sum(cluto)) == "66612cb11227254f9ef299bddbc10474")
  dir_cluto <- tempfile()
  pool <- list()
  for (corpus in c("new3", "reviews", "cranmed", "k1a", "la2", "ohscal", "tr31",
                   "tr11", "re0", "fbis")) {

    # Read the corpus as a sparse documents-by-terms matrix
    untar(cluto, files = paste0(corpus, ".mat"), exdir = dir_cluto)
    A <- slam::read_stm_CLUTO(file.path(dir_cluto, paste0(corpus, ".mat")))
    pool[[corpus]] <- docs(X = Matrix::sparseMatrix(i = A$i, j = A$j, x = A$v,
                                                    dims = c(A$nrow, A$ncol)))

  }

  ## Compositions with zeros

  # atlas1006: HITChip genus-level signals, the first sample (time 0) of each
  # subject.
  x <- get_data(pkg = "microbiome", item = "atlas1006")
  Y <- t(as(phyloseq::otu_table(x), "matrix"))
  pool$atlas1006 <- comp(Y = Y[phyloseq::sample_data(x)$time == 0, ])

  # enterotype: relative abundances of genera. The genera detected depend on the
  # sequencing technology, so many zeros are structural.
  x <- get_data(pkg = "phyloseq", item = "enterotype")
  pool$enterotype <- comp(Y = t(as(phyloseq::otu_table(x), "matrix")))

  # swapdiat: diatom percentages
  pool$swapdiat <- comp(Y = get_data(pkg = "analogue", item = "swapdiat"))

  # HIV: counts of genera
  pool$HIV <- comp(Y = get_data(pkg = "coda4microbiome", item = "HIV",
                                obj = "x_HIV"))

  # perthi02: counts of amino acids in proteins
  pool$perthi02 <- comp(Y = get_data(pkg = "ade4", item = "perthi02")$tab)

  # LPdata: concentrations (ug/g), zero below the detection limit
  pool$LPdata <- comp(Y = get_data(pkg = "zCompositions", item = "LPdata"))

  # piosphere: plant cover (%) of species around watering points; the covers
  # overlap, so closing them gives relative covers.
  pool$piosphere <- comp(Y = get_data(pkg = "ade4", item = "piosphere")$veg)

  # efsa: consumption of food groups (g/day) by country and population group
  x <- get_data(pkg = "robCompositions", item = "efsa")
  pool$efsa <- comp(Y = subset(x, select = -c(Country, Pop.Class)))

  # oribatid: counts of oribatid mite species in peat-moss cores
  pool$oribatid <- comp(Y = get_data(pkg = "ade4", item = "oribatid")$fau)

  # oliveoil: fatty acids of Italian olive oils (0.01%), zero below the
  # sensitivity of the instrument.
  x <- get_data(pkg = "pdfCluster", item = "oliveoil")
  pool$oliveoil <- comp(Y = subset(x, select = -c(macro.area, region)))

  ## Zero-free compositions

  # Hydrochem: concentrations (mg/L), H given as 10^(-pH) mol/L and converted
  # with its molar mass, 1.008 g/mol.
  x <- get_data(pkg = "compositions", item = "Hydrochem")
  x$H <- 1008 * x$H
  pool$Hydrochem <- comp(Y = subset(x, select = -c(Code, Site, Location, River,
                                                   Date)))

  # humus: elements of the Kola O horizon (mg/kg)
  x <- get_data(pkg = "mvoutlier", item = "humus")
  pool$humus <- comp(Y = subset(x, select = -c(ID, XCOO, YCOO, C, H, N, LOI, pH,
                                               Cond)))

  # gemas: elements of agricultural soils (mg/kg)
  x <- get_data(pkg = "robCompositions", item = "gemas")
  pool$gemas <- comp(Y = subset(x, select = -c(COUNTRY, longitude, latitude,
                                               Xcoord, Ycoord, MeanTemp,
                                               AnnPrec, soilclass, sand, silt,
                                               clay, LOI)))

  # Aar: major oxides (wt%) of the grain-size fractions of 11 sediment samples
  x <- get_data(pkg = "compositions", item = "Aar")
  pool$Aar <- comp(Y = subset(x, select = c(SiO2, TiO2, Al2O3, MnO, MgO, CaO,
                                            Na2O, K2O, P2O5, Fe2O3t)))

  # presid2002: votes for the candidates, first round, by department
  x <- get_data(pkg = "ade4", item = "presid2002")$tour1
  pool$presid2002 <- comp(Y = subset(x, select = -c(inscrits, abstentions,
                                                    votants, exprimes)))

  # trondelagO: aqua-regia elements of O-horizon soils (mg/kg), without the
  # Pb isotopes, which add up to Pb, and their ratios.
  x <- get_data(pkg = "robCompositions", item = "trondelagO")
  x <- subset(x, select = -c(X.Loc_ID, LITHO, longitude, latitude, E32wgs,
                             N32wgs, X.Medium, Alt_masl, LOI_480, pH, Pb204,
                             Pb206, Pb207, Pb208, X6_7Pb, X6_8Pb, X7_8Pb,
                             X6_4Pb, X7_4Pb, X8_4Pb))
  pool$trondelagO <- comp(Y = x)

  # Crohn: counts of genera plus one, which coda4microbiome adds to every count
  # (the selbal original has 28.8% zeros).
  pool$Crohn <- comp(Y = get_data(pkg = "coda4microbiome", item = "Crohn",
                                  obj = "x_Crohn"))

  # RBGlass2: trace elements of Romano-British glass (ppm)
  x <- get_data(pkg = "archdata", item = "RBGlass2")
  pool$RBGlass2 <- comp(Y = subset(x, select = -Site))

  # ageCatWorld: population (%) aged 0-14, 15-60, and over 60, by country
  x <- get_data(pkg = "robCompositions", item = "ageCatWorld")
  pool$ageCatWorld <- comp(Y = subset(x, select = -country))

  # mafragh: clay, silt, and sand fractions of marsh soils
  x <- get_data(pkg = "ade4", item = "mafragh")$env
  pool$mafragh <- comp(Y = x[, c("Clay", "Silt", "Sand")])

  ## Final check and order

  # Sparse storage and unit rows
  pool <- lapply(pool, function(d) {

    d$X <- Matrix::drop0(as(as(d$X, "CsparseMatrix"), "generalMatrix"))
    stopifnot(all(abs(sqrt(Matrix::rowSums(d$X^2)) - 1) < 1e-12))
    d

  })

  # Order of Table 1: by group, then by n
  grp <- match(sapply(pool, `[[`, "group"), c("text", "compzero", "compfree"))
  pool[order(grp, sapply(pool, function(d) nrow(d$X)))]

}

}

## Fits and regions
{

# Run fits and regions, if results file is not present. Takes about 1 h 50 min
# for the 30 datasets on 13 cores in an Apple M4 Pro. Each dataset is split at
# random into training (60%), calibration (20%), and test (20%) parts. The K
# scan and the fits use training data only. Coverage uses calibration and test
# parts.
results <- "folded-vmf-applications.rds"
if (file.exists(results)) {

  app <- readRDS(results)

} else {

  pool <- build_pool()
  datasets <- names(pool)
  app <- list()
  for (d in datasets) {

    # Start time of the dataset
    t0 <- Sys.time()

    # Split into training, calibration, and test parts
    X <- pool[[d]]$X
    n <- nrow(X)
    m <- round(0.2 * n)
    set.seed(20261016)
    id <- sample(n)
    i_tr <- sort(id[-seq_len(2 * m)])
    i_ca <- sort(id[seq_len(m)])
    i_te <- sort(id[m + seq_len(m)])
    X_tr <- X[i_tr, , drop = FALSE]
    X_ca <- X[i_ca, , drop = FALSE]
    X_te <- X[i_te, , drop = FALSE]

    # K scan and the MvMF, MrvMF, and MvMF+ fits at K_hat
    scan <- scan_K(X = X_tr, cores = cores)
    trio <- trio_K(X = X_tr, vmf_fit = scan$fit)

    # Plug-in and conformal regions: coverages and log-volumes
    regions <- regions_split(X_tr = X_tr, X_ca = X_ca, X_te = X_te, fit = trio,
                             cores = cores)

    # Group, descriptors, K_hat, log-likelihoods, diagnostics, and regions
    app[[d]] <- c(list(group = pool[[d]]$group),
                  list(n_all = n, zero_all = 1 - length(X@x) / (n * ncol(X)),
                       n = nrow(X_tr), p = ncol(X)),
                  scan[c("K_hat", "health")],
                  trio[c("loglik", "forced_gain")],
                  list(health_vmfp = trio$health, regions = regions))

    # Show progress
    cat(sprintf("Dataset %s done (%d of %d). Time taken = %s\n", d,
                match(d, datasets), length(datasets),
                format(round(Sys.time() - t0, 1))))

  }
  saveRDS(app, file = results)

}

}

## Tables 1 and 2, and the numbers of Section 5
{

# Table 1: n, p, and percentage of zero entries
say("Table 1 (n, p, zeros):")
say_df(x = data.frame(data = names(app), n = sapply(app, `[[`, "n_all"),
                      p = sapply(app, `[[`, "p"),
                      zeros = sprintf("%.1f",
                                      100 * sapply(app, `[[`, "zero_all"))))

# Table 2: K_hat; log-likelihood/n of MvMF+ minus those of MvMF and MrvMF, on
# the training part; and, for the plug-in and conformal regions, the
# coverage of C (= C+) and R and the log10 volume ratios of C and C+ to R.
one_digit <- function(x) sprintf("%.1f", round(x, 1) + 0)
delta_loglik <- function(a, m) (a$loglik[["vmfp"]] - a$loglik[[m]]) / a$n
region_cover_vol <- function(a, what, r, conformal = FALSE) {

  y <- if (conformal) a$regions$conformal else a$regions$plugin
  if (what == "cover") 100 * y$cover[[r]] else
    (y$lv[[r]] - y$lv[["R"]]) / log(10)

}
say(paste0("\nTable 2 (K, log-likelihood/n gain of MvMF+, plug-in and ",
           "conformal regions):"))
tab <- data.frame(data = names(app), K = sapply(app, `[[`, "K_hat"),
d_MvMF = one_digit(x = sapply(app, delta_loglik, m = "vmf")),
d_MrvMF = one_digit(x = sapply(app, delta_loglik, m = "rvmf")))
for (cf in c(FALSE, TRUE)) {

  pre <- if (cf) "cf_" else ""
  for (r in c("C", "R")) {

    tab[[paste0(pre, r)]] <- one_digit(x = sapply(app, region_cover_vol,
                                                  what = "cover", r = r,
                                                  conformal = cf))

  }
  for (r in c("C", "Cp")) {

    tab[[paste0(pre, r, "_R")]] <- one_digit(x = sapply(app, region_cover_vol,
                                                        what = "log10_ratio",
                                                        r = r, conformal = cf))

  }

}
say_df(x = tab)

# Dataset groups
grp <- sapply(app, `[[`, "group")
tx <- grp == "text"
zeros <- grp == "compzero"
zero_free <- grp == "compfree"

# Section 5.3: safeguards, movMF over every K scanned and MvMF+ at K_hat
H <- do.call(cbind, lapply(app, `[[`, "health"))
V <- sapply(app, function(a) a$health_vmfp[["capped_vmfp"]])
say(paste0("movMF runs excluded: %d of %d, %d ending in an error and %d with ",
           "kappa > 1e10; fits capped at 2000 iterations: movMF %d, MvMF+ %d"),
    sum(H[c("failed_vmf", "degenerate_vmf"), ]), sum(H["runs", ]),
    sum(H["failed_vmf", ]), sum(H["degenerate_vmf", ]),
    sum(H["capped_vmf", ]), sum(V))

# Section 5.4: log-likelihood/n gains of MvMF+ over MrvMF at K_hat and the
# datasets where the bound forces the gain.
gain_rvmf <- sapply(app, delta_loglik, m = "rvmf")
forced_gain <- sapply(app, `[[`, "forced_gain")
say("\nText: gain over MrvMF %.0f to %.0f; forced_gain: %d of %d",
    min(gain_rvmf[tx]), max(gain_rvmf[tx]), sum(forced_gain[tx]), sum(tx))
say("Zeros: gain %.1f to %.1f; forced_gain: %d of %d", min(gain_rvmf[zeros]),
    max(gain_rvmf[zeros]), sum(forced_gain[zeros]), sum(zeros))
say("Zero-free: gains %s; forced_gain anywhere: %s",
    paste(names(app)[zero_free], one_digit(x = gain_rvmf[zero_free]),
          collapse = ", "),
    any(forced_gain[zero_free]))

# Section 5.4: the 2^p orthants of S^{p-1} for the text corpora, in log10
p_text <- sapply(app[tx], `[[`, "p")
say("Text: p log10(2) from %.0f to %.0f", min(p_text) * log10(2),
    max(p_text) * log10(2))

# Section 5.4: regions. For each calibration, the coverage of C (= C+) and R,
# and the log10 volume ratios to R; |log10 C+/R| < 0.1 is about equal.
for (conformal in c(FALSE, TRUE)) {

  label <- if (conformal) "Conformal" else "Plug-in"
  coverage <- sapply(c(C = "C", R = "R"), function(r) {
    sapply(app, region_cover_vol, what = "cover", r = r, conformal = conformal)
  })
  log10_ratio <- sapply(c(C = "C", Cp = "Cp"), function(r) {
    sapply(app, region_cover_vol, what = "log10_ratio", r = r,
           conformal = conformal)
  })
  near <- !tx & abs(log10_ratio[, "Cp"]) < 0.1
  say(paste0("%s: coverage %.1f to %.1f in the text corpora, %.1f to %.1f in ",
             "the compositions; text log10 C/R %.1f to %.1f, C+/R %.1f to ",
             "%.1f; compositions largest C+/R %.1f (%s), |C+/R| < 0.1 in %s ",
             "(within %.0f%%)"),
      label,
      min(coverage[tx, ]), max(coverage[tx, ]),
      min(coverage[!tx, ]), max(coverage[!tx, ]),
      min(log10_ratio[tx, "C"]), max(log10_ratio[tx, "C"]),
      min(log10_ratio[tx, "Cp"]), max(log10_ratio[tx, "Cp"]),
      max(log10_ratio[!tx, "Cp"]),
      names(app)[!tx][which.max(log10_ratio[!tx, "Cp"])],
      paste(names(app)[near], sprintf("%.2f", log10_ratio[near, "Cp"]),
            collapse = ", "),
      100 * (10^max(abs(log10_ratio[near, "Cp"])) - 1))

  # A volume ratio rho^(p - 1) in S^{p-1} is a factor rho in linear size
  rho <- 10^(log10_ratio[tx, "Cp"] / (p_text - 1))
  say("%s: text C+/R per dimension, rho = %.4f to %.4f", label, min(rho),
      max(rho))

}

}

}
