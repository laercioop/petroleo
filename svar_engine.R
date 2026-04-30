## SVAR engine - implementacao em R base
## VAR reduzido por OLS, selecao AIC/BIC, identificacao por restricoes de sinal
## (rejection sampler de Rubio-Ramirez/Waggoner/Zha), IRFs e FEVD.

# -----------------------------------------------------------------------------
# 1. VAR(p) reduzido via OLS (equacao por equacao)
# -----------------------------------------------------------------------------
estimate_var <- function(Y, p, const = TRUE) {
  Y <- as.matrix(Y)
  Tn <- nrow(Y); k <- ncol(Y)
  if (p < 1) stop("p>=1")
  Yt  <- Y[(p+1):Tn, , drop = FALSE]
  X <- NULL
  for (i in 1:p) X <- cbind(X, Y[(p+1-i):(Tn-i), , drop = FALSE])
  if (const) X <- cbind(1, X)
  B <- solve(crossprod(X)) %*% crossprod(X, Yt)        # (1+kp) x k
  U <- Yt - X %*% B                                     # residuos
  Tr <- nrow(Yt)
  Sigma <- crossprod(U) / (Tr - ncol(X))                # covariancia residual (n-k ajuste)
  list(B = B, U = U, Sigma = Sigma, p = p, k = k, const = const,
       Tr = Tr, X = X, Y = Yt)
}

# Forma companheira para iterar IRFs
companion_form <- function(varfit) {
  k <- varfit$k; p <- varfit$p
  B <- varfit$B
  off <- if (varfit$const) 1 else 0
  A <- t(B[(off + 1):nrow(B), , drop = FALSE])           # k x kp
  if (p == 1) {
    F <- A
  } else {
    F <- rbind(A, cbind(diag(k * (p - 1)), matrix(0, k * (p - 1), k)))
  }
  F
}

# -----------------------------------------------------------------------------
# 2. Selecao de defasagens (AIC, BIC, HQC) - mesma amostra para comparabilidade
# -----------------------------------------------------------------------------
lag_selection <- function(Y, pmax = 12, const = TRUE) {
  Y <- as.matrix(Y); Tn <- nrow(Y); k <- ncol(Y)
  Tcom <- Tn - pmax                                      # amostra comum
  out <- matrix(NA_real_, pmax, 3,
                dimnames = list(paste0("p=", 1:pmax), c("AIC", "BIC", "HQC")))
  for (p in 1:pmax) {
    Yt <- Y[(pmax + 1):Tn, , drop = FALSE]
    X <- NULL
    for (i in 1:p) X <- cbind(X, Y[(pmax + 1 - i):(Tn - i), , drop = FALSE])
    if (const) X <- cbind(1, X)
    B <- solve(crossprod(X)) %*% crossprod(X, Yt)
    U <- Yt - X %*% B
    Sigma_ml <- crossprod(U) / Tcom                      # MLE (sem ajuste)
    ld <- determinant(Sigma_ml, logarithm = TRUE)$modulus
    npar <- k * p + ifelse(const, 1, 0)
    out[p, "AIC"] <- ld + 2 * npar * k / Tcom
    out[p, "BIC"] <- ld + log(Tcom) * npar * k / Tcom
    out[p, "HQC"] <- ld + 2 * log(log(Tcom)) * npar * k / Tcom
  }
  best <- apply(out, 2, which.min)
  list(table = out, best = best)
}

# -----------------------------------------------------------------------------
# 3. IRFs reduzidas (resposta a choques nas inovacoes ortogonalizadas via B0)
#    Dado choque estrutural em coluna j de B0 (impacto), iterar via companheira.
# -----------------------------------------------------------------------------
irf_struct <- function(varfit, B0, h = 24) {
  k <- varfit$k; p <- varfit$p
  Fc <- companion_form(varfit)                            # kp x kp
  J <- cbind(diag(k), matrix(0, k, k * (p - 1)))          # k x kp
  IRF <- array(0, dim = c(h + 1, k, k))                   # [horizonte, var, choque]
  Phi_h <- diag(k * p)
  for (hh in 0:h) {
    Psi <- J %*% Phi_h %*% t(J)                           # multiplicador de Wold em h
    IRF[hh + 1, , ] <- Psi %*% B0
    Phi_h <- Phi_h %*% Fc
  }
  IRF
}

# -----------------------------------------------------------------------------
# 4. FEVD - decomposicao da variancia do erro de previsao por choque estrutural
# -----------------------------------------------------------------------------
fevd_struct <- function(IRF) {
  hp <- dim(IRF)[1]; k <- dim(IRF)[2]
  FE <- array(0, c(hp, k, k))                              # [h+1, var, choque]
  cum_var <- matrix(0, k, k)
  for (hh in 1:hp) {
    cum_var <- cum_var + (IRF[hh, , ])^2
    tot <- rowSums(cum_var)
    FE[hh, , ] <- cum_var / tot                            # cada linha = var; coluna = choque
  }
  FE
}

# -----------------------------------------------------------------------------
# 5. Identificacao por restricoes de sinal (rejection sampler)
#    sign_mat: matriz k x k, com entradas +1 / -1 / 0 (sem restricao)
#              linha = variavel, coluna = choque
#    h_rest:   numero de horizontes (incluindo 0) em que a restricao deve valer
#    n_draws:  numero de rotacoes ortogonais aceitas a coletar
# -----------------------------------------------------------------------------
sign_restrict <- function(varfit, sign_mat, h_rest = 6, n_draws = 500,
                          max_tries = 200000, seed = 1) {
  set.seed(seed)
  k <- varfit$k
  P <- t(chol(varfit$Sigma))                               # Σ = P P'
  acc <- list(); tries <- 0
  while (length(acc) < n_draws && tries < max_tries) {
    tries <- tries + 1
    Z <- matrix(rnorm(k * k), k, k)
    qr_z <- qr(Z); Q <- qr.Q(qr_z); R <- qr.R(qr_z)
    sgn <- sign(diag(R)); sgn[sgn == 0] <- 1
    Q <- Q %*% diag(sgn)                                   # uniforme em O(k)
    B0 <- P %*% Q
    # Para cada coluna j (choque candidato), checar se ha sinal viavel no impacto
    B0_cand <- matrix(NA_real_, k, k); ok_cols <- TRUE
    for (j in 1:k) {
      restr_col <- sign_mat[, j]
      if (all(restr_col == 0)) {
        B0_cand[, j] <- B0[, j]
      } else {
        idx_r <- which(restr_col != 0)
        v <- B0[, j]
        if (all(sign(v[idx_r]) == restr_col[idx_r])) {
          B0_cand[, j] <- v
        } else if (all(sign(-v[idx_r]) == restr_col[idx_r])) {
          B0_cand[, j] <- -v
        } else { ok_cols <- FALSE; break }
      }
    }
    if (!ok_cols) next
    # Verificar restricoes de sinal em h=0..h_rest-1 para todos os choques
    IRF <- irf_struct(varfit, B0_cand, h = h_rest - 1)
    ok <- TRUE
    for (j in 1:k) {
      restr <- sign_mat[, j]
      if (all(restr == 0)) next
      idx <- which(restr != 0)
      for (hh in 1:h_rest) {
        sg <- sign(IRF[hh, idx, j])
        if (any(sg * restr[idx] <= 0)) { ok <- FALSE; break }
      }
      if (!ok) break
    }
    if (ok) acc[[length(acc) + 1]] <- B0_cand
  }
  list(draws = acc, tries = tries, accepted = length(acc))
}

# -----------------------------------------------------------------------------
# 6. Quantis de IRF / FEVD a partir de uma lista de B0 aceitos
# -----------------------------------------------------------------------------
irf_bands <- function(varfit, draws, h = 24, probs = c(0.16, 0.5, 0.84)) {
  n <- length(draws); k <- varfit$k
  arr <- array(NA_real_, c(n, h + 1, k, k))
  for (i in 1:n) arr[i, , , ] <- irf_struct(varfit, draws[[i]], h = h)
  out <- array(NA_real_, c(length(probs), h + 1, k, k))
  for (a in 1:(h + 1)) for (b in 1:k) for (c in 1:k)
    out[, a, b, c] <- quantile(arr[, a, b, c], probs = probs, na.rm = TRUE)
  list(quantiles = out, probs = probs, all = arr)
}

fevd_bands <- function(varfit, draws, h = 24, probs = c(0.16, 0.5, 0.84)) {
  n <- length(draws); k <- varfit$k
  arr <- array(NA_real_, c(n, h + 1, k, k))
  for (i in 1:n) {
    irf_i <- irf_struct(varfit, draws[[i]], h = h)
    arr[i, , , ] <- fevd_struct(irf_i)
  }
  out <- array(NA_real_, c(length(probs), h + 1, k, k))
  for (a in 1:(h + 1)) for (b in 1:k) for (c in 1:k)
    out[, a, b, c] <- quantile(arr[, a, b, c], probs = probs, na.rm = TRUE)
  list(quantiles = out, probs = probs, all = arr)
}
