## Smoke test do svar_engine: simular VAR(2) com choque de petroleo,
## estimar SVAR identificado por restricoes de sinal e checar IRFs/FEVD.

source("/home/user/petroleo/svar_engine.R")
set.seed(42)

# Variaveis (na ordem):
#   1: brent          (preco do Brent, log-diff)
#   2: prod           (producao global de petroleo, log-diff)
#   3: atividade      (atividade brasileira, log-diff)
#   4: cambio         (BRL/USD, log-diff)
#   5: ipca_combust   (variacao % m/m do IPCA combustiveis)
#   6: ipca_servicos  (variacao % m/m do IPCA servicos) - VARIAVEL DE INTERESSE

k <- 6
varnames <- c("brent","prod","atividade","cambio","ipca_combust","ipca_servicos")

# Definir B0 verdadeiro: efeito de impacto dos 6 choques sobre as 6 variaveis.
# Coluna 1 = choque de petroleo, com:  brent ↑, prod ↓, ipca_combust ↑.
# Demais choques irrelevantes para o teste, geramos com Choleski + ruido.
B0_true <- matrix(0, k, k)
B0_true[, 1] <- c(0.060, -0.008, -0.003, 0.012, 0.80, 0.05)   # choque de petroleo
B0_true[, 2] <- c(0.000,  0.006, 0.002, -0.001, 0.02, 0.00)   # choque de oferta
B0_true[, 3] <- c(0.005,  0.000, 0.008,  0.000, 0.05, 0.02)   # demanda global
B0_true[, 4] <- c(0.000,  0.000, 0.001, 0.030, 0.10, 0.05)    # cambial
B0_true[, 5] <- c(0.000,  0.000, 0.000,  0.000, 0.40, 0.05)   # combustiveis idiossincratico
B0_true[, 6] <- c(0.000,  0.000, 0.001,  0.001, 0.00, 0.30)   # servicos idiossincratico

# Coeficientes VAR(2): persistencia moderada
A1 <- diag(c(0.40, 0.55, 0.60, 0.45, 0.50, 0.65))
A1[5, 1] <- 0.15   # combust depende de Brent defasado (pass-through)
A1[6, 5] <- 0.05   # servicos: efeito pequeno de combust defasado
A1[6, 4] <- 0.08   # servicos sente cambio defasado
A1[6, 3] <- 0.08   # servicos depende de atividade
A2 <- diag(c(0.10, 0.15, 0.10, 0.05, 0.10, 0.10))
A2[5, 1] <- 0.10
A2[6, 5] <- 0.04

# Simular Tn observacoes
Tn <- 240
Y <- matrix(0, Tn + 50, k)            # burn-in 50
for (t in 3:nrow(Y)) {
  eps <- rnorm(k)
  u   <- B0_true %*% eps              # inovacoes reduzidas
  Y[t, ] <- A1 %*% Y[t - 1, ] + A2 %*% Y[t - 2, ] + u
}
Y <- Y[(50 + 1):nrow(Y), ]
colnames(Y) <- varnames

cat("=== Lag selection (AIC/BIC/HQC) ===\n")
ls_out <- lag_selection(Y, pmax = 8)
print(round(ls_out$table, 3))
cat("Best:\n"); print(ls_out$best)

p_use <- as.integer(ls_out$best["BIC"])
cat("Usando p =", p_use, "(criterio BIC)\n")

cat("\n=== Estimando VAR reduzido ===\n")
v <- estimate_var(Y, p = p_use)
cat("dim(B) =", dim(v$B), "  T =", v$Tr, "\n")
cat("Sigma residual:\n"); print(round(v$Sigma, 5))

cat("\n=== Identificacao por restricoes de sinal ===\n")
# Restricoes durante h = 0..5 (6 meses):
#  choque de petroleo (col 1): brent ↑, prod ↓, ipca_combust ↑
#  ipca_servicos: livre
sign_mat <- matrix(0, k, k)
sign_mat[1, 1] <-  1     # brent up
sign_mat[2, 1] <- -1     # prod down
sign_mat[5, 1] <-  1     # ipca_combust up
# Demais colunas: sem restricoes (sao identificadas como "outras" mas nao usamos)

sr <- sign_restrict(v, sign_mat, h_rest = 6, n_draws = 300, seed = 7)
cat("Aceitas:", sr$accepted, " | tentativas:", sr$tries, "\n")

cat("\n=== IRFs (mediana) por 24 meses ===\n")
ib <- irf_bands(v, sr$draws, h = 24)
med <- ib$quantiles[2, , , 1]   # respostas ao choque 1 (petroleo)
colnames(med) <- varnames
cat("Resposta de cada variavel ao choque de petroleo (mediana, h=0,3,6,12,24):\n")
print(round(med[c(1, 4, 7, 13, 25), ], 4))

cat("\n=== FEVD do IPCA servicos (% explicado pelo choque de petroleo) ===\n")
fb <- fevd_bands(v, sr$draws, h = 24)
fev_serv_oil <- fb$quantiles[2, , 6, 1]
cat("Mediana em h=1,6,12,24:", round(fev_serv_oil[c(2, 7, 13, 25)] * 100, 2), "%\n")

cat("\nTeste OK.\n")
