###############################################################################
# sidx_simulation.R -- simulation with KNOWN net merit H = w'g (answer to review point M2).
# Genotypes: real indica rice markers (G = Z Z'/p, n = 327). Four traits: main trait (tau) + three auxiliary traits.
# Genetic values vec(U) ~ N(0, Sigma_g (x) G); residuals vec(E) ~ N(0, Sigma_e (x) I); Y = U + E.
# Scenarios: genetic correlation main-auxiliary r_g in {0.2, 0.7} x residual correlation main-auxiliary r_e in {0.0, 0.6}.
# h2 = 0.4 (main) and 0.6 (auxiliary); aux-aux correlations: genetic 0.3, residual 0.2.
# Economic weights: 0.5 (main) and 1/6 for each auxiliary trait. CV2: main trait of 20% of the lines masked.
# Methods: SH (parameters from BGLR), SH-true (true parameters), ST-G, GSH-G (all models fitted with BGLR).
# Criteria (testing lines): accuracy = cor(H_hat, H_true); predictive ability = cor(H_hat, w'y);
#   the same two for the main trait: cor(g_hat_tau, g_tau) and cor(g_hat_tau, y_tau).
# Usage: Rscript sidx_simulation.R <data_dir> <out_dir> <n_reps> <nIter> <burnIn>
###############################################################################
suppressPackageStartupMessages(library(BGLR))
args <- commandArgs(TRUE); DD <- args[1]; OUT <- file.path(args[2], "simulation"); dir.create(OUT, FALSE, TRUE)
NR <- as.integer(args[3]); nIter <- as.integer(args[4]); burnIn <- as.integer(args[5]); tmp <- tempfile("sim_"); dir.create(tmp)
X <- as.matrix(read.csv(file.path(DD, "X_Indica.csv"), row.names = 1, check.names = FALSE))
f <- colMeans(X) / 2; keep <- pmin(f, 1 - f) > 0.05 & apply(X, 2, sd) > 0
Z <- scale(X[, keep], TRUE, TRUE); G <- tcrossprod(Z) / ncol(Z); n <- nrow(G)
EG <- eigen(G, symmetric = TRUE); LG <- EG$vectors %*% diag(sqrt(pmax(EG$values, 0)))
t <- 4; tg <- 1; aux <- 2:4; w <- c(0.5, 1/6, 1/6, 1/6); h2 <- c(0.4, 0.6, 0.6, 0.6)
cormat <- function(r_main, r_aux) { C <- matrix(r_aux, t, t); C[1, ] <- C[, 1] <- r_main; diag(C) <- 1; C }
sh_index <- function(Yna, rows, Sg, Se) {
  P <- Sg + Se; out <- matrix(NA, nrow(Yna), ncol(Yna))
  for (i in rows) { o <- which(!is.na(Yna[i, ])); out[i, ] <- Sg[, o, drop = FALSE] %*% solve(P[o, o], Yna[i, o]) }
  out
}
scen <- expand.grid(r_g = c(0.2, 0.7), r_e = c(0.0, 0.6))
res <- list()
for (s in seq_len(nrow(scen))) {
  Cg <- cormat(scen$r_g[s], 0.3); Ce <- cormat(scen$r_e[s], 0.2)
  Sg <- diag(sqrt(h2)) %*% Cg %*% diag(sqrt(h2)); Se <- diag(sqrt(1 - h2)) %*% Ce %*% diag(sqrt(1 - h2))
  stopifnot(min(eigen(Sg)$values) > 0, min(eigen(Se)$values) > 0)
  for (rep in seq_len(NR)) {
    f1 <- file.path(OUT, sprintf("sim_s%d_rep%02d.csv", s, rep)); if (file.exists(f1)) { res[[length(res) + 1]] <- read.csv(f1); next }
    set.seed(7000 + 100 * s + rep)
    U <- LG %*% matrix(rnorm(n * t), n, t) %*% chol(Sg)                    # vec(U) ~ N(0, Sg (x) G); chol(Sg) = R with R'R = Sg
    E <- matrix(rnorm(n * t), n, t) %*% chol(Se); Y <- U + E
    te <- sort(sample(n, round(0.2 * n))); Yna <- Y; Yna[te, tg] <- NA
    fm <- Multitrait(y = Yna, ETA = list(list(K = G, model = "RKHS")), nIter = nIter, burnIn = burnIn, verbose = FALSE,
                     saveAt = file.path(tmp, "mt_"))
    Sg_hat <- fm$ETA[[1]]$Cov$Omega; Se_hat <- fm$resCov$R
    Uh <- list("GSH-G" = fm$ETA[[1]]$u,
               "ST-G" = sapply(1:t, function(k) BGLR(y = Yna[, k], ETA = list(list(V = EG$vectors, d = EG$values, model = "RKHS")),
                                                     nIter = nIter, burnIn = burnIn, verbose = FALSE, saveAt = file.path(tmp, "st_"))$ETA[[1]]$u),
               "SH" = sh_index(Yna, te, Sg_hat, Se_hat), "SH-true" = sh_index(Yna, te, Sg, Se))
    Htrue <- as.numeric(U[te, ] %*% w); Hobs <- as.numeric(Y[te, ] %*% w)
    out <- do.call(rbind, lapply(names(Uh), function(m) { Hh <- as.numeric(Uh[[m]][te, ] %*% w); gh <- Uh[[m]][te, tg]
      data.frame(scenario = s, r_g = scen$r_g[s], r_e = scen$r_e[s], rep = rep, method = m,
                 acc_H = cor(Hh, Htrue), pa_H = cor(Hh, Hobs), acc_main = cor(gh, U[te, tg]), pa_main = cor(gh, Y[te, tg])) }))
    write.csv(out, f1, row.names = FALSE); res[[length(res) + 1]] <- out
    message(sprintf("scenario %d rep %d", s, rep))
  }
}
d <- do.call(rbind, res)
agg <- aggregate(cbind(acc_H, pa_H, acc_main, pa_main) ~ scenario + r_g + r_e + method, d, mean)
agg_sd <- aggregate(cbind(acc_H, pa_H, acc_main, pa_main) ~ scenario + r_g + r_e + method, d, function(x) sd(x) / sqrt(length(x)))
names(agg_sd)[5:8] <- paste0(names(agg_sd)[5:8], "_se"); agg <- merge(agg, agg_sd)
write.csv(agg[order(agg$scenario, agg$method), ], file.path(OUT, "simulation_summary.csv"), row.names = FALSE)
print(agg[order(agg$scenario, agg$method), ], digits = 3)
