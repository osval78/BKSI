###############################################################################
# sidx_example.R (version 6) -- step-by-step computation of the classical Smith-Hazel (SH) index and of the
# Generalized Smith-Hazel (GSH) index for one candidate (indica rice, partition 1, GY masked in the testing lines).
# Changes with respect to version 5:
#   * traits standardized with the training lines (as in sidx_index.R)
#   * the candidate is chosen by an objective rule: the testing line whose genomic self-relationship G_ii is closest
#     to 1 among the testing lines with all own records within 3 SD (a "typical" candidate)
#   * the GSH index is computed with Gamma = G and with Gamma = K, and its contributions are split by the sign of Gamma_ij
#   * reliabilities use Var(H_i) = w' Sigma_g w Gamma_ii (Equation for the reliability, corrected)
# Usage: Rscript sidx_example.R <data_dir> <out_dir>
###############################################################################
suppressPackageStartupMessages(library(BGLR))
args <- commandArgs(TRUE); DD <- args[1]; OUT <- file.path(args[2], "example"); dir.create(OUT, FALSE, TRUE)
X <- as.matrix(read.csv(file.path(DD, "X_Indica.csv"), row.names = 1, check.names = FALSE))
Yd <- read.csv(file.path(DD, "Y_Indica.csv"), check.names = FALSE); Y0 <- as.matrix(Yd[, c("GC", "GY", "PH", "PHR")])
lines <- as.character(Yd$Line); n <- nrow(Y0); t <- ncol(Y0); tg <- 2; aux <- c(1, 3, 4); tn <- colnames(Y0)
w <- c(GC = -0.15, GY = 0.50, PH = -0.15, PHR = 0.20)
f <- colMeans(X) / 2; keep <- pmin(f, 1 - f) > 0.05 & apply(X, 2, sd) > 0
Z <- scale(X[, keep], TRUE, TRUE); G <- tcrossprod(Z) / ncol(Z)
D2 <- outer(diag(G), diag(G), "+") - 2 * G; K <- exp(-D2 / median(D2[upper.tri(D2)]))
set.seed(3001); te <- sort(sample(n, round(0.2 * n))); trn <- setdiff(seq_len(n), te)
Y <- sweep(sweep(Y0, 2, colMeans(Y0[trn, ])), 2, apply(Y0[trn, ], 2, sd), "/")      # standardized with training lines
Yna <- Y; Yna[te, tg] <- NA
# objective choice of the candidate
okc <- te[apply(abs(Y[te, ]) < 3, 1, all)]; i <- okc[which.min(abs(diag(G)[okc] - 1))]
writeLines(c(sprintf("candidate=%s", lines[i]), sprintf("rule=testing line with G_ii closest to 1 among %d testing lines with |z|<3 in all traits", length(okc)),
             sprintf("G_ii=%.4f", G[i, i]), sprintf("median diag(G)=%.4f; range diag(G)=%.3f-%.3f", median(diag(G)), min(diag(G)), max(diag(G))),
             sprintf("number of |z|>4 records in indica: %d", sum(abs(Y) > 4))), file.path(OUT, "candidate_choice.txt"))
run_case <- function(Gam, lab) {
  fm <- Multitrait(y = Yna, ETA = list(list(K = Gam, model = "RKHS")), nIter = 12000, burnIn = 2000, verbose = FALSE,
                   saveAt = file.path(tempdir(), paste0("ex_", lab, "_")))
  Sg <- fm$ETA[[1]]$Cov$Omega; Se <- fm$resCov$R; P <- Sg + Se; mu <- colMeans(Yna, na.rm = TRUE)
  dimnames(Sg) <- dimnames(Se) <- dimnames(P) <- list(tn, tn)
  sgH <- as.numeric(Sg %*% w); names(sgH) <- tn; s2H <- as.numeric(t(w) %*% Sg %*% w); Gii <- Gam[i, i]
  # ---- classical index (own auxiliary records)
  Pinv <- solve(P[aux, aux]); b <- as.numeric(Pinv %*% sgH[aux]); dev <- Yna[i, aux] - mu[aux]; contrib <- b * dev
  bSgw <- sum(b * sgH[aux]); rel_SH_textbook <- bSgw / s2H
  relSH <- function(g) as.numeric((bSgw * g)^2 / ((t(b) %*% (Sg[aux, aux] * g + Se[aux, aux]) %*% b) * s2H * g))
  # ---- GSH index, closed form
  obs <- which(!is.na(as.vector(Yna))); V <- (kronecker(Sg, Gam) + kronecker(Se, diag(n)))[obs, obs]
  yo <- as.vector(Yna)[obs] - rep(mu, each = n)[obs]; CH <- kronecker(matrix(sgH, ncol = 1), Gam)[obs, , drop = FALSE]
  Vi <- chol2inv(chol(V)); ci <- CH[, i]; beta <- as.numeric(Vi %*% ci); contr <- beta * yo; H_GSH <- sum(contr)
  rel_GSH <- sum(ci * beta) / (s2H * Gii)
  cl <- rep(seq_len(n), t)[obs]; ct <- rep(seq_len(t), each = n)[obs]
  grp <- ifelse(cl == i, "Own records (auxiliary traits)", ifelse(ct == tg, "Other lines: GY records", "Other lines: auxiliary records"))
  sgn <- ifelse(cl == i, "own", ifelse(Gam[i, cl] > 0, "Gamma_ij > 0", "Gamma_ij <= 0"))
  G1 <- aggregate(contr ~ grp, FUN = sum); G1$n_records <- as.numeric(table(grp)[G1$grp])
  G2 <- aggregate(contr ~ sgn, FUN = sum); G2$n_records <- as.numeric(table(sgn)[G2$sgn])
  byl <- aggregate(contr ~ cl, FUN = sum); byl <- byl[byl$cl != i, ]; byl$Gamma_ij <- Gam[i, byl$cl]
  byl <- byl[order(-abs(byl$contr)), ]
  top <- head(byl, 6); top$line <- lines[top$cl]; top$y_GY <- Yna[top$cl, tg]
  top$beta_GY <- sapply(top$cl, function(j) { k <- which(cl == j & ct == tg); if (length(k)) beta[k] else NA })
  top20_neg <- sum(head(byl, 20)$Gamma_ij <= 0)
  j <- top$cl[1]; recs <- which(cl %in% c(i, j))
  det <- data.frame(line = lines[cl[recs]], role = ifelse(cl[recs] == i, "candidate", "other line"), trait = tn[ct[recs]], y_dev = yo[recs],
                    Sgw_k = sgH[ct[recs]], Gamma_ij = Gam[i, cl[recs]], c_entry = ci[recs], beta = beta[recs], contribution = contr[recs])
  k1 <- which(cl == i & ct == aux[1]); k2 <- which(cl == j & ct == tg)
  vb <- data.frame(entry = c("V[own GC, own GC]", "V[own GC, other GY]", "V[other GY, other GY]"), value = c(V[k1, k1], V[k1, k2], V[k2, k2]),
                   Gamma_jj = Gam[j, j], Gamma_ij = Gam[i, j])
  # ---- all testing lines: agreement with BGLR, predictive ability, reliabilities
  H_closed <- as.numeric(crossprod(CH[, te], Vi %*% yo)); H_bglr <- as.numeric(fm$ETA[[1]]$u[te, ] %*% w)
  SHall <- sapply(te, function(jj) sum(b * (Yna[jj, aux] - mu[aux]))); Hobs <- as.numeric(Y[te, ] %*% w)
  relGSH_te <- colSums(CH[, te] * (Vi %*% CH[, te])) / (s2H * diag(Gam)[te]); relSH_te <- sapply(diag(Gam)[te], relSH)
  summ <- data.frame(kernel = lab, candidate = lines[i], Gamma_ii = Gii, H_obs = Hobs[match(i, te)], y_GY_masked = Y[i, tg],
                     I_SH = sum(contrib), H_GSH_closed = H_GSH, H_GSH_BGLR = H_bglr[match(i, te)], rel_SH_textbook = rel_SH_textbook,
                     rel_SH = relSH(Gii), rel_GSH = rel_GSH, n_obs = length(obs), cor_closed_BGLR = cor(H_closed, H_bglr),
                     r_SH_test = cor(SHall, Hobs), r_GSH_closed_test = cor(H_closed, Hobs), r_GSH_BGLR_test = cor(H_bglr, Hobs),
                     mean_rel_GSH_test = mean(relGSH_te), mean_rel_SH_test = mean(relSH_te), top20_nonpositive_Gamma = top20_neg,
                     rank_SH = rank(-SHall)[match(i, te)], rank_GSH = rank(-H_closed)[match(i, te)], rank_obs = rank(-Hobs)[match(i, te)])
  pre <- function(x) file.path(OUT, sprintf("%s_%s.csv", lab, x))
  write.csv(round(cbind(Sg = Sg, Se = Se, P = P), 4), pre("step1_parameters"))
  write.csv(data.frame(trait = tn, w = w, Sg_w = sgH, h2 = diag(Sg) / diag(P)), pre("step1_weights"), row.names = FALSE)
  write.csv(data.frame(trait = tn[aux], Paa_row = apply(round(P[aux, aux], 4), 1, paste, collapse = "; "), sigma_gH = sgH[aux], b = b,
                       y_i = Yna[i, aux], mu = mu[aux], deviation = dev, contribution = contrib), pre("step2_SH"), row.names = FALSE)
  write.csv(round(cbind(P_aa = P[aux, aux], P_aa_inv = Pinv), 5), pre("step2_Paa_inverse"))
  write.csv(G1, pre("step3_groups"), row.names = FALSE); write.csv(G2, pre("step3_by_sign"), row.names = FALSE)
  write.csv(top[, c("line", "Gamma_ij", "y_GY", "beta_GY", "contr")], pre("step3_top_lines"), row.names = FALSE)
  write.csv(det, pre("step3_records_detail"), row.names = FALSE); write.csv(vb, pre("step3_V_block"), row.names = FALSE)
  write.csv(summ, pre("step4_summary"), row.names = FALSE)
  print(summ); print(G1); print(G2); print(top)
}
run_case(G, "G"); run_case(K, "K")
