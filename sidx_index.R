###############################################################################
# sidx_index.R (version 6) -- Classical Smith-Hazel (SH) index versus the Generalized Smith-Hazel (GSH) index
# (multi-trait kernel BLUP of net merit H = w'g). All models are fitted with BGLR.
#
# Changes with respect to version 5 (in response to the internal review):
#   * traits are standardized with the mean and SD of the TRAINING lines of each partition (no use of testing data)
#   * 10 partitions in every dataset (maize included)
#   * new methods: SH-r (classical index with Sigma_g, Sigma_e rescaled to the observed phenotypic variances) and
#     ST-I (single-trait GBLUP fitted to the phenotypic index w'y of the training lines)
#   * scenario CV1 (all traits of the testing lines masked) in addition to CV2 (only the main trait masked)
#   * predictive abilities are also stored for an alternative set of economic weights (equal absolute weights)
#   * running time of every partition is stored; MCMC samples of Sigma_g and Sigma_e (partition 1, GSH-G) are kept
#
# Datasets: Groudnut (main trait YPH), Indica and Japonica (GY), LC maize (each year a dataset; GY_Optimal)
# Kernels:  G = Z Z'/p (standardized markers, MAF > 0.05); K_ij = exp(-d_ij^2 / median(d^2)), d_ij^2 = G_ii + G_jj - 2 G_ij
# Methods CV2: SH, SH-r, ST-G, ST-I, GSH-G, GSH-K, GSH-GK.   Methods CV1: ST-G, ST-I, GSH-G, GSH-K.
# Usage: Rscript sidx_index.R <data_dir> <out_dir> <scenario CV2|CV1> <n_reps> <nIter> <burnIn> <dataset> [dataset ...]
###############################################################################
suppressPackageStartupMessages(library(BGLR))
args <- commandArgs(TRUE); DD <- args[1]; OUT <- args[2]; SC <- args[3]; R <- as.integer(args[4])
nIter <- as.integer(args[5]); burnIn <- as.integer(args[6]); DS <- args[-(1:6)]
RAW <- file.path(OUT, paste0("raw_", SC)); dir.create(RAW, FALSE, TRUE)
dir.create(file.path(OUT, "partitions"), FALSE, TRUE); dir.create(file.path(OUT, "mcmc"), FALSE, TRUE)
tmp <- tempfile("bglr_"); dir.create(tmp)

kernels <- function(X) {
  X <- as.matrix(X); f <- if (max(X, na.rm = TRUE) <= 1) colMeans(X) else colMeans(X) / 2   # allele frequency (0/1 or 0/1/2 coding)
  keep <- pmin(f, 1 - f) > 0.05 & apply(X, 2, sd) > 0
  Z <- scale(X[, keep], TRUE, TRUE); G <- tcrossprod(Z) / ncol(Z)
  D2 <- outer(diag(G), diag(G), "+") - 2 * G
  K <- exp(-D2 / median(D2[upper.tri(D2)]))
  list(G = G, K = K, p = sum(keep))
}
weights_for <- function(tr) {
  w <- c(YPH = .5, NPP = 1/6, PYPP = 1/6, SYPP = 1/6,
         GY = .5, PHR = .2, GC = -.15, PH = -.15,
         GY_Optimal = .4, GY_Drought = .2, PH_Optimal = -.1, PH_Drought = -.1, AD_Optimal = -.1, AD_Drought = -.1)
  w[tr]
}
equal_weights <- function(w) sign(w) / length(w)                       # alternative objective: equal absolute weights
load_ds <- function(ds) {
  if (startsWith(ds, "LC")) {
    yr <- as.integer(sub("LC", "", ds))
    Yd <- read.csv(file.path(DD, "Y_LC.csv"), check.names = FALSE); rows <- which(Yd$Year == yr)
    X <- read.csv(file.path(DD, "X_LC.csv"), row.names = 1, check.names = FALSE)[rows, ]
    Y <- as.matrix(Yd[rows, c("GY_Optimal", "PH_Optimal", "AD_Optimal", "GY_Drought", "PH_Drought", "AD_Drought")])
    ids <- as.character(Yd$Hybrid_name[rows]); target <- "GY_Optimal"
  } else {
    X <- read.csv(file.path(DD, sprintf("X_%s.csv", ds)), row.names = 1, check.names = FALSE)
    Yd <- read.csv(file.path(DD, sprintf("Y_%s.csv", ds)), check.names = FALSE)
    tr <- intersect(names(Yd), c("GC", "GY", "PH", "PHR", "NPP", "PYPP", "SYPP", "YPH")); Y <- as.matrix(Yd[, tr]); ids <- as.character(Yd$Line)
    target <- if ("GY" %in% tr) "GY" else "YPH"
  }
  ok <- rowSums(!is.na(Y)) > 0; Y <- Y[ok, , drop = FALSE]; X <- X[ok, ]; ids <- ids[ok]
  c(list(Y = Y, ids = ids, target = target), kernels(X))
}
mt_fit <- function(Yna, Ks, keep_mcmc = NULL) {
  pre <- file.path(tmp, "mt_")
  fm <- Multitrait(y = Yna, ETA = lapply(Ks, function(K) list(K = K, model = "RKHS")), nIter = nIter,
                   burnIn = burnIn, verbose = FALSE, saveAt = pre)
  if (!is.null(keep_mcmc)) { file.copy(paste0(pre, "Omega_1.dat"), paste0(keep_mcmc, "_Omega.dat"), TRUE)
                             file.copy(paste0(pre, "R.dat"), paste0(keep_mcmc, "_R.dat"), TRUE) }
  list(U = Reduce(`+`, lapply(fm$ETA, function(e) e$u)), Omega = Reduce(`+`, lapply(fm$ETA, function(e) e$Cov$Omega)),
       R = fm$resCov$R)
}
st_fit <- function(y, EV) BGLR(y = y, ETA = list(list(V = EV$vectors, d = EV$values, model = "RKHS")), nIter = nIter,
                               burnIn = burnIn, verbose = FALSE, saveAt = file.path(tmp, "st_"))$ETA[[1]]$u
sh_index <- function(Yna, rows, Sg, Se) {                                # E[g_i | own observed records], unrelated individuals
  P <- Sg + Se; mu <- colMeans(Yna, na.rm = TRUE); out <- matrix(NA, nrow(Yna), ncol(Yna))
  for (i in rows) { o <- which(!is.na(Yna[i, ]))
    out[i, ] <- if (length(o)) Sg[, o, drop = FALSE] %*% solve(P[o, o], Yna[i, o] - mu[o]) else 0 }
  out
}
for (ds in DS) {
  D <- load_ds(ds); n <- nrow(D$Y); t <- ncol(D$Y); tg <- which(colnames(D$Y) == D$target); aux <- setdiff(seq_len(t), tg)
  w <- weights_for(colnames(D$Y)); w2 <- equal_weights(w); EVG <- eigen(D$G, symmetric = TRUE)
  message(sprintf("[%s %s] n=%d traits=%s target=%s markers=%d", ds, SC, n, paste(colnames(D$Y), collapse = ","), D$target, D$p))
  part <- data.frame(line = D$ids)
  for (rep in seq_len(R)) {
    set.seed(3000 + rep)
    te <- sort(sample(n, round(0.2 * n))); tr <- setdiff(seq_len(n), te)        # identical partitions in CV2 and CV1
    part[[sprintf("rep%02d", rep)]] <- ifelse(seq_len(n) %in% te, "testing", "training")
    f <- file.path(RAW, sprintf("%s_rep%02d.csv", ds, rep)); if (file.exists(f)) next
    t0 <- Sys.time()
    # standardization with the training lines only
    m_tr <- colMeans(D$Y[tr, , drop = FALSE], na.rm = TRUE); s_tr <- apply(D$Y[tr, , drop = FALSE], 2, sd, na.rm = TRUE)
    Ys <- sweep(sweep(D$Y, 2, m_tr), 2, s_tr, "/")
    Yna <- Ys; if (SC == "CV2") Yna[te, tg] <- NA else Yna[te, ] <- NA
    kp <- if (rep == 1) file.path(OUT, "mcmc", sprintf("%s_%s_GSH-G", ds, SC)) else NULL
    fG <- mt_fit(Yna, list(D$G), kp); fK <- mt_fit(Yna, list(D$K))
    U <- list("GSH-G" = fG$U, "GSH-K" = fK$U, "ST-G" = sapply(seq_len(t), function(k) st_fit(Yna[, k], EVG)))
    if (SC == "CV2" && Sys.getenv("SIDX_NO_GK") != "1") U[["GSH-GK"]] <- mt_fit(Yna, list(D$G, D$K))$U   # GSH-GK skipped in maize (cost)
    Sg <- fG$Omega; Se <- fG$R; dimnames(Sg) <- dimnames(Se) <- list(colnames(Ys), colnames(Ys))
    # rescaled parameters: same heritabilities and correlations, variances equal to the observed phenotypic variances
    vobs <- apply(Yna[tr, , drop = FALSE], 2, var, na.rm = TRUE); S <- diag(sqrt(vobs / (diag(Sg) + diag(Se))))
    Sg_r <- S %*% Sg %*% S; Se_r <- S %*% Se %*% S
    if (SC == "CV2") { U[["SH"]] <- sh_index(Yna, te, Sg, Se); U[["SH-r"]] <- sh_index(Yna, te, Sg_r, Se_r) }
    # ST-I: single-trait GBLUP fitted to the phenotypic index of the training lines (w and w2)
    yI <- as.numeric(Yna %*% w); yI2 <- as.numeric(Yna %*% w2)                # NA for testing lines (main trait masked)
    uI <- st_fit(yI, EVG); uI2 <- rep(NA, n)                                # ST-I not refitted for the alternative weights
    Hobs <- as.numeric(Ys[te, ] %*% w); Hobs2 <- as.numeric(Ys[te, ] %*% w2); yo <- Ys[te, tg]; res <- list()
    for (m in c(names(U), "ST-I")) {
      if (m == "ST-I") { Hh <- uI[te]; Hh2 <- uI2[te]; gt <- rep(NA, length(te)) } else {
        Hh <- as.numeric(U[[m]][te, , drop = FALSE] %*% w); Hh2 <- as.numeric(U[[m]][te, , drop = FALSE] %*% w2); gt <- U[[m]][te, tg] }
      okk <- !is.na(Hobs) & !is.na(Hh); H1 <- Hobs[okk]; top <- order(-Hh[okk])[1:round(0.2 * sum(okk))]
      ok2 <- !is.na(Hobs2) & !is.na(Hh2)
      res[[m]] <- data.frame(dataset = ds, scenario = SC, target = D$target, rep = rep, method = m, n = n, n_train = length(tr),
                             n_test = sum(okk), r = cor(Hh[okk], H1), r_target = if (all(is.na(gt))) NA else cor(gt, yo, use = "complete.obs"),
                             sel_diff = (mean(H1[top]) - mean(H1)) / sd(H1), r_w2 = if (sum(ok2) > 2) cor(Hh2[ok2], Hobs2[ok2]) else NA)
    }
    out <- do.call(rbind, res)
    out$Sigma_g <- paste(round(Sg[upper.tri(Sg, TRUE)], 5), collapse = ";"); out$Sigma_e <- paste(round(Se[upper.tri(Se, TRUE)], 5), collapse = ";")
    out$var_obs <- paste(round(vobs, 5), collapse = ";"); out$traits <- paste(colnames(Ys), collapse = ";")
    out$b_SH <- paste(round(solve((Sg + Se)[aux, aux], Sg[aux, ] %*% w), 5), collapse = ";")
    out$b_SHr <- paste(round(solve((Sg_r + Se_r)[aux, aux], Sg_r[aux, ] %*% w), 5), collapse = ";")
    out$weights <- paste(round(w, 5), collapse = ";"); out$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    write.csv(out, f, row.names = FALSE)
    message(sprintf("[%s %s] rep %d done in %.0f s", ds, SC, rep, out$seconds[1]))
  }
  write.csv(part, file.path(OUT, "partitions", sprintf("%s_partitions.csv", ds)), row.names = FALSE)
}
