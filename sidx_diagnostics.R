###############################################################################
# sidx_diagnostics.R -- MCMC diagnostics (answer to review point M8).
# (1) Effective sample size (coda) and Geweke z-scores of every element of Sigma_g and Sigma_e from the chains of
#     the BKSI-G model in partition 1 (CV2) of every dataset (samples saved by sidx_index.R; BGLR thin = 10).
# (2) Long-chain check: the BKSI-G model of partition 1 is refitted with 4 times more iterations in Indica (and LC2020 if SIDX_NLC=2), and the posterior means (Sigma_g, predicted index of the testing lines) are compared with the standard run.
# Usage: Rscript sidx_diagnostics.R <data_dir> <out_dir>
###############################################################################
suppressPackageStartupMessages({ library(BGLR); library(coda) })
args <- commandArgs(TRUE); DD <- args[1]; OUT <- args[2]; MC <- file.path(OUT, "mcmc"); dir.create(file.path(OUT, "tables"), FALSE, TRUE)
std <- list(Groudnut = c(6000, 1000), Indica = c(6000, 1000), Japonica = c(6000, 1000),
            LC2017 = c(4000, 800), LC2018 = c(4000, 800), LC2019 = c(4000, 800), LC2020 = c(4000, 800))
thin <- 10; rows <- list()
for (ds in names(std)) for (par in c("Omega", "R")) {
  fl <- file.path(MC, sprintf("%s_CV2_BKSI-G_%s.dat", ds, par)); if (!file.exists(fl)) next
  S <- as.matrix(read.table(fl)); S <- S[-seq_len(std[[ds]][2] / thin), , drop = FALSE]
  ess <- effectiveSize(mcmc(S)); gz <- geweke.diag(mcmc(S))$z
  rows[[length(rows) + 1]] <- data.frame(dataset = ds, matrix = ifelse(par == "Omega", "Sigma_g", "Sigma_e"), n_samples = nrow(S),
    ESS_min = min(ess), ESS_median = median(ess), geweke_max_abs_z = max(abs(gz), na.rm = TRUE),
    prop_abs_z_gt_1.96 = mean(abs(gz) > 1.96, na.rm = TRUE))
}
diag_tab <- do.call(rbind, rows); write.csv(diag_tab, file.path(OUT, "tables", "mcmc_diagnostics.csv"), row.names = FALSE); print(diag_tab)
# ---------------- long chains
src <- new.env(); sys.source(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))), "sidx_common.R"), envir = src)
lc <- list()
for (ds in c("Indica", "LC2020")[seq_len(as.integer(Sys.getenv("SIDX_NLC", "2")))]) {
  D <- src$load_ds(DD, ds); n <- nrow(D$Y); tg <- which(colnames(D$Y) == D$target); w <- src$weights_for(colnames(D$Y))
  set.seed(3001); te <- sort(sample(n, round(0.2 * n))); tr <- setdiff(seq_len(n), te)
  Ys <- sweep(sweep(D$Y, 2, colMeans(D$Y[tr, ], na.rm = TRUE)), 2, apply(D$Y[tr, ], 2, sd, na.rm = TRUE), "/"); Yna <- Ys; Yna[te, tg] <- NA
  fits <- lapply(c(1, 4), function(m) { t0 <- Sys.time()
    fm <- Multitrait(y = Yna, ETA = list(list(K = D$G, model = "RKHS")), nIter = m * std[[ds]][1], burnIn = m * std[[ds]][2],
                     verbose = FALSE, saveAt = file.path(tempdir(), sprintf("lc%d_", m)))
    list(Sg = fm$ETA[[1]]$Cov$Omega, H = as.numeric(fm$ETA[[1]]$u[te, ] %*% w), sec = as.numeric(difftime(Sys.time(), t0, units = "secs"))) })
  lc[[ds]] <- data.frame(dataset = ds, iter_standard = std[[ds]][1], iter_long = 4 * std[[ds]][1],
    max_abs_diff_Sigma_g = max(abs(fits[[1]]$Sg - fits[[2]]$Sg)), cor_Sigma_g = cor(as.vector(fits[[1]]$Sg), as.vector(fits[[2]]$Sg)),
    cor_index_testing = cor(fits[[1]]$H, fits[[2]]$H), seconds_standard = fits[[1]]$sec, seconds_long = fits[[2]]$sec)
  print(lc[[ds]])
}
write.csv(do.call(rbind, lc), file.path(OUT, "tables", "mcmc_long_chain.csv"), row.names = FALSE)
