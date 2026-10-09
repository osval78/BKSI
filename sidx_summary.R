###############################################################################
# sidx_summary.R (version 6) -- tables and figures (base R). Usage: Rscript sidx_summary.R <data_dir> <out_dir>
# Paired comparisons use the corrected resampled t-test of Nadeau and Bengio (2003):
#   Var(d_bar) = (1/R + n_test/n_train) s_d^2,  CI = d_bar +/- t_{0.975,R-1} sqrt(Var),  p = 2 P(T_{R-1} > |d_bar|/sqrt(Var)).
###############################################################################
args <- commandArgs(TRUE); DD <- args[1]; OUT <- args[2]
TD <- file.path(OUT, "tables"); FD <- file.path(OUT, "figures"); dir.create(TD, FALSE, TRUE); dir.create(FD, FALSE, TRUE)
DS0 <- c("Groudnut", "Indica", "Japonica", "LC2017", "LC2018", "LC2019", "LC2020")
LAB <- c(Groudnut = "Groundnut", Indica = "Indica", Japonica = "Japonica", LC2017 = "LC2017", LC2018 = "LC2018", LC2019 = "LC2019", LC2020 = "LC2020")
M2 <- c("SH", "SH-r", "ST-G", "ST-I", "GSH-G", "GSH-K", "GSH-GK"); M1 <- c("ST-G", "ST-I", "GSH-G", "GSH-K")
cols <- c("SH" = "#7F7F7F", "SH-r" = "#BDBDBD", "ST-G" = "#E69F00", "ST-I" = "#CC79A7", "GSH-G" = "#0072B2", "GSH-K" = "#56B4E9", "GSH-GK" = "#009E73")
ml <- c("SH" = "Classical Smith-Hazel (SH)", "SH-r" = "SH, rescaled parameters (SH-r)", "ST-G" = "Weighted single-trait GBLUP (ST-G)",
        "ST-I" = "GBLUP of the phenotypic index (ST-I)", "GSH-G" = "GSH index, G (GSH-G)", "GSH-K" = "GSH index, K (GSH-K)",
        "GSH-GK" = "GSH index, G + K (GSH-GK)")
rd <- function(sc) { fl <- list.files(file.path(OUT, paste0("raw_", sc)), "_rep.*csv$", full.names = TRUE)
  if (!length(fl)) return(NULL); do.call(rbind, lapply(fl, read.csv, check.names = FALSE)) }
ci <- function(x) { x <- x[!is.na(x)]; m <- mean(x); if (length(x) < 2) return(c(m, NA, NA)); h <- qt(.975, length(x) - 1) * sd(x) / sqrt(length(x)); c(m, m - h, m + h) }
nb <- function(d, ratio) { R <- length(d); m <- mean(d); v <- (1 / R + ratio) * var(d); se <- sqrt(v); h <- qt(.975, R - 1) * se
  c(gain = m, lo = m - h, hi = m + h, p = if (R > 1 && se > 0) 2 * pt(-abs(m / se), R - 1) else NA) }
summarize <- function(d, M, sc, bases) {
  DS <- intersect(DS0, unique(d$dataset)); rows <- list(); gains <- list()
  for (ds in DS) { s <- d[d$dataset == ds, ]; ratio <- s$n_test[1] / s$n_train[1]
    for (m in intersect(M, unique(s$method))) { x <- s[s$method == m, ]; a <- ci(x$r); g <- ci(x$r_target)
      rows[[length(rows) + 1]] <- data.frame(scenario = sc, dataset = ds, method = m, reps = nrow(x), r = a[1], r_lo = a[2], r_hi = a[3],
        r_target = g[1], r_target_lo = g[2], r_target_hi = g[3], sel_diff = mean(x$sel_diff), r_w2 = mean(x$r_w2), seconds = mean(x$seconds)) }
    for (base in bases) for (m in setdiff(intersect(M, unique(s$method)), base)) {
      x <- merge(s[s$method == m, c("rep", "r")], s[s$method == base, c("rep", "r")], by = "rep"); df <- x$r.x - x$r.y; z <- nb(df, ratio)
      gains[[length(gains) + 1]] <- data.frame(scenario = sc, dataset = ds, baseline = base, method = m, reps = nrow(x), gain = z[["gain"]],
        lo = z[["lo"]], hi = z[["hi"]], rel = 100 * mean(df) / mean(x$r.y), win = mean(df > 0), p = z[["p"]]) } }
  list(TA = do.call(rbind, rows), TB = do.call(rbind, gains), DS = DS)
}
across_tab <- function(TA, TB, M, base) {
  W <- reshape(TA[, c("dataset", "method", "r")], idvar = "dataset", timevar = "method", direction = "wide"); names(W) <- sub("^r\\.", "", names(W))
  WT <- reshape(TA[, c("dataset", "method", "r_target")], idvar = "dataset", timevar = "method", direction = "wide"); names(WT) <- sub("^r_target\\.", "", names(WT))
  M <- intersect(M, names(W)); best <- table(factor(M[apply(W[, M], 1, which.max)], levels = M))
  data.frame(method = M, n_datasets = sapply(M, function(m) sum(!is.na(W[[m]]))), mean_r = sapply(M, function(m) mean(W[[m]], na.rm = TRUE)), mean_r_target = sapply(M, function(m) mean(WT[[m]], na.rm = TRUE)),
    gain_vs_base = sapply(M, function(m) mean(W[[m]] - W[[base]], na.rm = TRUE)), rel_vs_base = sapply(M, function(m) 100 * mean((W[[m]] - W[[base]]) / W[[base]], na.rm = TRUE)),
    n_better = sapply(M, function(m) sum(W[[m]] > W[[base]], na.rm = TRUE)),
    n_signif_better = sapply(M, function(m) if (m == base) NA else sum(TB$baseline == base & TB$method == m & TB$gain > 0 & TB$p < 0.05, na.rm = TRUE)),
    n_signif_worse = sapply(M, function(m) if (m == base) NA else sum(TB$baseline == base & TB$method == m & TB$gain < 0 & TB$p < 0.05, na.rm = TRUE)),
    n_best = as.numeric(best[M]))
}
d2 <- rd("CV2"); S2 <- summarize(d2, M2, "CV2", c("SH", "ST-G")); DS <- S2$DS
write.csv(S2$TA, file.path(TD, "accuracy_CV2.csv"), row.names = FALSE); write.csv(S2$TB, file.path(TD, "gains_CV2.csv"), row.names = FALSE)
write.csv(across_tab(S2$TA, S2$TB, M2, "SH"), file.path(TD, "across_CV2.csv"), row.names = FALSE)
d1 <- rd("CV1")
if (!is.null(d1)) { S1 <- summarize(d1, M1, "CV1", "ST-G")
  write.csv(S1$TA, file.path(TD, "accuracy_CV1.csv"), row.names = FALSE); write.csv(S1$TB, file.path(TD, "gains_CV1.csv"), row.names = FALSE)
  write.csv(across_tab(S1$TA, S1$TB, M1, "ST-G"), file.path(TD, "across_CV1.csv"), row.names = FALSE) }
# ---------------- dataset description
info <- list()
for (ds in DS) {
  if (startsWith(ds, "LC")) {
    Yd <- read.csv(file.path(DD, "Y_LC.csv"), check.names = FALSE); rows <- which(Yd$Year == as.integer(sub("LC", "", ds)))
    X <- as.matrix(read.csv(file.path(DD, "X_LC.csv"), row.names = 1, check.names = FALSE)[rows, ])
    Y <- Yd[rows, c("GY_Optimal", "PH_Optimal", "AD_Optimal", "GY_Drought", "PH_Drought", "AD_Drought")]; tg <- "GY_Optimal"
  } else {
    X <- as.matrix(read.csv(file.path(DD, sprintf("X_%s.csv", ds)), row.names = 1, check.names = FALSE))
    Yd <- read.csv(file.path(DD, sprintf("Y_%s.csv", ds)), check.names = FALSE)
    Y <- Yd[, intersect(names(Yd), c("GC", "GY", "PH", "PHR", "NPP", "PYPP", "SYPP", "YPH"))]; tg <- if ("GY" %in% names(Y)) "GY" else "YPH"
  }
  ok <- rowSums(!is.na(Y)) > 0; X <- X[ok, ]; Y <- Y[ok, ]
  f <- if (max(X, na.rm = TRUE) <= 1) colMeans(X) else colMeans(X) / 2; p <- sum(pmin(f, 1 - f) > 0.05 & apply(X, 2, sd) > 0)
  Z <- scale(X[, pmin(f, 1 - f) > 0.05 & apply(X, 2, sd) > 0]); dG <- rowSums(Z^2) / ncol(Z)
  rc <- cor(Y, use = "pairwise.complete.obs")[tg, ]; aux <- setdiff(names(Y), tg); n <- nrow(Y); nv <- round(0.2 * n)
  info[[ds]] <- data.frame(dataset = ds, n = n, markers_total = ncol(X), markers = p, target = tg,
    auxiliary = paste(sprintf("%s (%.2f)", aux, rc[aux]), collapse = "; "), n_train = n - nv, n_test = nv,
    partitions = length(unique(d2$rep[d2$dataset == ds])), partitions_CV1 = if (is.null(d1)) 0 else length(unique(d1$rep[d1$dataset == ds])),
    diagG_mean = mean(dG), diagG_min = min(dG), diagG_max = max(dG), n_abs_z_gt4 = sum(abs(scale(Y)) > 4, na.rm = TRUE), n_missing = sum(is.na(Y)))
}
write.csv(do.call(rbind, info), file.path(TD, "datasets.csv"), row.names = FALSE)
# ---------------- genetic parameters (multi-trait G model, averaged over partitions) and SH coefficients
pars <- list()
for (ds in DS) {
  s <- d2[d2$dataset == ds & !duplicated(d2[, c("dataset", "rep")]), ]; tr <- strsplit(s$traits[1], ";")[[1]]; k <- length(tr)
  tomat <- function(v) { v <- as.numeric(strsplit(v, ";")[[1]]); Mx <- matrix(0, k, k); Mx[upper.tri(Mx, TRUE)] <- v; Mx[lower.tri(Mx)] <- t(Mx)[lower.tri(Mx)]; Mx }
  Sg <- Reduce(`+`, lapply(s$Sigma_g, tomat)) / nrow(s); Se <- Reduce(`+`, lapply(s$Sigma_e, tomat)) / nrow(s)
  vob <- colMeans(do.call(rbind, lapply(s$var_obs, function(v) as.numeric(strsplit(v, ";")[[1]]))))
  bb <- colMeans(do.call(rbind, lapply(s$b_SH, function(v) as.numeric(strsplit(v, ";")[[1]]))))
  bbr <- colMeans(do.call(rbind, lapply(s$b_SHr, function(v) as.numeric(strsplit(v, ";")[[1]])))); tg <- which(tr == s$target[1]); aux <- setdiff(seq_len(k), tg)
  wv <- as.numeric(strsplit(s$weights[1], ";")[[1]])
  pars[[ds]] <- data.frame(dataset = ds, trait = tr, target = tr == s$target[1], w = wv, h2 = diag(Sg) / (diag(Sg) + diag(Se)),
    var_model = diag(Sg) + diag(Se), var_obs = vob, rg_target = Sg[, tg] / sqrt(diag(Sg) * Sg[tg, tg]), re_target = Se[, tg] / sqrt(diag(Se) * Se[tg, tg]),
    b_SH = replace(rep(NA, k), aux, bb), b_SHr = replace(rep(NA, k), aux, bbr))
}
write.csv(do.call(rbind, pars), file.path(TD, "parameters.csv"), row.names = FALSE)
# ---------------- figures
dotplot <- function(file, TA, M, col, lo, hi, ylab, ylim = c(0, 1), sep = 3.5, legpos = "bottomleft", h = 2884) {
  png(file.path(FD, file), width = 5538, height = h, res = 600); par(mar = c(3, 4.5, .8, .8)); off <- seq(-.4, .4, length.out = length(M))
  DSx <- intersect(DS0, unique(TA$dataset))
  plot(NA, xlim = c(.5, length(DSx) + .5), ylim = ylim, xaxt = "n", xlab = "", ylab = ylab, las = 1, cex.lab = .9, cex.axis = .8)
  axis(1, seq_along(DSx), LAB[DSx], cex.axis = .8)
  for (j in seq_along(M)) for (k in seq_along(DSx)) { s <- TA[TA$dataset == DSx[k] & TA$method == M[j], ]; if (!nrow(s) || is.na(s[[col]])) next
    if (!is.na(s[[lo]])) arrows(k + off[j], s[[lo]], k + off[j], s[[hi]], angle = 90, code = 3, length = .012, col = cols[M[j]])
    points(k + off[j], s[[col]], pch = 19, col = cols[M[j]], cex = .75) }
  if (!is.null(sep)) abline(v = sep, lty = 2, col = "grey70"); if (any(ylim < 0)) abline(h = 0)
  legend(legpos, ml[M], col = cols[M], pch = 19, bty = "n", cex = .6, ncol = 3); dev.off()
}
dotplot("Fig2_predictive_ability_index_CV2.png", S2$TA, M2, "r", "r_lo", "r_hi", "Predictive ability of the index", c(0, 1))
dotplot("Fig4_predictive_ability_main_trait_CV2.png", S2$TA[S2$TA$method != "ST-I", ], setdiff(M2, "ST-I"), "r_target", "r_target_lo", "r_target_hi",
        "Predictive ability, masked main trait", c(0, 1))
g <- S2$TB[S2$TB$baseline == "SH", ]; names(g)[names(g) == "gain"] <- "r"
dotplot("Fig3_gain_vs_SH_CV2.png", g, setdiff(M2, "SH"), "r", "lo", "hi", "Gain in predictive ability over SH", range(c(g$lo, g$hi, 0), na.rm = TRUE), legpos = "topleft")
if (!is.null(d1)) dotplot("Fig6_predictive_ability_index_CV1.png", S1$TA, M1, "r", "r_lo", "r_hi", "Predictive ability of the index (CV1)", c(0, 1))
# simulation figure
sf <- file.path(OUT, "simulation", "simulation_summary.csv")
if (file.exists(sf)) {
  sm <- read.csv(sf); MS <- c("SH-true", "SH", "ST-G", "GSH-G"); cs <- c("SH-true" = "#4D4D4D", "SH" = "#7F7F7F", "ST-G" = "#E69F00", "GSH-G" = "#0072B2")
  sl <- sprintf("r_g = %.1f\nr_e = %.1f", sm$r_g[match(1:4, sm$scenario)], sm$r_e[match(1:4, sm$scenario)])
  png(file.path(FD, "Fig5_simulation.png"), width = 5538, height = 2653, res = 600); par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2, .5))
  for (pan in list(c("acc_H", "pa_H", "Net merit H"), c("acc_main", "pa_main", "Main trait"))) {
    plot(NA, xlim = c(.5, 4.5), ylim = c(-.42, 1), xaxt = "n", yaxt = "n", xlab = "", ylab = "Correlation", main = pan[3], cex.main = .9)
    axis(2, seq(0, 1, .2), las = 1, cex.axis = .75); abline(h = -.05, col = "grey80")
    axis(1, 1:4, sl, cex.axis = .62, padj = .5); off <- seq(-.3, .3, length.out = 4)
    for (j in seq_along(MS)) for (s in 1:4) { x <- sm[sm$scenario == s & sm$method == MS[j], ]
      points(s + off[j], x[[pan[1]]], pch = 19, col = cs[j]); points(s + off[j], x[[pan[2]]], pch = 1, col = cs[j]) }
    legend("bottom", c(MS, "filled: accuracy", "open: predictive ability"), col = c(cs, NA, NA), pch = c(rep(19, 4), NA, NA),
           bty = "n", cex = .55, ncol = 3) }
  dev.off()
}
# scheme of one partition
png(file.path(FD, "Fig1_scheme.png"), width = 5538, height = 2653, res = 600)
par(mar = c(1, 1, 1, 1)); plot(NA, xlim = c(0, 24), ylim = c(0, 11), axes = FALSE, xlab = "", ylab = "")
tr <- c("Aux 1", "Aux 2", "Aux 3", "Main"); ln <- c(paste("Training line", 1:4), "Testing (CV2)", "Testing (CV1)")
for (r in 1:6) for (c in 1:4) { masked <- (r == 5 && c == 4) || r == 6
  rect(c * 1.6, 10 - r * 1.4, c * 1.6 + 1.5, 11.3 - r * 1.4, col = if (masked) "#F4CCCC" else if (r > 4) "#D9EAD3" else "#CFE2F3", border = "white")
  text(c * 1.6 + .75, 10.65 - r * 1.4, if (masked) "NA" else "y", cex = .8) }
text(1.5, 10.65 - (1:6) * 1.4, ln, adj = 1, cex = .6); text((1:4) * 1.6 + .75, 10, tr, cex = .75)
text(4.8, 0.3, "Phenotypes passed to BGLR (Y with NA)", cex = .75, font = 2)
xs <- 9.5; txt <- c("SH, SH-r: own auxiliary records of the testing line only (CV2)",
  "ST-G: each trait separately (all observed records of that trait) + G",
  "ST-I: GBLUP of the phenotypic index w'y of the training lines + G",
  "GSH-G / GSH-K / GSH-GK: all observed records of all lines + kernel(s)",
  "Predictive ability: cor(predicted index, observed index w'y), testing lines")
for (i in 1:5) { rect(xs, 10.2 - i * 1.85, 23.8, 11.5 - i * 1.85, col = c("#EEEEEE", "#FCE5CD", "#F3D9E8", "#CFE2F3", "#F4CCCC")[i], border = NA)
  text(xs + .3, 10.85 - i * 1.85, txt[i], adj = 0, cex = .68) }
dev.off()
print(S2$TA[, 1:9], digits = 3); print(read.csv(file.path(TD, "across_CV2.csv")), digits = 3)
