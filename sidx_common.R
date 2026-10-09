###############################################################################
# sidx_common.R -- helper functions shared by sidx_diagnostics.R (identical to those in sidx_index.R)
###############################################################################
kernels <- function(X) {
  X <- as.matrix(X); f <- if (max(X, na.rm = TRUE) <= 1) colMeans(X) else colMeans(X) / 2
  keep <- pmin(f, 1 - f) > 0.05 & apply(X, 2, sd) > 0
  Z <- scale(X[, keep], TRUE, TRUE); G <- tcrossprod(Z) / ncol(Z)
  D2 <- outer(diag(G), diag(G), "+") - 2 * G
  list(G = G, K = exp(-D2 / median(D2[upper.tri(D2)])), p = sum(keep))
}
weights_for <- function(tr) {
  w <- c(YPH = .5, NPP = 1/6, PYPP = 1/6, SYPP = 1/6, GY = .5, PHR = .2, GC = -.15, PH = -.15,
         GY_Optimal = .4, GY_Drought = .2, PH_Optimal = -.1, PH_Drought = -.1, AD_Optimal = -.1, AD_Drought = -.1)
  w[tr]
}
load_ds <- function(DD, ds) {
  if (startsWith(ds, "LC")) {
    yr <- as.integer(sub("LC", "", ds)); Yd <- read.csv(file.path(DD, "Y_LC.csv"), check.names = FALSE); rows <- which(Yd$Year == yr)
    X <- read.csv(file.path(DD, "X_LC.csv"), row.names = 1, check.names = FALSE)[rows, ]
    Y <- as.matrix(Yd[rows, c("GY_Optimal", "PH_Optimal", "AD_Optimal", "GY_Drought", "PH_Drought", "AD_Drought")])
    ids <- as.character(Yd$Hybrid_name[rows]); target <- "GY_Optimal"
  } else {
    X <- read.csv(file.path(DD, sprintf("X_%s.csv", ds)), row.names = 1, check.names = FALSE)
    Yd <- read.csv(file.path(DD, sprintf("Y_%s.csv", ds)), check.names = FALSE)
    tr <- intersect(names(Yd), c("GC", "GY", "PH", "PHR", "NPP", "PYPP", "SYPP", "YPH")); Y <- as.matrix(Yd[, tr]); ids <- as.character(Yd$Line)
    target <- if ("GY" %in% tr) "GY" else "YPH"
  }
  ok <- rowSums(!is.na(Y)) > 0
  c(list(Y = Y[ok, , drop = FALSE], ids = ids[ok], target = target), kernels(X[ok, ]))
}
