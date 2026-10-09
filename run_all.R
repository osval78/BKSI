###############################################################################
# run_all.R (version 6) -- reproduces the complete study (all models fitted with BGLR).
# Usage (from this folder):  Rscript run_all.R            # everything (about 3 h on one processor)
#                            Rscript run_all.R summary    # only tables and figures from stored raw results
# Requirements: R >= 4.0, packages BGLR and coda.
# Optional: the Word manuscript is assembled by manuscript_build/sidx_build.py (Python 3, pandoc, python-docx, pandas);
#           this step is not part of the statistical analysis.
# Partitions are stored one file per partition, so an interrupted run resumes where it stopped.
###############################################################################
args <- commandArgs(TRUE); only_summary <- length(args) && args[1] == "summary"
root <- normalizePath("."); code <- file.path(root, "code"); data <- file.path(root, "data"); out <- file.path(root, "outputs")
run <- function(script, ..., env = character()) {
  cmd <- paste(paste(env, collapse = " "), shQuote(file.path(R.home("bin"), "Rscript")), shQuote(file.path(code, script)), paste(shQuote(c(...)), collapse = " "))
  cat("\n==>", script, ..., "\n"); stopifnot(system(cmd) == 0)
}
if (!only_summary) {
  # CV2, groundnut and rice: 10 partitions, 6,000 iterations (burn-in 1,000), all methods
  run("sidx_index.R", data, out, "CV2", 10, 6000, 1000, "Groudnut", "Indica", "Japonica")
  # CV2, maize (one dataset per year): 5 partitions (3 were completed for LC2018), 4,000 iterations (burn-in 800), BKSI-GK skipped
  run("sidx_index.R", data, out, "CV2", 5, 4000, 800, "LC2020", "LC2019", "LC2017", env = "SIDX_NO_GK=1")
  run("sidx_index.R", data, out, "CV2", 3, 4000, 800, "LC2018", env = "SIDX_NO_GK=1")
  # simulation with known net merit
  run("sidx_simulation.R", data, out, 10, 5000, 1000)
  # step-by-step illustration (indica, partition 1)
  run("sidx_example.R", data, out)
  # CV1 (first 5 partitions)
  run("sidx_index.R", data, out, "CV1", 5, 6000, 1000, "Groudnut", "Indica", "Japonica")
  run("sidx_index.R", data, out, "CV1", 5, 4000, 800, "LC2020")
  # MCMC diagnostics
  run("sidx_diagnostics.R", data, out, env = "SIDX_NLC=1")
}
run("sidx_summary.R", data, out)
if (nzchar(Sys.which("python3")) && nzchar(Sys.which("pandoc")))
  system(paste("python3", shQuote(file.path(root, "manuscript_build", "sidx_build.py")), shQuote(file.path(root, "manuscript_src")),
               shQuote(out), shQuote(file.path(root, "manuscript"))))
cat("\nFINISHED\n")
