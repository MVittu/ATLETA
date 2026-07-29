# Run this from a SECOND R session while analysis_hsgp.R (or any rstan fit
# using sample_file) is running in the first one, e.g.:
#
#   Rscript experiments/multilevel/watch_progress.R \
#     experiments/multilevel/runs/07_hsgp/01_primary/association/models/chain.csv 4 2000
#
# Stan writes each chain's draws to its CSV file one row per saved iteration
# as sampling proceeds (this is the same mechanism CmdStan users rely on for
# `tail -f`), so counting rows gives a live, accurate progress readout even
# when rstan's own `refresh` console text isn't visible (e.g. parallel chains
# on Windows, run via Rscript in the background).

count_stan_rows <- function(path) {
  if (!file.exists(path)) return(0L)
  lines <- suppressWarnings(readLines(path, warn = FALSE))
  max(0L, sum(nzchar(lines) & !startsWith(lines, "#")) - 1L)  # minus header
}

watch_stan_progress <- function(sample_file_stub, chains, iter, refresh_seconds = 3) {
  stub <- sub("\\.csv$", "", sample_file_stub)
  files <- if (chains == 1L) paste0(stub, ".csv") else paste0(stub, "_", seq_len(chains), ".csv")
  repeat {
    counts <- pmax(vapply(files, count_stan_rows, integer(1)), 0L)
    pct <- pmin(100L, round(100 * counts / iter))
    line <- paste(sprintf("chain %d: %4d/%d (%3d%%)", seq_along(files), counts, iter, pct),
                  collapse = "   ")
    cat("\r", line, sep = "")
    utils::flush.console()
    if (all(counts >= iter)) break
    Sys.sleep(refresh_seconds)
  }
  cat("\n")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 3) {
  watch_stan_progress(args[1], as.integer(args[2]), as.integer(args[3]))
} else if (length(args) > 0) {
  stop("usage: Rscript watch_progress.R <sample_file_stub.csv> <chains> <iter>")
}
