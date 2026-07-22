configure_cpp_toolchain <- function() {
  if (.Platform$OS.type != "windows" || nzchar(Sys.getenv("R_MAKEVARS_USER"))) {
    return(invisible())
  }

  rtools_home <- Sys.getenv("RTOOLS45_HOME")
  if (!nzchar(rtools_home)) {
    rtools_home <- file.path(Sys.getenv("SystemDrive", unset = "C:"), "rtools45")
  }
  bin_dirs <- file.path(rtools_home, c("usr/bin", "x86_64-w64-mingw32.static.posix/bin"))
  if (!all(dir.exists(bin_dirs))) {
    stop("Rtools 4.5 was not found; install it or set RTOOLS45_HOME.")
  }

  Sys.setenv(PATH = paste(c(bin_dirs, Sys.getenv("PATH")), collapse = .Platform$path.sep))
  rtools_home <- normalizePath(rtools_home, winslash = "/")
  makevars <- tempfile("Makevars-rtools45-")
  writeLines(c(
    sprintf("SHELL = %s/usr/bin/sh.exe", rtools_home),
    sprintf("BINPREF = %s/x86_64-w64-mingw32.static.posix/bin/", rtools_home),
    "export PATH := /usr/bin:/bin:/x86_64-w64-mingw32.static.posix/bin"
  ), makevars)
  Sys.setenv(R_MAKEVARS_USER = makevars)
  invisible()
}

configure_cpp_toolchain()
rm(configure_cpp_toolchain)
