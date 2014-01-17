## register the parallel back-end
.onAttach <- function(libname, pkgname) {
  registerDoMC()
  cat("Currently using", getDoParWorkers(), "cores. ")
  cat("Set 'options(cores=<n.cores>)' to change.\n")
}
