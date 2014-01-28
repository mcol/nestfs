## register the parallel back-end
.onAttach <- function(libname, pkgname) {
  registerDoMC()
  packageStartupMessage("nestfs: currently using ", getDoParWorkers(),
                        " cores, set 'options(cores=<n.cores>)' to change.")
}
