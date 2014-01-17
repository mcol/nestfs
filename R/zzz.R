## register the parallel back-end
.onAttach <- function(libname, pkgname) {
  registerDoMC()
  packageStartupMessage("Currently using ", getDoParWorkers(), " cores. ",
                        "Set 'options(cores=<n.cores>)' to change.")
}
