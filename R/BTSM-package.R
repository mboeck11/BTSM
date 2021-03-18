#' @name BTSM
#' @title BTSM: Bayesian Time Series Models
#' @description This packages is a compendium of Bayesian Time Series Models...
#' @docType package
#' @importFrom stats rexp var
#' @useDynLib BTSM, .registration=TRUE
NULL

.onAttach <- function(lib, pkg) {
  if(interactive() || getOption("verbose")){
    packageStartupMessage(sprintf("Package %s %s attached. To cite, see citation(\"%s\").", pkg, utils::packageDescription(pkg)$Version, pkg))
  }
}

.onUnload <- function (libpath) {
  library.dynam.unload("BTSM", libpath)
}
