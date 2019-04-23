#' Calculate bootstrap confidence intervals for the HDR-based synchrony metric
#'
#' @param x A vector of individual completion times for the focal activity.
#' @param ci A number indicating the desired confidence level. The default is 95.
#' @param p A number indicating the percentage of individuals that have completed the focal activity. The default is 80.
#' @param nsim The number of simulations used to esimtate the HDR. The default is 5000.
#' @param nboot The number of bootstrap replicates used to estimate the CIs. The default is 5000.
#' @param ncores The number of processor cores over which to distribute the \code{nboot} replicates. The default is 1.
#' @return A named numerical vector containing the lower and upper confidence limits for the \code{hdrs} metric.
#' @examples
#' x <- rnorm(25)
#' hdrs.conf(x)
#' @export
hdrs.conf <- function(x, ci=95, p=80, nsim=5000, nboot=5000, ncores=1){

  doParallel::registerDoParallel(ncores)

  boot <- foreach::"%dopar%"(foreach::foreach(i=1:nboot, .combine='rbind'), {

    samp <- sample(x, length(x), replace=TRUE)
    hdrs(samp, p=p, nsim=nsim)

  })

  #This re-registers the sequential backend after the foreach loop is done.
  #this effectively gets rid of the doParallel backend registered above.
  foreach::registerDoSEQ()

  ql <- (1-ci/100)/2

  return(quantile(boot, c(ql, 1-ql)))

}
