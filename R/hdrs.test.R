#' A two-sample permutation test for differences in \code{hdrs} values.
#'
#' @param x A vector of individual completion times for the focal activity.
#' @param y A vector of individual completion times for the focal activity.
#' @param alternative A character string specifying the alternative hypothesis, must be one of "two.sided" (default), "more.sync" or "less.sync".
#' @param ci A number indicating the desired confidence level. The default is 95.
#' @param p A number indicating the percentage of individuals that have completed the focal activity. The default is 80.
#' @param nsim The number of simulations used to esimtate the HDR. The default is 5000.
#' @param nboot The number of bootstrap replicates used to estimate the CIs. The default is 5000.
#' @param ncores The number of processor cores over which to distribute the \code{nboot} replicates. The default is 1.
#' @return A list of class "hdrstest" containing the following values.
#' @note The \code{hdrs} metric quantifies the amount of time it took for \code{p} percent of the sampled individuals to
#' complete the specified activity (e.g., mating, calving, etc...). Therefore, a smaller \code{hdrs} value indicates a \emph{more}
#' synchronous population, while a larger \code{hdrs} value indicates \emph{less} synchrony. Consequently, the alternative
#' hypothesis that \code{x} is more synchronous than \code{y} ("more.sync") is equivalent to \code{hdrs(x) < hdrs(y)}, while
#' the "less.sync" alternative implies \code{hdrs(x) > hdrs(y)}. The null hypothesis is always \code{hdrs(x) = hdrs(y)}, and
#' the default "two.tailed" alternative is simply \code{hdrs(x) != hdrs(y)}.
#'
#' The \code{ncores} argument specifies the number of processor cores over which the \code{foreach} loop
#' will distribute permutation re-samples. This should be set no higher than the number of hyperthreaded
#' cores on your machine, which is usually twice the number of physical cores.
#' @examples
#' x <- rexp(25)
#' y <- rexp(25, rate=0.5)
#' hdrs.conf(x, y, alternative="more.sync", nperm=250)
#' @export
hdrs.test <- function(x, y, alternative="two.tailed", p=80, nsim=5000, nperm=5000, ncores=1){

  doParallel::registerDoParallel(ncores)

  #Get per group sample sizes
  nX <- length(x)
  nY <- length(y)

  #Total sample size
  n <- nX + nY

  #Join x and y datasets into one dataset
  dataAll <- c(x, y)

  permdiff <- foreach::"%dopar%"(foreach::foreach(i=1:nperm, .combine='rbind'), {

    #Randomly choose, w/o replacement, nx obs from the joint dataset
    inds <- sample(1:n, size=nX, replace=FALSE)

    #Assign the randomly chosen obs to the "x" perm group
    permX <- dataAll[inds]

    #Assign the remaining n - nx obs to the "y" perm group
    permY <- dataAll[-inds]

    #Estimate the hdr for both x and y perm groups and
    #store interval lengths
    hdrLenX <- hdrs(permX, p=p, nsim=nsim)
    hdrLenY <- hdrs(permY, p=p, nsim=nsim)

    #Return the test statistic, the diff in interval length
    hdrLenX - hdrLenY

  })

  #This re-registers the sequential backend after the foreach loop is done, which
  #effectively gets rid of the doParallel backend registered above.
  foreach::registerDoSEQ()

  #Calculate the observed difference in HDR synchrony between x and y.
  obsdiff <- hdrs(x, p=p, nsim=nsim) - hdrs(y, p=p, nsim=nsim)

  #Calculate the acheived significance level of the test, for the
  #appropriate alternative hypothsis
  asl <- switch(alternative,
                "more.sync" = length(permdiff[permdiff < obsdiff]) / nperm,
                "less.sync" = length(permdiff[permdiff > obsdiff]) / nperm,
                "two.tailed" = length(permdiff[abs(permdiff) > abs(obsdiff)]) / nperm
  )

  names(asl) <- "Acheived significance level"

  hdrprob <- p
  names(hdrprob) <- "HDRS percentile"

  names(obsdiff) <- "Observed difference in HDRS values"

  alt <- alternative
  names(alt) <- switch(alternative,
                       "more.sync"="x is more synchronous than y, implying hdrs(x) < hdrs(y)",
                       "less.sync"="x is less synchronous than y, implying hdrs(x) > hdrs(y)",
                       "two.tailed"="x is either more or less synchronous than y, implying hdrs(x) != hdrs(y)"
  )

  names(nperm) <- "Number of permutations"

  out <- list(
    method="Two-sample permutation test based on HDR Synchrony (HDRS)",
    hdr.prob=hdrprob,
    statistic=obsdiff,
    p.value=asl,
    alternative=alt,
    n.perm=nperm,
    distribution=permdiff
  )

  class(out)<- "hdrstest"

  return(out)

}
