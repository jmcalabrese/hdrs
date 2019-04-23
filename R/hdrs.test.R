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
#' @examples
#' x <- rnorm(25)
#' y <- rnorm(25, sd=1.5)
#' hdrs.conf(x, y, alternative="more.sync")
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
