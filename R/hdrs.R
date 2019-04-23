#' Calculate the HDR-based synchrony metric
#'
#' @param p A number indicating the percentage of individuals that have completed the focal activity. The default is 80.
#' @param nsim The number of simulations used to esimtate the HDR. The default is 5000.
#' @return A number representing the length of time it takes \code{p} percent of indivdiuals to complete the focal activity.
#' @note The \code{hdrs} metric quantifies the amount of time it took for \code{p} percent of the sampled individuals to
#' complete the specified activity (e.g., mating, calving, etc...). Therefore, a smaller \code{hdrs} value indicates a \emph{more}
#' synchronous population, while a larger \code{hdrs} value indicates \emph{less} synchrony.
#' @examples
#' x <- rexp(25)
#' hdrs(x)
#' @export
hdrs <- function(x, p=80, nsim=5000){

  #Get the hdr, which may consist of multiple subintervals if the data for the
  #focal dataset are multimodal.
  ints <- hdrcde::hdr(x, prob=p, nn=nsim)

  #Indexes for the vector of hdr intervals
  idx <- 1:length(ints$hdr)

  #Get the odd and even indexes
  odd <- idx %% 2 == 1
  evn <- idx %% 2 == 0

  #The even position elements of the hdr intervals are the end points,
  #while odd position elements are start points, so subtracting odd
  #elements from even elements results in a vector of the lengths of subintervals
  #that together make up the hdr. Summing this vector gives the hdr length
  return(sum(ints$hdr[evn] - ints$hdr[odd]))

}
