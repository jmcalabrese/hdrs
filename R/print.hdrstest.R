##' @method print hdrstest
##' @export
print.hdrstest <- function(x, digits = getOption("digits"), prefix = "\t", ...)  {

  cat("\n")
  cat(strwrap(x$method, prefix = prefix), sep = "\n")
  cat("\n")
  guts <- character()
  if (!is.null(x$hdr.prob))
    guts <- c(guts, paste(names(x$hdr.prob), "=", format(signif(x$hdr.prob,
                                                                max(1L, digits - 2L)))), "\n")
  if (!is.null(x$statistic))
    guts <- c(guts, paste(names(x$statistic), "=", format(signif(x$statistic,
                                                                 max(1L, digits - 2L)))), "\n")
  if (!is.null(x$p.value)) {
    min.asl <- 1/as.numeric(x$n.perm)
    fp <- format(x$p.value, digits = max(1L, digits - 3L))
    #guts <- c(guts, paste(names(x$p.value), "=", fp), "\n")}
    if(as.numeric(x$p.value)>=min.asl){
      guts <- c(guts, paste(names(x$p.value), "=", fp), "\n")
    }else{
      guts <- c(guts, paste(names(x$p.value), "<", min.asl), "\n")
    }
  }

  if (!is.null(x$n.perm))
    guts <- c(guts, paste(names(x$n.perm), "=", format(signif(x$n.perm,
                                                              max(1L, digits - 2L)))), "\n")
  cat(strwrap(guts), sep = "\n")
  if (!is.null(x$alternative)) {
    cat("Alternative hypothesis:", "\n")
    cat(strwrap(names(x$alternative), prefix = prefix), sep = "\n")
  }

}
