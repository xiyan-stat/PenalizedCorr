######
#' @title penarprint
#'
#' @description Print of objects from "\code{penacf}" and "\code{penpacf}"
#'
#' @param x an object of "acf.out" or "pacf.out".
#' @param digits minimal number of significant digits.
#' @param ... further arguments passed to or from other methods.
#'
#' @return list of results for penalized ar
#'
#'
#' @examples
#' \dontrun{
#' data <- arima.sim(n=100, model=list(ar=0.5))
#'
#' ac.res <- penar(data)
#' penarprint(ac.res)
#' }
#' @export
#####


penarprint <- function(x, digits = 3L, ...){
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  cat("Coefficients:\n")
  print(x$ar, digits = digits)
  cat("\nOrder selected:", x$order)
  cat("\n")
  invisible(x)
}
