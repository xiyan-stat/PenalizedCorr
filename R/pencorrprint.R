######
#' @title pencorrprint
#'
#' @description Print of objects from "\code{penacf}" and "\code{penpacf}"
#'
#' @param x an object of "acf.out" or "pacf.out".
#' @param digits minimal number of significant digits.
#' @param ... further arguments passed to or from other methods.
#'
#' @return list of results for penalized acf or pacf
#'
#'
#' @examples
#' \dontrun{
#' data <- arima.sim(n=100, model=list(ar=0.5))
#' pencorrprint(penacf(data))
#' }
#' @export
#####


pencorrprint <- function(x, digits = 3L, ...){
  type <- match(x$type, c("pencorrelation", "penpartial"))
  msg <- c("Penalized Autocorrelations", "Penalized Partial autocorrelations")
  cat("\n", msg[type]," of series ", sQuote(x$series), ", by lag\n\n",
      sep = "")
  nser <- dim(x$lag)[3]
  x$acf <- round(x$acf, digits)
  if (nser == 1){
    acfs <- setNames(drop(x$acf), format(drop(x$lag), digits = 3L))
    print(acfs, digits = digits, ...)
  }else{
    acfs <- format(x$acf, ...)
    lags <- format(x$lag, digits = 3L)
    acfs <- array(paste0(acfs, "(", lags, ")"), dim=dim(x$acf))
    dimnames(acfs)  <- list(rep("", nrow(x$lag)),"", x$snames)
    print(acfs, quote = FALSE, ...)
  }
  invisible(x)
}
