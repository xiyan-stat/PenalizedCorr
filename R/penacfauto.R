######
#' @title penacfauto
#'
#' @description An automatic algorithm to estimate the ACF and PACF: (1) select the order p from 'Penar'; (2) if p less than 3, the penalized ACF caluculated by inverting the penalized PACF estimator (posdef="DL"), otherwise use the original penalized ACF estimator.
#'
#' @param x a univariate numeric time series or a numeric vector.
#' @param lag.max maximum lag at which to calculate the penalized ACF. Defaults to the smaller of \eqn{N-1} and \eqn{10log_{10}(N/nser)} where N is the number of non-missing observations and nser is the number of series.
#' @param na.action function to be called to handle missing values. 'na.pass' can be used.
#' @param ... additional arguments for penalized ACF estimation.
#'
#' @return An object of penalized ACF estimation, which is a list with the following elements:
#' \describe{
#' \item{\code{orders}}{The vector of selected orders from \code{Penar}.}
#' \item{\code{penacf}}{The penalized ACF calculated by auto algorithm.}
#' }
#'
#' @examples
#' \dontrun{
#' data <- arima.sim(n=100, model=list(ar=0.5))
#'
#' penacfauto(data)
#' }
#' @export
#####

penacfauto <- function(x, lag.max = NULL, na.action = na.fail, ...){
  x <- na.action(as.ts(x))
  x.freq <- frequency(x)
  x <- as.matrix(x)
  if(!is.numeric(x))
    stop("'x' must be numeric")
  n <- as.integer(nrow(x))
  nser <- as.integer(ncol(x))
  if (is.na(n) | is.na(nser))
    stop("'n' and 'nser' must be integer")
  penacf_auto <- list()
  orders <- NULL
  for (i in 1:nser){
    p <- penar(x, order.max = lag.max, ...)$order[i]
    if (p <= 3){
      penacf_auto[[i]] <- penacf(x[,i], lag.max = lag.max, posdef = 'DL', plot=FALSE, print=FALSE, ...)$acf[-1]
    }else{
      penacf_auto[[i]] <- penacf(x[,i], lag.max = lag.max, plot=FALSE, print=FALSE, ...)$acf[-1]
    }
    orders = c(orders, p)
  }
  res <- list(order = orders, penacf = penacf_auto)
  res
}
