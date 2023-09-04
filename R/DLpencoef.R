######
#' @title DLpencoef
#'
#' @description Compute coefficients via the Durbin Levinson algorithm using the penalized partial autocorrelation function (PACF) estimation to obtain an penalized autocorrelation function (ACF) estimation which is positive definite.
#'
#' @param x a univariate or multivariate numeric time series.
#' @param lag.max maximum lag at which to calculate the coefficients. Defaults to the smaller of \eqn{N-1} and \eqn{10log_{10}(N/nser)} where N is the number of non-missing observations and nser is the number of series.
#' @param ... additional arguments for penalized PACF estimation.
#'
#' @return A coefficients vector to estimate a penalized autocorrelation function which is positive definite.
#'
#'
#' @examples
#' \dontrun{
#' data <- arima.sim(n=100, model=list(ar=0.5))
#'
#' DLpencoef(data)
#' DLpencoef(data, lag.max=10)
#' }
#' @export
#####

DLpencoef <- function(x, lag.max = NULL, ...){
  ppacf <- penpacf(x, lag.max = lag.max, plot = FALSE, print = FALSE, ...)
  penpacf0 <- ppacf$acf
  lag <- length(ppacf$lag)
  series <- ppacf$series
  coef <- 1L:lag
  coef[1] <- penpacf0[1]
  if(lag == 1) return(coef[1])
  for(i in 2:lag){
    penpacf <- penpacf(x, lag.max = i, plot = FALSE, print = FALSE, ...)$acf
    coef[i] <- penpacf[i]
    coef[1:(i-1)] <- coef[1:(i-1)] - penpacf[i] * coef[(i-1):(1)]
  }
  coef
}
