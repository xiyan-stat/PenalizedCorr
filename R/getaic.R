######
#' @title getaic
#'
#' @description Calculate the value of the Akaike Information Criterion (AIC) based on penalized PACF, which is used to choose the order of the autoregressive model.
#'
#' @param x a univariate or multivariate numeric time series.
#' @param lag.max maximum lag at which to calculate the penalized PACF. Defaults to the smaller of \eqn{N-1} and \eqn{10log_{10}(N/nser)} where N is the number of non-missing observations and nser is the number of series.
#' @param na.action function to be called to handle missing values. 'na.pass' can be used.
#' @param ... additional arguments for penalized PACF estimation.
#'
#' @return A vector of AIC calaulated using penalized partial autocorrelation function estimation.
#'
#' @examples
#' \dontrun{
#' data <- arima.sim(n=100, model=list(ar=0.5))
#'
#' getaic(data)
#' }
#' @export
#####

getaic <- function(x, lag.max = NULL, na.action = na.fail, ...){
  x <- na.action(as.ts(x))
  x <- as.matrix(x)
  if(!is.numeric(x))
    stop("'x' must be numeric")
  n <- as.integer(nrow(x))
  nser <- as.integer(ncol(x))
  if (is.na(n) | is.na(nser))
    stop("'n' and 'nser' must be integer")
  vaic <- NULL
  res <- penpacf(x, lag.max = lag.max, plot = FALSE, print = FALSE, ...)
  for (j in 1L:nser){
    penpacf <- res$acf[ , , j]
    lag <- length(res$lag[, , j])
    tmp <- n * log(var(x[ , j]) * prod(1 - penpacf^2)) + 2 * lag
    vaic <- c(vaic, tmp)
  }
  vaic
}
