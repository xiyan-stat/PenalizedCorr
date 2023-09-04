######
#' @title NNDpenacf
#'
#' @description Calculte Nonnegative definifite (NND) penalized estimator of autocorrelation function (ACF).
#'
#' @param x a univariate numeric time series or a numeric vector.
#' @param penacf vector of penalized ACF with posdef = "orginal".
#' @param lag.max maximum lag at which to calculate the penalized ACF. Defaults to the smaller of \eqn{N-1} and \eqn{10log_{10}(N/nser)} where N is the number of non-missing observations and nser is the number of series.
#' @param ... additional arguments for penalized ACF/PACF estimation.
#'
#' @return An object of penalized positive definite ACF estimation, which is a list with the following elements:
#' \describe{
#' \item{\code{rhotilde.NND}}{A vector containing the estimated penalized ACF which is Nonnegative definite.}
#' \item{\code{lag.max}}{The value of the 'lag.max' argument.}
#' }
#'
#' @examples
#' \dontrun{
#' data <- arima.sim(n=100, model=list(ar=0.5))
#'
#' penacf1 <- penacf(data, posdef="origin")
#' NNDpenacf(data, penacf1)
#' }
#' @export
#####

NNDpenacf <- function(x, penacf, lag.max = NULL, ...){
  x <- as.matrix(x)
  n <- as.integer(nrow(x))
  nser <- as.integer(ncol(x))
  if (nser > 1)
    stop("x must be univariate time series")
  lag.max <- if (is.null(lag.max))
    floor(10 * (log10(n) - log10(nser))) else round(lag.max)
  lag.max <- as.integer(min(lag.max, n - 1L))
  if (is.na(lag.max) || lag.max < 1L)
    stop("'lag.max' must be >= 1")
  rhotilde <- penacf
  rhohat <-  stats::acf(x, plot = FALSE, lag.max = lag.max)$acf[1:(lag.max+1)]
  R1 <- toeplitz(rhohat)
  R2 <- toeplitz(rhotilde)
  lam1 <- min(eigen(R1)$values)
  lam2 <- min(eigen(R2)$values)
  if(lam2 > 0){
    rhotilde.NND <- R2[1L, ][1L : (lag.max+1)]
  }else{
    alpha <- abs(lam2)/(lam1 + abs(lam2))
    tmp <- (alpha * R1 + (1L - alpha) * R2)
    rhotilde.NND <- tmp[1L, ][1L : (lag.max+1)]
  }
  res <- list(rhotilde.NND = rhotilde.NND, lag.max = lag.max)
  res
}
