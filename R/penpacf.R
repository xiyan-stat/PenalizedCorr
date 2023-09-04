######
#' @title penpacf
#'
#' @description The Penalized Partial Autocorrelation Function (PACF) Estimation
#'
#' @param x a univariate or multivariate numeric time series.
#' @param lag.max maximum lag at which to calculate the penalized/sample PACF. Defaults to the smaller of \eqn{N-1} and \eqn{10log_{10}(N/nser)} where N is the number of non-missing observations and nser is the number of series.
#' @param C a positive real number for l_h
#' @param penalized 'logical'. If 'TRUE' (the default) the penalized PACF is computed; if 'FALSE' the sample PACF is computed.
#' @param print  'logical'. If 'TRUE' (the default) the penalized/sample PACF is printed.
#' @param plot  'logical'. If 'TRUE' (the default) the penalized/sample PACF is plotted.
#' @param series names for the series. Defaults to 'deparse(substitute(x))'.
#' @param na.action function to be called to handle missing values. 'na.pass' can be used.
#' @param demean 'logical'. Should a mean be estimated during estimating.
#' @param ... additional arguments for specific methods.
#'
#' @return An object of penalized/sample PACF estimation, which is a list with the following elements:
#' \describe{
#' \item{\code{acf}}{An array containing the estimated penalized/sample PACF.}
#' \item{\code{lag}}{An array containing the lags at which the PACF is estimated.}
#' \item{\code{n.used}}{The number of observation in the time series.}
#' \item{\code{series}}{The name of the time series.}
#' \item{\code{snames}}{The series names for a multivariate time series.}
#' }
#'
#'
#' @examples
#' \dontrun{
#' data <- arima.sim(n=100, model=list(ar=0.5))
#'
#' # Examples for penalized PACF and sample PACF
#' penpacf(data)
#' penpacf(data, penalized = FALSE)
#' }
#' @export
#####

penpacf <- function(x, lag.max = NULL, C = NULL,
                    penalized = TRUE, print = TRUE,
                    plot = TRUE, na.action = na.fail,
                    demean = TRUE, series = NULL, ...){
  if (is.null(series)){
    series <- deparse(substitute(x))
  }
  x <- na.action(as.ts(x))
  xfreq <- frequency(x)
  x <- as.matrix(x)
  if(!is.numeric(x))
    stop("'x' must be numeric")
  n <- as.integer(nrow(x))
  nser <- as.integer(ncol(x))
  if (is.na(n) || is.na(nser))
    stop("'n' and 'nser' must be integer")
  lag.max <- if (is.null(lag.max))
    floor(10 * (log10(n) - log10(nser))) else round(lag.max)
  lag.max <- as.integer(min(lag.max, n - 1L))
  if (is.na(lag.max) || lag.max < 1L)
    stop("'lag.max' must be >= 1")
  m <- as.integer(min(lag.max, floor(sqrt(n))))
  if(demean){
    xm <- colMeans(x, na.rm = TRUE)
    x <- sweep(x, 2L, colMeans(x, na.rm = TRUE), check.margin = FALSE)
  }

  if(penalized){
    if(anyNA(x)) stop("'NA in 'x'")
    pacf <- array(NA,dim = c(lag.max, 1L, nser))
    for (j in 1L:nser){
      k <- 1L : lag.max
      rhat <- stats::pacf(x[ ,j], lag.max = lag.max, plot = FALSE, na.action = na.action)$acf[1L:lag.max]
      rhat <- rhat[k]
      r <- abs(rhat)
      if (is.null(C)){
        l <-  sqrt(log(n)/n) * sqrt((m - 1 + k)/m)
      }else{
        if(C <= 0){
          stop("C must be positive real number")
        }else{
          l <-  (C/sqrt(n)) * sqrt((m - 1 + k)/m)
        }
      }
      bias <- (k + 2) * sign(rhat) * (r - l)/n +
        1/n * (rhat + 1 + floor(k/2) - floor((k - 1)/2))
      rtar <- (1 + (r - l)/l) * rhat * (sign(l - r) + 1)/2 +
        (rhat + bias) * (sign(r - l) + 1)/2
      rtar <- sign(rtar) * pmin(abs(rtar),(1 - 1/n))
      lambda <- l/r^2 * (sign(l - r) + 1)/2 +
        sqrt(n) * k * (r - l) * (1 - l)/(1-r) * (sign(r - l) + 1)/2
      w <- lambda/(1 + lambda)
      pi.tilde <- w*rtar + (1 - w)*rhat
      tmp <- pi.tilde
      pacf[ , ,j] <- c(tmp)
    }
    lag <- matrix(1L, 1L, nser)
    lag[lower.tri(lag)] <- -1
    lag <- outer(1L:lag.max, lag/xfreq)

    acf.out <- structure(.Data = list(acf = pacf, type = "penpartial", n.used = n,
                                      lag = lag, series = series, snames = colnames(x)))
    if(print){
      if (plot) {
        pencorrplot(acf.out, ...)
        pencorrprint(acf.out, ...)
      } else pencorrprint(acf.out, ...)
    }else{
      if (plot) {
        pencorrplot(acf.out, ...)
        invisible(acf.out)
      } else invisible(acf.out)
    }
  }else{
    tmp <- stats::pacf(x, ...)

    if(print){
      tmp
    }else{
      invisible(tmp)
    }
  }
}
