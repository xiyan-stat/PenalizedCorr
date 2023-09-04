######
#' @title penacf
#'
#' @description The Penalized (Partial) Autocorrelation Function (ACF/PACF) Estimation
#'
#' @param x a univariate or multivariate numeric time series.
#' @param lag.max maximum lag at which to calculate the penalized/sample ACF/PACF. Defaults to the smaller of \eqn{N-1} and \eqn{10log_{10}(N/nser)} where N is the number of non-missing observations and nser is the number of series.
#' @param C a positive real number for l_h
#' @param type type of acf to be computed. Allow values are "correlation" (the default) or "partial".
#' @param penalized 'logical'. If 'TRUE' (the default) the penalized ACF/PACF is computed; if 'FALSE' the sample ACF/PACF is computed.
#' @param posdef the positive definite (two methods) penalized ACF is computed. (default is "origin") Allow values are "origin", "DL" or "NND".
#' @param print  'logical'. If 'TRUE' (the default) the penalized/sample ACF/PACF is printed.
#' @param plot  'logical'. If 'TRUE' (the default) the penalized/sample ACF/PACF is plotted.
#' @param series names for the series. Defaults to 'deparse(substitute(x))'.
#' @param na.action function to be called to handle missing values. 'na.pass' can be used.
#' @param demean 'logical'. Should a mean be estimated during estimating.
#' @param ... additional arguments for specific methods.
#'
#' @return An object of penalized/sample ACF/PACF estimation with the following elements:
#' \describe{
#' \item{\code{acf}}{An array containing the estimated penalized/sample ACF/PACF.}
#' \item{\code{lag}}{An array containing the lags at which the ACF/PACF is estimated.}
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
#' # Examples for penalized ACF/PACF and sample ACF/PACF
#' penacf(data)
#' penacf(data, penalized = FALSE)
#' penacf(data, type ="partial")
#' penacf(data, type ="partial", penalized = FALSE)
#'
#' x1 <- arima.sim(n=100, model=list(ar=0.5))
#' x2 <- arima.sim(n=100, model=list(ar=0.1))
#' x3 <- arima.sim(n=100, model=list(ar=0.9))
#' x <- cbind(x1, x2, x3)
#'
#' penacf(x)
#' penacf(x, penalized = FALSE)
#' penacf(x, type ="partial")
#' penacf(x, type ="partial", penalized = FALSE)
#' }
#' @export
#####


penacf <- function(x, lag.max = NULL, C = NULL,
                   type = c("correlation", "partial"),
                   penalized = TRUE,
                   posdef = c("origin","DL", "NND"),
                   print = TRUE, plot = TRUE, series = NULL,
                   na.action = na.fail, demean = TRUE, ...){
  type <- match.arg(type)
  if(type == "partial"){
    pacf.out <- penpacf(x, lag.max = lag.max,
                        penalized = penalized,
                        plot = FALSE, print = FALSE,
                        na.action = na.action, demean = demean, ...)
    if(print){
      if (plot) {
        pencorrplot(pacf.out, ...)
        pencorrprint(pacf.out, ...)
      } else pencorrprint(pacf.out, ...)
    }else{
      if (plot) {
        pencorrplot(pacf.out, ...)
        invisible(pacf.out)
      } else invisible(pacf.out)
    }
  }else{
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
    if (is.na(lag.max) || lag.max < 0)
      stop("'lag.max' must be >= 0")
    m <- as.integer(min(lag.max, floor(sqrt(n))))
    if(demean){
      xm <- colMeans(x, na.rm = TRUE)
      x <- sweep(x, 2L, colMeans(x, na.rm = TRUE), check.margin = FALSE)
    }

    if(penalized){
      if(anyNA(x)) stop("'NA in 'x'")
      acf <- array(NA, dim = c((lag.max+1), 1L, nser))
      for (j in 1:nser){
        k <- 1 : lag.max
        rhat <- stats::acf(x[ ,j], lag.max = lag.max, plot = FALSE, na.action = na.action)$acf[2:(lag.max+1)]
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
        rho.tilde <- w*rtar + (1 - w)*rhat
        tmp <- c(1, rho.tilde)

        posdef <- match.arg(posdef)

        if (posdef == "DL"){
          pacffit <- penpacf(x[ ,j], lag.max = lag.max, plot = FALSE, print=FALSE, ...)$acf
          if (lag.max == 1)
            rho.tilde <- pacffit
          rho.tilde <- 1:lag.max
          rho.tilde[1] <- pacffit[1]
          coefmat <- sapply((1:(lag.max-1)), DLpencoef, x = x[ ,j], ...)
          for (s in 1:(lag.max-1)){
            coef <- as.vector(coefmat[[s]])
            rho.tilde[(s+1)] <- cumprod((1 - pacffit^2))[s] * pacffit[(s+1)] + sum(coef %*% rho.tilde[s:1])
          }
          tmp <- c(1, rho.tilde)
        }
        if (posdef == "NND"){
          tmp <- NNDpenacf(x[ ,j], tmp, lag.max = lag.max, ...)$rhotilde.NND
        }
        acf[ , ,j] <- c(tmp)
      }

      lag <- matrix(1, 1, nser)
      lag[lower.tri(lag)] <- -1
      lag <- outer(0:lag.max, lag/xfreq)
      acf.out <- structure(list(acf = acf, type = "pencorrelation", n.used = n,
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
      tmp <- stats::acf(x, ...)
      if(print){
        tmp
      }else{
        invisible(tmp)
      }
    }
  }
}
