######
#' @title ar.penyw.default
#'
#' @description Fit an autoregressive time series model to the data using penalized correlation estimation, by default selecting the order by AIC.
#'
#' @param x a univariate numeric time series.
#' @param order.max maximum lag at which to calculate the penalized/sample PACF. Defaults to the smaller of \eqn{N-1} and \eqn{10log_{10}(N/nser)} where N is the number of non-missing observations and nser is the number of series.
#' @param penalized 'logical'. If 'TRUE' (the default) the penalized autoregressive model is fitted; if 'FALSE' the \code{ar} is fitted and \code{method} can be selected.
#' @param aic 'logical', If 'TRUE' then the Akaike Information Criterion is used to choose the order of the autoregressive model. If FALSE, the model of order order.max is fitted.
#' @param na.action function to be called to handle missing values. 'na.pass' can be used.
#' @param ... additional arguments for \code{ar} function.
#'
#' @return An object of penalized ar model fit with the following elements:
#' \describe{
#' \item{\code{order}}{The order of the fitted model. This is chosen by minimizing the AIC if 'AIC=TRUE', otherwise it is 'order.max'.}
#' \item{\code{ar}}{Estimated penalized autorregression coefficients for the fitted model.}
#' \item{\code{aic}}{The vector of AIC values.}
#' \item{\code{n.used}}{The number of observations in the time series.}
#' \item{\code{order.max}}{The value of the 'order.max' argument.}
#' \item{\code{call}}{The matched call.}
#' }
#'
#'
#' @examples
#' \dontrun{
#' data <- arima.sim(n=100, model=list(ar=0.5))
#'
#' # Example for penalized ar model fit
#' penar(data)
#' }
#' @export
#####

ar.penyw.default <- function(x, order.max = NULL, aic = TRUE,
                  na.action = na.fail, demean = TRUE,
                  series = NULL, ...){
  if(is.null(series))
    series <- deparse(substitute(x))
  x <- na.action(as.ts(x))
  if(is.ts(x))
    xtsp <- tsp(x)
  xfreq <- frequency(x)
  x <- as.matrix(x)
  if(!is.numeric(x))
    stop("'x' must be numeric")

  nser <- as.integer(ncol(x))
  n <- as.integer(nrow(x))
  if(is.na(n) | is.na(nser))
    stop("'n' and 'nser' must be integer")

  if(demean){
    xm <- colMeans(x)
    x <- sweep(x, 2L, xm, check.margin=FALSE)
  }else{
    xm <- rep.int(0, n)
  }
  if(is.null(order.max))
    order.max <- floor(10 * (log10(n)))
  order.max <- as.integer(min(order.max, n - 1L))
  if(is.na(order.max) || order.max < 1L)
    stop("'order.max' must be >= 1")
  else if (order.max >= n)
    stop("'order.max' must be less than 'n'")
  xacf0 <- acf(x, type = "covariance", lag.max = order.max, plot = FALSE,
              demean = demean)$acf[1]
  penacf <- penacf(x, lag.max = order.max, plot = FALSE, print = FALSE)$acf
  xacf <- penacf * xacf0

  if(nser > 1)
    stop("'x' must be multivariate")

  if(nser == 1){
    if (xacf[1L] == 0) stop("zero-variance series")
    ppacf <- penpacf(x, lag.max = order.max, plot = FALSE, print = FALSE, ...)$acf
    #univariate case
    r <- as.double(drop(xacf))
    v <- r[1L]
    d <- r[2L]
    a[1L] <- 1
    f[1L, 1L] <- r[2L] / r[1L]
    q <- f[1L, 1L] * r[2L]
    vars[1L] <- (1 - f[1L, 1L]^2) * r[1L]
    if (order.max == 1)
      pred.var = vars[1L]
    for (l in 2L:order.max) {
      a[l] <- -d / v
      if (l > 2) {
        l1 <- (l - 2) / 2
        l2 <- l1 + 1
        for (j in 2L:l2) {
          hold <- a[j]
          k <- l - j + 1
          a[j] <- a[j] + a[l] * a[k]
          a[k] <- a[k] + a[l] * hold
        }
        if (2 * l1 != l - 2) a[l2 + 1] <- a[l2 + 1] * (1 + a[l])
      }
      v <- v + a[l] * d
      f[l, l] <- (g[l + 1] - q) / v
      vars[l] <- vars[l - 1] * (1 - f[l, l]^2)
      d <- 0
      q <- 0
      for (i in 1:l) {
        k <- l - i + 2
        d <- d + a[i] * r[k]
        q <- q + f[l, i] * r[k]
      }
    }
    partialacf <- ppacf
    var.pred <- c(r[1L], vars)
    AICpen <- sapply(1:order.max, getaic, x=x)
    aic0 <- n * log(var(x))
    xaic <- c(aic0, AICpen)
    #xaic <- n * log(var.pred) + 2 * (0L:order.max) + 2 * demean
    maic <- min(aic)
    xaic <- setNames(if(is.finite(maic)) xaic - min(xaic) else
      ifelse(xaic == maic, 0, Inf), 0L:order.max)
    order <- if (aic) (0L:order.max)[xaic == 0L] else order.max
    ar <- if (order) DLpencoef(x, lag.max = order, ...) else numeric()
    var.pred <- var.pred[order + 1L]
    var.pred <- var.pred * n/(n - (order + 1L))
    resid <- if(order) c(rep.int(NA, order), embed(x, order + 1L) %*% c(1, -ar))
    else as.vector(x)
    if(is.ts(x)) {
      attr(resid, "tsp") <- xtsp
      attr(resid, "class") <- "ts"
    }
  }

  res <- list(order = order, ar = ar, var.pred = var.pred, x.mean  =  drop(xm),
              aic  =  xaic, n.used = n, order.max = order.max,
              partialacf = partialacf, resid = resid, method = "Penalized Yule-Walker",
              series = series, frequency = xfreq, call = match.call())
  if(nser == 1L && order)
    res$asy.var.coef <-
    solve(toeplitz(drop(xacf)[seq_len(order)]))*var.pred/n
  class(res) <- "ar"
  res
}
