######
#' @title penar
#'
#' @description Fit an autoregressive time series model to the data using penalized correlation estimation, by default selecting the order by AIC.
#'
#' @param x a univariate numeric time series.
#' @param order.max maximum lag at which to calculate the penalized/sample PACF. Defaults to the smaller of \eqn{N-1} and \eqn{10log_{10}(N/nser)} where N is the number of non-missing observations and nser is the number of series.
#' @param penalized 'logical'. If 'TRUE' (the default) the penalized autoregressive model is fitted; if 'FALSE' the \code{ar} is fitted and \code{method} can be selected.
#' @param AIC 'logical', If 'TRUE' then the Akaike Information Criterion is used to choose the order of the autoregressive model. If FALSE, the model of order order.max is fitted.
#' @param na.action function to be called to handle missing values. 'na.pass' can be used.
#' @param arprint  'logical'. If 'TRUE' (the default) the fitted penalized ar model is printed.
#' @param method character string specifying the method to fit the regular autoregressive time series model. Default to "yule-walker".
#' @param ... additional arguments for \code{ar} function.
#'
#' @return An object of penalized ar model fit with the following elements:
#' \describe{
#' \item{\code{order}}{The order of the fitted model. This is chosen by minimizing the AIC if 'AIC=TRUE', otherwise it is 'order.max'.}
#' \item{\code{penar}}{Estimated penalized autorregression coefficients for the fitted model.}
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

penar <- function(x, order.max = NULL,
                  penalized = TRUE, AIC = TRUE,
                  na.action = na.fail, arprint = TRUE,
                  method = c("yule-walker","burg", "ols", "mle", "yw"),
                  series = deparse(substitute(x)), ...){
  if(!penalized){
    ar.res <- stats::ar(x, order.max = order.max, aic = AIC, na.action = na.action, method = method, series = series,...)
    return(ar.res)
  }

  x <- na.action(as.ts(x))
  x <- as.matrix(x)
  if(!is.numeric(x))
    stop("'x' must be numeric")
  n <- as.integer(nrow(x))
  nser <- as.integer(ncol(x))
  if(is.na(n) | is.na(nser))
    stop("'n' and 'nser' must be integer")

  if(is.null(order.max))
    order.max <- floor(10 * (log10(n)))
  order.max <- as.integer(min(order.max, n - 1L))
  if(is.na(order.max) || order.max < 1L)
    stop("'order.max' must be >= 1")
  else if (order.max >= n)
    stop("'order.max' must be less than 'n'")

  if(is.null(series)) series <- deparse(substitute(x))

  if(!AIC){
    pencoef <- DLpencoef(x, lag.max = order.max, ...)
    penorder <- order.max
    AICpen <- NULL
  }

  AICpen <- sapply(1:order.max, getaic, x=x)
  if(nser == 1L){
    AICpen <- t(as.matrix(AICpen))
  }else{
    stop("'x' must be univariate")
  }
  penorder <- NULL
  pencoef <- list()
  for (j in 1L:nser){
    tmp.order <- order(AICpen[j, ])[1]
    tmp.coef <- DLpencoef(x[1L : n, j], lag.max = tmp.order, ...)
    aic0 <- n * log(var(x[, j]))
    if(aic0 < min(AICpen[j, ])){
      tmp.order = 0
      tmp.coef = 0
    }
    penorder = c(penorder, tmp.order)
    pencoef[[j]] =  tmp.coef
  }
  res <- list(order = penorder, ar = pencoef, aic = AICpen,
              n.used = n, order.max = order.max, call = match.call())
  if (arprint){
    penarprint(res)
  }else return(res)
}
