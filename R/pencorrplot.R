######
#' @title pencorrplot
#'
#' @description Plot method for "Penacf" and "Penpacf"
#'
#' @param x an object of "acf.out" or "pacf.out".
#' @param lag.max maximum lag at which to calculate the coefficients. Defaults to the smaller of \eqn{N-1} and \eqn{10log_{10}(N/nser)} where N is the number of non-missing observations and nser is the number of series.
#' @param ... additional arguments for penalized PACF estimation.
#'
#' @return plot of the penalized acf or pacf estimation
#'
#'
#' @examples
#' \dontrun{
#' data <- arima.sim(n=100, model=list(ar=0.5))
#'
#' pencorrplot(penacf(data))
#'
#' }
#' @export
#####

pencorrplot <- function (x, ci = 0.95, ci.type = c("white", "ma"),
                         type = "h", xlab = "Lag", ylab = NULL,
                         ylim = NULL, main = NULL, ci.col="blue",
                         ...){
  nser <- dim(x$lag)[3]
  snames <- x$snames
  if(nser < 1L)
    stop("x$lag must be >= 1")
  if (is.null(snames))
    snames <- paste("Series ", if (nser == 1L) x$series else 1L:nser)
  if (is.null(ylab))
    ylab <- switch(x$type,
                   pencorrelation = "Penalized ACF",
                   penpartial = "Penalized Partial ACF")

  ci.type <- match.arg(ci.type)
  with.ci <- ci > 0
  with.ci.ma <- with.ci && ci.type == "ma" && x$type == "pencorrelation"
  if (with.ci) clim0 <- qnorm((1 + ci)/2)/sqrt(x$n.used) else clim0 <- c(0, 0)
  if(with.ci.ma && x$lag[1L, 1L, 1L] != 0L) {
    warning("if first lag is 0, only can use ci.type=\"ma\"")
    with.ci.ma <- FALSE
  }

  for(i in 1:nser){
    if (is.null(ylim)) {
      ylim <- range(x$acf[, 1L, i], na.rm = TRUE)
      if (with.ci) ylim <- range(c(-clim0, clim0, ylim))
      if (with.ci.ma) {
        clim <- clim0 * sqrt(cumsum(c(1, 2*x$acf[-1, 1L, i]^2)))
        ylim <- range(c(-clim, clim, ylim))
      }else{
        clim <- if(x$type == "pencorrelation") rep(clim0, length(x$acf[, 1L, i])) else rep(clim0, length(x$acf[, 1L, i])+1)
      }
    }

    if(nser > 1){
      par(ask=TRUE)
    }else{
      par(ask=FALSE)
    }
    plot(x$lag[, 1L, i], x$acf[, 1L, i], type = type, xlab = xlab,
         ylab = ylab, ylim = ylim, ...)
    abline(h = 0)
    if (with.ci && ci.type == "white"){
      abline(h = c(clim, -clim), col = ci.col, lty = 2)
    } else if (with.ci.ma){
      lines(x$lag[, 1L, i], clim, col = ci.col, lty=2)
      lines(x$lag[, 1L, i], -clim, col = ci.col, lty=2)
    }
    title(if (!is.null(main)) main else snames[i], line = 2)
    if(nser > 1){
      xmax1 <- max(x$lag[, 1L, i]) - 1
      mtext(paste0("Page", i, "(", nser,")"), side=1, line=3.5, at=xmax1)
    }
  }
  invisible()
}
