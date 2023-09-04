# R Package: PenalizedCorr
# Penalized Correlation Function Estimation for Time Series Analysis
#
# This package includes 11 functions to estimate penalized/sample acf/pacf
# and fit AR model using the penalized estimators.
#
# Some examples are included in the help file for each function.

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

DLpencoef <- function(x, lag.max = NULL, ...){
  n <- as.integer(nrow(x))
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

penar <- function(x, order.max = NULL, AIC = TRUE, na.action = na.fail, arprint = TRUE, ...){
  x <- na.action(as.ts(x))
  x <- as.matrix(x)
  if(!is.numeric(x))
    stop("'x' must be numeric")
  n <- as.integer(nrow(x))
  nser <- as.integer(ncol(x))
  if (is.na(n) | is.na(nser))
    stop("'n' and 'nser' must be integer")
  if(!AIC){
    coef <- DLpencoef(x, lag.max = order.max, ...)
    res <- list(order = oder.max, coef = coef)
    res
  }
  if (is.null(order.max))
    order.max <- floor(10 * (log10(n)))
  order.max <- as.integer(min(order.max, n - 1L))
  if (is.na(order.max) || order.max < 1L)
    stop("'order.max' must be >= 1")
  else if (order.max >= n)
    stop("'order.max' must be less than 'n'")
  AICpen <- sapply(1:order.max, getAIC, x=x)
  if (nser == 1L)
    AICpen <- t(as.matrix(AICpen))
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
  res <- list(order = penorder, penar = pencoef, AIC = AICpen,
              n.used = n, order.max = order.max, call = match.call())
  if (arprint){
    penarprint(res)
  }else res
}

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

penarprint <- function(x, digits = 3L, ...){
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  cat("Coefficients:\n")
  print(x$penar, digits = digits)
  cat("\nOrder selected:", x$order)
  cat("\n")
  invisible(x)
}

#penarpredict function to fit the autoregressive model from penar
penarpredict <- function(x, order.max = NULL, ...){

}


