#' Compute Disjoint or Sliding Block Maxima of a Univariate Time Series
#'
#' @description For a given univariate time series of length n and block size
#' parameter r (\eqn{1 \le r \le n}), `blockmax` computes either the time
#' series of disjoint (non-overlapping) or sliding (overlapping) block maxima.
#'
#' @param data a vector of length n representing the time series
#' @param r the block length parameter
#' @param method either `disjoint` or `sliding`
#'
#' @return The time series of disjoint block maxima
#' (a vector of length \eqn{\lfloor n/r \rfloor})
#' or sliding block maxima
#' (a vector of length \eqn{n-r+1})
#' @export
#'
#' @examples
#' x <- runif(200)
#' blockmax(x, 20, "disjoint")
#' blockmax(x, 20, "sliding")
blockmax <- function(data, r, method=c("disjoint", "sliding")) {
  method <- match.arg(method)
  if (!is.vector(data, mode = "numeric"))
    stop("data must be a numeric vector")
  n <- length(data)
  if ( !(assertthat::is.count(r)) || r > n )
    stop("r must be a whole number between 1 and length(data)")
  result <- switch(method,
                   "disjoint" = apply(X = matrix(c(data,
                                rep(NA, ifelse(((floor(n/r) +1)*r -n) == r , 0,
                                             (floor(n/r) +1)*r -n))), nrow = r),
                                      MARGIN = 2,
                                      FUN = max, na.rm = TRUE),
                   "sliding" = RcppRoll::roll_max(data, r, na.rm = TRUE))

  result[result == -Inf] <- NA
  return(result)
}

