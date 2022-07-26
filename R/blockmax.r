#' Compute Disjoint or Sliding Block Maxima of a Univariate Time Series
#'
#' @description For a given univariate time series of length n and block size
#' parameter r (\eqn{1 \le r \le n}), `blockmax` computes either the time
#' series of disjoint (non-overlapping) or sliding (overlapping) block maxima.
#'
#' @param data a vector of length n representing the time series
#' @param r the block length parameter
#' @param method either `disjoint` or `sliding`
#' @param looplastblock logical; applies only when method is 'sliding'. If TRUE (not the default)
#' the first disjoint block of observations is concatenated to the last disjoint block of observations
#' in order to give each observation the same
#' chance to appear repeatedly in the sliding block maxima sample.
#' @details If the length of `data` is not a multiple of the block length `r`,
#' the last block maximum returned when using the `disjoint`-method is the maximum
#' of the remaining \eqn{(\lfloor n/r\rfloor) +1) \cdot r -n)} observations, where \eqn{n}
#' is the length of the data.
#'
#' @return The time series of disjoint block maxima
#' (a vector of length \eqn{n/r} or \eqn{\lfloor n/r \rfloor} +1, see details)
#' or sliding block maxima
#' (a vector of length \eqn{n-r+1})
#' @export
#'
#' @examples
#' x <- runif(200)
#' blockmax(x, 20, "disjoint")
#' blockmax(x, 20, "sliding")
blockmax <- function(data, r, method=c("disjoint", "sliding"), looplastblock = FALSE) {
  method <- match.arg(method)
  if (!is.vector(data, mode = "numeric"))
    stop("data must be a numeric vector")
  n <- length(data)
  if ( !(assertthat::is.count(r)) || r > n )
    stop("r must be a whole number between 1 and length(data)")
  if(method == "disjoint") {
    result <- apply(X = matrix(c(data,
                                 rep(NA, ifelse(((floor(n/r) +1)*r -n) == r , 0,
                                             (floor(n/r) +1)*r -n))), nrow = r),
                                 MARGIN = 2,
                                 FUN = max, na.rm = TRUE)
  } else {
    if(looplastblock) {
      result <- RcppRoll::roll_max(c(data, data[1:(r-1)]), r, na.rm = TRUE)
    }
    else {
      result <- RcppRoll::roll_max(data, r, na.rm = TRUE)
    }
  }

  result[result == -Inf] <- NA
  return(result)
}

