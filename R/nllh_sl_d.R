#' Negative log-Likelihood for sliding BM of ID process
#'
#' @param parvec Vector of parameter values
#' @param bm_uni Dataframe containing the unqiue values of sliding block maxima along with
#'  their frequency of appearance, for possibly several values of the duration $d$, as
#'  specified in column 'duration'.
#' @param sigma0link The link function for \eqn{\sigma_0}. Default is `make.link("log")`.
#' @param mult_sc Whether or not you pass a multi-scaling parameter.
#' @param dur_offset Whether or not you pass a duration offset parameter.
#' @param int_offset Whether or not you pass an intensity offset parameter.
#'
#' @return The value of the negative log-Likelihood at the given parameters.
#' @export
#'
#' @examples
nllh.sl.d <- function(parvec, bm_uni, sigma0link = make.link("log"), mult_sc = TRUE,
                                  dur_offset = FALSE , int_offset = FALSE){

  mut <- parvec[["mut"]]
  sigma0 <-  sigma0link$linkinv(parvec[["sigma0"]])

  xi <- parvec[["shape"]]

  if(!(dur_offset) ) {
    theta <- 0
  } else{
    theta <- parvec[["theta"]]
  }
  eta <- parvec[["eta"]]
  if(mult_sc){
    eta2 <- parvec[["eta2"]]
  }  else{
    eta2 <- 0
  }
  if(int_offset){
    tau <- parvec[["tau"]]
  } else{
    tau <- 0
  }

  d <- bm_uni$duration

  dpt <- d + theta
  sigma.d <- sigma0/(dpt^(eta + eta2)) + tau
  mu.d <- mut*(sigma0*(d + theta)^(-eta) + tau)
  y <- (bm_uni$sldata - mu.d)/sigma.d
  y <- 1 + xi * y

  if (dur_offset) {
    if ( eta <= 0 ||  theta < 0  || any(sigma.d <= 0) || any(y <= 0) || theta > 1 )
      return(10^6)
  } else {
    if ( eta + eta2 <= 0  || any(sigma.d <= 0) || any(y <=  0)  ||  tau < 0 ){
      return(10^6)
    }
  }


  nllh.value <- 1/sum(bm_uni$n)*
    as.numeric(t( bm_uni$n) %*% ( log(sigma.d) + y^(-1/xi) + log(y) * (1/xi +
                                                                         1)))


  return( nllh.value )

}
