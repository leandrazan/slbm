#' Check function for evaluating an estimated quantile
#'
#' @param u Difference of observed maximum in a test set and estimated p*100% quantile
#' @param probability
#'
#' @return value of check function
#' @export
#'
#' @examples
check_fun <- function( u, p){
  if( p < 0 | p > 1){stop("p is a probability (i.e. between 0 and 1")}
  ifelse( u  > 0 , p*u, (p-1)*u)
}

#' Compute the quantile score for an estimated quantile
#'
#' @param obs vector of observations from a test set
#' @param quant.est value of estimated p*100% quantile
#' @param p probability
#'
#' @return a numeric value giving the quantile score
#' @export
#'
#' @examples set.seed(1)
#'           train <- evd::rgev(100)
#'           estpars <- evd::fgev(train[1:97])$estimate
#'           quanthat <- evd::qgev(0.99, loc = estpars[1], scale = estpars[2],
#'            shape = estpars[3])
#'            compute_qs(train[98:100], quanthat, 0.99 )
#'
compute_qs <- function( obs, quant.est, p  ){
  mean(sapply( obs - quant.est, check_fun, p = p ))
}



#' Compute the Quantile Skill Index (QSI) of new model vs reference model
#'
#' @param qs.new.model quantile score of the new model
#' @param qs.ref.model quantile score of the reference model
#'
#' @return numeric value between -1 and 1. Values smaller than 0 indicate that
#' the reference model performs better, while values larger than 0 indicate that
#' the new model performs better
#' @export
#'
#' @examples
#' set.seed(1)
#' train <- evd::rgev(500)
#' test <- evd::rgev(30)
#' estpars1 <- evd::fgev(train)$estimate
#' quanthat1 <- evd::qgev(0.99, loc = estpars1[1], scale = estpars1[2],
#'                      shape = estpars1[3])
#' estpars2 <- evd::fgev(train, loc = 0 )$estimate # fix location to true value
#' quanthat2 <- evd::qgev(0.99, loc = 0, scale = estpars2[1],
#'                      shape = estpars2[2])
#' qs1 <- compute_qs(test, quanthat1, 0.99 )
#' qs2 <- compute_qs(test, quanthat2, 0.99 )
#' compute_qsi(qs2, qs1)

compute_qsi <- function(qs.new.model, qs.ref.model){
  if(anyNA(c(qs.new.model, qs.ref.model))){return(NA)}
  if( qs.new.model < qs.ref.model) {
    qsi <- 1 - qs.new.model/qs.ref.model
  }
  else{
    qsi <- qs.ref.model/qs.new.model - 1
  }
  return(qsi)
}

