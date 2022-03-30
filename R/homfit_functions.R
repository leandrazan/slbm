## ML estimation under homogeneity constraint
nllh_hom <- function(param , x){
  d <- ncol(x)
  delta  <- param[1]
  sigmas <- param[2:(d+1)]
  names(sigmas) <- paste0("sigma_", 1:d)
  shape <- param[d+2]
  nx <- nrow(x)

  nvec <- nx - colSums(is.na(x))
  # sigma0 <-param$sigma0
  # sigma1 <- param
  vals <- numeric(d)
  if( abs(shape) < 1e-10){

    for( i in 1:d){
      y <-  x[, i]/sigmas[i] - delta
      if( sigmas[i] <= 0 ) { vals[i] <- 1e+08} else{
       vals[i] <- nvec[i]* log(sigmas[i]) + sum(y, na.rm = TRUE)+ sum(exp(-y), na.rm = TRUE)
      }
    }
  }
  else{
    for( i in 1:d){
      y <- 1+ shape*( x[, i]/sigmas[i] - delta )
      if( any(y < 0, na.rm = TRUE) | sigmas[i] <= 0 ){vals[i] <- 1e+08}
      else{
        vals[i] <- nvec[i]* log(sigmas[i]) +  (1/shape +1)*sum(log(y), na.rm = TRUE) + sum(y^(-1/shape), na.rm =TRUE)
      }
    }
  }
  return(sum(vals))
}

gr_homgev_link <- function( param, x, scale.link = make.link("log")){
  xgr <- x

  d <- ncol(xgr)
  delta  <- param[1]
  sigmas <- scale.link$linkinv(param[2:(d+1)])
  # names(sigmas) <- paste0("sigma_", 1:d)
  shape <- param[d+2]
  nx <- nrow(xgr)

  nvec <- nx - colSums(is.na(xgr))
  # sigma0 <-param$sigma0
  # sigma1 <- param

  gr_delta <- numeric(d)
  gr_scale <-  numeric(d)
  gr_shape <- numeric(d)

  if( abs(shape) < 1e-10){

    for( i in 1:d){
      y <-  xgr[, i]/sigmas[i] - delta
      gr_delta[i] <- sum(-1+exp(-y), na.rm = TRUE)
      gr_scale[i] <- sum(1/sigmas[i] - xgr[, i]/sigmas[i]^2 + exp(-y)/sigmas[i]^2 , na.rm = TRUE)

    }
    gr_delta <- sum(gr_delta)
    gr_shape <- 0
  }
  else{
    for( i in 1:d){
      y <- 1+ shape*( xgr[, i]/sigmas[i] - delta )
      if( any(y < 0, na.rm = TRUE)){gr_shape[i] <- Inf}
      else{
        gr_delta[i] <- sum(- shape*(1/shape +1)/y + y^(-1/shape -1), na.rm = TRUE)
        gr_shape[i] <- sum( -1/shape^2 *log(y) + 1/y*(y-1)/shape*(1/shape+ 1) +
                              ( 1/shape^2 *log(y) -1/shape*( (y-1)/shape/y)  ) *y^(-1/shape), na.rm = TRUE)
        gr_scale[i] <- sum( 1/sigmas[i] - shape*xgr[, i]/sigmas[i]^2*(1/shape +1)/y + y^(-1/shape -1)*xgr[,i]/sigmas[i]^2 ,na.rm = TRUE)
        # vals[i] <- nvec[i]* log(sigmas[i]) +  (1/shape +1)*sum(log(y), na.rm = TRUE) + sum(y^(-1/shape), na.rm =TRUE)
      }
    }
    gr_delta <- sum(gr_delta)
    gr_shape <- sum(gr_shape)
  }
  return(c(gr_delta, gr_scale, gr_shape))

}

nllh_hom_link <- function(param , x, scale.link = make.link("log")){
  d <- ncol(x)
  delta  <- param[1]
  sigmas <- param[2:(d+1)]
  names(sigmas) <- paste0("sigma_", 1:d)
  shape <- param[d+2]
  nx <- nrow(x)

  nvec <- nx - colSums(is.na(x))
  # sigma0 <-param$sigma0
  # sigma1 <- param
  vals <- numeric(d)
  if( abs(shape) < 1e-10){

    for( i in 1:d){
      y <-  x[, i]/scale.link$linkinv(sigmas[i]) - delta
      # map2( list(apply(x,2, list),  sigmas),
      #            .f = function(.x, .sigma){
      #             .x/.sigma - delta  })
      vals[i] <- nvec[i]* log(scale.link$linkinv(sigmas[i])) + sum(y, na.rm = TRUE)+ sum(exp(-y), na.rm = TRUE)
    }
  }
  else{
    for( i in 1:d){
      y <- 1+ shape*( x[, i]/scale.link$linkinv(sigmas[i]) - delta )
      if( any(y < 0, na.rm = TRUE)){vals[i] <- Inf}
      else{
        vals[i] <- nvec[i]* log(scale.link$linkinv(sigmas[i])) +  (1/shape +1)*sum(log(y), na.rm = TRUE) + sum(y^(-1/shape), na.rm =TRUE)
      }
    }
  }
  return(sum(vals))
}


#' Maximum-Likelihood fitting of the GEV distribution under homogeneity constraint
#'
#' @description Maximum-Likelihood fitting of the GEV distribution under homogeneity
#' constraint, also known as the 'Indexflood'-assumption.
#' For (approximately) GEV distributed observations from several locations in space
#' it is assumed that the shape parameter as well as the ratio
#' of location and scale parameters (the \strong{dispersion} parameter) is constant
#' in space. To be precise, let
#' \eqn{X_{i,j}, i = 1, \ldots, n; j = 1, \ldots,d } denote the observations
#' from \eqn{d} stations. It is assumed that
#' \eqn{ X_{i,j} \sim GEV(\mu(j), \sigma(j), \gamma)},
#' where \eqn{ \mu(j) / \sigma(j) = \delta }, which does not depend on \eqn{d}.
#'
#' @param dat A matrix or data frame. Each column contains values for one location.
#' May also contain missing values.
#' @param start_vals Optional: a vector containing initial values for the log Likelihood
#' which is to be optimised. If omitted, initial values are estimated.
#' See below for further details.
#' @param scale.link The link function that is used for modelling the scale parameter.
#' Only used during optimisation.
#' @param method Optimisation method that is passed to \code{\link[stats]{optim}}.
#' @param ... Further arguments that are passed to the `control` argument of \code{\link[stats]{optim}}.
#'
#' @return Returns a list with the following components:
#'  \describe{
#'   \item{mle}{A \eqn{3 x d} dimensional matrix containing the estimated location,
#'   scale and shape parameter at each station.}
#'   \item{mledisp}{The estimated dispersion parameter, i.e. the ratio of the estimated
#'   location and scale parameters.}
#'   \item{conv}{The convergence code as obtained from \code{optim}.}
#'   \item{counts}{A two-element integer vector giving the number of calls to
#'   the negative log-Likelihood and its gradient, respectively, as obtained from
#'   \code{optim}.}
#'   \item{nllh}{Value of the minimised negative log-Likelihood.}
#'   }
#' @export
#'
#' @examples
#' xx <- sapply(1:5, function(x) evd::rgev(30, loc = 2*x, scale = x, shape = 0.2))
#' fgev_hom(xx)
#'
#' @details If no starting values are provided, the initial values are computed as
#' follows: First, GEV parameters are estimated stationwise.
#' The initial dispersion parameter
#' is set to the mean of stationwise dispersion parameters, the initial shape
#' parameter is set to the mean of stationwise shape parameters and the initial scale
#' parameters are the stationwise estimates.
#'
#' The scale parameters that are returned are the actual scale parameters.
#' E.g., when using \eqn{link(\sigma) = \tilde{\sigma}} as link function,
#'  \eqn{\sigma =  linkinv{ \tilde{\sigma}}} is returned.
#'
#'
fgev_hom <- function(dat, start_vals = NULL, scale.link = make.link("log"),
                          method = "BFGS", ...){
  if(!is.matrix(dat)) { dat <- as.matrix(dat) }

  add.args <- list(...)
  d <- ncol(dat)

  if(is.null(start_vals)) {
    start_st <- apply(dat, 2, function(x) evd::fgev(x, std.err = FALSE)$estimate)
    start_vals <- c(mean(start_st[1, ]/start_st[2, ]),
                    scale.link$linkfun(start_st[2, ]), mean(start_st[3, ]))
    #
    # start_st <-  evd::fgev(as.vector(dat), std.err = FALSE)$estimate
    # start_vals <- c(start_st[1]/start_st[2], rep(scale.link$linkfun(start_st[2]), d), start_st[3])

  }
  names(start_vals) <- c("delta", paste0("sigma_", 1:d), "shape")


  if(scale.link$name == "identity") {
    mlest <- optim(start_vals, fn = nllh_hom,
                   x = dat, method = method, control = add.args)
  } else {
   mlest <- optim(start_vals, fn = nllh_hom_link, gr = gr_homgev_link,
                   x = dat, method = method, control = add.args)
  }
  if(!(mlest$convergence == 0) ){warning("Optimization didn't succeed.")}
  mle <- mlest$par
  mle <- matrix(c(mle[1]*scale.link$linkinv( mle[2:(d+1)]), scale.link$linkinv(mle[2:(d+1)]),
                    rep(mle[d+2], d)),
                  byrow = TRUE, nrow = 3)
  rownames(mle) <- c("loc", "scale", "shape")
  return(list(mle = mle, mledisp = mle[1,1]/mle[2,1],
              conv = mlest$convergence, counts = mlest$counts, nllh = mlest$value))
}

