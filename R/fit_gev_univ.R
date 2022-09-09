
# params is a vector of length 4
nll_univ <- function(params,
                     data, type, rel_trend = TRUE){

  if(!(type %in% c("scale", "shift", "stationary"))) { stop("Type must be one of 'shift',
                                                            'scale' or 'stationary'.")}
  n.dat <- nrow(data)

  if(type == "scale" & !rel_trend) {

    mu0 <- c("mu" = params[1])
    sigma0 <- c("sigma" = params[2])
    gamma0 <- c("shape" = params[3])
    alpha0 <- c("alpha" = params[4])

    if(sigma0 <= 0) {return(1e+10)}
    else {
      mut <- mu0*exp(alpha0*data$temp_cvrt/mu0)
      sigmat <-  sigma0*exp(alpha0*data$temp_cvrt/mu0)


      if(abs(gamma0) < 1e-8){
        zt <- exp( -(data$slbm - mut)/sigmat)

        loglik <-  sum( data$n * (log(sigmat) + log(zt) + zt) , na.rm = TRUE)
      }
      else{
        zt <- 1 + gamma0*(data$slbm - mut)/sigmat
        if(any(zt < 0, na.rm = TRUE)){ loglik <- 1e+10}
        else{
          loglik <- sum(data$n*( log(sigmat) + (1/gamma0 +1)*log(zt) + zt^(-1/gamma0)),
                        na.rm = TRUE)
        }
      }
    }
  }
  if(type == "scale" & rel_trend) {

    mu0 <- c("mu" = params[1])
    sigma0 <- c("sigma" = params[2])
    gamma0 <- c("shape" = params[3])
    alpha0 <- c("alpha" = params[4])

    if(sigma0 <= 0) {return(1e+10)}
    else {
      mut <- mu0*exp(alpha0*data$temp_cvrt)
      sigmat <-  sigma0*exp(alpha0*data$temp_cvrt)


      if(abs(gamma0) < 1e-8){
        zt <- exp( -(data$slbm - mut)/sigmat)

        loglik <-  sum( data$n * (log(sigmat) + log(zt) + zt) , na.rm = TRUE)
      }
      else{
        zt <- 1 + gamma0*(data$slbm - mut)/sigmat
        if(any(zt < 0, na.rm = TRUE)){ loglik <- 1e+10}
        else{
          loglik <- sum(data$n*( log(sigmat) + (1/gamma0 +1)*log(zt) + zt^(-1/gamma0)),
                        na.rm = TRUE)
        }
      }
    }
  }
  if(type == "shift") {

    mu0 <- params[1]
    sigma <- params[2]
    xi <- params[3]
    mu1 <- params[4]

    if(sigma <= 0) {return(1e+10)}
    else {
      mut <- mu0 + mu1*data$temp_cvrt

      if(abs(xi) < 1e-8){
        zt <- exp( -(data$slbm - mut)/sigma)

        loglik <-  sum( data$n * (log(sigma) + log(zt) + zt) , na.rm = TRUE)
      }
      else{
        zt <- 1 + xi*(data$slbm - mut)/sigma
        if(any(zt < 0, na.rm = TRUE)){ loglik <- 1e+10}
        else{
          loglik <- sum(data$n*( log(sigma) + (1/xi +1)*log(zt) + zt^(-1/xi)),
                        na.rm = TRUE)
        }
      }
    }
  }
  if(type == "stationary") {

    mu <- params[1]
    sigma <- params[2]
    xi <- params[3]

    if(sigma <= 0) {return(1e+10)}
    else {

      if(abs(xi) < 1e-8){
        z <- exp( -(data$slbm - mu)/sigma)

        loglik <-  sum( data$n * (log(sigma) + log(z) + z) , na.rm = TRUE)
      }
      else{
        z <- 1 + xi*(data$slbm - mu)/sigma
        if(any(z < 0, na.rm = TRUE)){ loglik <- 1e+10}
        else{
          loglik <- sum(data$n*( log(sigma) + (1/xi +1)*log(z) + z^(-1/xi)),
                        na.rm = TRUE)
        }
      }
    }
  }
  loglik
}

#' GEV fit with trend
#' @description Fit a model that either shifts or scales with a temporal covariate
#' to univariate data
#' @param data A tibble containing values of the unique sliding BM along with the
#' corresponding value of the temporal covariate and the frequency of occurence of
#' the respective tupel. Can be obtained by applying \code{\link[slbm]{get_uniq_bm}}.
#' @param method The method used during optimisation; passed to optim.
#' @param maxiter Passed to optim.
#' @param hessian logical; whether to return the hessian matrix.
#' @param type One of 'scale' or 'shift', see details.
#' @return A list containing the parameter estimates, the value of the negative log-Likelihood,
#' the convergence code (ouput from \code{optim()}, 0 means everything was ok)
#' and the hessian (if \code{hessian = TRUE}).
#' @export
#'
#' @details The argument given in type determines whether the observations shift or scale
#' with time. For a temporal covariate \eqn{(X_t)_t}, shifting corresponds to a shift
#' in the location parameter as follows:
#' \deqn{ \mu(t) = \mu + \alpha X_t , \sigma(t) = \sigma, \gamma(t) = \gamma}
#' while scaling corresponds to the model where
#' \deqn{ \mu(t) = \mu \exp(\alpha X_t /\mu),  \sigma(t) = \sigma  \exp(\alpha X_t /\mu),
#' \gamma(t) = \gamma,
#' }
#' as inspired by the Clausius-Clapeyron relation.
#'
#' @examples
#' #' ##### generate some data #####
#' set.seed(1)
#' blcksz <- 90
#' xx <- evd::rgpd(100*90, shape = 0.2)
#'
#' # define a temporal covariate that is constant over a block of length blcksz
#' temp_cvrt <- rep(1:100/100, each = blcksz)[1:(99*blcksz + 1)]
#'
#' bms <- get_uniq_bm(xx, blcksz, temp_cvrt = temp_cvrt, looplastblock = FALSE)
#' bms
#'###############################
#' fit_gev_univ(data = bms, hessian = TRUE, type = "shift")
#' fit_gev_univ(data = bms, hessian = TRUE, type = "stationary")

fit_gev_univ <- function(data, method = "BFGS", maxiter = 100,
                         hessian = FALSE, type, ...) {

  add.args <- list(...)
  if(type == "scale" & is.null(add.args$rel_trend)) {stop("Please specify parametrization of scale model.")}

  start_st <-  evd::fgev(as.vector(data$slbm), std.err = FALSE)$estimate
  if(type == "stationary") {
    start_vals <- c(start_st[1], start_st[2], start_st[3])
    names(start_vals) <- c("loc", "scale", "shape")
  }
  else {
    start_vals <- c(start_st[1], start_st[2], start_st[3], 0)
    names(start_vals) <- c("mu", "sigma", "shape", "alpha")
  }
  #print(start_vals)
  mlest <- optim(start_vals, fn = nll_univ,
                 data = data,
                 method = method, control = list(maxit = maxiter),
                 hessian = hessian, type = type, ... = ...)
  if(!(mlest$convergence == 0) ){print("Optimization didn't succeed.")}

  if(!hessian) {
    return(list(mle = mlest$par, nll_expsc = mlest$value, conv = mlest$convergence))
  }
  else {
    return(list(mle = mlest$par, nll_expsc = mlest$value, conv = mlest$convergence,
                hessian = mlest$hessian))

  }

}

