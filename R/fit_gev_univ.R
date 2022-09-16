
# params is a vector of length 4
nll_univ <- function(params,
                     data, type, rel_trend = TRUE){

  if(!(type %in% c("scale", "shift", "stationary"))) { stop("Type must be one of 'shift',
                                                            'scale' or 'stationary'.")}
  n.dat <- nrow(data)

  if(type == "scale") {
    if(rel_trend) {

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
    else {

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
#'
#' @param data A tibble containing values of the unique sliding BM along with the
#' corresponding value of the temporal covariate and the frequency of occurence of
#' the respective tupel. Can be obtained by applying \code{\link[slbm]{get_uniq_bm}}.
#' @param method The method used during optimisation; passed to optim.
#' @param maxiter Passed to optim.
#' @param hessian logical; whether to return the hessian matrix.
#' @param return_cov logical; whether or not to return an estimate of the covariance matrix. If TRUE,
#' some further parameters/data need to be passed, see '...'
#' @param type One of   'stationary', 'scale' or 'shift', see details.
#' @param ... Several additional arguments which are only needed in some cases, e.g. in the scale model.
#' See `details` for further information.
#'
#' @return A list containing the parameter estimates, the value of the negative log-Likelihood,
#' the convergence code (ouput from \code{optim()}, 0 means everything was ok)
#' and the hessian (if \code{hessian = TRUE}).
#' @export
#'
#' @details # Details on type
#'
#' The argument given in type determines whether the observations shift or scale
#' with time. For a temporal covariate \eqn{(X_t)_t}, shifting corresponds to a shift
#' in the location parameter as follows:
#' \deqn{ \mu(t) = \mu + \alpha X_t , \sigma(t) = \sigma, \gamma(t) = \gamma}
#' while scaling corresponds to the model where
#' \deqn{ \mu(t) = \mu \exp(\alpha X_t /\mu),  \sigma(t) = \sigma  \exp(\alpha X_t /\mu),
#' \gamma(t) = \gamma,
#' }
#' as inspired by the Clausius-Clapeyron relation.
#' The latter model can also be parametrised as
#' \deqn{ \mu(t) = \mu \exp(\alpha X_t),  \sigma(t) = \sigma  \exp(\alpha X_t),
#' \gamma(t) = \gamma,
#' }
#' which is why one needs to specify the argument `rel_trend` whenever `type = scale`
#' @details # Detials on additional arguments
#' Additional arguments that need to be passed in some cases are
#'
#' * rel_trend : logical; specifies the parametrisation of the scale model,
#' i.e. only relevant when when 'type = scale'. When `TRUE`, the trend parameter \eqn{\alpha} is
#' seen relative to the location parameter \eqn{\mu}.
#'
#' When an estimate of the covariance matrix is required (`return_cov = TRUE`),
#' one needs to pass the following arguments and data:
#' * chain: logical; whether to use covariance matrix estimation based on chain rule
#' * varmeth: one of 'V', 'V2', 'both' (see documentation of \code{\link{est_var_univ}} for further details)
#' * blcksz The blocksize parameter.
#' * orig_slbm: the original sample of sliding block maxima (as generated by \code{\link{blockmax}})
#' * orig_cvrt: the temporal covariate for the original sliding block maxima.
#'
#' @examples
#' ##### generate some data #####
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
#' ### shift-fit with covariance matrix estimation:
#' fit_gev_univ(data = bms, hessian = TRUE, type = "shift", return_cov = TRUE,
#' varmeth = "V2", chain = TRUE, orig_slbm = blockmax(xx, 90, "sliding"),
#' orig_cvrt = temp_cvrt, blcksz = 90)

fit_gev_univ <- function(data, method = "BFGS", maxiter = 100,
                         hessian = FALSE, type, return_cov = FALSE, ...) {

  if(return_cov & !hessian) {
    hessian <- TRUE
  }
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
                 hessian = hessian, type = type, rel_trend = add.args$rel_trend)
  if(!(mlest$convergence == 0) ){print("Optimization didn't succeed.")}

  if(!return_cov) {
    if(!hessian) {
      return(list(mle = mlest$par, nll_expsc = mlest$value, conv = mlest$convergence,
                  fitted_model = paste(type, "model")))
    }
    else {
      return(list(mle = mlest$par, nll_expsc = mlest$value, conv = mlest$convergence,
                  hessian = mlest$hessian, fitted_model = paste(type, "model")))
    }
  }

  else {

    mllist <- list(mle = mlest$par, nll_expsc = mlest$value, conv = mlest$convergence,
                   hessian = mlest$hessian)

    varmeth <- add.args$varmeth
    chain <- add.args$chain

    orig_slbm <- add.args$orig_slbm
    orig_cvrt <- add.args$orig_cvrt

    blcksz <- add.args$blcksz

    if(is.null(varmeth) | is.null(chain) | is.null(blcksz)) {
      stop("Both variables 'varmeth' and 'chain' need to be
                specified for covariance matrix estimation.
           Further, the blocksize parameter is needed.")
    }
    if(is.null(orig_slbm)) {
      if( !(type == "stationary") & is.null(orig_cvrt)) {
        stop("Please pass the original sliding block maxima and the corresponding temporal covariate
             for covariance matrix estimation.")
      }
    }

    covhat <- slbm::est_var_univ(orig_slbm = orig_slbm, est_par = mllist, blcksz = blcksz,
                                 temp.cov = orig_cvrt, type = type, chain = chain, varmeth = varmeth,
                                 rel_trend = add.args$rel_trend)

    return(append(append(mllist, covhat), list(fitted_model = paste(type, "model"))))

  }

}

