nll_univ_dj <- function(params,
                        data, temp_cvrt, type, rel_trend = FALSE){

  if(!(type %in% c("scale", "shift", "stationary"))) { stop("Type must be one of 'shift',
                                                            'scale' or 'stationary'.")}
  n.dat <- length(data)

  if(type == "scale") {

    mu <- c("mu" = params[1])
    sigma <- c("sigma" = params[2])
    shape <- c("shape" = params[3])
    alpha <- c("alpha" = params[4])

    if(sigma <= 0) {return(1e+10)}
    else {
      if(rel_trend) {
        mut <- mu*exp(alpha*temp_cvrt/mu)
        sigmat <-  sigma*exp(alpha*temp_cvrt/mu)
      }
      else {
        mut <- mu*exp(alpha*temp_cvrt)
        sigmat <-  sigma*exp(alpha*temp_cvrt)
      }

      if(abs(shape) < 1e-8){
        zt <- exp( -(data - mut)/sigmat)

        loglik <- mean( (log(sigmat) + log(zt) + zt) , na.rm = TRUE)
      }
      else{
        zt <- 1 + shape*(data - mut)/sigmat
        if(any(zt < 0, na.rm = TRUE)){ loglik <- 1e+10}
        else{
          loglik <- mean(( log(sigmat) + (1/shape +1)*log(zt) + zt^(-1/shape)),
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
      mut <- mu0 + mu1*temp_cvrt

      if(abs(xi) < 1e-8){
        zt <- exp( -(data - mut)/sigma)

        loglik <-  mean(  (log(sigma) + log(zt) + zt) , na.rm = TRUE)
      }
      else{
        zt <- 1 + xi*(data - mut)/sigma
        if(any(zt < 0, na.rm = TRUE)){ loglik <- 1e+10}
        else{
          loglik <- mean(( log(sigma) + (1/xi +1)*log(zt) + zt^(-1/xi)),
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
        z <- exp( -(data - mu)/sigma)

        loglik <-  mean((log(sigma) + log(z) + z) , na.rm = TRUE)
      }
      else{
        z <- 1 + xi*(data - mu)/sigma
        if(any(z < 0, na.rm = TRUE)){ loglik <- 1e+10}
        else{
          loglik <- mean(( log(sigma) + (1/xi +1)*log(z) + z^(-1/xi)),
                         na.rm = TRUE)
        }
      }
    }
  }
  loglik
}



#' GEV fit to disjoint block maxima
#'
#' @param data vector of block maxima
#' @param temp_cvrt vector of temporal covariate
#' @param type one of 'stationary', 'shift', 'scale'
#' @param method passed to optim
#' @param maxiter passed to optim
#' @param hessian logical; whether to return hessian
#' @param rel_trend refers to parametrisation of scale model, therefore only
#' relevant when 'type == "shift"'
#' @param return_cov logical; whether to return a covariance matrix
#'
#' @return
#' @export
#'
#' @examples
#' # bla
fit_gev_univ_dj <- function(data, temp_cvrt, type, method = "BFGS", maxiter = 200,
                            hessian = FALSE, rel_trend = FALSE, return_cov = FALSE) {

  start_st <-  evd::fgev(as.vector(data), std.err = FALSE)$estimate
  if(type == "stationary") {
    start_vals <- c(start_st[1], start_st[2], start_st[3])
    names(start_vals) <- c("loc", "scale", "shape")
  }
  else {
    start_vals <- c(start_st[1], start_st[2], start_st[3], 0)
    names(start_vals) <- c("mu", "sigma", "shape", "alpha")
  }
  #print(start_vals)
  mlest <- optim(start_vals, fn = nll_univ_dj,
                 data = data, temp_cvrt = temp_cvrt,
                 hessian = hessian, type = type,
                 method = method, control = list(maxit = maxiter), rel_trend = rel_trend
  )
  if(!(mlest$convergence == 0) ){print("Optimization didn't succeed.")}

  if(!hessian) {
    return(list(mle = mlest$par, nll_expsc = mlest$value, conv = mlest$convergence))
  }
  else {
    if(!return_cov) {

      return(list(mle = mlest$par, nll_expsc = mlest$value, conv = mlest$convergence,
                  hessian = mlest$hessian))
    }
    else {
      mlfit <- list(mle = mlest$par, nll_expsc = mlest$value, conv = mlest$convergence,
                    hessian = mlest$hessian)

      zbm <- (data - mlfit$mle[1] - mlfit$mle[4]*temp_cvrt)/mlfit$mle[2]

      scorez <- score.function(zbm, theta = c(0,1, mlfit$mle[3]))

      covzsc <- cov(t(scorez))

      hessinv <- solve(mlfit$hessian)

      meanct <- mean(temp_cvrt, na.rm = TRUE)
      meanctsq <- mean(temp_cvrt^2, na.rm = TRUE)

      cov1 <- covzsc
      cov1 <- cbind(cov1, meanct *covzsc[1:3, 1])
      cov1 <- rbind(cov1, c(meanct *covzsc[1:3, 1], meanctsq*covzsc[1,1] ))

      cov1 <- hessinv %*%  diag(c(1/mlfit$mle[2], 1/mlfit$mle[2], 1, 1/mlfit$mle[2])) %*%
        cov1 %*%  diag(c(1/mlfit$mle[2], 1/mlfit$mle[2], 1, 1/mlfit$mle[2])) %*% hessinv
      cov1 <- cov1/length(data)

      return(append(mlfit, list(covhat = cov1)))

    }

  }

}
