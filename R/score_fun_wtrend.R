### score functions

# compute the score function at a given parameter vector for either the stationary model,
# a linear shift or a trend that scales the observations


# passt auch für xi = 0
score.function_univ <- function(x, theta, temp.cov = NULL, type = "shift") {

  if(!(type %in% c("shift", "scale", "stationary"))){
    stop("Type must be either 'shift', 'scale' or 'stationary'.")}
  if(!(type == "stationary") & is.null(temp.cov)) {stop("Temporal covariate must be provided.")}
  if(type == "shift") {

    mu0 <- theta["loc0"]
    mu1 <- theta["tempLoc1"]
    sigma <- theta["scale0"]
    xi <- theta["shape"]

    mu <- mu0 + mu1*temp.cov

    z <- (x-mu)/sigma
    in_supp <- which( (1+ xi*z) > 0)
    z <- z[in_supp]
    if(abs(xi) < 1e-08) {
      u <- exp(-z)

    }
    else {
      u <- (1+xi *z)^(-1/xi)
    }

    scoreloc0 <- (xi + 1 - u)/(sigma*(1+xi*z))
    scoretemploc  <- (xi + 1 - u)/(sigma*(1+xi*z))*temp.cov
    scorescale <- ((1-u)*z -1)/(sigma*(1+xi*z))

    if(abs(xi) < 1e-08){
      scoreshape <- (1-u)*z^2/2- z
    }
    else{
      scoreshape <- (1-u)*1/xi*(1/xi*log(1+xi*z)-z/(1+xi*z)) - z/(1+xi*z)
    }

    return(matrix(c(scoreloc0, scorescale, scoreshape, scoretemploc),
                  nrow = 4, byrow = TRUE))
  }

  if(type == "scale") {

    mu0 <- theta["mu0"]
    alpha0 <- theta["alpha"]
    sigma0 <- theta["sigma0"]
    xi <- theta["gamma"]


    mut <- mu0*exp(alpha0*temp.cov/mu0)
    sigmat <-  sigma0*exp(alpha0*temp.cov/mu0)

    zt <- (x - mut)/sigmat
    in_supp <- which(1+ xi*zt > 0)
    zt <- zt[in_supp]

    if(abs(xi) < 1e-08) {
      ut <- exp(-zt)
    }
    else {
      ut <- (1+xi *zt)^(-1/xi)
    }

    expo <- exp(alpha0*temp.cov/mu0)

    dmu0ut <- -ut^(xi +1)*(alpha0 *temp.cov *x/(mu0^2*sigmat) - 1/sigma0)

      # ut^(xi +1)/sigmat* (expo*(1 - alpha0*temp.cov/mu0) - (x - mut)*alpha0*temp.cov/mu0)
    dsigma0ut <- - ut^(xi +1)/sigma0*(-x/sigmat + mu0/sigma0)

    dalpha0ut <- ut^(xi +1)/sigmat* x*temp.cov/mu0

    scoreloc0 <- alpha0*temp.cov/mu0^2 + dmu0ut*( (xi +1)/ut -1)

    scorescale <- -1/sigma0 + dsigma0ut*( (xi +1)/ut -1)

    scorealpha <- -temp.cov/mu0  + dalpha0ut*( (xi +1)/ut  -1)

    if(abs(xi) < 1e-08){
      scoreshape <- (1-ut)*zt^2/2- zt
    }
    else{
      scoreshape <- (1-ut)*1/xi*(1/xi*log(1+xi*zt)-zt/(1+xi*zt)) - zt/(1+xi*zt)
    }

    return(matrix(c(scoreloc0, scorescale, scoreshape, scorealpha),
                  nrow = 4, byrow = TRUE))



  }

  if(type == "stationary") {
    mu <- theta[1]
    sigma <- theta[2]
    xi <- theta[3]
    z <- (x-mu)/sigma
    in_supp <- which( (1+ xi*z) > 0)
    z <- z[in_supp]
    if(abs(xi) < 1e-08) {
      u <- exp(-z)
    }
    else {
      u <- (1+xi *z)^(-1/xi)
    }

    scoreloc <- (xi + 1 - u)/(sigma*(1+xi*z))
    scorescale <- ((1-u)*z -1)/(sigma*(1+xi*z))

    if(abs(xi) < 1e-08){
      scoreshape <- (1-u)*z^2/2- z
    }
    else{
      scoreshape <- (1-u)*1/xi*(1/xi*log(1+xi*z)-z/(1+xi*z)) - z/(1+xi*z)
    }

    return(matrix(c(scoreloc, scorescale, scoreshape),
                  nrow = 3, byrow = TRUE))

  }


}

AllDiags <- function(inmat, sorted = TRUE) {
  Range <- ncol(inmat) - 1
  Range <- 0:Range
  if (isTRUE(sorted)) Range <- Range[order(abs(Range))]
  lapply(Range, function(x) {
    inmat[row(inmat) == (col(inmat) - x)]
  })
}

#' GEV parameter covariance matrix for sliding blocks ML estimator
#'
#' @param orig_slbm Full sample of sliding block maxima
#' @param est_par Output of \code{\link{fit_gev_univ}}.
#' @param blcksz The blocksize used.
#' @param temp.cov The temporal covariate (must be same length as \code{orig_slbm}).
#' @param type Either 'shift', 'scale' or 'stationary', depending if one assumes
#' that obervations are
#'  shifted or scaled w.r.t. the covariate or stationary (for more details, see
#'  \code{\link[slbm]{fit_gev_univ}}).
#' @param varmeth One of 'V', 'V2' or 'both'. Determines the method used for estimating
#' the covariance matrix.
#'
#' @return
#' @export
#'
#' @examples
#' ### simulate some data with a linear trend
#' xx <- evd::rgpd(90*100, shape = 0.2) + 2*rep(1:100/100, each = 90)
#' ## temporal covariate for sliding BM
#' temp_cvrt <- rep(1:100/100, each = 90)[1:8911]
#'
#' ## compute sample of uniuqe sliding BM
#' bms <- get_uniq_bm(xx, 90, temp_cvrt = temp_cvrt, looplastblock = FALSE)
#'
#' ## full sample of sliding for estimaing the covariance matrix
#' slbm <- blockmax(xx, r = 90, "sliding")
#' estim <- fit_gev_univ(data = bms, type = "shift", hessian = TRUE)
#' est_var_univ(slbm, est_par = estim, blcksz = 90, temp.cov = temp_cvrt,
#' varmeth = "both", type = "shift")

est_var_univ <- function(orig_slbm, est_par, blcksz,  temp.cov =  NULL, type = "shift",
                           varmeth = "both"){
  nsl <- length(orig_slbm)
  k <-  nsl /blcksz
  k <- ifelse(k == floor(k), k, floor(k))

  if(!(type %in% c("scale", "shift", "stationary"))) {
    stop("Type must be one of 'scale', 'shift', 'stationary'.")
  }
  if(is.null(est_par$hessian)) { stop("Hessian of log likelihood is missing.")}


  score.bdata <- score.function_univ(orig_slbm, theta =  est_par$mle, temp.cov = temp.cov,
                                       type = type)

  fishest <- est_par$hessian/nsl

  fishestinv <- solve(fishest)

  Y <- (fishestinv %*% score.bdata)

  useobs <- k*blcksz

  if(!(type == "stationary")) {

    Yloc <- Y[1, 1:useobs]
    Yscale <- Y[2, 1:useobs ]
    Yshape <- Y[3, 1:useobs]
    Ytrend <- Y[4, 1:useobs]

    m  <- t(matrix(Yshape, nrow = blcksz, ncol = k))
    cm <- cov(m)
    msc  <- t(matrix(Yscale, nrow = blcksz, ncol = k))
    cmsc <- cov(msc)
    mloc  <- t(matrix(Yloc, nrow = blcksz, ncol = k))
    cmloc <- cov(mloc)
    mtrend  <- t(matrix(Ytrend, nrow = blcksz, ncol = k))
    cmtrend <- cov(mtrend)
    cmshsc <- cov(m, msc)
    cmshloc <- cov(m, mloc)
    cmshtrend <- cov(m, mtrend)
    cmscloc <- cov(msc, mloc)
    cmsctrend <- cov(msc, mtrend)
    cmloctrend <- cov(mloc, mtrend)

    if(varmeth %in% c("both", "V")) {

      vshape <- 2*mean(cm[1,]) / k
      vscale <-  2*mean(cmsc[1,]) / k
      vloc <-  2*mean(cmloc[1,]) / k
      vtrend <-  2*mean(cmtrend[1,]) / k
      vshsc <- 2*mean(cmshsc[1,]) / k
      vshloc <- 2*mean(cmshloc[1,]) / k
      vshtrend <- 2*mean(cmshtrend[1,]) / k
      vscloc <- 2*mean(cmscloc[1,]) / k
      vsctrend <- 2*mean(cmsctrend[1,]) / k
      vloctrend <- 2*mean(cmloctrend[1,]) / k

      V <- matrix( c(vloc, vscloc, vshloc, vloctrend,
                     vscloc, vscale, vshsc, vsctrend,
                     vshloc, vshsc, vshape, vshtrend,
                     vloctrend, vsctrend, vshtrend, vtrend),
                   nrow = 4)
      rownames(V) <- colnames(V) <- c("loc", "scale", "shape", "trendpar")

    }

    if(varmeth %in% c("both", "V2")) {

      v2shape <- 2* mean(as.numeric(lapply(AllDiags(cm), mean)) )/ k
      v2scale <- 2* mean(as.numeric(lapply(AllDiags(cmsc), mean)) )/ k
      v2loc <- 2* mean(as.numeric(lapply(AllDiags(cmloc), mean)) )/ k
      v2trend <- 2* mean(as.numeric(lapply(AllDiags(cmtrend), mean)) )/ k
      v2shsc <- 2* mean(as.numeric(lapply(AllDiags(cmshsc), mean)) )/ k
      v2shloc <- 2* mean(as.numeric(lapply(AllDiags(cmshloc), mean)) )/ k
      v2shtrend <- 2* mean(as.numeric(lapply(AllDiags(cmshtrend), mean)) )/ k

      v2scloc <- 2* mean(as.numeric(lapply(AllDiags(cmscloc), mean)) )/ k
      v2sctrend <- 2* mean(as.numeric(lapply(AllDiags(cmsctrend), mean)) )/ k
      v2loctrend <- 2* mean(as.numeric(lapply(AllDiags(cmloctrend), mean)) )/ k


      V2 <- matrix( c(v2loc, v2scloc, v2shloc, v2loctrend,
                      v2scloc, v2scale, v2shsc, v2sctrend,
                      v2shloc, v2shsc, v2shape, v2shtrend,
                      v2loctrend, v2sctrend, v2shtrend, v2trend),
                    nrow = 4)
      rownames(V2) <- colnames(V2) <- c("loc", "scale", "shape", "trendpar")
    }
  }
  if(type == "stationary") {

    Yloc <- Y[1, 1:useobs]
    Yscale <- Y[2, 1:useobs ]
    Yshape <- Y[3, 1:useobs]

    m  <- t(matrix(Yshape, nrow = blcksz, ncol = k))
    cm <- cov(m)
    msc  <- t(matrix(Yscale, nrow = blcksz, ncol = k))
    cmsc <- cov(msc)
    mloc  <- t(matrix(Yloc, nrow = blcksz, ncol = k))
    cmloc <- cov(mloc)
    cmshsc <- cov(m, msc)
    cmshloc <- cov(m, mloc)
    cmscloc <- cov(msc, mloc)

    if(varmeth %in% c("both", "V")) {

      vshape <- 2*mean(cm[1,]) / k
      vscale <-  2*mean(cmsc[1,]) / k
      vloc <-  2*mean(cmloc[1,]) / k
      vshsc <- 2*mean(cmshsc[1,]) / k
      vshloc <- 2*mean(cmshloc[1,]) / k
      vscloc <- 2*mean(cmscloc[1,]) / k

      V <- matrix( c(vloc, vscloc, vshloc,
                     vscloc, vscale, vshsc,
                     vshloc, vshsc, vshape),
                   nrow = 3)
      rownames(V) <- colnames(V) <- c("loc", "scale", "shape")

    }

    if(varmeth %in% c("both", "V2")) {

      v2shape <- 2* mean(as.numeric(lapply(AllDiags(cm), mean)) )/ k
      v2scale <- 2* mean(as.numeric(lapply(AllDiags(cmsc), mean)) )/ k
      v2loc <- 2* mean(as.numeric(lapply(AllDiags(cmloc), mean)) )/ k
      v2shsc <- 2* mean(as.numeric(lapply(AllDiags(cmshsc), mean)) )/ k
      v2shloc <- 2* mean(as.numeric(lapply(AllDiags(cmshloc), mean)) )/ k
      v2scloc <- 2* mean(as.numeric(lapply(AllDiags(cmscloc), mean)) )/ k


      V2 <- matrix( c(v2loc, v2scloc, v2shloc,
                      v2scloc, v2scale, v2shsc,
                      v2shloc, v2shsc, v2shape),
                    nrow = 3)
      rownames(V2) <- colnames(V2) <- c("loc", "scale", "shape")
    }

  }

  if(varmeth == "both") {return(list(V = V, V2 = V2))}
  if(varmeth == "V"){return(list(V = V))}
  if(varmeth == "V2"){return(list(V2 = V2))}

}

qdelta_rl <- function(theta, Tyrl, type, ref_gmst = NULL) {

  ct <- -log(1-1/Tyrl)

  if(type == "shift") {
    if(is.null(ref_gmst)) { stop("Reference covariate value must be provided.")}
    if(is.null(names(theta))) { names(theta) <- c("loc0", "scale0", "shape", "TempLoc1")}
    mu0 <- theta["loc0"]
    mu1 <- theta["tempLoc1"]
    sigma <- theta["scale0"]
    xi <- theta["shape"]

    dmu0 <- 1
    if(abs(xi) < 1e-08) {
      dsigma <- -log(ct)
      dxi <- sigma*log(ct)^2/2
    }
    else {
      dsigma <- (ct^(-xi) -1)/xi
      dxi <- sigma/xi^2*(1-ct^(-xi)*(xi*log(ct) +1))
    }
    dmu1 <- ref_gmst

    return(
      c("mu0" = dmu0, "sigma0" = unname(dsigma), "shape" = unname(dxi), "loctrend" = dmu1))

  }
  if(type == "scale") {
    if(is.null(ref_gmst)) { stop("Reference covariate value must be provided.")}
    if(is.null(names(theta))) { names(theta) <- c("mu0", "sigma0", "gamma", "alpha")}

    mu0 <- theta["mu0"]
    alpha0 <- theta["alpha"]
    sigma0 <- theta["sigma0"]
    xi <- theta["gamma"]

    expo <- exp(alpha0*ref_gmst/mu0)

    if(abs(xi) < 1e-08) {
      dmu0 <- expo*( log(ct) *alpha0*ref_gmst/mu0^2*sigma0 + 1 - alpha0*ref_gmst/mu0)
      dsigma0 <- expo*log(ct)
      dxi <- sigma0*expo*log(ct)^2/2
      dalpha <- expo*ref_gmst*( sigma0/mu0 * (-log(ct)) +1)
    }
    else {

      dmu0 <- expo*( - (ct^(-xi) -1)/xi *alpha0*ref_gmst/mu0^2*sigma0 + 1 - alpha0*ref_gmst/mu0)

      dsigma0 <- expo*(ct^(-xi) -1)/xi

      dxi <- sigma0*expo/xi^2*(1-ct^(-xi)*(xi*log(ct) +1))

      dalpha <- expo*ref_gmst*( sigma0/mu0 * (ct^(-xi) -1)/xi +1)
    }

    return(
      c("mu0" = unname(dmu0), "sigma0" = unname(dsigma0), "shape" = unname(dxi),
      "alpha0" = unname(dalpha))
    )
  }
  if(type == "stationary"){

    mu <- theta["loc"]
    sigma <- theta["scale"]
    xi <- theta["shape"]

    if(abs(xi) < 1e-08) {
      return(c(1, -log(ct),  log(ct)^2/2 ))
    } else {
      return(c(1, (ct^(-xi) -1)/xi, sigma/xi^2*(1-ct^(-xi)*(xi*log(ct) +1) )))
    }
  }
}



#' Estimate Variance of RL estimation
#'
#' @param theta A named of parameter estimates in the order:
#' constant part of location parameter, constant part of scale parameter, shape parameter,
#' trend parameter.
#' @param Tyrl The period for which the variance of the corresponding Return Level is to be
#'  estimated.
#' @param type One of 'scale', 'shift', or 'stationary', depending whether the model assumes
#'  a scale, a shift or stationary behaviour of observations with respect to the temporal covariate.
#' @param ref_gmst A numeric vector giving reference values of the temporal covariate
#' at which the variance of the corresponding RL is estimated, or \code{NULL} when stationarity
#' is assumed.
#' @param Covmat The (estimated) covariance matrix of the parameter vector \code{theta}.
#'
#' @return
#' @export
#'
#' @examples
#' ### simulate some data with a linear trend
#' xx <- evd::rgpd(90*100, shape = 0.2) + 2*rep(1:100/100, each = 90)
#' ## temporal covariate for sliding BM
#' temp_cvrt <- rep(1:100/100, each = 90)[1:8911]
#'
#' ## compute sample of unique sliding BM
#' bms <- get_uniq_bm(xx, 90, temp_cvrt = temp_cvrt, looplastblock = FALSE)
#'
#' ## full sample of sliding for estimaing the covariance matrix
#' slbm <- blockmax(xx, r = 90, "sliding")
#' estim <- fit_gev_univ(data = bms, type = "shift", hessian = TRUE)
#' covest <- est_var_univ(slbm, est_par = estim, blcksz = 90, temp.cov = temp_cvrt,
#' varmeth = "V2", type = "shift")
#'
#' ## define reference value of temporal covariate
#' reft <- 0.8
#'
#' ## location parameter and 100-year RL in that climate
#' loctr <- estim$mle[1] + estim$mle[4]*reft
#' rlhat <- evd::qgev(1-1/100, loc = loctr, scale = estim$mle[2], shape = estim$mle[3])
#' rlhat
#' estimate_var_rl(estim$mle, Tyrl = 100, type = "shift", ref_gmst = reft, Covmat = covest$V)
#'
#'
#'#### assuming stationarity
#' xx <- evd::rgpd(90*100, shape = 0.2)
#'
#' ## compute sample of unique sliding BM
#' bms <- get_uniq_bm(xx, 90, looplastblock = TRUE)
#'
#' ## full sample of sliding for estimaing the covariance matrix
#' slbm <- blockmax(xx, r = 90, "sliding")
#' estim <- fit_gev_univ(data = bms, type = "stationary", hessian = TRUE)
#' covest <- est_var_univ(slbm, est_par = estim, blcksz = 90, varmeth = "V2", type = "stationary")
#'
#' # estimated 100-year RL
#' rlhat <- evd::qgev(1-1/100, loc = estim$mle[1], scale = estim$mle[2], shape = estim$mle[3])
#' rlhat
#' estimate_var_rl(estim$mle, Tyrl = 100, type = "stationary", Covmat = covest$V)
#'
estimate_var_rl <- function(theta, Tyrl, type, ref_gmst = NULL, Covmat) {

  if(type == "stationary") {
    ref_gmst <- 1
  }


  qdelt <- purrr::map(ref_gmst, ~
                        qdelta_rl(theta = theta, Tyrl = Tyrl, type = type, ref_gmst = .x)
  )

  res <- purrr::map_dbl( qdelt, ~ .x %*% Covmat %*% .x)
  if(!(type == "stationary")) {
    names(res) <- paste("refGMST", ref_gmst)
  }
  res[res < 0] <- NA
  res
}

#
# xx <- evd::rgpd(100*90, shape = 0.2) + rep(1:100/100, each = 90)
# df.xx <- data.frame(Obs = xx, Station = "X1")
# unbm <- get_uniq_bm(data = df.xx, blcksz = 90, temp_cvrt =rep(1:100/100, each = 90)[1:8911])
# est_sl <- fit_spgev_sl(unbm, loc.sp.form = ~1, scale.sp.form = ~1, loc.temp.form = ~ GMST,
#                     spat.cov = data.frame(lon = 1),
#                     st_val_meth = "LeastSqTemp",
#                     datastart = blockmax(xx, r = 90, "disjoint"),
#                     temp.cov = 1:100/100, return_hess = TRUE,
#                     scale.link = make.link("identity")
#                     )
#
# est_sl
# aa <- score.function_wtrend(x = 3.64, theta = est_sl$mle,
#                       temp.cov = 0.15)
# aa
#
# rowMeans(aa)
# unbm1 <- unbm %>% mutate(uniq_data = map(uniq_data, ~ .x[7, ]))
#
# prepdata <- prep4spatml(loc.sp.form = ~ 1,
#                scale.sp.form = ~ 1,
#                loc.temp.form = ~ GMST, data = unbm1, spat.cov = data.frame(lon= 1))
#
# numDeriv::grad( function(x){ - nll_spat_temp_sl_prep(params = c("loc0" = x,
#                                                                 "scale0" = est_sl$mle[2],
#                                                                 "shape" = est_sl$mle[3],
#                                                                 "tempLoc1" = est_sl$mle[4]),
#                                                      loc.sp.form = ~ 1,
#                                                      scale.sp.form = ~1,
#                                                      loc.temp.form = ~ GMST,
#                                                      dataprep = prepdata,
#                                                      spat.cov = data.frame(lon = 1),
#                                                      scale.link = make.link("identity"))} ,
#   x =  est_sl$mle[1])
#
#
#' set.seed(3)
#' blcksz <- 90
#' xx <- evd::rgpd(100*90, shape = 0.2)
#' #'
#' yy <- data.frame(Station = "X1", Obs = xx)
#' #'
#' # define a temporal covariate that is constant over a block of length blcksz
#' temp_cvrt <- rep(1:100/100, each = blcksz)[1:(99*blcksz + 1)]
#' bms <- get_uniq_bm(yy, blcksz, temp_cvrt = temp_cvrt, looplastblock = FALSE)
#' bms <- bms$uniq_data[[1]]
#' bms
#' slbm <- blockmax(xx, r = 90, "sliding")
#' bms[11, ]
#' #'###############################
#' #' fit_gev_expsc(data = bms)
#' est_sl <- fit_gev_expsc(data = bms, hessian = TRUE)
#' est_sl
#'
#' aa <- score.function_wtrend(x =  4.519812, theta = est_sl$mle,
#'                         temp.cov = 0.45, type = "scale")
#' aa
#'
#' numDeriv::grad( function(x){ - nll_expsc(params = c("loc0" = x,  #est_sl$mle[1],
#'                                                                 "sigma0" = est_sl$mle[2],
#'                                                                 "gamma" = est_sl$mle[3],
#'                                                                 "alpha" = est_sl$mle[4]),
#'                                                      data = bms[11, ] )} ,
#'                 x =  est_sl$mle[1])
#'
#' numDeriv::grad( function(x){ - nll_expsc(params = c("loc0" = est_sl$mle[1],
#'                                                     "sigma0" = x, # est_sl$mle[2],
#'                                                     "gamma" = est_sl$mle[3],
#'                                                     "alpha" = est_sl$mle[4]),
#'                                          data = bms[11, ] )} ,
#'                 x =  est_sl$mle[2])
#'
#' numDeriv::grad( function(x){ - nll_expsc(params = c("loc0" = est_sl$mle[1],
#'                                                     "sigma0" = est_sl$mle[2],
#'                                                     "gamma" = x, #est_sl$mle[3],
#'                                                     "alpha" = est_sl$mle[4]),
#'                                          data = bms[11, ] )} ,
#'                 x =  est_sl$mle[3])
#'
#' numDeriv::grad( function(x){ - nll_expsc(params = c("loc0" = est_sl$mle[1],
#'                                                     "sigma0" = est_sl$mle[2],
#'                                                     "gamma" = est_sl$mle[3],
#'                                                     "alpha" = x ) ,# est_sl$mle[4]),
#'                                          data = bms[11, ] )} ,
#'                 x =  est_sl$mle[4])
#'
#' est_var_wtrend(slbm, est_par = est_sl, blcksz = 90, temp.cov = rep(1:100/100, each = 90)[1:8911],
#'                type = "shift")
#'
