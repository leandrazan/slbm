# ... additional arguments:
# *rel_trend : logical, whether in scale model one estimates the
# trend parameter relative to the location parameter

# xx <- evd::rgpd(90*100, shape = 0.2) + 2*rep(1:100/100, each = 90)
# # ## temporal covariate for sliding BM
#  temp_cvrt <- rep(1:100/100, each = 90)[1:(90*100 - 90 +1)]
# #
# # ## compute sample of uniuqe sliding BM
#  bms <- get_uniq_bm(xx, 90, temp_cvrt = temp_cvrt, looplastblock = FALSE)
# #
# # ## full sample of sliding for estimaing the covariance matrix
#  slbm <- blockmax(xx, r = 90, "sliding")
#  estim <- fit_gev_univ(data = bms, type = "shift", hessian = TRUE)
# # aa <- score.fun(slbm, theta = estim$mle, temp.cov = temp_cvrt, type = "shift", chain = FALSE)
# # bb <- score.fun(slbm, theta = estim$mle, temp.cov = temp_cvrt, type = "shift", chain = TRUE)
# # bb <- rbind(bb,  bb[1,]*temp_cvrt)
# # identical(aa, bb)
#
# tempc <- rep(as.numeric(scale(1:100/100)), each = 90)
# xx <- evd::rgpd(90*100, shape = 0.2)*2*exp(tempc*1.1)
# ## temporal covariate for sliding BM
# temp_cvrt <- tempc[1:(90*100 - 90 +1)]
#
# ## compute sample of uniuqe sliding BM
# bms <- get_uniq_bm(xx, 90, temp_cvrt = temp_cvrt, looplastblock = FALSE)
#
# # full sample of sliding for estimaing the covariance matrix
#  slbm <- blockmax(xx, r = 90, "sliding")
# estim <- fit_gev_univ(data = bms, type = "scale", hessian = TRUE, rel_trend = FALSE)
# aa <- score.fun(slbm, theta = estim$mle, temp.cov = temp_cvrt, type = "scale", chain = FALSE, rel_trend = FALSE)
#
# bb <- score.fun(slbm, theta = estim$mle, temp.cov = temp_cvrt, type = "scale", chain = TRUE, rel_trend = FALSE)

# bb[1, ] <- exp(estim$mle[4]*temp_cvrt)*bb[1, ]
# bb[2, ] <- exp(estim$mle[4]*temp_cvrt)*bb[2, ]
# bb <- rbind(bb, estim$mle[1]*temp_cvrt*bb[1, ] + estim$mle[2]*temp_cvrt*bb[2, ] )
# all.equal(aa, bb)
# #
# estim <- fit_gev_univ(data = bms, type = "scale", hessian = TRUE, rel_trend = TRUE)
# aa <- score.fun(slbm, theta = estim$mle, temp.cov = temp_cvrt, type = "scale", chain = FALSE, rel_trend = TRUE)
# bb <- score.fun(slbm, theta = estim$mle, temp.cov = temp_cvrt, type = "scale", chain = TRUE, rel_trend = TRUE)
# bb2 <- bb
# muh <- estim$mle[1]
# sigmah <- estim$mle[2]
# alphh <- estim$mle[4]
# bb2[1, ] <- exp(alphh/muh*temp_cvrt)* (
#   (1-alphh*temp_cvrt/muh)* bb[1, ]  - sigmah*alphh*temp_cvrt/muh^2*bb[2, ])
#
# bb2[2, ] <- exp(alphh/muh*temp_cvrt)*bb[2, ]
#
# bb2 <- rbind(bb2, exp(alphh/muh*temp_cvrt)*
#                ( temp_cvrt*bb[1, ] + sigmah/muh*temp_cvrt*bb[2, ] ))
# all.equal(aa, bb2)

## if chain = TRUE, the three dimensional scorefunction of the stationary model evaluated
## in the respective parameters is returned. This is not the score function of the
## non-stationary model, since it needs to be multiplied with the T_sigma or B_c() matrix
# score.fun <- function(x, theta, temp.cov = NULL, type = "shift", chain = TRUE,
#                                 ...) {
#
#   add.args <- list(...)
#   if(!(type %in% c("shift", "scale", "stationary"))){
#     stop("Type must be either 'shift', 'scale' or 'stationary'.")}
#   if(!(type == "stationary") & is.null(temp.cov)) {stop("Temporal covariate must be provided.")}
#
#   if(chain) {
#     if(type == "shift") {
#
#       mu0 <- theta["mu"]
#       mu1 <- theta["alpha"]
#       sigma <- theta["sigma"]
#       xi <- theta["shape"]
#
#       mu <- mu0 + mu1*temp.cov
#
#       z <- (x-mu)/sigma
#       in_supp <- which( (1+ xi*z) > 0)
#       z[!in_supp] <- NA
#
#       if(abs(xi) < 1e-08) {
#         u <- exp(-z)
#       }
#       else {
#         u <- (1+xi *z)^(-1/xi)
#       }
#
#       scoreloc <- (xi + 1 - u)/(sigma*(1+xi*z))
#       scorescale <- ((1-u)*z -1)/(sigma*(1+xi*z))
#
#       if(abs(xi) < 1e-08){
#         scoreshape <- (1-u)*z^2/2- z
#       }
#       else{
#         scoreshape <- (1-u)*1/xi*(1/xi*log(1+xi*z)-z/(1+xi*z)) - z/(1+xi*z)
#       }
#
#       scoremat <- matrix(c(scoreloc, scorescale, scoreshape),
#                          nrow = 3, byrow = TRUE)
#
#
#       return(scoremat)
#
#     }
#     if(type == "scale") {
#       # whether trend parameter is set relative to location parameter
#       rel_trend <- add.args$rel_trend
#       if(is.null(rel_trend)) {stop("Please specify the parametrisation of the scale model vie 'rel_trend'.")}
#       mu0 <- theta["mu"]
#       alpha0 <- theta["alpha"]
#       sigma0 <- theta["sigma"]
#       xi <- theta["shape"]
#
#       if(rel_trend) {
#         mut <- mu0*exp(alpha0*temp.cov/mu0)
#         sigmat <-  sigma0*exp(alpha0*temp.cov/mu0)
#       }
#       else {
#         mut <- mu0*exp(alpha0*temp.cov)
#         sigmat <- sigma0*exp(alpha0*temp.cov)
#       }
#
#       zt <- (x - mut)/sigmat
#       in_supp <- which(1+ xi*zt > 0)
#       zt[!in_supp] <- NA
#
#       if(abs(xi) < 1e-08) {
#         ut <- exp(-zt)
#       }
#       else {
#         ut <- (1+xi *zt)^(-1/xi)
#       }
#
#       scoreloc <- (xi + 1 - ut)/(sigmat*(1+xi*zt))
#       scorescale <- ((1-ut)*zt -1)/(sigmat*(1+xi*zt))
#
#       if(abs(xi) < 1e-08){
#         scoreshape <- (1-ut)*zt^2/2- zt
#       }
#       else{
#         scoreshape <- (1-ut)*1/xi*(1/xi*log(1+xi*zt)-zt/(1+xi*zt)) - zt/(1+xi*zt)
#       }
#
#       scoremat <- matrix(c(scoreloc, scorescale, scoreshape),
#                          nrow = 3, byrow = TRUE)
#
#       return(scoremat)
#
#     }
#   }
#   else {
#     if(type == "shift") {
#
#       mu0 <- theta["mu"]
#       mu1 <- theta["alpha"]
#       sigma <- theta["sigma"]
#       xi <- theta["shape"]
#
#       mu <- mu0 + mu1*temp.cov
#
#       z <- (x-mu)/sigma
#       in_supp <- which( (1+ xi*z) > 0)
#       z[!in_supp] <- NA
#       if(abs(xi) < 1e-08) {
#         u <- exp(-z)
#
#       }
#       else {
#         u <- (1+xi *z)^(-1/xi)
#       }
#
#       scoreloc0 <- (xi + 1 - u)/(sigma*(1+xi*z))
#       scoretemploc  <- scoreloc0*temp.cov
#       scorescale <- ((1-u)*z -1)/(sigma*(1+xi*z))
#
#       if(abs(xi) < 1e-08){
#         scoreshape <- (1-u)*z^2/2- z
#       }
#       else{
#         scoreshape <- (1-u)*1/xi*(1/xi*log(1+xi*z)-z/(1+xi*z)) - z/(1+xi*z)
#       }
#
#       return(matrix(c(scoreloc0, scorescale, scoreshape, scoretemploc),
#                     nrow = 4, byrow = TRUE))
#     }
#
#     if(type == "scale") {
#
#       mu0 <- theta["mu"]
#       alpha0 <- theta["alpha"]
#       sigma0 <- theta["sigma"]
#       xi <- theta["shape"]
#
#       if(!add.args$rel_trend) {
#
#         mut <- mu0*exp(alpha0*temp.cov)
#         sigmat <-  sigma0*exp(alpha0*temp.cov)
#
#         zt <- (x - mut)/sigmat
#         in_supp <- which(1+ xi*zt > 0)
#         zt[!in_supp] <- NA
#
#         if(abs(xi) < 1e-08) {
#           ut <- exp(-zt)
#         }
#         else {
#           ut <- (1+xi *zt)^(-1/xi)
#         }
#
#         # expo <- exp(alpha0*temp.cov/mu0)
#
#         dmu0ut <- ut^(xi +1)/sigma0
#
#         # ut^(xi +1)/sigmat* (expo*(1 - alpha0*temp.cov/mu0) - (x - mut)*alpha0*temp.cov/mu0)
#         dsigma0ut <- ut^(xi +1)*zt/sigma0
#
#         dalpha0ut <- ut^(xi +1)*(zt + mu0/sigma0)*temp.cov
#
#         scoreloc0 <-  dmu0ut*( (xi +1)/ut -1)
#
#         scorescale <- -1/sigma0 + dsigma0ut*( (xi +1)/ut -1)
#
#         scorealpha <- -temp.cov + dalpha0ut*( (xi +1)/ut  -1)
#
#         if(abs(xi) < 1e-08){
#           scoreshape <- (1-ut)*zt^2/2- zt
#         }
#         else{
#           scoreshape <- (1-ut)*1/xi*(1/xi*log(1+xi*zt)-zt/(1+xi*zt)) - zt/(1+xi*zt)
#         }
#
#         return(matrix(c(scoreloc0, scorescale, scoreshape, scorealpha),
#                       nrow = 4, byrow = TRUE))
#
#       }
#       else {
#         mut <- mu0*exp(alpha0*temp.cov/mu0)
#         sigmat <-  sigma0*exp(alpha0*temp.cov/mu0)
#
#         zt <- (x - mut)/sigmat
#         in_supp <- which(1+ xi*zt > 0)
#         zt[!in_supp] <- NA
#
#         if(abs(xi) < 1e-08) {
#           ut <- exp(-zt)
#         }
#         else {
#           ut <- (1+xi *zt)^(-1/xi)
#         }
#
#         # expo <- exp(alpha0*temp.cov/mu0)
#
#         dmu0ut <- -ut^(xi +1)*(alpha0 *temp.cov *x/(mu0^2*sigmat) - 1/sigma0)
#
#         # ut^(xi +1)/sigmat* (expo*(1 - alpha0*temp.cov/mu0) - (x - mut)*alpha0*temp.cov/mu0)
#         dsigma0ut <- - ut^(xi +1)/sigma0*(-x/sigmat + mu0/sigma0)
#
#         dalpha0ut <- ut^(xi +1)/sigmat* x*temp.cov/mu0
#
#         scoreloc0 <- alpha0*temp.cov/mu0^2 + dmu0ut*( (xi +1)/ut -1)
#
#         scorescale <- -1/sigma0 + dsigma0ut*( (xi +1)/ut -1)
#
#         scorealpha <- -temp.cov/mu0  + dalpha0ut*( (xi +1)/ut  -1)
#
#         if(abs(xi) < 1e-08){
#           scoreshape <- (1-ut)*zt^2/2- zt
#         }
#         else{
#           scoreshape <- (1-ut)*1/xi*(1/xi*log(1+xi*zt)-zt/(1+xi*zt)) - zt/(1+xi*zt)
#         }
#       }
#
#       return(matrix(c(scoreloc0, scorescale, scoreshape, scorealpha),
#                     nrow = 4, byrow = TRUE))
#
#
#
#     }
#
#     if(type == "stationary") {
#       mu <- theta[1]
#       sigma <- theta[2]
#       xi <- theta[3]
#       z <- (x-mu)/sigma
#       in_supp <- which( (1+ xi*z) > 0)
#       z[!in_supp] <- NA
#       if(abs(xi) < 1e-08) {
#         u <- exp(-z)
#       }
#       else {
#         u <- (1+xi *z)^(-1/xi)
#       }
#
#       scoreloc <- (xi + 1 - u)/(sigma*(1+xi*z))
#       scorescale <- ((1-u)*z -1)/(sigma*(1+xi*z))
#
#       if(abs(xi) < 1e-08){
#         scoreshape <- (1-u)*z^2/2- z
#       }
#       else{
#         scoreshape <- (1-u)*1/xi*(1/xi*log(1+xi*z)-z/(1+xi*z)) - z/(1+xi*z)
#       }
#
#       return(matrix(c(scoreloc, scorescale, scoreshape),
#                     nrow = 3, byrow = TRUE))
#
#     }
#   }
#
#
# }


#
# blcksz <- 183
# xx <- dailyPrec %>% filter(lubridate::month(Date) %in% 4:9,
#                            lubridate::year(Date) >= 1950)
# gmst <- GMST %>% filter(Year >= 1950, Year <= 2010)
# bms <- get_uniq_bm(xx$cumPrec, blcksz, temp_cvrt = rep(gmst$smoothedGMST, each = blcksz) , looplastblock = FALSE)
#
# # full sample of sliding for estimaing the covariance matrix
# slbm <- blockmax(xx$cumPrec, r = blcksz, "sliding")
# estim <- fit_gev_univ(data = bms, type = "scale", hessian = TRUE, rel_trend = FALSE)
#
# estim
#
# bb <- score.fun(slbm, theta = estim$mle, temp.cov = rep(gmst$smoothedGMST, each = blcksz)[1:length(slbm)],
#                 type = "scale", chain = FALSE, rel_trend = FALSE)
#
# est_var_chain(orig_slbm = slbm, est_par =  estim, blcksz =  blcksz,
#               temp.cov = rep(gmst$smoothedGMST, each = blcksz)[1:length(slbm)],
#               type = "scale", varmeth = "both", rel_trend = FALSE)

#
# est_var_univ(slbm, est_par = estim, blcksz = 90, temp.cov = temp_cvrt,
#  varmeth = "both", type = "scale")
#
# est_var_univ_shift(slbm, est_par = estim, blcksz = blcksz, temp.cov = temp_cvrt, temp.cov.dj = temp_cvrt)
# est_var_chain(slbm, estim, blcksz, temp_cvrt, "shift", "both")
# est_var_chain(slbm, estim, blcksz, temp_cvrt, "scale", "both", rel_trend = FALSE)


#' GEV parameter covariance matrix for sliding blocks ML estimator in shift model
#' @inheritParams est_var_univ
#' @param temp.cov.dj Optional: values of the temporal covariate one would use for
#' disjoint blocks.
#' @return A list with components `V` and/or `V2` containing estimated covariance matrices.
#' @export
#'
#' @examples
#' ### simulate some data with a linear trend
#' xx <- evd::rgpd(90*100, shape = 0.2) + 2*rep(1:100/100, each = 90)
#' ## temporal covariate for sliding BM
#' temp_cvrt <- rep(1:100/100, each = 90)[1:(90*100 - 90 +1)]
#'
#' ## compute sample of uniuqe sliding BM
#' bms <- get_uniq_bm(xx, 90, temp_cvrt = temp_cvrt, looplastblock = FALSE)
#'
#' ## full sample of sliding for estimaing the covariance matrix
#' slbm <- blockmax(xx, r = 90, "sliding")
#' estim <- fit_gev_univ(data = bms, type = "shift", hessian = TRUE)
#' est_var_univ_shift(slbm, est_par = estim, blcksz = 90, temp.cov = temp_cvrt,
#' varmeth = "both")
est_var_univ_shift <- function(orig_slbm, est_par, blcksz,  temp.cov =  NULL,
                               temp.cov.dj = NULL,
                               varmeth = "both"){
  nsl <- length(orig_slbm)
  k <-  nsl /blcksz
  k <- ifelse(k == floor(k), k, floor(k))


  if(is.null(est_par$hessian)) { stop("Hessian of log likelihood is missing.")}

  temp_cvrt_sl <- temp.cov[1:length(orig_slbm)]

  score.bdata <- score.function_univ_tilde(orig_slbm,
                                           theta =  est_par$mle, temp.cov = temp_cvrt_sl,
                                           type = "shift")

  temp_cvrt_sl <- temp.cov[1:length(orig_slbm)]


  fishest <- est_par$hessian/nsl

  fishestinv <- solve(fishest)

  useobs <- k*blcksz
  dimtheta <- 4

  if(!is.null(temp.cov.dj)) {
    meanbt <- mean(temp.cov.dj[1:k])

    meanbtsquare <- mean(temp.cov.dj[1:k]^2)
  }
  else {
    meanbt <- mean(temp_cvrt_sl)

    meanbtsquare <- mean(temp_cvrt_sl^2)
  }

  if(varmeth %in% c("V", "both")) {
    Y <- t(score.bdata[, 1:useobs])
    A2 <- array(Y, dim = c(blcksz, k, dimtheta))

    # result of purrr::map(..):
    # list of length blcksz. j-th list element is a matrix of dimension 4 x 4,
    # the empirical cross-covariance matrix of (score.fun(M_{r,1}^{(t)})_t and
    # (score.fun(M_{r,j}^{(t)})_t
    B2 <- 2*Reduce( "+" , purrr::map(1:blcksz, ~ cov(A2[1, , ], A2[.x, , ])) )/blcksz/k

    B2 <- (B2 + t(B2))/2


    B2[1:3, 4] <- B2[4,1:3 ] <-  meanbt*B2[1:3, 4]
    B2[4, 4] <- meanbtsquare*B2[4,4]

    V <- fishestinv %*% B2 %*% fishestinv
  }
  if(varmeth %in% c("V2", "both")) {

    Y <- score.bdata[, 1:useobs]
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

    v2shape <- 2* mean(as.numeric(lapply(AllDiags(cm), mean)) )/ k
    v2scale <- 2* mean(as.numeric(lapply(AllDiags(cmsc), mean)) )/ k
    v2loc <- 2* mean(as.numeric(lapply(AllDiags(cmloc), mean)) )/ k
    v2trend <- 2* mean(as.numeric(lapply(AllDiags(cmtrend), mean)) )/ k

    # use AllDiags2 for the mixed terms to make it symmetric
    v2shsc <- 2* mean(as.numeric(lapply(AllDiags2(cmshsc), mean)) )/ k
    v2shloc <- 2* mean(as.numeric(lapply(AllDiags2(cmshloc), mean)) )/ k
    v2shtrend <- 2* mean(as.numeric(lapply(AllDiags2(cmshtrend), mean)) )/ k

    v2scloc <- 2* mean(as.numeric(lapply(AllDiags2(cmscloc), mean)) )/ k
    v2sctrend <- 2* mean(as.numeric(lapply(AllDiags2(cmsctrend), mean)) )/ k
    v2loctrend <- 2* mean(as.numeric(lapply(AllDiags2(cmloctrend), mean)) )/ k


    B2 <- matrix( c(v2loc, v2scloc, v2shloc, v2loctrend,
                    v2scloc, v2scale, v2shsc, v2sctrend,
                    v2shloc, v2shsc, v2shape, v2shtrend,
                    v2loctrend, v2sctrend, v2shtrend, v2trend),
                  nrow = 4)

    B2[1:3, 4] <- B2[4,1:3 ] <-  meanbt*B2[1:3, 4]
    B2[4, 4] <- meanbtsquare*B2[4,4]

    V2 <- fishestinv %*% B2 %*% fishestinv

    rownames(V2) <- colnames(V2) <- c("mu", "sigma", "shape", "alpha")


  }


  if(varmeth == "both") {return(list(V = V, V2 = V2))}
  if(varmeth == "V"){return(list(V = V))}
  if(varmeth == "V2"){return(list(V2 = V2))}

}

# schneller als est_var_univ. Aktuell nur V variante
est_var_univ_faster <- function(orig_slbm, est_par, blcksz,  temp.cov =  NULL,
                                type = "shift",
                                varmeth = "both"){
  nsl <- length(orig_slbm)
  k <-  nsl /blcksz
  k <- ifelse(k == floor(k), k, floor(k))

  if(!(type %in% c("scale", "shift", "stationary"))) {
    stop("Type must be one of 'scale', 'shift', 'stationary'.")
  }
  if(is.null(est_par$hessian)) { stop("Hessian of log likelihood is missing.")}


  temp_cvrt_sl <- temp.cov[1:length(orig_slbm)]

  score.bdata <- score.function_univ(orig_slbm, theta =  est_par$mle, temp.cov = temp_cvrt_sl,
                                     type = type)

  fishest <- est_par$hessian/nsl

  fishestinv <- solve(fishest)

  useobs <- k*blcksz
  dimtheta <- ifelse(type == "stationary", 3, 4)

  Y <- t(score.bdata[, 1:useobs])

  A2 <- array(Y, dim = c(blcksz, k, dimtheta))


  B2 <- 2*Reduce( "+" , purrr::map(1:blcksz, ~ cov(A2[1, , ], A2[.x, , ])) )/blcksz/k

  B2 <- (B2 + t(B2))/2
  V <- fishestinv %*% B2 %*% fishestinv

  if(varmeth == "both") {return(list(V = V, V2 = V2))}
  if(varmeth == "V"){return(list(V = V))}
  if(varmeth == "V2"){return(list(V2 = V2))}

}


score.function_univ <- function(x, theta, temp.cov = NULL, type = "shift") {

  if(!(type %in% c("shift", "scale", "stationary"))){
    stop("Type must be either 'shift', 'scale' or 'stationary'.")}
  if(!(type == "stationary") & is.null(temp.cov)) {stop("Temporal covariate must be provided.")}
  if(type == "shift") {

    mu0 <- theta["mu"]
    mu1 <- theta["alpha"]
    sigma <- theta["sigma"]
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
    scoretemploc  <- scoreloc0*temp.cov
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

    mu0 <- theta["mu"]
    alpha0 <- theta["alpha"]
    sigma0 <- theta["sigma"]
    xi <- theta["shape"]

    if(alpha0 == 0 & mu0 == 0) {
      mut <- 0
      sigmat <- sigma0
    } else {
      mut <- mu0*exp(alpha0*temp.cov/mu0)
      sigmat <-  sigma0*exp(alpha0*temp.cov/mu0)
    }


    zt <- (x - mut)/sigmat
    in_supp <- which(1+ xi*zt > 0)
    zt <- zt[in_supp]

    if(abs(xi) < 1e-08) {
      ut <- exp(-zt)
    }
    else {
      ut <- (1+xi *zt)^(-1/xi)
    }

    # expo <- exp(alpha0*temp.cov/mu0)

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


score.function_univ_tilde <- function(x, theta, temp.cov = NULL, type = "shift") {

  if(!(type %in% c("shift", "scale", "stationary"))){
    stop("Type must be either 'shift', 'scale' or 'stationary'.")}
  if(!(type == "stationary") & is.null(temp.cov)) {stop("Temporal covariate must be provided.")}
  if(type == "shift") {

    mu0 <- theta["mu"]
    mu1 <- theta["alpha"]
    sigma <- theta["sigma"]
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

    scorescale <- ((1-u)*z -1)/(sigma*(1+xi*z))

    if(abs(xi) < 1e-08){
      scoreshape <- (1-u)*z^2/2- z
    }
    else{
      scoreshape <- (1-u)*1/xi*(1/xi*log(1+xi*z)-z/(1+xi*z)) - z/(1+xi*z)
    }

    return(matrix(c(scoreloc0, scorescale, scoreshape, scoreloc0),
                  nrow = 4, byrow = TRUE))
  }


}


