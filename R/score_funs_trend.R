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
score.fun <- function(x, theta, temp.cov = NULL, type = "shift", chain = TRUE,
                                ...) {

  add.args <- list(...)
  if(!(type %in% c("shift", "scale", "stationary"))){
    stop("Type must be either 'shift', 'scale' or 'stationary'.")}
  if(!(type == "stationary") & is.null(temp.cov)) {stop("Temporal covariate must be provided.")}

  if(chain) {
    if(type == "shift") {

      mu0 <- theta["mu"]
      mu1 <- theta["alpha"]
      sigma <- theta["sigma"]
      xi <- theta["shape"]

      mu <- mu0 + mu1*temp.cov

      z <- (x-mu)/sigma
      in_supp <- which( (1+ xi*z) > 0)
      z[!in_supp] <- NA

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

      scoremat <- matrix(c(scoreloc, scorescale, scoreshape),
                         nrow = 3, byrow = TRUE)


      return(scoremat)

    }
    if(type == "scale") {
      # whether trend parameter is set relative to location parameter
      rel_trend <- add.args$rel_trend

      mu0 <- theta["mu"]
      alpha0 <- theta["alpha"]
      sigma0 <- theta["sigma"]
      xi <- theta["shape"]

      if(rel_trend) {
        mut <- mu0*exp(alpha0*temp.cov/mu0)
        sigmat <-  sigma0*exp(alpha0*temp.cov/mu0)
      }
      else {
        mut <- mu0*exp(alpha0*temp.cov)
        sigmat <- sigma0*exp(alpha0*temp.cov)
      }

      zt <- (x - mut)/sigmat
      in_supp <- which(1+ xi*zt > 0)
      zt[!in_supp] <- NA

      if(abs(xi) < 1e-08) {
        ut <- exp(-zt)
      }
      else {
        ut <- (1+xi *zt)^(-1/xi)
      }

      scoreloc <- (xi + 1 - ut)/(sigmat*(1+xi*zt))
      scorescale <- ((1-ut)*zt -1)/(sigmat*(1+xi*zt))

      if(abs(xi) < 1e-08){
        scoreshape <- (1-ut)*zt^2/2- zt
      }
      else{
        scoreshape <- (1-ut)*1/xi*(1/xi*log(1+xi*zt)-zt/(1+xi*zt)) - zt/(1+xi*zt)
      }

      scoremat <- matrix(c(scoreloc, scorescale, scoreshape),
                         nrow = 3, byrow = TRUE)

      return(scoremat)

    }
  }
  if(type == "shift") {

    mu0 <- theta["mu"]
    mu1 <- theta["alpha"]
    sigma <- theta["sigma"]
    xi <- theta["shape"]

    mu <- mu0 + mu1*temp.cov

    z <- (x-mu)/sigma
    in_supp <- which( (1+ xi*z) > 0)
    z[!in_supp] <- NA
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

    if(!add.args$rel_trend) {

      mut <- mu0*exp(alpha0*temp.cov)
      sigmat <-  sigma0*exp(alpha0*temp.cov)

      zt <- (x - mut)/sigmat
      in_supp <- which(1+ xi*zt > 0)
      zt[!in_supp] <- NA

      if(abs(xi) < 1e-08) {
        ut <- exp(-zt)
      }
      else {
        ut <- (1+xi *zt)^(-1/xi)
      }

      # expo <- exp(alpha0*temp.cov/mu0)

      dmu0ut <- ut^(xi +1)/sigma0

      # ut^(xi +1)/sigmat* (expo*(1 - alpha0*temp.cov/mu0) - (x - mut)*alpha0*temp.cov/mu0)
      dsigma0ut <- ut^(xi +1)*zt/sigma0

      dalpha0ut <- ut^(xi +1)*(zt + mu0/sigma0)*temp.cov

      scoreloc0 <-  dmu0ut*( (xi +1)/ut -1)

      scorescale <- -1/sigma0 + dsigma0ut*( (xi +1)/ut -1)

      scorealpha <- -temp.cov + dalpha0ut*( (xi +1)/ut  -1)

      if(abs(xi) < 1e-08){
        scoreshape <- (1-ut)*zt^2/2- zt
      }
      else{
        scoreshape <- (1-ut)*1/xi*(1/xi*log(1+xi*zt)-zt/(1+xi*zt)) - zt/(1+xi*zt)
      }

      return(matrix(c(scoreloc0, scorescale, scoreshape, scorealpha),
                    nrow = 4, byrow = TRUE))

    }
    else {
      mut <- mu0*exp(alpha0*temp.cov/mu0)
      sigmat <-  sigma0*exp(alpha0*temp.cov/mu0)

      zt <- (x - mut)/sigmat
      in_supp <- which(1+ xi*zt > 0)
      zt[!in_supp] <- NA

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
    z[!in_supp] <- NA
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


compute_cov_stat <- function(Y, varmeth, k, useobs, blcksz) {

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

    # variances
    vshape <- 2*mean(cm[1,]) / k
    vscale <-  2*mean(cmsc[1,]) / k
    vloc <-  2*mean(cmloc[1,]) / k

    # covariances
    vshsc <- 2*mean(cmshsc[1,]) / k
    vshloc <- 2*mean(cmshloc[1,]) / k
    vscloc <- 2*mean(cmscloc[1,]) / k
    vscsh <- 2*mean(cmshsc[,1]) / k
    vlocsh <- 2*mean(cmshloc[,1]) / k
    vlocsc <- 2*mean(cmscloc[,1]) / k


    V <- matrix( c(vloc, vscloc, vshloc,
                   vlocsc, vscale, vshsc,
                   vlocsc, vscsh, vshape),
                 nrow = 3)
    V <- (V + t(V))/2
    rownames(V) <- colnames(V) <- c("loc", "scale", "shape")

  }

  if(varmeth %in% c("both", "V2")) {
    # variances. Use AllDiags
    v2shape <- 2* mean(as.numeric(lapply(AllDiags(cm), mean)) )/ k
    v2scale <- 2* mean(as.numeric(lapply(AllDiags(cmsc), mean)) )/ k
    v2loc <- 2* mean(as.numeric(lapply(AllDiags(cmloc), mean)) )/ k
    # covariances. Use AllDiags2
    v2shsc <- 2* mean(as.numeric(lapply(AllDiags2(cmshsc), mean)) )/ k
    v2shloc <- 2* mean(as.numeric(lapply(AllDiags2(cmshloc), mean)) )/ k
    v2scloc <- 2* mean(as.numeric(lapply(AllDiags2(cmscloc), mean)) )/ k


    V2 <- matrix( c(v2loc, v2scloc, v2shloc,
                    v2scloc, v2scale, v2shsc,
                    v2shloc, v2shsc, v2shape),
                  nrow = 3)
    rownames(V2) <- colnames(V2) <- c("loc", "scale", "shape")
  }

  if(varmeth == "both") {return(list(V = V, V2 = V2))}
  if(varmeth == "V"){return(list(V = V))}
  if(varmeth == "V2"){return(list(V2 = V2))}

}


compute_cov_nonstat <- function(covstat, type, temp_cvrt, Jinv, ...) {

  add.args <- list(...)

  if(type == "shift") {
    meanct <- mean(temp_cvrt, na.rm = TRUE)
    meanctsq <- mean(temp_cvrt^2, na.rm = TRUE)


    covstat2 <- covstat
    covstat2 <- cbind(covstat, meanct*covstat[1,])
    covstat2 <- rbind(covstat2, c(meanct*covstat[1,], meanctsq*covstat[1,1]))

    return(Jinv %*% covstat2 %*% Jinv)
  }
  if(type == "scale") {
    rel_trend <- add.args$rel_trend
    param <- add.args$param

    if(!rel_trend) {

      alpha0 <- param["alpha"]
      sigma0 <- param["sigma"]
      mu0 <- param["mu"]

      meanexp_ac <- mean(exp(alpha0*temp_cvrt), na.rm = TRUE)
      meanexp_ac_c <- mean(temp_cvrt*exp(alpha0*temp_cvrt), na.rm = TRUE)
      meanexp_ac_sq <- mean(exp(alpha0*temp_cvrt)^2, na.rm = TRUE)
      meanexp_ac_c_sq <- mean(temp_cvrt*exp(alpha0*temp_cvrt)^2, na.rm = TRUE)
      meanexp_ac_csq_sq <- mean(temp_cvrt^2*exp(alpha0*temp_cvrt)^2, na.rm = TRUE)


      covstat2 <- covstat
      covstat2[1:2, 1:2] <- meanexp_ac_sq*covstat[1:2, 1:2]
      covstat2[3, 1:2] <- meanexp_ac*covstat[3, 1:2]
      covstat2[1:2, 3] <- meanexp_ac*covstat[1:2, 3]

      alpha_ab <- c(meanexp_ac_c_sq*(mu0*covstat[1,1] + sigma0*covstat[1,2]),
                    meanexp_ac_c_sq*(mu0*covstat[1,2] + sigma0*covstat[2,2]),
                    meanexp_ac_c*(mu0*covstat[1,3] + sigma0*covstat[2,3])
                    )

      sig44 <- meanexp_ac_csq_sq*(mu0*(mu0*covstat[1,1] + sigma0*covstat[1,2]) +
                                  sigma0*(mu0*covstat[1,1] + sigma0*covstat[2,2]))

      covstat2 <- rbind(covstat2, alpha_ab)
      covstat2 <- cbind(covstat2, c(alpha_ab,sig44))

      return( (Jinv %*% covstat2 %*% Jinv) )
    }
    else{
      print("to do")}

  }

}

# xx <- dailyPrec %>% filter(lubridate::month(Date) %in% 10:12, lubridate::year(Date) >= 1950)
# gmst <- GMST %>% filter(Year >= 1950, Year <= 2010)
# bms <- get_uniq_bm(xx$cumPrec, 92, temp_cvrt = rep(gmst$smoothedGMST, each = 92) , looplastblock = FALSE)
#
# # # full sample of sliding for estimaing the covariance matrix
# slbm <- blockmax(xx$cumPrec, r = 92, "sliding")
# estim <- fit_gev_univ(data = bms, type = "scale", hessian = TRUE, rel_trend = FALSE)
#
# estim
#
# bb <- score.fun(slbm, theta = estim$mle, temp.cov = rep(gmst$smoothedGMST, each = 92)[1:length(slbm)],
#                 type = "scale", chain = TRUE, rel_trend = FALSE)
# est_var_chain(orig_slbm = slbm, est_par =  estim, blcksz =  blcksz,
#               temp.cov = rep(gmst$smoothedGMST, each = 92)[1:length(slbm)],
#               type = "scale", varmeth = "both", rel_trend = FALSE)

#
# est_var_univ(slbm, est_par = estim, blcksz = 90, temp.cov = temp_cvrt,
#  varmeth = "both", type = "scale")
#
# est_var_univ_shift(slbm, est_par = estim, blcksz = blcksz, temp.cov = temp_cvrt, temp.cov.dj = temp_cvrt)
# est_var_chain(slbm, estim, blcksz, temp_cvrt, "shift", "both")
# est_var_chain(slbm, estim, blcksz, temp_cvrt, "scale", "both", rel_trend = FALSE)


est_var_chain <- function(orig_slbm, est_par, blcksz,  temp.cov =  NULL,
                         type = "shift",
                         varmeth = "both", ...){

  add.args <- list(...)

  nsl <- length(orig_slbm)
  k <-  nsl /blcksz
  k <- ifelse(k == floor(k), k, floor(k))
  useobs <- k*blcksz

  if(!(type %in% c("scale", "shift"))) {
    stop("Type must be one of 'scale', 'shift'.")
  }
  if(is.null(est_par$hessian)) { stop("Hessian of log likelihood is missing.")}


  temp_cvrt_sl <- temp.cov[1:length(orig_slbm)]

  Y <- score.fun(orig_slbm, theta =  est_par$mle, temp.cov = temp_cvrt_sl,
                           chain = TRUE,
                                     type = type, rel_trend = add.args$rel_trend)

  fishest <- est_par$hessian/nsl

  fishestinv <- solve(fishest)

  # list with covariance matrices
  Gammahat <- compute_cov_stat(Y = Y, varmeth = varmeth, k = k, useobs = useobs, blcksz = blcksz)

  purrr::map(Gammahat, ~ compute_cov_nonstat(.x, type = type, temp_cvrt = temp_cvrt_sl,
                                             Jinv = fishestinv,
                                             rel_trend = add.args$rel_trend,
                                             param = est_par$mle))



}
