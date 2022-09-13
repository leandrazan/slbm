### score functions

# compute the score function at a given parameter vector for either the stationary model,
# a linear shift or a trend that scales the observations


# passt auch für xi = 0


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
      if(is.null(rel_trend)) {stop("Please specify the parametrisation of the scale model via 'rel_trend'.")}
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
  else {
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


      if(is.null(add.args$rel_trend)) {
        stop("Please specify the parametrisation of the scale model via 'rel_trend'.")
        }


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


}


AllDiags <- function(inmat, sorted = TRUE) {
  Range <- ncol(inmat) - 1
  Range <- 0:Range
  if (isTRUE(sorted)) Range <- Range[order(abs(Range))]
  lapply(Range, function(x) {
    inmat[row(inmat) == (col(inmat) - x)]
  })
}
AllDiags2 <- function(inmat, sorted = TRUE) {
  Range <- ncol(inmat) - 1
  Range <- 0:Range
  if (isTRUE(sorted)) Range <- Range[order(abs(Range))]
  lapply(Range, function(x) {
    c(inmat[row(inmat) == (col(inmat) - x)], inmat[col(inmat) == (row(inmat) - x)])
  })
}

compute_cov_stat_4chain <- function(Y, varmeth, k, useobs, blcksz) {

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


compute_cov_nonstat_chain <- function(covstat, type, temp_cvrt, Jinv, ...) {

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
                                    sigma0*(mu0*covstat[1,2] + sigma0*covstat[2,2]))

      covstat2 <- rbind(covstat2, alpha_ab)
      covstat2 <- cbind(covstat2, c(alpha_ab,sig44))

      return( (Jinv %*% covstat2 %*% Jinv) )
    }
    else{

      alpha0 <- param["alpha"]
      sigma0 <- param["sigma"]
      mu0 <- param["mu"]

      meanexp_ac <- mean(exp(alpha0/mu0*temp_cvrt), na.rm = TRUE)
      meanexp_ac_c <- mean(temp_cvrt*exp(alpha0/mu0*temp_cvrt), na.rm = TRUE)
      meanexp_ac_sq <- mean(exp(alpha0/mu0*temp_cvrt)^2, na.rm = TRUE)
      meanexp_ac_c_sq <- mean(temp_cvrt*exp(alpha0/mu0*temp_cvrt)^2, na.rm = TRUE)
      meanexp_ac_csq_sq <- mean(temp_cvrt^2*exp(alpha0/mu0*temp_cvrt)^2, na.rm = TRUE)

      a1 <-  (1-alpha0/mu0*temp_cvrt)

      a2 <- sigma0/mu0

      sigma11 <- exp(alpha0/mu0*temp_cvrt)^2 *
        (a1* (a1*covstat[1,1] - a2*alpha0*temp_cvrt/mu0*covstat[1,2]) -
           a2*alpha0*temp_cvrt/mu0 * (a1*covstat[1,2] - a2*alpha0*temp_cvrt/mu0 *covstat[2,2]))

      sigma11 <- mean(sigma11, na.rm = TRUE)

      sigma12 <- meanexp_ac_sq *covstat[1,2] - meanexp_ac_c_sq *(alpha0/mu0*covstat[1,2] +
                                                                   sigma0*alpha0/mu0^2*covstat[2,2])

      sigma13 <- meanexp_ac * covstat[1,3] - meanexp_ac_c*(a2*covstat[1,3] - a2*alpha0/mu0*covstat[2,3])

      sigma14 <- meanexp_ac_c_sq *(covstat[1,1] + a2*covstat[1,2]) -
        meanexp_ac_csq_sq*alpha0/mu0*( covstat[1,1] + 2*a2*covstat[1,2] + a2*covstat[2,2])

      sigma22 <- meanexp_ac_sq*covstat[2,2]

      sigma23 <- meanexp_ac*covstat[2,3]

      sigma24 <- meanexp_ac_c_sq *(covstat[1,2] + a2*covstat[2,2])

      sigma33 <- covstat[3,3]

      sigma34 <- meanexp_ac_c*(covstat[1,3] + a2*covstat[2,3] )

      sigma44 <- meanexp_ac_csq_sq *(covstat[1,1] + 2*a2*covstat[1,2] +a2^2*covstat[2,2])

      covstat2 <- matrix(c(sigma11, sigma12, sigma13, sigma14,
                           sigma12, sigma22, sigma23, sigma24,
                           sigma13, sigma23, sigma33, sigma34,
                           sigma14, sigma24, sigma34, sigma44), ncol = 4)

      return( (Jinv %*% covstat2 %*% Jinv) )

    }

  }

}


est_var_univ_nochain <- function(orig_slbm, est_par, blcksz,  temp.cov =  NULL,
                         type = "shift",
                           varmeth = "both", ...){


  nsl <- length(orig_slbm)
  k <-  nsl /blcksz
  k <- ifelse(k == floor(k), k, floor(k))

  if(!(type %in% c("scale", "shift", "stationary"))) {
    stop("Type must be one of 'scale', 'shift', 'stationary'.")
  }
  if(is.null(est_par$hessian)) { stop("Hessian of log likelihood is missing.")}


  temp_cvrt_sl <- temp.cov[1:length(orig_slbm)]

  score.bdata <- score.fun(orig_slbm, theta =  est_par$mle, temp.cov = temp_cvrt_sl,
                                        type = type, chain = FALSE, ... = ...)
  # score.bdata <- score.function_univ(orig_slbm, theta =  est_par$mle, temp.cov = temp_cvrt_sl,
  #                                      type = type)

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

      vscsh <- 2*mean(cmshsc[,1]) / k
      vlocsh <- 2*mean(cmshloc[,1]) / k
      vtrendsh <- 2*mean(cmshtrend[,1]) / k
      vlocsc <- 2*mean(cmscloc[,1]) / k
      vtrendsc <- 2*mean(cmsctrend[,1]) / k
      vtrendloc <- 2*mean(cmloctrend[,1]) / k

      V <- matrix( c(vloc, vscloc, vshloc, vloctrend,
                     vlocsc, vscale, vshsc, vsctrend,
                     vlocsh, vscsh, vshape, vshtrend,
                     vtrendloc, vtrendsc, vtrendsh, vtrend),
                   nrow = 4)
      V <- (V + t(V))/2
      # V <- matrix( c(vloc, vscloc, vshloc, vloctrend,
      #                vscloc, vscale, vshsc, vsctrend,
      #                vshloc, vshsc, vshape, vshtrend,
      #                vloctrend, vsctrend, vshtrend, vtrend),
      #              nrow = 4)
      rownames(V) <- colnames(V) <- c("mu", "sigma", "shape", "alpha")

    }

    if(varmeth %in% c("both", "V2")) {

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


      V2 <- matrix( c(v2loc, v2scloc, v2shloc, v2loctrend,
                      v2scloc, v2scale, v2shsc, v2sctrend,
                      v2shloc, v2shsc, v2shape, v2shtrend,
                      v2loctrend, v2sctrend, v2shtrend, v2trend),
                    nrow = 4)
      rownames(V2) <- colnames(V2) <- c("mu", "sigma", "shape", "alpha")
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

  }

  if(varmeth == "both") {return(list(V = V, V2 = V2))}
  if(varmeth == "V"){return(list(V = V))}
  if(varmeth == "V2"){return(list(V2 = V2))}

}


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
  Gammahat <- compute_cov_stat_4chain(Y = Y, varmeth = varmeth, k = k, useobs = useobs, blcksz = blcksz)

  purrr::map(Gammahat, ~ compute_cov_nonstat_chain(.x, type = type, temp_cvrt = temp_cvrt_sl,
                                                   Jinv = fishestinv,
                                                   rel_trend = add.args$rel_trend,
                                                   param = est_par$mle))



}

#' GEV parameter covariance matrix for sliding blocks ML estimator
#'
#' @param orig_slbm Full sample of sliding block maxima
#' @param est_par Output of \code{\link{fit_gev_univ}}.
#' @param blcksz The blocksize parameter.
#' @param temp.cov The temporal covariate (must be same length as \code{orig_slbm}).
#' @param type Either 'shift', 'scale' or 'stationary', depending if one assumes
#' that obervations shift, scale or do not change w.r.t. the covariate or stationary
#' (for more details, see
#'  \code{\link[slbm]{fit_gev_univ}}).
#' @param varmeth One of 'V', 'V2' or 'both'. Determines the method used for estimating
#' the covariance matrix.
#' @param chain logical; whether the estimator for the covariance matrix should be one obtained
#'  from the chain rule. Only relvant for the non-stationary models.
#'
#' @param ... A further argument that needs to be passed when `type = scale`, specifying the
#' parametrisation of the model.
#' @return A list with components `V` and/or `V2` containing estimated covariance matrices.
#' @export
#'
#' @examples
#' ### simulate some data with a linear trend
# xx <- evd::rgpd(90*100, shape = 0.2) + 2*rep(1:100/100, each = 90)
# ## temporal covariate for sliding BM
# temp_cvrt <- rep(1:100/100, each = 90)[1:(90*100 - 90 +1)]
#
# ## compute sample of uniuqe sliding BM
# bms <- get_uniq_bm(xx, 90, temp_cvrt = temp_cvrt, looplastblock = FALSE)
#
# ## full sample of sliding for estimaing the covariance matrix
# slbm <- blockmax(xx, r = 90, "sliding")
# estim <- fit_gev_univ(data = bms, type = "shift", hessian = TRUE)
# est_var_univ(slbm, est_par = estim, blcksz = 90, temp.cov = temp_cvrt,
# varmeth = "both", type = "shift", chain = TRUE)
est_var_univ <-  function(orig_slbm, est_par, blcksz,  temp.cov =  NULL,
                          type = "shift",
                          varmeth = "both", chain = TRUE, ...) {
  if(type == "stationary" | !chain) {
    est_var_univ_nochain(orig_slbm = orig_slbm, est_par = est_par,
                         blcksz = blcksz, temp.cov = temp.cov, type = type, varmeth = varmeth, ... = ...)
  }
  else {
    est_var_chain(orig_slbm = orig_slbm, est_par = est_par,
                  blcksz = blcksz, temp.cov = temp.cov, type = type, varmeth = varmeth,
                  ... = ...)
  }
}



qdelta_rl <- function(theta, Tyrl, type, ref_gmst = NULL, rel_trend = TRUE) {

  ct <- -log(1-1/Tyrl)

  if(type == "shift") {
    if(is.null(ref_gmst)) { stop("Reference covariate value must be provided.")}
    if(is.null(names(theta))) { names(theta) <- c("mu", "sigma", "shape", "alpha")}
    mu0 <- theta["mu"]
    mu1 <- theta["alpha"]
    sigma <- theta["sigma"]
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
      c("mu" = dmu0, "sigma" = unname(dsigma), "shape" = unname(dxi), "alpha" = dmu1))

  }
  if(type == "scale") {
    if(is.null(ref_gmst)) { stop("Reference covariate value must be provided.")}
    if(is.null(names(theta))) { names(theta) <- c("mu", "sigma", "shape", "alpha")}

    if(rel_trend) {
      mu0 <- theta["mu"]
      alpha0 <- theta["alpha"]
      sigma0 <- theta["sigma"]
      xi <- theta["shape"]

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
    }
    else {
      mu0 <- theta["mu"]
      alpha0 <- theta["alpha"]
      sigma0 <- theta["sigma"]
      xi <- theta["shape"]

      expo <- exp(alpha0*ref_gmst)

      if(abs(xi) < 1e-08) {
        dmu0 <- expo
        dsigma0 <- expo*log(ct)
        dxi <- sigma0*expo*log(ct)^2/2
        dalpha <- expo*ref_gmst*( sigma0 * (-log(ct)) +1)
      }
      else {

        dmu0 <- expo

        dsigma0 <- expo*(ct^(-xi) -1)/xi

        dxi <- sigma0*expo/xi^2*(1-ct^(-xi)*(xi*log(ct) +1))

        dalpha <- expo*ref_gmst*( sigma0 * (ct^(-xi) -1)/xi +1)
      }

    }

    return(
      c("mu" = unname(dmu0), "sigma" = unname(dsigma0), "shape" = unname(dxi),
      "alpha" = unname(dalpha))
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
#' @param rel_trend logical; only relevant when `type = "scale"`. Specifies the parametrisation
#' in the scale model.
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
estimate_var_rl <- function(theta, Tyrl, type, ref_gmst = NULL, Covmat,
                            rel_trend = TRUE) {

  if(type == "stationary") {
    ref_gmst <- 1
  }


  qdelt <- purrr::map(ref_gmst, ~
                        qdelta_rl(theta = theta, Tyrl = Tyrl, type = type, ref_gmst = .x,
                                  rel_trend = rel_trend)
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
