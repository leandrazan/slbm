# origslbm <- apply(ExampleData, 2, blockmax, r = 90, "sliding")
# origslbm <- reshape2::melt(origslbm)
# origslbm <- rename(origslbm, "Obs" = "value", "Station" = "Var2") %>% select(-Var1)
# params <-  fit_spgev_sl(data = bmuniq, loc.sp.form = ~ lon + lat,
#                         scale.sp.form = ~ lon +lat,
#                         loc.temp.form = ~ GMST, spat.cov = spatial_cvrt,
#                         datastart = djbm, st_val_meth = "LeastSqTemp",
#                         temp.cov = tempcv$smoothedGMST, return_hess = TRUE)$mle
#
# scorefun_spat(params= params, loc.sp.form = ~ lon + lat,scale.sp.form = ~ lon +lat,
#               spat.cov = spatial_cvrt, origslbm = origslbm)

scorefun_spat <- function(params, loc.sp.form, scale.sp.form,
                                  loc.temp.form = NULL, scale.temp.form = NULL,
                                  origslbm, spat.cov, scale.link = make.link("identity")){


  params.sp <- list(loc = params[startsWith(names(params),"loc")],
                    scale = params[startsWith(names(params),"scale")],
                    shape = params[startsWith(names(params),"shape")])

  params.temp <- list(loc = params[startsWith(names(params),"tempLoc")],
                      scale = params[startsWith(names(params),"tempScale")])

  params.temp <- purrr::map(params.temp, function(.x) ifelse(length(.x) == 0, 0, .x) )

  # model matrices
  model.loc.sp <- model.matrix(loc.sp.form, model.frame(loc.sp.form,
                                                        data = spat.cov, na.action = na.pass))

  model.scale.sp <- model.matrix(scale.sp.form, model.frame(scale.sp.form,
                                                            data = spat.cov, na.action = na.pass))

  ### to be implemented
   if(!is.null(loc.temp.form)) {
    loc.temp.form <- update(loc.temp.form,  ~ . + 0)

    data <- data %>%
      dplyr::mutate(MatLocTemp =  purrr::map(uniq_data, function(.x){

        datatemp <-  data.frame(.x$temp_cvrt)
        colnames(datatemp) <- all.vars(loc.temp.form)
        model.matrix(loc.temp.form,
                     model.frame(loc.temp.form,
                                 data = datatemp,
                                 na.action = na.pass))

      }
      ))

  }

  if(!is.null(scale.temp.form)) {
    scale.temp.form <- update(scale.temp.form,  ~ . + 0)


    data <- data %>%
      dplyr::mutate(MatScaleTemp =  purrr::map(uniq_data, function(.x){
        datatemp <-  data.frame(.x$temp_cvrt)
        colnames(datatemp) <- all.vars(scale.temp.form)
        model.matrix(scale.temp.form,
                     model.frame(scale.temp.form,
                                 data = datatemp,
                                 na.action = na.pass))

      }
      ))
  }


   # values of paramters
  loc_sp <- data.frame(loc_sp = model.loc.sp %*% params.sp$loc)

  scale_sp <-  data.frame(scale_sp = scale.link$linkinv(model.scale.sp %*% params.sp$scale))

  n.loc.sp <- ncol(model.loc.sp)
  n.scale.sp <- ncol(model.scale.sp)

  xx <- origslbm
  xx <- xx %>% group_by(Station) %>% nest(data = Obs) %>% mutate(data = purrr::map(data, ~ .x$Obs)) %>%
    ungroup()

  xx <- xx %>% bind_cols(loc_sp, scale_sp)

  # mutate value of (1+gamma*(x - mud)/sigmad)^(-1/gamma)
  xx <- xx %>% mutate(zd = purrr::pmap(list(data, loc_sp, scale_sp),
                                       function(data, loc_sp, scale_sp) {
                                           zz <- (data - loc_sp)/scale_sp
                                           zz[ (1+params.sp$shape*zz) <=  0] <- NA
                                           zz
                                         }),
                      ud = purrr::map(zd,
                                 function(x) {
                                   if(abs(params.sp$shape) < 1e-08) {
                                     exp( - x)
                                   }
                                   else {
                                    (1+ params.sp$shape*x)^(-1/params.sp$shape)

                                   }
                                 }))

  colnames(model.loc.sp) <- paste0("loc", 1:n.loc.sp)

  yy <- xx %>% dplyr::bind_cols(model.loc.sp)
  # mutate values of derivative of d/dmu_i loglik evaluated at observations with corresponding
  #parameters
  yy <- yy %>% dplyr::mutate_at(vars(paste0("loc", 1:n.loc.sp)),
                     ~ purrr::pmap(list(ud,  scale_sp, .x),
                                    function(ud,  scale_sp, .x) {
                                      ud^(params.sp$shape +1) *.x /scale_sp *
                         ( (params.sp$shape + 1)/ud -1) } )
    )
  # finally the sum over the d stations of the derivatives w.r.t. location
  locabl <- (apply(yy[ ,paste0("loc", 1:n.loc.sp) ], 2, function(x) Reduce("+", x )))

  # check:
  # aa <- xx[1, ]$ud[[1]]^(params.sp$shape +1)*1/
  #   xx[1,]$scale_sp[[1]]*( (params.sp$shape +1)/xx[1,]$ud[[1]] -1)

  ## same for scale derivative
   ## use xx from before
  colnames(model.scale.sp) <- paste0("scale", 1:n.scale.sp)

  yy <- xx %>% dplyr::bind_cols(model.scale.sp)

  yy <- yy %>% dplyr::mutate_at(vars(paste0("scale", 1:n.scale.sp)),
                                ~ purrr::pmap(list(zd, ud, loc_sp, scale_sp, .x),
                                              function(zd, ud, loc_sp,  scale_sp, .x) {
                                                ud^(params.sp$shape +1) *.x /scale_sp * zd*
                                                  ( (params.sp$shape + 1)/ud -1 - 1/scale_sp) } ))


## check:
# aa <- xx[1, ]$ud[[1]]^(params.sp$shape +1)*model.scale.sp[1, 2]/
#   xx[1,]$scale_sp[[1]]^2*( xx[1, ]$data[[1]] - xx[1,]$loc_sp[[1]]) *
#  ((params.sp$shape +1)/xx[1,]$ud[[1]] -1 - 1/xx[1,]$scale_sp[[1]])

  # finally the sum over the d stations of the derivatives w.r.t. scale
  scaleabl <- (apply(yy[ ,paste0("scale", 1:n.scale.sp) ], 2, function(x) Reduce("+", x )))


  ## for derivative wrt shape
  yy <- xx %>% dplyr::mutate(shapeabl = purrr::pmap(list(zd, ud, scale_sp, loc_sp),
                                              function(zd, ud, scale_sp, loc_sp) {
                                                xi <- params.sp$shape

                                                if(abs(params.sp$shape) < 1e-08) {
                                                  (1-ud)*zd^2/2- zd
                                                }
                                                else {
                                                  (1-ud)*1/xi*(1/xi*log(1+xi*zd)-zd/(1+xi*zd)) - zd/(1+xi*zd)
                                                }
                                              }))


 shapeabl <- Reduce("+", yy$shapeabl)

 dermat <- (cbind(locabl, scaleabl, shapeabl))/nrow(xx)

 dermat

}




compute_sl_cov <- function(Y, column1, column2, varmeth, blcksz) {

  nsl <- ncol(Y)
  k <-  nsl /blcksz
  k <- ifelse(k == floor(k), k, floor(k))

  useobs <- k*blcksz

  Y1 <- Y[column1, 1:useobs ]
  Y2 <- Y[column2, 1:useobs ]

  m1 <- t(matrix(Y1, nrow = blcksz, ncol = k))
  m2 <- t(matrix(Y2, nrow = blcksz, ncol = k))

  cm1m2 <- cov(m1, m2)
  if(varmeth == "V2") {
    return( 2* mean(as.numeric(lapply(AllDiags(cm1m2), mean, na.rm = TRUE)) )/ k)
  }
  else {
    return(2*mean(cm1m2[1, ], na.rm = TRUE)/k)
  }
}

### jetzt auf alle Kombination von Spalten die sliding Kovarianz schätzung anwenden !

compute_sl_cov_spat <- function(origslbm, params, hessmat, blcksz, varmeth,
                                loc.sp.form, scale.sp.form, spatial_cvrt) {

  scorevec <- scorefun_spat(params = params, loc.sp.form = loc.sp.form, scale.sp.form = scale.sp.form,
                             spat.cov = spatial_cvrt, origslbm = origslbm)

   hessmat <- hessmat/nrow(origslbm)
   fishestinv <- solve(hessmat)

   Y <- (fishestinv %*% t(scorevec))

   nrowy <- nrow(Y)
   covsl <- array(dim = c(nrowy, nrowy))
   for(j in 1:nrowy) {
     for( l in j:nrowy) {
       covsl[j, l] <- covsl[l,j] <- compute_sl_cov(Y, j, l, blcksz = 90, varmeth = "V2")
     }
   }

   colnames(covsl) <- rownames(covsl) <- colnames(fishestinv)
   covsl
 }


# compute_sl_cov_spat(aa, hessmat = hessmat$hessian,
#                     blcksz = 90, varmeth = "V2")

#
# hessmat <-  fit_spgev_sl(data = bmuniq, loc.sp.form = ~ lon + lat,
#                          scale.sp.form = ~ lon +lat,
#                          loc.temp.form = NULL, spat.cov = spatial_cvrt,
#                          datastart = djbm, st_val_meth = "LeastSq",
#                          temp.cov = tempcv$smoothedGMST, return_hess = TRUE)
#
#
# compute_sl_cov_spat(aa, hessmat = hessmat$hessian, k = 40, blcksz = 90, varmeth = "V2")
#
