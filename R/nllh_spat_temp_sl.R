
#' Negative (pseudo) log-Likelihood for spatial GEV model with covariates
#'
#' Computes the negative (pseudo) log-Likelihood for the spatial GEV model with
#' parameters
#' dependending on spatial (and/or temporal) covariates.
#'
#' @param params Named vector of parameters for which the negative log-Likelihood
#' is to be computed. See details.
#' @param loc.sp.form R formula definining the spatial model for the location parameter.
#' @param scale.sp.form R formula definining the spatial model for the scale parameter.
#' @param loc.temp.form R formula definining the temporal trend for the location parameter.
#' @param scale.temp.form R formula definining the temporal trend for the scale parameter.
#' @param dataprep Prepared data of sliding block maxima and design matrices
#' as obtained from \code{\link[slbm]{prep4spatml}}.
#' @param spat.cov A data frame containing the spatial covariates used in the formulations
#' of the spatial formulas.
#' @param scale.link The link function that is used for modelling the scale parameter.
#'
#' @return The value of the negative log-Likelihood (a single numeric value).
#' @export
#'
#' @examples
#'
#' data("ExampleData")
#' data("GMST")
#' tempcv <- GMST %>% dplyr::filter(Year %in% c(1980:2019))
#' tempcv <- rep(tempcv$smoothedGMST, each = 90)[1:(39*90 +1)]
#' tempcv <- data.frame(GMST = tempcv)
#' set.seed(3)
#' spatial_cvrt <- data.frame(lat  = seq(0, 8, length = 8),
#'    lon = runif(8), ele = runif(8))
#' yy <- data.frame(ExampleData) %>%
#'       tidyr::pivot_longer( 1:8, names_to = "Station", values_to = "Obs")
#' ##############################
#' bms <- get_uniq_bm(yy, 90, temp_cvrt = tempcv$GMST)
#' prepdata <- prep4spatml(loc.sp.form = ~ lon + lat + ele,
#'               scale.sp.form = ~ lon +lat,
#'               loc.temp.form = ~ GMST, data = bms, spat.cov = spatial_cvrt)
#'
#' params  <- c("loc0" = 18, "loc1" = -1, "loc2" = 0.4, "loc3" = 0.02, "scale0" =0.2,
#'                "scale1" = 0.1, "scale2" = 0.1, "shape" = 0.3, "tempLoc1" = 1)
#'
#' nll_spat_temp_sl_prep(params = params, loc.sp.form = ~ lon + lat + ele,
#'           scale.sp.form = ~ lon + lat, loc.temp.form = ~ GMST,
#'           dataprep = prepdata, spat.cov = spatial_cvrt,
#'           scale.link = make.link("log"))
#' @importFrom stats lm model.frame model.matrix na.pass reformulate update
#'
nll_spat_temp_sl_prep <- function(params, loc.sp.form, scale.sp.form,
                                  loc.temp.form = NULL, scale.temp.form = NULL,
                                  dataprep, spat.cov, scale.link = make.link("log")){


  params.sp <- list(loc = params[startsWith(names(params),"loc")],
                    scale = params[startsWith(names(params),"scale")],
                    shape = params[startsWith(names(params),"shape")])

  params.temp <- list(loc = params[startsWith(names(params),"tempLoc")],
                      scale = params[startsWith(names(params),"tempScale")])

  params.temp <- purrr::map(params.temp, function(.x) ifelse(length(.x) == 0, 0, .x) )

  model.loc.sp <- dataprep$MatLocSp

  model.scale.sp <- dataprep$MatScaleSp

  loc_sp <- data.frame(loc_sp = model.loc.sp %*% params.sp$loc)

  scale_sp <-  data.frame(scale_sp = model.scale.sp %*% params.sp$scale)

  dataprep1 <- dataprep$data

  dataprep1 <- dataprep1 %>% dplyr::bind_cols( loc_sp , scale_sp )



  if(!is.null(loc.temp.form)) {
    loc.temp.form <- update(loc.temp.form,  ~ . + 0)
    dataprep1 <- dataprep1 %>%
      dplyr::mutate(LocTemp =  purrr::map(MatLocTemp,
                                          ~ as.numeric( .x %*% params.temp$loc) )
      ) %>%
      dplyr::select( - MatLocTemp)

    locpars <- purrr::map2(.x = dataprep1$LocTemp, .y = dataprep1$loc_sp, ~ .x + .y)


  } else {
    locpars <- dataprep1$loc_sp
  }


  if(!is.null(scale.temp.form)) {
    scale.temp.form <- update(scale.temp.form,  ~ . + 0)

    dataprep1 <- dataprep1 %>%
      dplyr::mutate(ScaleTemp =  purrr::map(MatScaleTemp,
                                            ~ as.numeric( .x %*% params.temp$scale) )
      ) %>%
      dplyr::select( - MatScaleTemp)

    scalepars <- purrr::map2(.x = dataprep1$ScaleTemp, .y = dataprep1$scale_sp, ~ .x + .y)


  } else {
    scalepars <- dataprep1$scale_sp
  }

  yy <- purrr::map(dataprep1$uniq_data, ~ .x$slbm)

  link.scalepars <- scale.link$linkinv(scalepars)



  if(abs(params.sp$shape ) < 1e-8) {
    yy <-  purrr::pmap(list(.x = yy, .y = locpars, .z = link.scalepars),
                       function(.x, .y, .z) {
                         exp(- (.x - .y)/.z)
                       } )
  } else {
    yy <- purrr::pmap(list(.x = yy, .y = locpars, .z = link.scalepars),
                      function(.x, .y, .z) {
                        (1 + params.sp$shape*( (.x - .y)/.z))^(-1/params.sp$shape)
                      })
  }

  stopyn <-  (any(unlist(purrr::map(yy, ~ any(.x < 0, na.rm = TRUE)))) |
                any(unlist(link.scalepars ) <= 0) |
                (any(unlist(purrr::map(yy, ~ is.nan(.x))))))

  if(stopyn) {
    return( 1e+10)
  } else {


    yy <- purrr::pmap( list(.x = yy, .y = link.scalepars , .z = dataprep1$uniq_data),
                       .f = function(.x, .y, .z) {
                         u <- (log(.y) - (params.sp$shape +1)*log(.x) + .x)
                         u*.z$n
                       }

    )

    sum(unlist(yy), na.rm = TRUE)

  }




}

#' Negative (pseudo) log-Likelihood for spatial GEV model with covariates under
#' homogeneity constraint.
#'
#' Computes the negative (pseudo) log-Likelihood for the spatial GEV model with
#' parameters dependending on spatial (and/or temporal) covariates under
#' the assumption of constant dispersion (location/scale) parameter.
#'
#' @inheritParams nll_spat_temp_sl_prep
#'
#' @return The value of the negative log-Likelihood (a single numeric value).
#' @export
#'
#' @examples
#' data("ExampleData")
#' yy <- data.frame(ExampleData) %>%
#'       tidyr::pivot_longer( 1:8, names_to = "Station", values_to = "Obs")
#' set.seed(3)
#' spatial_cvrt <- data.frame(lat  = seq(0, 8, length = 8),
#'    lon = runif(8), ele = runif(8))
#' ##############################
#' bms <- get_uniq_bm(yy, 90)
#' prepdata <- prep4spatml(loc.sp.form = ~ lon + lat,
#'               scale.sp.form = ~ lon +lat, data = bms, spat.cov = spatial_cvrt)
#'
#' params  <- c("disp" = 10, "scale0" = 1,
#'                "scale1" = -0.5, "scale2" = 0.1,  "shape" = -0.1)
#'
#' nll_spat_temp_sl_hom(params = params,
#'           scale.sp.form = ~ lon + lat + ele,
#'           dataprep = prepdata, spat.cov = spatial_cvrt,
#'            scale.link = make.link("identity"))
#'
nll_spat_temp_sl_hom <- function(params, scale.sp.form,
                                 scale.temp.form = NULL,
                                 dataprep, spat.cov, scale.link = make.link("log")){


  params.sp <- list(disp = params[startsWith(names(params),"disp")],
                    scale = params[startsWith(names(params),"scale")],
                    shape = params[startsWith(names(params),"shape")])

  params.temp <- list(scale = params[startsWith(names(params),"tempScale")])

  params.temp <- purrr::map(params.temp, function(.x) ifelse(length(.x) == 0, 0, .x) )

  model.scale.sp <- dataprep$MatScaleSp

  #  oc_sp <- data.frame(loc_sp = model.loc.sp %*% params.sp$loc)

  scale_sp <-  data.frame(scale_sp = model.scale.sp %*% params.sp$scale)

  loc_sp <- params.sp$disp*scale_sp
  names(loc_sp) <- "loc_sp"
  dataprep <- dataprep$data

  dataprep <- dataprep %>% dplyr::bind_cols( loc_sp, scale_sp )



  if(!is.null(scale.temp.form)) {
    scale.temp.form <- update(scale.temp.form,  ~ . + 0)

    dataprep <- dataprep %>%
      dplyr::mutate(ScaleTemp =  purrr::map(MatScaleTemp,
                                            ~ as.numeric( .x %*% params.temp$scale) )
      ) %>%
      dplyr::select( - MatScaleTemp)

    scalepars <- purrr::map2(.x = dataprep$ScaleTemp, .y = dataprep$scale_sp, ~ .x + .y)


  } else {
    scalepars <- dataprep$scale_sp
  }

  yy <- purrr::map(dataprep$uniq_data, ~ .x$slbm)

  link.scalepars <- scale.link$linkinv(scalepars)



  if(abs(params.sp$shape ) < 1e-8) {
    yy <-  purrr::pmap(list(.x = yy, .z = link.scalepars),
                       function(.x,  .z) {
                         exp(- (.x /.z -  params.sp$disp))
                       } )
  } else {
    yy <- purrr::pmap(list(.x = yy, .z = link.scalepars),
                      function(.x, .z) {
                        (1 + params.sp$shape*( .x /.z -  params.sp$disp))^(-1/params.sp$shape)
                      })
  }

  stopyn <- (any(unlist(purrr::map(yy, ~ any(.x < 0, na.rm = TRUE)))) | any(unlist(link.scalepars )<= 0) )
  if(stopyn) {
    return( 1e+10 )
  } else {

    yy <- purrr::pmap( list(.x = yy, .y = link.scalepars , .z = dataprep$uniq_data),
                       .f = function(.x, .y, .z) {
                         u <- (log(.y) - (params.sp$shape +1)*log(.x) + .x)
                         u*.z$n
                       }

    )

    sum(unlist(yy))
  }




}



