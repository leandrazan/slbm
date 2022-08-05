pareto <- function(x, shape){
  if( any(x < 0) | any(x >1)){stop("x must be between 0 and 1 ")}
  return((1-x)^(-shape))
}


#' Simulate spatially dependent data with margins in MDA(GEV)
#'
#' @param n Number of observations at each location
#' @param locations Matrix with named columns giving the coordinates of the locations
#' (lon and lat)
#' @param margins A list with elements
#' \describe{
#' \item{dist}{The marginal distribution of the resulting data. Must be one of
#' "pareto", "gpd", "norm", "gev".}
#' \item{shape}{The shape parameter of the chosen marginal distribution.}
#' }
#' @param ms A list with components that are passed to \code{rmaxstab} from the
#'  \code{SpatialExtremes}-package:
#' \describe{
#' \item{kovmod}{The covariance model for the max-stable process.}
#' \item{nugget, range, smooth}{The nugget, range and smooth parameters.}
#' }
#'
#' @return Returns a matrix with numeric values of dimension \eqn{ n x d},
#'  where \eqn{d} is the number of
#' locations passed in `coords`. The dependence structure of the returned data is
#' that  of the chosen max-stable model, while the marginal distributions are as
#' chosen.
#'
#' @export
#'
#' @examples
#' coords <- cbind(runif(2), runif(2))
#' colnames(coords ) <- c("lat", "lon")
#' mdax <- MDAdata(100, locations = coords, margins = list(distr = "gpd", shape = -0.2),
#'     ms = list(kovmod = "whitmat", nugget = 0, range = 3, smooth = 0.44))
#' plot.ts(mdax[, 1])
#' lines(mdax[, 2], col = 2)
#'
#' @importFrom SpatialExtremes rmaxstab
MDAdata <- function(n, locations, margins = list(distr = "gpd", shape = 0.2),
                    ms = list(kovmod = "whitmat", nugget = 0, range = 3, smooth = 1)){

  # generate max stable process with unit frechet margins
  frechdata <- switch(ms$kovmod,
                      "whitmat" = SpatialExtremes::rmaxstab(n = n, coord = locations, cov.mod = ms$kovmod,
                                           range = ms$range, smooth = ms$smooth, nugget = ms$nugget),
                      "powexp" = SpatialExtremes::rmaxstab(n = n, coord = locations, cov.mod = ms$kovmod, nugget = ms$nugget,
                                          range = ms$range, smooth = ms$smooth),
                      "gauss" = SpatialExtremes::rmaxstab(n = n, coord = locations, cov.mod = ms$kovmod, cov11 = ms$cov11,
                                         cov12 = ms$cov12,cov22 = ms$cov22)
  )

  # select inverse distribution function to desired marginal distribution
  qdistr_func <- switch(margins$distr, "Pareto" = qpareto,
                        "gpd" = evd::qgpd,
                        "norm" = qnorm,
                        "gev" = evd::qgev)
  mdadata <- qdistr_func(exp(-1/frechdata), shape = margins$shape)

  return(mdadata)

}


#' Theoretical coefficients of linear GEV parameter model
#'
#' Computes approximate theoretical values of the coefficients of the linear models
#' of GEV parameters for the case where block maxima are of the form
#' \eqn{\max{\sigma(t,d)X_1 + \mu(t,d), ..., \sigma(t,d)X_r + \mu(t,d) }},
#' with \eqn{\sigma(t,d)} and \eqn{\mu(t,d)} linear functions of spatial and
#' temporal covariates, and \eqn{X_i} GPD distributed random variables.
#'
#' @param blcksz The blocksize $r$ for which block maxima are computed.
#' @param loc0 The constant coefficient in \eqn{\mu(t,d)}.
#' @param scale0 The constant coefficient in \eqn{\sigma(t,d)}.
#' @param loclat The coefficient in \eqn{\mu(t,d)} belonging to latitude.
#' @param scalelat The coefficient in \eqn{\sigma(t,d)} belonging to latitude.
#' @param loclon The coefficient in \eqn{\mu(t,d)} belonging to longitude.
#' @param scalelon The coefficient in \eqn{\sigma(t,d)} belonging to longitude.
#' @param locele The coefficient in \eqn{\mu(t,d)} belonging to elevation.
#' @param scaleele The coefficient in \eqn{\sigma(t,d)} belonging to elevaiton.
#' @param loctemp The coefficient in \eqn{\mu(t,d)} belonging to the temporal covariate.
#' @param scaletemp The coefficient in \eqn{\sigma(t,d)} belonging to the temporal covariate.
#' @param shape The shape parameter.
#'
#' @return A vector with approximate coefficient values (only those which are not 0).
#' @export
#'
#' @examples
#' theo_pars(90, loc0 = 10, scale0 = 1, scalelon = 0.2,  shape = -0.2)
theo_pars <- function(blcksz, loc0 = 0, scale0 = 1, loclat = 0, scalelat = 0, loclon = 0,
                      scalelon = 0,
                      locele = 0, scaleele = 0,
                      loctemp = 0, scaletemp = 0, shape = 0.2
                      ) {

  ar <- blcksz^shape
  br <- (ar -1)/shape

  locconst <- loc0 + br*scale0
  loclat <- loclat + br*scalelat
  loclon <- loclon + br*scalelon
  locele <- locele + br*scaleele

  scaleconst <- ar*scale0
  scalelat <- ar *scalelat
  scalelon <- ar*scalelon
  scaleele <- ar*scaleele

  loctemp <- loctemp + br*scaletemp
  scaletemp <- ar*scaletemp

  parvec <- c(locconst = locconst, loclat = loclat, loclon = loclon, locele = locele,
          scaleconst = scaleconst, scalelat = scalelat, scalelon = scalelon,
          scaleele = scaleele, loctemp = loctemp, scaletemp = scaletemp,
          shape = shape)

  parvec[!(parvec == 0)]

}

#' Simulate data from the MDA of a spatial GEV model
#'
#' Simulate data from the MDA of a spatial-temporal GEV model with the dependence
#' structure of a max-stable process.
#' The input data \eqn{X_1, X_2 ..., } is transformed through
#'
#' \eqn{ \sigma(t,d) X_i + \mu(t,d)},
#' where \eqn{\sigma(t,d), \mu(t,d)} are linear functions of spatial and temporal
#' covariates. Theoretical values of the resulting GEV parameter coefficients
#' in case the \eqn{X_i} are GPD distributed can be computed with
#' \code{\link{theo_pars}}.
#'
#' @param mdadat Data matrix with observations in the domain of attraction of the GEV
#' distribution, as generated from \code{\link{MDAdata}}.
#' @param params The vector of coefficients for the linear models of covariates.
#' @param spat.cov A data frame containing the spatial covariates used in the formulations
#' of the spatial formulas.
#' @param temp.cov A named matrix containing the temporal covariates used in the formulations
#' of the temporal formulas.
#' @param loc.sp.form R formula definining the spatial model for the location parameter.
#' @param scale.sp.form R formula definining the spatial model for the scale parameter.
#' @param loc.temp.form R formula definining the temporal trend for the location parameter.
#' @param scale.temp.form R formula definining the temporal trend for the scale parameter.

#' @return A matrix with same dimension \eqn{n x d} as input data,
#' transformed with the affine linear transformation
#' \eqn{ \sigma(t,d) X_i + \mu(t,d)}.
#'
#' @export
#'
#' @examples
#' data("GMST")
#' ny <- 100
#' blcksz <- 90
#' coords <- cbind(spatial_cvrt$lat, spatial_cvrt$lon)
#' colnames(coords ) <- c("lat", "lon")
#' mdax <- MDAdata(ny*blcksz, locations = coords, margins = list(distr = "gpd", shape = -0.2),
#'                 ms = list(kovmod = "whitmat", nugget = 0, range = 3, smooth = 0.44))
#'
#' tempcov1 <- data.frame(GMST = rep(GMST$smoothedGMST[1:100], each = 90))
#' xx <- sim_spatial(mdax, params = c("loc0" = 10, "loc1" = 2, "loc2" = 0.05,
#'                                    "scale0" = 2, "scale1" = 0.2, "shape" = -0.2,
#'                                    "tempLoc" = 2),
#'                   loc.sp.form = ~ lat + lon, scale.sp.form = ~ lon,
#'                   spat.cov = spatial_cvrt, temp.cov = tempcov1,
#'                   loc.temp.form = ~ GMST)
#'
#' ### plot the resulting disjoint block maxima
#' bmx <- apply(xx, 2, blockmax, r = 90,  "disjoint")
#' plot.ts(bmx)
#' # fit spatio-temporal GEV model
#' \dontrun{
#' tempcovdj <- matrix(GMST$smoothedGMST[1:100])
#' colnames(tempcovdj) <- "GMST"
#' SpatialExtremes::fitspatgev(bmx, covariables = as.matrix(spatial_cvrt),
#' loc.form = ~ lat + lon, scale.form = ~lon,
#' shape.form = ~ 1, temp.form.loc = ~ GMST,
#' temp.cov = tempcovdj)$fitted.values
#' }
#'
#' theo_pars(90, loc0 = 10, scale0 = 2, loclat = 2, scalelat = 0, loclon = 0.05,
#'  scalelon = 0.2,
#' loctemp = 2, shape = -0.2)
#'
sim_spatial <- function(mdadat, params, spat.cov, temp.cov = NULL,
                        loc.sp.form = ~1, scale.sp.form = ~1, loc.temp.form = NULL,
                        scale.temp.form = NULL) {

  params.sp <- list(loc = params[startsWith(names(params),"loc")],
                    scale = params[startsWith(names(params),"scale")],
                    shape = params[startsWith(names(params),"shape")])

  params.temp <- list(loc =  params[startsWith(names(params),"tempLoc")],
                      scale = params[startsWith(names(params),"tempScale")])

  n.loc.temp <- length(all.vars(loc.temp.form))

  n.scale.temp <- length(all.vars(scale.temp.form))

  model.loc.sp <- model.matrix(loc.sp.form, model.frame(loc.sp.form,
                                                         data = spat.cov, na.action = na.pass))

  model.scale.sp <- model.matrix(scale.sp.form, model.frame(scale.sp.form,
                                                            data = spat.cov, na.action = na.pass))

  scale_param <- (model.scale.sp %*% params.sp$scale)

  scale_param <- matrix(rep(scale_param, each = nrow(mdadat)), nrow = nrow(mdadat))

  loc_param <- model.loc.sp %*% params.sp$loc

  loc_param <-  matrix(rep(loc_param, each = nrow(mdadat)), nrow = nrow(mdadat))

  if(n.loc.temp > 0 ) {
    loc.temp.form <- update(loc.temp.form,  ~ . + 0)
    model.loc.temp <-  model.matrix(loc.temp.form, model.frame(loc.temp.form,
                                                           data = temp.cov, na.action = na.pass))

    loc_param.temp <- model.loc.temp %*% params.temp$loc
    loc_param.temp <-  matrix(rep(loc_param.temp, ncol(mdadat)), nrow = nrow(mdadat))

    loc_param <- loc_param + loc_param.temp

  }
  if(n.scale.temp > 0 ) {
   scale.temp.form <- update(scale.temp.form,  ~ . + 0)
   model.scale.temp <-  model.matrix(scale.temp.form, model.frame(scale.temp.form,
                                                             data = temp.cov, na.action = na.pass))
   scale_param.temp <- (model.scale.temp %*% params.temp$scale)
   scale_param.temp <-  matrix(rep(scale_param.temp, ncol(mdadat)), nrow = nrow(mdadat))

   scale_param <- scale_param + scale_param.temp
  }



  #browser()

  mdadat <- mdadat *scale_param + loc_param

}




# tcvrt <- tempcovdj
# cvrt_sl <- c(rep(tcvrt[1], blcksz/2), rep(tcvrt[2:(ny)], each = blcksz))
# cvrt_sl <- c(cvrt_sl, rep(tcvrt[ny +1], ny*blcksz - length(cvrt_sl) ))
#
#
# yy <- data.frame(xx) %>% pivot_longer( 1:8, names_to = "Station", values_to = "Obs")
# bms <- get_uniq_bm(yy, blcksz, cvrt_sl[1: (nrow(xx) - blcksz +1) ])
#
# system.time({
#   estim1 <- fit_spgev_sl(data = bms, loc.sp.form = ~ lon  + lat ,  scale.sp.form = ~ lon,
#                          loc.temp.form = ~ GMST, scale.temp.form = NULL,
#                          spat.cov = spatial_cvrt, use_gr = FALSE,
#                          method = "BFGS", st_val_meth = "LeastSqTemp",
#                          datastart = apply(xx, 2, blockmax, r = blcksz, "disjoint"),
#                          temp.cov  = tempcovdj,
#                          print_start_vals = TRUE,
#                          maxit  = 2000, scale.link = make.link("identity"))
# })
#
#
# estim1
