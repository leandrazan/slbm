#' Compute initial values for fitting a spatio-temporal GEV model.
#'
#' Compute initial values for fitting a spatio-temporal GEV model.
#'
#'
#' @param data The data matrix from which initial GEV parameters are estimated.
#' @param loc.sp.form R formula definining the spatial model for the location parameter.
#' @param scale.sp.form R formula definining the spatial model for the scale parameter.
#' @param loc.temp.form R formula definining the temporal trend for the location parameter.
#' @param scale.temp.form R formula definining the temporal trend for the scale parameter.
#' @param spat.cov A data frame containing the spatial covariates used in the formulations
#' of the spatial formulas.
#' @param temp.cov A data frame (or matrix) with named columns containing the
#' temporal covariates used in the formulations
#' of the temporal formulas.
#' @param method The method used for computing initial values: One of "LeastSq",
#' "LeastSqTemp", "spatialGev", "fgevMean". See details for more information.
#'
#' @param print_start_vals Logical; whether to print the initial parameters to the console.
#' @param scale.link The link function that is used for modelling the scale parameter.
#' @param ... Further agruments that can be passed. Atm, the only supported further
#' argument is \code{type = "IF"} for computing initial values under homogeneity
#' constraint (Indexflood assumption). See details for further information.
#'
#' @details
#'
#'
#' The different methods supplied in `methods` work as follows.
#' \describe{
#' \item{LeastSq}{The three GEV parameters are estimated stationwise (with the
#' function \code{\link[evd]{fgev}} of the `evd` package).
#' Then, linear models with model equations as given by the \emph{spatial} formulas
#' are fitted for location and scale. The estimated model coefficients are returned
#' as respective initial values. For the initial vale of the shape parameter, the average
#' over stationwise shape parameters is returned.
#' If any formula for a trend in location and/or scale parameter is given, the
#' corresponding initial coefficients are set to 0.
#'
#' The link function passed in `scale.link` is used for modelling
#' the scale parameter.
#' }
#' \item{LeastSqTemp}{This method computes initial values for the temporal location and
#' the spatial
#' coefficients, so the `loc.temp.form`
#' and the `temp.cov`
#' arguments should not be `NULL`.
#' GEV parameters and the location trend coefficient are estimated stationwise,
#' with the function \code{\link[evd]{fgev}} of the `evd` package and by passing
#' the temporal covariate to its `nsloc`-argument).
#'
#' The initial spatial coefficients and shape are then computed as described in the `LeastSq`-method,
#' while
#' the mean of the estimated stationwise location trend coefficients
#' is taken as the initial location trend coefficient. For the scale parameter,
#' the link function passed in `scale.link` is used for modelling.
#' Currently, with this method, a trend in the scale parameter is not considered.
#' }
#' \item{spatialGEV}{
#' This method computes initial valus with the function
#' \code{fitspatgev} of the `spatialExtremes` package.
#' With this method, initial coefficients for all spatial and temporal covariates
#' can be computed.
#' Since \code{fitspatgev} does not support link functions,
#' the use of any other link function than the identitiy function
#' when the scale parameter depends on any covariate is not recommended.
#' }
#' \item{fgevMean}{This is the most basic and least accurate method of the possible
#' methods.
#' The initial values for the constant parts of the parameters (the intercepts)
#' are the averages of stationwise GEV parameter estimates. Coefficients
#' concerning any of the covariates are set to 0.
#' }
#' }
#'
#' If the additional argument `type = IF` is set,
#' initial parameters for the GEV model under homogeneity constraint are computed.
#' That is, the ratio of location and scale parameter is assumed to be constant in
#' space and time. In this case, the linear model of the location parameter
#' must be a multiple \eqn{\delta} of the scale parameter's linear model.
#'  This
#' \eqn{\delta} is a parameter of the model.
#' Initial values are computed for \eqn{\delta} as well as the coefficients of the
#' formula passed in `scale.form` (`loc.form` can
#' be `NULL` or anything else, since it will not be used).
#' First, parameter values \eqn{(\delta_i, \sigma_i, \gamma_i)} are estimated stationwise.
#' Then, the average values of
#' \eqn{\delta_i} and \eqn{\gamma_i}
#' are returned as initial values of \eqn{\delta} and \eqn{\gamma}, respectivelym
#' while the estimated coefficients
#' of the linear model fitted to
#' the estimated scale parameter \eqn{\sigma_i} are returned as initial values for
#' the spatial coefficients.
#'
#' @return A vector containing a set of initial parameter values.
#' @export
#'
#' @examples
#' ### prepare data
#' data("ExampleData")
#' data("GMST")
#'
#' djbm <- apply(ExampleData, 2, blockmax, r = 90, "disjoint")
#'
#' set.seed(3)
#' spatial_cvrt <- data.frame(lat  = seq(0, 8, length = 8),
#'    lon = runif(8), ele = runif(8))
#' tempcov <- GMST$smoothedGMST[101:140]
#'
#' ### compute initial values based on disjoint block maxima
#' get_start_vals(djbm, loc.sp.form = ~ lon + lat, scale.sp.form = ~ lon + lat,
#' loc.temp.form = ~ GMST, scale.link = make.link("identity"),
#' spat.cov = spatial_cvrt, temp.cov = tempcov, method = "spatialGev")
#'
get_start_vals <- function(data, loc.sp.form = ~ 1, scale.sp.form = ~ 1,
                           loc.temp.form = NULL, scale.temp.form = NULL,
                           spat.cov, temp.cov = NULL, method = "LeastSq",
                           print_start_vals = TRUE, scale.link = make.link("log"),
                           ...){


  if(!is.matrix(data)) {data <- as.matrix(data)}

  if( !(method %in% c("LeastSq", "LeastSqTemp", "spatialGev", "fgevMean"))) {
    stop("This is not a valid method for computing starting values.")
  }

  n.loc.temp <- length(all.vars(loc.temp.form))
  n.scale.temp <- length(all.vars(scale.temp.form))
  n.loc.sp <- length(all.vars(loc.sp.form)) + 1
  n.scale.sp <- length(all.vars(scale.sp.form)) + 1

  n.param <- n.loc.sp + n.scale.sp + 1

  IFargs <- list(...)

  if(exists("type", where = IFargs) & !is.null(IFargs$type)) {

    #browser()
    YY <- fgev_hom(data, method = "BFGS", scale.link = scale.link)$mle

    .datfr <- data.frame(t(YY)) %>% dplyr::bind_cols(spat.cov)

    lm.sigma <- reformulate(deparse(scale.sp.form[[2]]), paste0(scale.link$name, "(scale)"))

    st_disp <- mean(.datfr$loc/.datfr$scale)
    names(st_disp) <- "disp"
    st_sigmas <- lm(lm.sigma, data = .datfr)$coefficients
    names(st_sigmas) <-  paste0("scale", 0:(length(st_sigmas)-1))
    st_shape <- mean(.datfr$shape)
    names(st_shape) <- "shape"

    if(is.null(scale.temp.form)){
      start_vals <- c(st_disp, st_sigmas, st_shape)

    }
    else{

      n.scale.temp <- length(all.vars(scale.temp.form))

      st_temp_scale <- rep(0, n.scale.temp)
      names(st_temp_scale) <- paste0("tempScale", 1:n.scale.temp)

      start_vals <- c(st_disp, st_sigmas, st_shape, st_temp_scale)
    }

  }


  # sets coefficients of temporal covariates to 0
  if(method == "LeastSq") {

      YY <- apply(data, 2,  function(.x) evd::fgev(.x, std.err = F)$estimate)

      .datfr <- data.frame(t(YY)) %>% dplyr::bind_cols(spat.cov)


      lm.mu <- reformulate(deparse(loc.sp.form[[2]]), "loc")
      lm.sigma <- reformulate(deparse(scale.sp.form[[2]]), paste0(scale.link$name, "(scale)"))

      st_mus <- lm(lm.mu, data = .datfr )$coefficients
      names(st_mus) <- paste0("loc", 0:(length(st_mus)-1))

      st_sigmas <- lm(lm.sigma, data = .datfr)$coefficients
      names(st_sigmas) <-  paste0("scale", 0:(length(st_sigmas)-1))

      st_shape <- mean(.datfr$shape)
      names(st_shape) <- "shape"

      if(n.loc.temp == 0 & n.scale.temp == 0) {
        start_vals <- c(st_mus, st_sigmas, st_shape)
      }
      else { # initialise zeros for temporal coefficients

        st_temp_loc <- rep(0, n.loc.temp)
        if(n.loc.temp > 0){
          names(st_temp_loc) <- paste0("tempLoc", 1:n.loc.temp)
        }

        st_temp_scale <- rep(0, n.scale.temp)
        if(n.scale.temp > 0 ){
          names(st_temp_scale) <- paste0("tempScale", 1:n.scale.temp)
        }

        start_vals <- c(st_mus, st_sigmas, st_shape, st_temp_loc, st_temp_scale)
      }

  }

  # estimates location trend based on mean of stationwise location trends
  # (one covariate only)
  if(method == "LeastSqTemp") {

    if(is.null(temp.cov)) {
      stop("Please provide a temporal covariate of the same length as the data.")
    }
    else{
      temp.cov <- data.frame(temp.cov)
      if(!(nrow(temp.cov) == nrow(data)) | ncol(temp.cov) > 1) {
        stop("Please provide a temporal one-dimensional covariate of the same
             length as the data.")
      }

      YY <- apply(data, 2,
                  function(.x) evd::fgev(.x, nsloc = temp.cov,
                                         std.err = F)$estimate)

      .datfr <- data.frame(t(YY)) %>% dplyr::bind_cols(spat.cov)


      lm.mu <- reformulate(deparse(loc.sp.form[[2]]), "loc")
      # lm.sigma <- reformulate(deparse(scale.sp.form[[2]]), "log(scale)")
      lm.sigma <- reformulate(deparse(scale.sp.form[[2]]), paste0(scale.link$name, "(scale)"))

      st_mus <- lm(lm.mu, data = .datfr )$coefficients
      names(st_mus) <- paste0("loc", 0:(length(st_mus)-1))

      st_sigmas <- lm(lm.sigma, data = .datfr)$coefficients
      names(st_sigmas) <-  paste0("scale", 0:(length(st_sigmas)-1))

      st_shape <- mean(.datfr$shape)
      names(st_shape) <- "shape"

      # initialise zeros for temporal coefficients

      st_temp_loc <- mean(.datfr[, 2])
      names(st_temp_loc) <- "tempLoc1"


      st_temp_scale <- rep(0, n.scale.temp)
      if(n.scale.temp > 0 ){
        names(st_temp_scale) <- paste0("tempScale", 1:n.scale.temp)
      }

      start_vals <- c(st_mus, st_sigmas, st_shape, st_temp_loc, st_temp_scale)


    }
  }

  # can estimate arbitrary amount of (trend) coefficients.
  # No link function for scale parameter
  if(method == "spatialGev"){

    if (is.null(temp.cov) & any( c(n.loc.temp, n.scale.temp) > 0)) {

      warning("No temporal covariate was supplied to compute starting values with
             this method. Initial Temporal parameters are set to 0.")

      start_params <- SpatialExtremes::fitspatgev(data = data,
                                                  covariables = as.matrix(spat.cov),
                                                  loc.form =  loc.sp.form,
                                                  scale.form = scale.sp.form,
                                                  shape.form = ~ 1)$fitted.values
    } else {
      temp.cov <- as.matrix(temp.cov)
      if(is.null(colnames(temp.cov))) {
        colnames(temp.cov) <- all.vars(loc.temp.form)
      }

      start_params <- SpatialExtremes::fitspatgev(data = data,
                                                  covariables = as.matrix(spat.cov),
                                                  loc.form =  loc.sp.form,
                                                  scale.form = scale.sp.form,
                                                  shape.form = ~ 1,
                                                  temp.form.loc = loc.temp.form,
                                                  temp.form.scale = scale.temp.form,
                                                  temp.cov = temp.cov)$fitted.values
    }

    st_param <- list(loc = start_params[startsWith(names(start_params),"loc")],
                     scale = start_params[startsWith(names(start_params),"scale")],
                     shape = start_params[startsWith(names(start_params),"shape")])


    start_vals <- numeric(n.param)
    start_vals[1:n.loc.sp] <- st_param$loc

    # caution: SpatialExtremes doesn't use link functionv for scale parameter,
    # method best to use without link function (identity) or when scale is constant
    start_vals[(n.loc.sp +1):(n.param -1)] <- scale.link$linkfun(st_param$scale)

    start_vals[n.param] <- st_param$shape
    names(start_vals) <- c(paste0("loc", 0:(n.loc.sp - 1)),
                           paste0("scale", 0:(n.scale.sp - 1)), "shape")

    if(n.loc.temp  > 0) {
      sv_temp_loc <- numeric(n.loc.temp)
      if(!(is.null(temp.cov))) {
        sv_temp_loc <- start_params[startsWith(names(start_params),"tempCoeffLoc")]
      }
      names(sv_temp_loc) <- paste0("tempLoc", 1:n.loc.temp)
      start_vals <- c(start_vals, sv_temp_loc)
    }
    if(n.scale.temp  > 0) {
      sv_temp_scale <- numeric(n.scale.temp)
      if(!(is.null(temp.cov))) {
        sv_temp_scale <- start_params[startsWith(names(start_params),"tempCoeffScale")]
      }
      names(sv_temp_scale) <- paste0("tempScale", 1:n.scale.temp)
      start_vals <- c(start_vals, sv_temp_scale)
    }
  }

  # most basic method. Estimates only stationary parameters, other coefficients are 0
  if(method == "fgevMean"){

    start_st <-  apply(data,2, function(.x){evd::fgev(.x[!is.na(.x)],  std.err = F)})
    start_st <- purrr::map(purrr::map_dfr(start_st, function(x) x$estimate), mean)
    start_vals <- numeric(n.param)
    start_vals[1] <- start_st$loc
    start_vals[n.loc.sp + 1] <- scale.link$linkfun(start_st$scale)
    start_vals[n.param] <- start_st$shape
    names(start_vals) <- c(paste0("loc", 0:(n.loc.sp-1)),
                           paste0("scale", 0:(n.scale.sp -1)), "shape")


    if(n.loc.temp  > 0) {
      sv_temp_loc <- rep(0, n.loc.temp)
      names(sv_temp_loc) <- paste0("tempLoc", 1:n.loc.temp)
      start_vals <- c(start_vals, sv_temp_loc)
    }
    if(n.scale.temp  > 0) {
      sv_temp_scale <- rep(0, n.scale.temp)
      names(sv_temp_scale) <- paste0("tempScale", 1:n.scale.temp)
      start_vals <- c(start_vals, sv_temp_scale)
    }
  }

  if(print_start_vals){
    message( "start values are", "\n")
    print(c(  round(start_vals,3) ))
  }

  start_vals

}


