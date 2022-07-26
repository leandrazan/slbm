# data is tibble of unique pairs of sliding bm and temporal covariate with
# respective frequency of appearance, as computed by get_unique_bm()

## use st_val_meth = sp_gev only when scale doesnt depend on covariate, since then, start values
## are computed with SpatialExtremes package that does not use a link function

## ... can be the following arguments: maxit or similar for optim,
## temp.cov for computing starting values with SpatialExtremes when depending on temporal cov.
#  type = IF for homogeneity assumption

#' Fit spatial-temporal GEV model based on sliding block maxima
#'
#' @param data A nested tibble containing the unique sliding block maxima and their
#' frequency for the considered locations as obtained from the function \link{get_uniq_bm}.
#' @param loc.sp.form R formula definining the spatial model for the location parameter.
#' @param scale.sp.form R formula definining the spatial model for the scale parameter.
#' @param loc.temp.form R formula definining the temporal trend for the location parameter.
#' @param scale.temp.form R formula definining the temporal trend for the scale parameter.
#' @param spat.cov A data frame containing the spatial covariates used in the formulations
#' of the spatial formulas.
#' @param start_vals Optional; a vector containing initial values for the log Likelihood
#' which is to be optimised. If omitted, initial values are computed with \code{\link{get_start_vals()}}.
#' @param datastart The data matrix on which starting values will be computed if
#' none are passed in `start_vals`.
#' @param use_gr Logical; whether to use theoretical gradient function when optimising.
#' Not working yet.
#' @param method The optimisation method that is passed to `optim()`. Defaults to "BFGS".
#' @param st_val_meth Method passed to `get_start_vals()` for cpmputing starting values if none
#' are provided.
#' @param print_start_vals Logical; whether to print out the initial parameter values.
#' @param scale.link The link function that is used to model the scale parameter.
#' Defaults to \eqn{log(\sigma) = \sigma_0 + \sigma_1 + ... }.
#' @param ... Further arguments. These can be
#' \describe{
#' \item{temp.cov}{A numeric vector or matrix containing temporal covariates for the
#' data provided in `datastart` when needed to compute initial values. For more details,
#' see the documentation of \code{\link{get_start_vals}}.}
#' \item{type}{
#' If you set `type  = "IF"`, coefficients will be fitted under the assumption
#' of constant dispersion parameter. In this case, any formulas passed concerning the location
#' parameter will be ignored and only the formulas concerning the scale parameter will be considered. }
#' \item{maxit, reltol ... }{Furhter components that are passed as a list to the
#' `control` argument of `optim()`.}
#' }
#'
#' @return Returns a list containing the following components:
#' \describe{
#' \item{mle}{A numeric vector containing the estimated parameter values.}
#' \item{nllh}{The negative log-Likelihood evaluated at the estimated parameters.}
#' \item{conv}{The convergence code, inherited from `optim`.}
#' \item{counts}{A two-element integer vector giving the number of calls to the negative
#' log-Likelihood and its gradient, respectively, as obtained from optim.}
#' }
#' @export
#'
#' @examples
#'
#' ######## prepare data ##########################################
#' data("ExampleData")
#' data("GMST")
#' tempcv <- GMST %>% dplyr::filter(Year %in% c(1980:2019))
#' tempcvsl <- rep(tempcv$smoothedGMST, each = 90)[1:(39*90 +1)]
#' tempcvsl <- data.frame(GMST = tempcvsl)
#' set.seed(3)
#' spatial_cvrt <- data.frame(lat  = seq(0, 8, length = 8),
#'    lon = runif(8), ele = runif(8))
#' yy <- data.frame(ExampleData) %>%
#'       tidyr::pivot_longer( 1:8, names_to = "Station", values_to = "Obs")
#' ##############################
#' bmuniq <- get_uniq_bm(yy, 90, temp_cvrt = tempcvsl$GMST, looplastblock = FALSE)
#' # disjoint block maxima are used for computing starting values:
#' djbm <- apply(ExampleData, 2, blockmax, r = 90, "disjoint")
#' fit_spgev_sl(data = bmuniq, loc.sp.form = ~ lon + lat,
#'               scale.sp.form = ~ lon +lat,
#'               loc.temp.form = ~ GMST, spat.cov = spatial_cvrt,
#'               datastart = djbm, st_val_meth = "LeastSqTemp",
#'               temp.cov = tempcv$smoothedGMST, return_hess = TRUE)
#'
#'
fit_spgev_sl <- function(data, loc.sp.form = ~ 1, scale.sp.form = ~ 1,
                         loc.temp.form = NULL, scale.temp.form = NULL,
                         spat.cov, start_vals = NULL, datastart = NULL, use_gr = FALSE,
                         method = "BFGS", st_val_meth = "LeastSq", print_start_vals = TRUE,
                         scale.link = make.link("log"),
                         return_hessian = FALSE, ...){

  add.args <- list(...)

  d.dat <- nrow(data)

  if(d.dat == 1 & missing(spat.cov)) {
    spat.cov <- data.frame(lon = 1)
  }

  n.spat <- nrow(spat.cov)

  if (d.dat != n.spat) {
    stop("'data' and 'covariates' doesn't match")
  }

  if (is.null(start_vals) & is.null(datastart)) {
    stop("Please provide either starting values or a dataset where starting values can
         be computed on.")
  }

  dataprep <- prep4spatml(loc.sp.form = loc.sp.form, scale.sp.form = scale.sp.form,
                          loc.temp.form = loc.temp.form, scale.temp.form = scale.temp.form,
                          data = data, spat.cov = spat.cov)




  if(is.null(start_vals)){
    if(!exists("temp.cov", where = add.args)) { add.args$temp.cov <- NULL }

    start_vals <-  get_start_vals(data = datastart,
                                  loc.sp.form = loc.sp.form,
                                  scale.sp.form = scale.sp.form,
                                  loc.temp.form = loc.temp.form,
                                  scale.temp.form = scale.temp.form,
                                  spat.cov = spat.cov,
                                  temp.cov = add.args$temp.cov,
                                  method = st_val_meth,
                                  print_start_vals = print_start_vals,
                                  scale.link = scale.link,
                                  type = add.args$type)


    add.args <- add.args[!(names(add.args) == "temp.cov")]

  }


  if(use_gr){
    mlest <- optim(start_vals, fn = nll_spat, gr =  gr_spat_gev,
                   loc.sp.form = loc.sp.form,
                   scale.sp.form = scale.sp.form,
                   loc.temp.form = loc.temp.form, scale.temp.form = scale.temp.form,
                   data = data, spat.cov = spat.cov, temp.cov = temp.cov,
                   method = method, control = add.args)
  }
  if(exists("type", where = add.args)) {

    add.args <- add.args[!(names(add.args) == "type")]

   mlest <- optim(start_vals, fn = nll_spat_temp_sl_hom,
                 scale.sp.form = scale.sp.form,
                 scale.temp.form = scale.temp.form,
                 dataprep = dataprep, spat.cov = spat.cov,
                 scale.link = scale.link,
                 method = method, control = add.args)
  }
  else{
    mlest <- optim(start_vals, fn = nll_spat_temp_sl_prep,
                   loc.sp.form = loc.sp.form,
                   scale.sp.form = scale.sp.form,
                   loc.temp.form = loc.temp.form, scale.temp.form = scale.temp.form,
                   dataprep = dataprep, spat.cov = spat.cov,
                   scale.link = scale.link,
                   method = method, control = add.args, hessian = return_hessian)
  }

  if(!(mlest$convergence == 0)) {
    warning("Optimization did not succeed. Try increasing the maximum number of iterations or try a different method for finding starting values.")
  }

  if(!return_hessian) {
    list(mle = mlest$par, nllh = mlest$value, conv = mlest$convergence, counts = mlest$counts)
  }
  else {
    list(mle = mlest$par, nllh = mlest$value, conv = mlest$convergence, counts = mlest$counts,
         hessian = mlest$hessian)
  }

}



# data("ExampleData")
# data("GMST")
# tempcv <- GMST %>% dplyr::filter(Year %in% c(1980:2019))
# tempcvsl <- rep(tempcv$smoothedGMST, each = 90)[1:(39*90 +1)]
# tempcvsl <- data.frame(GMST = tempcvsl, Datetime = rep(1:(90*40),
#  each = 40)[1:length(tempcvsl)])
# set.seed(3)
# spatial_cvrt <- data.frame(lat  = seq(0, 8, length = 8),
#    lon = runif(8), ele = runif(8))
# dates <- rep(1:(40), each = 90)
# yy <- data.frame(Datetime = dates) %>% bind_cols(data.frame(ExampleData))
#
# bmuniq <- get_uniq_bm(yy, 90, temp.cov = tempcvsl, tempvar = GMST)
#
# ##### disjoint block maxima are used for computing starting values:
# djbm <- apply(ExampleData, 2, blockmax, r = 90, "disjoint")
# fit_spgev_sl(data = bmuniq, loc.sp.form = ~ lon + lat,
#               scale.sp.form = ~ lon +lat,
#               loc.temp.form = ~ GMST, spat.cov = spatial_cvrt,
#               datastart = djbm, st_val_meth = "LeastSqTemp",
#               temp.cov = tempcv$smoothedGMST)
