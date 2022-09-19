#' Compute stationwise unique sliding block maxima and their frequency
#'
#' Compute stationwise unique sliding block maxima and their respective
#' frequency from observations from several locations in space over time.
#' If the distribution is assumed to be non-stationary, the frequencies where
#' block maxima and temporal covariate coincide are computed.
#'
#' @param data A tibble with columns `Station` and `Obs`, where `Obs` contains
#' measurements of the variable of interest at the Station referenced by `Station`.
#' If dealing with observations from one station only (univariate case), it can also be a
#' vector containing the observations.
#' @param blcksz The blocksize for which block maxima are computed.
#' @param temp_cvrt Optional; a numeric vector of length `nrow(data) - blcksz +1`
#' containing values of a temporal
#' covariate. The specific length is necessary in order to match each
#' sliding block maximum to a value of the temproal
#' covariate. If too long, the vector will be cut off accordingly.
#'  To be provided if you assume the distribution of your
#'  data to be non-stationary and want to fit parameters accordingly.
#' @param looplastblock logical; whether or not the first disjoint block of observations is
#'  concatenated to the last disjoint block of observations when computing the sample of sliding
#'  block maxima. Because of stationarity issues, this should only be used when
#'  sample is assumed to be stationary.
#'
#' @return Returns a nested tibble with columns `Station` and `uniq_data`:
#' for each station, `uniq_data` contains a tibble with columns `slbm`, `temp_cvrt`
#' (if considered)
#' and `n`, with the unqiue values of sliding block maxima (and temporal covariate)
#' along with the respective frequency.
#' In the univariate case, a tibble containing the unique values and their frequency
#' (and, if considered, the corresponding value of the temporal covariate) is returned.
#' @export
#'
#' @examples
#' ##### generate some data #####
#' blcksz <- 90
#' xx <- sapply(c(0.1, 0.4), function(x) evd::rgpd(10*blcksz, shape = x))
#'
#' yy <- data.frame(xx) %>% tidyr::pivot_longer( 1:2, names_to = "Station", values_to = "Obs")
#' ##############################
#'
#' # define a temporal covariate that is constant over a block of length blcksz
#' temp_cvrt <- rep(1:10/10, each = blcksz)[1:(9*blcksz + 1)]
#'
#' ##### apply function #####
#' bms <- get_uniq_bm(yy, blcksz, temp_cvrt = temp_cvrt )
### TO DO : treat the case when not all stations have the same record lengths.
### fill with NA in observational data
get_uniq_bm <- function(data, blcksz, temp_cvrt = NULL, looplastblock = FALSE){

  if(is.vector(data)) {
    data <- data.frame(Obs = data, Station = "X1")
    # if temporal covariate is used
    if(!is.null(temp_cvrt)) {
      nsl <- nrow(data) - blcksz +1
      if(!length(temp_cvrt) == nsl) { temp_cvrt <- temp_cvrt[1:nsl] }

      (data %>% dplyr::group_by(Station) %>%
        tidyr::nest() %>%
        dplyr::mutate( uniq_data  = purrr::map( .x = data, .f = function(.x){
          bmx <- blockmax(.x$Obs, r = blcksz, "sliding", looplastblock = looplastblock)
          bmx <- data.frame(slbm = bmx, temp_cvrt = temp_cvrt)
          bmx %>%  dplyr::group_by(slbm, temp_cvrt) %>%
            dplyr::summarise(n  = dplyr::n(), .groups = "drop") } )) %>%
        dplyr::ungroup() %>%
        dplyr::select(-data))$uniq_data[[1]]
    } else {
      (data %>% dplyr::group_by(Station) %>%
        tidyr::nest() %>%
        dplyr::mutate( uniq_data  = purrr::map( .x =data, .f = function(.x){
          bmx <- blockmax(.x$Obs, r = blcksz, "sliding", looplastblock = looplastblock)
          bmx <- data.frame(slbm = bmx)
          bmx %>%  dplyr::group_by(slbm) %>%
            dplyr::summarise(n  = dplyr::n()) } )) %>%
        dplyr::ungroup() %>%
        dplyr::select(-data))$uniq_data[[1]]
    }
    }

  # if temporal covariate is used
  else {
    if(!is.null(temp_cvrt)) {

      data <- data %>% dplyr::group_by(Station) %>%
        tidyr::nest()
      nsl <- nrow(data[1, ]$data[[1]]) - blcksz +1
      if(!length(temp_cvrt) == nsl) { temp_cvrt <- temp_cvrt[1:nsl] }

      data  %>%
        dplyr::mutate( uniq_data  = purrr::map( .x = data, .f = function(.x){
          bmx <- blockmax(.x$Obs, r = blcksz, "sliding", looplastblock = looplastblock)
          bmx <- data.frame(slbm = bmx, temp_cvrt = temp_cvrt)
          bmx %>%  dplyr::group_by(slbm, temp_cvrt) %>%
            dplyr::summarise(n  = dplyr::n(), .groups = "drop") } )) %>%
        dplyr::ungroup() %>%
        dplyr::select(-data)
    } else {
      data %>% dplyr::group_by(Station) %>%
        tidyr::nest() %>%
        dplyr::mutate( uniq_data  = purrr::map( .x =data, .f = function(.x){
          bmx <- blockmax(.x$Obs, r = blcksz, "sliding", looplastblock = looplastblock)
          bmx <- data.frame(slbm = bmx)
          bmx %>%  dplyr::group_by(slbm) %>%
            dplyr::summarise(n  = dplyr::n()) } )) %>%
        dplyr::ungroup() %>%
        dplyr::select(-data)
    }
  }
}



#' Prepare sliding block maxima data for ML fitting
#'
#' This is a helper function that prepares data of sliding block maxima as well
#' as covariates and response surfaces that are considered such that optimisation
#' gets more efficient.
#'
#' @param loc.sp.form R formula definining the spatial model for the location parameter.
#' @param scale.sp.form R formula definining the spatial model for the scale parameter.
#' @param loc.temp.form R formula definining the temporal trend for the location parameter.
#' @param scale.temp.form R formula definining the temporal trend for the scale parameter.
#' @param data A nested tibble as obtained from \code{\link[slbm]{get_uniq_bm}}.
#' @param spat.cov A data frame containing the spatial covariates used in the formulations
#' of the spatial formulas.
#'
#' @return A list containing the following components:
#' \describe{
#'  \item{data}{The original data that was put in. If there is a trend in either
#'   location and/or scale parameter (i.e. when the corrsesponding formulas are not NULL),
#'   data contains an additional column `MatLocTemp` or `MatScaleTemp` (or both)
#'   containing the design matrices for the temporal
#'   trends as generated by `model.matrix()` (these are vectors of the same length
#'    as the matrices in `uniq_data`).}
#'   \item{MatLocSp}{The design matrix for the spatial location model.}
#'   \item{MatScaleSp}{The design matrix for the spatial scale model.}
#' }
#' @export
#'
#' @examples
#' ##### generate some data #####
#' blcksz <- 90
#' xx <- sapply(c(0.1, 0.4), function(x) evd::rgpd(10*blcksz, shape = x))
#'
#' yy <- data.frame(xx) %>% tidyr::pivot_longer( 1:2, names_to = "Station", values_to = "Obs")
#'
#' # define a temporal covariate that is constant over a block of length blcksz
#' temp_cvrt <- rep(1:10/10, each = blcksz)[1:(9*blcksz + 1)]
#'
#' bms <- get_uniq_bm(yy, blcksz, temp_cvrt = temp_cvrt, looplastblock = FALSE)
#'###############################
#'
#' prep4spatml(loc.sp.form = ~ lon, scale.sp.form = ~ lon +lat,
#' loc.temp.form = ~ GMST, data = bms, spat.cov = data.frame( lon = c(1,2), lat = c(2,4)))
#'
prep4spatml <- function(loc.sp.form, scale.sp.form,
                        loc.temp.form = NULL, scale.temp.form = NULL,
                        data, spat.cov) {

  model.loc.sp <- model.matrix(loc.sp.form, model.frame(loc.sp.form,
                                                        data = spat.cov, na.action = na.pass))

  model.scale.sp <- model.matrix(scale.sp.form, model.frame(scale.sp.form,
                                                            data = spat.cov, na.action = na.pass))


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



  list( data = data, MatLocSp = model.loc.sp, MatScaleSp = model.scale.sp)
}

