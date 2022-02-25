
#' Compute concatenated block maxima
#' @description Compute sliding block maxima for all possible compositions of years
#' arising from cross validation with test set of size given in 'nlo'
#'
#' @param dailydata dataframe with hourly/daily ... observations in column named 'obs',
#' and an index for each observation (e.g. the year/season in which the observation was
#' made) in column named 'Index'
#' @param nlo number of years that make up the test set
#' @param blcksz blocksize for computing sliding block maxima
#'
#' @return tibble with column 'ind1' giving the first index of observations making up
#' the concatenated time series; column 'ind2' giving the second index of observations
#' making up the concatenated time series; column 'newslbm': a list
#' of vectors containing the sliding block maxima of the concatenated time series
#'
#' Note that 'ind1' is allowed to be larger than 'ind2'
#' @export
#'
#' @examples
#' df <- data.frame(obs = evd::rgpd(30*50), Index = rep(1:50, each = 30))
#' compute_conc_bm(df, nlo = 3, blcksz = 30)
compute_conc_bm <- function(dailydata, nlo = 3, blcksz){


  years_obs <- 1:ceiling(nrow(dailydata)/blcksz)
  ny <- length(years_obs)

  start_year <- tibble::tibble(ind1 = rep(years_obs, each = nlo), index = rep(2:(nlo+1), ny))
  start_year <- start_year %>% dplyr::mutate(ind2 = ind1 + index)
  start_year <- start_year %>% dplyr::mutate(ind2 = ifelse(ind2 <= max(years_obs), ind2, ind2 - ny))
  start_year <- start_year %>% dplyr::select(-index)

  df_conc_slbms <-  start_year %>%
    dplyr::mutate( newslbm = purrr::map2(.x = ind1, .y = ind2 ,
                                         .f = function(.x, .y){
                                  bm_concat_years(dailydata = dailydata, ind1 = .x, ind2 = .y,
                                                                    blcksz = blcksz)
                                         }) )

  return(df_conc_slbms)
}
