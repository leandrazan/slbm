#' Compute sliding block maxima of the time series obtained by concatenating
#'  observations from two arbitrary years.
#'
#' @param dailydata dataframe with hourly/daily ... observations in column named 'obs',
#' and the Index for the observation (e.g. the year/season in which the observation was made)
#'in column named 'Index'
#' @param ind1 Observations with that index make up the first part of the
#' concatenated series
#' @param ind2 Observations with that index make up the first part of the
#' concatenated series
#' @param blcksz The blocksize for computing sliding block maxima.
#'
#' @return a numeric vector containing the sliding block maxima
#'
#' @export
#'
#' @examples
#' df <- data.frame(obs = evd::rgpd(30*5), Index = rep(1:5, each = 30))
#' bm_concat_years(df, 2, 4, blcksz = 30)
bm_concat_years <- function(dailydata, ind1, ind2, blcksz = 90){
  conc_ts <- dailydata %>% dplyr::filter(Index %in% c(ind1, ind2))

  blockmax(conc_ts$obs, r = blcksz, "sliding")
}
