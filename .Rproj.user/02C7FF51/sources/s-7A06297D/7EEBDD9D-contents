#' Compute disjoint and sliding block maxima sample of ID process
#'
#' @param data a data frame or tibble of precipitation data, with a column named
#' * datetime: containing the date-time of the observation as POSIXct, POSIXlt object
#' * prec: containing the observed precipitation amount
#' @param ds the durations for which you want to compute the ID process and its annual
#' maxima
#'
#' @param resolution the resolution of the observed data
#' @param seasonlength The length of the season of the input data. E.g., if observations
#' are from the whole year, it is 365; if observations are from JJA only, it is 92.
#'
#' @return a tibble with columns duration, blcksz (the blocksize for the sliding block maxima sample,
#' depending on the duration), djbm (list of tibbles of disjont block maxima),
#' slbm ( list of tibbles og sliding block maxima) (one tibble for each duration)
#' @export
#'
#' @examples
#' dates <- seq(as.POSIXct("2000-01-01 00:00:00"),
#' as.POSIXct("2001-12-31 23:00:00"),by = 'hour')
#' prec <- rgamma(length(dates), shape = 0.1)
#' ExampleData <- data.frame(datetime = dates, prec = prec)
#'
#' get_agg_bm(ExampleData, ds = c(1,2,4,8,16, 24, 48))

get_agg_bm <- function(data, ds, resolution = "hourly", seasonlength = 365){

  agg_df <-  fun_aggregate2df(data, ds, resolution = resolution,
                              seasonlength = seasonlength) # dataframe of ID-process

  n_daily <- ifelse( resolution == "daily", 1, 24)  # number of daily observations

  r <- seasonlength*n_daily    # blocksize

  agg_df <- agg_df %>%
    dplyr::mutate(slbm = purrr::map(.x = aggs, function(.x){
      bmx <- blockmax(c(.x$agg.sum , .x$agg.sum[1: (r-1)]), r = r, "sliding")
      years <- .x$Year[1:length(bmx)]
      return(data.frame(sldata = bmx, Year = years))
      }) ) %>%
    dplyr::select(- aggs)


  return(agg_df)
}
