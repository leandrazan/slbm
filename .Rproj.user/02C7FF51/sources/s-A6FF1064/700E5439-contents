#' Compute intensity duration process and block maxima for several durations $d$
#'
#' @param data a data frame or tibble of precipitation data, with a column named
#' * datetime: containing the date-time of the observation as POSIXct, POSIXlt object
#' * prec: containing the observed precipitation amount
#' @param ds the durations for which you want to compute the intensity process.
#' For each $d$, the intensity process is computed, which is the moving average over d observations
#' of precipitation, computed within distinct years.
#'
#' @param resolution The resolution of the data, must be one of "hourly" or "daily"
#'
#' @return returns a tibble with colums
#' * aggs: list of tibbles of aggregated observations, Year they belong to and the duration d.
#' One tibble for each duration $d$
#' * djbm: list of tibbles of annual maxima of the intensity duration process, along with
#'         corresponding year and duration. One tibble for each duration $d$.
#' @export
#'
#' @examples dates <- seq(as.POSIXct("2000-01-01 00:00:00"),
#' as.POSIXct("2001-12-31 23:00:00"),by = 'hour')
#' prec <- rgamma(length(dates), shape = 0.1)
#' ExampleData <- data.frame(datetime = dates, prec = prec)
#'
#' fun_aggregate2df(ExampleData, ds = c(1,2,4,8,16, 24, 48))
fun_aggregate2df <- function(data, ds, resolution = "hourly"){
  purrr::map_dfr(ds, .f = function(ds){ return( fun_aggregate(data = data, d = ds,
                                                              resolution = resolution))}) %>%
    dplyr::bind_cols(duration = ds)

}
