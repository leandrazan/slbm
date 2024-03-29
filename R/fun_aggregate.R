#' Compute Intensity duration process and yearly block maxima for one duration \eqn{d}
#'
#' @param data a data frame or tibble of precipitation data, with a column named
#' * datetime: containing the date-time of the observation as POSIXct, POSIXlt object
#' * prec: containing the observed precipitation amount
#'
#' @param d the duration for which you want to compute the intensity process.
#' It is the number \eqn{d} of observations over which the the moving average
#' of precipitation is computed, computed within distinct years.
#'
#' @param resolution The resolution of the data, must be one of "hourly" or "daily"
#'
#' @param seasonlength The length of the season of the input data. E.g., if observations
#' are from the whole year, it is 365; if observations are from JJA only, it is 92.
#' @return returns a tibble with colum
#' * aggs: tibble of aggregated observations, Year they belong to and the duration d
#' * djbm: tibble of annual maxima of the intensity duration process, along with
#'         corresponding year and duration
#' @export
#'
#' @examples dates <- seq(as.POSIXct("2000-01-01 00:00:00"),
#' as.POSIXct("2001-12-31 23:00:00"),by = 'hour')
#' prec <- rgamma(length(dates), shape =0.1)
#' ExampleData <- data.frame(datetime = dates, prec = prec)
#'
#' fun_aggregate(ExampleData, 120)
#'
fun_aggregate <- function( data, d, resolution = "hourly", seasonlength  = 365){

  d1 <- d
  if(resolution == "daily"){ d1 <- d/24}

  n_daily <- ifelse( resolution == "daily", 1, 24)  # number of daily observations

  r <- seasonlength*n_daily    # blocksize

  data <- data %>% dplyr::mutate(Year = lubridate::year(datetime)) %>%
    dplyr::select(-c(datetime))
  ny <- length(unique(data$Year))
  agg.sum <-  RcppRoll::roll_mean(data$prec, n = d1, align = "right", na.rm = TRUE)

  aggs <- data[ d1 :nrow(data), ] %>% dplyr::select(Year) %>%
    dplyr::bind_cols(agg.sum = agg.sum) %>%
    dplyr::mutate( duration = d)

  djbm <- blockmax(aggs$agg.sum, r = r, "disjoint")[1: ny]
  DJBM <-  data.frame( Year = unique(data$Year), djdata = djbm, duration = d)

  return(tibble::tibble(aggs = list(aggs), djbm = list(DJBM)))
}
