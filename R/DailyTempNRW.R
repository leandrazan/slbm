#' Daily Temperature Series at 9 stations in NRW.
#'
#'
#' @format A list with components
#' \describe{
#' \item{TimeSeries}{ tibble (8,188 x 10): the first column `Datetime` contains
#' the date of the observation, the columns 2 to 9 the observed daily maximum temperature
#' at the respective station (see column names) as made during summer months (June, July, August)
#' }
#' \item{SpatialCov}{Spatial covariates such as Station ID, Station name, elevation,
#' latitude, longitude.}
#' }
#' @details
#' The data was downloaded from the DWD.
#' @source \url{https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/daily/kl/historical/}
#'
"DailyTempNRW"

