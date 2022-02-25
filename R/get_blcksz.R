#' Compute blocksize for duration d
#'
#' @param d the duration
#' @param resolution resolution of measurements, must be one of "hourly" or "daily"
#' @return returns the number of observations from an intensity duration process
#' that fall into one year, when one has hourly observations on 365 days
#' @export
#'
#' @examples get.blcksz(4)
get_blcksz <- function(d, resolution = "hourly"){
  switch( resolution, "hourly"  = 365*24- d +1,
                      "daily" = 365 - d +1 )
}
