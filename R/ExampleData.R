#' Synthetic example dataset of `daily` observations at 8 stations.
#'
#'
#' @format A matrix with 3600 rows and 8 columns.
#' Each column contains measurements of a station.
#'
#' @details
#' The example data was generated using the \code{\link[slbm]{sim_spatial}}
#' function as follows:
#'
#' @examples
#' \dontrun{
#' ### dataset was generated as follows:
#' data("GMST")
#' set.seed(1)
#' ny <- 40
#' blcksz <- 90
#' coords <- cbind(spatial_cvrt$lat, spatial_cvrt$lon)
#' colnames(coords ) <- c("lat", "lon")
#' mdax <- MDAdata(ny*blcksz, locations = coords, margins = list(distr = "gpd", shape = -0.2),
#'                 ms = list(kovmod = "whitmat", nugget = 0, range = 3, smooth = 0.44))
#'
#' tempcov1 <- data.frame(GMST = rep(GMST$smoothedGMST[101:140], each = 90))
#' xx <- sim_spatial(mdax, params = c("loc0" = 10, "loc1" = 2, "loc2" = 0.05,
#'                                    "scale0" = 2, "scale1" = 0.2, "shape" = -0.2,
#'                                    "tempLoc" = 2),
#'                   loc.sp.form = ~ lon + lat, scale.sp.form = ~ lon,
#'                 spat.cov = spatial_cvrt, temp.cov = tempcov1,
#'                   loc.temp.form = ~ GMST)
#'
#' plot.ts(apply(xx, 2, blockmax, r = 90, "disjoint"))
#'}
#'
"ExampleData"


