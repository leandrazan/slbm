
#' Get dataframes for use in quantile score functions
#'
#' @param data Dataframe or tibble with variables named 'obs' and 'Index' containing
#' the observations and the year/season the observation was made in.
#' @param blcksz The blocksize to compute blockmaxima
#' @param return_sldata logical; whether to return dataframe containing sliding data
#' @param return_djdata logical; whether to return dataframe containing disjoint data
#'
#' @return list with elements djbm and slbm
#' @export
#'
#' @examples
#'
#' get_bm_data_qs(data.frame( obs = rnorm(50*30), Index = rep(1:50, each = 30)),
#' blcksz = 30)
get_bm_data_qs <- function(data, blcksz, return_sldata = TRUE, return_djdata = TRUE){

  stopifnot("The dataframe must contain variables named 'obs' and 'Index' " =
              c("Index", "obs") %in% colnames(data) )

  dj.bm <- blockmax(data$obs, r = blcksz, "disjoint")
  n.years <- length(dj.bm)
  dj.bm <- data.frame( obs = dj.bm, Index = 1:n.years)

  conc_bm <- compute_conc_bm(data, blcksz = blcksz)
  sl.bm <- tibble::tibble( sldata = blockmax(data$obs, r = blcksz, "sliding"))
  sl.bm$Index <- rep(1:n.years, each = blcksz)[1:nrow(sl.bm)]
  if( return_sldata & return_djdata){
   return(list(djbm = dj.bm, slbm = sl.bm))
  } else if (return_sldata & !(return_djdata) ) {
    return(list(slbm = sl.bm))
  }  else (
    return(list(djbm = dj.bm))
  )
}
