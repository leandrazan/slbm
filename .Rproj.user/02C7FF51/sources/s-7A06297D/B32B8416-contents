
bm_concat_years_id <- function(agg_df, year1, year2, blcksz = 8760){
  agg_df <- agg_df %>% dplyr::mutate(conc_ts = purrr::map(aggs, function(.x){
    .x %>% dplyr::filter(Year %in% c(year1, year2))
  }))
  agg_df %>%
    dplyr::mutate(conc_slbm = purrr::map(conc_ts, ~ blockmax(.x$agg.sum, r = blcksz, "sliding"))) %>%
    dplyr::select(conc_slbm, duration)
}


#' Compute concatenated block maxima
#' @description Compute sliding block maxima of the intensity duration process
#' for all possible compositions of two consecutive years after removing observations
#' from a testset.
#'
#' For performing cross-validation, one splits the data into several test and a training
#' data sets. In the sliding block case, one cannot only remove the observations with
#' an index belonging to the test set, because the sliding blocks 'look into the future',
#' which means that a sliding block starting in year 2 also looks at observations from
#' year 3 (if it's not the block corresponding to the disjoint block).
#' Because it takes quite long to remove the data of the years that are chosen as
#' test set from the original data and then computing sliding block maxima, this
#' functions computes the sliding block maxima of only those years that can be consecutive
#' after removing observations belonging to a testset. E.g., when the testset consists of
#' years 2002, 2003, 2004 and we have data from 2000 to 2010, we only need to compute
#' sliding block maxima of the concatenated series with data from 2002 and 2005.
#'
#'
#' @param agg_df Dataframe with values of intensity duration process; output of
#' function 'fun_aggregate2df'
#' @param nlo number of Indices (years/seasons ...) that make up the test set
#' @param resolution Resolution of the data
#' @param seasonlength The length of the season of the input data. E.g., if observations
#' are from the whole year, it is 365; if observations are from JJA only, it is 92.
#' @param testset Whether to compute only those concatenated block maxima that are
#' needed when cross-validation sets are obtained by dividing the data into disjoint
#' subsets of consecutive years. Defaults to "consec" (which is this reduced setting),
#' anything else computes block maxima for all possibilities of first and second year.
#' @return A tibble with
#' * ind1: The year from which
#' * newslbm: sliding block maxima as obtained when observations with the i
#' @export
#'
#' @examples
#' dates <- seq(as.POSIXct("2000-01-01 00:00:00"),
#' as.POSIXct("2005-12-31 23:00:00"),by = 'hour')
#' prec <- rgamma(length(dates), shape = 0.1)
#' ExampleData <- data.frame(datetime = dates, prec = prec)
#'
#' aggdat <-  fun_aggregate2df(ExampleData, ds = c(1,5,8))
#' compute_conc_bm_id(aggdat, nlo = 2, testset = "consec")
#'
compute_conc_bm_id <- function(agg_df, nlo = 3, resolution = "hourly", seasonlength = 365,
                               testset = "consec"){

  n_daily <- ifelse( resolution == "daily", 1, 24)  # number of daily observations

  blcksz <- seasonlength*n_daily    # blocksize

  years_obs <- agg_df$djbm[[1]]$Year
  ny <- length(years_obs)

  start_year <- tibble::tibble(ind1 = rep(years_obs, each = nlo), index = rep(2:(nlo+1), ny))
  start_year <- start_year %>% dplyr::mutate(ind2 = ind1 + index)
  start_year <- start_year %>%
    dplyr::mutate(ind2= ifelse(ind2 <= max(years_obs), ind2, ind2 - ny))
  start_year <- start_year %>% dplyr::select(-index)

  if(testset == "consec"){
    styseq <- seq(years_obs[(nlo +1)], years_obs[ny], nlo)  # sequence with first
    # years of testsets
    start_year <- start_year %>% dplyr::mutate(Diff = abs(ind1 - ind2) ) %>%
      dplyr::filter(Diff == (nlo+1), ind1 %in% (styseq -1)) %>%
      dplyr::select(-Diff)
  }

  df_conc_slbms <-  start_year %>%
    dplyr::mutate( newslbm = purrr::map2(.x = ind1, .y = ind2  ,
                           .f = function(.x, .y){
                             tryCatch(bm_concat_years_id(agg_df = agg_df, year1 = .x, year2 = .y, blcksz = blcksz),
                                      error = function(e) NA)
                           }) )

  return(df_conc_slbms)
}

