#' d-GEV parameter estimation and quantile score evaluation on a test set
#'
#' @param agg_bm Tibble as obtained from the output of [get_agg_bm()()].
#' @param conc_bm Tibble of sliding block maxima for all possible combinations of
#' years/seasons that are consecutive after removing observations from a year/season
#' for the test set
#' @param leave_out The index (i.e. years/seasons ...) that make up the test set.
#' These years are left out in the fitting procedure.
#' @param ds Vector of durations for which quantile of the d-GEV is estimated.
#' @param quants Vector of quantiles to estimate based on the d-GEV fit.
#' @param mult_sc logical. Whether to fit a multi-scale parameter eta2
#' @param dur_offset logical. Whether to fit a duration offset parameter theta
#' @param int_offset logical. Whether to fit an intensity offset parameter tau
#' @param optimMethod Optimisation method used in [optim()].
#' @param Maxit Passed to the control argument of [optim()]:
#' maximum number of iterations performed during the optimisation.
#' @param ... Further arguments that can be passed to the control argument of [optim()].
#'
#' @return Returns a tibble containing the estimated parameters, the estimation method,
#' a code for the convergence (zero indicates successful convergence), the values of \eqn{p} for
#' which quantiles are estimated, the values of the durations \eqn{d} for which
#' quantiles are estimated, the test set (data that is used for computing the quantile score)
#' and the quantile score. Quantile Scores can only be computed for those durations which are also
#' present in the testdata.
#'
#' @export
#'
#' @examples dates <- seq(as.POSIXct("2000-01-01 00:00:00"),
#' as.POSIXct("2020-12-31 23:00:00"),by = 'hour')
#' prec <- rgamma(length(dates), shape = 0.1)
#' example_data <- data.frame(datetime = dates, prec = prec)
#'
#' agbm <- get_agg_bm(example_data, ds = c(1,2,4,8,16, 24, 48))
#' agdf <- fun_aggregate2df( example_data, ds = c(1,2,4,8,16, 24, 48) )
#' concbm <- compute_conc_bm_id(agdf)
#' dgev_cv_fit(agbm, concbm, ds = c(1,2) , leave_out = c(2004, 2008, 2017),
#'  quants = c(.5, .9, .99) )
#'
dgev_cv_fit <- function(agg_bm, conc_bm, leave_out , ds, quants, mult_sc = FALSE,  dur_offset = FALSE,
                       int_offset = FALSE, optimMethod = "Nelder-Mead", Maxit = 1500, ...){

  agg_bm_test <- agg_bm %>% dplyr::select(djbm) %>%
    tidyr::unnest(cols = djbm) %>% dplyr::filter(Year %in% leave_out)

  dj_train <- agg_bm %>% dplyr::select(-slbm) %>%
    dplyr::mutate(djbm = purrr::map(djbm, ~ .x %>% dplyr::filter(!(Year %in% leave_out))))

  minmaxn <- do.call( function(.x){data.frame(Min = min(.x), Max = max(.x),
                                              Nobs = length(.x))} , list(agg_bm[1,]$djbm[[1]]$Year))


  sl_lo <- data.frame(year1 = leave_out, year2 = leave_out -1)
  sl_lo <- sl_lo %>% dplyr::mutate(year2 = ifelse(year2 >= minmaxn$Min, year2, year2 + minmaxn$Nobs))
  sl_train1 <- agg_bm %>% dplyr::select( - djbm) %>%
    dplyr::mutate(slbm = purrr::map(slbm, ~ .x %>% dplyr::filter(!(Year %in% unlist(sl_lo)))))

  conc_bm <- conc_bm %>% dplyr::filter(ind1 %in% sl_lo$year2)  ## filter for possible new starting years

  # neither the new starting year nor the following year can be in the new sliding block
  conc_bm <- conc_bm %>%
    dplyr::filter(!(ind1 %in% leave_out) ) %>%
    dplyr::filter(!(ind2 %in% leave_out))

  conc_bm <- conc_bm %>% dplyr::mutate(Diff = ind2 - ind1) %>%
    dplyr::mutate(Diff = ifelse(Diff >= 0, Diff, Diff + minmaxn$Nobs)) %>%
    dplyr::group_by(ind1) %>%
    dplyr::filter(Diff == min(Diff)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-Diff)

  conc_bm <- conc_bm %>% tidyr::unnest(cols = newslbm) %>%
    dplyr::select(-c( ind2)) %>%
    tidyr::unnest(cols = conc_slbm) %>%
    dplyr::group_by(duration) %>%
    tidyr::nest(cols = c( ind1 , conc_slbm))  %>% dplyr::ungroup()


  newsample <- conc_bm %>% dplyr::left_join(sl_train1, by = "duration") %>%
    dplyr::mutate(slbm = purrr::map2(.x = cols, .y = slbm , function(.x, .y){
      .x <- dplyr::rename(.x, "Year" = ind1, "sldata" = conc_slbm)
      dplyr::bind_rows(.x, .y)

    })) %>% dplyr::select(-cols)


  agg_bm <- dplyr::right_join(dj_train, newsample, by = "duration")



  fitb <- gev.d.fit.sl(agg_bm = agg_bm, mult_sc = mult_sc,
                                 dur_offset = dur_offset, int_offset = int_offset, method = optimMethod,
                       Maxit = Maxit, ...)

  mlest <- dplyr::bind_cols(fitb$mle, conv =  fitb$conv)
  rmcol <- character()
  if(!mult_sc) {
    mlest <- mlest %>% dplyr::bind_cols( eta2 = 0)
    rmcol <- c( rmcol, "eta2")
  }
  if(!int_offset) {
    mlest <- mlest %>% dplyr::bind_cols( tau = 0)
    rmcol <- c(rmcol, "tau")
  }
  if(!dur_offset) {
    mlest <- mlest %>% dplyr::bind_cols( theta = 0)
    rmcol <- c(rmcol, "theta")
  }
  table_quants <- mlest %>% tidyr::expand_grid(ps = quants,
                                               ds = ds)


  table_quants <- table_quants %>%
    dplyr::mutate( est.quants = purrr::pmap_dbl( list(mut, sigma0, shape, theta,  eta, eta2, tau, ps, ds) ,
                                                 .f = function( .mut, .sigma0, .shape, .theta,
                                                                .eta, .eta2, .tau, .ps, .d){
                                                   IDF::qgev.d(.ps, mut = .mut, sigma0 = .sigma0, xi = .shape ,
                                                               eta = .eta , eta2 = (.eta2 + .eta), tau = .tau ,
                                                               theta = .theta, d = .d)

                                                 }
    ))

  table_quants <- table_quants[ , !(colnames(table_quants) %in% rmcol)]

  zns  <- agg_bm_test %>%
    dplyr::select(-Year)

  zns <- dplyr::rename(zns, "ds"= "duration")
  zns <- zns %>% dplyr::group_by(ds) %>% tidyr::nest(Testdata = djdata)
  table_quants <- dplyr::left_join(table_quants, zns, by = "ds")
  table_quants <- table_quants %>%
    dplyr::mutate( quant_score = purrr::pmap_dbl( .l = list( .ps = ps , .q = est.quants, .obs = Testdata),
                                                  .f = function(.ps, .q, .obs){
                                                    if(is.null(.obs)){ NA} else{
                                                    compute_qs(obs = unlist(.obs) , quant.est = .q,
                                                               p = .ps) }
                                                  })) %>%
    dplyr::arrange(ps, ds)

  return(table_quants)
}
