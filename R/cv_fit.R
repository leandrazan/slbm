
#' cross-validated d-GEV fit based on disjoint and sliding BM
#'
#' @param agg_bm Output of get_agg_bm, i.e. a tibble containing the sample of
#'  disjoint and sliding block maxima of the ID process, for different durations $d$
#' @param leave_out the years to use as test data; observations of these years are
#' excluded in training data used for the fit, fit is evaluated on these years
#' @param ds durations for which the quantile skill/ score/ index should be computed
#' @param quants quantiles for which the quantile skill/ score/ index should be computed
#' @param dur_offset whether the duration offset parameter should be fitted
#'
#' @return a tibble with estimated parameters and quantile scores, computed on
#' disjoint and sliding block maxima
#' @export
#'
#' @examples  dates <- seq(as.POSIXct("2000-01-01 00:00:00"),
#' as.POSIXct("2020-12-31 23:00:00"),by = 'hour')
#' prec <- rgamma(length(dates))
#' example_data <- data.frame(MESS_DATUM = dates, RR = prec)
#'
#' agbm <- get_agg_bm(example_data, ds = c(1,2,4,8,16, 24, 48))
#' cv_fit(agbm, leave_out = 2001:2003, ds =   c(1,2,4,8,16, 24, 48) , quants = c(.5, .9, .99) )
cv_fit <- function(agg_bm, leave_out , ds, quants, mult_sc = FALSE, dur_offset = F){


  fitb <- gev.d.fit.sl.ms(agg_bm = agg_bm, leave_out = leave_out , mult_sc = mult_sc,
                          dur_offset = dur_offset)
  mlest <- fitb$mle
  if(!mult_sc) {
    mlest <- mlest %>% bind_cols( eta2 = 0)
  }
  table_quants <- mlest %>% tidyr::expand_grid(ps = quants,
                                               ds = ds)


  table_quants <- table_quants %>%
    dplyr::mutate( est.quants = purrr::pmap_dbl( list(mut, sigma0, shape,  eta, eta2,  ps, ds) ,
                                                 .f = function( .mut, .sigma0, .shape,  .eta, .eta2,  .ps, .d){
                                                   IDF::qgev.d(.ps, mut = .mut, sigma0 = .sigma0, xi = .shape ,
                                                               eta = .eta , eta2 = (.eta2 + .eta), theta = 0,  d = .d)

                                                 }
    ))

  zns  <- agg_bm %>% dplyr::select( djbm) %>%
    tidyr::unnest(cols = djbm) %>%
    dplyr::filter(Year %in% leave_out) %>%
    select(-Year)

  zns <- dplyr::rename(zns, "ds"= "duration")
  zns <- zns %>% dplyr::group_by(ds) %>% tidyr::nest(djdat = djdata)
  table_quants <- dplyr::left_join(table_quants, zns, by = "ds")
  table_quants <- table_quants %>%
    dplyr::mutate( quant_score = purrr::pmap_dbl( .l = list( .ps = ps , .q = est.quants, .obs = djdat),
                                                  .f = function(.ps, .q, .obs){
                                                    compute_qs(obs = unlist(.obs) , quant.est = .q,
                                                               p = .ps)
                                                  })) %>%
    dplyr::arrange(ps, ds)

  return(table_quants)
}



cv_fit_slnew  <- function(agg_df, leave_out , ds, quants, mult_sc = FALSE,  dur_offset = F){

  agg_bm_test <- agg_df %>% select(djbm) %>% unnest(cols = djbm) %>% filter(Year %in% leave_out)

  agg_df <- agg_df %>% mutate(aggs = map(aggs, ~ .x %>% filter(!(Year %in% leave_out))),
                              djbm = map(djbm, ~ .x %>% filter(!(Year %in% leave_out))))

  agg_df <-  agg_df %>%
    dplyr::mutate(duration = ds) %>%
    dplyr::mutate(blcksz = purrr::map_dbl(duration, get.blcksz))


  agg_bm <- agg_df %>%
    dplyr::mutate( slbm = purrr::map2( .x  = aggs, .y = blcksz,
                                       .f = function(.x, .y){
                                         sldf <- data.frame(sldata = blockmax(.x$agg.sum, r = .y, "sliding"))
                                         return(dplyr::bind_cols(sldf, Year = .x$Year[1:nrow(sldf)] ))
                                       }) ) %>%
    select( - aggs)

  fitb <- gev.d.fit.sl.slbmnew(agg_bm = agg_bm, leave_out = 0 , mult_sc = mult_sc,
                          dur_offset = dur_offset)
  mlest <- fitb$mle
  if(!mult_sc) {
    mlest <- mlest %>% bind_cols( eta2 = 0)
  }
  table_quants <- mlest %>% tidyr::expand_grid(ps = quants,
                                               ds = ds)


  table_quants <- table_quants %>%
    dplyr::mutate( est.quants = purrr::pmap_dbl( list(mut, sigma0, shape,  eta, eta2,  ps, ds) ,
                                                 .f = function( .mut, .sigma0, .shape,  .eta, .eta2,  .ps, .d){
                                                   IDF::qgev.d(.ps, mut = .mut, sigma0 = .sigma0, xi = .shape ,
                                                               eta = .eta , eta2 = (.eta2 + .eta), theta = 0,  d = .d)

                                                 }
    ))

  zns  <- agg_bm_test %>%
    select(-Year)

  zns <- dplyr::rename(zns, "ds"= "duration")
  zns <- zns %>% dplyr::group_by(ds) %>% tidyr::nest(djdat = djdata)
  table_quants <- dplyr::left_join(table_quants, zns, by = "ds")
  table_quants <- table_quants %>%
    dplyr::mutate( quant_score = purrr::pmap_dbl( .l = list( .ps = ps , .q = est.quants, .obs = djdat),
                                                  .f = function(.ps, .q, .obs){
                                                    compute_qs(obs = unlist(.obs) , quant.est = .q,
                                                               p = .ps)
                                                  })) %>%
    dplyr::arrange(ps, ds)

  return(table_quants)
}


cv_fit_cont <- function(agg_df, leave_out , ds, quants, mult_sc = FALSE,  dur_offset = F){

  agg_bm_test <- agg_df %>% select(djbm) %>% unnest(cols = djbm) %>% filter(Year %in% leave_out)

  agg_df <- agg_df %>% dplyr::mutate(aggs = map(aggs, ~ .x %>% filter(!(Year %in% leave_out))),
                              djbm = map(djbm, ~ .x %>% filter(!(Year %in% leave_out)))) %>%
    dplyr::mutate(duration = ds)


  agg_bm <- agg_df %>% mutate(slbm = map(.x = aggs, function(.x){
    bmx <- blockmax(.x$agg.sum , r = 8760, "sliding")
    # nbmx <- length(bmx)
    years <- .x$Year[1:length(bmx)]
    return(data.frame(sldata = bmx, Year = years, duration = .x$duration[1]))
  }) ) %>%
    select(- aggs)



  fitb <- gev.d.fit.sl.slbmnew(agg_bm = agg_bm, leave_out = 0 , mult_sc = mult_sc,
                               dur_offset = dur_offset)
  mlest <- fitb$mle
  if(!mult_sc) {
    mlest <- mlest %>% bind_cols( eta2 = 0)
  }
  table_quants <- mlest %>% tidyr::expand_grid(ps = quants,
                                               ds = ds)


  table_quants <- table_quants %>%
    dplyr::mutate( est.quants = purrr::pmap_dbl( list(mut, sigma0, shape,  eta, eta2,  ps, ds) ,
                                                 .f = function( .mut, .sigma0, .shape,  .eta, .eta2,  .ps, .d){
                                                   IDF::qgev.d(.ps, mut = .mut, sigma0 = .sigma0, xi = .shape ,
                                                               eta = .eta , eta2 = (.eta2 + .eta), theta = 0,  d = .d)

                                                 }
    ))

  zns  <- agg_bm_test %>%
    select(-Year)

  zns <- dplyr::rename(zns, "ds"= "duration")
  zns <- zns %>% dplyr::group_by(ds) %>% tidyr::nest(djdat = djdata)
  table_quants <- dplyr::left_join(table_quants, zns, by = "ds")
  table_quants <- table_quants %>%
    dplyr::mutate( quant_score = purrr::pmap_dbl( .l = list( .ps = ps , .q = est.quants, .obs = djdat),
                                                  .f = function(.ps, .q, .obs){
                                                    compute_qs(obs = unlist(.obs) , quant.est = .q,
                                                               p = .ps)
                                                  })) %>%
    dplyr::arrange(ps, ds)

  return(table_quants)
}

cv_fit_eff <- function(agg_bm, conc_bm, leave_out , ds, quants, mult_sc = FALSE,  dur_offset = FALSE,
                       int_offset = FALSE, Maxit = 1500){

  agg_bm_test <- agg_bm %>% select(djbm) %>% unnest(cols = djbm) %>% filter(Year %in% leave_out)

  dj_train <- agg_bm %>% select(-slbm) %>%
    dplyr::mutate(djbm = map(djbm, ~ .x %>% filter(!(Year %in% leave_out))))

  minmaxn <- do.call( function(.x){data.frame(Min = min(.x), Max = max(.x),
                                              Nobs = length(.x))} , list(agg_bm[1,]$djbm[[1]]$Year))


  sl_lo <- data.frame(year1 = leave_out, year2 = leave_out -1)
  sl_lo <- sl_lo %>% mutate(year2 = ifelse(year2 >= minmaxn$Min, year2, year2 + minmaxn$Nobs))


  sl_train1 <- agg_bm %>% select( - djbm) %>%
    dplyr::mutate(slbm = map(slbm, ~ .x %>% filter(!(Year %in% unlist(sl_lo)))))

  conc_bm <- conc_bm %>% filter(sty %in% sl_lo$year2)  ## filter for possible new starting years

  # neither the new starting year nor the following year can be in the new sliding block
  conc_bm <- conc_bm %>%
    filter(!(sty %in% leave_out) ) %>%
    filter(!(year_f %in% leave_out))

  conc_bm <- conc_bm %>% mutate(Diff = year_f - sty) %>%
    mutate(Diff = ifelse(Diff >= 0, Diff, Diff + minmaxn$Nobs)) %>%
    group_by(sty) %>% filter(Diff == min(Diff)) %>% ungroup() %>% select(-Diff)

  conc_bm <- conc_bm %>% unnest(cols = newslbm) %>%
    select(-c( year_f)) %>%
    unnest(cols = conc_slbm) %>%
    group_by(duration) %>%
    nest(cols = c( sty , conc_slbm))  %>% ungroup()


  newsample <- conc_bm %>% left_join(sl_train1, by = "duration") %>%
    mutate(slbm = map2(.x = cols, .y = slbm , function(.x, .y){
    .x <- rename(.x, "Year" = sty, "sldata" = conc_slbm)
    bind_rows(.x, .y)

    })) %>% select(-cols)


  agg_bm <- right_join(dj_train, newsample, by = "duration")



  fitb <- gev.d.fit.sl.msfl.curv(agg_bm = agg_bm, leave_out = 0 , mult_sc = mult_sc,
                               dur_offset = dur_offset, int_offset = int_offset, Maxit = Maxit)

  mlest <- bind_cols(fitb$mle, conv = c(NA, fitb$conv))

  if(!mult_sc) {
    mlest <- mlest %>% bind_cols( eta2 = 0)
  }
  if(!int_offset) {
    mlest <- mlest %>% bind_cols( tau = 0)
  }
  if(!dur_offset) {
    mlest <- mlest %>% bind_cols( theta = 0)
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

  zns  <- agg_bm_test %>%
    select(-Year)

  zns <- dplyr::rename(zns, "ds"= "duration")
  zns <- zns %>% dplyr::group_by(ds) %>% tidyr::nest(djdat = djdata)
  table_quants <- dplyr::left_join(table_quants, zns, by = "ds")
  table_quants <- table_quants %>%
    dplyr::mutate( quant_score = purrr::pmap_dbl( .l = list( .ps = ps , .q = est.quants, .obs = djdat),
                                                  .f = function(.ps, .q, .obs){
                                                    compute_qs(obs = unlist(.obs) , quant.est = .q,
                                                               p = .ps)
                                                  })) %>%
    dplyr::arrange(ps, ds)

  return(table_quants)
}


# agg_df <- fun.aggregate4sl2df(data1, ds = 2^(0:5))
#  agg_bm <- get_agg_bm4sl(agg_df)
#  conc_bm <- compute_conc_bm(agg_df, nlo = 5)
# cv_fit_eff(agg_bm = agg_bm, conc_bm = conc_bm, leave_out = 1999:2003, ds = 2^(0:5), quants = .99,
#            int_offset = TRUE)



cv_fit_sleval <- function(agg_bm, agg_df,  conc_bm, leave_out , ds, quants, mult_sc = FALSE,  dur_offset = FALSE,
                       int_offset = FALSE, Maxit = 1500){

  agg_bm_test <- agg_bm %>% select(djbm) %>% unnest(cols = djbm) %>% filter(Year %in% leave_out)

  dj_train <- agg_bm %>% select(-slbm) %>%
    dplyr::mutate(djbm = map(djbm, ~ .x %>% filter(!(Year %in% leave_out))))

  minmaxn <- do.call( function(.x){data.frame(Min = min(.x), Max = max(.x),
                                              Nobs = length(.x))} , list(agg_bm[1,]$djbm[[1]]$Year))


  sl_lo <- data.frame(year1 = leave_out, year2 = leave_out -1)
  sl_lo <- sl_lo %>% mutate(year2 = ifelse(year2 >= minmaxn$Min, year2, year2 + minmaxn$Nobs))
  sl_train1 <- agg_bm %>% select( - djbm) %>%
    dplyr::mutate(slbm = map(slbm, ~ .x %>% filter(!(Year %in% unlist(sl_lo)))))

  sl_train <- agg_df %>% select(aggs, duration) %>%
    dplyr::mutate(aggs = map(aggs, ~ .x %>% filter(Year %in% leave_out))) %>%
    mutate( BMtest = map(aggs, ~ unique(blockmax(c(.x$agg.sum, .x$agg.sum[1:8759]), r = 8760 , "sliding")))) %>%
    select(BMtest, duration)

  conc_bm <- conc_bm %>% filter(sty %in% sl_lo$year2)  ## filter for possible new starting years

  # neither the new starting year nor the following year can be in the new sliding block
  conc_bm <- conc_bm %>%
    filter(!(sty %in% leave_out) ) %>%
    filter(!(year_f %in% leave_out))

  conc_bm <- conc_bm %>% mutate(Diff = year_f - sty) %>%
    mutate(Diff = ifelse(Diff >= 0, Diff, Diff + minmaxn$Nobs)) %>%
    group_by(sty) %>% filter(Diff == min(Diff)) %>% ungroup() %>% select(-Diff)

  conc_bm <- conc_bm %>% unnest(cols = newslbm) %>%
    select(-c( year_f)) %>%
    unnest(cols = conc_slbm) %>%
    group_by(duration) %>%
    nest(cols = c( sty , conc_slbm))  %>% ungroup()


  newsample <- conc_bm %>% left_join(sl_train1, by = "duration") %>%
    mutate(slbm = map2(.x = cols, .y = slbm , function(.x, .y){
      .x <- rename(.x, "Year" = sty, "sldata" = conc_slbm)
      bind_rows(.x, .y)

    })) %>% select(-cols)


  agg_bm <- right_join(dj_train, newsample, by = "duration")

  agg_bm_sleval <- right_join(sl_train, newsample, by = "duration")

  fitb <- gev.d.fit.sl.msfl.curv(agg_bm = agg_bm, leave_out = 0 , mult_sc = mult_sc,
                                 dur_offset = dur_offset, int_offset = int_offset, Maxit = Maxit)

  mlest <- bind_cols(fitb$mle, conv = c(NA, fitb$conv))

  if(!mult_sc) {
    mlest <- mlest %>% bind_cols( eta2 = 0)
  }
  if(!int_offset) {
    mlest <- mlest %>% bind_cols( tau = 0)
  }
  if(!dur_offset) {
    mlest <- mlest %>% bind_cols( theta = 0)
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

  zns  <- agg_bm_test %>%
    select(-Year)


  zns <- dplyr::rename(zns, "ds"= "duration")
  zns <- zns %>% dplyr::group_by(ds) %>% tidyr::nest(djdat = djdata)

  sl_train <- dplyr::rename(sl_train, "ds"= "duration")

  table_quants  <- dplyr::left_join(table_quants, zns, by = "ds")
  table_quants <- dplyr::left_join(table_quants, sl_train, by = "ds")

  table_quants <- table_quants %>%
    dplyr::mutate( quant_score_djeval = purrr::pmap_dbl( .l = list( .ps = ps , .q = est.quants, .obsdj = djdat),
                                                  .f = function(.ps, .q, .obsdj){
                                                    compute_qs(obs = unlist(.obsdj) , quant.est = .q,
                                                               p = .ps)
                                                  })) %>%
    dplyr::mutate( quant_score_sleval = purrr::pmap_dbl( .l = list( .ps = ps , .q = est.quants, .obssl = BMtest),
                                                         .f = function(.ps, .q, .obssl){
                                                           compute_qs(obs = .obssl , quant.est = .q,
                                                                      p = .ps)
                                                         })) %>%

    dplyr::arrange(ps, ds)

  return(table_quants)
}


cv_fit_eff <- function(agg_bm, conc_bm, leave_out , ds, quants, mult_sc = FALSE,  dur_offset = FALSE,
                       int_offset = FALSE, Maxit = 1500){

  agg_bm_test <- agg_bm %>% select(djbm) %>% unnest(cols = djbm) %>% filter(Year %in% leave_out)

  dj_train <- agg_bm %>% select(-slbm) %>%
    dplyr::mutate(djbm = map(djbm, ~ .x %>% filter(!(Year %in% leave_out))))

  minmaxn <- do.call( function(.x){data.frame(Min = min(.x), Max = max(.x),
                                              Nobs = length(.x))} , list(agg_bm[1,]$djbm[[1]]$Year))


  sl_lo <- data.frame(year1 = leave_out, year2 = leave_out -1)
  sl_lo <- sl_lo %>% mutate(year2 = ifelse(year2 >= minmaxn$Min, year2, year2 + minmaxn$Nobs))
  sl_train1 <- agg_bm %>% select( - djbm) %>%
    dplyr::mutate(slbm = map(slbm, ~ .x %>% filter(!(Year %in% unlist(sl_lo)))))

  conc_bm <- conc_bm %>% filter(sty %in% sl_lo$year2)  ## filter for possible new starting years

  # neither the new starting year nor the following year can be in the new sliding block
  conc_bm <- conc_bm %>%
    filter(!(sty %in% leave_out) ) %>%
    filter(!(year_f %in% leave_out))

  conc_bm <- conc_bm %>% mutate(Diff = year_f - sty) %>%
    mutate(Diff = ifelse(Diff >= 0, Diff, Diff + minmaxn$Nobs)) %>%
    group_by(sty) %>% filter(Diff == min(Diff)) %>% ungroup() %>% select(-Diff)

  conc_bm <- conc_bm %>% unnest(cols = newslbm) %>%
    select(-c( year_f)) %>%
    unnest(cols = conc_slbm) %>%
    group_by(duration) %>%
    nest(cols = c( sty , conc_slbm))  %>% ungroup()


  newsample <- conc_bm %>% left_join(sl_train1, by = "duration") %>%
    mutate(slbm = map2(.x = cols, .y = slbm , function(.x, .y){
      .x <- rename(.x, "Year" = sty, "sldata" = conc_slbm)
      bind_rows(.x, .y)

    })) %>% select(-cols)


  agg_bm <- right_join(dj_train, newsample, by = "duration")



  fitb <- gev.d.fit.sl.msfl.curv(agg_bm = agg_bm, leave_out = 0 , mult_sc = mult_sc,
                                 dur_offset = dur_offset, int_offset = int_offset, Maxit = Maxit)

  mlest <- bind_cols(fitb$mle, conv = c(NA, fitb$conv))

  if(!mult_sc) {
    mlest <- mlest %>% bind_cols( eta2 = 0)
  }
  if(!int_offset) {
    mlest <- mlest %>% bind_cols( tau = 0)
  }
  if(!dur_offset) {
    mlest <- mlest %>% bind_cols( theta = 0)
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

  zns  <- agg_bm_test %>%
    select(-Year)

  zns <- dplyr::rename(zns, "ds"= "duration")
  zns <- zns %>% dplyr::group_by(ds) %>% tidyr::nest(djdat = djdata)
  table_quants <- dplyr::left_join(table_quants, zns, by = "ds")
  table_quants <- table_quants %>%
    dplyr::mutate( quant_score = purrr::pmap_dbl( .l = list( .ps = ps , .q = est.quants, .obs = djdat),
                                                  .f = function(.ps, .q, .obs){
                                                    compute_qs(obs = unlist(.obs) , quant.est = .q,
                                                               p = .ps)
                                                  })) %>%
    dplyr::arrange(ps, ds)

  return(table_quants)
}
