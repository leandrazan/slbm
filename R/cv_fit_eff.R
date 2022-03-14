#' Compute parameter estimates and quantile score on a test set
#' @description  Estimate quantiles based on sliding and/or disjoint blocks and
#' efficiently evaluate estimates based on the quantile score of a subset of the data
#'
#' @param djbm dataframe with columns 'obs' and 'Index' containing observed disjoint block
#' maxima and the year they were observed in, respectively. Needs to be passed
#' also when only using sliding blocks method to obtain the test set for cross-validation.
#' @param slbm dataframe  with columns 'sldata' and 'Index' containing observed sliding block
#' maxima and an Index for the block maxima, e.g. the year/season the first observation
#' of the respective sliding block was observerd in. With this Index, the
#' observations that are left out during cross-validation are specified. For example, if you want
#' to leave out all observations that were made in Years 1999, 2000, 2001, when computing
#' sliding block maxima, the column 'Index' must contain the year in which the first observation of the
#' sliding block was observed. When passing the years 1999, 2000, 2001 to 'leave_out', then
#' all observations from these years are left out, and sliding block maxima for the now consecutive
#' years (1998 is followd by 2002) are assembled from 'conc_bm'
#' @param conc_bm tibble containing sliding block maxima obtained from concatenating
#' observations from years/seasons... that are usually not consecutive, as computed from the function
#' 'compute_conc_bm'
#' @param leave_out The indices you want to leave out from the training data; observations
#' with these indices make up the test data
#' @param quants The quantiles which you want to estimate and compute the quantile score for
#' @param method Specifies which estimation method to use: must be
#'  one of "disjoint", "sliding" or "both" (default).
#' @param fixpar Optional: when applying the sliding block method, you can
#' fix one of the parameters to the respective estimated value obtained from the disjoint block method.
#' For the default (NULL) no parameter is fixed. Other options are "loc", "scale" or c("loc", "scale")
#' for fixing location, scale or both location and scale parameter. This overrides the choice of 'method'
#' and returns values of the disjoint method and the sliding method with fixed values if not NULL.
#' @return tibble containing the estimated parameters, the estimation method, the
#' estimated quantiles, the test set as well as the corresponding quantile scores.
#' @export
#'
#' @examples
#' df <- data.frame(obs = evd::rgpd(30*50), Index = rep(1:50, each = 30))
#' cv_fit_eff(data.frame( obs = blockmax(df$obs, r = 30, "disjoint"),  Index = 1:50),
#' data.frame( sldata = blockmax(df$obs, r = 30, "sliding"),
#' Index = rep(1:50, each = 30)[1:(49*30 +1)]),
#' conc_bm = compute_conc_bm(df, nlo = 3, blcksz = 30), leave_out = c(4, 35, 42))
#'
cv_fit_eff <- function(djbm, slbm, conc_bm, leave_out ,
                       quants = 1-1/c(2,5,10, 20, 50, 100, 200 ), method = "both",
                       fixpar = NULL ){

  bm_test <- djbm %>% dplyr::filter(Index %in% leave_out)

  dj_train <- djbm  %>% dplyr::filter(!(Index %in% leave_out))

  if(is.null(fixpar) & method == "disjoint"){
    fit.dj <- evd::fgev(dj_train$obs)
    est.par <- data.frame(t(fit.dj$estimate)) %>%  dplyr::bind_cols(estimator = "disjoint")
  } else {



    minmaxn <- do.call( function(.x){data.frame(Min = min(.x), Max = max(.x),
                                                Nobs = length(.x))} , list(djbm$Index))


    sl_lo <- data.frame(year1 = leave_out, year2 = leave_out -1)
    sl_lo <- sl_lo %>% dplyr::mutate(year2 = ifelse(year2 >= minmaxn$Min, year2, year2 + minmaxn$Nobs))


    sl_train1 <- slbm  %>%  dplyr::filter(!(Index %in% unlist(sl_lo)))

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
      dplyr::select(-c(ind2))   %>% dplyr::ungroup()

    conc_bm <- dplyr::rename(conc_bm, "Index" = ind1, "sldata" = newslbm)

    newsamp_sl <- dplyr::bind_rows(conc_bm, sl_train1)

    fit.dj <- evd::fgev(dj_train$obs)
    fit.sl <- evd::fgev(newsamp_sl$sldata)

    est.par <- data.frame(t(fit.dj$estimate)) %>%
      dplyr::bind_rows(data.frame(t(fit.sl$estimate)) )%>%
      dplyr::bind_cols(estimator = c("disjoint", "sliding"))

    if(!is.null(fixpar)){
      djparam <- fit.dj$estimate
      if( identical(fixpar , "loc")){
        fit.sl.fixloc <- evd::fgev(newsamp_sl$sldata, loc = djparam[1])
        fit.sl.fixloc <- data.frame( loc = djparam[1], t(fit.sl.fixloc$estimate))

        fix.par <- data.frame(fit.sl.fixloc) %>%
          dplyr::bind_cols(estimator =  "sl.locfix")

      }
      if( identical(fixpar, "scale")){
        fit.sl.fixscale<- evd::fgev(newsamp_sl$sldata, scale = djparam[2])
        fit.sl.fixscale <- data.frame( scale = djparam[2], t(fit.sl.fixscale$estimate))

        fix.par <- data.frame(fit.sl.fixscale) %>%
          dplyr::bind_cols(estimator = "sl.scalefix")
      }
      if( identical(fixpar, c("loc", "scale") )){
        fit.sl.fixscale <- evd::fgev(newsamp_sl$sldata, scale = djparam[2])
        fit.sl.fixloc <- evd::fgev(newsamp_sl$sldata, loc = djparam[1])
        fit.sl.fixloc <- data.frame( loc = djparam[1], t(fit.sl.fixloc$estimate))
        fit.sl.fixscale <- data.frame( scale = djparam[2], t(fit.sl.fixscale$estimate))

        fix.par <- data.frame(fit.sl.fixscale) %>%
          dplyr::bind_rows(data.frame(fit.sl.fixloc) )%>%
          dplyr::bind_cols(estimator = c("sl.scalefix", "sl.locfix"))
      }

      est.par <- est.par %>% dplyr::bind_rows(fix.par)

    }
  }

  est.par <- est.par %>% tidyr::expand_grid(ps = quants)

  est.par <- est.par %>%
    dplyr::mutate( est.quants = purrr::pmap_dbl( list(loc, scale, shape, ps) ,
                                                 .f = function( .loc, .scale, .shape, .ps){
                                                   evd::qgev(.ps, loc = .loc, scale = .scale, shape = .shape)
                                                 }
    ))

  zns  <- bm_test

  est.par$djdat <- rep(list(zns), nrow(est.par))

  table_quants <- est.par %>%
    dplyr::mutate( quant_score = purrr::pmap_dbl( .l = list( .ps = ps , .q = est.quants),
                                                  .f = function(.ps, .q){
                                                    compute_qs(obs = zns$obs, quant.est = .q,
                                                               p = .ps)
                                                  })) %>%
    dplyr::arrange(ps)


  if(method == "sliding" & is.null(fixpar)){
    table_quants %>% dplyr::filter( estimator == "sliding")
  } else {
    table_quants
  }




}
