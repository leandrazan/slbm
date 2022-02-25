#' Estimate quantiles and compute quantile score based on cross-validation
#'
#' @description This function split the data into a training and a test set,
#' fits a GEV distribution to block maxima of the training set, estimates quantiles
#' based on the fit and then computes the quantile score on the test set.
#' GEV distributions are fitted based on both disjoint and sliding block maxima.
#'
#' @param data The daily (hourly, subhourly ...) data from which block maxima are computed
#' @param djbm dataframe with columns 'obs' and 'Index' containing observed disjoint block
#' maxima and and Index (e.g. the year/season in which the observation was made)
#' @param leave_out The indices you want to leave out from the training data; observations
#' with these indices make up the test data
#' @param blcksz The blocksize for computing block maxima
#' @param ps the quantiles which you want to estimate and compute the quantile score for
#'
#' @return tibble containing the estimated parameters, the estimation method, the
#' estimated quantiles, the test set as well as the corresponding quantile scores.
#' For a large number of cross-validations, this function is very slow.
#' Use 'cv_fit_eff' in that case, which is more efficient.
#' @export
#'
#' @examples
#' df <- data.frame(obs = evd::rgpd(30*50), Index = rep(1:50, each = 30))
#' fit_cv_iter(df, data.frame( obs = blockmax(df$obs, r = 30, "disjoint"), Index = 1:50),
#' leave_out = c(4, 35, 42),  blcksz = 30)
fit_cv_slow <- function(data, djbm, leave_out, blcksz,
                        ps =  1-1/c(2,5,10, 20, 50, 100, 200 )){

  data.train <- data %>% dplyr::filter(! (Index %in% leave_out))



  bm.dj.train <- (djbm %>% dplyr::filter(!(Index %in% leave_out)))$obs
  bm.dj.test <- (djbm %>% dplyr::filter(Index %in% leave_out))$obs
  bm.sl.train <- blockmax(data.train$obs, r = blcksz, "sliding")

  fit.dj <- evd::fgev(bm.dj.train)
  fit.sl <- evd::fgev(bm.sl.train)

  est.par <- data.frame(t(fit.dj$estimate)) %>%
    dplyr::bind_rows(data.frame(t(fit.sl$estimate)) )%>%
    dplyr::bind_cols(estimator = c("disjoint", "sliding"))

  est.par <- est.par %>% tidyr::expand_grid(ps = ps)

  est.par <- est.par %>%
    dplyr::mutate( est.quants = purrr::pmap_dbl( list(loc, scale, shape, ps) ,
                                                 .f = function( .loc, .scale, .shape, .ps){
                                                   evd::qgev(.ps, loc = .loc, scale = .scale, shape = .shape)
                                                 }
    ))

  zns  <- bm.dj.test
  est.par$djdat <- rep(list(zns), nrow(est.par))

  table_quants <- est.par %>%
    dplyr::mutate( quant_score = purrr::pmap_dbl( .l = list( .ps = ps , .q = est.quants),
                                                  .f = function(.ps, .q){
                                                    compute_qs(obs = zns, quant.est = .q,
                                                               p = .ps)
                                                  })) %>%
    dplyr::arrange(ps)


  table_quants


}
