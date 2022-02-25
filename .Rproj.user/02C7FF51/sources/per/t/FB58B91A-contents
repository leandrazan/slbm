


#' Cross-validation for quantile estimation
#' @description  Fit a GEV distribution to block maxima of a univariate time series,
#'  estimate quantiles and compute quantile scores based on cross-validation
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
#' @param method Specifies which estimation method to use: must be
#'  one of "disjoint", "sliding" or "both" (default).
#' @param testset Specifies how to construct the test set: either "consec", when test data
#' consists of maxima from consecutive years/seasons, such that the test sets are disjoint;
#' or 'random', when from all possible combinations of indices that can be left out, a random subset
#' of length 'n.cv' is sampled.
#' @param quants The quantiles which you want to estimate and compute the quantile score for.
#' @param n.lo The amount of indices (i.e. how many years/seasons) you want to use to construct one test set
#' @param n.cv Only to be specified when 'testset = "random" ': specifies how many random test sets
#' to choose for cross-validation.
#' @param fixpar Optional: when applying the sliding block method, you can
#' fix one of the parameters to the respective estimated value obtained from the disjoint block method.
#' For the default (NULL) no parameter is fixed. Other options are "loc", "scale" or c("loc", "scale")
#' for fixing location, scale or both location and scale parameter. This overrides the choice of 'method'
#' and returns values of the disjoint method and the sliding method with fixed values if not NULL.
#' @return Tibble containing the estimated parameters, the estimation method, the
#' estimated quantiles, the test set as well as the corresponding quantile scores for every
#' cross-validation test set.
#' @example
#' df <- data.frame(obs = evd::rgpd(30*50), Index = rep(1:50, each = 30))
#' concbm <- compute_conc_bm(df, nlo = 3, blcksz = 30)
#' djbm <- data.frame( obs = blockmax(df$obs, r = 30, "disjoint"),  Index = 1:50)
#' slbm <- data.frame( sldata = blockmax(df$obs, r = 30, "sliding"),
#'                           Index = rep(1:50, each = 30)[1:(49*30 +1)])
#' cvfits <- fit_gev_cv_qs(djbm = djbm, slbm = slbm, conc_bm = concbm, testset = "consec" ,
#' quants = 1-1/c(2,5,10, 20, 50, 100, 200 ), n.cv = NULL,
#' fixpar = NULL )
#'
#' # Now you can compute the mean quantile scores of the cross-validated values:
#' cvfits %>% dplyr::select(-c(djdat, est.quants)) %>%
#' dplyr::group_by(estimator, ps) %>%
#' dplyr::summarise( MeanQS = mean(quant_score),  .groups = "drop")

#'
fit_gev_cv_qs <- function(djbm, slbm, conc_bm, method = "both", testset = "consec" ,
                       quants = 1-1/c(2,5,10, 20, 50, 100, 200 ),n.lo = 3,  n.cv = NULL,
                       fixpar = NULL){
  n.years <- nrow(djbm)
#
  if(testset == "consec") {
    cvstart <- seq(1, n.years, n.lo)
    sbst <- purrr::map(cvstart, ~ .x + (0:(n.lo -1)))
    sbst[[length(sbst)]] <-  sbst[[length(sbst)]][ sbst[[length(sbst)]] <= n.years]

    n.cv <- length(sbst)
  } else if( testset == "random") {
    assertthat::assert_that(is.numeric(n.cv))

    all.subsets <- combn(n.years, n.lo, simplify = F)
    sbst <- sample(all.subsets, min(n.cv, nrow(all.subsets)) )
    n.cv <- length(sbst)
  } else {
    stop("This is not a valid method to construct the test set. Choose either 'consec' or 'random'")
  }

  results_cv <- purrr::map_dfr(sbst, ~
                                 tryCatch(
                                   cv_fit_eff(djbm = djbm, slbm = slbm,
                                              conc_bm = conc_bm, leave_out = .x, fixpar = fixpar,
                                              method = method),
                                   error = function(egal){tibble::tibble()}))


    results_cv

}
