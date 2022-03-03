
#' Cross-validation for d-GEV quantile estimation
#'
#' @param agg_bm Tibble as obtained from the output of [get_agg_bm()()].
#' @param conc_bm Tibble of sliding block maxima for all possible combinations of
#' years/seasons that are consecutive after removing observations from a year/season
#' for the test set
#' @param ds Vector of durations for which quantile of the d-GEV is estimated.
#' @param quants Vector of quantiles to estimate based on the d-GEV fit.
#' @param mult_sc logical. Whether to fit a multi-scale parameter eta2
#' @param dur_offset logical. Whether to fit a duration offset parameter theta
#' @param int_offset logical. Whether to fit an intensity offset parameter tau
#' @param testset Specifies how to construct the test set: either "consec", when test data
#' consists of maxima from consecutive years/seasons, such that the test sets are disjoint;
#' or 'random', when from all possible combinations of indices that can be left out, a random subset
#' of length 'n.cv' is sampled.
#' @param n.lo The amount of indices (i.e. how many years/seasons) you want to use to construct one test set
#' @param n.cv Only to be specified when 'testset = "random" ': specifies how many random test sets
#' to choose for cross-validation.
#' @param method  Specifies which estimation method to use: must be
#'  one of "disjoint", "sliding" or "both" (default).
#' @param returnLO Whether or not (default) the output should have an additional column containing the indices
#' of the years/seasons that the respective test set consists of. Note that the output already
#' has a colum with containing the test set (the values only, no indices).
#' @param optimMethod Optimisation method used in [optim()].
#' @param Maxit Passed to the control argument of [optim()]:
#' maximum number of iterations performed during the optimisation.
#' @param ... Further arguments that can be passed to the control argument of [optim()].
##'
#' @return Returns a tibble of dimension
#' \code{length(quants)}\eqn{ \times } \code{length(ds)} \eqn{\times} \eqn{n_{cv} \times 2 },
#' where \eqn{n_cv} is either the amount of possible testsets when the method for constructing the
#' test sets is `consec`, or as specified in the argument `n.cv` when `testset = "random" `.
#'  The tibble contains the estimated parameters, the estimation method,
#' a code for the convergence (zero indicates successful convergence), the values of \eqn{p} for
#' which quantiles are estimated, the values of the durations \eqn{d} for which
#' quantiles are estimated, the test set (data that is used for computing the quantile score)
#' and the quantile score. Quantile Scores can only be computed for those durations which are also
#' present in the testdata.
#'
#' @export
#'
#' @examples
#' dates <- seq(as.POSIXct("2000-01-01 00:00:00"),
#' as.POSIXct("2020-12-31 23:00:00"),by = 'hour')
#' prec <- rgamma(length(dates), shape = 0.1)
#' example_data <- data.frame(datetime = dates, prec = prec)
#'
#' agbm <- get_agg_bm(example_data, ds = c(1,2,4,8,16, 24, 48))
#' agdf <- fun_aggregate2df( example_data, ds = c(1,2,4,8,16, 24, 48) )
#' concbm <- compute_conc_bm_id(agdf)
#' cvres <- fit_dgev_cv_qs(agbm, concbm,  ds = c(1,4,8, 24))
#' # compute mean quantile scores based on cross-validation:
#' cvres %>% dplyr::select( - c( mut, sigma0, shape, eta, conv, Testdata )) %>%
#'           dplyr::group_by(ds, ps, estimator) %>%
#'           dplyr::summarise(MQS = mean(quant_score), .groups = "drop")
fit_dgev_cv_qs <- function(agg_bm, conc_bm, ds,
                           quants =  1-1/c(50, 100 ), mult_sc = FALSE,
                           dur_offset = FALSE, int_offset = FALSE, testset = "consec",
                           n.lo = 3, n.cv = NULL, method = "both", returnLO = FALSE,
                           optimMethod = "Nelder-Mead", Maxit = 1500, ... ){


  years_obs <- agg_bm[1, ]$djbm[[1]]$Year
  n.years <- length(years_obs)

  if(testset == "consec") {
    cvstart <- seq(years_obs[1], years_obs[n.years], n.lo)
    sbst <- purrr::map(cvstart, ~ .x + (0:(n.lo -1)))
    sbst[[length(sbst)]] <-  sbst[[length(sbst)]][ sbst[[length(sbst)]] <= years_obs[n.years]]
    sbst <- sbst[map(sbst, ~ length(.x)) == n.lo]

    n.cv <- length(sbst)
  } else if( testset == "random") {
    assertthat::assert_that(is.numeric(n.cv))

    all.subsets <- combn(years_obs, n.lo, simplify = F)
    sbst <- sample(all.subsets, min(n.cv, nrow(all.subsets)) )
    n.cv <- length(sbst)
  } else {
    stop("This is not a valid method to construct the test set. Choose either 'consec' or 'random'")
  }


   tibres <- tibble::tibble(left_out = sbst)  %>%
      dplyr::mutate( Result = purrr::map( left_out, .f = function(.x){
        tryCatch(
          dgev_cv_fit(agg_bm = agg_bm, conc_bm = conc_bm, leave_out = .x, ds =  ds,
                 quants = quants, mult_sc = mult_sc, dur_offset = dur_offset,
                 int_offset = int_offset, optimMethod = optimMethod, Maxit = Maxit, ... )  ,
          error = function(egal){
            tibble::tibble( mut = NA, sigma0 = NA, shape = NA,
                                         eta = NA,
                                         estimator = NA,  ps = NA, ds  = NA,
                                         est.quants = NA, Testdata = NA, quant_score = NA)

        } )
        }))

   if(!returnLO){
    tibres <-   tibres %>% dplyr::select(-left_out)
   }
   if(!(method == "both")){
     tibres <- tibres %>%
       dplyr::mutate( Result = purrr::map(Result, ~ .x %>% dplyr::filter( estimator == method)))
   }

   tibres %>% tidyr::unnest(cols = Result)

}
