#' Fit d-GEV based on sliding BM
#'
#' @param agg_bm Output of get_agg_bm, i.e. a tibble containing the sample of disjoint
#'  and sliding block maxima of the ID process,
#' for different durations $d$
#' @param mult_sc logical. Whether to fit a multi-scale parameter eta2
#' @param dur_offset logical. Whether to fit a duration offset parameter theta
#' @param int_offset logical. Whether to fit an intensity offset parameter tau
#' @param sigma0link The link function for \eqn{\sigma_0}. Default is \code{ make.link("log")}.
#' @param method The optimisation method used in \code{\link{optim}}.
#' @param Maxit Passed to the control argument of \code{\link{optim}}:
#' maximum number of iterations performed during the optimisation.
#' @param ... Further arguments that can be passed to the control argument of \code{\link{optim}}.
#' @return Returns a list with the following arguments:
#' * mle : The vector of estimated parameters
#' * conv: logical; whether optimisation succeeded
#' * nllh: value of the negative log-Likelihoods (disjoint and sliding)
#' * counts: Number of calls to the negative log-Likelihood
#' * dur_offset: logical; whether the duration offset was fitted
#' * mult_sc: logical; whether the duration offset was fitted
#' @export
#'
#' @examples
gev.d.fit.sl <- function(agg_bm,  mult_sc = TRUE, dur_offset = FALSE,
                                   int_offset = FALSE, sigma0link = make.link("log"),
                                    method = "Nelder-Mead", Maxit = 1000, ...){


  bm_dj  <- agg_bm %>% dplyr::select(djbm) %>% tidyr::unnest(cols = djbm)

  # delete dj bm from tibble for saving memory
  agg_bm_wy <- agg_bm %>% dplyr::select(-djbm)

  # compute initial parameters from dj bm
  optim_dj <- IDF::gev.d.fit(xdat  = bm_dj$djdata[!is.na(bm_dj$djdata)],
                             ds = bm_dj$duration[!is.na(bm_dj$djdata)],
                             sigma0link = sigma0link,
                             theta_zero = !dur_offset  ,
                             eta2_zero = !mult_sc,
                             tau_zero = !int_offset,
                             show = FALSE)
  par_init <- optim_dj$mle

  if( !dur_offset &  !mult_sc & !int_offset) {
    names( par_init) <- c("mut", "sigma0", "shape",  "eta")
  }
  if( !dur_offset &  !mult_sc & int_offset) {
    names( par_init) <- c("mut", "sigma0", "shape",  "eta", "tau")
  }
  if(!dur_offset  &  mult_sc & !int_offset){
    names( par_init) <- c("mut", "sigma0", "shape",  "eta", "eta2")
    par_init[["eta2"]] <- par_init[["eta2"]]  - par_init[["eta"]]
  }
  if(!dur_offset &  mult_sc & int_offset){
    names( par_init) <- c("mut", "sigma0", "shape",  "eta", "eta2", "tau")
    par_init[["eta2"]] <- par_init[["eta2"]]  - par_init[["eta"]]
  }
  if( dur_offset &  !mult_sc & !int_offset) {
    names( par_init) <- c("mut", "sigma0", "shape", "theta",  "eta")
  }
  if( dur_offset &  !mult_sc & int_offset) {
    names( par_init) <- c("mut", "sigma0", "shape", "theta",  "eta", "tau")
  }
  if(dur_offset  &  mult_sc & !int_offset){
    names( par_init) <- c("mut", "sigma0", "shape", "theta",  "eta", "eta2")
    par_init[["eta2"]] <- par_init[["eta2"]]  - par_init[["eta"]]
  }
  if(dur_offset & mult_sc & int_offset) {
    names( par_init) <- c("mut", "sigma0", "shape", "theta", "eta", "eta2", "tau")
    par_init[["eta2"]] <- par_init[["eta2"]]  - par_init[["eta"]]
  }



  # get tibble with unique values and their frequency

  bm_unique <- agg_bm_wy   %>%
    dplyr::mutate( uniq_data =
                     purrr::map(.x = slbm,
                                .f = function(.x){
                                  .df <- .x %>%
                                    dplyr::group_by(sldata) %>%
                                    dplyr::summarise( n = dplyr::n(), .groups = "drop")
                                }
                     )) %>%  dplyr::select(- slbm)


  bm_unique <- bm_unique %>% tidyr::unnest(cols = uniq_data)

  optim_sl <-  optim(par_init, fn = nllh.sl.d, bm_uni = bm_unique,
                     mult_sc = mult_sc,
                     dur_offset = dur_offset, int_offset = int_offset, sigma0link = sigma0link,
                     control = list(maxit = Maxit, ...),
                     method = method)

  par_init[["sigma0"]] <- sigma0link$linkinv(par_init[["sigma0"]])
  optim_sl$par[["sigma0"]] <- sigma0link$linkinv(optim_sl$par[["sigma0"]])

  par_est <- data.frame( t(par_init) ) %>%
    dplyr::bind_rows(data.frame(t(optim_sl$par))) %>%
    dplyr::bind_cols(estimator = c("disjoint", "sliding"))

  return(list( mle = par_est, conv = c(dj = optim_dj$conv, sl = optim_sl$convergence),
               nllh = data.frame(dj = optim_dj$nllh/length(bm_dj$djdata),
                                 sl = optim_sl$value),
               counts = data.frame(SL = optim_sl$counts),
               dur_offset = dur_offset,
               mult_sc = mult_sc))

}
