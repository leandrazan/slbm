#' Compute unique block maxima for bootstrap
#'
#' @param data tibble with columns named "Station" and "Obs", where the latter
#' contains daily observations
#' @param blcksz The blocksize parameter
#' @param indexblock A matrix with columns `blockind` and `obsind`, where `obsind` is the
#' consecutive counter of the observations (from to `nrow(data)`) and `blockind` is the
#' index of the K-block the observation belongs to (here, K is the blockbootstrap parameter).
#' @param temp_cvrt either NULL or the temporal covariate
#' @param looplastblock logical; whether within the blocks of size `K` times `blcksz` the
#' oobservations are looped (first block concatenated to last block)
#' @param returnfullsamp logical; whether the full (looped) sliding block maxima sample shall
#' be returned. Needed for parametric covariance estimation.
#' @param K Number of disjoint blocks that make up
#'
#' @return Returns a tibble with columns `Station`, `Kblockind`,`uniq_data` and, if `returnfullsamp  = TRUE`,
#' `full_data`. These are the station reference, the block index of the blocks of size `K` times `blcksz`,
#' the sample of unique sliding block maxima values and their frequency and the looped sliding block maxima.
#' @export
#'
#' @examples
#' blcksz <- 90
#' ny <- 50
#' xx <- evd::rgpd(ny*blcksz, shape = 0.2)
#' df.xx <- data.frame(Station = "X1", Obs = xx)
#' k <- 4
#' ndata <- blcksz * ny
#' .nKblocks <- ceiling(ndata/(k*blcksz))
#' indexblock <- data.frame(blockind = c(rep(1:(.nKblocks-1), each = k*blcksz),
#'                                       rep(.nKblocks, ndata - k*blcksz*(.nKblocks-1))),
#'                          obsind = 1:ndata)
#'
#' sluniq_wb <- get_uniq_bm_boot(df.xx, blcksz = blcksz, K = k, indexblock = indexblock,
#'                               temp_cvrt = NULL, looplastblock = TRUE,
#'                               returnfullsamp = TRUE)
#'
get_uniq_bm_boot <- function (data, blcksz, K,  indexblock, temp_cvrt = NULL,
                              looplastblock = TRUE, returnfullsamp = FALSE) {


  ndat <- nrow(data)
  # fill up last Kblock if it contains less observations than the other blocks
  if(!(ndat/(K*blcksz)  == floor(ndat/(K*blcksz)))) {
    m <- ceiling(ndat/blcksz)
    mk <- floor(m/K)
    diffmk <- K - m + mk*K
    data <- data %>% dplyr::bind_rows(data[(ndat - diffmk*blcksz +1):ndat,  ])
    indexblock <- indexblock %>% dplyr::bind_rows(data.frame(blockind = rep(mk +1, diffmk*blcksz ), obsind = ndat + (1:(diffmk*blcksz))))
  }
  data$Kblockind <- indexblock$blockind

  if (!is.null(temp_cvrt)) {
    n.cvrt <- length(temp_cvrt)
    # add column containing the temporal covariate. If temporal covariate is given for
    # sliding bm only, repeat the last value to obtain temporal covariate for each day
    data$tempcvrt <- c(temp_cvrt, rep(temp_cvrt[n.cvrt], (nrow(data) - n.cvrt)))

    if(returnfullsamp) {
      # return the sliding BM sample itself and the weighted sliding BM sample

      uu <- data %>% dplyr::group_by(Station, Kblockind) %>% tidyr::nest() %>%
        dplyr::mutate(full_data = purrr::map(.x = data, .f = function(.x) {
          bmx <- slbm::blockmax(.x$Obs, r = blcksz, method = "sliding",
                                looplastblock = looplastblock)

          bmx <- data.frame(slbm = bmx, temp_cvrt = .x$tempcvrt[1:length(bmx)])
          bmx}
        ))

      uu1 <- uu %>%
        dplyr::mutate(uniq_data = purrr::map(full_data,  .f = function(.x) {
          bmx <- .x %>% dplyr::group_by(slbm, temp_cvrt) %>%
            dplyr::summarise(n = dplyr::n(), .groups = "drop")
          bmx
        })) %>% dplyr::ungroup() %>% dplyr::select(-data)

      return(uu1)

    }
    else {
      # return only the weighted sliding BM sample
      data %>% dplyr::group_by(Station, Kblockind) %>% tidyr::nest() %>%
        dplyr::mutate(uniq_data = purrr::map(.x = data, .f = function(.x) {

          bmx <- slbm::blockmax(.x$Obs, r = blcksz, method = "sliding",
                                looplastblock = looplastblock)

          bmx <- data.frame(slbm = bmx, temp_cvrt = .x$tempcvrt[1:length(bmx)])
          bmx <- bmx %>% dplyr::group_by(slbm, temp_cvrt) %>%
            dplyr::summarise(n = dplyr::n(), .groups = "drop")
          # bmx$n[which.max(bmx$slbm)] <- blcksz
          bmx
        })) %>% dplyr::ungroup() %>% dplyr::select(-data)
    }
  }

  else {
    if(returnfullsamp) {

      uu <- data %>% dplyr::group_by(Station, Kblockind) %>% tidyr::nest() %>%
        dplyr::mutate(full_data = purrr::map(.x = data, .f = function(.x) {

          bmx <- slbm::blockmax(.x$Obs, r = blcksz, "sliding", looplastblock = looplastblock)

          bmx <- data.frame(slbm = bmx)
          bmx}))
      uu1 <- uu %>%
        dplyr::mutate(uniq_data = purrr::map(full_data,  .f = function(.x) {
          bmx <- .x %>% dplyr::group_by(slbm) %>%
            dplyr::summarise(n = dplyr::n(), .groups = "drop")
          bmx
        })) %>% dplyr::ungroup() %>% dplyr::select(-data)

      return(uu1)


    } else {
      data %>% dplyr::group_by(Station, Kblockind) %>% tidyr::nest() %>%
        dplyr::mutate(uniq_data = purrr::map(.x = data, .f = function(.x) {

          bmx <- slbm::blockmax(.x$Obs, r = blcksz, "sliding", looplastblock = looplastblock)

          bmx <- data.frame(slbm = bmx)
          bmx <- bmx %>% dplyr::group_by(slbm) %>%
            dplyr::summarise(n = dplyr::n(), .groups = "drop")
          bmx
        })) %>% dplyr::ungroup() %>% dplyr::select(-data)
    }
  }
}

#' Compute RL of (fitted) GEV distribution
#'
#' @param theta Estimated GEV parameters
#' @param Tyrl The return periods for which the RLs are computed
#' @param type The model: 'stationary'. 'shift' or 'scale'
#' @param ref_gmst The reference value of the temporal covariate for which the
#' RL is computed. Ignored when \code{type = "stationary"}.
#' @param ... For additional argument `rel_trend` for specifying the parametrisation which
#'  needs to be passed when `type = "scale"`
#' @return A tibble containg the estimated RL for each combination of \code{Tyrl} and
#' \code{ref_gmst}.
#' @export
#'
#' @examples
#' compute_rl(theta = c(loc0 = 7, scale0 = 2, shape = 0.2, TempLoc1 = 2),
#' Tyrl = c(50, 100), type = "shift", ref_gmst = c(0.5, 0.95))
#'
compute_rl <- function(theta, Tyrl, type, ref_gmst = NULL, ...) {
  add.args <- list(...)

  if(type == "stationary") {
    rl <- evd::qgev(1-1/Tyrl, loc = theta[1], scale = theta[2], shape = theta[3])
    return(dplyr::tibble(rl = rl, Year = Tyrl ))
  }
  if(type == "scale") {
    if(is.null(add.args$rel_trend)) {stop("Please specify the parametrisation of the scale model.")}
    if(add.args$rel_trend) {
      mut <- theta[1]*exp(theta[4]/theta[1]*ref_gmst)
      sigmat <- theta[2]*exp(theta[4]/theta[1]*ref_gmst)
    }
    else {
      mut <- theta[1]*exp(theta[4]*ref_gmst)
      sigmat <- theta[2]*exp(theta[4]*ref_gmst)
    }
    return(purrr::map_dfr(Tyrl, ~ {
      rl <- evd::qgev(1-1/.x, loc = mut, scale = sigmat, shape = theta[3])
      dplyr::tibble(rl = rl, refGMST = ref_gmst, Year = .x)
    }))
  }
  if(type == "shift") {
    mut <- theta[1] + theta[4]*ref_gmst
    return(purrr::map_dfr(Tyrl, ~ {
      rl <- evd::qgev(1-1/.x, loc = mut, scale = theta[2], shape = theta[3])
      dplyr::tibble(rl = rl, refGMST = ref_gmst, Year = .x)
    }))
  }

}

#'  Studentized bootstrap repetition
#'
#' @param nKblocks Indices of blocks that make up the bootstrap sample
#' @param sluniq Tibble with column `Station` and `uniq_data`, where the latter
#' contains a tibble with columns `Kblockind`, `slbm`, `n` giving the index of the
#' K-block, the unqiue sliding block maxima within that block and their frequency, respectively.
#' @param slorig The complete sliding block maxima sample for the looped blocks.
#' @param blcksz The blocksize parameter.
#' @param start_vals Starting values.
#' @param scale.link The link function for the scale parameter.
#' @param estimate_RL Must be TRUE. GEV parameter Estimation not implemented atm.
#' @param seed A seed for the sampling of blocks.
#' @param reltol Passed to optim.
#' @param Tyrl The return periods for which RLs are estimated
#' @param ref_gmst Vector of reference values of the temporal covariate for which
#' RLs are estimated. Ignored when type is 'stationary'.
#' @param varmeth The method used to estimate the covariance matrix of GEV parameters: 'V' or 'V2'.
#' @param type The type of the model to fit: one of `stationary`, `shift` or `scale`
#' @param chain logical and only relevant for non-stationary models: Whether the method for
#'  covariance estimation should be based on chain rule
#' @param ... Further arguments passed to `fit_gev_univ`. When `type = scale`, one must
#'  pass the logical argument `rel_trend` which specifies the parametrisation of the scale model.
#'
#' @return A list containing
#'
#' * `est_boot` the parameter estimates of the bootstrap sample
#' and one or both of the following:
#' * `estvarV` the estimated covariance matrix of the parameter estimates (method 1)
#' * `estvarV2` the estimated covariance matrix of the parameter estimates (method 1)
#'
#' when \code{stimate_RL = FALSE}. If \code{estimate_RL = TRUE}, a list with the following
#' components is returned:
#'
#' * `rlboot` the estimated RLs
#' and one or both of the following:
#' * `rlvarV` the estimated variances of the RLs (method 1)
#' * `rlvarV2` the estimated variances of the RLs (method 2)
#' @export
#'
#' @examples
#' # none
sample_boot_sl <- function(nKblocks, sluniq, slorig,
                           blcksz,
                           start_vals = NULL,
                           type = "stationary",
                           scale.link = make.link("identity"),
                           seed, reltol = 1e-08, estimate_RL = TRUE,
                           Tyrl = c(50, 100),
                           ref_gmst = NULL,
                           varmeth = "V2", chain = TRUE,
                           shape_constraint = -0.9,  ...) {


  set.seed(seed)

  add.args <- list(...)
  # the nKblocks blocks of size Kblock that make up the bootstrap sample
  chosenblocks <- sample(1:nKblocks, nKblocks, replace = TRUE)

  # the bootstrap sample. sluniq consists of columns containing the Kblockindex,
  # the unique values of sliding block maxima and their resp. frequency.
  # observations with Kblockind in chosenblocks are filtered
  boot_samp <- sluniq[chosenblocks, ] %>% tidyr::unnest(cols = uniq_data)

  slorig <- purrr::map_dfr(slorig[chosenblocks, ]$full_data, ~ .x)
  # slorig <- unlist((slorig  %>% dplyr::filter(Kblockind %in% chosenblocks))$full_data)


  # the ML estimator based on weighted log-Likelihood computed on the bootstrap sample along with estimate of covariance matrix
  est_boot <-  slbm::fit_gev_univ(data = boot_samp, hessian = TRUE, type  = type, return_cov = TRUE, varmeth = varmeth,
                                  chain = chain, rel_trend = add.args$rel_trend,
                                  orig_slbm = slorig$slbm, orig_cvrt = slorig$temp_cvrt, blcksz = blcksz)

  if(est_boot$mle["shape"] <= shape_constraint) { return(list(rl_boot = tibble::tibble(rl = NA, refGMST = NA, Year = NA),
                                                              rlvarV2 = tibble::tibble(var_hat = NA, refGMST = NA, Year = NA),
                                                              mle = c("mu" = NA, "sigma" = NA, "shape" = NA, "alpha" = NA),
                                                              conv = NA,
                                                              hessian = array(dim = c(4,4)), V2 = array(dim = c(4,4))))}
  else {
    if(estimate_RL == FALSE) {
      return(est_boot)
     # return(append( est_boot, estvar))
    }
    else {

      rlhat <- slbm::compute_rl(theta = est_boot$mle, Tyrl = Tyrl, type = type, ref_gmst = ref_gmst, rel_trend = add.args$rel_trend)
      if(varmeth %in% c("both", "V2")) {
        varrl_V2 <- purrr::map_dfr(Tyrl, ~ {
          dplyr::tibble(var_hat =
                   slbm::estimate_var_rl(theta = est_boot$mle, Tyrl = .x, type = type,
                                         ref_gmst = ref_gmst, Covmat = est_boot$V2, rel_trend = add.args$rel_trend),
                 refGMST = ref_gmst, Year = .x) })
      }
      if(varmeth %in% c("both", "V")) {
        varrl_V <- purrr::map_dfr(Tyrl, ~ {
          dplyr::tibble(var_hat =
                   slbm::estimate_var_rl(theta = est_boot$mle, Tyrl = .x, type = type,
                                         ref_gmst = ref_gmst, Covmat = est_boot$V, rel_trend = add.args$rel_trend),
                 refGMST = ref_gmst, Year = .x) })
      }
      if(isTRUE(estimate_RL)) {
        if(varmeth == "V") {
          return(list(rl_boot = rlhat,  rlvarV = varrl_V))
        }
        if(varmeth == "V2") { return(list(rl_boot = rlhat, rlvarV2 = varrl_V2)) }
        else {return(list(rl_boot = rlhat, rlvarV = varrl_V, rlvarV2 = varrl_V2)) }
      }
      else if(estimate_RL == "both") {

          if(varmeth == "V") {
            return(append(list(rl_boot = rlhat,  rlvarV = varrl_V), est_boot))
          }
          if(varmeth == "V2") { return(append(list(rl_boot = rlhat, rlvarV2 = varrl_V2), est_boot)) }
          else {return(append(list(rl_boot = rlhat, rlvarV = varrl_V, rlvarV2 = varrl_V2), est_boot)) }

      }

    }
  }

}


#' Studentized block bootstrap for sliding block maxima
#'
#' @param sluniq Sample of unique sliding block maxima with repsective frequency and
#' corresping covariate value (if type is not stationary)
#' @param slorig Sample of original sliding block maxima with
#' corresping covariate value (if type is not stationary)
#' @param nKblocks Number of Kblocks from which to bootstrap
#' @param blcksz The blocksize parameter
#' @param B The number of bootstrap repetitions
#' @param start_vals Starting values for the optimisation of the negative log-Likelihood of
#' the bootstrap sample.
#' @param type One of 'stationary', 'shift' or 'scale'.
#' @param scale.link Link function. Not implemented atm.
#' @param reltol Passed to optim.
#' @param estimate_RL One of \code{TRUE, FALSE} or \code{"both"}. If \code{FALSE},
#' only GEV parameter estimates are computed. If \code{TRUE}, only RL estimates are computed.
#' If \code{"both"}, both parameter and RL estimates are computed.
#' @param Tyrl The return periods for which RLs are estimated
#' @param ref_gmst Vector of reference values of the temporal covariate for which
#' RLs are estimated. Ignored when type is 'stationary'.
#' @param varmeth The method used to estimate the covariance matrix of GEV parameters: 'V', 'V2', or 'both'.
#' @param chain logical and only relevant for non-stationary models: Whether the method for
#'  covariance estimation should be based on chain rule
#' @param ... Further arguments passed to `fit_gev_univ`. When `type = scale`, one must
#'  pass the logical argument `rel_trend` which specifies the parametrisation of the scale model.
#' @return
#' @export
#'
#' @examples
#' blcksz <- 90
#' ny <- 100
#' xx <-  evd::rgpd(ny*blcksz, shape = -0.2) + 2*rep(1:ny/ny, each = blcksz)
#' df.xx <- data.frame(Station = "X1", Obs = xx)
#' k <- 3
#' ndata <- blcksz * ny
#' .nKblocks <- ceiling(ndata/(k*blcksz))
#' indexblock <- data.frame(blockind = c(rep(1:(.nKblocks-1), each = k*blcksz),
#'                                       rep(.nKblocks, ndata - k*blcksz*(.nKblocks-1))),
#'                          obsind = 1:ndata)
#'
#' sluniq_wb <- get_uniq_bm_boot(df.xx, blcksz = blcksz, K = k, indexblock = indexblock,
#'                               temp_cvrt = rep(1:ny/ny, each = blcksz), looplastblock = TRUE,
#'                               returnfullsamp = TRUE)
#'
#' full_slbm <- sluniq_wb %>% dplyr::select(-uniq_data)
#' sluniq_wb <- sluniq_wb %>% dplyr::select(-full_data)
#'
#'
#' estim_lo <- slbm::fit_gev_univ(tidyr::unnest(sluniq_wb, cols = uniq_data),
#'  type = "shift", hessian = TRUE)
#'
#' boot_sl(sluniq = sluniq_wb, slorig = full_slbm, blcksz = blcksz, nKblocks = .nKblocks,
#' B = 10, type = "shift", ref_gmst = (c(0.8, 0.9)), Tyrl = c(50, 100), start_vals = estim_lo$mle,
#' scale.link = make.link("identity"), reltol = 1e-09, estimate_RL = TRUE)
boot_sl <- function(sluniq, slorig, nKblocks,
                    blcksz, B = 50,
                    start_vals = NULL,
                    type = "stationary",
                    scale.link = scale.link, reltol = 1e-08,
                    estimate_RL = TRUE,
                    Tyrl = c(50, 100), ref_gmst = 0.8, varmeth = "V2", chain = TRUE, ...) {

  add.args <- list(...)

  # perform bootstrap B times
  boot_estimators <- purrr::map(1:B,
                                ~ slbm::sample_boot_sl(nKblocks = nKblocks,
                                                 sluniq = sluniq,
                                                 slorig = slorig,
                                                 blcksz = blcksz,
                                                 start_vals = start_vals,
                                                 type = type,
                                                 scale.link = scale.link,
                                                 seed = .x,
                                                 # reltol = 1e-08,
                                                 estimate_RL = estimate_RL,
                                                 Tyrl = Tyrl,
                                                 ref_gmst = ref_gmst,
                                                 varmeth = varmeth, chain = chain,
                                                 rel_trend = add.args$rel_trend))

  browser()
  if(estimate_RL == FALSE) {
   if(varmeth == "both") {
      boot_parest <- list(resV =
                            purrr::map_dfr(boot_estimators, ~ dplyr::tibble(
                              data.frame(t(.x$mle), conv = .x$conv),
                              CovestV = list(.x$V))),
                          resV2 = purrr::map_dfr(boot_estimators, ~ dplyr::tibble(
                            data.frame(t(.x$mle), conv = .x$conv),
                            CovestV = list(.x$V2))))
   } else if(varmeth == "V") {
     boot_parest <- list(resV = purrr::map_dfr(boot_estimators, ~ dplyr::tibble(
       data.frame(t(.x$mle), conv = .x$conv),
       CovestV = list(.x$V))))
   } else {
     boot_parest <- list(resV2 = purrr::map_dfr(boot_estimators, ~ dplyr::tibble(
       data.frame(t(.x$mle), conv = .x$conv),
       CovestV = list(.x$V2))))
   }

  }
  else if(estimate_RL == "both") {

    if(varmeth %in%  c("both", "V2")) {
      if(type == "stationary") {
        boot_parest_v2 <- purrr::map_dfr(boot_estimators, ~
                                           { resrl <- dplyr::left_join(.x$rl_boot, .x$rlvarV2,
                                                                       by = c("Year"))
                                           respars <- dplyr::tibble(
                                             data.frame(t(.x$mle), conv = .x$conv),
                                             CovestV = list(.x$V2))
                                           respars %>% dplyr::bind_cols(resrl) }
        )
      }
      else {
        boot_parest_v2 <- purrr::map_dfr(boot_estimators, ~
                                         { resrl <- dplyr::left_join(.x$rl_boot, .x$rlvarV2,
                                                                     by = c("refGMST", "Year"))
                                           respars <- dplyr::tibble(
                                            data.frame(t(.x$mle), conv = .x$conv),
                                           CovestV = list(.x$V2))
                                            respars %>% dplyr::bind_cols(resrl) }
                                         )
      }
    }
    if (varmeth %in%  c("both", "V")) {
      if(type == "stationary") {
        boot_parest_v <- purrr::map_dfr(boot_estimators, ~
                                           { resrl <- dplyr::left_join(.x$rl_boot, .x$rlvarV,
                                                                       by = c("Year"))
                                           respars <- dplyr::tibble(
                                             data.frame(t(.x$mle), conv = .x$conv),
                                             CovestV = list(.x$V))
                                           respars %>% dplyr::bind_cols(resrl) }
        )
      }
     else {
      boot_parest_v <- purrr::map_dfr(boot_estimators, ~
                                        {resrl <- dplyr::left_join(.x$rl_boot, .x$rlvarV,
                                                                    by = c("refGMST", "Year"))
                                          respars <- dplyr::tibble(
                                            data.frame(t(.x$mle), conv = .x$conv),
                                            CovestV = list(.x$V))
                                          respars %>% dplyr::bind_cols(resrl)
                                        }
                                        )
     }
    }
    if(varmeth == "both") {
      boot_parest <- list( resV = boot_parest_v, resV2 = boot_parest_v2)
    }
    else {
      boot_parest <- base::switch(varmeth, "V" = list(resV = boot_parest_v),
                            "V2" = list(resV2 = boot_parest_v2))
    }
  }
  else if(isTRUE(estimate_RL)) {

    if(varmeth %in%  c("both", "V2")) {
      if(type == "stationary") {
        boot_parest_v2 <- purrr::map_dfr(boot_estimators, ~
                                           { dplyr::left_join(.x$rl_boot, .x$rlvarV2,
                                                              by = c("Year"))
                                           }
        )
      }
      else {
        boot_parest_v2 <- purrr::map_dfr(boot_estimators, ~
                                           { dplyr::left_join(.x$rl_boot, .x$rlvarV2,
                                                                       by = c("refGMST", "Year"))
                                            }
        )
      }


    }
    if (varmeth %in%  c("both", "V")) {
      if(type == "stationary") {
        boot_parest_v <- purrr::map_dfr(boot_estimators, ~
                                           { dplyr::left_join(.x$rl_boot, .x$rlvarV,
                                                              by = c("Year"))
                                           }
        )
      }
      else {
        boot_parest_v <- purrr::map_dfr(boot_estimators, ~
                                          { dplyr::left_join(.x$rl_boot, .x$rlvarV,
                                                                     by = c("refGMST", "Year"))
                                           }
        )
      }
    }
    if(varmeth == "both") {
      boot_parest <- list( resV = boot_parest_v, resV2 = boot_parest_v2)
    } else {
      boot_parest <- base::switch(varmeth, "V" = list(resV = boot_parest_v),
                                  "V2" = list( resV2 = boot_parest_v2))
    }

  }



  boot_parest

}

errfct <- function(type) {
  if(type == "stationary") {
    list(mle = c(loc = NA, scale = NA, shape = NA),
         hessian = NA)
  } else {
    list(mle = c(mu = NA, sigma = NA, shape = NA,
                 alpha = NA),
         hessian = NA)
  }
}

#' Estimate RLs and CIs based on sliding blocks
#'
#' @description Estimate RLs and respective confidence intervals of GEV fits based on
#' sliding Block maxima, with possibility to consider a shift or a scaling of maxima
#' with respect to a temporal covariate.
#' Confidence intervals are obtained via block bootstrap.
#'
#' @param x Numeric vector of data with daily resolution
#' @param Kblock Value(s) of number of blocks making up one K-block within the block-bootstrap
#' @param B Number of bootstrap iterations
#' @param temp.cov temporal covariate, can be NULL
#' @param scale.link link function for the scale parameter
#' @param looplastblock logical; wether last block is looped
#' @param reltol Passed to optim
#' @param blcksz The blocksize parameter
#' @param type Either one of 'stationary' (no trend), 'shift' (linear trend in location parameter)
#' or 'scale' (trend in location and scale parameter such that their ratio is constant over time)
#' @param Tyrl The return periods for which RLs are estimated
#' @param ref_gmst Vector of reference values of the temporal covariate for which
#' RLs are estimated. Ignored when type is 'stationary'.
#' @param estimate_RL One of \code{TRUE, FALSE} or \code{"both"}. If \code{FALSE},
#' only GEV parameter estimates are computed. If \code{TRUE}, only RL estimates are computed.
#' If \code{"both"}, both parameter and RL estimates are computed.
#' @param varmeth The method used to estimate the covariance matrix of GEV parameters: 'V' or 'V2'.
#' @param conf.level A probability between 0 and 1; The confidence level for which confidence intervals are estimated.
#' @param chain logical and only relevant for non-stationary models: Whether the method for
#'  covariance estimation should be based on chain rule
#' @param ... Further arguments passed to `fit_gev_univ`. When `type = scale`, one must
#'  pass the logical argument `rel_trend` which specifies the parametrisation of the scale model.
#'
#' @return When \code{estimate_RL = TRUE}:
#'  A tibble containing for each combination of \code{Year}, \code{ref_gmst} and
#' \code{Kblocks} the following values:
#'  \describe{
#'   \item{rl}{The RL that was estimated on the concatenated Kblock-sample.
#'   Only differs from the RL that was estimated on the full sliding BM sample when
#'   \code{looplastblock = TRUE}.}
#'   \item{var_hat}{The estimated Variance of the RL estimation.}
#'   \item{qlow, qupp}{The \eqn{ 1-\alpha/2} and \eqn{\alpha/2}- quantiles of the bootstrapped studentized RLs,
#'   where \eqn{\alpha} is the chosen confidence level.}
#'   \item{lower, upper}{The lower and upper bound of the confidence interval.}
#'   \item{further parameters}{Denoted by  \code{mu, sigma0, gamma, alpha} for 'scale', and 'shift', and by
#'   \code{loc, scale, shape} for 'stationary', giving the estimated GEV parameters as estimated
#'   on the Kblock-sample.}
#'   }
#'
#'   When \code{estimate_RL = "both"}: A list containing the above tibble as list element
#'   named \code{res_rl} as well as a tibble with confidence interval bounds for each of the
#'   parameter vectors.
#' @export
#'
#' @examples
#' \dontrun{
#' blcksz <- 90
#' ny <- 50
#' ## with a shift in location
#' yy <- evd::rgpd(ny*blcksz, shape = 0.2) + 2*rep(1:ny/ny, each = blcksz)
#' ci_student_boot_sl(yy, blcksz = 90, temp.cov = rep(1:ny/ny, each = blcksz),
#' B = 100, type = "shift", ref_gmst = c(0.5, 0.95), Kblock = 4)
#'
#' ## observations are scaling w.r.t. GMST
#' yy <- evd::rgpd(ny*blcksz, shape = 0.2)*exp(0.5*rep(1:ny/ny, each = blcksz))
#' ci_student_boot_sl(yy, blcksz = 90, temp.cov = rep(1:ny/ny, each = blcksz), Kblock = 4,
#' B = 100, type = "scale", ref_gmst = c(0.5, 0.95), chain = TRUE, rel_trend = FALSE, estimate_RL = "both", rel_trend = FALSE)
#'
##### fix looplastblock = FALSE
ci_student_boot_sl <- function(x, blcksz,
                               Kblock = c(4,8,10),
                               B = 200, temp.cov = NULL,
                               type = "shift",
                               scale.link = make.link("identity"),
                               Tyrl = c(50, 100), ref_gmst = c(0.8, 0.9), varmeth = "V2",
                               looplastblock = TRUE, reltol = 1e-09,
                               estimate_RL = TRUE, conf.level = 0.05, chain = TRUE,
                               returnBootSamples = FALSE, ...) {

  add.args <- list(...)
  # lenght of time series (of daily observations given in x)
  ndata <- length(x)

  # the sample of disjoint block maxima

  djbm <- slbm::blockmax(x, r = blcksz, "disjoint")


  # compute the unique sliding BM and their frequency
  df.x <- data.frame(Station = "X1", Obs = x)

  sluniq_comp <- slbm::get_uniq_bm(x, blcksz = blcksz, temp_cvrt = temp.cov,
                                   looplastblock = FALSE)

  slbm_orig <- slbm::blockmax(x, r = blcksz, "sliding")

  temp_cvrt_orig <- temp.cov[1:length(slbm_orig)]
  # estimate on original block maxima sample
  est_sl <- tryCatch(slbm::fit_gev_univ(data = sluniq_comp, hessian = TRUE, type = type,
                                      #  return_cov = TRUE, varmeth = varmeth, chain = chain,
                                      rel_trend = add.args$rel_trend
                                      #  orig_slbm = slbm_orig, orig_cvrt = temp_cvrt_orig, blcksz = blcksz
                                      ),
                     error = errfct(type))

  rlest <- dplyr::tibble()

  paramest <- dplyr::tibble()

  est_sl_df <- data.frame(t(est_sl$mle))

  # pefrom block bootstrap for each block bootstrap parameter k
  for(k in Kblock) {
  # number of blocks of size K*blcksz that are observed
  .nKblocks <- ceiling(ndata/(k*blcksz))


  # for each observation: which Kblock does it belong to, and is it in the last
  # block of site blcksz within that larger Kblock
  indexblock <- data.frame(blockind = c(rep(1:(.nKblocks-1), each = k*blcksz),
                                        rep(.nKblocks, ndata - k*blcksz*(.nKblocks-1))),
                           obsind = 1:ndata)


  # compute the unique sliding BM within each of the
  # nKblocks of size Kblock, their frequency and their Kblockindex
  # (index of the bigger block of size K the sliding BM belongs to)
  sluniq_wb <- slbm::get_uniq_bm_boot(df.x, blcksz = blcksz, K = k, indexblock = indexblock,
                                temp_cvrt = temp.cov, looplastblock = looplastblock,
                                returnfullsamp = TRUE)

  full_slbm_Kblock <- sluniq_wb %>% dplyr::select(-uniq_data)
  sluniq_wb <- sluniq_wb %>% dplyr::select(-full_data)

 # estimate Variance of Kblock estimator
  slorig_Kblock <- purrr::map_dfr(full_slbm_Kblock$full_data, ~ .x)
  #  compute ML estimator
  estim_Kblock <- tryCatch(slbm::fit_gev_univ(tidyr::unnest(sluniq_wb, cols = uniq_data),
                                              type = type, hessian = TRUE, rel_trend = add.args$rel_trend,
                                              return_cov = TRUE,
                                              varmeth = varmeth, chain = chain,
                                              orig_slbm = slorig_Kblock$slbm, orig_cvrt = slorig_Kblock$temp_cvrt, blcksz = blcksz
                                              ), error = errfct(type))

  # bootstrap RLs and compute their parametric variance estimation
  bootres <- slbm::boot_sl(
    sluniq = sluniq_wb,
    slorig = full_slbm_Kblock,
    blcksz = blcksz,
    nKblocks = .nKblocks,
    B = B,
    type = type,
    ref_gmst = ref_gmst,
    Tyrl = Tyrl,
    start_vals = estim_Kblock$mle,
    scale.link = make.link("identity"),
    reltol = 1e-09, estimate_RL = estimate_RL, varmeth = varmeth,
    chain = chain, rel_trend = add.args$rel_trend)

    # put parameters into df
     df_estim_Kblock <- data.frame(t(estim_Kblock$mle))

      # compute the quantiles and confidence bounds from the bootstrap sample
     res_K <- get_quants_from_boot(bootres = bootres, estim_Kblock = estim_Kblock,
                                   ref_gmst = ref_gmst,
                                      Tyrl = Tyrl,
                                   estimate_RL = estimate_RL, varmeth = varmeth,
                                   type = type, conf.level = conf.level,
                                   rel_trend = add.args$rel_trend,
                                   var_crit = add.args$var_crit,
                                   var_crit_rl = add.args$var_crit_rl)

     if(estimate_RL == "both") {
       paramest <- paramest %>%
         dplyr::bind_rows(res_K$res_params %>% dplyr::mutate( Kblock = k) )

       rlest <- rlest %>%
         dplyr::bind_rows(res_K$res_rl %>% dplyr::mutate(df_estim_Kblock, Kblock = k))


     }
     else if(estimate_RL == FALSE) {
       paramest <- paramest %>%
         dplyr::bind_rows(res_K$res_params %>% dplyr::mutate( Kblock = k) )

     }
     else {
       rlest <- rlest %>%
         dplyr::bind_rows(res_K$res_rl %>% dplyr::mutate(df_estim_Kblock, Kblock = k) )
     }
   }



  if(estimate_RL %in% c(TRUE, "both")) {
    rl_hat_full <- slbm::compute_rl(theta = est_sl$mle, Tyrl = Tyrl, type = type, ref_gmst = ref_gmst,
                                    rel_trend = add.args$rel_trend)

    rl_hat_full <- dplyr::rename(rl_hat_full, "rlfull" = "rl")
    if(type == "stationary") {
      rlest <- rlest %>% dplyr::left_join(rl_hat_full, by = c("Year"))
    }
    else {
      rlest <- rlest %>% dplyr::left_join(rl_hat_full, by = c("refGMST", "Year"))
    }
  }

  if(!returnBootSamples) {
     if(isTRUE(estimate_RL)) { return( list(res_rl = rlest)) }
     if(estimate_RL == FALSE){ return( list(res_params = paramest)) }
     if(estimate_RL == "both"){ return( list(res_params  = paramest,
                                          res_rl  = rlest)) }
  }
  else {
    if(isTRUE(estimate_RL)) { return( list(res_rl = rlest, bootsamples = bootres)) }
    if(estimate_RL == FALSE){
      bootres <- bootres$resV2 %>% dplyr::mutate( purrr::map_dfr(CovestV, ~ { aa <- diag(.x)
      names(aa) <- paste0("var" , names(aa))
      as.data.frame(t(aa))

      }) ) %>% dplyr::select(- CovestV)
      bootres <- bootres$resV2 %>% tidyr::pivot_longer(cols  = 1:4, names_to = "parameter", values_to  = "boot_est")
      bootres <- bootres %>% mutate(varhat_boot = purrr::map2_dbl(CovestV, parameter, ~ { .x[.y,  .y]})) %>% select(-CovestV)


      return( list(res_params = paramest, bootsamples = bootres %>% left_join(paramest %>% select(parameter, par_est)))) }
    if(estimate_RL == "both"){ return( list(res_params  = paramest,
                                            res_rl  = rlest, bootsamples = bootres)) }
  }




}


compute_quants_boot <- function(bootres, Tyrl, estim_Kblock, rlhat_Kblock = NULL,
                                type, ref_gmst, Covmat, conf.level = 0.05, target = "RL", rel_trend = NULL,
                                var_crit = NULL, var_crit_rl = NULL) {

  if(target == "RL") {

    if(!is.null(var_crit_rl)) {
      bootres <- bootres %>% dplyr::filter(var_hat > var_crit_rl)
    }
    # compute variance of RL estimation on Kblock sample
    varrl <- purrr::map_dfr(Tyrl, ~ {
      dplyr::tibble(var_hat =
                      slbm::estimate_var_rl(theta = estim_Kblock$mle, Tyrl = .x, type = type,
                                            ref_gmst = ref_gmst, Covmat = Covmat, rel_trend = rel_trend),
                    refGMST = ref_gmst, Year = .x) })

    if(type == "stationary") {
      # compute quantiles of the bootstrapped studentized RLs
      tstar_quants <- bootres %>%
        dplyr::mutate(tstar = (rlboot - rl)/sqrt(var_hat)) %>%
        dplyr::group_by(Year) %>%
        dplyr::summarise(qupp = quantile(tstar, p = c(1-conf.level/2), na.rm = TRUE),
                         qlow = quantile(tstar, p = conf.level/2, na.rm = TRUE),
                         .groups = "drop")

      # join results and compute confidence bounds
      res <- varrl %>%
        dplyr::left_join(tstar_quants, by = c("Year")) %>%
        dplyr::left_join(rlhat_Kblock, by = c("Year")) %>%
        dplyr::mutate(lower = rl - sqrt(var_hat)*qupp, upper = rl - sqrt(var_hat)*qlow)
    }
    else {

      # compute quantiles of the bootstrapped studentized RLs
      tstar_quants <- bootres %>%
        dplyr::mutate(tstar = (rlboot - rl)/sqrt(var_hat)) %>%
        dplyr::group_by(refGMST, Year) %>%
        dplyr::summarise(qupp = quantile(tstar, p = c(1-conf.level/2), na.rm = TRUE),
                         qlow = quantile(tstar, p = conf.level/2, na.rm = TRUE),
                         .groups = "drop")

      # join results and compute confidence bounds
      res <- varrl %>%
        dplyr::left_join(tstar_quants, by = c("refGMST", "Year")) %>%
        dplyr::left_join(rlhat_Kblock, by = c("refGMST", "Year")) %>%
        dplyr::mutate(lower = rl - sqrt(var_hat)*qupp, upper = rl - sqrt(var_hat)*qlow)
    }
    return(res)
  }
  if(target == "Param") {
    # compute variance of RL estimation on Kblock sample
    par_vars <- diag(Covmat)
    parsK <- as.data.frame(t(estim_Kblock$mle))

    if(type == "stationary") {
      bootpars <-  bootres %>%  dplyr::select(loc, scale, shape, CovestV) %>% unique()

      bootpars <- bootpars %>% dplyr::mutate(
        purrr::map_dfr(CovestV, ~ { aa <- diag(.x)
        names(aa) <- paste0("var" , names(aa))
        as.data.frame(t(aa))}) ) %>%
        dplyr::select(- CovestV)

      funtstar <- function(vec, pars) {(vec[1:3] - pars)/ sqrt(vec[4:6])}
      tstar_pars <- apply(bootpars, 1, funtstar, pars = parsK)
      tstar_pars <- purrr::map_dfr(tstar_pars, ~.x)

      # tstar_pars <- bootpars %>% dplyr::mutate( star =
      #                                             purrr::pmap(list(loc, scale, shape, CovestV),
      #                                                         function(loc, scale, shape, CovestV) {
      #                                                           dplyr::tibble((loc - parsK["loc"])/sqrt(CovestV[1,1]),
      #                                                                         (scale - parsK["scale"])/sqrt(CovestV[2,2]),
      #                                                                         (shape - parsK["shape"])/sqrt(CovestV[3,3]))
      #                                                         }))
     n.par <-  3
    }
    else {
      bootpars <-  bootres %>%  dplyr::select(mu, sigma, shape, alpha, CovestV) %>% unique()

      bootpars <- bootpars %>% dplyr::mutate( purrr::map_dfr(CovestV, ~ { aa <- diag(.x)
        names(aa) <- paste0("var" , names(aa))
        as.data.frame(t(aa))

      }) ) %>% dplyr::select(- CovestV)

      if(!is.null(var_crit)) {
        bootpars <- bootpars %>% dplyr::filter(varshape > var_crit)
      }

      funtstar <- function(vec, pars) {(vec[1:4] - pars)/ sqrt(vec[5:8])}
      tstar_pars <- apply(bootpars, 1, funtstar, pars = parsK)
      tstar_pars <- purrr::map_dfr(tstar_pars, ~.x)

      # tstar_pars <- bootpars %>% dplyr::mutate( star =
      #                                                purrr::pmap(list(mu, sigma, shape, alpha, CovestV),
      #                                                            function(mu, sigma, shape, alpha, CovestV) {
      #                                                              dplyr::tibble((mu - parsK["mu"])/sqrt(CovestV[1,1]),
      #                                                                            (sigma - parsK["sigma"])/sqrt(CovestV[2,2]),
      #                                                                            (shape - parsK["shape"])/sqrt(CovestV[3,3]),
      #                                                                            (alpha - parsK["alpha"])/sqrt(CovestV[4,4]))
      #                                                            }))
      n.par <- 4
    }

    # compute quantiles of the bootstrapped studentized parameter estimates
    # tstar_pars_quants <- tstar_pars %>% dplyr::select(star) %>% tidyr::unnest(cols = star) %>%
    #   dplyr::summarise_at(1:n.par,  ~ quantile(.x, p = c( 1 - conf.level/2, conf.level/2), na.rm = TRUE)) %>%
    #   dplyr::mutate(quant = paste0("q", c( 1 - conf.level/2, conf.level/2)))

    tstar_pars_quants <- tstar_pars %>%
      dplyr::summarise_at(1:n.par,  ~ quantile(.x, p = c( 1 - conf.level/2, conf.level/2), na.rm = TRUE)) %>%
      dplyr::mutate(quant = paste0("q", c( 1 - conf.level/2, conf.level/2)))

    lower_bounds <- (parsK - sqrt(par_vars) * tstar_pars_quants[1, 1:n.par]) %>%
      tidyr::pivot_longer(cols = 1:n.par,
                          names_to = "parameter", values_to = "lower")

    upper_bounds <- (parsK - sqrt(par_vars) * tstar_pars_quants[2, 1:n.par]) %>%
      tidyr::pivot_longer(cols = 1:n.par,
                          names_to = "parameter", values_to = "upper")

    parsK <- tidyr::pivot_longer(parsK, cols = 1:n.par, names_to = "parameter", values_to = "par_est")

    res_pars <- parsK %>% dplyr::left_join(lower_bounds, by = "parameter") %>%
      dplyr::left_join(upper_bounds, by = "parameter") %>%
      dplyr::mutate(varest = par_vars)

    return(res_pars)
  }
}


get_quants_from_boot <- function(bootres, estim_Kblock, ref_gmst,
                                 Tyrl,  estimate_RL, varmeth, type, conf.level = 0.05,
                                 var_crit = NULL, var_crit_rl = NULL, ...) {

  add.args <- list(...)
  Covlist <- switch(varmeth, "both" = list(estim_Kblock$V, estim_Kblock$V2),
                    "V" = list(estim_Kblock$V),
                    "V2" = list(estim_Kblock$V2))

  if(estimate_RL %in% c(FALSE, "both")) {
    res_pars <- purrr::map2(.x = bootres, .y = Covlist,
                ~ compute_quants_boot(bootres = .x,  Tyrl = Tyrl, estim_Kblock = estim_Kblock,
                                      type = type, ref_gmst = ref_gmst, Covmat = .y,
                                      conf.level = conf.level, target = "Param", var_crit = var_crit))

    res_pars <- purrr::map2_dfr(res_pars, names(res_pars),  ~ .x %>% dplyr::mutate(varmeth = .y))
  }


  if(estimate_RL %in% c(TRUE, "both")) {

    # RL estimated on kblock sample
    rlhat_Kblock <- slbm::compute_rl(theta = estim_Kblock$mle,
                                     Tyrl = Tyrl, type = type, ref_gmst = ref_gmst, rel_trend = add.args$rel_trend)


    bootres <- purrr::map(bootres, ~ dplyr::rename(.x, "rlboot" = "rl"))
    if(type == "stationary") {
      bootres <- purrr::map(bootres, ~ .x %>% dplyr::left_join(rlhat_Kblock,
                                                               by = "Year"))
    }
    else {
      bootres <- purrr::map(bootres, ~ .x %>% dplyr::left_join(rlhat_Kblock,
                                                               by = c("refGMST", "Year")))
    }

    res <- purrr::map2(.x = bootres, .y = Covlist,
                  ~ compute_quants_boot(bootres = .x,  Tyrl = Tyrl, estim_Kblock = estim_Kblock,
                                        rlhat_Kblock = rlhat_Kblock,
                          type = type, ref_gmst = ref_gmst, Covmat = .y, conf.level = conf.level, target = "RL",
                          rel_trend = add.args$rel_trend, var_crit_rl = var_crit_rl))

    res <-  purrr::map2_dfr(res, names(res),  ~ .x %>% dplyr::mutate(varmeth = .y))
  }
  if(isTRUE(estimate_RL)) { return(list(res_rl = res)) }
  if(estimate_RL == FALSE){ return(list(res_params = res_pars)) }
  if(estimate_RL == "both"){ return( list(res_params  = res_pars,
                                          res_rl  = res)) }

}




ci_student_boot_samples <- function(x, blcksz,
                               Kblock = c(4,8,10),
                               B = 200, temp.cov = NULL,
                               type = "shift",
                               scale.link = make.link("identity"),
                               Tyrl = c(50, 100), ref_gmst = c(0.8, 0.9), varmeth = "V2",
                               looplastblock = TRUE, reltol = 1e-09,
                               estimate_RL = TRUE, conf.level = 0.05, chain = TRUE,
                               returnBootSamples = FALSE, ...) {

  add.args <- list(...)
  # lenght of time series (of daily observations given in x)
  ndata <- length(x)

  # the sample of disjoint block maxima

  djbm <- slbm::blockmax(x, r = blcksz, "disjoint")


  # compute the unique sliding BM and their frequency
  df.x <- data.frame(Station = "X1", Obs = x)

  sluniq_comp <- slbm::get_uniq_bm(x, blcksz = blcksz, temp_cvrt = temp.cov,
                                   looplastblock = FALSE)

  slbm_orig <- slbm::blockmax(x, r = blcksz, "sliding")

  temp_cvrt_orig <- temp.cov[1:length(slbm_orig)]
  # estimate on original block maxima sample
  est_sl <- tryCatch(slbm::fit_gev_univ(data = sluniq_comp, hessian = TRUE, type = type,
                                        #  return_cov = TRUE, varmeth = varmeth, chain = chain,
                                        rel_trend = add.args$rel_trend
                                        #  orig_slbm = slbm_orig, orig_cvrt = temp_cvrt_orig, blcksz = blcksz
  ),  error = errfct(type))

  rlest <- dplyr::tibble()

  paramest <- dplyr::tibble()

  est_sl_df <- data.frame(t(est_sl$mle))

  # pefrom block bootstrap for each block bootstrap parameter k
  for(k in Kblock) {
    # number of blocks of size K*blcksz that are observed
    .nKblocks <- ceiling(ndata/(k*blcksz))


    # for each observation: which Kblock does it belong to, and is it in the last
    # block of site blcksz within that larger Kblock
    indexblock <- data.frame(blockind = c(rep(1:(.nKblocks-1), each = k*blcksz),
                                          rep(.nKblocks, ndata - k*blcksz*(.nKblocks-1))),
                             obsind = 1:ndata)


    # compute the unique sliding BM within each of the
    # nKblocks of size Kblock, their frequency and their Kblockindex
    # (index of the bigger block of size K the sliding BM belongs to)
    sluniq_wb <- slbm::get_uniq_bm_boot(df.x, blcksz = blcksz, indexblock = indexblock, K = k,
                                        temp_cvrt = temp.cov, looplastblock = looplastblock,
                                        returnfullsamp = TRUE)

    full_slbm_Kblock <- sluniq_wb %>% dplyr::select(-uniq_data)
    sluniq_wb <- sluniq_wb %>% dplyr::select(-full_data)

    # estimate Variance of Kblock estimator
    slorig_Kblock <- purrr::map_dfr(full_slbm_Kblock$full_data, ~ .x)
    #  compute ML estimator
    estim_Kblock <- tryCatch(slbm::fit_gev_univ(tidyr::unnest(sluniq_wb, cols = uniq_data),
                                                type = type, hessian = TRUE, rel_trend = add.args$rel_trend,
                                                return_cov = TRUE,
                                                varmeth = varmeth, chain = chain,
                                                orig_slbm = slorig_Kblock$slbm, orig_cvrt = slorig_Kblock$temp_cvrt, blcksz = blcksz
    ), error = errfct(type))

    # bootstrap RLs and compute their parametric variance estimation
    bootres <- slbm::boot_sl(
      sluniq = sluniq_wb,
      slorig = full_slbm_Kblock,
      blcksz = blcksz,
      nKblocks = .nKblocks,
      B = B,
      type = type,
      ref_gmst = ref_gmst,
      Tyrl = Tyrl,
      start_vals = estim_Kblock$mle,
      scale.link = make.link("identity"),
      reltol = 1e-09, estimate_RL = estimate_RL, varmeth = varmeth,
      chain = chain, rel_trend = add.args$rel_trend)

    #browser()
    # put parameters into df
    df_estim_Kblock <- data.frame(parameter = names(estim_Kblock$mle),
                                  par_est = unname(estim_Kblock$mle),
                                  var_hat = diag(estim_Kblock$V2))






      bootres_rl <- bootres$resV2 %>% select(rl, refGMST, Year, var_hat) %>%
        rename("rlboot" = "rl", "var_hat_boot" = "var_hat")

      bootres <- bootres$resV2 %>% select(mu, sigma, shape, alpha, CovestV) %>%
        tidyr::pivot_longer(cols  = 1:4, names_to = "parameter", values_to  = "boot_est") %>%
        mutate(varhat_boot = purrr::map2_dbl(CovestV, parameter, ~ { .x[.y,  .y]})) %>% select(-CovestV)


      rlhat <- slbm::compute_rl(estim_Kblock$mle, Tyrl = Tyrl, type = type, ref_gmst = ref_gmst)

      varrl <- purrr::map_dfr(Tyrl, ~ {
        dplyr::tibble(var_hat =
                        slbm::estimate_var_rl(theta = estim_Kblock$mle, Tyrl = .x, type = type,
                                              ref_gmst = ref_gmst, Covmat = estim_Kblock$V2, rel_trend = rel_trend),
                      refGMST = ref_gmst, Year = .x) })

      bootres_rl <- bootres_rl %>% left_join(varrl,  by = c("refGMST", "Year")) %>%
        left_join(rlhat,  by = c("refGMST", "Year"))

      bootres <- bootres %>% left_join(df_estim_Kblock)


      rlest <- rlest %>% bind_rows(tibble::tibble(boot_samples = list(bootres_rl), Kblock = k))
      paramest <- paramest %>% bind_rows(tibble::tibble(boot_samples = list(bootres), Kblock = k))
  }



  return(list(res_rl = rlest, res_par = paramest))




}
