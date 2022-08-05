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
#'
#' @return Returns a tibble with columns `Station`, `Kblockind`,`uniq_data` and, if `returnfullsamp  = TRUE`,
#' `ful_data`. These are the station reference, the block index of the blocks of size `K` times `blcksz`,
#' the sample of unique sliding block maxima values and their frequency and the looped sliding block maxima.
#' @export
#'
#' @examples
#' blcksz <- 90
#' ny <- 50
#' yy <- evd::rgpd(ny*blcksz, shape = 0.2)
#' df.yy <- data.frame(Station = "X1", Obs = yy)
#' k <- 3
#' ndata <- blcksz * ny
#' .nKblocks <- ceiling(ndata/(k*blcksz))
#' indexblock <- data.frame(blockind = c(rep(1:(.nKblocks-1), each = k*blcksz),
#'                                       rep(.nKblocks, ndata - k*blcksz*(.nKblocks-1))),
#'                          obsind = 1:ndata)
#'
#' sluniq_wb <- get_uniq_bm_boot(df.yy, blcksz = blcksz, indexblock = indexblock,
#'                               temp_cvrt = NULL, looplastblock = TRUE,
#'                               returnfullsamp = TRUE)
#'
# bei looplastblock = TRUE mit Kovariable
get_uniq_bm_boot <- function (data, blcksz, indexblock, temp_cvrt = NULL,
                              looplastblock = TRUE, returnfullsamp = FALSE) {

  data$Kblockind <- indexblock$blockind
  if (!is.null(temp_cvrt)) {
    n.cvrt <- length(temp_cvrt)
    data$tempcvrt <- c(temp_cvrt, rep(temp_cvrt[n.cvrt], (nrow(data) - n.cvrt)))

    if(returnfullsamp) {

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

#' Compute RL of fitted GEV distribution
#'
#' @param theta Estimated GEV parameters
#' @param Tyrl The return periods for which the RLs are computed
#' @param type The model: 'stationary'. 'shift' or 'scale'
#' @param ref_gmst The reference value of the temporal covariate for which the
#' RL is computed. Ignored when \code{type = "stationary"}.
#'
#' @return A tibble containg the estimated RL for each combination of \code{Tyrl} and
#' \code{ref_gmst}.
#' @export
#'
#' @examples
#' compute_rl(theta = c(loc0 = 7, scale0 = 2, shape = 0.2, TempLoc1 = 2),
#' Tyrl = c(50, 100), type = "shift", ref_gmst = c(0.5, 0.95))
#'
compute_rl <- function(theta, Tyrl, type, ref_gmst = NULL) {

  if(type == "stationary") {
    rl <- evd::qgev(1-1/Tyrl, loc = theta[1], scale = theta[2], shape = theta[3])
    return(dplyr::tibble(rl = rl, Year = Tyrl ))
  }
  if(type == "scale") {
    mut <- theta[1]*exp(theta[4]/theta[1]*ref_gmst)
    sigmat <- theta[2]*exp(theta[4]/theta[1]*ref_gmst)
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
#'
#' @return A list containing
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
                           varmeth = "V2") {
  set.seed(seed)
  # the nKblocks blocks of size Kblock that make up the bootstrap sample
  chosenblocks <- sample(1:nKblocks, nKblocks, replace = TRUE)

  # the bootstrap sample. sluniq consists of columns containing the Kblockindex,
  # the unique values of sliding block maxima and their resp. frequency.
  # observations with Kblockind in chosenblocks are filtered
  boot_samp <- sluniq[chosenblocks, ] %>% tidyr::unnest(cols = uniq_data)

  slorig <- purrr::map_dfr(slorig[chosenblocks, ]$full_data, ~ .x)
  # slorig <- unlist((slorig  %>% dplyr::filter(Kblockind %in% chosenblocks))$full_data)


  # the ML estimator based on weighted log-Likelihood computed on the bootstrap sample
  est_boot <-  slbm::fit_gev_univ(data = boot_samp, hessian = TRUE, type  = type )


  estvar <- slbm::est_var_univ(orig_slbm = slorig$slbm, est_par = est_boot, blcksz = blcksz,
                               temp.cov = slorig$temp_cvrt, type = type,  varmeth = varmeth)

  if(!estimate_RL) {
    return(list(est_boot = est_boot, estvarV = estvar$V, estvarV2 = estvar$V2))
  }
  else {

    rlhat <- slbm::compute_rl(theta = est_boot$mle, Tyrl = Tyrl, type = type, ref_gmst = ref_gmst)
    if(varmeth %in% c("both", "V2")) {
      varrl_V2 <- purrr::map_dfr(Tyrl, ~ {
        dplyr::tibble(var_hat =
                 slbm::estimate_var_rl(theta = est_boot$mle, Tyrl = .x, type = type,
                                       ref_gmst = ref_gmst, Covmat = estvar$V2),
               refGMST = ref_gmst, Year = .x) })
    }
    if(varmeth %in% c("both", "V")) {
      varrl_V <- purrr::map_dfr(Tyrl, ~ {
        dplyr::tibble(var_hat =
                 slbm::estimate_var_rl(theta = est_boot$mle, Tyrl = .x, type = type,
                                       ref_gmst = ref_gmst, Covmat = estvar$V),
               refGMST = ref_gmst, Year = .x) })
    }
    if(varmeth == "V") {
      return(list(rl_boot = rlhat,  rlvarV = varrl_V))
    }
    if(varmeth == "V2") { return(list(rl_boot = rlhat, rlvarV2 = varrl_V2)) }
    else {return(list(rl_boot = rlhat, rlvarV = varrl_V, rlvarV2 = varrl_V2)) }

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
#' @param estimate_RL Must be TRUE. Not implemented atm (when interested in GEV
#' parameters rather than RLs)
#' @param Tyrl The return periods for which RLs are estimated
#' @param ref_gmst Vector of reference values of the temporal covariate for which
#' RLs are estimated. Ignored when type is 'stationary'.
#' @param varmeth The method used to estimate the covariance matrix of GEV parameters: 'V' or 'V2'.
#'
#' @return
#' @export
#'
#' @examples
#' blcksz <- 90
#' ny <- 100
#' yy <-  evd::rgpd(ny*blcksz, shape = 0.2)*2*rep(1:ny/ny, each = blcksz)
#' df.yy <- data.frame(Station = "X1", Obs = yy)
#' k <- 3
#' ndata <- blcksz * ny
#' .nKblocks <- ceiling(ndata/(k*blcksz))
#' indexblock <- data.frame(blockind = c(rep(1:(.nKblocks-1), each = k*blcksz),
#'                                       rep(.nKblocks, ndata - k*blcksz*(.nKblocks-1))),
#'                          obsind = 1:ndata)
#'
#' sluniq_wb <- get_uniq_bm_boot(df.yy, blcksz = blcksz, indexblock = indexblock,
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
#' boot_sl(
#' sluniq = sluniq_wb,
#' slorig = full_slbm,
#' blcksz = blcksz,
#' nKblocks = .nKblocks,
#' B = 10,
#' type = "shift",
#' ref_gmst = (c(0.8, 0.9)),
#' Tyrl = c(50, 100),
#' start_vals = estim_lo$mle,
#' scale.link = make.link("identity"),
#' reltol = 1e-09, estimate_RL = TRUE)
boot_sl <- function(sluniq, slorig, nKblocks,
                    blcksz, B = 50,
                    start_vals = NULL,
                    type = "stationary",
                    scale.link = scale.link, reltol = 1e-08,
                    estimate_RL = TRUE,
                    Tyrl = c(50, 100), ref_gmst = 0.8, varmeth = "V2") {


  # perform bootstrap B times
  boot_estimators <- purrr::map(1:B,
                                ~ slbm::sample_boot_sl(nKblocks = nKblocks,
                                                 sluniq = sluniq,
                                                 slorig = slorig,
                                                 blcksz = blcksz,
                                                 start_vals = start_vals,
                                                 type = type,
                                                 scale.link = scale.link,
                                                 .x,
                                                 # reltol = 1e-08,
                                                 estimate_RL = estimate_RL,
                                                 Tyrl = Tyrl,
                                                 ref_gmst = ref_gmst,
                                                 varmeth = varmeth))

  if(!estimate_RL) {
    boot_parest <- purrr::map_dfr(boot_estimators, ~ tibble(
      data.frame(t(.x$est_boot$mle), conv = .x$est_boot$conv),
      CovestV = list(.x$estvarV), CovestV2 = list(.x$estvarV2)))
  }
  else {
    if(varmeth %in%  c("both", "V2")) {
      boot_parest_v2 <- purrr::map_dfr(boot_estimators, ~
                                         dplyr::left_join(.x$rl_boot, .x$rlvarV2, by = c("refGMST", "Year")))
    }
    if (varmeth %in%  c("both", "V")) {
      boot_parest_v <- purrr::map_dfr(boot_estimators, ~
                                        dplyr::left_join(.x$rl_boot, .x$rlvarV, by = c("refGMST", "Year")))
    }
    if(varmeth == "both") {
      boot_parest <- boot_parest_v2 %>% dplyr::mutate(varmeth = "V2") %>%
        dplyr::bind_rows(boot_parest_v %>% dplyr::mutate(varmeth = "V"))
    } else {
      boot_parest <- base::switch(varmeth, "V" = boot_parest_v,
                            "V2" = boot_parest_v2)
    }
  }



  boot_parest

}

errfct <- function(type) {
  if(type == "stationary") {
    list(mle = c(loc = NA, scale = NA, shape = NA),
         hessian = NA)
  } else if (type == "shift") {
    list(mle = c(loc0 = NA, scale0 = NA, shape = NA,
                 tempLoc1 = NA),
         hessian = NA)
  } else {
    list(mle = c(mu0 = NA, sigma0 = NA, gamma = NA,
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
#' @param yy Numeric vector of data with daily resolution
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
#' @param varmeth The method used to estimate the covariance matrix of GEV parameters: 'V' or 'V2'.
#'
#' @return A tibble containing for each combination of \code{Year}, \code{ref_gmst} and
#' \code{Kblocks} the following values:
#'  \describe{
#'   \item{rl}{The RL that was estimated on the concatenated Kblock-sample.
#'   Only differs from the RL that was estimated on the full sliding BM sample when
#'   \code{looplastblock = TRUE}.}
#'   \item{var_hat}{The estimated Variance of the RL estimation.}
#'   \item{q0975, q0025}{The 97.5%- and 2.5%- quantiles of the bootstrapped studentized RLs.}
#'   \item{lower, upper}{The lower and upper bound of the confidence interval.}
#'   \item{further parameters}{Denoted by  \code{mu0, sigma0, gamma, alpha} for 'scale',
#'   \code{loc0, scale0, shape, TempLoc1} for 'shift' and by
#'   \code{loc, scale, shape} for 'stationary', giving the estimated GEV parameters as estimated
#'   on the Kblock-sample.}
#'   }
#' @export
#'
#' @examples
#' blcksz <- 90
#' ny <- 50
#' ## with a shift in location
#' yy <- evd::rgpd(ny*blcksz, shape = 0.2) + 2*rep(1:ny/ny, each = blcksz)
#' ci_student_boot_sl(yy, blcksz = 90, temp.cov = rep(1:ny/ny, each = blcksz),
#' B = 100, type = "shift", ref_gmst = c(0.5, 0.95))
#'
#' ## observations are scaling w.r.t. GMST
#' yy <- evd::rgpd(ny*blcksz, shape = 0.2)*exp(2/7.3*rep(1:ny/ny, each = blcksz))
#' ci_student_boot_sl(yy, blcksz = 90, temp.cov = rep(1:ny/ny, each = blcksz),
#' B = 100, type = "scale", ref_gmst = c(0.5, 0.95))
#'
##### fix looplastblock = FALSE
ci_student_boot_sl <- function(yy, blcksz,
                               Kblock = c(4,8,10),
                               B = 200, temp.cov = NULL,
                               type = "shift",
                               scale.link = make.link("identity"),
                               Tyrl = c(50, 100), ref_gmst = c(0.8, 0.9), varmeth = "V2",
                               looplastblock = TRUE, reltol = 1e-09) {

  # lenght of time series (of daily observations given in yy)
  ndata <- length(yy)

  # the sample of disjoint block maxima

  djbm <- slbm::blockmax(yy, r = blcksz, "disjoint")


  # compute the unique sliding BM and their frequency
  df.yy <- data.frame(Station = "X1", Obs = yy)

  sluniq_comp <- slbm::get_uniq_bm(yy, blcksz = blcksz, temp_cvrt = temp.cov,
                                   looplastblock = FALSE)

  slbm_orig <- slbm::blockmax(yy, r = blcksz, "sliding")

  est_sl <- tryCatch(slbm::fit_gev_univ(data = sluniq_comp,
                                        hessian = TRUE,
                                        type = type), error = errfct(type))


  rl_hat_full <- slbm::compute_rl(theta = est_sl$mle, Tyrl = Tyrl, type = type, ref_gmst = ref_gmst)

  rl_hat_full <- dplyr::rename(rl_hat_full, "rlfull" = "rl")


  parest <- dplyr::tibble()
  est_sl_df <- data.frame(t(est_sl$mle))

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
    sluniq_wb <- slbm::get_uniq_bm_boot(df.yy, blcksz = blcksz, indexblock = indexblock,
                                  temp_cvrt = temp.cov, looplastblock = looplastblock,
                                  returnfullsamp = TRUE)

    full_slbm_Kblock <- sluniq_wb %>% dplyr::select(-uniq_data)
    sluniq_wb <- sluniq_wb %>% dplyr::select(-full_data)

    #  compute ML estimator
    estim_Kblock <- tryCatch(slbm::fit_gev_univ(tidyr::unnest(sluniq_wb, cols = uniq_data),
                                                type = type, hessian = TRUE), error = errfct(type))

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
      reltol = 1e-09, estimate_RL = TRUE)


    # estimate Variance of Kblock estimator
    slorig_Kblock <- purrr::map_dfr(full_slbm_Kblock$full_data, ~ .x)


    estvar_Kblock <- slbm::est_var_univ(orig_slbm = slorig_Kblock$slbm,
                                        est_par = estim_Kblock, blcksz = blcksz,
                                        temp.cov = slorig_Kblock$temp_cvrt,
                                        type = type,  varmeth = varmeth)

    # RL estimated on kblock sample
    rlhat_Kblock <- slbm::compute_rl(theta = estim_Kblock$mle,
                               Tyrl = Tyrl, type = type, ref_gmst = ref_gmst)
    bootres <- dplyr::rename(bootres, "rlboot" = "rl")
    bootres <- bootres %>% dplyr::left_join(rlhat_Kblock, by = c("refGMST", "Year"))

    if(varmeth  ==  "V2") {
      # compute variance of RL estimation on Kblock sample
      varrl_V2 <- purrr::map_dfr(Tyrl, ~ {
        dplyr::tibble(var_hat =
                        slbm::estimate_var_rl(theta = estim_Kblock$mle, Tyrl = .x, type = type,
                                              ref_gmst = ref_gmst, Covmat = estvar_Kblock$V2),
                      refGMST = ref_gmst, Year = .x) })

      # compute quantiles of the bootstrapped studentized RLs
      tstar_quants_V2 <- bootres %>%
        dplyr::mutate(tstar = (rlboot - rl)/sqrt(var_hat)) %>%
        dplyr::group_by(refGMST, Year) %>%
        dplyr::summarise(q0975 = quantile(tstar, p = c(0.975), na.rm = TRUE),
                         q0025 = quantile(tstar, p = 0.025, na.rm = TRUE),
                         .groups = "drop")

      # join results and compute confidence bounds
      res <- varrl_V2 %>%
        dplyr::left_join(tstar_quants_V2, by = c("refGMST", "Year")) %>%
        dplyr::left_join(rlhat_Kblock, by = c("refGMST", "Year")) %>%
        dplyr::mutate(lower = rl - sqrt(var_hat)*q0975, upper = rl - sqrt(var_hat)*q0025)

    }

    if(varmeth  ==  "V") {
      # compute variance of RL estimation on Kblock sample
      varrl_V <- purrr::map_dfr(Tyrl, ~ {
        dplyr::tibble(var_hat =
                        slbm::estimate_var_rl(theta = estim_Kblock$mle, Tyrl = .x, type = type,
                                              ref_gmst = ref_gmst, Covmat = estvar_Kblock$V),
                      refGMST = ref_gmst, Year = .x) })

      # compute quantiles of the bootstrapped studentized RLs
      tstar_quants_V <- bootres %>%
        dplyr::mutate(tstar = (rlboot - rl)/sqrt(var_hat)) %>%
        dplyr::group_by(refGMST, Year) %>%
        dplyr::summarise(q0975 = quantile(tstar, p = c(0.975), na.rm = TRUE),
                         q0025 = quantile(tstar, p = 0.025, na.rm = TRUE),
                         .groups = "drop")

      # join results and compute confidence bounds
      res <- varrl_V %>%
        dplyr::left_join(tstar_quants_V, by = c("refGMST", "Year")) %>%
        dplyr::left_join(rlhat_Kblock, by = c("refGMST", "Year")) %>%
        dplyr::mutate(lower = rl - sqrt(var_hat)*q0975, upper = rl - sqrt(var_hat)*q0025)

    }

    # put parameters into df
    estim_Kblock$mle <- data.frame(t(estim_Kblock$mle))

    parest <- parest %>%
      dplyr::bind_rows(res %>% dplyr::mutate(estim_Kblock$mle, Kblock = k) )


  }





  parest %>% dplyr::left_join(rl_hat_full, by = c("refGMST", "Year"))


}
