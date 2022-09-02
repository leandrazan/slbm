

data("GMST")
set.seed(1)
ny <- 50
blcksz <- 90
spatial_cvrt <- data.frame(lat  = seq(0, 8, length = 8),
  lon = runif(8))


iteration_spat_cov <- function(blcksz, n.years, spatial_cvrt,
                               loc0, loc1 = NULL, loc2 = NULL, scale0 , scale1 = NULL, scale2 = NULL,
                               shape, loc.sp.form, scale.sp.form) {


  coords <- cbind(spatial_cvrt$lon, spatial_cvrt$lat)
  mdax <- MDAdata(n.years*blcksz, locations = coords,
                  margins = list(distr = "gpd", shape = shape),
                  ms = list(kovmod = "whitmat", nugget = 0, range = 3, smooth = 0.44))

  #tempcov1 <- data.frame(GMST = rep(GMST$smoothedGMST[91:140], each = 90))
  xx <- sim_spatial(mdax, params = c("loc0" = loc0, "loc1" = loc1, "loc2" = loc2,
                                     "scale0" = scale0, "scale1" = scale1, "shape" = shape),
                                    # "tempLoc" = 2),
                    loc.sp.form = loc.sp.form, scale.sp.form = scale.sp.form,
                    spat.cov = spatial_cvrt )#, temp.cov = tempcov1,
                   # loc.temp.form = ~ GMST)

  origslbm <- apply(xx, 2, blockmax, r = blcksz, "sliding")
  origslbm <- reshape2::melt(origslbm)
  origslbm <- rename(origslbm, "Obs" = "value", "Station" = "Var2") %>% select(-Var1)

  yy <- data.frame(xx) %>%
       tidyr::pivot_longer( 1:ncol(xx), names_to = "Station", values_to = "Obs")
 ##############################
  bmuniq <- get_uniq_bm(yy, blcksz, looplastblock = FALSE)
  # disjoint block maxima are used for computing starting values:
  djbm <- apply(xx, 2, blockmax, r = 90, "disjoint")


  params <-  fit_spgev_sl(data = bmuniq, loc.sp.form = loc.sp.form,
                          scale.sp.form = scale.sp.form,
                          spat.cov = spatial_cvrt,
                          datastart = djbm, st_val_meth = "LeastSq",
                          return_hess = TRUE, print_start_vals = FALSE,
                          scale.link = make.link("identity"))

  covest <- compute_sl_cov_spat(origslbm, params = params$mle, hessmat = params$hessian,
                       blcksz = blcksz, varmeth = "V2", loc.sp.form = loc.sp.form,
                       scale.sp.form = scale.sp.form, spatial_cvrt = spatial_cvrt)


 as_tibble(t(params$mle)) %>% mutate(Cov = list(covest))
}


#
# Res <- purrr::map_dfr(1:100,
#                        ~ iteration_spat_cov(90, 50, spatial_cvrt = spatial_cvrt, loc0 = 10,
#                                             loc1= 2, loc2 = 1.2, scale0 = 2,
#                                             scale1 = 1.2,
#                    shape = 0.2, loc.sp.form = ~ lon + lat, scale.sp.form = ~ lon))
#
# Res
# cov(Res[ , 1:6])
#
# Res %>% mutate( purrr::map_dfr(Cov, ~ data.frame(
#   varloc0 = .x[1, 1], loc0loc1 = .x[1, 2], loc0loc2 = .x[1,3],
#   loc0scale0 = .x[1, 4], loc0scale1 = .x[1, 5], loc0shape = .x[1, 6],
#   varloc1 = .x[2, 2], loc1loc2 = .x[2,3], loc1scale0 = .x[2,4], loc1scale1 = .x[2,5],
#   loc1shape = .x[2,6],
#   varloc2 = .x[3,3], loc2scale0 = .x[3,4], loc2scale1 = .x[3,5], loc2shape = .x[3,6],
#   varscale0 = .x[4,4], scale0scale1 = .x[4,5], scale0shape = .x[4,6],
#   varscale1 = .x[5,5], scale1shape = .x[5,6],
#   varshape = .x[6,6] ))) %>%
#   summarise_at(8:28, ~mean(.x))

#
iteration_shift_cov <- function(blcksz, n.years, temp_cvrt,
                               mu, alpha, sigma,
                               shape, type = "shift", method = "neu", varmeth = "V") {

  ### simulate some data with a linear trend
  xx <- evd::rgpd(blcksz*n.years, shape = shape) + 2*rep(temp_cvrt, each = blcksz)
  ## temporal covariate for sliding BM
  temp_cvrt_dj <- temp_cvrt
  temp_cvrt <- rep(temp_cvrt, each = blcksz)[1:((n.years -1)*blcksz +1)]

  ## compute sample of uniuqe sliding BM
  bms <- get_uniq_bm(xx, blcksz, temp_cvrt = temp_cvrt, looplastblock = FALSE)

  ## full sample of sliding for estimaing the covariance matrix
  slbm <- blockmax(xx, r = blcksz, "sliding")
  estim <- fit_gev_univ(data = bms, type = "shift", hessian = TRUE)

  if(method == "neu") {
    covest <- est_var_univ_shift(slbm, est_par = estim, blcksz = blcksz, temp.cov = temp_cvrt,
                                 temp.cov.dj = temp_cvrt_dj,
      varmeth = varmeth)
  }
  else {
    covest <- est_var_univ(slbm, est_par = estim, blcksz = blcksz, temp.cov = temp_cvrt,
                                 type = "shift",
                                 varmeth = varmeth)
  }



  as_tibble(t(estim$mle)) %>% mutate(Cov = list(covest))
}

# set.seed(1)
# Res <- purrr::map_dfr(1:500, ~ iteration_shift_cov(90, 100, GMST[43:142, 3],
#                                                  mu = 30, alpha = 3, sigma = 1.5,
#                                                  shape = 0.2, method = "alt", varmeth = "V2"))
#
#
# methalt <- Res %>% mutate(
#   purrr::map_dfr(Cov, ~ data.frame(varmu = .x$V2[1,1], musigma = .x$V2[1,2], mushape = .x$V2[1,3],
#                                    mualpha = .x$V2[1,4],
#                                    varsigma = .x$V2[2,2], sigmashape = .x$V2[2,3],
#                                    sigmaalpha = .x$V2[2,4],
#                                    varshape  = .x$V2[3,3], shapealpha = .x$V2[3,4],
#                                    varalpha = .x$V2[4,4]))
# )%>% summarise_at(6:15, ~ median(.x, na.rm = TRUE))
#
# set.seed(1)
# Res <- purrr::map_dfr(1:500, ~ iteration_shift_cov(90, 100, GMST[43:142, 3],
#                                                    mu = 30, alpha = 3, sigma = 1.5,
#                                                    shape = 0.2, method = "neu", varmeth = "V2"))
#
#
# methneu <- Res %>% mutate(
#   purrr::map_dfr(Cov, ~ data.frame(varmu = .x$V2[1,1], musigma = .x$V2[1,2], mushape = .x$V2[1,3],
#                                    mualpha = .x$V2[1,4],
#                                    varsigma = .x$V2[2,2], sigmashape = .x$V2[2,3],
#                                    sigmaalpha = .x$V2[2,4],
#                                    varshape  = .x$V2[3,3], shapealpha = .x$V2[3,4],
#                                    varalpha = .x$V2[4,4]))
# )%>% summarise_at(6:15, ~ median(.x, na.rm = TRUE))
#
# A1 <- cov(Res[, 1:4])
#  methneu
#   A1
# methalt
#  goals <- A1[upper.tri(A1, diag = TRUE)]
#  methneu - goals
#  methalt - goals
# abs(methneu - goals) < abs(methalt - goals)
#  round(methneu - methalt, 10)
# A1
