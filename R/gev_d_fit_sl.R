

nllh.sl.d <- function( parvec, bm_uni, sigma0link = make.link("log"), dur_offset = F ){

  mut <- parvec[["mut"]]
  sigma0 <-  parvec[["sigma0"]]  #  sigma0link$linkfun(parvec[["sigma0"]]) # parvec contains the value tilde(sigma0) = exp(sigma0) > 0
  xi <- parvec[["shape"]]
  if(!(dur_offset) ) {
    theta <- 0
  } else{
    theta <- parvec[["theta"]]
  }
  eta <- parvec[["eta"]]

  d <- bm_uni$duration

  dpt <- d + theta
  sigma.d <- sigma0/(dpt^eta)  # sigma0link$linkinv(sigma0)/(dpt^eta)
  y <- bm_uni$sldata/sigma.d - mut
  y <- 1 + xi * y

  if (dur_offset) {
    if ( eta <= 0 ||  theta < 0  || any(sigma.d <= 0) || any(y <= 0) || theta > 1 )
      return(10^6)
  } else {
    if ( eta <= 0  || any(sigma.d <= 0) || any(y <=
                                               0))
      return(10^6)
  }


  nllh.value <- 1/sum(bm_uni$n)*
    as.numeric(t( bm_uni$n) %*% ( log(sigma.d) + y^(-1/xi) + log(y) * (1/xi +
                                                                         1)))


  return( nllh.value )

}


#' Maximum Likelihood fitting of the d-GEV distribution to the sample of disjoint and
#' sliding block maxima
#'
#' @param agg_bm Output of get_agg_bm, i.e. a tibble containing the sample of disjoint
#'  and sliding block maxima of the ID process,
#' for different durations $d$
#' @param leave_out optional. Vector containing the years to exclude in the fit.
#' Useful for cross-validation.
#' @param dur_offset logical. Whether to fit a duration offset parameter theta
#'
#' @return list of 4.
#' * mle: dataframe of ML estimates based on disjoint and sliding BM sample
#' * nllh: values of the negative log-Likelihood evaluated at ML estimates
#' * counts: number of function evaluations when fitting based on sliding BM (passsed from optim)
#' * dur_offset: whether a duration offset parameter was fitted
#' @export
#'
#' @examples dates <- seq(as.POSIXct("2000-01-01 00:00:00"),
#' as.POSIXct("2001-12-31 23:00:00"),by = 'hour')
#' prec <- rgamma(length(dates))
#' example_data <- data.frame(MESS_DATUM = dates, RR = prec)
#'
#' agbm <- get_agg_bm(example_data, ds = c(1,2,4,8,16, 24, 48))
#' gev.d.fit.sl(agbm)
gev.d.fit.sl <- function(agg_bm, leave_out = 0 , dur_offset = F){

  # select disjoint bm to calculate starting values
  bm_dj <- agg_bm %>% dplyr::select(djbm) %>% tidyr::unnest(cols = djbm)

  bm_dj_train <- bm_dj %>% dplyr::filter(!(Year %in% leave_out))

  # delete dj bm from tibble for saving memory
  agg_bm_wy <- agg_bm %>% dplyr::select(-djbm)

  # compute initial parameters from dj bm
  optim_dj <- IDF::gev.d.fit(xdat  = bm_dj_train$djdata[!is.na(bm_dj_train$djdata)],
                             ds = bm_dj_train$duration[!is.na(bm_dj_train$djdata)],
                             sigma0link = make.link('log'),
                             theta_zero = !dur_offset  ,
                             show = FALSE)
  par_init <- optim_dj$mle

  if(!dur_offset){
    names( par_init) <- c("mut", "sigma0", "shape",  "eta")
  } else{
    names( par_init) <- c("mut", "sigma0", "shape", "theta", "eta")
  }

  par_init[["sigma0"]] <- exp(par_init[["sigma0"]])

  # get tibble with unique values and their frequency
  bm_unique <- agg_bm_wy   %>%
    dplyr::mutate( uniq_data =
                     purrr::map(.x = slbm,
                                .f = function(.x){
                                  .df <- .x %>% dplyr::filter(!(Year %in% leave_out)) %>%
                                    dplyr::group_by(sldata) %>%
                                    dplyr::summarise( n = n(), .groups = "drop")
                                }
                     )) %>%  dplyr::select(- slbm)

  bm_unique <- bm_unique %>% tidyr::unnest(cols = uniq_data)

  optim_sl <-  optim(par_init, fn = nllh.sl.d, bm_uni = bm_unique,
                     dur_offset = dur_offset )

  par_est <- data.frame( t(par_init) ) %>%
    dplyr::bind_rows(data.frame(t(optim_sl$par))) %>%
    dplyr::bind_cols(estimator = c("disjoint", "sliding"))

  return(list( mle = par_est, conv = optim_sl$convergence,
               nllh = data.frame(dj = optim_dj$nllh/length(bm_dj$djdata),
                                 sl = optim_sl$value),
               counts = data.frame(SL = optim_sl$counts),
               dur_offset = dur_offset))

}




nllh.sl.d.ms <- function( parvec, bm_uni, sigma0link = make.link("log"), mult_sc = TRUE,
                          dur_offset = F ){

  mut <- parvec[["mut"]]
  sigma0 <-  sigma0link$linkinv(parvec[["sigma0"]])  #  sigma0link$linkfun(parvec[["sigma0"]]) # parvec contains the value tilde(sigma0) = exp(sigma0) > 0
  xi <- parvec[["shape"]]

  if(!(dur_offset) ) {
    theta <- 0
  }  else{
    theta <- parvec[["theta"]]
  }
  eta <- parvec[["eta"]]
  if(mult_sc){
    eta2 <- parvec[["eta2"]]
  } else{
    eta2 <- 0
  }
  d <- bm_uni$duration

  dpt <- d + theta
  sigma.d <- sigma0/(dpt^(eta + eta2))  # sigma0link$linkinv(sigma0)/(dpt^eta)
  y <- bm_uni$sldata/sigma.d - mut*d^eta2
  y <- 1 + xi * y

  if (dur_offset) {
    if ( eta <= 0 ||  theta < 0  || any(sigma.d <= 0) || any(y <= 0) || theta > 1 )
      return(10^6)
  } else {
     if ( eta + eta2 <= 0  || any(sigma.d <= 0) || any(y <=
                                                      0))
      return(10^6)
  }


  nllh.value <- 1/sum(bm_uni$n)*
    as.numeric(t( bm_uni$n) %*% ( log(sigma.d) + y^(-1/xi) + log(y) * (1/xi +
                                                                         1)))


  return( nllh.value )

}



nllh.sl.d.ms.fl <- function( parvec, bm_uni, sigma0link = make.link("log"), mult_sc = TRUE,
                             dur_offset = F , int_offset = F){

  mut <- parvec[["mut"]]
  sigma0 <-  sigma0link$linkinv(parvec[["sigma0"]])  #  sigma0link$linkfun(parvec[["sigma0"]]) # parvec contains the value tilde(sigma0) = exp(sigma0) > 0
  xi <- parvec[["shape"]]

  if(!(dur_offset) ) {
    theta <- 0
  } else{
    theta <- parvec[["theta"]]
  }
  eta <- parvec[["eta"]]
  if(mult_sc){
    eta2 <- parvec[["eta2"]]
  }  else{
    eta2 <- 0
  }
  if(int_offset){
    tau <- parvec[["tau"]]
  } else{
    tau <- 0
  }

  d <- bm_uni$duration

  dpt <- d + theta
  sigma.d <- sigma0/(dpt^(eta + eta2)) + tau # sigma0link$linkinv(sigma0)/(dpt^eta)
  mu.d <- mut*(sigma0*(d + theta)^(-eta) + tau)
  y <- (bm_uni$sldata - mu.d)/sigma.d
  y <- 1 + xi * y

  if (dur_offset) {
    if ( eta <= 0 ||  theta < 0  || any(sigma.d <= 0) || any(y <= 0) || theta > 1 )
      return(10^6)
  } else {
    if ( eta + eta2 <= 0  || any(sigma.d <= 0) || any(y <=  0)  ||  tau < 0 ){
      return(10^6)
    }
  }


  nllh.value <- 1/sum(bm_uni$n)*
    as.numeric(t( bm_uni$n) %*% ( log(sigma.d) + y^(-1/xi) + log(y) * (1/xi +
                                                                         1)))


  return( nllh.value )

}


nllh.sl.d.ms.fl.curv <- function( parvec, bm_uni, sigma0link = make.link("log"), mult_sc = TRUE,
                             dur_offset = FALSE , int_offset = FALSE){

  mut <- parvec[["mut"]]
  sigma0 <-  sigma0link$linkinv(parvec[["sigma0"]])  #  sigma0link$linkfun(parvec[["sigma0"]]) # parvec contains the value tilde(sigma0) = exp(sigma0) > 0
  xi <- parvec[["shape"]]

  if(!(dur_offset) ) {
    theta <- 0
  } else{
    theta <- parvec[["theta"]]
  }
  eta <- parvec[["eta"]]
  if(mult_sc){
    eta2 <- parvec[["eta2"]]
  }  else{
    eta2 <- 0
  }
  if(int_offset){
    tau <- parvec[["tau"]]
  } else{
    tau <- 0
  }

  d <- bm_uni$duration

  dpt <- d + theta
  sigma.d <- sigma0/(dpt^(eta + eta2)) + tau # sigma0link$linkinv(sigma0)/(dpt^eta)
  mu.d <- mut*(sigma0*(d + theta)^(-eta) + tau)
  y <- (bm_uni$sldata - mu.d)/sigma.d
  y <- 1 + xi * y

  if (dur_offset) {
    if ( eta <= 0 ||  theta < 0  || any(sigma.d <= 0) || any(y <= 0) || theta > 1 )
      return(10^6)
  } else {
    if ( eta + eta2 <= 0  || any(sigma.d <= 0) || any(y <=  0)  ||  tau < 0 ){
      return(10^6)
    }
  }


  nllh.value <- 1/sum(bm_uni$n)*
    as.numeric(t( bm_uni$n) %*% ( log(sigma.d) + y^(-1/xi) + log(y) * (1/xi +
                                                                         1)))


  return( nllh.value )

}



gev.d.fit.sl.ms <- function(agg_bm, leave_out = 0 , mult_sc = TRUE, dur_offset = FALSE,
                            sigma0link = make.link("log")){

  # select disjoint bm to calculate starting values
  bm_dj <- agg_bm %>% dplyr::select(djbm) %>% tidyr::unnest(cols = djbm)

  bm_dj_train <- bm_dj %>% dplyr::filter(!(Year %in% leave_out))

  # delete dj bm from tibble for saving memory
  agg_bm_wy <- agg_bm %>% dplyr::select(-djbm)

  # compute initial parameters from dj bm
  optim_dj <- IDF::gev.d.fit(xdat  = bm_dj_train$djdata[!is.na(bm_dj_train$djdata)],
                             ds = bm_dj_train$duration[!is.na(bm_dj_train$djdata)],
                             sigma0link = sigma0link,
                             theta_zero = !dur_offset  ,
                             eta2_zero = !mult_sc,
                             show = FALSE)
  par_init <- optim_dj$mle

  if(!dur_offset &  !mult_sc){
    names( par_init) <- c("mut", "sigma0", "shape",  "eta")
  }
  if(!dur_offset &  mult_sc){
    names( par_init) <- c("mut", "sigma0", "shape",  "eta", "eta2")
    par_init[["eta2"]] <- par_init[["eta2"]]  - par_init[["eta"]]
  }
  if(dur_offset & mult_sc) {
    names( par_init) <- c("mut", "sigma0", "shape", "theta", "eta", "eta2")
    par_init[["eta2"]] <- par_init[["eta2"]]  - par_init[["eta"]]
  }




  # get tibble with unique values and their frequency

  inc_weight <- bm_dj_train %>% filter( Year %in% (leave_out -1) )
  inc_weight <- inc_weight %>% group_by(duration) %>% nest(col_inc = c(Year, djdata)) %>% ungroup()

  bm_unique <- agg_bm_wy   %>%
    dplyr::mutate( uniq_data =
                     purrr::map(.x = slbm,
                                .f = function(.x){
                                  .df <- .x %>% dplyr::filter(!(Year %in% leave_out)) %>%
                                    filter( !((Year %in% (leave_out - 1) & Index > 1)) ) %>%
                                    dplyr::group_by(sldata) %>%
                                    dplyr::summarise( n = n(), .groups = "drop")
                                }
                     )) %>%  dplyr::select(- slbm)
  bm_unique <- bm_unique %>% left_join(inc_weight, by = "duration")
  bm_unique <- bm_unique %>% mutate( uniq_data =  pmap(list(blcksz = blcksz, uniq_data = uniq_data, col_inc =  col_inc),
                             .f = function(blcksz, uniq_data, col_inc){
                               uniq_data$n[uniq_data$sldata %in% col_inc$djdata] <- uniq_data$n[uniq_data$sldata %in% col_inc$djdata] +
                                 blcksz
                               return(uniq_data)
                             }))
  bm_unique <- bm_unique %>% select(-col_inc) %>% tidyr::unnest(cols = uniq_data)

  optim_sl <-  optim(par_init, fn = nllh.sl.d.ms, bm_uni = bm_unique,
                     mult_sc = mult_sc,
                     dur_offset = dur_offset, sigma0link= sigma0link )

  par_init[["sigma0"]] <- sigma0link$linkinv(par_init[["sigma0"]])
  optim_sl$par[["sigma0"]] <- sigma0link$linkinv(optim_sl$par[["sigma0"]])

  par_est <- data.frame( t(par_init) ) %>%
    dplyr::bind_rows(data.frame(t(optim_sl$par))) %>%
    dplyr::bind_cols(estimator = c("disjoint", "sliding"))

  return(list( mle = par_est, conv = optim_sl$convergence,
               nllh = data.frame(dj = optim_dj$nllh/length(bm_dj$djdata),
                                 sl = optim_sl$value),
               counts = data.frame(SL = optim_sl$counts),
               dur_offset = dur_offset,
               mult_sc = mult_sc))

}



gev.d.fit.sl.slbmnew <- function(agg_bm, leave_out = 0 , mult_sc = TRUE, dur_offset = FALSE,
                                 sigma0link = make.link("log")){

 # bm_dj_train <- bm_dj %>% dplyr::filter(!(Year %in% leave_out))
  bm_dj_train <- agg_bm %>% dplyr::select(djbm) %>% tidyr::unnest(cols = djbm) %>%
    dplyr::filter(!(Year %in% leave_out))
  bm_dj <- bm_dj_train
  # delete dj bm from tibble for saving memory
  agg_bm_wy <- agg_bm %>% dplyr::select(-djbm)

  # compute initial parameters from dj bm
  optim_dj <- IDF::gev.d.fit(xdat  = bm_dj_train$djdata[!is.na(bm_dj_train$djdata)],
                             ds = bm_dj_train$duration[!is.na(bm_dj_train$djdata)],
                             sigma0link = sigma0link,
                             theta_zero = !dur_offset  ,
                             eta2_zero = !mult_sc,
                             show = FALSE)
  par_init <- optim_dj$mle

  if(!dur_offset &  !mult_sc){
    names( par_init) <- c("mut", "sigma0", "shape",  "eta")
  }
  if(!dur_offset &  mult_sc){
    names( par_init) <- c("mut", "sigma0", "shape",  "eta", "eta2")
    par_init[["eta2"]] <- par_init[["eta2"]]  - par_init[["eta"]]
  }
  if(dur_offset & mult_sc) {
    names( par_init) <- c("mut", "sigma0", "shape", "theta", "eta", "eta2")
    par_init[["eta2"]] <- par_init[["eta2"]]  - par_init[["eta"]]
  }




  # get tibble with unique values and their frequency

  bm_unique <- agg_bm_wy   %>%
    dplyr::mutate( uniq_data =
                     purrr::map(.x = slbm,
                                .f = function(.x){
                                  .df <- .x %>%
                                    dplyr::group_by(sldata) %>%
                                    dplyr::summarise( n = n(), .groups = "drop")
                                }
                     )) %>%  dplyr::select(- slbm)


  bm_unique <- bm_unique %>% tidyr::unnest(cols = uniq_data)

  optim_sl <-  optim(par_init, fn = nllh.sl.d.ms, bm_uni = bm_unique,
                     mult_sc = mult_sc,
                     dur_offset = dur_offset, sigma0link = sigma0link )

  par_init[["sigma0"]] <- sigma0link$linkinv(par_init[["sigma0"]])
  optim_sl$par[["sigma0"]] <- sigma0link$linkinv(optim_sl$par[["sigma0"]])

  par_est <- data.frame( t(par_init) ) %>%
    dplyr::bind_rows(data.frame(t(optim_sl$par))) %>%
    dplyr::bind_cols(estimator = c("disjoint", "sliding"))

  return(list( mle = par_est, conv = optim_sl$convergence,
               nllh = data.frame(dj = optim_dj$nllh/length(bm_dj$djdata),
                                 sl = optim_sl$value),
               counts = data.frame(SL = optim_sl$counts),
               dur_offset = dur_offset,
               mult_sc = mult_sc))

}

gev.d.fit.sl.msfl <- function(agg_bm, leave_out = 0 , mult_sc = TRUE, dur_offset = FALSE,
                              int_offset = FALSE, sigma0link = make.link("log")){

  # bm_dj_train <- bm_dj %>% dplyr::filter(!(Year %in% leave_out))
  bm_dj_train <- agg_bm %>% dplyr::select(djbm) %>% tidyr::unnest(cols = djbm) %>%
    dplyr::filter(!(Year %in% leave_out))
  bm_dj <- bm_dj_train
  # delete dj bm from tibble for saving memory
  agg_bm_wy <- agg_bm %>% dplyr::select(-djbm)

  # compute initial parameters from dj bm
  optim_dj <- IDF::gev.d.fit(xdat  = bm_dj_train$djdata[!is.na(bm_dj_train$djdata)],
                             ds = bm_dj_train$duration[!is.na(bm_dj_train$djdata)],
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
                                    dplyr::summarise( n = n(), .groups = "drop")
                                }
                     )) %>%  dplyr::select(- slbm)


  bm_unique <- bm_unique %>% tidyr::unnest(cols = uniq_data)

  optim_sl <-  optim(par_init, fn = nllh.sl.d.ms.fl, bm_uni = bm_unique,
                     mult_sc = mult_sc,
                     dur_offset = dur_offset, int_offset = int_offset, sigma0link = sigma0link )

  par_init[["sigma0"]] <- sigma0link$linkinv(par_init[["sigma0"]])
  optim_sl$par[["sigma0"]] <- sigma0link$linkinv(optim_sl$par[["sigma0"]])

  par_est <- data.frame( t(par_init) ) %>%
    dplyr::bind_rows(data.frame(t(optim_sl$par))) %>%
    dplyr::bind_cols(estimator = c("disjoint", "sliding"))

  return(list( mle = par_est, conv = optim_sl$convergence,
               nllh = data.frame(dj = optim_dj$nllh/length(bm_dj$djdata),
                                 sl = optim_sl$value),
               counts = data.frame(SL = optim_sl$counts),
               dur_offset = dur_offset,
               mult_sc = mult_sc))

}

gev.d.fit.sl.msfl.curv <- function(agg_bm, leave_out = 0 , mult_sc = TRUE, dur_offset = FALSE,
                              int_offset = FALSE, sigma0link = make.link("log"), Maxit = 1000){

  # bm_dj_train <- bm_dj %>% dplyr::filter(!(Year %in% leave_out))
  bm_dj_train <- agg_bm %>% dplyr::select(djbm) %>% tidyr::unnest(cols = djbm) %>%
    dplyr::filter(!(Year %in% leave_out))
  bm_dj <- bm_dj_train
  # delete dj bm from tibble for saving memory
  agg_bm_wy <- agg_bm %>% dplyr::select(-djbm)

  # compute initial parameters from dj bm
  optim_dj <- IDF::gev.d.fit(xdat  = bm_dj_train$djdata[!is.na(bm_dj_train$djdata)],
                             ds = bm_dj_train$duration[!is.na(bm_dj_train$djdata)],
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
                                    dplyr::summarise( n = n(), .groups = "drop")
                                }
                     )) %>%  dplyr::select(- slbm)


  bm_unique <- bm_unique %>% tidyr::unnest(cols = uniq_data)

  optim_sl <-  optim(par_init, fn = nllh.sl.d.ms.fl.curv, bm_uni = bm_unique,
                     mult_sc = mult_sc,
                     dur_offset = dur_offset, int_offset = int_offset, sigma0link = sigma0link,
                     control = list(maxit = Maxit))

  par_init[["sigma0"]] <- sigma0link$linkinv(par_init[["sigma0"]])
  optim_sl$par[["sigma0"]] <- sigma0link$linkinv(optim_sl$par[["sigma0"]])

  par_est <- data.frame( t(par_init) ) %>%
    dplyr::bind_rows(data.frame(t(optim_sl$par))) %>%
    dplyr::bind_cols(estimator = c("disjoint", "sliding"))

  return(list( mle = par_est, conv = optim_sl$convergence,
               nllh = data.frame(dj = optim_dj$nllh/length(bm_dj$djdata),
                                 sl = optim_sl$value),
               counts = data.frame(SL = optim_sl$counts),
               dur_offset = dur_offset,
               mult_sc = mult_sc))

}
