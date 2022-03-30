## Likelihood of spatial block maxima with covariates
#library(tidyverse)
#rm(list = ls())


### make one version that operates on sliding block maxima and can handle
## various temporal covariates
### the other version optimises the weighted log Likelihood, only one temporal covariate
## possible
### also implement under IF assumption


## uses exp() as link function for scale parameter
nll_spat <- function(params, loc.sp.form, scale.sp.form,
                     loc.temp.form = NULL, scale.temp.form = NULL,
                     data, spat.cov, temp.cov = NULL, hom = FALSE ){


  n.dat <- nrow(data)
  n.site <- ncol(data)
  n.spat <- nrow(spat.cov)
  n.temp <- ncol(temp.cov)


  params.sp <- list(loc = params[startsWith(names(params),"loc")],
                    scale = params[startsWith(names(params),"scale")],
                    shape = params[startsWith(names(params),"shape")])

  params.temp <- list(loc = params[startsWith(names(params),"tempLoc")],
                      scale = params[startsWith(names(params),"tempScale")])

  params.temp <- purrr::map(params.temp, function(.x) ifelse(length(.x) == 0, 0, .x) )

  model.loc.sp <- model.matrix(loc.sp.form, model.frame(loc.sp.form,
                                                        data = spat.cov, na.action = na.pass))

  model.scale.sp <- model.matrix(scale.sp.form, model.frame(scale.sp.form,
                                                            data = spat.cov, na.action = na.pass))




  if(!is.null(loc.temp.form)) {
    loc.temp.form <- update(loc.temp.form,  ~ . + 0)

    model.loc.temp <- model.matrix(loc.temp.form, model.frame(loc.temp.form,
                                                              data = data.frame(temp.cov),
                                                              na.action = na.pass))
  } else {
    model.loc.temp <- matrix(rep(0, n.dat), ncol = 1)
  }

  if(!is.null(scale.temp.form)) {
    scale.temp.form <- update(scale.temp.form,  ~ . + 0)
    scale.temp.form <- update(scale.temp.form,  ~ . + 0)
    model.scale.temp <- model.matrix(scale.temp.form,
                                     model.frame(scale.temp.form,
                                     data =  data.frame(temp.cov), na.action = na.pass))
  } else {
    model.scale.temp <- matrix(rep(0, n.dat), ncol = 1)
  }

  param.loc <- model.loc.sp %*% params.sp$loc
  param.loc <- matrix(rep(param.loc, n.dat), ncol = n.site, byrow = TRUE) +
    matrix(rep(model.loc.temp %*% params.temp$loc, n.site), ncol = n.site, byrow =F)

  param.scale <-  model.scale.sp %*% params.sp$scale

  param.scale <- matrix(rep( param.scale, n.dat), ncol = n.site, byrow = TRUE) +
    matrix(rep(model.scale.temp %*% params.temp$scale, n.site), ncol = n.site, byrow = F) # exp(model.scale.sp %*% paramH0$scale)

  # param.shape <- matrix(rep(params.sp$shape, n.site*n.obs), ncol = n.site)



  #sigma.mat <- matrix(rep( model.scale.sp %*% param$scale, n.dat),
  #                   ncol = d.dat, byrow = TRUE)

  #loc.mat <- matrix(rep(model.loc.sp %*% param$loc, n.dat),
  #  ncol = d.dat, byrow = TRUE)

  if(abs(params.sp$shape) < 1e-8){
    yy <- data/exp(param.scale) - param.loc/exp(param.scale)

    likh <-  sum( param.scale  + yy + exp(-yy), na.rm = TRUE)
  }
  else{
    yy <- 1+ params.sp$shape*(data/exp(param.scale) - param.loc/exp(param.scale))

    likh <- sum( param.scale + (1/params.sp$shape +1)*log(yy) + yy^(-1/params.sp$shape),
                 na.rm = TRUE)
    if(any(yy <0, na.rm = TRUE)){ likh <- 1e+10}

  }

  return(likh)
}


#yyuniq <- get_uniq_bm(yy, blcksz = 90)

#bm_unique <- bm_unique %>% tidyr::unnest(cols = uniq_data)

# data consists of unique values of sliding bm and resp. frequency of occurence
# implement without temporal covariates first
# later: unique block maxima grouped by station and covariate value
nll_spat_sl <- function(params, loc.sp.form, scale.sp.form,
         loc.temp.form = NULL, scale.temp.form = NULL,
         data, spat.cov, temp.cov = NULL, hom = FALSE, scale.link = make.link("log")){


  n.site <- nrow(data)
  n.spat <- nrow(spat.cov)
  # n.temp <- ncol(temp.cov)


  params.sp <- list(loc = params[startsWith(names(params),"loc")],
                    scale = params[startsWith(names(params),"scale")],
                    shape = params[startsWith(names(params),"shape")])

  # params.temp <- list(loc = params[startsWith(names(params),"tempLoc")],
  #                     scale = params[startsWith(names(params),"tempScale")])
  #
  # params.temp <- purrr::map(params.temp, function(.x) ifelse(length(.x) == 0, 0, .x) )

  model.loc.sp <- model.matrix(loc.sp.form, model.frame(loc.sp.form,
                                                        data = spat.cov, na.action = na.pass))

  model.scale.sp <- model.matrix(scale.sp.form, model.frame(scale.sp.form,
                                                            data = spat.cov, na.action = na.pass))

  param.loc <- data.frame(loc = model.loc.sp %*% params.sp$loc)



  param.scale <-  data.frame(scale = model.scale.sp %*% params.sp$scale)

  param_new <- dplyr::bind_cols(param.loc, param.scale) %>% dplyr::mutate( shape = params.sp$shape)


  data <- data %>% dplyr::bind_cols(param_new)

  if(abs(params.sp$shape) < 1e-8){
    data <- data %>%
      dplyr::mutate( likh = purrr::pmap_dbl( list( .ud = uniq_data, .loc = loc, .scale = scale, .shape = shape),
                                             .f = function(.ud, .loc, .scale, .shape){
                                               yy <- (.ud$slbm - .loc)/exp(.scale)
                                               t(.ud$n) %*% (.scale + yy + exp(-yy))
                                             }))

  }
  else{
    data <- data %>%
      dplyr::mutate( likh = purrr::pmap_dbl( list( .ud = uniq_data, .loc =loc, .scale = scale, .shape = shape),
                                             .f = function(.ud, .loc, .scale, .shape){
                                               yy <- 1 + .shape*(.ud$slbm - .loc)/exp(.scale)
                                               .ud$n %*% (.scale + (1/.shape +1)* log(yy) + yy^(-1/.shape) )
                                             }))

    # yy <- 1+ params.sp$shape*(data/exp(param.scale) - param.loc/exp(param.scale))
    #
    # likh <- sum( param.scale + (1/params.sp$shape +1)*log(yy) + yy^(-1/params.sp$shape),
    #              na.rm = TRUE)
    if(any(yy <0, na.rm = TRUE)){ likh <- 1e+10}

  }
  # data$likh
  sum(data$likh)
  #return(likh)
}


log_dens <- function(obs, weight,  loc_temp, loc_sp, scale_temp, scale_sp, shape) {

  if(is.null(loc_temp)) { loc_temp <- 0}
  if(is.null(scale_temp)) { scale_temp <- 0}
  if(is.null(loc_sp)) { loc_sp <- 0}
  if(is.null(scale_sp)) { scale_sp <- 0}

  purrr::pmap_dbl( list(.obs = obs, .weight = weight, .loc_temp = loc_temp, .scale_temp = scale_temp),
               .f = function(.obs, .weight, .loc_temp, .scale_temp) {

                 -.weight*evd::dgev(x = .obs, loc = .loc_temp + loc_sp,
                                    scale = exp(.scale_temp + scale_sp),
                                   shape = shape, log = TRUE)

               })
  # -weight*evd::dgev(x = obs, loc = loc_temp + loc_sp, scale = scale_temp + scale_sp,
  #           shape = shape, log = TRUE)
}





# use link function for scale parameter: to be completed
nll_spat_temp_sl_prep <- function(params, loc.sp.form, scale.sp.form,
                             loc.temp.form = NULL, scale.temp.form = NULL,
                             dataprep, spat.cov, scale.link = make.link("log")){


  params.sp <- list(loc = params[startsWith(names(params),"loc")],
                    scale = params[startsWith(names(params),"scale")],
                    shape = params[startsWith(names(params),"shape")])

  params.temp <- list(loc = params[startsWith(names(params),"tempLoc")],
                      scale = params[startsWith(names(params),"tempScale")])

  params.temp <- purrr::map(params.temp, function(.x) ifelse(length(.x) == 0, 0, .x) )

  model.loc.sp <- dataprep$MatLocSp

  model.scale.sp <- dataprep$MatScaleSp

  loc_sp <- data.frame(loc_sp = model.loc.sp %*% params.sp$loc)

  scale_sp <-  data.frame(scale_sp = model.scale.sp %*% params.sp$scale)

  dataprep <- dataprep$data

  dataprep <- dataprep %>% dplyr::bind_cols( loc_sp , scale_sp )



  if(!is.null(loc.temp.form)) {
    loc.temp.form <- update(loc.temp.form,  ~ . + 0)
    dataprep <- dataprep %>%
      dplyr::mutate(LocTemp =  purrr::map(MatLocTemp,
        ~ as.numeric( .x %*% params.temp$loc) )
      ) %>%
      dplyr::select( - MatLocTemp)

    locpars <- purrr::map2(.x = dataprep$LocTemp, .y = dataprep$loc_sp, ~ .x + .y)


  } else {
    locpars <- dataprep$loc_sp
  }


  if(!is.null(scale.temp.form)) {
    scale.temp.form <- update(scale.temp.form,  ~ . + 0)

    dataprep <- dataprep %>%
      dplyr::mutate(ScaleTemp =  purrr::map(MatScaleTemp,
                                          ~ as.numeric( .x %*% params.temp$scale) )
      ) %>%
      dplyr::select( - MatScaleTemp)

    scalepars <- purrr::map2(.x = dataprep$ScaleTemp, .y = dataprep$scale_sp, ~ .x + .y)


  } else {
    scalepars <- dataprep$scale_sp
  }

  yy <- purrr::map(dataprep$uniq_data, ~ .x$slbm)

  link.scalepars <- scale.link$linkinv(scalepars)



  if(abs(params.sp$shape ) < 1e-8) {
    yy <-  purrr::pmap(list(.x = yy, .y = locpars, .z = link.scalepars),
                       function(.x, .y, .z) {
                         exp(- (.x - .y)/.z)
                       } )
  } else {
    yy <- purrr::pmap(list(.x = yy, .y = locpars, .z = link.scalepars),
                function(.x, .y, .z) {
                  (1 + params.sp$shape*( (.x - .y)/.z))^(-1/params.sp$shape)
                })
  }

  stopyn <- (any(unlist(purrr::map(yy, ~ any(.x < 0, na.rm = TRUE)))) | any(unlist(link.scalepars )<= 0) )
  if(stopyn) {
    return( 1e+10 )
    } else {

      yy <- purrr::pmap( list(.x = yy, .y = link.scalepars , .z = dataprep$uniq_data),
                         .f = function(.x, .y, .z) {
                           u <- (log(.y) - (params.sp$shape +1)*log(.x) + .x)
                           u*.z$n
                         }

      )

    sum(unlist(yy))
    }




}

# only formulas for scale are given; location is const*scale
nll_spat_temp_sl_hom <- function(params, scale.sp.form,
                                  scale.temp.form = NULL,
                                  dataprep, spat.cov, scale.link = make.link("log")){


  params.sp <- list(disp = params[startsWith(names(params),"disp")],
                    scale = params[startsWith(names(params),"scale")],
                    shape = params[startsWith(names(params),"shape")])

  params.temp <- list(scale = params[startsWith(names(params),"tempScale")])

  params.temp <- purrr::map(params.temp, function(.x) ifelse(length(.x) == 0, 0, .x) )

  model.scale.sp <- dataprep$MatScaleSp

  #  oc_sp <- data.frame(loc_sp = model.loc.sp %*% params.sp$loc)

  scale_sp <-  data.frame(scale_sp = model.scale.sp %*% params.sp$scale)

  loc_sp <- params.sp$disp*scale_sp
  names(loc_sp) <- "loc_sp"
  dataprep <- dataprep$data

  dataprep <- dataprep %>% dplyr::bind_cols( loc_sp, scale_sp )



  if(!is.null(scale.temp.form)) {
    scale.temp.form <- update(scale.temp.form,  ~ . + 0)

    dataprep <- dataprep %>%
      dplyr::mutate(ScaleTemp =  purrr::map(MatScaleTemp,
                                            ~ as.numeric( .x %*% params.temp$scale) )
      ) %>%
      dplyr::select( - MatScaleTemp)

    scalepars <- purrr::map2(.x = dataprep$ScaleTemp, .y = dataprep$scale_sp, ~ .x + .y)


  } else {
    scalepars <- dataprep$scale_sp
  }

  yy <- purrr::map(dataprep$uniq_data, ~ .x$slbm)

  link.scalepars <- scale.link$linkinv(scalepars)



  if(abs(params.sp$shape ) < 1e-8) {
    yy <-  purrr::pmap(list(.x = yy, .z = link.scalepars),
                       function(.x,  .z) {
                         exp(- (.x /.z -  params.sp$disp))
                       } )
  } else {
    yy <- purrr::pmap(list(.x = yy, .z = link.scalepars),
                      function(.x, .z) {
                        (1 + params.sp$shape*( .x /.z -  params.sp$disp))^(-1/params.sp$shape)
                      })
  }

  stopyn <- (any(unlist(purrr::map(yy, ~ any(.x < 0, na.rm = TRUE)))) | any(unlist(link.scalepars )<= 0) )
  if(stopyn) {
    return( 1e+10 )
  } else {

    yy <- purrr::pmap( list(.x = yy, .y = link.scalepars , .z = dataprep$uniq_data),
                       .f = function(.x, .y, .z) {
                         u <- (log(.y) - (params.sp$shape +1)*log(.x) + .x)
                         u*.z$n
                       }

    )

    sum(unlist(yy))
  }




}

nll_spat_IF <- function(params ,scale.sp.form,
                        scale.temp.form = NULL,
                        data, spat.cov, temp.cov = NULL){

  n.dat <- nrow(data)
  n.site <- ncol(data)
  n.spat <- nrow(spat.cov)
  n.temp <- ncol(temp.cov)


  params.sp <- list(disp = params[startsWith(names(params),"disp")],
                    scale = params[startsWith(names(params),"scale")],
                    shape = params[startsWith(names(params),"shape")])

  params.temp <- list( scale = params[startsWith(names(params),"tempScale")])

  params.temp <- purrr::map(params.temp, function(.x) ifelse(length(.x) == 0, 0, .x) )

  model.scale.sp <- model.matrix(scale.sp.form, model.frame(scale.sp.form,
                                                            data = spat.cov, na.action = na.pass))


  if( !is.null(scale.temp.form)){
    scale.temp.form <- update(scale.temp.form,  ~ . + 0)

    model.scale.temp <- model.matrix(scale.temp.form, model.frame(scale.temp.form,
                                                                  data =  data.frame(temp.cov), na.action = na.pass))
  }
  else{
    model.scale.temp <- matrix(rep(0, n.dat), ncol = 1)
  }


  param.scale <-  model.scale.sp %*% params.sp$scale

  param.scale <- matrix(rep( param.scale, n.dat), ncol = n.site, byrow = TRUE) +
    matrix(rep(model.scale.temp %*% params.temp$scale, n.site), ncol = n.site, byrow = F)


  if(abs(params.sp$shape) < 1e-8){
    yy <- data/exp(param.scale) - params.sp$disp
    likh <-  sum(param.scale  + yy + exp(-yy), na.rm = TRUE)
  }
  else{
    yy <- 1+ params.sp$shape*(data/exp(param.scale) - params.sp$disp)

    if(any(yy <0, na.rm = TRUE)){ likh <- 1e+10} else {
      likh <- sum(param.scale  + (1/params.sp$shape +1)*log(yy) + yy^(-1/params.sp$shape),
                  na.rm = TRUE)
    }
  }

  return(likh)
}

### Gradients
# to do : case shape = 0, case IF assumption
#################################################################################
gr_spat_gev <- function(params, loc.sp.form , scale.sp.form,
                        loc.temp.form = NULL, scale.temp.form = NULL,
                        data, spat.cov, temp.cov = NULL ){
  #
  # param <- list(loc = params[startsWith(names(params),"loc")],
  #               scale = params[startsWith(names(params),"scale")],
  #               shape = params[startsWith(names(params),"shape")])

  params.sp <- list(loc = params[startsWith(names(params),"loc")],
                    scale = params[startsWith(names(params),"scale")],
                    shape = params[startsWith(names(params),"shape")])

  params.temp <- list(loc = params[startsWith(names(params),"tempLoc")],
                      scale = params[startsWith(names(params),"tempScale")])

  params.temp <- purrr::map(params.temp, function(.x) ifelse(length(.x) == 0, 0, .x) )

  n.dat <- nrow(data)
  d.dat <- n.site <- ncol(data)
  n.loc.sp <- length(all.vars(loc.sp.form)) +1
  n.scale.sp <- length(all.vars(scale.sp.form)) +1

  n.loc.temp <- length(all.vars(loc.temp.form))
  n.scale.temp <- length(all.vars(scale.temp.form))

  n.param <- n.loc.sp+ n.scale.sp + n.loc.temp + n.scale.temp + 1

  model.loc.sp <- model.matrix(loc.sp.form, model.frame(loc.sp.form,
                                                        data = spat.cov, na.action = na.pass))

  model.scale.sp <- model.matrix(scale.sp.form, model.frame(scale.sp.form,
                                                            data = spat.cov, na.action = na.pass))

  if( !is.null(loc.temp.form)){
    loc.temp.form <- update(loc.temp.form,  ~ . + 0)

    model.loc.temp <- model.matrix(loc.temp.form, model.frame(loc.temp.form,
                                                              data = data.frame(temp.cov),
                                                              na.action = na.pass))
  }
  else{
    model.loc.temp <- matrix(rep(0, n.dat), ncol = 1)
  }

  if( !is.null(scale.temp.form)){
    scale.temp.form <- update(scale.temp.form,  ~ . + 0)
    scale.temp.form <- update(scale.temp.form,  ~ . + 0)
    model.scale.temp <- model.matrix(scale.temp.form, model.frame(scale.temp.form,
                                                                  data =  data.frame(temp.cov), na.action = na.pass))
  }
  else{
    model.scale.temp <- matrix(rep(0, n.dat), ncol = 1)
  }

  param.loc <- model.loc.sp %*% params.sp$loc
  param.loc <- matrix(rep(param.loc, n.dat), ncol = n.site, byrow = TRUE) +
    matrix(rep(model.loc.temp %*% params.temp$loc, n.site), ncol = n.site, byrow =F)

  param.scale <-  model.scale.sp %*% params.sp$scale

  param.scale <- matrix(rep( param.scale, n.dat), ncol = n.site, byrow = TRUE) +
    matrix(rep(model.scale.temp %*% params.temp$scale, n.site), ncol = n.site, byrow = F) # exp(model.scale.sp %*% paramH0$scale)


  gr_loc <- numeric(n.loc.sp)
  gr_scale <- numeric(n.scale.sp)

  gr_loc_temp <- numeric(n.loc.temp)
  gr_scale_temp <- numeric(n.scale.temp)

  # sigma.mat <- matrix(rep( model.scale.sp %*% param$scale, n.dat),
  #                     ncol = d.dat, byrow = TRUE)
  # loc.mat <- matrix(rep(model.loc.sp %*% param$loc, n.dat),
  #                   ncol = d.dat, byrow = TRUE)


  if(abs(params.sp$shape) < 1e-8){
    print("no")
    yy <- data/exp(sigma.mat) - loc.mat/exp(sigma.mat)
    for( i.loc in 1:n.loc.sp){
      mat.cvrt <- matrix(rep(model.loc.sp[, i.loc], n.dat), ncol = d.dat, byrow = TRUE)
      gr_loc[i.loc] <-  sum( (1/params.sp$shape +1)*yy *(-params.sp$shape)/exp(sigma.mat)*mat.cvrt +
                               yy^(-1/params.sp$shape -1)*exp(-sigma.mat)*mat.cvrt, na.rm = TRUE)
    }
    for(i.scale in 1:n.scale.sp){
      gr_scale[i.scale] <- NA # sum( (1/param$shape +1)*yy *(-param$shape)/exp(sigma.mat)*mat.cvrt +
      # yy^(-1/shape -1)*exp(-sigma.mat)*mat.cvrt, na.rm = TRUE)
    }

  }
  else{

    yy <- 1+ params.sp$shape*(data/exp(param.scale) - param.loc/exp(param.scale))

    for( i.loc in 1:n.loc.sp){
      mat.cvrt <- matrix(rep(model.loc.sp[, i.loc], n.dat), ncol = d.dat, byrow = TRUE)
      gr_loc[i.loc] <-  sum( (1/params.sp$shape +1)*yy^(-1) *(-params.sp$shape)/exp(param.scale)*mat.cvrt +
                               yy^(-1/params.sp$shape -1)*exp(-param.scale)*mat.cvrt, na.rm = TRUE)

    }
    for( i.scale in 1:n.scale.sp){
      mat.cvrt <- matrix(rep(model.scale.sp[, i.scale], n.dat), ncol = d.dat, byrow = TRUE)
      gr_scale[i.scale] <- sum( mat.cvrt +
                                  params.sp$shape*mat.cvrt/exp(param.scale)*(-data + param.loc)*
                                  ( (1/params.sp$shape +1)*yy^(-1) -1/params.sp$shape*yy^(-1/params.sp$shape -1)), na.rm = TRUE)

    }
    if(n.scale.temp > 0 ){
      for( i.temp.scale in 1:n.scale.temp){
        mat.cvrt <- matrix(rep(model.scale.temp[, i.temp.scale], d.dat), ncol = d.dat, byrow = F)
        gr_scale_temp[i.temp.scale] <- sum( mat.cvrt +
                                              params.sp$shape*mat.cvrt/exp(param.scale)*(-data + param.loc)*
                                              ( (1/params.sp$shape +1)*yy^(-1) -1/params.sp$shape*yy^(-1/params.sp$shape -1)), na.rm = TRUE)

      }}
    if(n.loc.temp > 0 ){
      for( i.temp.loc in 1:n.loc.temp){
        mat.cvrt <- matrix(rep(model.loc.temp[, i.temp.loc], d.dat), ncol = d.dat, byrow = F)
        gr_loc_temp[i.temp.loc] <-  sum( (1/params.sp$shape +1)*yy^(-1) *(-params.sp$shape)/exp(param.scale)*mat.cvrt +
                                           yy^(-1/params.sp$shape -1)*exp(-param.scale)*mat.cvrt, na.rm = TRUE)

      }}

    gr_shape <- sum(-1/params.sp$shape^2*log(yy) + (1/params.sp$shape +1)*(yy -1)/params.sp$shape/yy +
                      yy^(-1/params.sp$shape)*(1/params.sp$shape^2 *log(yy) -
                                                 1/params.sp$shape *(yy-1)/params.sp$shape/yy),
                    na.rm = TRUE)
  }
  # print(c(gr_loc, gr_scale, gr_shape))
  return(c(gr_loc, gr_scale, gr_shape, gr_loc_temp, gr_scale_temp))
}


gr_spat_gev_IF <- function(params,  scale.sp.form,
                           scale.temp.form = NULL,
                           data, spat.cov, temp.cov = NULL ){


  params.sp <- list(disp = params[startsWith(names(params),"disp")],
                    scale = params[startsWith(names(params),"scale")],
                    shape = params[startsWith(names(params),"shape")])

  params.temp <- list(
    scale = params[startsWith(names(params),"tempScale")])

  params.temp <- purrr::map(params.temp, function(.x) ifelse(length(.x) == 0, 0, .x) )

  n.dat <- nrow(data)
  d.dat <- n.site <-  ncol(data)

  n.scale.sp <- length(all.vars(scale.sp.form)) +1


  n.scale.temp <- length(all.vars(scale.temp.form))

  n.param <-  1 + n.scale.sp + n.scale.temp + 1

  model.scale.sp <- model.matrix(scale.sp.form, model.frame(scale.sp.form,
                                                            data = spat.cov, na.action = na.pass))

  if( !is.null(scale.temp.form)){
    scale.temp.form <- update(scale.temp.form,  ~ . + 0)
    scale.temp.form <- update(scale.temp.form,  ~ . + 0)
    model.scale.temp <- model.matrix(scale.temp.form, model.frame(scale.temp.form,
                                                                  data =  data.frame(temp.cov), na.action = na.pass))
  }
  else{
    model.scale.temp <- matrix(rep(0, n.dat), ncol = 1)
  }

  param.scale <-  model.scale.sp %*% params.sp$scale

  param.scale <- matrix(rep( param.scale, n.dat), ncol = n.site, byrow = TRUE) +
    matrix(rep(model.scale.temp %*% params.temp$scale, n.site), ncol = n.site, byrow = F) # exp(model.scale.sp %*% paramH0$scale)

  gr_scale <- numeric(n.scale.sp)

  gr_scale_temp <- numeric(n.scale.temp)

  if(abs(params.sp$shape) < 1e-8){
    print("no")
    yy <- data/exp(sigma.mat) - param$disp

    gr_disp <-  sum( (1/param$shape +1)*yy *(-param$shape)/exp(sigma.mat)*mat.cvrt +
                       yy^(-1/param$shape -1)*exp(-sigma.mat)*mat.cvrt, na.rm = TRUE)

    for(i.scale in 1:n.scale.sp){
      gr_scale[i.scale] <- NA # sum( (1/param$shape +1)*yy *(-param$shape)/exp(sigma.mat)*mat.cvrt +
      # yy^(-1/shape -1)*exp(-sigma.mat)*mat.cvrt, na.rm = TRUE)
    }

  }
  else{

    yy <- 1+ params.sp$shape*(data/exp(param.scale) - params.sp$disp)

    gr_disp <-  sum( (1/params.sp$shape +1)*(-params.sp$shape)*yy^(-1) +
                       yy^(-1/params.sp$shape -1),  na.rm = TRUE)


    for( i.scale in 1:n.scale.sp){
      mat.cvrt <- matrix(rep(model.scale.sp[, i.scale], n.dat), ncol = d.dat, byrow = TRUE)
      gr_scale[i.scale] <- sum( mat.cvrt +
                                  params.sp$shape*mat.cvrt/exp(param.scale)*(-data )*
                                  ( (1/params.sp$shape +1)*yy^(-1)
                                    -1/params.sp$shape*yy^(-1/params.sp$shape -1)), na.rm = TRUE)

    }
    if(n.scale.temp > 0 ){
      for( i.temp.scale in 1:n.scale.temp){
        mat.cvrt <- matrix(rep(model.scale.temp[, i.temp.scale], d.dat), ncol = d.dat, byrow = F)
        gr_scale_temp[i.temp.scale] <- sum( mat.cvrt +
                                              params.sp$shape*mat.cvrt/exp(param.scale)*(-data)*
                                              ( (1/params.sp$shape +1)*yy^(-1)
                                                -1/params.sp$shape*yy^(-1/params.sp$shape -1)), na.rm = TRUE)

      }}
    gr_shape <- sum(-1/params.sp$shape^2*log(yy) + (1/params.sp$shape +1)*(yy -1)/params.sp$shape/yy +
                      yy^(-1/params.sp$shape)*(1/params.sp$shape^2 *log(yy)
                                               -1/params.sp$shape *(yy-1)/params.sp$shape/yy),
                    na.rm = TRUE)
  }
  # print(c(gr_loc, gr_scale, gr_shape))
  return(c(gr_disp, gr_scale, gr_shape, gr_scale_temp))
}

gr_spat_gev_sl <- function(params, loc.sp.form , scale.sp.form,
                        loc.temp.form = NULL, scale.temp.form = NULL,
                        data, spat.cov, temp.cov = NULL ){

  params.sp <- list(loc = params[startsWith(names(params),"loc")],
                    scale = params[startsWith(names(params),"scale")],
                    shape = params[startsWith(names(params),"shape")])

  params.temp <- list(loc = params[startsWith(names(params),"tempLoc")],
                      scale = params[startsWith(names(params),"tempScale")])

  params.temp <- purrr::map(params.temp, function(.x) ifelse(length(.x) == 0, 0, .x) )

  model.loc.sp <- dataprep$MatLocSp

  model.scale.sp <- dataprep$MatScaleSp

  loc_sp <- data.frame(loc_sp = model.loc.sp %*% params.sp$loc)

  scale_sp <-  data.frame(scale_sp = model.scale.sp %*% params.sp$scale)

  dataprep <- dataprep$data

  dataprep <- dataprep %>% dplyr::bind_cols(loc_sp, scale_sp)



  if(!is.null(loc.temp.form)) {
    loc.temp.form <- update(loc.temp.form,  ~ . + 0)
    dataprep <- dataprep %>%
      dplyr::mutate(LocTemp =  purrr::map(MatLocTemp,
                                          ~ as.numeric( .x %*% params.temp$loc) )
      ) %>%
      dplyr::select( - MatLocTemp)

    locpars <- purrr::map2(.x = dataprep$LocTemp, .y = dataprep$loc_sp, ~ .x + .y)


  } else {
    locpars <- dataprep$loc_sp
  }


  if(!is.null(scale.temp.form)) {
    scale.temp.form <- update(scale.temp.form,  ~ . + 0)

    dataprep <- dataprep %>%
      dplyr::mutate(ScaleTemp =  purrr::map(MatScaleTemp,
                                            ~ as.numeric( .x %*% params.temp$scale) )
      ) %>%
      dplyr::select( - MatScaleTemp)

    scalepars <- purrr::map2(.x = dataprep$ScaleTemp, .y = dataprep$scale_sp, ~ .x + .y)


  } else {
    scalepars <- dataprep$scale_sp
  }

  yy <- purrr::map(dataprep$uniq_data, ~ .x$slbm)

  link.scalepars <- scale.link$linkinv(scalepars)


  gr_loc <- numeric(n.loc.sp)
  gr_scale <- numeric(n.scale.sp)

  gr_loc_temp <- numeric(n.loc.temp)
  gr_scale_temp <- numeric(n.scale.temp)

  # sigma.mat <- matrix(rep( model.scale.sp %*% param$scale, n.dat),
  #                     ncol = d.dat, byrow = TRUE)
  # loc.mat <- matrix(rep(model.loc.sp %*% param$loc, n.dat),
  #                   ncol = d.dat, byrow = TRUE)


  if(abs(params.sp$shape) < 1e-8){
    print("no")
    yy <- data/exp(sigma.mat) - loc.mat/exp(sigma.mat)
    for( i.loc in 1:n.loc.sp){
      mat.cvrt <- matrix(rep(model.loc.sp[, i.loc], n.dat), ncol = d.dat, byrow = TRUE)
      gr_loc[i.loc] <-  sum( (1/params.sp$shape +1)*yy *(-params.sp$shape)/exp(sigma.mat)*mat.cvrt +
                               yy^(-1/params.sp$shape -1)*exp(-sigma.mat)*mat.cvrt, na.rm = TRUE)
    }
    for(i.scale in 1:n.scale.sp){
      gr_scale[i.scale] <- NA # sum( (1/param$shape +1)*yy *(-param$shape)/exp(sigma.mat)*mat.cvrt +
      # yy^(-1/shape -1)*exp(-sigma.mat)*mat.cvrt, na.rm = TRUE)
    }

  }
  else{

    yy <- 1+ params.sp$shape*(data/exp(param.scale) - param.loc/exp(param.scale))

    for( i.loc in 1:n.loc.sp){
      mat.cvrt <- matrix(rep(model.loc.sp[, i.loc], n.dat), ncol = d.dat, byrow = TRUE)
      gr_loc[i.loc] <-  sum( (1/params.sp$shape +1)*yy^(-1) *(-params.sp$shape)/exp(param.scale)*mat.cvrt +
                               yy^(-1/params.sp$shape -1)*exp(-param.scale)*mat.cvrt, na.rm = TRUE)

    }
    for( i.scale in 1:n.scale.sp){
      mat.cvrt <- matrix(rep(model.scale.sp[, i.scale], n.dat), ncol = d.dat, byrow = TRUE)
      gr_scale[i.scale] <- sum( mat.cvrt +
                                  params.sp$shape*mat.cvrt/exp(param.scale)*(-data + param.loc)*
                                  ( (1/params.sp$shape +1)*yy^(-1) -1/params.sp$shape*yy^(-1/params.sp$shape -1)), na.rm = TRUE)

    }
    if(n.scale.temp > 0 ){
      for( i.temp.scale in 1:n.scale.temp){
        mat.cvrt <- matrix(rep(model.scale.temp[, i.temp.scale], d.dat), ncol = d.dat, byrow = F)
        gr_scale_temp[i.temp.scale] <- sum( mat.cvrt +
                                              params.sp$shape*mat.cvrt/exp(param.scale)*(-data + param.loc)*
                                              ( (1/params.sp$shape +1)*yy^(-1) -1/params.sp$shape*yy^(-1/params.sp$shape -1)), na.rm = TRUE)

      }}
    if(n.loc.temp > 0 ){
      for( i.temp.loc in 1:n.loc.temp){
        mat.cvrt <- matrix(rep(model.loc.temp[, i.temp.loc], d.dat), ncol = d.dat, byrow = F)
        gr_loc_temp[i.temp.loc] <-  sum( (1/params.sp$shape +1)*yy^(-1) *(-params.sp$shape)/exp(param.scale)*mat.cvrt +
                                           yy^(-1/params.sp$shape -1)*exp(-param.scale)*mat.cvrt, na.rm = TRUE)

      }}

    gr_shape <- sum(-1/params.sp$shape^2*log(yy) + (1/params.sp$shape +1)*(yy -1)/params.sp$shape/yy +
                      yy^(-1/params.sp$shape)*(1/params.sp$shape^2 *log(yy) -
                                                 1/params.sp$shape *(yy-1)/params.sp$shape/yy),
                    na.rm = TRUE)
  }
  # print(c(gr_loc, gr_scale, gr_shape))
  return(c(gr_loc, gr_scale, gr_shape, gr_loc_temp, gr_scale_temp))
}

#################################################################################
# the estimated scale parameters are the sigma_i of
# log(sigma(s)) = sigma_0 + sigma_1*cvrt1 + sigma_2*scrt2 + ...

# ... contains arguments for starting values under IF assumption
get_start_vals <- function(data, loc.sp.form = ~ 1, scale.sp.form = ~ 1,
                           loc.temp.form = NULL, scale.temp.form = NULL,
                           spat.cov, temp.cov = NULL, method = "LeastSq",
                           print_start_vals = TRUE, scale.link = make.link("log"),
                           ...){

  if( !(method %in% c("LeastSq", "LeastSqTemp", "spatialGev", "fgevMean"))) {
    stop("This is not a valid method for computing starting values.")
  }

  n.loc.temp <- length(all.vars(loc.temp.form))
  n.scale.temp <- length(all.vars(scale.temp.form))
  n.loc.sp <- length(all.vars(loc.sp.form)) + 1
  n.scale.sp <- length(all.vars(scale.sp.form)) + 1

  n.param <- n.loc.sp + n.scale.sp + 1

  IFargs <- list(...)

  # sets coefficients of temporal covariates to 0
  if(method == "LeastSq") {

    if(exists("type", where = IFargs) & !is.null(IFargs$type)) {

      #browser()
      YY <- fgev_hom(data, method = "BFGS", scale.link = scale.link)$mle

      .datfr <- data.frame(t(YY)) %>% dplyr::bind_cols(spat.cov)

      lm.sigma <- reformulate(deparse(scale.sp.form[[2]]), paste0(scale.link$name, "(scale)"))

      st_disp <- mean(.datfr$loc/.datfr$scale)
      names(st_disp) <- "disp"
      st_sigmas <- lm(lm.sigma, data = .datfr)$coefficients
      names(st_sigmas) <-  paste0("scale", 0:(length(st_sigmas)-1))
      st_shape <- mean(.datfr$shape)
      names(st_shape) <- "shape"

      if(is.null(scale.temp.form)){
        start_vals <- c(st_disp, st_sigmas, st_shape)

      }
      else{

        n.scale.temp <- length(all.vars(scale.temp.form))

        st_temp_scale <- rep(0, n.scale.temp)
        names(st_temp_scale) <- paste0("tempScale", 1:n.scale.temp)

        start_vals <- c(st_disp, st_sigmas, st_shape, st_temp_scale)
      }

    }

    else {

        YY <- apply(data, 2,  function(.x) evd::fgev(.x, std.err = F)$estimate)

        .datfr <- data.frame(t(YY)) %>% dplyr::bind_cols(spat.cov)


        lm.mu <- reformulate(deparse(loc.sp.form[[2]]), "loc")
       # lm.sigma <- reformulate(deparse(scale.sp.form[[2]]), "log(scale)")
        lm.sigma <- reformulate(deparse(scale.sp.form[[2]]), paste0(scale.link$name, "(scale)"))

        st_mus <- lm(lm.mu, data = .datfr )$coefficients
        names(st_mus) <- paste0("loc", 0:(length(st_mus)-1))

        st_sigmas <- lm(lm.sigma, data = .datfr)$coefficients
        names(st_sigmas) <-  paste0("scale", 0:(length(st_sigmas)-1))

        st_shape <- mean(.datfr$shape)
        names(st_shape) <- "shape"

        if(n.loc.temp == 0 & n.scale.temp == 0) {
          start_vals <- c(st_mus, st_sigmas, st_shape)
        }
        else { # initialise zeros for temporal coefficients

          st_temp_loc <- rep(0, n.loc.temp)
          if(n.loc.temp > 0){
            names(st_temp_loc) <- paste0("tempLoc", 1:n.loc.temp)
          }

          st_temp_scale <- rep(0, n.scale.temp)
          if(n.scale.temp > 0 ){
            names(st_temp_scale) <- paste0("tempScale", 1:n.scale.temp)
          }

          start_vals <- c(st_mus, st_sigmas, st_shape, st_temp_loc, st_temp_scale)
        }
    }

  }

  # estimates location trend based (one covariate only)
  if(method == "LeastSqTemp") {

    YY <- apply(data, 2,
                function(.x) evd::fgev(.x, nsloc = data.frame(temp.cov), std.err = F)$estimate)

    .datfr <- data.frame(t(YY)) %>% dplyr::bind_cols(spat.cov)


    lm.mu <- reformulate(deparse(loc.sp.form[[2]]), "loc")
      # lm.sigma <- reformulate(deparse(scale.sp.form[[2]]), "log(scale)")
    lm.sigma <- reformulate(deparse(scale.sp.form[[2]]), paste0(scale.link$name, "(scale)"))

    st_mus <- lm(lm.mu, data = .datfr )$coefficients
    names(st_mus) <- paste0("loc", 0:(length(st_mus)-1))

    st_sigmas <- lm(lm.sigma, data = .datfr)$coefficients
    names(st_sigmas) <-  paste0("scale", 0:(length(st_sigmas)-1))

    st_shape <- mean(.datfr$shape)
    names(st_shape) <- "shape"

   # initialise zeros for temporal coefficients

    st_temp_loc <- mean(.datfr[, 2])
    names(st_temp_loc) <- "tempLoc1"


    st_temp_scale <- rep(0, n.scale.temp)
    if(n.scale.temp > 0 ){
      names(st_temp_scale) <- paste0("tempScale", 1:n.scale.temp)
    }

    start_vals <- c(st_mus, st_sigmas, st_shape, st_temp_loc, st_temp_scale)


  }

  # can estimate arbitrary amount of (trend) coefficients.
  # No link function for scale parameter
  if(method == "spatialGev"){

      if (is.null(temp.cov) & any( c(n.loc.temp, n.scale.temp) > 0)) {

        warning("No temporal covariate was supplied to compute starting values with
             this method. Initial Temporal parameters are set to 0.")

        start_params <- SpatialExtremes::fitspatgev(data = data,
                                                    covariables = as.matrix(spat.cov),
                                                    loc.form =  loc.sp.form,
                                                    scale.form = scale.sp.form,
                                                    shape.form = ~ 1)$fitted.values
      } else {

        start_params <- SpatialExtremes::fitspatgev(data = data,
                                                    covariables = as.matrix(spat.cov),
                                                    loc.form =  loc.sp.form,
                                                    scale.form = scale.sp.form,
                                                    shape.form = ~ 1,
                                                    temp.form.loc = loc.temp.form,
                                                    temp.form.scale = scale.temp.form,
                                                    temp.cov = as.matrix(temp.cov))$fitted.values
      }

      st_param <- list(loc = start_params[startsWith(names(start_params),"loc")],
                       scale = start_params[startsWith(names(start_params),"scale")],
                       shape = start_params[startsWith(names(start_params),"shape")])


      start_vals <- numeric(n.param)
      start_vals[1:n.loc.sp] <- st_param$loc

      # caution: SpatialExtremes doesn't use link functionv for scale parameter,
      # method best to use without link function (identity) or when scale is constant
      start_vals[(n.loc.sp +1):(n.param -1)] <- scale.link$linkfun(st_param$scale)

      start_vals[n.param] <- st_param$shape
      names(start_vals) <- c(paste0("loc", 0:(n.loc.sp - 1)),
                              paste0("scale", 0:(n.scale.sp - 1)), "shape")

      if(n.loc.temp  > 0) {
        sv_temp_loc <- numeric(n.loc.temp)
        if(!(is.null(temp.cov))) {
          sv_temp_loc <- start_params[startsWith(names(start_params),"tempCoeffLoc")]
        }
        names(sv_temp_loc) <- paste0("tempLoc", 1:n.loc.temp)
        start_vals <- c(start_vals, sv_temp_loc)
      }
      if(n.scale.temp  > 0) {
        sv_temp_scale <- numeric(n.scale.temp)
        if(!(is.null(temp.cov))) {
          sv_temp_scale <- start_params[startsWith(names(start_params),"tempCoeffScale")]
        }
        names(sv_temp_scale) <- paste0("tempScale", 1:n.scale.temp)
        start_vals <- c(start_vals, sv_temp_scale)
      }
    }

  # most basic method. Estimates only stationary parameters, other coefficients are 0
  if(method == "fgevMean"){

      start_st <-  apply(data,2, function(.x){evd::fgev(.x[!is.na(.x)],  std.err = F)})
      start_st <- purrr::map(purrr::map_dfr(start_st, function(x) x$estimate), mean)
      start_vals <- numeric(n.param)
      start_vals[1] <- start_st$loc
      start_vals[n.loc.sp + 1] <- scale.link$linkfun(start_st$scale)
      start_vals[n.param] <- start_st$shape
      names(start_vals) <- c(paste0("loc", 0:(n.loc.sp-1)),
                              paste0("scale", 0:(n.scale.sp -1)), "shape")


      if(n.loc.temp  > 0) {
        sv_temp_loc <- rep(0, n.loc.temp)
        names(sv_temp_loc) <- paste0("tempLoc", 1:n.loc.temp)
        start_vals <- c(start_vals, sv_temp_loc)
      }
      if(n.scale.temp  > 0) {
        sv_temp_scale <- rep(0, n.scale.temp)
        names(sv_temp_scale) <- paste0("tempScale", 1:n.scale.temp)
        start_vals <- c(start_vals, sv_temp_scale)
      }
  }

  if(print_start_vals){
    message( "start values are", "\n")
    print(c(  round(start_vals,3) ))
  }

  start_vals

}

fit_spgev <- function(data, loc.sp.form, scale.sp.form,
                      loc.temp.form = NULL, scale.temp.form = NULL,
                      spat.cov, temp.cov = NULL, start_vals = NULL, use_gr = TRUE,
                      method = "BFGS", st_val_meth = "LeastSq", print_start_vals = TRUE, ...){

  n.dat <- nrow(data)
  d.dat <- ncol(data)
  n.spat <- nrow(spat.cov)
  n.temp <- nrow(temp.cov)

  if (d.dat != n.spat)
    stop("'data' and 'covariates' doesn't match")
  use.temp.cov <- c(!is.null(loc.temp.form), !is.null(scale.temp.form))
  if (any(use.temp.cov) && (n.dat != n.temp))
    stop("'data' and 'temp.cov' doesn't match")
  if (any(use.temp.cov) && is.null(temp.cov))
    stop("'temp.cov' must be supplied if at least one temporal formula is given")

  n.loc.sp <- length(all.vars(loc.sp.form)) +1
  n.scale.sp <- length(all.vars(scale.sp.form)) +1

  n.param <- n.loc.sp + n.scale.sp + 1


  if(is.null(start_vals)){

    if(st_val_meth == "LeastSq"){
      start_vals <-  get_start_vals(data = data, loc.sp.form = loc.sp.form, scale.sp.form = scale.sp.form,
                                    loc.temp.form = loc.temp.form, scale.temp.form = scale.temp.form,
                                    spat.cov = spat.cov, temp.cov = temp.cov,
                                    method = "LeastSq", ...)
      if(print_start_vals){cat("start values are", start_vals, "\n") }
    }
    if(st_val_meth == "sp_gev"){
      start_params <- fitspatgev(data = data, covariables = as.matrix(spat.cov),
                                 loc.form =  loc.sp.form, scale.form = scale.sp.form,
                                 shape.form = ~ 1)$fitted.values

      st_param <- list(loc= start_params[startsWith(names(start_params),"loc")],
                       scale = start_params[startsWith(names(start_params),"scale")],
                       shape = start_params[startsWith(names(start_params),"shape")])


      start_vals <- numeric(n.param)
      start_vals[1:n.loc.sp] <-st_param$loc
      start_vals[(n.loc.sp +1):(n.param -1)] <- log(st_param$scale)
      start_vals[n.param] <- st_param$shape
      names(start_vals ) <- c(paste0("loc", 0:(n.loc.sp-1)),
                              paste0("scale", 0:(n.scale.sp -1)), "shape")

      if(print_start_vals){
        cat("start values are ", start_vals , "\n")
      }
    }
    if(st_val_meth == "fgev"){

      start_st <-  apply(data,2, function(.x){evd::fgev(.x[!is.na(.x)],  std.err = F)})
      start_st <- purrr::map(purrr::map_dfr(start_st, function(x) x$estimate), mean)
      start_vals <- numeric(n.param)
      start_vals[1] <- start_st$loc
      start_vals[n.loc.sp +1] <- log(start_st$scale)
      start_vals[n.param] <- start_st$shape
      names(start_vals ) <- c(paste0("loc", 0:(n.loc.sp-1)),
                              paste0("scale", 0:(n.scale.sp -1)), "shape")
      if(print_start_vals){
        cat("start values are ", start_vals)
      }
    }
  }
  else( cat("start values are ", start_vals))
  if(use_gr){
    mlest <- optim(start_vals, fn = nll_spat, gr =  gr_spat_gev,
                   loc.sp.form = loc.sp.form,
                   scale.sp.form = scale.sp.form,
                   loc.temp.form = loc.temp.form, scale.temp.form = scale.temp.form,
                   data = data, spat.cov = spat.cov, temp.cov = temp.cov,
                   method = method,control = list( maxit = 1000))
  }
  else{
    mlest <- optim(start_vals, fn = nll_spat,
                   loc.sp.form = loc.sp.form,
                   scale.sp.form = scale.sp.form,
                   loc.temp.form = loc.temp.form, scale.temp.form = scale.temp.form,
                   data = data, spat.cov = spat.cov, temp.cov = temp.cov,
                   method = method,control = list( maxit = 1000))
  }
  return(list(mle = mlest$par, nllh = mlest$value, conv = mlest$convergence))
}

fit_spgev_IF <- function(data,  scale.sp.form,
                         scale.temp.form = NULL,
                         spat.cov, temp.cov = NULL, start_vals = NULL,
                         use_gr = TRUE, method = "BFGS", st_val_meth = "LeastSq", start_methIF = "Lmom"){


  n.dat <- nrow(data)
  d.dat <- ncol(data)
  n.spat <- nrow(spat.cov)
  n.temp <- nrow(temp.cov)

  if (d.dat != n.spat)
    stop("'data' and 'covariates' doesn't match")
  use.temp.cov <- c(!is.null(scale.temp.form))
  if (any(use.temp.cov) && (n.dat != n.temp))
    stop("'data' and 'temp.cov' doesn't match")
  if (any(use.temp.cov) && is.null(temp.cov))
    stop("'temp.cov' must be supplied if at least one temporal formula is given")


  n.scale.sp <- length(all.vars(scale.sp.form)) +1

  n.param <- 1+ n.scale.sp + 1


  if(is.null(start_vals)){
    if(st_val_meth == "LeastSq"){
      start_vals <-  get_start_vals(data = data, scale.sp.form = scale.sp.form,
                                    loc.temp.form = loc.temp.form, scale.temp.form = scale.temp.form,
                                    spat.cov = spat.cov, temp.cov = temp.cov, type = "IF",
                                    start_methIF = start_methIF )
      cat("start values are", start_vals, "\n")
    }
    if(st_val_meth == "sp_gev"){
      start_params <- fitspatgev(data = data, covariables = as.matrix(spat.cov),
                                 loc.form =  scale.sp.form, scale.form = scale.sp.form,
                                 shape.form = ~ 1)$fitted.values

      st_param <- list(loc= start_params[startsWith(names(start_params),"loc")],
                       scale = start_params[startsWith(names(start_params),"scale")],
                       shape = start_params[startsWith(names(start_params),"shape")])


      model.scale.sp <- model.matrix(scale.sp.form, model.frame(scale.sp.form,
                                                                data = spat.cov, na.action = na.pass))

      sigma.mat <- model.scale.sp %*% st_param$scale
      loc.mat <- model.scale.sp %*% st_param$loc
      disp_st <- mean(loc.mat/sigma.mat)

      start_vals <- numeric(n.param)
      start_vals[1] <- disp_st
      start_vals[2:(n.param -1)] <- log(st_param$scale)
      start_vals[n.param] <- st_param$shape
      names(start_vals ) <- c("disp",
                              paste0("scale", 0:(n.scale.sp -1)), "shape")
      cat("start values are ", start_vals , "\n")
    }
    if(st_val_meth == "fgev"){

      start_st <-  apply(data,2, evd::fgev,  std.err = F)
      start_st <- purrr::map(purrr::map_dfr(start_st, function(x) x$estimate), mean)
      start_vals <- numeric(n.param)
      start_vals[1] <- start_st$loc/start_st$scale
      start_vals[2:(n.param -1)] <- log(start_st$scale)
      start_vals[n.param] <- start_st$shape
      names(start_vals ) <- c("disp",
                              paste0("scale", 0:(n.scale.sp -1)), "shape")
      cat("start values are ", start_vals , "\n")
    }
  }
  if(use_gr){
    mlest <- optim(start_vals, fn = nll_spat_IF, gr = gr_spat_gev_IF,

                   scale.sp.form = scale.sp.form,
                   scale.temp.form = scale.temp.form,
                   data = data, spat.cov = spat.cov, temp.cov = temp.cov,
                   method = method)
  }
  else{
    mlest <- optim(start_vals, fn = nll_spat_IF,

                   scale.sp.form = scale.sp.form,
                   scale.temp.form = scale.temp.form,
                   data = data, spat.cov = spat.cov, temp.cov = temp.cov,
                   method = method)
  }
  return(list(mle = mlest$par, nllh = mlest$value, conv = mlest$convergence))
}


# data is tibble of unique pairs of sliding bm and temporal covariate with
# respective frequency of appearance, as computed by get_unique_bm()

## use st_val_meth = sp_gev only when scale doesnt depend on covariate, since then, start values
## are computed with SpatialExtremes package that does not use a link function

## ... can be the following arguments: maxit or similar for optim,
## temp.cov for computing starting values with SpatialExtremes when depending on temporal cov.
#  type = IF for homogeneity assumption
fit_spgev_sl <- function(data, loc.sp.form, scale.sp.form,
                      loc.temp.form = NULL, scale.temp.form = NULL,
                      spat.cov, start_vals = NULL, datastart = NULL, use_gr = FALSE,
                      method = "BFGS", st_val_meth = "LeastSq", print_start_vals = TRUE,
                      scale.link = make.link("log"), ...){

  add.args <- list(...)

  d.dat <- nrow(data)
  n.spat <- nrow(spat.cov)

  if (d.dat != n.spat) {
    stop("'data' and 'covariates' doesn't match")
  }

  if (is.null(start_vals) & is.null(datastart)) {
    stop("Please provide either starting values or a dataset where starting values can
         be computed on.")
  }

  dataprep <- prep4spatml(loc.sp.form = loc.sp.form, scale.sp.form = scale.sp.form,
                          loc.temp.form = loc.temp.form, scale.temp.form = scale.temp.form,
                          data = data, spat.cov = spat.cov)


  if(is.null(start_vals)){
    if(!exists("temp.cov", where = add.args)) { add.args$temp.cov <- NULL }

    start_vals <-  get_start_vals(data = datastart,
                                  loc.sp.form = loc.sp.form,
                                  scale.sp.form = scale.sp.form,
                                  loc.temp.form = loc.temp.form,
                                  scale.temp.form = scale.temp.form,
                                  spat.cov = spat.cov,
                                  temp.cov = add.args$temp.cov,
                                  method = st_val_meth,
                                  print_start_vals = print_start_vals,
                                  scale.link = scale.link,
                                  type = add.args$type)

    add.args <- add.args[!(names(add.args) == "temp.cov")]

  }


  if(use_gr){
    mlest <- optim(start_vals, fn = nll_spat, gr =  gr_spat_gev,
                   loc.sp.form = loc.sp.form,
                   scale.sp.form = scale.sp.form,
                   loc.temp.form = loc.temp.form, scale.temp.form = scale.temp.form,
                   data = data, spat.cov = spat.cov, temp.cov = temp.cov,
                   method = method, control = add.args)
  }
  if(exists("type", where = add.args)) {

    add.args <- add.args[!(names(add.args) == "type")]
    mlest <- optim(start_vals, fn = nll_spat_temp_sl_hom,
                   scale.sp.form = scale.sp.form,
                   scale.temp.form = scale.temp.form,
                   dataprep = dataprep, spat.cov = spat.cov,
                   scale.link = scale.link,
                   method = method, control = add.args)
  }
  else{
    mlest <- optim(start_vals, fn = nll_spat_temp_sl_prep,
                   loc.sp.form = loc.sp.form,
                   scale.sp.form = scale.sp.form,
                   loc.temp.form = loc.temp.form, scale.temp.form = scale.temp.form,
                   dataprep = dataprep, spat.cov = spat.cov,
                   scale.link = scale.link,
                   method = method, control = add.args)
  }

  if(!(mlest$convergence == 0)) {
    warning("Optimization did not succeed. Try increasing the maximum number of iterations or try a different method for finding starting values.")
  }

  list(mle = mlest$par, nllh = mlest$value, conv = mlest$convergence, counts = mlest$counts)

}




# fit_spgev_sl(data = bms, loc.sp.form = ~ lon + ele, scale.sp.form =  ~ lon,
#              loc.temp.form = ~ GMST, scale.temp.form = NULL,
#              spat.cov = spatial_cvrt, start_vals = params, use_gr = FALSE,
#              method = "Nelder-Mead", st_val_meth = "LeastSq", print_start_vals = TRUE)
