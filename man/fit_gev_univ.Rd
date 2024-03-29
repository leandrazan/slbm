% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_gev_univ.R
\name{fit_gev_univ}
\alias{fit_gev_univ}
\title{GEV fit with trend}
\usage{
fit_gev_univ(
  data,
  method = "BFGS",
  maxiter = 100,
  reltol = 1e-08,
  hessian = FALSE,
  type,
  return_cov = FALSE,
  ...
)
}
\arguments{
\item{data}{A tibble containing values of the unique sliding BM along with the
corresponding value of the temporal covariate and the frequency of occurence of
the respective tupel. Can be obtained by applying \code{\link[slbm]{get_uniq_bm}}.}

\item{method}{The method used during optimisation; passed to optim.}

\item{maxiter}{Passed to optim.}

\item{hessian}{logical; whether to return the hessian matrix.}

\item{type}{One of   'stationary', 'scale' or 'shift', see details.}

\item{return_cov}{logical; whether or not to return an estimate of the covariance matrix. If TRUE,
some further parameters/data need to be passed, see '...'}

\item{...}{Several additional arguments which are only needed in some cases, e.g. in the scale model.
See \code{details} for further information.}
}
\value{
A list containing the parameter estimates, the value of the negative log-Likelihood,
the convergence code (ouput from \code{optim()}, 0 means everything was ok)
and the hessian (if \code{hessian = TRUE}).
}
\description{
Fit a GEV model that either shifts or scales with a temporal covariate
to univariate data
}
\section{Details on type}{
The argument given in type determines whether the observations shift or scale
with time. For a temporal covariate \eqn{(c_t)_t}, shifting corresponds to a shift
in the location parameter as follows:
\deqn{ \mu(t) = \mu + \alpha c_t , \sigma(t) = \sigma, \gamma(t) = \gamma}
while scaling corresponds to the model where
\deqn{ \mu(t) = \mu \exp(\alpha c_t /\mu),  \sigma(t) = \sigma  \exp(\alpha c_t /\mu),
\gamma(t) = \gamma,
}
as inspired by the Clausius-Clapeyron relation.
The latter model can also be parametrised as
\deqn{ \mu(t) = \mu \exp(\alpha c_t),  \sigma(t) = \sigma  \exp(\alpha c_t),
\gamma(t) = \gamma,
}
which is why one needs to specify the argument \code{rel_trend} whenever \code{type = scale}
}

\section{Detials on additional arguments}{
Additional arguments that need to be passed in some cases are
\itemize{
\item rel_trend : logical; specifies the parametrisation of the scale model,
i.e. only relevant when when 'type = scale'. When \code{TRUE}, the trend parameter \eqn{\alpha} is
seen relative to the location parameter \eqn{\mu}.
}

When an estimate of the covariance matrix is required (\code{return_cov = TRUE}),
one needs to pass the following arguments and data:
\itemize{
\item chain: logical; whether to use covariance matrix estimation based on chain rule
\item varmeth: one of 'V', 'V2', 'both' (see documentation of \code{\link{est_var_univ}} for further details)
\item blcksz The blocksize parameter.
\item orig_slbm: the original sample of sliding block maxima (as generated by \code{\link{blockmax}})
\item orig_cvrt: the temporal covariate for the original sliding block maxima.
}
}

\examples{
##### generate some data #####
set.seed(1)
blcksz <- 90
xx <- evd::rgpd(100*90, shape = 0.2)

# define a temporal covariate that is constant over a block of length blcksz
temp_cvrt <- rep(1:100/100, each = blcksz)[1:(99*blcksz + 1)]

bms <- get_uniq_bm(xx, blcksz, temp_cvrt = temp_cvrt, looplastblock = FALSE)
bms
###############################
fit_gev_univ(data = bms, hessian = TRUE, type = "shift")
### shift-fit with covariance matrix estimation:
fit_gev_univ(data = bms, hessian = TRUE, type = "shift", return_cov = TRUE,
varmeth = "V2", chain = TRUE, orig_slbm = blockmax(xx, 90, "sliding"),
orig_cvrt = temp_cvrt, blcksz = 90)
}
