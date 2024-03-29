% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_mda_data.R
\name{MDAdata}
\alias{MDAdata}
\title{Simulate spatially dependent data with margins in MDA(GEV)}
\usage{
MDAdata(
  n,
  locations,
  margins = list(distr = "gpd", shape = 0.2),
  ms = list(kovmod = "whitmat", nugget = 0, range = 3, smooth = 1)
)
}
\arguments{
\item{n}{Number of observations at each location}

\item{locations}{Matrix with named columns giving the coordinates of the locations
(lon and lat)}

\item{margins}{A list with elements
\describe{
\item{dist}{The marginal distribution of the resulting data. Must be one of
"pareto", "gpd", "norm", "gev".}
\item{shape}{The shape parameter of the chosen marginal distribution.}
}}

\item{ms}{A list with components that are passed to \code{rmaxstab} from the
\code{SpatialExtremes}-package:
\describe{
\item{kovmod}{The covariance model for the max-stable process.}
\item{nugget, range, smooth}{The nugget, range and smooth parameters.}
}}
}
\value{
Returns a matrix with numeric values of dimension \eqn{ n x d},
where \eqn{d} is the number of
locations passed in \code{coords}. The dependence structure of the returned data is
that  of the chosen max-stable model, while the marginal distributions are as
chosen.
}
\description{
Simulate spatially dependent data with margins in MDA(GEV)
}
\examples{
coords <- cbind(runif(2), runif(2))
colnames(coords ) <- c("lat", "lon")
mdax <- MDAdata(100, locations = coords, margins = list(distr = "gpd", shape = -0.2),
    ms = list(kovmod = "whitmat", nugget = 0, range = 3, smooth = 0.44))
plot.ts(mdax[, 1])
lines(mdax[, 2], col = 2)

}
