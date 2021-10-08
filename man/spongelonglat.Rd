\name{spongelonglat}
\alias{spongelonglat}
\docType{data}
\title{A dataset of sponge species richness in the Timor Sea region, northern Australia marine margin}
\description{
This dataset contains 77 samples of 7 predictive variables including longitude, latitude, bathy, backscatter and their derived variables. It is the sponge dataset in `spm` package, but with long and lat instead of easting and northing.}
\usage{data("spongelonglat")}
\format{
  A data frame with 77 observations on the following 8 variables.
  \describe{
    \item{\code{sponge}}{a numeric vector}
    \item{\code{tpi3}}{a numeric vector}
    \item{\code{var7}}{a numeric vector}
    \item{\code{entro7}}{a numeric vector}
    \item{\code{bs34}}{a numeric vector}
    \item{\code{bs11}}{a numeric vector}
    \item{\code{long}}{a numeric vector}
    \item{\code{lat}}{a numeric vector}
  }
}
\details{
For details, please see sponge dataset in library(spm). Where the long and lat were projected to easting and northing.}
\source{sponge dataset in library(spm)}
\references{
Li, J., B. Alvarez, J. Siwabessy, M. Tran, Z. Huang, R. Przeslawski, L. Radke, F. Howard, and S. Nichol. 2017. Application of random forest, generalised linear model and their hybrid methods with geostatistical techniques to count data: Predicting sponge species richness. Environmental Modelling & Software, 97: 112-129.}
\examples{
data(spongelonglat)
## maybe str(spongelonglat) ; plot(spongelonglat) ...
}
\keyword{datasets}
