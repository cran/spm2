\name{sponge2}
\alias{sponge2}
\docType{data}
\title{A dataset of sponge species richness in the Timor Sea region, northern Australia marine margin}
\description{This dataset contains 77 samples of 81 variables including easting (longitude), northing (latitude), bathy, backscatter and their derived variables.}
\usage{data("sponge2")}
\format{
  A data frame with 77 observations on the following 89 variables.
  \describe{
    \item{\code{easting}}{a numeric vector, m}
    \item{\code{northing}}{a numeric vector, m}
    \item{\code{species.richness}}{a numeric vector, no unit}
    \item{\code{mud}}{a numeric vector, percentage}
    \item{\code{sand}}{a numeric vector, percentage}
    \item{\code{gravel}}{a numeric vector, percentage}
    \item{\code{bathy}}{a numeric vector, m}
    \item{\code{bs25}}{a numeric vector, dB}
    \item{\code{bs10}}{a numeric vector, dB}
    \item{\code{bs11}}{a numeric vector, dB}
    \item{\code{bs12}}{a numeric vector, dB}
    \item{\code{bs13}}{a numeric vector, dB}
    \item{\code{bs14}}{a numeric vector, dB}
    \item{\code{bs15}}{a numeric vector, dB}
    \item{\code{bs16}}{a numeric vector, dB}
    \item{\code{bs17}}{a numeric vector, dB}
    \item{\code{bs18}}{a numeric vector, dB}
    \item{\code{bs19}}{a numeric vector, dB}
    \item{\code{bs20}}{a numeric vector, dB}
    \item{\code{bs21}}{a numeric vector, dB}
    \item{\code{bs22}}{a numeric vector, dB}
    \item{\code{bs23}}{a numeric vector, dB}
    \item{\code{bs24}}{a numeric vector, dB}
    \item{\code{bs26}}{a numeric vector, dB}
    \item{\code{bs27}}{a numeric vector, dB}
    \item{\code{bs28}}{a numeric vector, dB}
    \item{\code{bs29}}{a numeric vector, dB}
    \item{\code{bs30}}{a numeric vector, dB}
    \item{\code{bs31}}{a numeric vector, dB}
    \item{\code{bs32}}{a numeric vector, dB}
    \item{\code{bs33}}{a numeric vector, dB}
    \item{\code{bs34}}{a numeric vector, dB}
    \item{\code{bs35}}{a numeric vector, dB}
    \item{\code{bs36}}{a numeric vector, dB}
    \item{\code{bs_o}}{a numeric vector, dB}
    \item{\code{bs_homo_o}}{a numeric vector}
    \item{\code{bs_entro_o}}{a numeric vector, no unit}
    \item{\code{bs_var_o}}{a numeric vector, dB^2}
    \item{\code{bs_lmi_o}}{a numeric vector}
    \item{\code{bathy_o}}{a numeric vector, m}
    \item{\code{bathy_lmi_o}}{a numeric vector}
    \item{\code{tpi_o}}{a numeric vector, no unit}
    \item{\code{slope_o}}{a numeric vector}
    \item{\code{plan_cur_o}}{a numeric vector}
    \item{\code{prof_cur_o}}{a numeric vector}
    \item{\code{relief_o}}{a numeric vector}
    \item{\code{rugosity_o}}{a numeric vector}
    \item{\code{dist.coast}}{a numeric vector, m}
    \item{\code{rugosity3}}{a numeric vector}
    \item{\code{rugosity5}}{a numeric vector}
    \item{\code{rugosity7}}{a numeric vector}
    \item{\code{tpi3}}{a numeric vector, no unit}
    \item{\code{tpi5}}{a numeric vector, no unit}
    \item{\code{tpi7}}{a numeric vector, no unit}
    \item{\code{bathy_lmi3}}{a numeric vector}
    \item{\code{bathy_lmi5}}{a numeric vector}
    \item{\code{bathy_lmi7}}{a numeric vector}
    \item{\code{plan_curv3}}{a numeric vector}
    \item{\code{plan_curv5}}{a numeric vector}
    \item{\code{plan_curv7}}{a numeric vector}
    \item{\code{relief_3}}{a numeric vector}
    \item{\code{relief_5}}{a numeric vector}
    \item{\code{relief_7}}{a numeric vector}
    \item{\code{slope3}}{a numeric vector}
    \item{\code{slope5}}{a numeric vector}
    \item{\code{slope7}}{a numeric vector}
    \item{\code{prof_cur3}}{a numeric vector}
    \item{\code{prof_cur5}}{a numeric vector}
    \item{\code{prof_cur7}}{a numeric vector}
    \item{\code{entro3}}{a numeric vector, no unit}
    \item{\code{entro5}}{a numeric vector, no unit}
    \item{\code{entro7}}{a numeric vector, no unit}
    \item{\code{homo3}}{a numeric vector}
    \item{\code{homo5}}{a numeric vector}
    \item{\code{homo7}}{a numeric vector}
    \item{\code{var3}}{a numeric vector, dB^2}
    \item{\code{var5}}{a numeric vector, dB^2}
    \item{\code{var7}}{a numeric vector, dB^2}
    \item{\code{bs_lmi3}}{a numeric vector}
    \item{\code{bs_lmi5}}{a numeric vector}
    \item{\code{bs_lmi7}}{a numeric vector}
  }
}
\details{For details, please see the source. This dataset was published as an appendix of the paper listed in the source. Where the long and lat were reprojected to easting and northing.}
\source{see Appendix A-D. Supplementary data at: "http://dx.doi.org/10.1016/j.envsoft.2017.07.016."}
\references{
Li, J., B. Alvarez, J. Siwabessy, M. Tran, Z. Huang, R. Przeslawski, L. Radke, F. Howard, and S. Nichol. 2017. Application of random forest, generalised linear model and their hybrid methods with geostatistical techniques to count data: Predicting sponge species richness. Environmental Modelling & Software, 97: 112-129.
}
\keyword{datasets}
