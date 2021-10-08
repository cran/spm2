\name{bees}
\alias{bees}
\docType{data}
\title{A dataset of bees count data and relevant information in oilseed Brassica fields in an Australian temperate landscape.
}
\description{
This dataset contains 212 samples of 61 variables including three bee species, inflorescence, temperature, wid speed and various derived landscape variables.
}
\usage{data("bees")}
\format{
  A data frame with 212 observations on the following 61 variables.
  \describe{
    \item{\code{transid}}{a factor with levels \code{G1} \code{G10} \code{G11} \code{G12} \code{G13} \code{G14} \code{G15} \code{G16} \code{G1D1} \code{G1D10} \code{G1D11} \code{G1D12} \code{G1D13} \code{G1D14} \code{G1D15} \code{G1D16} \code{G1D17} \code{G1D2} \code{G1D3} \code{G1D4} \code{G1D5} \code{G1D6} \code{G1D7} \code{G1D8} \code{G1D9} \code{G1G1} \code{G1G10} \code{G1G11} \code{G1G12} \code{G1G13} \code{G1G14} \code{G1G15} \code{G1G16} \code{G1G17} \code{G1G2} \code{G1G3} \code{G1G4} \code{G1G5} \code{G1G6} \code{G1G7} \code{G1G8} \code{G1G9} \code{G2} \code{G2H1} \code{G2H10} \code{G2H11} \code{G2H2} \code{G2H3} \code{G2H4} \code{G2H5} \code{G2H6} \code{G2H7} \code{G2H8} \code{G2H9} \code{G2S1} \code{G2S10} \code{G2S11} \code{G2S2} \code{G2S3} \code{G2S4} \code{G2S5} \code{G2S6} \code{G2S7} \code{G2S8} \code{G2S9} \code{G3} \code{G4} \code{G5} \code{G6} \code{G7} \code{G8} \code{G9} \code{GCH1} \code{GCH10} \code{GCH11} \code{GCH12} \code{GCH13} \code{GCH14} \code{GCH15} \code{GCH16} \code{GCH17} \code{GCH18} \code{GCH2} \code{GCH3} \code{GCH4} \code{GCH5} \code{GCH6} \code{GCH7} \code{GCH8} \code{GCH9} \code{GCS1} \code{GCS10} \code{GCS11} \code{GCS12} \code{GCS13} \code{GCS14} \code{GCS15} \code{GCS16} \code{GCS17} \code{GCS18} \code{GCS2} \code{GCS3} \code{GCS4} \code{GCS5} \code{GCS6} \code{GCS7} \code{GCS8} \code{GCS9} \code{H1} \code{H10} \code{H11} \code{H12} \code{H13} \code{H14} \code{H15} \code{H16} \code{H2} \code{H3} \code{H4} \code{H5} \code{H6} \code{H7} \code{H8} \code{H9} \code{MC1G1} \code{MC1G2} \code{MC1G3} \code{MC1G4} \code{MC1G5} \code{MC1G6} \code{MC1G7} \code{MC1G8} \code{MC1G9} \code{MC1H1} \code{MC1H2} \code{MC1H3} \code{MC1H4} \code{MC1H5} \code{MC1H6} \code{MC1H7} \code{MC1H8} \code{MC1H9} \code{MC2AH1} \code{MC2AH10} \code{MC2AH11} \code{MC2AH12} \code{MC2AH13} \code{MC2AH14} \code{MC2AH15} \code{MC2AH16} \code{MC2AH17} \code{MC2AH18} \code{MC2AH2} \code{MC2AH3} \code{MC2AH4} \code{MC2AH5} \code{MC2AH6} \code{MC2AH7} \code{MC2AH8} \code{MC2AH9} \code{MC2AS1} \code{MC2AS10} \code{MC2AS11} \code{MC2AS12} \code{MC2AS13} \code{MC2AS14} \code{MC2AS15} \code{MC2AS16} \code{MC2AS17} \code{MC2AS18} \code{MC2AS2} \code{MC2AS3} \code{MC2AS4} \code{MC2AS5} \code{MC2AS6} \code{MC2AS7} \code{MC2AS8} \code{MC2AS9} \code{MC2BB1} \code{MC2BB10} \code{MC2BB11} \code{MC2BB12} \code{MC2BB13} \code{MC2BB14} \code{MC2BB15} \code{MC2BB16} \code{MC2BB17} \code{MC2BB2} \code{MC2BB3} \code{MC2BB4} \code{MC2BB5} \code{MC2BB6} \code{MC2BB7} \code{MC2BB8} \code{MC2BB9} \code{MC2BG1} \code{MC2BG10} \code{MC2BG11} \code{MC2BG12} \code{MC2BG13} \code{MC2BG14} \code{MC2BG15} \code{MC2BG16} \code{MC2BG17} \code{MC2BG2} \code{MC2BG3} \code{MC2BG4} \code{MC2BG5} \code{MC2BG6} \code{MC2BG7} \code{MC2BG8} \code{MC2BG9}}
    \item{\code{transsurv}}{a numeric vector}
    \item{\code{plotsurv}}{a numeric vector}
    \item{\code{paddock}}{a numeric vector}
    \item{\code{plot}}{a factor with levels \code{G1-1} \code{G1-10} \code{G1-11} \code{G1-12} \code{G1-13} \code{G1-14} \code{G1-15} \code{G1-16} \code{G1-17} \code{G1-2} \code{G1-3} \code{G1-4} \code{G1-5} \code{G1-6} \code{G1-7} \code{G1-8} \code{G1-9} \code{G21} \code{G210} \code{G211} \code{G22} \code{G23} \code{G24} \code{G25} \code{G26} \code{G27} \code{G28} \code{G29} \code{GC1} \code{GC10} \code{GC11} \code{GC12} \code{GC13} \code{GC14} \code{GC15} \code{GC16} \code{GC17} \code{GC18} \code{GC2} \code{GC3} \code{GC4} \code{GC5} \code{GC6} \code{GC7} \code{GC8} \code{GC9} \code{MC-1} \code{MC-10} \code{MC-11} \code{MC-12} \code{MC-13} \code{MC-14} \code{MC-15} \code{MC-16} \code{MC-2} \code{MC-3} \code{MC-4} \code{MC-5} \code{MC-6} \code{MC-7} \code{MC-8} \code{MC-9} \code{MC1-1} \code{MC1-2} \code{MC1-3} \code{MC1-4} \code{MC1-5} \code{MC1-6} \code{MC1-7} \code{MC1-8} \code{MC1-9} \code{MC2-A1} \code{MC2-A10} \code{MC2-A11} \code{MC2-A12} \code{MC2-A13} \code{MC2-A14} \code{MC2-A15} \code{MC2-A16} \code{MC2-A17} \code{MC2-A18} \code{MC2-A2} \code{MC2-A3} \code{MC2-A4} \code{MC2-A5} \code{MC2-A6} \code{MC2-A7} \code{MC2-A8} \code{MC2-A9} \code{MC2-B1} \code{MC2-B10} \code{MC2-B11} \code{MC2-B12} \code{MC2-B13} \code{MC2-B14} \code{MC2-B15} \code{MC2-B16} \code{MC2-B17} \code{MC2-B2} \code{MC2-B3} \code{MC2-B4} \code{MC2-B5} \code{MC2-B6} \code{MC2-B7} \code{MC2-B8} \code{MC2-B9}}
    \item{\code{obs}}{a factor with levels \code{Andrew} \code{Barbara} \code{Micah} \code{Sonia} \code{Steve}}
    \item{\code{hbee}}{a numeric vector}
    \item{\code{nbee}}{a numeric vector}
    \item{\code{hover}}{a numeric vector}
    \item{\code{date}}{a numeric vector}
    \item{\code{sx}}{a numeric vector}
    \item{\code{fx}}{a numeric vector}
    \item{\code{sy}}{a numeric vector}
    \item{\code{fy}}{a numeric vector}
    \item{\code{loc}}{a factor with levels \code{100} \code{150} \code{200} \code{400} \code{450} \code{edge}}
    \item{\code{pair}}{a factor with levels \code{I} \code{O}}
    \item{\code{inf}}{a numeric vector}
    \item{\code{rankinf}}{a numeric vector}
    \item{\code{dupl}}{a numeric vector}
    \item{\code{temp}}{a numeric vector}
    \item{\code{windspeed}}{a numeric vector}
    \item{\code{winddir}}{a factor with levels \code{N} \code{NE} \code{NNE} \code{NNW} \code{NW}}
    \item{\code{cloudc}}{a numeric vector}
    \item{\code{disttoedgecalc}}{a numeric vector}
    \item{\code{disttoedgemeasured}}{a numeric vector}
    \item{\code{w100}}{a numeric vector}
    \item{\code{w200}}{a numeric vector}
    \item{\code{w300}}{a numeric vector}
    \item{\code{w400}}{a numeric vector}
    \item{\code{w500}}{a numeric vector}
    \item{\code{w600}}{a numeric vector}
    \item{\code{w700}}{a numeric vector}
    \item{\code{w800}}{a numeric vector}
    \item{\code{w900}}{a numeric vector}
    \item{\code{w1000}}{a numeric vector}
    \item{\code{w1500}}{a numeric vector}
    \item{\code{w2000}}{a numeric vector}
    \item{\code{c100}}{a numeric vector}
    \item{\code{c200}}{a numeric vector}
    \item{\code{c300}}{a numeric vector}
    \item{\code{c400}}{a numeric vector}
    \item{\code{c500}}{a numeric vector}
    \item{\code{c1000}}{a numeric vector}
    \item{\code{c1500}}{a numeric vector}
    \item{\code{area}}{a numeric vector}
    \item{\code{perimeter}}{a numeric vector}
    \item{\code{gyration}}{a numeric vector}
    \item{\code{paratio}}{a numeric vector}
    \item{\code{shape}}{a numeric vector}
    \item{\code{fractaldimention}}{a numeric vector}
    \item{\code{circumscircle}}{a numeric vector}
    \item{\code{contiguity}}{a numeric vector}
    \item{\code{links100}}{a numeric vector}
    \item{\code{links200}}{a numeric vector}
    \item{\code{links300}}{a numeric vector}
    \item{\code{links400}}{a numeric vector}
    \item{\code{links500}}{a numeric vector}
    \item{\code{links1000}}{a numeric vector}
    \item{\code{links1500}}{a numeric vector}
    \item{\code{links2000}}{a numeric vector}
    \item{\code{windspeed2}}{a numeric vector}
  }
}
\details{
For details, please see the source. This dataset was published as an appendix of the paper listed in the source. Where the long and lat were reprojected to easting and northing
}
\source{
The data source is [https://doi.org/10.25919/5f17b34638cca] or [https://data.csiro.au/collections/collection/CIcsiro:45533], which provides bees count data and relevant predictive variables along with a brief description of the data. The detailed descriptions of the data are available in: "Arthur, A. D., et al. (2010). "Influence of woody vegetation on pollinator densities in oilseed Brassica fields in an Australian temperate landscape." Basic and Applied Ecology 11: 406-414."
}
\references{
Arthur, A. D., Li, J., Henry, S., Cunningham, S.A. (2010). "Influence of woody vegetation on pollinator densities in oilseed Brassica fields in an Australian temperate landscape." Basic and Applied Ecology 11: 406-414.
}
\keyword{datasets}
