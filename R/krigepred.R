#' @title Generate spatial predictions using kriging methods ('krige')
#'
#' @description This function is to make spatial predictions using kriging methods ('krige').
#'
#' @param trainx a dataframe contains longitude (long), latitude (lat) and
#' predictive variables of point samples. The location information must be named
#' as 'long' and 'lat'.
#' @param trainy a vector of response, must have length equal to the number of
#' rows in trainx.
#' @param trainx2 a dataframe contains longitude (long), latitude (lat) and
#' predictive variables of point locations (i.e., the centres of grids) to
#' be predicted. The location information must be named as 'long' and 'lat' and
#' in the first two columns respectively..
#' @param nmax for local kriging: the number of nearest observations that
#' should be used for a kriging prediction or simulation, where nearest is
#' defined in terms of the space of the spatial locations. By default, 12
#' observations are used.
#' @param transformation transform the response variable to normalise the data;
#' can be "sqrt" for square root, "arcsine" for arcsine, "log" or "none"
#' for non transformation. By default, "none" is used.
#' @param delta numeric; to avoid log(0) in the log transformation.
#' @param formula formula defining the response vector and (possible) regressor.
#' an object (i.e., 'variogram.formula') for 'variogram' or a formula for
#' 'krige'. see 'variogram' and 'krige' in gstat for details.
#' @param vgm.args arguments for vgm, e.g. variogram model of response
#' variable and anisotropy parameters. see notes vgm in gstat for details.
#' By default, "Sph" is used.
#' @param anis anisotropy parameters: see notes vgm in gstat for details.
#' @param alpha direction in plane (x,y). see variogram in gstat for details.
#' @param block block size. see krige in gstat for details.
#' @param beta for simple kriging. see krige in gstat for details.
#' @param ... other arguments passed on to gstat.
#'
#' @return A dataframe of longitude, latitude, predictions and variances.
#'
#' @note The variances in the output are not transformed back when a transformation
#' is used. This is because kriging variances are not measuring the uncertainty of
#' predictions but they are indicator of the spatial distribution of sample density.
#' The variances are exported only for interested users; and if needed,
#' they can be transformed back from the output.
#'
#' @references Pebesma, E.J., 2004. Multivariable geostatistics in S: the gstat
#' package. Computers & Geosciences, 30: 683-691.
#'
#' @author Jin Li
#' @examples
#' \donttest{
#' library(sp)
#' library(spm)
#' data(swmud)
#' data(sw)
#' okpred1 <- krigepred(swmud[, c(1,2)], swmud[, 3], sw, nmax = 7, transformation =
#' "arcsine", vgm.args = ("Sph"))
#' names(okpred1)
#' }
#'
#' @export
krigepred <- function (trainx, trainy, trainx2, nmax = 12, transformation =
  "none", delta = 1, formula = var1 ~ 1, vgm.args = ("Sph"), anis = c(0, 1), alpha = 0, block = 0,  beta, ...) {
  if (transformation == "none") {trainy1 = trainy} else (
  if (transformation == "sqrt") {trainy1 = sqrt(trainy)} else (
  if (transformation == "arcsine") {trainy1 = asin(sqrt(trainy / 100))} else (
  if (transformation == "log") {trainy1 = log(trainy + delta)} else (
  stop ("This transfromation is not supported in this version!")))))

  data.dev <- trainx
  data.pred <- trainx2
  data.dev$var1 <- trainy1
  sp::coordinates(data.dev) = ~ long + lat
  vgm1 <- gstat::variogram(object = formula, data.dev, alpha = alpha)
  # model.1 <- gstat::fit.variogram(vgm1, gstat::vgm(vgm.args, anis = anis))
  model.1 <- gstat::fit.variogram(vgm1, gstat::vgm(mean(vgm1$gamma), vgm.args, mean(vgm1$dist), min(vgm1$gamma)/10, anis = anis))

  if (model.1$range[2] <= 0) {model.1$range[2] <- min(vgm1$dist)} # a negative range is not allowed. To keep it working, we replace it with a minimum dist.

  sp::coordinates(data.pred) = ~ long + lat
  ok.pred1 <- gstat::krige(formula = formula, data.dev, data.pred,
    model = model.1, nmax = nmax, block = block, beta = beta)
  if (transformation == "none") {ok.pred2 = ok.pred1$var1.pred}
  if (transformation == "sqrt") {ok.pred2 = ok.pred1$var1.pred ^ 2}
  if (transformation == "arcsine") {ok.pred2 = (sin(ok.pred1$var1.pred)) ^ 2 * 100}
  if (transformation == "log") {ok.pred2 = exp(ok.pred1$var1.pred)-delta}
  var1.pred = ok.pred2
  var1.var = ok.pred1$var1.var
  ok.pred <- cbind(trainx2[, c(1:2)], var1.pred, var1.var)
  ok.pred
}
