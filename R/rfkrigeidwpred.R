#' @title Generate spatial predictions using the hybrid methods of 'random forest' ('RF'),
#' 'kriging' and inverse distance weighted ('IDW').
#'
#' @description This function is for generating spatial predictions using the
#' hybrid methods of 'RF', 'kriging' and 'IDW', including all methods implemented
#' in 'rfkrigeidwcv'.
#'
#' @param longlat	a dataframe contains longitude and latitude of point samples.
#' @param trainx a dataframe contains longitude (long), latitude (lat),
#' predictive variables and the response variable of point samples. That is,
#' the location information must be named as 'long' and 'lat'.
#' @param predx	a dataframe or matrix contains columns of predictive variables
#' for the grids to be predicted.
#' @param trainy a vector of the response variable in the formula, that is, the left
#' part of the formula.
#' @param longlatpredx	a dataframe contains longitude and latitude of point locations
#' (i.e., the centers of grids) to be predicted.
#' @param mtry a function of number of remaining predictor variables to use as
#' the 'mtry' parameter in the 'randomForest' call.
#' @param ntree number of trees to grow. This should not be set to too small a
#' number, to ensure that every input row gets predicted at least a few times.
#' By default, 500 is used.
#' #' @param longlatpredx	a dataframe contains longitude and latitude of point locations
#' (i.e., the centers of grids) to be predicted. The location information must be
#' named as 'long' and 'lat'.
#' @param transformation transform the residuals of 'rf' to normalise the data;
#' can be "sqrt" for square root, "arcsine" for arcsine, "log" or "none"
#' for non transformation. By default, "none" is used.
#' @param delta numeric; to avoid log(0) in the log transformation. The default is 1.
#' @param formula.krige formula defining the response vector and (possible) regressor.
#' an object (i.e., 'variogram.formula') for 'variogram' or a formula for
#' 'krige'. see 'variogram' and 'krige' in 'gstat' for details.
#' @param vgm.args arguments for 'vgm', e.g. variogram model of response
#' variable and anisotropy parameters. see 'vgm' in 'gstat' for details.
#' By default, "Sph" is used.
#' @param anis anisotropy parameters: see notes 'vgm' in 'gstat' for details.
#' @param alpha direction in plane (x,y). see variogram in 'gstat' for details.
#' @param block block size. see 'krige' in 'gstat' for details.
#' @param beta for simple kriging. see 'krige' in 'gstat' for details.
#' @param nmaxkrige for a local predicting: the number of nearest observations that
#'  should be used for a prediction or simulation, where nearest is defined in
#'  terms of the space of the spatial locations. By default, 12 observations
#'  are used.
#' @param idp	 a numeric number specifying the inverse distance weighting power.
#' @param nmaxidw for a local predicting: the number of nearest observations that
#'  should be used for a prediction or simulation, where nearest is defined in
#'  terms of the space of the spatial locations. By default, 12 observations
#'  are used.
#' @param hybrid.parameter the default is 2 that is for 'rfkrigerfidw';
#' for 'rfrfkrigerfidw', it needs to be 3.
#' @param lambda, ranging from 0 to 2; the default is 1 for 'rfkrigerfidw'
#' and 'rfrfkrigerfidw'; and if it is < 1, more weight is placed on 'krige',
#' otherwise more weight is placed on 'idw'; and if it is 0, 'idw' is not
#' considered and the resultant methods is 'rfkrige' when the default
#' 'hybrid.parameter' is used; and if it is 2, then the resultant method is
#' 'rfidw' when the default 'hybrid.parameter' is used.
#' @param ... other arguments passed on to 'rf', 'krige' and 'gstat'.
#'
#' @return A dataframe of longitude, latitude, and predictions.
#'
#' @references Li, J., Potter, A., Huang, Z., and Heap, A. (2012). Predicting Seabed
#' Sand Content across the Australian Margin Using Machine Learning and Geostatistical
#'  Methods, Geoscience Australia, Record 2012/48, 115pp.
#'
#' Li, J., Heap, A., Potter, A., and Danilel, J.J. (2011). Predicting Seabed Mud Content
#' across the Australian Margin II: Performance of Machine Learning Methods and Their
#' Combination with Ordinary Kriging and Inverse Distance Squared, Geoscience Australia,
#' Record 2011/07, 69pp.
#'
#' Liaw, A. and M. Wiener (2002). Classification and Regression by
#' randomForest. R News 2(3), 18-22.
#'
#' Pebesma, E.J., 2004. Multivariable geostatistics in S: the gstat package.
#' Computers & Geosciences, 30: 683-691.
#'
#' @author Jin Li
#' @examples
#' \donttest{
#' library(spm)
#'
#' data(sponge)
#' data(sponge.grid)
#' longlat <- sponge[, 1:2]
#' y = sponge[, 3]
#' trainx = sponge[, -3]
#'
#' set.seed(1234)
#'
#' rfkrigeidwpred1 <- rfkrigeidwpred(longlat = longlat, trainx =  trainx,
#' predx = sponge.grid, trainy = y, longlatpredx = sponge.grid[, c(1:2)],
#' formula.krige = res1 ~ 1, vgm.args = "Sph", nmaxkrige = 12, idp = 2, nmaxidw = 12)
#'
#' names(rfkrigeidwpred1)
#'
#' range(rfkrigeidwpred1$predictions)
#'}
#'
#' @export
rfkrigeidwpred <- function (longlat, trainx, predx, trainy, longlatpredx, mtry = function(p) max(1, floor(sqrt(p))), ntree = 500, transformation = "none", delta = 1, formula.krige = res1 ~ 1, vgm.args = c("Sph"), anis = c(0, 1), alpha = 0, block = 0, beta, nmaxkrige = 12, idp = 2, nmaxidw = 12, hybrid.parameter = 2, lambda = 1, ...) {

  n <- nrow(trainx)
  p <- ncol(trainx)
  names(longlat) <- c("long", "lat")
  names(longlatpredx) <- c("long", "lat")

  # rf modeling
  rf1 <- randomForest::randomForest(trainx, trainy, mtry = mtry(p), ntree=ntree)

  # rf predictions
  pred.rf1 <- stats::predict(rf1, predx)

  # the residuals of rf for krige
  data.dev1 <- longlat
  data.pred1 <- longlatpredx
  dev.rf1 <- stats::predict(rf1, trainx)

  res1 <- trainy - dev.rf1
  data.dev1$res1 <- res1

  # idw of the residuals
  gstat1 <- gstat::gstat(id = "res1", formula = res1 ~ 1, locations = ~ long + lat, data = data.dev1, set = list(idp = idp), nmax = nmaxidw)

  # idw predictions
  pred.idw1 <- stats::predict(gstat1, data.pred1)

  res1 <- trainy - dev.rf1
  if (transformation == "none") {data.dev1$res1 = res1} else (
    if (transformation == "sqrt") {data.dev1$res1 = sqrt(res1 + abs(min(res1)))} else (
      if (transformation == "arcsine") {data.dev1$res1 = asin(sqrt((res1 + abs(min(res1))) / 100))} else (
        if (transformation == "log") {data.dev1$res1 = log(res1 + abs(min(res1)) + delta)} else (
          stop ("This transfromation is not supported in this version!")))))
    # The '+ abs(min(res1))' above is to set possible negative values to 0.

  # vgm of the residuals
  sp::coordinates(data.dev1) = ~ long + lat
  vgm1 <- gstat::variogram(object = formula.krige, data.dev1, alpha = alpha)
  model.1 <- gstat::fit.variogram(vgm1, gstat::vgm(mean(vgm1$gamma), vgm.args, mean(vgm1$dist), min(vgm1$gamma)/10, anis = anis))
  if (model.1$range[2] <= 0) (cat("A zero or negative range was fitted to variogram", "\n"))
  if (model.1$range[2] <= 0) (model.1$range[2] <- min(vgm1$dist))  # set negative range to be positive

  # krige predictions
  sp::coordinates(data.pred1) = ~long + lat
  pred.krige1 <- gstat::krige(formula = formula.krige, data.dev1, data.pred1, model = model.1, nmax=nmaxkrige, block = block, beta = beta)$var1.pred

  if (transformation == "none") {pred.krige = pred.krige1}
  if (transformation == "sqrt") {pred.krige = pred.krige1 ^ 2 - abs(min(res1))}
  if (transformation == "arcsine") {pred.krige = (sin(pred.krige1)) ^ 2 * 100 -  abs(min(res1))}
  if (transformation == "log") {pred.krige = exp(pred.krige1) - abs(min(res1)) - delta}

  predictions <- (pred.krige * (2 - lambda) + pred.idw1$res1.pred * lambda + pred.rf1 * hybrid.parameter) / hybrid.parameter
  rfkrigeidw.pred <- cbind(longlatpredx, predictions)
  rfkrigeidw.pred
}

