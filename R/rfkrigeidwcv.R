#' @title Cross validation, n-fold and leave-one-out for the hybrid methods of
#' 'random forest' ('RF'), 'kriging' and inverse distance weighted ('IDW')
#'
#' @description This function is a cross validation function for 38 hybrid
#'  methods of 'RF', 'kriging' and 'IDW', including the average of 'rfkrige'
#' and 'rfidw' ('rfkrigerfidw') and  the average of 'rf', 'rfkrige' and 'rfidw'
#' ('rfrfkrigerfidw'), where 'kriging' methods include ordinary kriging ('OK'),
#'  simple kriging ('SK'), block 'OK' ('BOK') and block 'SK'('BSK') and 'IDW'
#'   also covers 'NN' and 'KNN'.. The data splitting is based on a stratified
#' random sampling method (see the 'datasplit' function for details).
#'
#' @param longlat	a dataframe contains longitude and latitude of point samples.
#' @param trainx a dataframe or matrix contains columns of predictive variables.
#' @param trainy a vector of the response variable.
#' @param mtry a function of number of remaining predictor variables to use as
#' the 'mtry' parameter in the 'randomForest' call.
#' @param ntree number of trees to grow. This should not be set to too small a
#' number, to ensure that every input row gets predicted at least a few times.
#' By default, 500 is used.
#' @param transformation transform the residuals of 'rf' to normalize the data for 'krige';
#' can be "sqrt" for square root, "arcsine" for arcsine, "log" or "none"
#' for non transformation. By default, "none" is used.
#' @param delta numeric; to avoid log(0) in the log transformation. The default is 1.
#' @param formula formula defining the response vector and (possible) regressor.
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
#' @param validation validation methods, include 'LOO': leave-one-out, and 'CV':
#' cross-validation.
#' @param cv.fold integer; number of folds in the cross-validation. if > 1,
#' then apply n-fold cross validation; the default is 10, i.e., 10-fold cross
#' validation that is recommended.
#' @param predacc can be either "VEcv" for vecv or "ALL" for all measures
#' in function pred.acc.
#' @param ... other arguments passed on to 'randomForest', 'krige' and 'gstat'.
#'
#' @return A list with the following components:
#'  me, rme, mae, rmae, mse, rmse, rrmse, vecv and e1; or vecv only
#' @note This function is largely based on 'rfcv' in 'randomForest', and 'krigecv'
#' in 'spm2'.
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
#' # rfrfokrfidw
#' data(sponge)
#' longlat <- sponge[, 1:2]
#' set.seed(1234)
#' rfrfkrigerfidwcv1 <- rfkrigeidwcv(longlat = longlat,
#' trainx = sponge[, -3], trainy = sponge[, 3], formula = res1 ~ 1, vgm.args = ("Sph"),
#' nmaxkrige = 12, idp = 2, nmaxidw = 12, hybrid.parameter = 3, validation = "CV",
#' predacc = "ALL")
#' rfrfkrigerfidwcv1
#'
#' # rfokrfidw for count data
#' data(sponge)
#' longlat <- sponge2[, 1:2]
#' y = sponge[, 3]
#' trainx = sponge[, -3]
#' set.seed(1234)
#' n <- 20 # number of iterations,60 to 100 is recommended.
#' VEcv <- NULL
#' for (i in 1:n) {
#'  rfkrigerfidwcv1 <- rfkrigeidwcv(longlat = longlat,
#'  trainx = trainx, trainy = y, formula = res1 ~ 1, vgm.args = ("Sph"),
#'  nmaxkrige = 12, idp = 2, nmaxidw = 12, validation = "CV",  predacc = "VEcv")
#'  VEcv [i] <- rfkrigerfidwcv1
#'  }
#'  plot(VEcv ~ c(1:n), xlab = "Iteration for rfokrfidw", ylab = "VEcv (%)")
#'  points(cumsum(VEcv) / c(1:n) ~ c(1:n), col = 2)
#'  abline(h = mean(VEcv), col = 'blue', lwd = 2)
#'}
#'
#' @export
rfkrigeidwcv <- function (longlat, trainx, trainy, mtry = function(p) max(1, floor(sqrt(p))), ntree = 500, transformation = "none", delta = 1, formula = res1 ~ 1, vgm.args = c("Sph"), anis = c(0, 1), alpha = 0, block = 0, beta, nmaxkrige = 12, idp = 2, nmaxidw = 12, hybrid.parameter = 2, lambda = 1, validation = "CV", cv.fold = 10, predacc = "VEcv", ...) {

  if (validation == "LOO") {idx <- 1:length(trainy)}
  if (validation == "CV")  {idx <- datasplit(trainy, k.fold = cv.fold)}

  names(longlat) <- c("long", "lat")

  # cross validation
  n <- nrow(trainx)
  p <- ncol(trainx)
  cv.pred <- NULL

  if (validation == "LOO") {
    for (i in 1 : length(trainy)) {
      data.dev <- trainx[idx != i, , drop = FALSE]
      data.pred <- trainx[idx == i, , drop = FALSE]

      # rf modelling
      rf1 <- randomForest::randomForest(data.dev, trainy[idx != i], mtry = mtry(p), ntree=ntree)
      # rf predictions
      pred.rf1 <- stats::predict(rf1, data.pred)


      # the residuals of rf for krige
      dev.rf1 <- stats::predict(rf1, data.dev)

      data.dev1 <- longlat[idx != i, , drop = FALSE] # for krige
      data.pred1 <- longlat[idx == i, , drop = FALSE] # for krige

      res1 <- trainy[idx != i] - dev.rf1

      data.dev1$res1 <- res1

      # idw of the residuals
      gstat1 <- gstat::gstat(id = "res1", formula = res1 ~ 1, locations = ~ long + lat, data = data.dev1, set = list(idp = idp), nmax = nmaxidw)

      # idw predictions
      pred.idw1 <- stats::predict(gstat1, data.pred1)

      # for krige
      if (transformation == "none") {data.dev1$res1 = res1} else (
        if (transformation == "sqrt") {data.dev1$res1 = sqrt(res1 + abs(min(res1)))} else (
          if (transformation == "arcsine") {data.dev1$res1 = asin(sqrt((res1 + abs(min(res1))) / 100))} else (
            if (transformation == "log") {data.dev1$res1 = log(res1 + abs(min(res1)) + delta)} else (
              stop ("This transfromation is not supported in this version!")))))
        # The '+ abs(min(res1))' above is to set possible negative values to 0.

      # vgm of the residuals
      sp::coordinates(data.dev1) = ~ long + lat
      vgm1 <- gstat::variogram(object = formula, data.dev1, alpha = alpha)
      model.1 <- gstat::fit.variogram(vgm1, gstat::vgm(mean(vgm1$gamma), vgm.args, mean(vgm1$dist), min(vgm1$gamma)/10, anis = anis))
      if (model.1$range[2] <= 0) (cat("A zero or negative range was fitted to variogram", "\n"))
      if (model.1$range[2] <= 0) (model.1$range[2] <- min(vgm1$dist))  # set negative range to be positive

      # krige predictions
      sp::coordinates(data.pred1) = ~long + lat
      pred.krige1 <- gstat::krige(formula = formula, data.dev1, data.pred1, model = model.1, nmax=nmaxkrige, block = block, beta = beta)$var1.pred

      if (transformation == "none") {pred.krige = pred.krige1}
      if (transformation == "sqrt") {pred.krige = pred.krige1 ^ 2 - abs(min(res1))}
      if (transformation == "arcsine") {pred.krige = (sin(pred.krige1)) ^ 2 * 100 -  abs(min(res1))}
      if (transformation == "log") {pred.krige = exp(pred.krige1) - abs(min(res1)) - delta}

      cv.pred[idx == i] <- (pred.krige * (2 - lambda) + pred.idw1$res1.pred * lambda + pred.rf1 * hybrid.parameter) / hybrid.parameter
    }
  }

  if (validation == "CV") {
  for (i in 1 : cv.fold) {
    data.dev <- trainx[idx != i, , drop = FALSE]
    data.pred <- trainx[idx == i, , drop = FALSE]

    # rf modelling
    rf1 <- randomForest::randomForest(data.dev, trainy[idx != i], mtry = mtry(p), ntree=ntree)

    # rf predictions
    pred.rf1 <- stats::predict(rf1, data.pred)

    # the residuals of rf for krige
    dev.rf1 <- stats::predict(rf1, data.dev)
    data.dev1 <- longlat[idx != i, , drop = FALSE] # for krige
    data.pred1 <- longlat[idx == i, , drop = FALSE] # for krige
    res1 <- trainy[idx != i] - dev.rf1

    data.dev1$res1 <- res1

    # idw of the residuals
    gstat1 <- gstat::gstat(id = "res1", formula = res1 ~ 1, locations = ~ long + lat, data = data.dev1, set = list(idp = idp), nmax = nmaxidw)

    # idw predictions
    pred.idw1<- stats::predict(gstat1, data.pred1)

    if (transformation == "none") {data.dev1$res1 = res1} else (
      if (transformation == "sqrt") {data.dev1$res1 = sqrt(res1 + abs(min(res1)))} else (
        if (transformation == "arcsine") {data.dev1$res1 = asin(sqrt((res1 + abs(min(res1))) / 100))} else (
          if (transformation == "log") {data.dev1$res1 = log(res1 + abs(min(res1)) + delta)} else (
            stop ("This transfromation is not supported in this version!")))))
    # The '+ abs(min(res1))' above is to set possible negative values to 0.

    # vgm of the residuals
    sp::coordinates(data.dev1) = ~ long + lat
    vgm1 <- gstat::variogram(object = formula, data.dev1, alpha = alpha)
    model.1 <- gstat::fit.variogram(vgm1, gstat::vgm(mean(vgm1$gamma), vgm.args, mean(vgm1$dist), min(vgm1$gamma)/10, anis = anis))
    if (model.1$range[2] <= 0) (cat("A zero or negative range was fitted to variogram", "\n"))
    if (model.1$range[2] <= 0) (model.1$range[2] <- min(vgm1$dist))  # set negative range to be positive

    # krige predictions
    sp::coordinates(data.pred1) = ~long + lat
    pred.krige1 <- gstat::krige(formula = formula, data.dev1, data.pred1, model = model.1, nmax=nmaxkrige, block = block, beta = beta)$var1.pred

    if (transformation == "none") {pred.krige = pred.krige1}
    if (transformation == "sqrt") {pred.krige = pred.krige1 ^ 2 - abs(min(res1))}
    if (transformation == "arcsine") {pred.krige = (sin(pred.krige1)) ^ 2 * 100 -  abs(min(res1))}
    if (transformation == "log") {pred.krige = exp(pred.krige1) - abs(min(res1)) - delta}

    cv.pred[idx == i] <- (pred.krige * (2 - lambda) + pred.idw1$res1.pred * lambda + pred.rf1* hybrid.parameter) / hybrid.parameter
    }
  }

  # predicitve error and accuracy assessment
  if (predacc == "VEcv") {predictive.accuracy = spm::vecv(trainy, cv.pred)} else (
  if (predacc == "ALL") {predictive.accuracy = spm::pred.acc(trainy, cv.pred)} else (
  stop ("This measure is not supported in this version!")))
  predictive.accuracy
}

