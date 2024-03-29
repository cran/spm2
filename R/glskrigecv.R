#' @title Cross validation, n-fold and leave-one-out for the hybrid method of
#'  generalized least squares ('gls') and kriging ('krige') ('glskrige')
#'
#' @description This function is a cross validation function
#' for the hybrid method  of 'gls' and 'krige' ('glskrige'), where the data
#'  splitting is based on a stratified random  sampling method (see the
#'   'datasplit' function for details)
#'
#' @param model a formula defining the response variable and predictive variables.
#' @param longlat	a dataframe contains longitude and latitude of point samples.
#' @param trainxy a dataframe contains longitude (long), latitude (lat),
#' predictive variables and the response variable of point samples. That is,
#' the location information must be names as 'long' and 'lat'.
#' @param y a vector of the response variable in the formula, that is, the left
#' part of the formula.
#' @param weights describing the within-group heteroscedasticity structure. Defaults
#'  to "NULL", corresponding to homoscedastic errors. See '?gls' in 'nlme'
#' for details.
#' @param transformation transform the residuals of 'gls' to normalize the data;
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
#' @param corr.args arguments for 'correlation' in 'gls'. See '?corClasses' in 'nlme'
#' for details. By default, "NULL" is used. When "NULL" is used,
#' then 'gls' is actually performing 'lm'.
#' @param validation validation methods, include 'LOO': leave-one-out, and 'CV':
#' cross-validation.
#' @param cv.fold integer; number of folds in the cross-validation. if > 1,
#' then apply n-fold cross validation; the default is 10, i.e., 10-fold cross
#' validation that is recommended.
#' @param predacc can be either "VEcv" for vecv or "ALL" for all measures
#' in function pred.acc.
#' @param ... other arguments passed on to 'gls' and 'krige'.
#'
#' @return A list with the following components:
#'  me, rme, mae, rmae, mse, rmse, rrmse, vecv and e1; or vecv only.
#'
#' @note This function is largely based on rfcv in 'randomForest', 'krigecv' in 'spm2'
#'  and 'gls' in 'library(nlme)'.
#'
#' @references Pinheiro, J. C. and D. M. Bates (2000). Mixed-Effects Models
#' in S and S-PLUS. New York, Springer.
#'
#' Pebesma, E.J., 2004. Multivariable geostatistics in S: the gstat package.
#' Computers & Geosciences, 30: 683-691.
#'
#' @author Jin Li
#' @examples
#' \donttest{
#' library(spm)
#' library(nlme)
#'
#' data(petrel)
#' gravel <- petrel[, c(1, 2, 6:9, 5)]
#' longlat <- petrel[, c(1, 2)]
#' range1 <- 0.8
#' nugget1 <- 0.5
#' model <- log(gravel + 1) ~  long + lat +  bathy + dist + I(long^2) + I(lat^2) +
#' I(lat^3) + I(bathy^2) + I(bathy^3) + I(dist^2) + I(dist^3) + I(relief^2) + I(relief^3)
#'
#' glskrigecv1 <- glskrigecv(model = model, longlat = longlat, trainxy = gravel,
#' y = log(gravel[, 7] +1), transformation = "none", formula.krige = res1 ~ 1,
#' vgm.args = "Sph", nmaxkrige = 12, validation = "CV",
#'  corr.args = corSpher(c(range1, nugget1), form = ~ lat + long, nugget = TRUE),
#'  predacc = "ALL")
#' glskrigecv1
#'
#' # For glskrige
#' set.seed(1234)
#' n <- 20 # number of iterations,60 to 100 is recommended.
#' VEcv <- NULL
#' for (i in 1:n) {
#' glskrigecv1 <- glskrigecv(model = model, longlat = longlat, trainxy = gravel,
#' y = log(gravel[, 7] +1), transformation = "none", formula.krige = res1 ~ 1,
#' vgm.args = "Sph", nmaxok = 12, validation = "CV",
#' corr.args = corSpher(c(range1, nugget1), form = ~ lat + long, nugget = TRUE),
#'  predacc = "VEcv")
#' VEcv [i] <- glskrigecv1
#' }
#' plot(VEcv ~ c(1:n), xlab = "Iteration for GLSOK", ylab = "VEcv (%)")
#' points(cumsum(VEcv) / c(1:n) ~ c(1:n), col = 2)
#' abline(h = mean(VEcv), col = 'blue', lwd = 2)
#'}
#'
#' @export
glskrigecv <- function (model = var1 ~ 1, longlat, trainxy, y, corr.args = NULL, weights = NULL, transformation = "none", delta = 1, formula.krige = res1 ~ 1, vgm.args = c("Sph"), anis = c(0, 1), alpha = 0, block = 0, beta, nmaxkrige = 12, validation = "CV", cv.fold = 10, predacc = "VEcv", ...) {

  if (validation == "LOO") {idx <- 1:length(y)}
  if (validation == "CV")  {idx <- datasplit(y, k.fold = cv.fold)}

  names(longlat) <- c("long", "lat")

  # cross validation
  n <- nrow(trainxy)
  p <- ncol(trainxy) - 1
  cv.pred <- NULL

  if (validation == "LOO") {
    for (i in 1 : length(y)) {
      data.dev <- trainxy[idx != i, , drop = FALSE]
      data.pred <- trainxy[idx == i, , drop = FALSE]

      # gls modeling
      gls1 <- nlme::gls(model, data.dev, correlation = corr.args, weights = weights)
      # gls predictions
      pred.gls1 <- stats::predict(gls1, data.pred, type = "response")

      # the residuals of gls for krige
      dev.gls1 <- stats::predict(gls1, data.dev, type="response")
      data.dev1 <- longlat[idx != i, , drop = FALSE] # for krige
      data.pred1 <- longlat[idx == i, , drop = FALSE] # for krige
      res1 <- y[idx != i] - dev.gls1

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

      cv.pred[idx == i] <- pred.krige + pred.gls1
    }
  }

  if (validation == "CV") {
  for (i in 1 : cv.fold) {
    data.dev <- trainxy[idx != i, , drop = FALSE]
    data.pred <- trainxy[idx == i, , drop = FALSE]

    # gls modeling
    gls1 <- nlme::gls(model, data.dev, correlation = corr.args, weights = weights)
    # gls predictions
    pred.gls1 <- stats::predict(gls1, data.pred, type = "response")

    # the residuals of gls for krige
    dev.gls1 <- stats::predict(gls1, data.dev, type="response")
    data.dev1 <- longlat[idx != i, , drop = FALSE] # for krige
    data.pred1 <- longlat[idx == i, , drop = FALSE] # for krige

    res1 <- y[idx != i] - dev.gls1

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

    cv.pred[idx == i] <- pred.krige + pred.gls1
  }
  }

  # predictive error and accuracy assessment
  if (predacc == "VEcv") {predictive.accuracy = spm::vecv(y, cv.pred)} else (
  if (predacc == "ALL") {predictive.accuracy = spm::pred.acc(y, cv.pred)} else (
  stop ("This measure is not supported in this version!")))
  predictive.accuracy
}
