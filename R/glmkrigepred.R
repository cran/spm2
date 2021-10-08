#' @title Generate spatial predictions using the hybrid method of generalised
#' linear models  ('glm') and 'krige'
#'
#' @description This function is for generating spatial predictions using the hybrid method of
#' 'glm' and 'krige', including all methods implemented
#' in 'glmkrigecv'. (see reference #1 for further info).
#'
#' @param formula.glm a formula defining the response variable and predictive variables for 'glm'.
#' @param longlat	a dataframe contains longitude and latitude of point samples.
#' @param trainxy a dataframe contains longitude (long), latitude (lat),
#' predictive variables and the response variable of point samples.
#' @param predx	a dataframe or matrix contains columns of predictive variables
#' for the grids to be predicted.
#' @param y a vector of the response variable in the formula, that is, the left
#' part of the formula.
#' @param longlatpredx	a dataframe contains longitude and latitude of point locations
#' (i.e., the centers of grids) to be predicted.
#' @param family a description of the error distribution and link function to
#' be used in the model. See '?glm' for details.
#' @param transformation transform the residuals of 'glm' to normalise the data;
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
#' @param ... other arguments passed on to 'glm' and 'krige'.
#'
#' @return A dataframe of longitude, latitude, and predictions.
#'
#' @references Li, J., Alvarez, B., Siwabessy, J., Tran, M., Huang, Z.,
#' Przeslawski, R., Radke, L., Howard, F. and Nichol, S. (2017). "Application
#' of random forest, generalised linear model and their hybrid methods with
#' geostatistical techniques to count data: Predicting sponge species richness."
#' Environmental Modelling & Software 97: 112-129.
#'
#' Pebesma, E.J., 2004. Multivariable geostatistics in S: the gstat package.
#' Computers & Geosciences, 30: 683-691.
#'
#' @author Jin Li
#' @examples
#' \donttest{
#' library(spm)
#'
#' data(petrel)
#' data(petrel.grid)
#'
#' gravel <- petrel[, c(1, 2, 6:9, 5)]
#' longlat <- petrel[, c(1, 2)]
#' model <- log(gravel + 1) ~  lat +  bathy + I(long^3) + I(lat^2) + I(lat^3)
#' y <- log(gravel[, 7] +1)
#'
#' glmkrigepred1 <- glmkrigepred(formula.glm = model, longlat = longlat, trainxy =  gravel,
#' predx = petrel.grid, y = y, longlatpredx = petrel.grid[, c(1:2)],
#' transformation = "none", formula.krige = res1 ~ 1, vgm.args = "Sph", nmaxkrige = 12)
#'  # Since the default 'family' is used, actually a 'lm' model is used.
#'
#' names(glmkrigepred1)
#'
#' # Back transform 'glmkrigepred$predictions' to generate the final predictions
#' glmkrige.predictions <- exp(glmkrigepred1$predictions) - 1
#' range(glmkrige.predictions)
#'}
#'
#' @export
glmkrigepred <- function (formula.glm = NULL, longlat, trainxy, predx, y, longlatpredx, family = "gaussian", transformation = "none", delta = 1, formula.krige = res1 ~ 1, vgm.args = c("Sph"), anis = c(0, 1), alpha = 0, block = 0, beta, nmaxkrige = 12, ...) {

  names(longlat) <- c("long", "lat")
  names(longlatpredx) <- c("long", "lat")

  # glm modeling
  glm1 <- stats::glm(formula.glm, trainxy, family = family)
  # glm predictions
  pred.glm1 <- stats::predict(glm1, predx, type = "response")

  # the residuals of glm for krige
  data.dev1 <- longlat
  data.pred1 <- longlatpredx

  dev.glm1 <- stats::predict(glm1, trainxy, type="response")
  res1 <- y - dev.glm1

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

  predictions <- pred.krige + pred.glm1
  glmkrige.pred <- cbind(longlatpredx, predictions)
  glmkrige.pred
}

