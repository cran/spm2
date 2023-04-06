#' @title Cross validation, n-fold and leave-one-out for the hybrid methods of
#' generalised linear models  ('glm'), 'kriging' and inverse distance weighted ('IDW').
#'
#' @description This function is a cross validation function
#' for 38 hybrid  methods of 'glm', 'kriging' and 'IDW', including the average
#' of 'glmkrige' and 'glmidw' ('glmkrigeglmidw') and  the average of 'glm',
#' 'glmkrige' and 'glmidw' ('glmglmkrigeglmidw'), where 'kriging' methods
#' include ordinary kriging  ('OK'), simple kriging ('SK'), block 'OK' ('BOK')
#' and block 'SK'('BSK') and 'IDW' also covers 'NN' and 'KNN' (for details, see
#' reference #1). This function can also be sued for 38 hybrid methods of 'lm',
#' 'kriging' and 'IDW'.
#'
#' @param formula.glm a formula defining the response variable and predictive variables for 'glm'.
#' @param longlat	a dataframe contains longitude and latitude of point samples.
#' @param trainxy a dataframe contains longitude (long), latitude (lat),
#' predictive variables and the response variable of point samples.
#' @param y a vector of the response variable in the formula, that is, the left
#' part of the formula.
#' @param family a description of the error distribution and link function to
#' be used in the model. See '?glm' for details.
#' @param transformation transform the residuals of 'glm' to normalise the data for 'krige';
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
#' @param hybrid.parameter the default is 2 that is for 'glmkrigeglmidw';
#' for 'glmglmkrigeglmidw', it needs to be 3.
#' @param lambda, ranging from 0 to 2; the default is 1 for 'glmkrigeglmidw'
#' and 'glmglmkrigeglmidw'; and if it is < 1, more weight is placed on 'krige',
#' otherwise more weight is placed on 'idw'; and if it is 0, 'idw' is not
#' considered and the resultant methods is 'glmkrige' when the default
#' 'hybrid.parameter' is used; and if it is 2, then the resultant method is
#' 'glmidw' when the default 'hybrid.parameter' is used.
#' @param validation validation methods, include 'LOO': leave-one-out, and 'CV':
#' cross-validation.
#' @param cv.fold integer; number of folds in the cross-validation. if > 1,
#' then apply n-fold cross validation; the default is 10, i.e., 10-fold cross
#' validation that is recommended.
#' @param predacc can be either "VEcv" for vecv or "ALL" for all measures
#' in function pred.acc.
#' @param ... other arguments passed on to 'glm', 'krige' and 'gstat'.
#'
#' @return A list with the following components:
#'  me, rme, mae, rmae, mse, rmse, rrmse, vecv and e1; or vecv only
#' @note This function is largely based on 'rfcv' in 'randomForest', 'krigecv'
#' in 'spm2'and 'glm' in 'stats'.
#'
#' @references Li, J. (2022). Spatial Predictive Modeling with R. Boca Raton,
#' Chapman and Hall/CRC.
#'
#' Li, J., Alvarez, B., Siwabessy, J., Tran, M., Huang, Z.,
#' Przeslawski, R., Radke, L., Howard, F. and Nichol, S. (2017). "Application
#' of random forest, generalised linear model and their hybrid methods with
#' geostatistical techniques to count data: Predicting sponge species richness."
#' Environmental Modelling & Software 97: 112-129.
#'
#' A. Liaw and M. Wiener (2002). Classification and Regression by
#' randomForest. R News 2(3), 18-22.
#'
#' Pebesma, E.J., 2004. Multivariable geostatistics in S: the gstat package.
#' Computers & Geosciences, 30: 683-691.
#'
#' @author Jin Li
#' @examples
#' \donttest{
#' library(spm)
#' # glmokglidw
#' data(petrel)
#' gravel <- petrel[, c(1, 2, 6:9, 5)]
#' longlat <- petrel[, c(1, 2)]
#' model <- log(gravel + 1) ~  lat +  bathy + I(long^3) + I(lat^2) + I(lat^3)
#' y <- log(gravel[, 7] +1)
#' set.seed(1234)
#' glmkrigeglmidwcv1 <- glmkrigeidwcv(formula.glm = model, longlat = longlat,
#' trainxy =  gravel, y = y, transformation = "none", formula.krige = res1 ~ 1,
#' vgm.args = "Sph", nmaxkrige = 12, idp = 2, nmaxidw = 12, validation = "CV",
#'  predacc = "ALL")
#' glmkrigeglmidwcv1 # Since the default 'family' is used, actually a 'lm' model is used.
#'
#' # glmokglmidw
#' data(spongelonglat)
#' longlat <- spongelonglat[, 7:8]
#' model <- sponge ~ long + I(long^2)
#' y = spongelonglat[, 1]
#' set.seed(1234)
#' glmkrigeglmidwcv1 <- glmkrigeidwcv(formula.glm = model, longlat = longlat,
#' trainxy = spongelonglat, y = y, family = poisson, transformation = "arcsine",
#' formula.krige = res1 ~ 1, vgm.args = ("Sph"), nmaxkrige = 12, idp = 2,
#' nmaxidw = 12, validation = "CV", predacc = "ALL")
#' glmkrigeglmidwcv1
#'
#' # glmglmokglmidw
#' data(spongelonglat)
#' longlat <- spongelonglat[, 7:8]
#' model <- sponge ~ long + I(long^2)
#' y = spongelonglat[, 1]
#' set.seed(1234)
#' glmglmkrigeglmidwcv1 <- glmkrigeidwcv(formula.glm = model, longlat = longlat,
#' trainxy = spongelonglat, y = y, family = poisson, transformation = "arcsine",
#' formula.krige = res1 ~ 1, vgm.args = ("Sph"), nmaxkrige = 12, idp = 2,
#' nmaxidw = 12, hybrid.parameter = 3, validation = "CV", predacc = "ALL")
#' glmglmkrigeglmidwcv1
#'
#' # glmokglidw for count data
#' data(spongelonglat)
#' longlat <- spongelonglat[, 7:8]
#' model <- sponge ~ . # use all predictive variables in the dataset
#' y = spongelonglat[, 1]
#' set.seed(1234)
#' n <- 20 # number of iterations,60 to 100 is recommended.
#' VEcv <- NULL
#' for (i in 1:n) {
#'  glmkrigeglmidwcv1 <- glmkrigeidwcv(formula.glm = model, longlat = longlat,
#'  trainxy = spongelonglat, y = y, family = poisson, formula.krige = res1 ~ 1,
#'  vgm.args = ("Sph"), nmaxkrige = 12, idp = 2, nmaxidw = 12, validation = "CV",
#'  predacc = "VEcv")
#'  VEcv [i] <- glmkrigeglmidwcv1
#'  }
#'  plot(VEcv ~ c(1:n), xlab = "Iteration for GLM", ylab = "VEcv (%)")
#'  points(cumsum(VEcv) / c(1:n) ~ c(1:n), col = 2)
#'  abline(h = mean(VEcv), col = 'blue', lwd = 2)
#'
#' # glmokglmidw for percentage data
#' longlat <- petrel[, c(1, 2)]
#' model <- gravel / 100 ~  lat +  bathy + I(long^3) + I(lat^2) + I(lat^3)
#' set.seed(1234)
#' n <- 20 # number of iterations,60 to 100 is recommended.
#' VEcv <- NULL
#' for (i in 1:n) {
#' glmkrigeglmidwcv1 <- glmkrigeidwcv(formula.glm = model, longlat = longlat,
#' trainxy = gravel, y = gravel[, 7] / 100, family = binomial(link=logit),
#' formula.krige = res1 ~ 1, vgm.args = ("Sph"), nmaxkrige = 12, idp = 2,
#' nmaxidw = 12, validation = "CV", predacc = "VEcv")
#' VEcv [i] <- glmkrigeglmidwcv1
#' }
#' plot(VEcv ~ c(1:n), xlab = "Iteration for GLM", ylab = "VEcv (%)")
#' points(cumsum(VEcv) / c(1:n) ~ c(1:n), col = 2)
#' abline(h = mean(VEcv), col = 'blue', lwd = 2)
#'}
#'
#' @export
glmkrigeidwcv <- function (formula.glm = NULL, longlat, trainxy, y, family = "gaussian", transformation = "none", delta = 1, formula.krige = res1 ~ 1, vgm.args = c("Sph"), anis = c(0, 1), alpha = 0, block = 0, beta, nmaxkrige = 12, idp = 2, nmaxidw = 12, hybrid.parameter = 2, lambda = 1, validation = "CV", cv.fold = 10, predacc = "VEcv", ...) {

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

      # glm modelling
      glm1 <- stats::glm(formula.glm, data.dev, family = family)
      # glm predictions
      pred.glm1 <- stats::predict(glm1, data.pred, type = "response")

      # the residuals of glm for krige
      data.dev1 <- longlat[idx != i, , drop = FALSE] # for krige
      data.pred1 <- longlat[idx == i, , drop = FALSE] # for krige

      dev.glm1 <- stats::predict(glm1, data.dev, type="response")
      res1 <- y[idx != i] - dev.glm1
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

      cv.pred[idx == i] <- (pred.krige * (2 - lambda) + pred.idw1$res1.pred * lambda + pred.glm1 * hybrid.parameter) / hybrid.parameter
    }
  }

  if (validation == "CV") {
  for (i in 1 : cv.fold) {
    data.dev <- trainxy[idx != i, , drop = FALSE]
    data.pred <- trainxy[idx == i, , drop = FALSE]

    # glm modelling
    glm1 <- stats::glm(formula.glm, data.dev, family = family)

    # glm predictions
    pred.glm1 <- stats::predict(glm1, data.pred, type = "response")

    # the residuals of glm for krige
    data.dev1 <- longlat[idx != i, , drop = FALSE] # for krige
    data.pred1 <- longlat[idx == i, , drop = FALSE] # for krige

    dev.glm1 <- stats::predict(glm1, data.dev, type="response")
    res1 <- y[idx != i] - dev.glm1
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

    cv.pred[idx == i] <- (pred.krige * (2 - lambda) + pred.idw1$res1.pred * lambda + pred.glm1* hybrid.parameter) / hybrid.parameter
    }
  }

  # predicitve error and accuracy assessment
  if (predacc == "VEcv") {predictive.accuracy = spm::vecv(y, cv.pred)} else (
  if (predacc == "ALL") {predictive.accuracy = spm::pred.acc(y, cv.pred)} else (
  stop ("This measure is not supported in this version!")))
  predictive.accuracy
}

