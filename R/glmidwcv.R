#' @title Cross validation, n-fold and leave-one-out for the hybrid method of
#' generalised linear models  ('glm') and inverse distance weighted ('IDW') ('glmidw')
#'
#' @description This function is a cross validation function
#' for the hybrid method of 'glm' and 'idw' using 'gstat' (glmidw) (see
#'  reference #1), where the  data splitting is based on a stratified random
#'  sampling method (see the  'datasplit' function for details).
#'
#' @param formula a formula defining the response variable and predictive variables
#'  for 'glm'.
#' @param longlat	a dataframe contains longitude and latitude of point samples.
#' @param trainxy a dataframe contains longitude (long), latitude (lat),
#' predictive variables and the response variable of point samples.
#' @param y a vector of the response variable in the formula, that is, the left
#' part of the formula.
#' @param family a description of the error distribution and link function to
#' be used in the model. See '?glm' for details.
#' @param idp	 a numeric number specifying the inverse distance weighting power.
#' @param nmaxidw for a local predicting: the number of nearest observations that
#'  should be used for a prediction or simulation, where nearest is defined in
#'  terms of the space of the spatial locations. By default, 12 observations
#'  are used.
#' @param validation validation methods, include 'LOO': leave-one-out, and 'CV':
#' cross-validation.
#' @param cv.fold integer; number of folds in the cross-validation. if > 1,
#' then apply n-fold cross validation; the default is 10, i.e., 10-fold cross
#' validation that is recommended.
#' @param predacc can be either "VEcv" for vecv or "ALL" for all measures
#' in function pred.acc.
#' @param ... other arguments passed on to 'glm' and 'gstat'.
#'
#' @return A list with the following components:
#'  me, rme, mae, rmae, mse, rmse, rrmse, vecv and e1; or vecv only.
#'
#' @note This function is largely based on 'rfcv' in 'randomForest', 'idwcv'
#' in 'spm'and 'glm' in 'stats'.
#'
#' @references Li, J., Alvarez, B., Siwabessy, J., Tran, M., Huang, Z.,
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
#'
#' data(petrel)
#' gravel <- petrel[, c(1, 2, 6:9, 5)]
#' longlat <- petrel[, c(1, 2)]
#' model <- log(gravel + 1) ~  lat +  bathy + I(long^3) + I(lat^2) + I(lat^3)
#' y <- log(gravel[, 7] +1)
#'
#' set.seed(1234)
#' glmidwcv1 <- glmidwcv(formula = model, longlat = longlat, trainxy =  gravel,
#' y = y, idp = 2, nmaxidw = 12, validation = "CV", predacc = "ALL")
#' glmidwcv1 # Since the default 'family' is used, actually a 'lm' model is used.
#'
#' data(spongelonglat)
#' longlat <- spongelonglat[, 7:8]
#' model <- sponge ~ long + I(long^2)
#' y = spongelonglat[, 1]
#' set.seed(1234)
#' glmidwcv1 <- glmidwcv(formula = model, longlat = longlat, trainxy = spongelonglat,
#' y = y, family = poisson, idp = 2, nmaxidw = 12, validation = "CV",
#' predacc = "ALL")
#' glmidwcv1
#'
#' # glmidw for count data
#' data(spongelonglat)
#' longlat <- spongelonglat[, 7:8]
#' model <- sponge ~ . # use all predictive variables in the dataset
#' y = spongelonglat[, 1]
#' set.seed(1234)
#' n <- 20 # number of iterations,60 to 100 is recommended.
#' VEcv <- NULL
#' for (i in 1:n) {
#'  glmidwcv1 <- glmidwcv(formula = model, longlat = longlat, trainxy = spongelonglat,
#'  y = y, family = poisson, idp = 2, nmaxidw = 12, validation = "CV",
#'  predacc = "VEcv")
#'  VEcv [i] <- glmidwcv1
#'  }
#'  plot(VEcv ~ c(1:n), xlab = "Iteration for GLM", ylab = "VEcv (%)")
#'  points(cumsum(VEcv) / c(1:n) ~ c(1:n), col = 2)
#'  abline(h = mean(VEcv), col = 'blue', lwd = 2)
#'
#' # glmidw for percentage data
#' library(MASS)
#' longlat <- petrel[, c(1, 2)]
#' model <- gravel / 100 ~  lat +  bathy + I(long^3) + I(lat^2) + I(lat^3)
#' set.seed(1234)
#' n <- 20 # number of iterations,60 to 100 is recommended.
#' VEcv <- NULL
#' for (i in 1:n) {
#' glmidwcv1 <- glmcv(formula = model, longlat = longlat, trainxy = gravel,
#' y = gravel[, 7] / 100, family = binomial(link=logit), idp = 2, nmaxidw = 12,
#' validation = "CV", predacc = "VEcv")
#' VEcv [i] <- glmidwcv1
#' }
#' plot(VEcv ~ c(1:n), xlab = "Iteration for GLM", ylab = "VEcv (%)")
#' points(cumsum(VEcv) / c(1:n) ~ c(1:n), col = 2)
#' abline(h = mean(VEcv), col = 'blue', lwd = 2)
#'}
#'
#' @export
glmidwcv <- function (formula = NULL, longlat, trainxy, y, family = "gaussian", idp = 2, nmaxidw = 12, validation = "CV", cv.fold = 10, predacc = "VEcv", ...) {

  if (validation == "LOO") {idx <- 1:length(y)}
  if (validation == "CV")  {idx <- datasplit(y, k.fold = cv.fold)}

  names(longlat) <- c("long", "lat")

  # cross validation
  cv.pred <- NULL

  if (validation == "LOO") {
    for (i in 1 : length(y)) {
      data.dev <- trainxy[idx != i, , drop = FALSE]
      data.pred <- trainxy[idx == i, , drop = FALSE]

      # glm modeling
      glm1 <- stats::glm(formula, data.dev, family = family)

      # glm predictions
      pred.glm1 <- stats::predict(glm1, data.pred, type = "response")

      # the residuals of glm for idw
      data.dev1 <- longlat[idx != i, , drop = FALSE] # for idw
      data.pred1 <- longlat[idx == i, , drop = FALSE] # for idw

      #data.dev1$res1 <- stats::predict(glm1, data.dev, type="response")$residuals
      dev.glm1 <- stats::predict(glm1, data.dev, type="response")
      data.dev1$res1 <- y[idx != i] - dev.glm1

      # idw of the residuals
      gstat1 <- gstat::gstat(id = "res1", formula = res1 ~ 1, locations = ~ long + lat, data = data.dev1, set = list(idp = idp), nmax = nmaxidw)

      # idw predictions
      pred.idw1<- stats::predict(gstat1, data.pred1)

      cv.pred[idx == i] <- pred.idw1$res1.pred + pred.glm1
    }
  }

  if (validation == "CV") {
  for (i in 1 : cv.fold) {
    data.dev <- trainxy[idx != i, , drop = FALSE]
    data.pred <- trainxy[idx == i, , drop = FALSE]

    # glm modeling
    glm1 <- stats::glm(formula, data.dev, family = family)

    # glm predictions
    pred.glm1 <- stats::predict(glm1, data.pred, type = "response")

    # the residuals of glm for idw

    data.dev1 <- longlat[idx != i, , drop = FALSE] # for idw
    data.pred1 <- longlat[idx == i, , drop = FALSE] # for idw

    #data.dev1$res1 <- stats::predict(glm1, data.dev, type="response")$residuals
    dev.glm1 <- stats::predict(glm1, data.dev, type="response")
    data.dev1$res1 <- y[idx != i] - dev.glm1

    # idw of the residuals
    gstat1 <- gstat::gstat(id = "res1", formula = res1 ~ 1, locations = ~ long + lat, data = data.dev1, set = list(idp = idp), nmax = nmaxidw)

    # idw predictions
    pred.idw1<- stats::predict(gstat1, data.pred1)

    cv.pred[idx == i] <- pred.idw1$res1.pred + pred.glm1
    }
  }

  # predictive error and accuracy assessment
  if (predacc == "VEcv") {predictive.accuracy = spm::vecv(y, cv.pred)} else (
  if (predacc == "ALL") {predictive.accuracy = spm::pred.acc(y, cv.pred)} else (
  stop ("This measure is not supported in this version!")))
  predictive.accuracy
}
