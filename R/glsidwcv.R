#' @title Cross validation, n-fold and leave-one-out for the hybrid method of
#'  generalized least squares ('gls') and inverse distance weighted ('idw')
#'  (glsidw)
#'
#' @description This function is a cross validation function
#' for the hybrid method  of 'gls' and 'idw', where the data splitting is based
#'  on a stratified random  sampling method (see the 'datasplit' function for details)
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
#' @param idp	 a numeric number specifying the inverse distance weighting power.
#' @param nmaxidw for a local predicting: the number of nearest observations that
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
#' @param ... other arguments passed on to 'gls' and 'gstat'.
#'
#' @return A list with the following components:
#'  me, rme, mae, rmae, mse, rmse, rrmse, vecv and e1; or vecv only.
#'
#' @note This function is largely based on rfcv in 'randomForest' and
#' 'gls' in 'library(nlme)'.
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
#' glsidwcv1 <- glsidwcv(model = model, longlat = longlat, trainxy = gravel,
#' y = log(gravel[, 7] +1), idp = 2, nmaxidw = 12, validation = "CV",
#'  corr.args = corSpher(c(range1, nugget1), form = ~ lat + long, nugget = T),
#'  predacc = "ALL")
#' glsidwcv1
#'
#' # For glsidw
#' set.seed(1234)
#' n <- 20 # number of iterations,60 to 100 is recommended.
#' VEcv <- NULL
#' for (i in 1:n) {
#' glsidwcv1 <- glsidwcv(model = model, longlat = longlat, trainxy = gravel,
#' y = log(gravel[, 7] +1), idp = 2, nmaxidw = 12, validation = "CV",
#' corr.args = corSpher(c(range1, nugget1), form = ~ lat + long, nugget = T),
#' predacc = "VEcv")
#' VEcv [i] <- glsidwcv1
#' }
#' plot(VEcv ~ c(1:n), xlab = "Iteration for GLSIDW", ylab = "VEcv (%)")
#' points(cumsum(VEcv) / c(1:n) ~ c(1:n), col = 2)
#' abline(h = mean(VEcv), col = 'blue', lwd = 2)
#'}
#'
#' @export
glsidwcv <- function (model = var1 ~ 1, longlat, trainxy, y, corr.args = NULL, weights = NULL, idp = 2, nmaxidw = 12, validation = "CV", cv.fold = 10, predacc = "VEcv", ...) {

  if (validation == "LOO") {idx <- 1:length(y)}
  if (validation == "CV")  {idx <- datasplit(y, k.fold = cv.fold)}

  names(longlat) <- c("long", "lat")

  # cross validation
  cv.pred <- NULL

  if (validation == "LOO") {
    for (i in 1 : length(y)) {
      data.dev <- trainxy[idx != i, , drop = FALSE]
      data.pred <- trainxy[idx == i, , drop = FALSE]

      # gls modeling
      gls1 <- nlme::gls(model, data.dev, correlation = corr.args, weights = weights)

      # gls predictions
      pred.gls1 <- stats::predict(gls1, data.pred, type = "response")

      # the residuals of gls for idw
      dev.gls1 <- stats::predict(gls1, data.dev, type="response")

      data.dev1 <- longlat[idx != i, , drop = FALSE] # for idw
      data.pred1 <- longlat[idx == i, , drop = FALSE] # for idw

      data.dev1$res1 <- y[idx != i] - dev.gls1

      # idw of the residuals
      gstat1 <- gstat::gstat(id = "res1", formula = res1 ~ 1, locations = ~ long + lat, data = data.dev1, set = list(idp = idp), nmax = nmaxidw)

      # idw predictions
      pred.idw1<- stats::predict(gstat1, data.pred1)

      cv.pred[idx == i] <- pred.idw1$res1.pred + pred.gls1
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

    # the residuals of gls for idw
    dev.gls1 <- stats::predict(gls1, data.dev, type="response")

    data.dev1 <- longlat[idx != i, , drop = FALSE] # for idw
    data.pred1 <- longlat[idx == i, , drop = FALSE] # for idw

    data.dev1$res1 <- y[idx != i] - dev.gls1

    # idw of the residuals
    gstat1 <- gstat::gstat(id = "res1", formula = res1 ~ 1, locations = ~ long + lat, data = data.dev1, set = list(idp = idp), nmax = nmaxidw)

    # idw predictions
    pred.idw1<- stats::predict(gstat1, data.pred1)

    cv.pred[idx == i] <- pred.idw1$res1.pred + pred.gls1
  }
  }

  # predictive error and accuracy assessment
  if (predacc == "VEcv") {predictive.accuracy = spm::vecv(y, cv.pred)} else (
  if (predacc == "ALL") {predictive.accuracy = spm::pred.acc(y, cv.pred)} else (
  stop ("This measure is not supported in this version!")))
  predictive.accuracy
}
