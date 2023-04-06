#' @title Cross validation, n-fold and leave-one-out for generalized least squares ('gls')
#'
#' @description This function is a cross validation function for 'gls' method
#' in 'nlme' package.
#'
#' @param model a formula defining the response variable and predictive variables.
#' @param trainxy a dataframe contains longitude (long), latitude (lat),
#' predictive variables and the response variable of point samples. That is,
#' the location information must be names as 'long' and 'lat'.
#' @param y a vector of the response variable in the formula, that is, the left
#' part of the formula.
#' @param weights describing the within-group heteroscedasticity structure. Defaults
#'  to "NULL", corresponding to homoscedastic errors. See '?gls' in 'nlme'
#' for details.
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
#' @param ... other arguments passed on to 'gls'.
#'
#' @return A list with the following components:
#'  me, rme, mae, rmae, mse, rmse, rrmse, vecv and e1; or vecv only
#' @note This function is largely based on rfcv in 'randomForest' and
#' 'gls' in 'library(nlme)'.
#'
#' @references Pinheiro, J. C. and D. M. Bates (2000). Mixed-Effects Models
#' in S and S-PLUS. New York, Springer.
#'
#' @author Jin Li
#' @examples
#' \donttest{
#' library(spm)
#' library(nlme)
#'
#' data(petrel)
#' gravel <- petrel[, c(1, 2, 6:9, 5)]
#' range1 <- 0.8
#' nugget1 <- 0.5
#'
#' model <- log(gravel + 1) ~  long + lat +  bathy + dist + I(long^2) + I(lat^2) +
#' I(lat^3) + I(bathy^2) + I(bathy^3) + I(dist^2) + I(dist^3) + I(relief^2) + I(relief^3)
#'
#'
#' glscv1 <- glscv(model = model, gravel, log(gravel[, 7] +1), validation = "CV",
#'  corr.args = corSpher(c(range1, nugget1), form = ~ long + lat, nugget = TRUE),
#'  predacc = "ALL")
#' glscv1
#'
#' #For gls
#' set.seed(1234)
#' n <- 20 # number of iterations,60 to 100 is recommended.
#' VEcv <- NULL
#' for (i in 1:n) {
#' glscv1 <- glscv(model = model, gravel, log(gravel[, 7] +1), validation = "CV",
#'           corr.args = corSpher(c(range1, nugget1), form = ~ long + lat,
#'           nugget = TRUE), predacc = "VEcv")
#' VEcv [i] <- glscv1
#' }
#' plot(VEcv ~ c(1:n), xlab = "Iteration for GLS", ylab = "VEcv (%)")
#' points(cumsum(VEcv) / c(1:n) ~ c(1:n), col = 2)
#' abline(h = mean(VEcv), col = 'blue', lwd = 2)
#'
#' # For lm, that is, gls with 'correlation = NULL'
#' n <- 20 # number of iterations,60 to 100 is recommended.
#' VEcv <- NULL
#' set.seed(1234)
#' for (i in 1:n) {
#' glscv1 <- glscv(model = model, gravel, log(gravel[, 7] +1),
#' validation = "CV", predacc = "VEcv")
#' VEcv [i] <- glscv1
#' }
#' plot(VEcv ~ c(1:n), xlab = "Iteration for GLS", ylab = "VEcv (%)")
#' points(cumsum(VEcv) / c(1:n) ~ c(1:n), col = 2)
#' abline(h = mean(VEcv), col = 'blue', lwd = 2)
#'}
#'
#' @export
glscv <- function (model = var1 ~ 1, trainxy, y, corr.args = NULL, weights = NULL, validation = "CV", cv.fold = 10, predacc = "VEcv", ...) {

  if (validation == "LOO") {idx <- 1:length(y)}
  if (validation == "CV")  {idx <- datasplit(y, k.fold = cv.fold)}

  # cross validation
  cv.pred <- NULL

  if (validation == "LOO") {
    for (i in 1 : length(y)) {
      data.dev <- trainxy[idx != i, , drop = FALSE]
      data.pred <- trainxy[idx == i, , drop = FALSE]
      gls1 <- nlme::gls(model, data.dev, correlation = corr.args, weights = weights)
      cv.pred[idx == i] <- stats::predict(gls1, data.pred, type = "response")
    }
  }

  if (validation == "CV") {
  for (i in 1 : cv.fold) {
    data.dev <- trainxy[idx != i, , drop = FALSE]
    data.pred <- trainxy[idx == i, , drop = FALSE]
    gls1 <- nlme::gls(model, data.dev, correlation = corr.args, weights = weights)
    cv.pred[idx == i] <- stats::predict(gls1, data.pred, type = "response")
  }
  }

  # predicitve error and accuracy assessment
  if (predacc == "VEcv") {predictive.accuracy = spm::vecv(y, cv.pred)} else (
  if (predacc == "ALL") {predictive.accuracy = spm::pred.acc(y, cv.pred)} else (
  stop ("This measure is not supported in this version!")))
  predictive.accuracy
}
