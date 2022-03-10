#' @title Cross validation, n-fold and leave-one-out for thin plate splines ('TPS')
#'
#' @description This function is a cross validation function for 'Tps' method
#' in 'fields' package.
#'
#' @param trainx a dataframe contains longitude (long), latitude (lat) and
#' predictive variables of point samples. That is, they must be names as
#'  'long' and 'lat'.
#' @param trainy a vector of response, must have length equal to the number
#'  of rows in trainx.
#' @param m A polynomial function of degree (m-1) will be included in the
#' model as the drift (or spatial trend) component. Default is 'm = NULL'
#' that is the value such that 2m-d is greater than zero where d is the
#' dimension of x.
#' @param p polynomial power for Wendland radial basis functions as in 'Tps'.
#' 'p = NULL' that leads to a default value of 2 for spatial predictive
#'  modelling based on 'x' containing only the location information.
#' @param theta the tapering range. 'theta = 3' degrees is a very generous
#' taper range. For spatial predictive modeling the taper should be large
#' enough to about 20 non zero nearest neighbors for every location.
#' @param lon.lat if 'TRUE' locations are interpreted as longitude and
#'  latitude and great circle distance is used to find distances among
#'  locations.
#' @param lambda smoothing parameter, the default is 'NULL'. See '?Tps' for further info.
#' @param validation validation methods, include 'LOO': leave-one-out, and
#'  'CV': cross-validation.
#' @param cv.fold integer; number of folds in the cross-validation. if > 1,
#' then apply n-fold cross validation; the default is 10, i.e., 10-fold cross
#' validation that is recommended.
#' @param predacc can be either "VEcv" for vecv or "ALL" for all measures
#' in function pred.acc.
#' @param ... other arguments passed on to 'gstat'.
#'
#' @return A list with the following components:
#'  me, rme, mae, rmae, mse, rmse, rrmse, vecv and e1; or vecv only
#' @note This function is largely based on 'krigecv' in this package and
#'  'Tps' and 'fastTpsMLE' in 'fields' package.
#'
#' @references Douglas Nychka and Reinhard Furrer and John Paige and Stephan
#' Sain and Florian Gerber and Matthew Iverson, 2020. fields: Tools for
#' Spatial Data, R package version 10.3
#' {https://CRAN.R-project.org/package=fields}.
#'
#' @author Jin Li
#' @examples
#' \donttest{
#' library(fields)
#' library(spm)
#' data(petrel)
#'
#' tpscv1 <- tpscv(petrel[, c(1,2)], petrel[, 5], cv.fold = 5, predacc = "VEcv")
#' tpscv1
#'
#' tpscv1 <- tpscv(petrel[, c(1,2)], petrel[, 5], lambda = 0.9, cv.fold = 5, predacc = "VEcv")
#' tpscv1
#'
#' tpscv1 <- tpscv(petrel[, c(1,2)], petrel[, 5], validation = "LOO", predacc = "VEcv")
#' tpscv1
#'
#' set.seed(1234)
#' n <- 20 # number of iterations,60 to 100 is recommended.
#' VEcv <- NULL
#' for (i in 1:n) {
#' tpscv1 <- tpscv(petrel[, c(1,2)], petrel[, 5], cv.fold = 10,  lambda = 0.13, predacc = "VEcv")
#' VEcv [i] <- tpscv1
#' }
#' plot(VEcv ~ c(1:n), xlab = "Iteration for TPS", ylab = "VEcv (%)")
#' points(cumsum(VEcv) / c(1:n) ~ c(1:n), col = 2)
#' abline(h = mean(VEcv), col = 'blue', lwd = 2)
#'
#' n <- 20 # number of iterations,60 to 100 is recommended.
#' # set.seed(1234)
#' VEcv <- NULL
#' for (i in 1:n) {
#' set.seed(1234 + i) # set random seed for each iteration. You can remove
#' # this line and use above set.seed(1234) and see what you can get.
#' tpscv1 <- tpscv(petrel[, c(1,2)], petrel[, 5], predacc = "VEcv")
#' VEcv [i] <- tpscv1
#' }
#' plot(VEcv ~ c(1:n), xlab = "Iteration for TPS", ylab = "VEcv (%)")
#' points(cumsum(VEcv) / c(1:n) ~ c(1:n), col = 2)
#' abline(h = mean(VEcv), col = 'blue', lwd = 2)
#'}
#' @export
tpscv <- function (trainx, trainy, m = NULL, p = NULL, theta = 3, lambda = NULL, lon.lat = TRUE, validation = "CV", cv.fold = 10, predacc = "VEcv", ...) {
  classtrainy <- is.factor(trainy)
  if (classtrainy) {
    stop ("This function is not for categorical response variable")
  } else (
  if (validation == "LOO") {idx <- 1:length(trainy)} else (
  if (validation == "CV")  {idx <- datasplit(trainy, k.fold = cv.fold)} else (stop ("This validation method is not supported in this version!"))))

  # cross validation
  cv.pred <- NULL
  if (validation == "LOO") {
    for (i in 1 : length(trainy)) {
      x <- trainx[idx != i, , drop = FALSE]
      y <- trainy[idx != i]
      x.pred <- trainx[idx == i, , drop = FALSE]

      tps1 <- fields::Tps(x, y, m = m, p = p, lon.lat = lon.lat, lambda = lambda)
      cv.pred[idx == i] <- as.vector(stats::predict(tps1, x.pred))
    }
  }

  if (validation == "CV") {
    for (i in 1 : cv.fold) {
      x <- trainx[idx != i, , drop = FALSE]
      y <- trainy[idx != i]
      x.pred <- trainx[idx == i, , drop = FALSE]

      tps1 <- fields::Tps(x, y, m = m, p = p, lon.lat = lon.lat, lambda = lambda)
      cv.pred[idx == i] <- as.vector(stats::predict(tps1, x.pred))
    }
  }

  if (predacc == "VEcv") {predictive.accuracy = spm::vecv(trainy, cv.pred)} else (
  if (predacc == "ALL") {predictive.accuracy = spm::pred.acc(trainy, cv.pred)} else (
  stop ("This measure is not supported in this version!")))
  predictive.accuracy
}
