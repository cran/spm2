#' @title Cross validation, n-fold and leave-one-out for support vector machine ('svm')
#'
#' @description This function is a cross validation function
#' for 'svm' regression in 'e1071' package.
#'
#' @param formula a formula defining the response variable and predictive variables.
#' @param trainxy a dataframe contains predictive variables and the response
#' variable of point samples. The location information, longitude (long),
#' latitude (lat), need to be included in the 'trainx' for spatial predictive
#'  modelling, need to be named as 'long' and 'lat'.
#' @param y a vector of the response variable in the formula, that is, the left
#' part of the formula.
#' @param scale A logical vector indicating the variables to be scaled (default: TRUE).
#' @param type the default setting is 'NULL'. See '?svm' for various options.
#' @param kernel the default setting is 'radial'. See '?svm' for other options.
#' @param degree a parameter needed for kernel of type polynomial (default: 3).
#' @param gamma a parameter needed for all 'kernels' except 'linear'
#' (default: 1/(data dimension)).
#' @param coef0 a parameter needed for kernels of type 'polynomial' and 'sigmoid'(default: 0).
#' @param cost cost of constraints violation (default: 1).
#' @param nu a parameter needed for 'nu-classification', 'nu-regression', and 'one-classification' (default: 0.5).
#' @param tolerance tolerance of termination criterion (default: 0.001).
#' @param epsilon 'epsilon' in the insensitive-loss function (default: 0.1).
#' @param validation validation methods, include 'LOO': leave-one-out, and 'CV':
#' cross-validation.
#' @param cv.fold integer; number of folds in the cross-validation. if > 1,
#' then apply n-fold cross validation; the default is 10, i.e., 10-fold cross
#' validation that is recommended.
#' @param predacc can be either "VEcv" for 'vecv' or "ALL" for all measures
#' in function pred.acc.
#' @param ... other arguments passed on to 'svm'.
#'
#' @return A list with the following components:
#'  me, rme, mae, rmae, mse, rmse, rrmse, vecv and e1; or vecv only
#' @note This function is largely based on 'rfcv' in 'randomForest' and
#' 'svm' in 'e1071'.
#'
#' @references David Meyer, Evgenia Dimitriadou, Kurt Hornik, Andreas Weingessel and Friedrich
#' Leisch (2020). e1071: Misc Functions of the Department of Statistics, Probability
#' Theory Group (Formerly: E1071), TU Wien. R package version 1.7-4.
#' https://CRAN.R-project.org/package=e1071.
#'
#' @author Jin Li
#' @examples
#' \donttest{
#' library(spm)
#'
#' data(petrel)
#' gravel <- petrel[, c(1, 2, 6:9, 5)]
#' model <- log(gravel + 1) ~ lat +  bathy + I(long^3) + I(lat^2) + I(lat^3)
#' set.seed(1234)
#' svmcv1 <- svmcv(formula = model, gravel, log(gravel[, 7] +1), validation = "CV",
#'  predacc = "ALL")
#' svmcv1
#'
#' data(sponge2)
#' model <- species.richness ~ .
#' set.seed(1234)
#' svmcv1 <- svmcv(formula = model, sponge2[, -4], sponge[, 3], gamma = 0.01, cost = 3.5,
#' scale = TRUE, validation = "CV",  predacc = "VEcv")
#' svmcv1
#'
#' # For svm
#' model <- species.richness ~ .
#' set.seed(1234)
#' n <- 20 # number of iterations,60 to 100 is recommended.
#' VEcv <- NULL
#' for (i in 1:n) {
#' svmcv1 <- svmcv(formula = model, sponge2[, -4], sponge[, 3], gamma = 0.01, cost = 3.5,
#' scale = TRUE, validation = "CV",  predacc = "VEcv")
#' VEcv [i] <- svmcv1
#' }
#' plot(VEcv ~ c(1:n), xlab = "Iteration for svm", ylab = "VEcv (%)")
#' points(cumsum(VEcv) / c(1:n) ~ c(1:n), col = 2)
#' abline(h = mean(VEcv), col = 'blue', lwd = 2)
#'}
#'
#' @export
svmcv <- function (formula = NULL, trainxy, y, scale = TRUE, type = NULL, kernel = "radial", degree = 3, gamma = if (is.vector(trainxy)) 1 else 1 / ncol(trainxy), coef0 = 0, cost = 1, nu = 0.5, tolerance = 0.001, epsilon = 0.1, validation = "CV", cv.fold = 10, predacc = "VEcv", ...) {

  if (validation == "LOO") {idx <- 1:length(y)}
  if (validation == "CV")  {idx <- datasplit(y, k.fold = cv.fold)}

  # cross validation
  cv.pred <- NULL

  if (validation == "LOO") {
    for (i in 1 : length(y)) {
      data.dev <- trainxy[idx != i, , drop = FALSE]
      data.pred <- trainxy[idx == i, , drop = FALSE]
      svm1 <- e1071::svm(formula, data.dev, scale = scale, type = type, kernel = kernel, degree = degree, gamma = gamma, coef0 = coef0, cost = cost, nu = nu, tolerance = tolerance, epsilon = epsilon)
      cv.pred[idx == i] <- stats::predict(svm1, data.pred)
    }
  }

  if (validation == "CV") {
  for (i in 1 : cv.fold) {
    data.dev <- trainxy[idx != i, , drop = FALSE]
    data.pred <- trainxy[idx == i, , drop = FALSE]
    svm1 <- e1071::svm(formula, data.dev, scale = scale, type = type, kernel = kernel, degree = degree, gamma = gamma, coef0 = coef0, cost = cost, nu = nu, tolerance = tolerance, epsilon = epsilon)
    cv.pred[idx == i] <- stats::predict(svm1, data.pred)
  }
  }

  # predictive error and accuracy assessment
  if (predacc == "VEcv") {predictive.accuracy = spm::vecv(y, cv.pred)} else (
  if (predacc == "ALL") {predictive.accuracy = spm::pred.acc(y, cv.pred)} else (
  stop ("This measure is not supported in this version!")))
  predictive.accuracy
}
