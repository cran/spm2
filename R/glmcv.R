#' @title Cross validation, n-fold and leave-one-out for generalised linear models  ('glm')
#'
#' @description This function is a cross validation function
#' for 'glm' method in 'stats' package.
#'
#' @param formula a formula defining the response variable and predictive variables.
#' @param trainxy a dataframe contains predictive variables and the response
#' variable of point samples. The location information, longitude (long),
#' latitude (lat), need to be included in the 'trainx' for spatial predictive
#'  modeling.
#' @param y a vector of the response variable in the formula, that is, the left
#' part of the formula.
#' @param family a description of the error distribution and link function to
#' be used in the model. See '?glm' for details.
#' @param validation validation methods, include 'LOO': leave-one-out, and 'CV':
#' cross-validation.
#' @param cv.fold integer; number of folds in the cross-validation. if > 1,
#' then apply n-fold cross validation; the default is 10, i.e., 10-fold cross
#' validation that is recommended.
#' @param predacc can be either "VEcv" for vecv or "ALL" for all measures
#' in function pred.acc.
#' @param ... other arguments passed on to 'glm'.
#'
#' @return A list with the following components:
#'  me, rme, mae, rmae, mse, rmse, rrmse, vecv and e1; or vecv only
#' @note This function is largely based on 'rfcv' in 'randomForest' and
#' 'glm' in 'stats'.
#'
#' @references A. Liaw and M. Wiener (2002). Classification and Regression by
#' randomForest. R News 2(3), 18-22.
#'
#' @author Jin Li
#' @examples
#' \donttest{
#' library(spm)
#'
#' data(petrel)
#' gravel <- petrel[, c(1, 2, 6:9, 5)]
#' model <- log(gravel + 1) ~  lat +  bathy + I(long^3) + I(lat^2) + I(lat^3)
#' set.seed(1234)
#' glmcv1 <- glmcv(formula = model, gravel, log(gravel[, 7] +1), validation = "CV",
#'  predacc = "ALL")
#' glmcv1 # Since the default 'family' is used, it is actually a 'lm' model.
#'
#' data(sponge)
#' model <- sponge ~ easting + I(easting^2)
#' set.seed(1234)
#' glmcv1 <- glmcv(formula = model, sponge, sponge[, 3], family = poisson,
#' validation = "CV",  predacc = "ALL")
#' glmcv1
#'
#' # For glm
#' model <- gravel / 100 ~  lat +  bathy + I(long^3) + I(lat^2) + I(lat^3)
#' set.seed(1234)
#' n <- 20 # number of iterations,60 to 100 is recommended.
#' VEcv <- NULL
#' for (i in 1:n) {
#' glmcv1 <- glmcv(formula = model, gravel, gravel[, 7] / 100, family =
#' binomial(link=logit), validation = "CV", predacc = "VEcv")
#' VEcv [i] <- glmcv1
#' }
#' plot(VEcv ~ c(1:n), xlab = "Iteration for GLM", ylab = "VEcv (%)")
#' points(cumsum(VEcv) / c(1:n) ~ c(1:n), col = 2)
#' abline(h = mean(VEcv), col = 'blue', lwd = 2)
#'}
#'
#' @export
glmcv <- function (formula = NULL, trainxy, y, family = "gaussian", validation = "CV", cv.fold = 10, predacc = "VEcv", ...) {

  if (validation == "LOO") {idx <- 1:length(y)}
  if (validation == "CV")  {idx <- datasplit(y, k.fold = cv.fold)}

  # cross validation
  cv.pred <- NULL

  if (validation == "LOO") {
    for (i in 1 : length(y)) {
      data.dev <- trainxy[idx != i, , drop = FALSE]
      data.pred <- trainxy[idx == i, , drop = FALSE]
      glm1 <- stats::glm(formula, data.dev, family = family)
      cv.pred[idx == i] <- stats::predict(glm1, data.pred, type = "response")
    }
  }

  if (validation == "CV") {
  for (i in 1 : cv.fold) {
    data.dev <- trainxy[idx != i, , drop = FALSE]
    data.pred <- trainxy[idx == i, , drop = FALSE]
    glm1 <- stats::glm(formula, data.dev, family = family)
    cv.pred[idx == i] <- stats::predict(glm1, data.pred, type = "response")
  }
  }

  # predictive error and accuracy assessment
  if (predacc == "VEcv") {predictive.accuracy = spm::vecv(y, cv.pred)} else (
  if (predacc == "ALL") {predictive.accuracy = spm::pred.acc(y, cv.pred)} else (
  stop ("This measure is not supported in this version!")))
  predictive.accuracy
}
