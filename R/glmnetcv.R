#' @title Cross validation, n-fold and leave-one-out, for 'glmnet' in 'glmnet' package
#'
#' @description This function is a cross validation function
#' for 'glmnet' method in 'glmnet' package.
#'
#' @param trainx a matrix contains predictive variables of point samples.
#' The location information, longitude (long), latitude (lat), need to be included
#' in the 'trainx' for spatial predictive modelling.
#' @param y a vector of the response variable in the formula, that is, the left
#' part of the formula.
#' @param family a description of the error distribution and link function to
#' be used in the model. See '?glmnet' for details.
#' @param alpha, an elasticnet mixing parameter, with $0 <= alpha <= 1$.
#' See '?glmnet' for details.
#' @param relax, if TRUE then for each active set in the path of solutions,
#' the model is refit without any regularization. See '?glmnet' for more information.
#' @param type.measure, loss to use for cross-validation. The default is
#' type.measure="mse". See '?cv.glmnet' for more information.
#' @param validation validation methods, include 'LOO': leave-one-out, and 'CV':
#' cross-validation.
#' @param cv.fold integer; number of folds in the cross-validation. if > 1,
#' then apply n-fold cross validation; the default is 10, i.e., 10-fold cross
#' validation that is recommended.
#' @param predacc can be either "VEcv" for vecv or "ALL" for all measures
#' in function pred.acc.
#' @param ... other arguments passed on to 'fields'.
#'
#' @return A list with the following components:
#'  me, rme, mae, rmae, mse, rmse, rrmse, vecv and e1; or vecv only
#' @note This function is largely based on 'glmcv' in this 'spm2' package.
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
#' x <- as.matrix(petrel[, c(1, 2, 6:9)])
#' y <- log(petrel[, 5] + 1)
#' set.seed(1234)
#' glmnetcv1 <- glmnetcv(x, y, validation = "CV",  predacc = "ALL")
#' glmnetcv1
#'
#' data(sponge)
#' x <- as.matrix(cbind(sponge$easting, sponge$easting^2))
#' set.seed(1234)
#' glmnetcv1 <- glmnetcv(x, sponge[, 3], family = poisson, validation = "CV",
#' predacc = "ALL")
#' glmnetcv1
#'
#' # For glmnet with gaussian
#' x <- as.matrix(petrel[, c(1, 2, 6:9)])
#' y <- log(petrel[, 5] + 1)
#' set.seed(1234)
#' n <- 20 # number of iterations,60 to 100 is recommended.
#' VEcv <- NULL
#' for (i in 1:n) {
#'   glmnetcv1 <- glmnetcv(x, y, validation = "CV", predacc = "VEcv")
#'   VEcv [i] <- glmnetcv1
#'  }
#' plot(VEcv ~ c(1:n), xlab = "Iteration for glmnet", ylab = "VEcv (%)")
#' points(cumsum(VEcv) / c(1:n) ~ c(1:n), col = 2)
#' abline(h = mean(VEcv), col = 'blue', lwd = 2)
#'
#' # For glmnet with binomial
#' x <- as.matrix(cbind(petrel[, c(2, 6)], petrel$long^3, petrel$lat^2, petrel$lat^3))
#' set.seed(1234)
#' n <- 20 # number of iterations,60 to 100 is recommended.
#' VEcv <- NULL
#' for (i in 1:n) {
#' glmnetcv1 <- glmnetcv(x, petrel[, 5] / 100, family = binomial(link=logit),
#' validation = "CV", predacc = "VEcv")
#' VEcv [i] <- glmnetcv1
#' }
#' plot(VEcv ~ c(1:n), xlab = "Iteration for glmnet", ylab = "VEcv (%)")
#' points(cumsum(VEcv) / c(1:n) ~ c(1:n), col = 2)
#' abline(h = mean(VEcv), col = 'blue', lwd = 2)
#'}
#'
#' @export
glmnetcv <- function (trainx, y, family = "gaussian", alpha = 0.5, relax = FALSE, type.measure = "mse", validation = "CV", cv.fold = 10, predacc = "VEcv", ...) {

  if (validation == "LOO") {idx <- 1:length(y)}
  if (validation == "CV")  {idx <- spm2::datasplit(y, k.fold = cv.fold)}

  # cross validation
  cv.pred <- NULL

  if (validation == "LOO") {
    for (i in 1 : length(y)) {
      data.dev <- trainx[idx != i, , drop = FALSE]
      data.pred <- trainx[idx == i, , drop = FALSE]
      y1 <- y[idx != i, drop = FALSE]
      glmnet1 <- glmnet::glmnet(data.dev, y1, family = family, alpha = alpha, relax = relax)
      enet.cv <- glmnet::cv.glmnet(data.dev, y1, alpha = alpha, type.measure = type.measure)
      cv.pred[idx == i] <- stats::predict(glmnet1, data.pred, s = enet.cv$lambda.min, type = "response")
    }
  }

  if (validation == "CV") {
  for (i in 1 : cv.fold) {
    data.dev <- trainx[idx != i, , drop = FALSE]
    data.pred <- trainx[idx == i, , drop = FALSE]
    y1 <- y[idx != i, drop = FALSE]
    glmnet1 <- glmnet::glmnet(data.dev, y1, family = family, alpha = alpha, relax = relax)
    enet.cv <- glmnet::cv.glmnet(data.dev, y1, alpha = alpha, type.measure = type.measure)
    cv.pred[idx == i] <- stats::predict(glmnet1, data.pred, s = enet.cv$lambda.min, type = "response")
  }
  }

  # predicitve error and accuracy assessment
  if (predacc == "VEcv") {predictive.accuracy = spm::vecv(y, cv.pred)} else (
  if (predacc == "ALL") {predictive.accuracy = spm::pred.acc(y, cv.pred)} else (
  stop ("This measure is not supported in this version!")))
  predictive.accuracy
}
