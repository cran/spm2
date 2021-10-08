#' @title Generate spatial predictions using generalised linear models  ('glm')
#'
#' @description This function is for generating spatial predictions using 'glm' method
#' in 'stats' package.
#'
#' @param formula a formula defining the response variable and predictive variables.
#' @param trainxy a dataframe contains predictive variables and the response
#' variable of point samples. The location information, longitude (long),
#' latitude (lat), need to be included in the 'trainx' for spatial predictive
#'  modeling, need to be named as 'long' and 'lat'.
#' @param family a description of the error distribution and link function to
#' be used in the model. See '?glm' for details.
#' @param longlatpredx	a dataframe contains longitude and latitude of point
#' locations (i.e., the centers of grids) to be predicted.
#' @param predx	a dataframe or matrix contains columns of predictive variables
#' for the grids to be predicted.
#' @param ... other arguments passed on to 'glm'.
#'
#' @return A dataframe of longitude, latitude and predictions.
#'
#' @author Jin Li
#' @examples
#' \donttest{
#' library(spm)
#' data(petrel)
#' data(petrel.grid)
#'
#' gravel <- petrel[, c(1, 2, 6:9, 5)]
#' model <- log(gravel + 1) ~  lat +  bathy + I(long^3) + I(lat^2) + I(lat^3)
#'
#' glmpred1 <- glmpred(formula = model, trainxy = gravel,
#' longlatpredx = petrel.grid[, c(1:2)], predx = petrel.grid)
#'
#' names(glmpred1)
#'
#' # Back transform 'glmpred1$pred.glm1' to generate the final predictions
#' glm.predictions <- exp(glmpred1$pred.glm1) - 1
#' range(glm.predictions)
#'}
#'
#' @export
glmpred <- function (formula = NULL, trainxy, longlatpredx, predx, family = "gaussian", ...) {
    glm1 <- stats::glm(formula, trainxy, family = family)
    pred.glm1 <- stats::predict(glm1, predx, type = "response")

    glm.pred <- cbind(longlatpredx, pred.glm1)
    glm.pred
}
