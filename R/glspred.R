#' @title Generate spatial predictions using generalized least squares ('gls')
#'
#' @description This function is for generating spatial predictions using  'gls' method
#' in 'nlme' package.
#'
#' @param model a formula defining the response variable and predictive variables.
#' @param trainxy a dataframe contains longitude (long), latitude (lat),
#' predictive variables and the response variable of point samples. That is,
#' the location information must be names as 'long' and 'lat'.
#' @param longlatpredx	a dataframe contains longitude and latitude of point
#' locations (i.e., the centers of grids) to be predicted, need to be named as
#'  'long' and 'lat'.
#' @param predx	a dataframe or matrix contains columns of predictive variables
#' for the grids to be predicted.
#' @param weights describing the within-group heteroscedasticity structure. Defaults
#'  to "NULL", corresponding to homoscedastic errors. See '?gls' in 'nlme'
#' for details.
#' @param corr.args arguments for 'correlation' in 'gls'. See '?corClasses' in 'nlme'
#' for details. By default, "NULL" is used. When "NULL" is used,
#' then 'gls' is actually performing 'lm'.
#' @param ... other arguments passed on to 'gls'.
#'
#' @return A dataframe of longitude, latitude and predictions.
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
#' data(petrel.grid)
#'
#' gravel <- petrel[, c(1, 2, 6:9, 5)]
#' range1 <- 0.8
#' nugget1 <- 0.5
#' model <- log(gravel + 1) ~  long + lat +  bathy + dist + I(long^2) + I(lat^2)
#' + I(lat^3) + I(bathy^2) + I(bathy^3) + I(dist^2) + I(dist^3) + I(relief^2) + I(relief^3)
#'
#' glspred1 <- glspred(model = model, trainxy = gravel,
#'  longlatpredx = petrel.grid[, c(1:2)], predx = petrel.grid,
#'  corr.args = corSpher(c(range1, nugget1), form = ~ lat + long, nugget = T))
#'
#' names(glspred1)
#'
#' # Back transform 'glspred1$predictions' to generate the final predictions
#' gls.predictions <- exp(glspred1$predictions) - 1
#' range(gls.predictions)
#'}
#'
#' @export
glspred <- function (model = var1 ~ 1, trainxy, longlatpredx, predx, corr.args = NULL, weights = NULL, ...) {
      gls1 <- nlme::gls(model, trainxy, correlation = corr.args, weights = weights)
      predictions <- stats::predict(gls1, predx, type = "response")

      gls.pred <- cbind(longlatpredx, predictions)
      gls.pred
}
