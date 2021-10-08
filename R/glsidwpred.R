#' @title Generate spatial predictions using the hybrid method of generalized
#'  least squares ('gls')  and inverse distance weighted ('IDW') ('glsidw')
#'
#' @description This function is for generating spatial predictions using the hybrid
#' method of 'gls' and 'idw' ('glsidw') (see reference #1).
#'
#' @param model a formula defining the response variable and predictive variables.
#' @param longlat	a dataframe contains longitude and latitude of point samples.
#' @param trainxy a dataframe contains longitude (long), latitude (lat),
#' predictive variables and the response variable of point samples. That is,
#' the location information must be names as 'long' and 'lat'.
#' @param y a vector of the response variable in the formula, that is, the left
#' part of the formula.
#' @param longlatpredx	a dataframe contains longitude and latitude of point locations
#' (i.e., the centers of grids) to be predicted. The location information must be
#' named as 'long' and 'lat'.
#' @param predx	a dataframe or matrix contains columns of predictive variables
#' for the grids to be predicted.
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
#' @param idp	 a numeric number specifying the inverse distance weighting power.
#' @param nmaxidw for a local predicting: the number of nearest observations that
#'  should be used for a prediction or simulation, where nearest is defined in
#'  terms of the space of the spatial locations. By default, 12 observations
#'  are used.
#' @param ... other arguments passed on to 'gls' and 'gstat'.
#'
#' @return A dataframe of longitude, latitude, and predictions.
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
#' data(petrel.grid)
#'
#' gravel <- petrel[, c(1, 2, 6:9, 5)]
#' longlat <- petrel[, c(1, 2)]
#' range1 <- 0.8
#' nugget1 <- 0.5
#' model <- log(gravel + 1) ~  long + lat +  bathy + dist + I(long^2) + I(lat^2) +
#' I(lat^3) + I(bathy^2) + I(bathy^3) + I(dist^2) + I(dist^3) + I(relief^2) + I(relief^3)
#' y <- log(gravel[, 7] +1)
#'
#' glsidwpred1 <- glsidwpred(model = model, longlat = longlat, trainxy = gravel,
#' y = y, longlatpredx = petrel.grid[, c(1:2)], predx = petrel.grid,
#'  idp = 2, nmaxidw = 12, corr.args = corSpher(c(range1, nugget1),
#'  form = ~ lat + long, nugget = T))
#'
#' names(glsidwpred1)
#'
#' # Back transform 'glsidwpred$predictions' to generate the final predictions
#' glsidw.predictions <- exp(glsidwpred1$predictions) - 1
#' range(glsidw.predictions)
#'}
#'
#' @export
glsidwpred <- function (model = var1 ~ 1, longlat, trainxy, y, longlatpredx, predx, corr.args = NULL, weights = NULL, idp = 2, nmaxidw = 12, ...) {
      # gls modeling
      gls1 <- nlme::gls(model, trainxy, correlation = corr.args, weights = weights)

      # gls predictions
      pred.gls1 <- stats::predict(gls1, predx, type = "response")

      # the residuals of gls for idw
      data.dev1 <- longlat
      data.pred1 <- longlatpredx
      dev.gls1 <- stats::predict(gls1, trainxy, type="response")
      data.dev1$res1 <- y - dev.gls1

      # idw of the residuals
      gstat1 <- gstat::gstat(id = "res1", formula = res1 ~ 1, locations = ~ long + lat, data = data.dev1, set = list(idp = idp), nmax = nmaxidw)

      # idw predictions
      pred.idw1<- stats::predict(gstat1, data.pred1)

      predictions <- pred.idw1$res1.pred + pred.gls1
      glsidw.pred <- cbind(longlatpredx, predictions)
      glsidw.pred
}
