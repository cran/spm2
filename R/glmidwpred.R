#' @title Generate spatial predictions using the hybrid method of generalised
#' linear models  ('glm') and inverse distance weighted ('IDW') ('glmidw')
#'
#' @description This function is for generating spatial predictions using the hybrid
#' method of 'glm' and 'idw' ('glmidw') (see reference #1).
#'
#' @param formula a formula defining the response variable and predictive variables
#'  for 'glm'.
#' @param longlat	a dataframe contains longitude and latitude of point samples. The
#'  location information must be named as 'long' and 'lat'.
#' @param trainxy a dataframe contains longitude (long), latitude (lat),
#' predictive variables and the response variable of point samples. That is,
#' the location information must be named as 'long' and 'lat'.
#' @param y a vector of the response variable in the formula, that is, the left
#' part of the formula.
#' @param longlatpredx	a dataframe contains longitude and latitude of point locations
#' (i.e., the centers of grids) to be predicted.
#' @param predx	a dataframe or matrix contains columns of predictive variables
#' for the grids to be predicted.
#' @param family a description of the error distribution and link function to
#' be used in the model. See '?glm' for details.
#' @param idp	 a numeric number specifying the inverse distance weighting power.
#' @param nmaxidw for a local predicting: the number of nearest observations that
#'  should be used for a prediction or simulation, where nearest is defined in
#'  terms of the space of the spatial locations. By default, 12 observations
#'  are used.
#' @param ... other arguments passed on to 'glm'.
#'
#' @return A dataframe of longitude, latitude, and predictions.
#'
#' @references Li, J., Alvarez, B., Siwabessy, J., Tran, M., Huang, Z.,
#' Przeslawski, R., Radke, L., Howard, F. and Nichol, S. (2017). "Application
#' of random forest, generalised linear model and their hybrid methods with
#' geostatistical techniques to count data: Predicting sponge species richness."
#' Environmental Modelling & Software 97: 112-129.
#'
#' Pebesma, E.J., 2004. Multivariable geostatistics in S: the gstat package.
#' Computers & Geosciences, 30: 683-691.
#'
#' @author Jin Li
#' @examples
#' \donttest{
#' library(spm)
#' data(petrel)
#' data(petrel.grid)
#'
#' gravel <- petrel[, c(1, 2, 6:9, 5)]
#' longlat <- petrel[, c(1, 2)]
#' model <- log(gravel + 1) ~  lat +  bathy + I(long^3) + I(lat^2) + I(lat^3)
#' y <- log(gravel[, 7] +1)
#'
#' glmidwpred1 <- glmidwpred(formula = model, longlat = longlat, trainxy =  gravel,
#' y = y, longlatpredx = petrel.grid[, c(1:2)], predx = petrel.grid, idp = 2,
#'  nmaxidw = 12)
#'  # Since the default 'family' is used, actually a 'lm' model is used.
#'
#' names(glmidwpred1)
#'
#' # Back transform 'glmidwpred$predictions' to generate the final predictions
#' glmidwpred1$predictions.bt <- exp(glmidwpred1$predictions) - 1
#' range(glmidwpred1$predictions.bt)
#'}
#'
#' @export
glmidwpred <- function (formula = NULL, longlat, trainxy, y, longlatpredx, predx, family = "gaussian", idp = 2, nmaxidw = 12, ...) {

    names(longlat) <- c("long", "lat")
    names(longlatpredx) <- c("long", "lat")

    # glm modeling
    glm1 <- stats::glm(formula, trainxy, family = family)

    # glm predictions
    pred.glm1 <- stats::predict(glm1, predx, type = "response")

    # the residuals of glm for idw
    data.dev1 <- longlat
    data.pred1 <- longlatpredx

    dev.glm1 <- stats::predict(glm1, trainxy, type="response")
    data.dev1$res1 <- y - dev.glm1

    # idw of the residuals
    gstat1 <- gstat::gstat(id = "res1", formula = res1 ~ 1, locations = ~ long + lat, data = data.dev1, set = list(idp = idp), nmax = nmaxidw)

    # idw predictions
    pred.idw1<- stats::predict(gstat1, data.pred1)

    predictions <- pred.idw1$res1.pred + pred.glm1
    glmidw.pred <- cbind(longlatpredx, predictions)
    glmidw.pred
}
