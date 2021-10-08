#' @title Generate spatial predictions using the hybrid method of support
#' vector machine ('svm') regression and 'krige' (svmkrige)
#'
#' @description This function is for generating spatial predictions using the hybrid method of
#' 'svm' and 'krige' (svmkrige).
#'
#' @param formula.svm a formula defining the response variable and predictive variables for 'svm'.
#' @param longlat	a dataframe contains longitude and latitude of point samples.
#' @param trainxy a dataframe contains longitude (long), latitude (lat),
#' predictive variables and the response variable of point samples.
#' @param predx	a dataframe or matrix contains columns of predictive variables
#' for the grids to be predicted.
#' @param y a vector of the response variable in the formula, that is, the left
#' part of the formula.
#' @param longlatpredx	a dataframe contains longitude and latitude of point locations
#' (i.e., the centers of grids) to be predicted.
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
#'  See '?svm' for details.
#' @param transformation transform the residuals of 'svm' to normalise the data;
#' can be "sqrt" for square root, "arcsine" for arcsine, "log" or "none"
#' for non transformation. By default, "none" is used.
#' @param delta numeric; to avoid log(0) in the log transformation. The default is 1.
#' @param formula.krige formula defining the response vector and (possible) regressor.
#' an object (i.e., 'variogram.formula') for 'variogram' or a formula for
#' 'krige'. see 'variogram' and 'krige' in 'gstat' for details.
#' @param vgm.args arguments for 'vgm', e.g. variogram model of response
#' variable and anisotropy parameters. see 'vgm' in 'gstat' for details.
#' By default, "Sph" is used.
#' @param anis anisotropy parameters: see notes 'vgm' in 'gstat' for details.
#' @param alpha direction in plane (x,y). see variogram in 'gstat' for details.
#' @param block block size. see 'krige' in 'gstat' for details.
#' @param beta for simple kriging. see 'krige' in 'gstat' for details.
#' @param nmaxkrige for a local predicting: the number of nearest observations that
#'  should be used for a prediction or simulation, where nearest is defined in
#'  terms of the space of the spatial locations. By default, 12 observations
#'  are used.
#' @param ... other arguments passed on to 'svm' and 'krige'.
#'
#' @return A dataframe of longitude, latitude, and predictions.
#'
#' @references Li, J., Potter, A., Huang, Z., and Heap, A. (2012). Predicting Seabed
#' Sand Content across the Australian Margin Using Machine Learning and Geostatistical
#'  Methods, Geoscience Australia, Record 2012/48, 115pp.
#'
#' Li, J., Heap, A., Potter, A., and Danilel, J.J. (2011). Predicting Seabed Mud Content
#' across the Australian Margin II: Performance of Machine Learning Methods and Their
#' Combination with Ordinary Kriging and Inverse Distance Squared, Geoscience Australia,
#' Record 2011/07, 69pp.
#'
#' David Meyer, Evgenia Dimitriadou, Kurt Hornik, Andreas Weingessel and Friedrich
#' Leisch (2020). e1071: Misc Functions of the Department of Statistics, Probability
#' Theory Group (Formerly: E1071), TU Wien. R package version 1.7-4.
#' https://CRAN.R-project.org/package=e1071.
#'
#' Pebesma, E.J., 2004. Multivariable geostatistics in S: the gstat package.
#' Computers & Geosciences, 30: 683-691.
#'
#' @author Jin Li
#' @examples
#' \donttest{
#' library(spm)
#'
#' data(petrel)
#' data(petrel.grid)
#'
#' gravel <- petrel[, c(1, 2, 6:9, 5)]
#' longlat <- petrel[, c(1, 2)]
#' model <- log(gravel + 1) ~  lat +  bathy + I(long^3) + I(lat^2) + I(lat^3)
#' y <- log(gravel[, 7] +1)
#'
#' svmkrigepred1 <- svmkrigepred(formula.svm = model, longlat = longlat, trainxy =  gravel,
#' predx = petrel.grid, y = y, longlatpredx = petrel.grid[, c(1:2)],
#' transformation = "none", formula.krige = res1 ~ 1, vgm.args = "Sph", nmaxkrige = 12)
#'
#' names(svmkrigepred1)
#'
#' # Back transform 'svmkrigepred$predictions' to generate the final predictions
#' svmkrige.predictions <- exp(svmkrigepred1$predictions) - 1
#' range(svmkrige.predictions)
#'}
#'
#' @export
svmkrigepred <- function (formula.svm = NULL, longlat, trainxy, predx, y, longlatpredx, scale = TRUE, type = NULL, kernel = "radial", degree = 3, gamma = if (is.vector(trainxy)) 1 else 1 / ncol(trainxy), coef0 = 0, cost = 1, nu = 0.5, tolerance = 0.001, epsilon = 0.1, transformation = "none", delta = 1, formula.krige = res1 ~ 1, vgm.args = c("Sph"), anis = c(0, 1), alpha = 0, block = 0, beta, nmaxkrige = 12, ...) {

  names(longlat) <- c("long", "lat")
  names(longlatpredx) <- c("long", "lat")

      # svm modeling
      svm1 <- e1071::svm(formula.svm, trainxy, scale = scale, type = type, kernel = kernel, degree = degree, gamma = gamma, coef0 = coef0, cost = cost, nu = nu, tolerance = tolerance, epsilon = epsilon)

      # svm predictions
      pred.svm1 <- stats::predict(svm1, predx, type = "response")

      # the residuals of svm for krige
      data.dev1 <- longlat
      data.pred1 <- longlatpredx
      dev.svm1 <- stats::predict(svm1, trainxy, type="response")
      res1 <- y - dev.svm1

      if (transformation == "none") {data.dev1$res1 = res1} else (
        if (transformation == "sqrt") {data.dev1$res1 = sqrt(res1 + abs(min(res1)))} else (
          if (transformation == "arcsine") {data.dev1$res1 = asin(sqrt((res1 + abs(min(res1))) / 100))} else (
            if (transformation == "log") {data.dev1$res1 = log(res1 + abs(min(res1)) + delta)} else (
              stop ("This transfromation is not supported in this version!")))))
        # The '+ abs(min(res1))' above is to set possible negative values to 0.

      # vgm of the residuals
      sp::coordinates(data.dev1) = ~ long + lat
      vgm1 <- gstat::variogram(object = formula.krige, data.dev1, alpha = alpha)
      model.1 <- gstat::fit.variogram(vgm1, gstat::vgm(mean(vgm1$gamma), vgm.args, mean(vgm1$dist), min(vgm1$gamma)/10, anis = anis))
      if (model.1$range[2] <= 0) (cat("A zero or negative range was fitted to variogram", "\n"))
      if (model.1$range[2] <= 0) (model.1$range[2] <- min(vgm1$dist))  # set negative range to be positive

      # krige predictions
      sp::coordinates(data.pred1) = ~long + lat
      pred.krige1 <- gstat::krige(formula = formula.krige, data.dev1, data.pred1, model = model.1, nmax=nmaxkrige, block = block, beta = beta)$var1.pred

      if (transformation == "none") {pred.krige = pred.krige1}
      if (transformation == "sqrt") {pred.krige = pred.krige1 ^ 2 - abs(min(res1))}
      if (transformation == "arcsine") {pred.krige = (sin(pred.krige1)) ^ 2 * 100 -  abs(min(res1))}
      if (transformation == "log") {pred.krige = exp(pred.krige1) - abs(min(res1)) - delta}

      predictions <- pred.krige + pred.svm1
      svmkrige.pred <- cbind(longlatpredx, predictions)
      svmkrige.pred
}

