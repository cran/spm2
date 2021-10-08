#' @title Cross validation, n-fold and leave-one-out for the hybrid method of support vector machine ('svm')
#'  regression and inverse distance weighted ('IDW') (svmidw)
#'
#' @description This function is a cross validation function for the hybrid
#' method of 'svm' regression and 'idw' using 'gstat' (svmidw), where the
#'  data splitting is based on a stratified random sampling method (see the
#'  'datasplit' function for details).
#'
#' @param formula a formula defining the response variable and predictive variables
#'  for 'svm'.
#' @param longlat	a dataframe contains longitude and latitude of point samples.
#' @param trainxy a dataframe contains longitude (long), latitude (lat),
#' predictive variables and the response variable of point samples. That is,
#' the location information must be named as 'long' and 'lat'.
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
#' @param idp	 a numeric number specifying the inverse distance weighting power.
#' @param nmaxidw for a local predicting: the number of nearest observations that
#'  should be used for a prediction or simulation, where nearest is defined in
#'  terms of the space of the spatial locations. By default, 12 observations
#'  are used.
#' @param validation validation methods, include 'LOO': leave-one-out, and 'CV':
#' cross-validation.
#' @param cv.fold integer; number of folds in the cross-validation. if > 1,
#' then apply n-fold cross validation; the default is 10, i.e., 10-fold cross
#' validation that is recommended.
#' @param predacc can be either "VEcv" for 'vecv' or "ALL" for all measures
#' in function pred.acc.
#' @param ... other arguments passed on to 'svm' and 'gstat'.
#'
#' @return A list with the following components:
#'  me, rme, mae, rmae, mse, rmse, rrmse, vecv and e1; or vecv only.
#'
#' @note This function is largely based on 'rfcv' in 'randomForest', 'idwcv'
#' in 'spm'and 'svm' in 'e1071'.
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
#' gravel <- petrel[, c(1, 2, 6:9, 5)]
#' longlat <- petrel[, c(1, 2)]
#' model <- log(gravel + 1) ~  lat +  bathy + I(long^3) + I(lat^2) + I(lat^3)
#' y <- log(gravel[, 7] +1)
#'
#' set.seed(1234)
#' svmidwcv1 <- svmidwcv(formula = model, longlat = longlat, trainxy =  gravel,
#' y = y, idp = 2, nmaxidw = 12, validation = "CV", predacc = "ALL")
#' svmidwcv1
#'
#' # svmidw for count data
#' data(sponge2)
#' model <- species.richness ~ . # use all predictive variables in the dataset
#' longlat <- sponge2[, 1:2]
#' set.seed(1234)
#' n <- 20 # number of iterations,60 to 100 is recommended.
#' VEcv <- NULL
#' for (i in 1:n) {
#'  svmidwcv1 <- svmidwcv(formula = model, longlat = longlat, trainxy = sponge2[, -4],
#'  y = sponge[, 3], gamma = 0.01,  cost = 3.5, scale = TRUE, idp = 2, nmaxidw = 12,
#'  validation = "CV", predacc = "VEcv")
#'  VEcv [i] <- svmidwcv1
#'  }
#'  plot(VEcv ~ c(1:n), xlab = "Iteration for svmidw", ylab = "VEcv (%)")
#'  points(cumsum(VEcv) / c(1:n) ~ c(1:n), col = 2)
#'  abline(h = mean(VEcv), col = 'blue', lwd = 2)
#'}
#'
#' @export
svmidwcv <- function (formula = NULL, longlat, trainxy, y, scale = TRUE, type = NULL, kernel = "radial", degree = 3, gamma = if (is.vector(trainxy)) 1 else 1 / ncol(trainxy), coef0 = 0, cost = 1, nu = 0.5, tolerance = 0.001, epsilon = 0.1, idp = 2, nmaxidw = 12, validation = "CV", cv.fold = 10, predacc = "VEcv", ...) {

  if (validation == "LOO") {idx <- 1:length(y)}
  if (validation == "CV")  {idx <- datasplit(y, k.fold = cv.fold)}

  names(longlat) <- c("long", "lat")

  # cross validation
  cv.pred <- NULL

  if (validation == "LOO") {
    for (i in 1 : length(y)) {
      data.dev <- trainxy[idx != i, , drop = FALSE]
      data.pred <- trainxy[idx == i, , drop = FALSE]

      # svm modeling
      svm1 <- e1071::svm(formula, data.dev, scale = scale, type = type, kernel = kernel, degree = degree, gamma = gamma, coef0 = coef0, cost = cost, nu = nu, tolerance = tolerance, epsilon = epsilon)

      # svm predictions
      pred.svm1 <- stats::predict(svm1, data.pred)

      # the residuals of svm for idw
      dev.svm1 <- stats::predict(svm1, data.dev)

      data.dev1 <- longlat[idx != i, , drop = FALSE] # for idw
      data.pred1 <- longlat[idx == i, , drop = FALSE] # for idw

      data.dev1$res1 <- y[idx != i] - dev.svm1

      # idw of the residuals
      gstat1 <- gstat::gstat(id = "res1", formula = res1 ~ 1, locations = ~ long + lat, data = data.dev1, set = list(idp = idp), nmax = nmaxidw)

      # idw predictions
      pred.idw1<- stats::predict(gstat1, data.pred1)

      cv.pred[idx == i] <- pred.idw1$res1.pred + pred.svm1
    }
  }

  if (validation == "CV") {
  for (i in 1 : cv.fold) {
    data.dev <- trainxy[idx != i, , drop = FALSE]
    data.pred <- trainxy[idx == i, , drop = FALSE]

    # svm modeling
    svm1 <- e1071::svm(formula, data.dev, scale = scale, type = type, kernel = kernel, degree = degree, gamma = gamma, coef0 = coef0, cost = cost, nu = nu, tolerance = tolerance, epsilon = epsilon)

    # svm predictions
    pred.svm1 <- stats::predict(svm1, data.pred)

    # the residuals of svm for idw
    dev.svm1 <- stats::predict(svm1, data.dev)

    data.dev1 <- longlat[idx != i, , drop = FALSE] # for idw
    data.pred1 <- longlat[idx == i, , drop = FALSE] # for idw

    data.dev1$res1 <- y[idx != i] - dev.svm1

    # idw of the residuals
    gstat1 <- gstat::gstat(id = "res1", formula = res1 ~ 1, locations = ~ long + lat, data = data.dev1, set = list(idp = idp), nmax = nmaxidw)

    # idw predictions
    pred.idw1<- stats::predict(gstat1, data.pred1)

    cv.pred[idx == i] <- pred.idw1$res1.pred + pred.svm1
    }
  }

  # predictive error and accuracy assessment
  if (predacc == "VEcv") {predictive.accuracy = spm::vecv(y, cv.pred)} else (
  if (predacc == "ALL") {predictive.accuracy = spm::pred.acc(y, cv.pred)} else (
  stop ("This measure is not supported in this version!")))
  predictive.accuracy
}
