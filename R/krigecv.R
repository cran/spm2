#' @title Cross validation, n-fold and leave-one-out for kriging methods ('krige')
#'
#' @description This function is a cross validation function
#' for kriging methods ('krige') in 'gstat'.
#'
#' @param longlat	a dataframe contains longitude and latitude of point samples.
#' @param trainy a vector of response, must have length equal to the number of
#' rows in trainx.
#' @param trainpredx a dataframe contains predictive variables of point samples.
#' If longitude and latitude are going to be used as predictive variables, they
#' should also be included but they should be named in names other than 'long' and 'lat'.
#' @param validation validation methods, include 'LOO': leave-one-out, and 'CV':
#' cross-validation.
#' @param cv.fold integer; number of folds in the cross-validation. if > 1,
#' then apply n-fold cross validation; the default is 10, i.e., 10-fold cross
#' validation that is recommended.
#' @param nmax for local kriging: the number of nearest observations that
#' should be used for a kriging prediction or simulation, where nearest is
#' defined in terms of the space of the spatial locations. By default, 12
#' observations are used.
#' @param transformation transform response variable to normalise the data;
#' can be "sqrt" for square root, "arcsine" for arcsine, "log" or "none"
#' for non transformation. By default, "none" is used.
#' @param delta numeric; to avoid 'log(0)' in "log" transformation. The default is 1.
#' @param formula formula defining response vector and (possible) regressor.
#' an object (i.e., 'variogram.formula') for 'variogram' or a formula for
#' 'krige'. see 'variogram' and 'krige' in the 'gstat' package for details.
#' @param vgm.args arguments for 'vgm', e.g. variogram model of response
#' variable and anisotropy parameters. see 'vgm' in the 'gstat' package for details.
#' By default, "Sph" is used.
#' @param anis anisotropy parameters: see notes 'vgm' in the 'gstat' package for details.
#' @param alpha direction in plane (x,y). see variogram in the 'gstat' package for details.
#' @param block block size. see 'krige' in the 'gstat' package for details.
#' @param beta for simple kriging. see 'krige' in the 'gstat' package for details.
#' @param predacc can be either "VEcv" for vecv or "ALL" for all measures
#' in function pred.acc.
#' @param ... other arguments passed on to the function 'gstat'.
#'
#' @return A list with the following components:
#'  me, rme, mae, rmae, mse, rmse, rrmse, vecv and e1; or vecv only
#' @note This function is largely based on rfcv in 'randomForest' and
#' some functions in 'library(gstat)'.  When 'A zero or negative range was fitted
#' to variogram' occurs, to allow 'gstat' running, the range was set to be positive by
#' using 'min(vgm1$dist)'. In this case, caution should be taken in applying this
#' method. If it still occurs for 'okpred' function, different method may need to be
#' used.
#'
#' @references Li, J., 2013. Predictive Modelling Using Random Forest and Its
#' Hybrid Methods with Geostatistical Techniques in Marine Environmental
#' Geosciences, In: Christen, P., Kennedy, P., Liu, L., Ong, K.-L., Stranieri,
#' A., Zhao, Y. (Eds.), The proceedings of the Eleventh Australasian Data
#' Mining Conference (AusDM 2013), Canberra, Australia, 13-15 November 2013.
#' Conferences in Research and Practice in Information Technology, Vol. 146.
#'
#' A. Liaw and M. Wiener (2002). Classification and Regression by
#' randomForest. R News 2(3), 18-22.
#'
#' Pebesma, E.J., 2004. Multivariable geostatistics in S: the gstat
#' package. Computers & Geosciences, 30: 683-691.
#'
#' @author Jin Li
#' @examples
#' \donttest{
#' library(sp)
#' library(spm)
#' data(swmud)
#' data(petrel)
#'
#' set.seed(1234)
#' okcv1 <- krigecv(longlat = swmud[, c(1,2)], trainy = swmud[, 3], nmax = 7, transformation =
#' "arcsine", vgm.args = ("Sph"), predacc = "VEcv")
#' okcv1
#'
#' set.seed(1234)
#' skcv1 <- krigecv(longlat = swmud[, c(1,2)], trainy = swmud[, 3], nmax = 7, transformation =
#' "arcsine", vgm.args = ("Sph"), predacc = "VEcv", beta = mean(swmud[, 3]))
#' skcv1
#'
#' set.seed(1234)
#' ukcv1 <- krigecv(longlat = swmud[, c(1,2)], trainy = swmud[, 3], nmax = 7, transformation =
#' "arcsine", formula = var1 ~ long + lat, vgm.args = ("Sph"), predacc = "VEcv")
#' ukcv1
#'
#' set.seed(1234)
#' okcv2 <- krigecv(longlat = swmud[, c(1,2)], trainy = swmud[, 3], validation = "LOO", nmax = 7,
#' transformation = "arcsine", vgm.args = ("Sph"), predacc = "ALL")
#' okcv2
#'
#' set.seed(1234)
#' n <- 20 # number of iterations,60 to 100 is recommended.
#' VEcv <- NULL
#' for (i in 1:n) {
#' okcv1 <- krigecv(longlat = petrel[, c(1,2)], trainy = petrel[, 5], nmax = 12,
#' transformation = "arcsine", predacc = "VEcv")
#' VEcv [i] <- okcv1
#' }
#' plot(VEcv ~ c(1:n), xlab = "Iteration for OK", ylab = "VEcv (%)")
#' points(cumsum(VEcv) / c(1:n) ~ c(1:n), col = 2)
#' abline(h = mean(VEcv), col = 'blue', lwd = 2)
#'
#' set.seed(1234)
#' n <- 20 # number of iterations, 60 to 100 is recommended.
#' measures <- NULL
#' for (i in 1:n) {
#' okcv1 <- krigecv(longlat = petrel[, c(1,2)], trainy = petrel[, 3], nmax = 12, transformation =
#' "arcsine", predacc = "ALL")
#' measures <- rbind(measures, okcv1$vecv)
#' }
#' plot(measures ~ c(1:n), xlab = "Iteration for OK", ylab = "VEcv (%)")
#' points(cumsum(measures) / c(1:n) ~ c(1:n), col = 2)
#' abline(h = mean(measures), col = 'blue', lwd = 2)
#' }
#'
#' @export
krigecv <- function (longlat, trainy, trainpredx = NULL, validation = "CV", cv.fold = 10, nmax = 12, transformation = "none", delta = 1, formula = var1 ~ 1, vgm.args = ("Sph"), anis = c(0, 1), alpha = 0, block = 0, beta, predacc = "VEcv", ...) {
  classRF <- is.factor(trainy)

  names(longlat) <- c("long", "lat")

  if (is.null(trainpredx)) {trainx <- longlat} else (trainx = cbind(longlat, trainpredx))

  if (transformation == "none") {trainy1 = trainy} else (
  if (transformation == "sqrt") {trainy1 = sqrt(trainy)} else (
  if (transformation == "arcsine") {trainy1 = asin(sqrt(trainy / 100))} else (
  if (transformation == "log") {trainy1 = log(trainy + delta)} else (
  stop ("This transfromation is not supported in this version!")))))
  if (classRF) {
    stop ("This function is not for categorical response variable")
  } else (
  if (validation == "LOO") {idx <- 1:length(trainy)} else (
  if (validation == "CV")  {idx <- datasplit(trainy, k.fold = cv.fold)} else (stop ("This validation method is not supported in this version!"))))

  # cross validation
  cv.pred <- NULL
  cv.pred1 <- NULL

  if (validation == "LOO") {
    for (i in 1 : length(trainy1)) {
      data.dev <- trainx[idx != i, , drop = FALSE]
      data.pred <- trainx[idx == i, , drop = FALSE]
      data.dev$var1 <- trainy1[idx != i]
      sp::coordinates(data.dev) = ~ long + lat
      vgm1 <- gstat::variogram(object = formula, data.dev, alpha = alpha)
      #model.1 <- gstat::fit.variogram(vgm1, gstat::vgm(vgm.args, anis = anis))
      model.1 <- gstat::fit.variogram(vgm1, gstat::vgm(mean(vgm1$gamma), vgm.args, mean(vgm1$dist), min(vgm1$gamma)/10, anis = anis))
      if (model.1$range[2] <= 0) (cat("A zero or negative range was fitted to variogram. To allow gstat running, the range was set to be positive by using min(vgm1$dist). ", "\n"))
      if (model.1$range[2] <= 0) (model.1$range[2] <- min(vgm1$dist)) # set negative range to be positive
      sp::coordinates(data.pred) = ~long + lat
      cv.pred1[idx == i] <- gstat::krige(formula = formula, data.dev, data.pred, model = model.1, nmax = nmax, block = block, beta = beta)$var1.pred
    }
  }

  if (validation == "CV") {
  for (i in 1 : cv.fold) {
    data.dev <- trainx[idx != i, , drop = FALSE]
    data.pred <- trainx[idx == i, , drop = FALSE]
    data.dev$var1 <- trainy1[idx != i]
    sp::coordinates(data.dev) = ~ long + lat
    vgm1 <- gstat::variogram(object = formula, data.dev, alpha = alpha)
    #model.1 <- gstat::fit.variogram(vgm1, gstat::vgm(vgm.args, anis = anis))
    model.1 <- gstat::fit.variogram(vgm1, gstat::vgm(mean(vgm1$gamma), vgm.args, mean(vgm1$dist), min(vgm1$gamma)/10, anis = anis))
    if (model.1$range[2] <= 0) (cat("A zero or negative range was fitted to variogram. To allow gstat running, the range was set to be positive by using min(vgm1$dist). ", "\n"))
    if (model.1$range[2] <= 0) (model.1$range[2] <- min(vgm1$dist)) # set negative range to be positive
    sp::coordinates(data.pred) = ~long + lat
    cv.pred1[idx == i] <- gstat::krige(formula = formula, data.dev, data.pred, model = model.1, nmax = nmax, block = block, beta = beta)$var1.pred
  }
  }

  if (transformation == "none") {cv.pred = cv.pred1}
  if (transformation == "sqrt") {cv.pred = cv.pred1 ^ 2}
  if (transformation == "arcsine") {cv.pred = (sin(cv.pred1)) ^ 2 * 100}
  if (transformation == "log") {cv.pred = exp(cv.pred1)-delta}
  # predicitve error and accuracy assessment
  if (predacc == "VEcv") {predictive.accuracy = spm::vecv(trainy, cv.pred)} else (
  if (predacc == "ALL") {predictive.accuracy = spm::pred.acc(trainy, cv.pred)} else (
  stop ("This measure is not supported in this version!")))
  predictive.accuracy
}
