#' @title Cross validation, n-fold and leave-one-out for the hybrid methods of
#'  generalized boosted regression modeling ('gbm'), 'kriging' and
#'  inverse distance weighted ('IDW').
#'
#' @description This function is a cross validation function for 38 hybrid
#'  methods of 'gbm', 'kriging' and 'IDW', including the average of 'gbmkrige'
#' and 'gbmidw' ('gbmkrigegbmidw') and  the average of 'gbm', 'gbmkrige' and 'gbmidw'
#' ('gbmgbmkrigegbmidw'), where 'kriging' methods include ordinary kriging
#'  ('OK'), simple kriging ('SK'), block 'OK' ('BOK') and block 'SK'('BSK')
#' and 'IDW' also covers 'NN' and 'KNN'. The data splitting is based on a stratified
#' random sampling method (see the 'datasplit' function for details).
#'
#' @param longlat	a dataframe contains longitude and latitude of point samples.
#' @param trainx a dataframe or matrix contains columns of predictive variables.
#' @param trainy a vector of the response variable.
#' @param var.monotone an optional vector, the same length as the number of
#' predictors, indicating which variables have a monotone increasing (+1),
#' decreasing (-1), or arbitrary (0) relationship with the outcome. By default,
#' a vector of 0 is used.
#' @param family either a character string specifying the name of the distribution to
#' use or a list with a component name specifying the distribution and any
#' additional parameters needed. See gbm for details. By default, "gaussian" is
#' used.
#' @param n.trees the total number of trees to fit. This is equivalent to the
#' number of iterations and the number of basis functions in the additive
#' expansion. By default, 3000 is used.
#' @param learning.rate a shrinkage parameter applied to each tree in the
#' expansion. Also known as step-size reduction.
#' @param interaction.depth the maximum depth of variable interactions.
#' 1 implies an additive model, 2 implies a model with up to 2-way
#' interactions, etc. By default, 2 is used.
#' @param bag.fraction the fraction of the training set observations randomly
#' selected to propose the next tree in the expansion. By default, 0.5 is used.
#' @param train.fraction The first train.fraction * nrows(data) observations
#' are used to fit the gbm and the remainder are used for computing
#' out-of-sample estimates of the loss function.
#' @param n.minobsinnode minimum number of observations in the trees terminal
#' nodes. Note that this is the actual number of observations not the total
#' weight. By default, 10 is used.
#' @param transformation transform the residuals of 'gbm' to normalize the data for 'krige';
#' can be "sqrt" for square root, "arcsine" for arcsine, "log" or "none"
#' for non transformation. By default, "none" is used.
#' @param weights an optional vector of weights to be used in the fitting
#' process. Must be positive but do not need to be normalized.
#' If keep.data = FALSE in the initial call to gbm then it is the user's
#' responsibility to resupply the weights to gbm.more. By default, a vector of
#' 1 is used.
#' @param keep.data a logical variable indicating whether to keep the data and
#' an index of the data stored with the object. Keeping the data and index
#' makes subsequent calls to gbm.more faster at the cost of storing an extra
#' copy of the dataset. By default, 'FALSE' is used.
#' @param verbose If TRUE, gbm will print out progress and performance
#' indicators. By default, 'TRUE' is used.
#' @param delta numeric; to avoid log(0) in the log transformation. The default is 1.
#' @param formula formula defining the response vector and (possible) regressor.
#' an object (i.e., 'variogram.formula') for 'variogram' or a formula for
#' 'krige'. see 'variogram' and 'krige' in 'gstat' for details. The default is
#' 'formula = res1 ~ 1'.
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
#' @param idp	 a numeric number specifying the inverse distance weighting power.
#' @param nmaxidw for a local predicting: the number of nearest observations that
#'  should be used for a prediction or simulation, where nearest is defined in
#'  terms of the space of the spatial locations. By default, 12 observations
#'  are used.
#' @param hybrid.parameter the default is 2 that is for 'gbmkrigegbmidw';
#' for 'gbmgbmkrigegbmidw', it needs to be 3.
#' @param lambda, ranging from 0 to 2; the default is 1 for 'gbmkrigegbmidw'
#' and 'gbmgbmkrigegbmidw'; and if it is < 1, more weight is placed on 'krige',
#' otherwise more weight is placed on 'idw'; and if it is 0, 'idw' is not
#' considered and the resultant methods is 'gbmkrige' when the default
#' 'hybrid.parameter' is used; and if it is 2, then the resultant method is
#' 'gbmidw' when the default 'hybrid.parameter' is used.
#' @param validation validation methods, include 'LOO': leave-one-out, and 'CV':
#' cross-validation.
#' @param cv.fold integer; number of folds in the cross-validation. if > 1,
#' then apply n-fold cross validation; the default is 10, i.e., 10-fold cross
#' validation that is recommended.
#' @param predacc can be either "VEcv" for vecv or "ALL" for all measures
#' in function pred.acc.
#' @param n.cores The number of CPU cores to use. See gbm for details. By
#' default, 6 is used.
#' @param ... other arguments passed on to 'randomForest', 'krige' and 'gstat'.
#'
#' @return A list with the following components:
#'  me, rme, mae, rmae, mse, rmse, rrmse, vecv and e1; or vecv only
#'
#' @note This function is largely based on 'gbmcv' in 'spm', and 'krigecv'
#' in 'spm2'.
#'
#' @references Li, J. (2022). Spatial Predictive Modeling with R. Boca Raton,
#' Chapman and Hall/CRC.
#'
#' Li, J., Potter, A., Huang, Z., and Heap, A. (2012). Predicting Seabed
#' Sand Content across the Australian Margin Using Machine Learning and Geostatistical
#'  Methods, Geoscience Australia, Record 2012/48, 115pp.
#'
#' Li, J., Heap, A., Potter, A., and Danilel, J.J. (2011). Predicting Seabed Mud Content
#' across the Australian Margin II: Performance of Machine Learning Methods and Their
#' Combination with Ordinary Kriging and Inverse Distance Squared, Geoscience Australia,
#' Record 2011/07, 69pp.
#'
#' Greg Ridgeway with contributions from others (2015). gbm: Generalized
#' Boosted Regression Models. R package version 2.1.1.
#' https://CRAN.R-project.org/package=gbm
#'
#' Pebesma, E.J., 2004. Multivariable geostatistics in S: the gstat package.
#' Computers & Geosciences, 30: 683-691.
#'
#' @author Jin Li
#' @examples
#'
#' \donttest{
#' library(spm)
#' # gbmgbmokgbmidw
#' data(sponge)
#' longlat <- sponge[, 1:2]
#' set.seed(1234)
#' gbmgbmkrigegbmidwcv1 <- gbmkrigeidwcv(longlat = longlat,
#' trainx = sponge[, -3], trainy = sponge[, 3], family = "poisson", interaction.depth = 3,
#' transformation = "none", formula = res1 ~ 1, vgm.args = "Sph",
#' nmaxkrige = 12, idp = 2, nmaxidw = 12, hybrid.parameter = 3, validation = "CV",
#' predacc = "ALL", n.cores = 2)
#' gbmgbmkrigegbmidwcv1
#'
#' # gbmokgbmidw for count data
#' data(sponge)
#' longlat <- sponge2[, 1:2]
#' y = sponge[, 3]
#' trainx = sponge[, -3]
#' set.seed(1234)
#' n <- 20 # number of iterations,60 to 100 is recommended.
#' VEcv <- NULL
#' for (i in 1:n) {
#'   gbmkrigegbmidwcv1 <- gbmkrigeidwcv(longlat = longlat,
#'   trainx = trainx, trainy = y, family = "poisson", interaction.depth = 3,
#'   transformation = "none", formula = res1 ~ 1, vgm.args = ("Sph"),
#'   nmaxkrige = 12, idp = 2, nmaxidw = 12, hybrid.parameter = 2, validation = "CV",
#'   predacc = "VEcv", n.cores = 2)
#'   VEcv [i] <- gbmkrigegbmidwcv1
#'  }
#'  plot(VEcv ~ c(1:n), xlab = "Iteration for gbmokgbmidw", ylab = "VEcv (%)")
#'  points(cumsum(VEcv) / c(1:n) ~ c(1:n), col = 2)
#'  abline(h = mean(VEcv), col = 'blue', lwd = 2)
#'}
#'
#' @export
gbmkrigeidwcv <- function (longlat, trainx, trainy,
                           var.monotone = rep(0, ncol(trainx)),
                           family = "gaussian",
                           n.trees = 3000,          # default number of trees
                           learning.rate = 0.001,
                           interaction.depth = 2,
                           bag.fraction = 0.5,
                           train.fraction = 1.0,
                           n.minobsinnode = 10,
                           transformation = "none",
                           weights = rep(1, nrow(trainx)),   # by default set equal
                           keep.data = FALSE,
                           verbose = TRUE,
                           delta = 1,
                           formula = res1 ~ 1,
                           vgm.args = "Sph",
                           anis = c(0, 1),
                           alpha = 0,
                           block = 0,
                           beta,
                           nmaxkrige = 12,
                           idp = 2,
                           nmaxidw = 12,
                           hybrid.parameter = 2,
                           lambda = 1,
                           validation = "CV",
                           cv.fold = 10,
                           predacc = "VEcv",
                           n.cores = 6, ...) {

  if (validation == "LOO") {idx <- 1:length(trainy)}
  if (validation == "CV")  {idx <- datasplit(trainy, k.fold = cv.fold)}

  names(longlat) <- c("long", "lat")

  # cross validation
  n <- nrow(trainx)
  p <- ncol(trainx)
  cv.pred <- NULL

  if (validation == "LOO") {
    for (i in 1 : length(trainy)) {
      data.dev <- trainx[idx != i, , drop = FALSE]
      data.pred <- trainx[idx == i, , drop = FALSE]

      # gbm modelling
      gbm1 <- gbm::gbm(trainy[idx != i] ~ ., data = data.dev,
                           var.monotone = var.monotone,
                           distribution = as.character(family),
                           n.trees = n.trees,
                           shrinkage = learning.rate,
                           interaction.depth = interaction.depth,
                           bag.fraction = bag.fraction,
                           train.fraction = train.fraction,
                           n.minobsinnode = n.minobsinnode,
                           weights = weights[idx != i],
                           cv.folds = cv.fold,
                           keep.data = keep.data,
                           verbose = verbose,
                           n.cores = n.cores)
      # gbm predictions
      best.iter <- gbm::gbm.perf(gbm1, method = "cv")
      print(best.iter)
      pred.gbm1 <- gbm::predict.gbm(gbm1, data.pred, n.trees = best.iter, type = "response")

      # the residuals of gbm for krige
      dev.gbm1 <- gbm::predict.gbm(gbm1, data.dev, n.trees = best.iter, type = "response")

      data.dev1 <- longlat[idx != i, , drop = FALSE] # for krige
      data.pred1 <- longlat[idx == i, , drop = FALSE] # for krige

      res1 <- trainy[idx != i] - dev.gbm1

      data.dev1$res1 <- res1

      # idw of the residuals
      gstat1 <- gstat::gstat(id = "res1", formula = res1 ~ 1, locations = ~ long + lat, data = data.dev1, set = list(idp = idp), nmax = nmaxidw)

      # idw predictions
      pred.idw1 <- stats::predict(gstat1, data.pred1)

      # for krige
      if (transformation == "none") {data.dev1$res1 = res1} else (
        if (transformation == "sqrt") {data.dev1$res1 = sqrt(res1 + abs(min(res1)))} else (
          if (transformation == "arcsine") {data.dev1$res1 = asin(sqrt((res1 + abs(min(res1))) / 100))} else (
            if (transformation == "log") {data.dev1$res1 = log(res1 + abs(min(res1)) + delta)} else (
              stop ("This transfromation is not supported in this version!")))))
        # The '+ abs(min(res1))' above is to set possible negative values to 0.

      # vgm of the residuals
      sp::coordinates(data.dev1) = ~ long + lat
      vgm1 <- gstat::variogram(object = formula, data.dev1, alpha = alpha)
      model.1 <- gstat::fit.variogram(vgm1, gstat::vgm(mean(vgm1$gamma), vgm.args, mean(vgm1$dist), min(vgm1$gamma)/10, anis = anis))
      if (model.1$range[2] <= 0) (cat("A zero or negative range was fitted to variogram", "\n"))
      if (model.1$range[2] <= 0) (model.1$range[2] <- min(vgm1$dist))  # set negative range to be positive

      # krige predictions
      sp::coordinates(data.pred1) = ~long + lat
      pred.krige1 <- gstat::krige(formula = formula, data.dev1, data.pred1, model = model.1, nmax=nmaxkrige, block = block, beta = beta)$var1.pred

      if (transformation == "none") {pred.krige = pred.krige1}
      if (transformation == "sqrt") {pred.krige = pred.krige1 ^ 2 - abs(min(res1))}
      if (transformation == "arcsine") {pred.krige = (sin(pred.krige1)) ^ 2 * 100 -  abs(min(res1))}
      if (transformation == "log") {pred.krige = exp(pred.krige1) - abs(min(res1)) - delta}

      cv.pred[idx == i] <- (pred.krige * (2 - lambda) + pred.idw1$res1.pred * lambda + pred.gbm1 * hybrid.parameter) / hybrid.parameter
    }
  }

  if (validation == "CV") {
  for (i in 1 : cv.fold) {
    data.dev <- trainx[idx != i, , drop = FALSE]
    data.pred <- trainx[idx == i, , drop = FALSE]

    # gbm modelling
    gbm1 <- gbm::gbm(trainy[idx != i] ~ ., data = data.dev,
                     var.monotone = var.monotone,
                     distribution = as.character(family),
                     n.trees = n.trees,
                     shrinkage = learning.rate,
                     interaction.depth = interaction.depth,
                     bag.fraction = bag.fraction,
                     train.fraction = train.fraction,
                     n.minobsinnode = n.minobsinnode,
                     weights = weights[idx != i],
                     cv.folds = cv.fold,
                     keep.data = keep.data,
                     verbose = verbose,
                     n.cores = n.cores)
    # gbm predictions
    best.iter <- gbm::gbm.perf(gbm1, method = "cv")
    print(best.iter)
    pred.gbm1 <- gbm::predict.gbm(gbm1, data.pred, n.trees = best.iter, type = "response")

    # the residuals of gbm for krige
    dev.gbm1 <- gbm::predict.gbm(gbm1, data.dev, n.trees = best.iter, type = "response")

    data.dev1 <- longlat[idx != i, , drop = FALSE] # for krige
    data.pred1 <- longlat[idx == i, , drop = FALSE] # for krige
    res1 <- trainy[idx != i] - dev.gbm1

    data.dev1$res1 <- res1

    # idw of the residuals
    gstat1 <- gstat::gstat(id = "res1", formula = res1 ~ 1, locations = ~ long + lat, data = data.dev1, set = list(idp = idp), nmax = nmaxidw)

    # idw predictions
    pred.idw1<- stats::predict(gstat1, data.pred1)

    if (transformation == "none") {data.dev1$res1 = res1} else (
      if (transformation == "sqrt") {data.dev1$res1 = sqrt(res1 + abs(min(res1)))} else (
        if (transformation == "arcsine") {data.dev1$res1 = asin(sqrt((res1 + abs(min(res1))) / 100))} else (
          if (transformation == "log") {data.dev1$res1 = log(res1 + abs(min(res1)) + delta)} else (
            stop ("This transfromation is not supported in this version!")))))
    # The '+ abs(min(res1))' above is to set possible negative values to 0.

    # vgm of the residuals
    sp::coordinates(data.dev1) = ~ long + lat
    vgm1 <- gstat::variogram(object = formula, data.dev1, alpha = alpha)
    model.1 <- gstat::fit.variogram(vgm1, gstat::vgm(mean(vgm1$gamma), vgm.args, mean(vgm1$dist), min(vgm1$gamma)/10, anis = anis))
    if (model.1$range[2] <= 0) (cat("A zero or negative range was fitted to variogram", "\n"))
    if (model.1$range[2] <= 0) (model.1$range[2] <- min(vgm1$dist))  # set negative range to be positive

    # krige predictions
    sp::coordinates(data.pred1) = ~long + lat
    pred.krige1 <- gstat::krige(formula = formula, data.dev1, data.pred1, model = model.1, nmax=nmaxkrige, block = block, beta = beta)$var1.pred

    if (transformation == "none") {pred.krige = pred.krige1}
    if (transformation == "sqrt") {pred.krige = pred.krige1 ^ 2 - abs(min(res1))}
    if (transformation == "arcsine") {pred.krige = (sin(pred.krige1)) ^ 2 * 100 -  abs(min(res1))}
    if (transformation == "log") {pred.krige = exp(pred.krige1) - abs(min(res1)) - delta}

    cv.pred[idx == i] <- (pred.krige * (2 - lambda) + pred.idw1$res1.pred * lambda + pred.gbm1* hybrid.parameter) / hybrid.parameter
    }
  }

  # predicitve error and accuracy assessment
  if (predacc == "VEcv") {predictive.accuracy = spm::vecv(trainy, cv.pred)} else (
  if (predacc == "ALL") {predictive.accuracy = spm::pred.acc(trainy, cv.pred)} else (
  stop ("This measure is not supported in this version!")))
  predictive.accuracy
}

