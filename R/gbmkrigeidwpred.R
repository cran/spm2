#' @title Generate spatial predictions using the hybrid methods of
#' generalized boosted regression modeling ('gbm'), 'kriging' and inverse distance weighted ('IDW').
#'
#' @description This function is for generating spatial predictions using the
#' hybrid methods of 'gbm', 'kriging' and 'IDW', including all methods implemented
#' in 'gbmkrigeidwcv'.
#'
#' @param longlat	a dataframe contains longitude and latitude of point samples.
#' @param trainx a dataframe contains longitude (long), latitude (lat),
#' predictive variables and the response variable of point samples.
#' @param predx	a dataframe or matrix contains columns of predictive variables
#' for the grids to be predicted.
#' @param trainy a vector of the response variable in the formula, that is, the left
#' part of the formula.
#' @param longlatpredx	a dataframe contains longitude and latitude of point locations
#' (i.e., the centers of grids) to be predicted.
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
#' @param train.fraction The first 'train.fraction * nrows(data)' observations
#' are used to fit the gbm and the remainder are used for computing
#' out-of-sample estimates of the loss function.
#' @param n.minobsinnode minimum number of observations in the trees terminal
#' nodes. Note that this is the actual number of observations not the total
#' weight. By default, 10 is used.
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
#' @param transformation transform the residuals of 'gbm' to normalise the data;
#' can be "sqrt" for square root, "arcsine" for arcsine, "log" or "none"
#' for non transformation. By default, "none" is used.
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
#' @param cv.fold integer; number of folds in the cross-validation. if > 1,
#' then apply n-fold cross validation; the default is 10, i.e., 10-fold cross
#' validation that is recommended.
#' @param n.cores The number of CPU cores to use. See gbm for details. By
#' default, 6 is used.
#' @param ... other arguments passed on to 'gbm', 'krige' and 'gstat'.
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
#' Greg Ridgeway with contributions from others (2015). gbm: Generalized
#' Boosted Regression Models. R package version 2.1.1.
#' https://CRAN.R-project.org/package=gbm
#'
#' Pebesma, E.J., 2004. Multivariable geostatistics in S: the gstat package.
#' Computers & Geosciences, 30: 683-691.
#'
#' @author Jin Li
#' @examples
#' \donttest{
#' library(spm)
#'
#' data(sponge)
#' data(sponge.grid)
#' longlat <- sponge[, 1:2]
#'
#' set.seed(1234)
#'
#' gbmkrigeidwpred1 <- gbmkrigeidwpred(longlat = longlat, trainx = sponge[, -3],
#' predx = sponge.grid, trainy = sponge[, 3], longlatpredx = sponge.grid[, c(1:2)],
#' family = "poisson", interaction.depth = 3, transformation = "none", formula = res1 ~ 1,
#' vgm.args = "Sph", nmaxkrige = 12, idp = 2, nmaxidw = 12, hybrid.parameter = 3,
#' n.cores = 8)
#'
#' names(gbmkrigeidwpred1)
#'
#' range(gbmkrigeidwpred1$predictions)
#'}
#'
#' @export
gbmkrigeidwpred <- function (longlat, trainx, predx, trainy, longlatpredx,
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
                             cv.fold = 10,
                             n.cores = 8,
                             ...) {

  n <- nrow(trainx)
  p <- ncol(trainx)
  names(longlat) <- c("long", "lat")
  names(longlatpredx) <- c("long", "lat")

  # gbm modeling
  gbm1 <- gbm::gbm(trainy ~ ., data = trainx,
                   var.monotone = var.monotone,
                   distribution = as.character(family),
                   n.trees = n.trees,
                   shrinkage = learning.rate,
                   interaction.depth = interaction.depth,
                   bag.fraction = bag.fraction,
                   train.fraction = train.fraction,
                   n.minobsinnode = n.minobsinnode,
                   weights = weights,
                   cv.folds = cv.fold,
                   keep.data = keep.data,
                   verbose = verbose,
                   n.cores = n.cores)

  # gbm predictions
  best.iter <- gbm::gbm.perf(gbm1, method = "cv")
  print(best.iter)
  pred.gbm1 <- gbm::predict.gbm(gbm1, predx, n.trees = best.iter, type = "response")

  # the residuals of gbm for krige
  data.dev1 <- longlat
  data.pred1 <- longlatpredx
  dev.gbm1 <- gbm::predict.gbm(gbm1, trainx, n.trees = best.iter, type = "response")

  res1 <- trainy - dev.gbm1
  data.dev1$res1 <- res1

  # idw of the residuals
  gstat1 <- gstat::gstat(id = "res1", formula = res1 ~ 1, locations = ~ long + lat, data = data.dev1, set = list(idp = idp), nmax = nmaxidw)

  # idw predictions
  pred.idw1 <- stats::predict(gstat1, data.pred1)

  res1 <- trainy - dev.gbm1
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

  predictions <- (pred.krige * (2 - lambda) + pred.idw1$res1.pred * lambda + pred.gbm1 * hybrid.parameter) / hybrid.parameter
  gbmkrigeidw.pred <- cbind(longlatpredx, predictions)
  gbmkrigeidw.pred
}

