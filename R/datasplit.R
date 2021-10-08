#' @title Split data for k-fold cross-validation
#'
#' @description This function is a data splitting function for k-fold cross-
#' validation and uses a stratified random sampling technique. It resamples the
#' training data based on sample quantiles.
#'
#' @param trainy a vector of response, must have a length equal to sample size.
#' @param k.fold integer; number of folds in the cross-validation. if > 1,
#' then apply k-fold cross validation; the default is 10, i.e., 10-fold cross
#' validation that is recommended.
#'
#' @return A list of samples each with an index of k-fold number.
#'
#' @note This function is largely based on rfcv in randomForest.
#'
#' @references A. Liaw and M. Wiener (2002). Classification and Regression
#' by randomForest. R News 2(3), 18-22.
#'
#' @author Jin Li
#' @examples
#' library(spm)
#' data(petrel)
#' idx1 <- datasplit(petrel[, 3], k.fold = 10)
#' table(idx1)
#'
#' @export
datasplit <- function (trainy, k.fold = 10) {
  classRF <- is.factor(trainy)
  n <- length(trainy)
  if (classRF) {
    stop ("This function is not for categorical response variable")
  } else {
    f <- cut(trainy, c(-Inf, stats::quantile(trainy, 1:4/5), Inf))
  }
  nlvl <- table(f)
  idx <- numeric(n)
  for (i in 1:length(nlvl)) {
    idx[which(f == levels(f)[i])] <- sample(rep(1:k.fold,
    length = nlvl[i]))
  }
  idx
}
