#' @title Correct classification rate for predictive models based on cross
#' -validation
#'
#' @description This function is to calculates correct classification
#' rate (ccr) for categorical data with the observed (obs) data specified
#' as factor. It based on the differences between the predicted values for
#' and the observed values of validation samples for cross-validation. For 0
#' and 1 data, the observed values need to be specified as factor in order
#' to use this accuracy measure. It is modified from the function 'pred.acc'
#'  in 'spm' package.
#'
#' @param obs a vector of observation values of validation samples.
#' @param pred a vector of prediction values of predictive models for validation samples.
#' @return A list with the following component:
#' ccr (correct classification rate) for categorical data.
#'
#' @references Jin Li (2019). spm: Spatial Predictive Modeling. R package
#' version 1.2.0. https://CRAN.R-project.org/package=spm.
#'
#' @author Jin Li
#' @examples
#' set.seed(1234)
#' x <- as.factor(sample(letters[1:2], 30, TRUE))
#' y <- sample(x, 30)
#' ccr(x, y)
#'
#' @export
ccr <- function (obs, pred) {
  if (is.factor(obs)) {
      data1 <- as.data.frame(cbind(pred, obs))
      sum(diag(table(data1))) / sum(table(data1)) * 100
    } else (
      stop ("This is for categorical data!"))
}
