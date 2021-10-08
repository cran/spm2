#' @title Digit number after decimal point for a numeric variable
#'
#' @description This function is to derive the digit number after decimal point for a
#' numeric variable (e.g., lat and long).
#'
#' @param x one or more decimal numbers.
#' @param dechar The character used to separate the decimal part of a number.
#' @param nint The number of characters to which the integer part of the numbers should be padded.
#' @param ndec The number of characters to which the decimal part of the numbers should be padded.
#' @param pad.left Whether the left (integer) side of the numbers should be padded as well as the right.
#'
#' @return A list of integer number to show digit number after decimal point of x.
#'
#' @note This function is modified from decimal.align in 'prettyR' package.
#'
#' @references Jim Lemon and Philippe Grosjean (2019). 'prettyR': Pretty Descriptive Stats. R package version 2.1.1.
#' https://CRAN.R-project.org/package=prettyR.
#'
#' @author Jin Li
#' @examples
#'
#' x<-c(0.1, 2.2, 3.03, 44.444, 555.0005, 6666.66666)
#' decimaldigit(x)
#'
#' @export
decimaldigit <- function (x, dechar = ".", nint = NA, ndec = NA, pad.left = TRUE)
{
  splitchar <- paste("[", dechar, "]", sep = "")
  splitlist <- strsplit(as.character(x), splitchar)

  ndec <- max(nchar(unlist(lapply(splitlist, "[", 2))))
  decs <- unlist(lapply(splitlist, "[", 2))
  decs[is.na(decs)] <- "0"
  digit.no <- nchar(decs)
  return(digit.no)
}
