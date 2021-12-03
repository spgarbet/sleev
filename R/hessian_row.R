#' Hessian Row
#' 
#' Multiply two variables and return colsums
#' TODO
#' 
#' @param x 
#' @param pm
hessian_row <- function(x, pm) {
  return(colSums(x * pm))
}
