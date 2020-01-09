#' Silly Printer function
#'
#' @param r What you want in the second column.
#' @param x What you want in the first column.
#'
#' @return A tibble
#' @export
#'
#' @importFrom tibble data_frame
#' @importFrom  utils head adist
#' @examples
#' sillyPrinter(x = rnorm(5), r = rnorm(5))
sillyPrinter = function(r,x){
  y = data_frame(x = x, r = r)
  print(head(y))
  return(y)
}
