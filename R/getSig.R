#' Title Obtaining significance through p-values.
#'
#' @param p p value
#'
#' @return Significance shown in "***","**","*"
#' @export
#' @examples
#' getSig(0.04)
getSig <- function(p) {
  if (p < 0.001) {"***"}
  else if(p < 0.01) {"**"}
  else if (p < 0.05) {"*"}
  else {""}
}
