# my_calc_sem ####

#' my_calc_sem
#' @import ggplot2
#' @import Luminescence
#' @export
my_calc_sem <- function(object) {
  if(sum(is.na(object)) > 0) warning("Careful of those NA's")
  sd(object) / sqrt(length(object))
}
