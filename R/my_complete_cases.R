# my_complete_cases ####

#' my complete cases
#'
#' @export
my_complete_cases <- function(data) {
  data <- do.call(data.frame, lapply(data, function(x) replace(x, is.infinite(x), NA)))
  complete.cases(data)
}
