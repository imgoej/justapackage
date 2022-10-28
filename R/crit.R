# crit ####

#' crit
#'
#' @export
crit <- function(criteria, error, signif = 4, limits, limittype = "<>") {
  criteria  <- round(criteria, signif)
  if(is.null(error)) {error = 0}
  error     <- round(error, signif)
  data      <- do.call(data.frame, lapply(data.frame(criteria, error), function(x) replace(x, is.infinite(x), NA)))
  criteria <- data[, 1]
  error <- data[, 2]
  if(limittype == "<>") {
    index <- (criteria + abs(error) >= limits[1]) & (criteria - abs(error) <= limits[2])
  }
  if(limittype == ">") {
    index <- (criteria - abs(error) > limits[1])
  }
  if(limittype == "<") {
    index <- (criteria - abs(error) < limits[1])
  }
  index
}
