# my_blanks_checker ####

#' my complete cases
#'
#' @export
my_blanks_checker <- function(x, signal.threshold, interval = 6:11) {
  my.data.frame <- data.frame("Signal" = sapply(1:length(x@DATA), function(y) { sum(x@DATA[[y]][interval]) }))
  my.data.frame <- cbind("Disc" = x@METADATA$POSITION, "Grain" = x@METADATA$GRAIN, my.data.frame)
  my.data.frame <- my.data.frame[my.data.frame$"Signal" > signal.threshold,]
  return(my.data.frame)
}
