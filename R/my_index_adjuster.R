# my_index_adjuster ####

#' my_index_adjuster
#' @export
my_index_adjuster <- function(index = NULL, add = NULL, remove = NULL) {
  if(!is.null(add)){
    for (i in 1:length(index)) {
      index[[i]][add] <- TRUE
    }
  }
  if(!is.null(remove)){
    for (i in 1:length(index)) {
      index[[i]][remove] <- FALSE
    }
  }
  index
}
