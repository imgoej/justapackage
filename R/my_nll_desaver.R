# my_nll_desaver ####

#' my_nll_desaver
#'
#' @export
my_nll_desaver <- function(x, bin.file) {
  # create matrix with rows of 6 WITHOUT duplicated values if length is not a match
  X <- matrix(x[1:(6 * ceiling(length(x)/6))], nrow = 6)

  # if length of x is not a multiple of 6, make a blank cell
  X[is.na(X)] <- ""

  # convert matrix to a data.frame
  X <- as.data.frame(X)

  # set binx-file name as column name
  colnames(X) <- rep(bin.file, ncol(X))

  # copy to clipboard
  write.table(X, "clipboard", sep = "\t", col.names = T, row.names = F)
}
