# my_criteria_function2 ####

#' my_criteria_function2
#' @import ggplot2
#' @import Luminescence
#' @export
my_criteria_function2 <- function(data = NULL, index, blanks = NULL, keep.alternatives = F) {
  if(!is.null(blanks)) {
    index.blanks <- sapply(1:nrow(blanks), function(x) {
      !(data$Disc == blanks[x, 1] & data$Grain == blanks[x, 2])
    })
    index.blanks <- (rowSums(index.blanks)-ncol(index.blanks)+1) == 1
    index <- cbind.data.frame(index, index.blanks)
  }
  ncol.index <- ncol(index)
  index$status <- rowSums(index) - ncol(index)
  index$status <- index$status == 0
  index$status[is.na(index$status)] <- FALSE
  index$status2[index$status == TRUE] <- "OK"
  index$status2[index$status == FALSE] <- "FAILED"
  index$status2 <- factor(index$status2, levels = c("OK", "FAILED"))
  prop <- data.frame("total" = length(index$status), "accepted" = sum(index$status))
  prop$rejected <- prop$total-prop$accepted
  prop$p <- prop$accepted / prop$total

  trial <- function(combo.col, keep.alt = F) {
    N <- as.numeric(index.combinations.unrepeated.matrix[combo.col, ])
    crit.data <- index[, as.numeric(N)]
    X <- rowSums(crit.data) - ncol(crit.data)
    X <- X == 0
    X[is.na(X)] <- FALSE
    if(isTRUE(keep.alt)){
      return(X)
    } else{
      return(sum(X, na.rm = T))
    }
  }

  index.combinations <- Map(combn, list(1:ncol.index), seq_along(1:ncol.index), simplify = FALSE)
  index.combinations.unrepeated <- unlist(index.combinations, recursive = FALSE)
  index.combinations.unrepeated.matrix <- suppressWarnings(do.call(rbind, index.combinations.unrepeated))

  index.comb.acceptance = sapply(1:nrow(index.combinations.unrepeated.matrix), trial)
  index.comb.names = index.combinations.unrepeated
  index.comb.prop = index.comb.acceptance / nrow(data)
  index.comb = cbind.data.frame(
    index.comb.acceptance,
    index.comb.prop)
  index.comb$crit.comb = index.comb.names
  index.comb = index.comb[, c(3,1,2)]
  index.comb$combination.index = 1:nrow(index.combinations.unrepeated.matrix)
  colnames(index.comb) = c("combination","n.accepted", "n.prop","combination.index")

  if(!is.null(data)) {
    p <- ggplot(data, aes(x = as.factor(Grain), y = as.factor(Disc))) +
      geom_tile(aes(fill = index$status2), alpha = 0.5, color = "white") +
      xlab("Grain") + ylab("Disc") +
      scale_fill_manual(values = c("green", "red3")) +
      theme_classic() + theme(legend.position = "none", axis.line = element_blank()) +
      coord_fixed(ratio=1)
    print(p)
  }
  cat("\n----------------------------------------------------------------------------------\n")
  cat(">> aliquot/grain acceptance level << ")
  cat("\n----------------------------------------------------------------------------------\n")
  print(prop, row.names = FALSE)
  cat("\n\n")
  cat("\n----------------------------------------------------------------------------------\n")
  cat(">> criteria impact overview << ")
  cat("\n----------------------------------------------------------------------------------\n")
  print(index.comb, row.names = FALSE, na.print = "", digits = 3)
  cat("\n\n")

  if(keep.alternatives == T) {
    index.alternatives = lapply(1:nrow(index.combinations.unrepeated.matrix), function(x) {trial(x, keep.alt = T)})
    return(index.alternatives)
  } else {
    return(index$status)
  }
}
