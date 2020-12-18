#' Cam per D0
#'
#' Calculates equivalent dose for specified D0 thresholds.
#' Estimates are then plotted against chosen D0 thresholds and the point of intersect is estimated using linear regression between the two points containing the intersect.
#'
#' @import ggplot2
#' @import Luminescence
#'@param data 3-column data frame containing (in order) De, De.error and D0 features. Alternatively you can specify colnames or column numbers of said features (see De/De.error/D0).
#'@param minD0 Vector of D0 thresholds.
#'@param method Method of estimating central tendency in equivalent dose analysis. Currently supported are "CAM" (default) and "ADM".
#'@param De/De.error/D0 Column names of relevant variables in the OSL dataset. Selects the first three columns by default.
#'@param ... Arguments passed to CAM or ADM. For example one might be interested in specifying "log" or "sigmab" in the CAM procedure or "sigma_m" in the ADM procedure.
#'@return Returns a dataframe of equivalent doses and D0 thresholds along with intersect coordinate(s) and error estimates.
#'@export
cam_per_d0 <- function(data, minD0, method = "CAM", De = 1, De.error = 2, D0 = 3, ...) {
  unsatDATA <- data[complete.cases(data[, c(De, De.error, D0)]), ]
  satDATA <- data[!complete.cases(data[, c(De, De.error)]), ]
  satDATA <- satDATA[!is.na(satDATA[, D0]), ]
  minD0 <- minD0[minD0 <= max(unsatDATA[, D0])]
  DATA <- lapply(minD0, function(x) {
    unsatDATA.tmp <- unsatDATA[unsatDATA[, D0] >= x,]
    n <- nrow(unsatDATA.tmp)
    nsaturated <- nrow(satDATA[satDATA[, D0] >= x, ])
    if(nrow(unsatDATA.tmp) <= 1) DATA.tmp <- NA
    if(x > max(unsatDATA[, D0])) DATA.tmp <- NA
    else {
      if(method == "CAM")
        DATA.tmp <- get_RLum(calc_CentralDose(
          unsatDATA.tmp[, c(De, De.error)],
          plot = FALSE,
          verbose = FALSE,
          ...))[1:4]
      if(method == "ADM")
        DATA.tmp <- get_RLum(calc_AverageDose(
          unsatDATA.tmp[, c(De, De.error)],
          plot = F,
          verbose = F,
          sigma_m = 0.2,
          ...))[1:2]
      DATA.tmp <- as.data.frame(DATA.tmp)
    }
    DATA.tmp <- cbind.data.frame(DATA.tmp, "n" = n, "nsaturated" = nsaturated, "nsaturatedPROP" = nsaturated/(nsaturated+n))
    return(DATA.tmp)
  })
  DATA <- do.call(rbind, DATA)
  DATA <- cbind.data.frame("minD0" = minD0, DATA)
  DATA <- cbind.data.frame(DATA, "method" = method)
  colnames(DATA)[2:3] <- c("de", "de.error")
  above <- DATA$de > DATA$minD0
  if(sum(abs(diff(above))) == 0) {
    intersect.results <- data.frame("intersect" = 0, "x-coordinates" = 0, "x-error" = 0, "y-coordinates" = 0, "y-error" = 0)
    intersect.x = 0
    intersect.y = 0
    intersect.points.start = NULL
    intersect.points.end = NULL
    segment_data <- data.frame(
      "x" = 0,
      "xend" = 2,
      "y" = 2,
      "yend" = 2
    )
    ggplotlist <- list()
  }
  else {
    intersect.points.start <- which(diff(above) != 0)
    intersect.points.end <- intersect.points.start+1
    slopes <- (DATA$de[intersect.points.end] - DATA$de[intersect.points.start]) / (DATA$minD0[intersect.points.end] - DATA$minD0[intersect.points.start])
    yintercepts <- DATA$de[intersect.points.end]-slopes*DATA$minD0[intersect.points.end]
    intersect.x = (yintercepts-0) / (1-slopes)
    intersect.y = yintercepts + intersect.x*slopes
    segment_data <- data.frame(
      "x" = 0,
      "xend" = intersect.x,
      "y" = intersect.x,
      "yend" = intersect.x
    )
    x1 = DATA$minD0[intersect.points.start]
    x2 = DATA$minD0[intersect.points.end]
    y1 = DATA$de[intersect.points.start]
    y2 = DATA$de[intersect.points.end]

    dx = sqrt( ((y1-x1)*(x2-x1)/(x2-x1-y2+y1)^2*DATA$de.error[intersect.points.end])^2 + ((x2^2-x2*x1-x2*y2+x1*y2)/(x2-x1-y2+y1)^2*DATA$de.error[intersect.points.start])^2)


    intersect.results <- data.frame("intersect" = 1:length(intersect.x), "x.y.coordinate" = intersect.x, "error" = dx)
    ggplotlist <- list(geom_segment(data = segment_data, aes(x = x, xend = xend, y = y, yend = yend), linetype = "dashed", color = "red"),
                       geom_segment(data = segment_data, aes(x = xend, xend = xend, y = x, yend = yend), linetype = "dashed", color  = "red"),
                       geom_segment(data = segment_data, aes(x = DATA$minD0[intersect.points.start], xend = DATA$minD0[intersect.points.end], y = DATA$de[intersect.points.start], yend = DATA$de[intersect.points.end]), linetype = "dashed"))
  }
  legenddata <- data.frame("label" = c("CAM", "1to1 line", "Regression line", "intersect line"), "dummyx" = 0, "dummyy" = 0)

  suppressMessages(
    p <- ggplot(DATA, aes(x = minD0, y = de)) + geom_point(size = 3, shape = 23, fill = "black") + geom_errorbar(aes(ymin = de-de.error, ymax = de+de.error, width = 0.1*max(de))) + xlim(0, max(DATA$de)*2) + ylim(0, max(minD0)*2) + scale_x_continuous(breaks = minD0) +
      ggplotlist +
      geom_abline(size = 0.2) + theme_classic() +
      ylab(bquote("D"["e"])) + xlab("D0 threshold") +
      geom_line(data = legenddata, aes(x = dummyx, y = dummyy, color = label, linetype = label)) + scale_linetype_manual(values = c("solid", "blank","dashed", "dashed")) + scale_color_manual(values = c("black", "black", "red", "black")) +
      geom_point(data = legenddata, aes(x = DATA$minD0[1], y = DATA$de[1], alpha = label), shape = 23, size = 3, fill = "black") + scale_alpha_manual(values = c(0,1,0,0)) +
      theme(legend.position = c(1, 1), legend.justification = c(1,1), legend.title = element_blank())
  )
  cat("\n----------------------------------------------------------------------------------\n")
  cat(">> Results << ")
  cat("\n----------------------------------------------------------------------------------\n")
  print.data.frame(DATA, print.gap = 2, digits = 4, row.names = F)
  cat("\n----------------------------------------------------------------------------------\n")
  cat(">> intersects << ")
  cat("\n----------------------------------------------------------------------------------\n")
  print.data.frame(intersect.results, row.names = F)
  if(sum(abs(diff(above))) == 0) print("I could not find any intersects")
  suppressMessages (print(p))
  return(list(DATA, intersect.results))
}
