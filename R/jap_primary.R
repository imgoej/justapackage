# my_complete_cases ####

#' my complete cases
#'
#' @export
my_complete_cases <- function(data) {
  data <- do.call(data.frame, lapply(data, function(x) replace(x, is.infinite(x), NA)))
  complete.cases(data)
}

# cam per D0 ####

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
cam_per_d0 <- function(data, minD0, method = "CAM", error.mode = "propagate", ...) {
  argg <- c(as.list(environment()), list(...))
  argg <- do.call(cbind, argg[lengths(argg) == 1])
  colnames(data) <- c("De", "De.error", "D0")
  unsatDATA <- data[my_complete_cases(data[, c("De", "De.error", "D0")]), ]
  satDATA <- data[!my_complete_cases(data[, c("De", "De.error")]), ]
  satDATA <- satDATA[!is.na(satDATA[, "D0"]), ]
  minD0 <- minD0[minD0 <= max(unsatDATA[, "D0"])]
  DATA.print <- lapply(minD0, function(x) {
    unsatDATA.tmp <- unsatDATA[unsatDATA[, "D0"] >= x,]
    n <- nrow(unsatDATA.tmp)
    nsaturated <- nrow(satDATA[satDATA[, "D0"] >= x, ])
    if(nrow(unsatDATA.tmp) > 1 & x < max(unsatDATA[, "D0"])) {
      if(method == "CAM")
        DATA.tmp <- get_RLum(calc_CentralDose(
          unsatDATA.tmp[, c("De", "De.error")],
          plot = FALSE,
          verbose = FALSE,
          ...))[1:4]
      if(method == "ADM")
        DATA.tmp <- get_RLum(calc_AverageDose(
          unsatDATA.tmp[, c("De", "De.error")],
          plot = F,
          verbose = F,
          ...))[1:2]
      DATA.tmp <- as.data.frame(DATA.tmp)
      DATA.tmp <- cbind.data.frame(DATA.tmp, "n" = n, "nsaturated" = nsaturated, "nsaturatedPROP" = nsaturated/(nsaturated+n), argg)
    } else {DATA.tmp <- NA
    warning("Some of the chosen D0 thresholds contain less than two observations and cannot be evaluated")}
   return(DATA.tmp)
  })
  DATA.print <- do.call(rbind, DATA.print)
  DATA.print <- cbind.data.frame("minD0" = minD0, DATA.print)
  colnames(DATA.print)[2:3] <- c("de", "de.error")
  DATA <- DATA.print[!is.na(DATA.print[, "de"]), ]
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
  } else {
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

    if(error.mode == "propagate") {
      dy1 = -(( (x2-x1)*(y2-x2) ) / (y1-y2+x2-x1)^2)
      dy2 = ( (x2-x1)*(y1-x1) ) / (y2-y1-x2+x1)^2
      dx = abs(dy1*DATA$de.error[intersect.points.start]) + abs(dy2*DATA$de.error[intersect.points.end])
    }
    if(error.mode == "w.error"){
      distance = abs(x2-x1)
      dratios = (intersect.x-x1)/distance
      dx = DATA$de.error[intersect.points.start]*(1-dratios)+DATA$de.error[intersect.points.end]*dratios
      }

    intersect.results = data.frame("intersect" = 1:nrow(DATA.print), "x.y.coordinate" = NA, "error" = NA)
    intersect.results[1:length(intersect.x), c("x.y.coordinate", "error")] <- cbind(intersect.x, dx)
    ggplotlist <- list(geom_segment(data = segment_data, aes(x = x, xend = xend, y = y, yend = yend), linetype = "dashed", color = "red"),
                       geom_segment(data = segment_data, aes(x = xend, xend = xend, y = x, yend = yend), linetype = "dashed", color  = "red"),
                       geom_segment(data = segment_data, aes(x = DATA$minD0[intersect.points.start], xend = DATA$minD0[intersect.points.end], y = DATA$de[intersect.points.start], yend = DATA$de[intersect.points.end]), linetype = "dashed"))
  }
  legenddata <- data.frame("label" = c(method, "1to1 line", "Regression line", "intersect line"), "dummyx" = 0, "dummyy" = 0)

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
  print.data.frame(DATA.print, print.gap = 2, digits = 4, row.names = F)
  cat("\n----------------------------------------------------------------------------------\n")
  cat(">> intersects << ")
  cat("\n----------------------------------------------------------------------------------\n")
  print.data.frame(intersect.results, row.names = F)
  if(sum(abs(diff(above))) == 0) print("I could not find any intersects")
  suppressMessages (print(p))
  return(list(DATA.print, intersect.results))
}

# my SAR analyzer ####

#' My SAR analyzer
#'
#' Does things.
#' And returns things.
#'
#' @import Luminescence
#' @export
my_SAR_analyzer2 <- function(Risoe.object, position = NULL, grain = NULL, run = NULL, set = NULL, ltype = "OSL", dtype = NULL, protocol = "unknown", keep.empty = TRUE, txtProgressBar = FALSE,
                             signal.integral = c(6:11), background.integral = c(350:400), rejection.criteria = list(
                               recycling.ratio = NA,
                               recuperation.rate = NA,
                               testdose.error = NA,
                               palaeodose.error = NA,
                               exceed.max.regpoint = NA),
                             sig0 = 0.025,
                             sigmab = NULL,
                             background.count.distribution =  "poisson",
                             fit.method = "EXP",
                             fit.force_through_origin = TRUE,
                             fit.weights = TRUE,
                             fit.includingRepeatedRegPoints = TRUE,
                             plot = FALSE,
                             plot.single = FALSE,
                             onlyLxTxTable = FALSE,
                             readerrate = 1,
                             drdose = 1,
                             label = "sample",
                             recuprecyc = TRUE
) {
  if(is.null(position)) { position <- unique(Risoe.object@METADATA$POSITION)}
  RLum.object <- Risoe.BINfileData2RLum.Analysis(Risoe.object,
                                                 pos = position,
                                                 grain = grain,
                                                 run = run,
                                                 set = set,
                                                 ltype = ltype,
                                                 dtype = dtype,
                                                 protocol = protocol,
                                                 keep.empty = keep.empty,
                                                 txtProgressBar = txtProgressBar)

  if(!is.list(RLum.object)) RLum.object <- list(RLum.object)

  results <- sapply(1:length(RLum.object), function(x) {
    print(c(RLum.object[[x]][[1]]$POSITION, RLum.object[[x]][[1]]$GRAIN))
    analyse_SAR.CWOSL(RLum.object[[x]],
                      signal.integral.min = min(signal.integral),
                      signal.integral.max = max(signal.integral),
                      background.integral.min = min(background.integral),
                      background.integral.max = max(background.integral),
                      rejection.criteria = rejection.criteria,
                      sig0 = sig0,
                      sigmab = sigmab,
                      background.count.distribution = background.count.distribution,
                      fit.method = fit.method,
                      fit.force_through_origin = fit.force_through_origin,
                      fit.weights = fit.weights,
                      fit.includingRepeatedRegPoints = fit.includingRepeatedRegPoints,
                      plot = plot,
                      plot.single = plot.single,
                      onlyLxTxTable = onlyLxTxTable
    )})

  data.tmp <- lapply(1:length(results), function(x) {
    results[[x]]$data
  })
  LnLxTnTx.table.tmp <- lapply(1:length(results), function(x) {
    results[[x]]$LnLxTnTx.table
  })

  formula <- lapply(1:length(results), function(x) {
    results[[x]]$Formula
  })

  results <- list(
    "data" = do.call(rbind.data.frame, data.tmp),
    "LnLxTnTx.table" = do.call(rbind.data.frame, LnLxTnTx.table.tmp),
    "formula" = formula
  )

  Ln <- results$LnLxTnTx.table$LnLx[results$LnLxTnTx.table$Name == "Natural"]
  Ln.BG <- results$LnLxTnTx.table$LnLx.BG[results$LnLxTnTx.table$Name == "Natural"]
  Tn <- results$LnLxTnTx.table$TnTx[results$LnLxTnTx.table$Name == "Natural"]
  Tn.BG <- results$LnLxTnTx.table$TnTx.BG[results$LnLxTnTx.table$Name == "Natural"]


  Disc <- sapply(1:length(RLum.object), function(x) {
    RLum.object[[x]][[1]]$POSITION
  })
  Grain <- sapply(1:length(RLum.object), function(x) {
    RLum.object[[x]][[1]]$GRAIN
  })

  data <- results$data
  PH <- rep(0, nrow(data))
  if(!is.null(position)) {
    PH <- sapply(position, function(x){
      Risoe.object@METADATA$AN_TEMP[Risoe.object@METADATA$POSITION == x & Risoe.object@METADATA$LTYPE == ltype][1]
    }) }
  if(PH[1] == 0) {
    PH <- sapply(position, function(x){
      Risoe.object@METADATA$TEMPERATURE[Risoe.object@METADATA$POSITION == x & Risoe.object@METADATA$LTYPE == ltype][1]
    })
  }

  if(recuprecyc == TRUE) {
    Recyc_Depl_Recup <- lapply(results$data$UID, function(x){
      B <- results$LnLxTnTx.table[results$LnLxTnTx.table$UID == x, ]
      repeated <- B[, c("Name", "Dose")] [B$Repeated == TRUE,]
      origin <- B[, c("Name", "Dose")] [B$Repeated == FALSE & B$Dose %in% repeated$Dose,]
      compare <- rbind(origin, repeated)
      split <- split(compare, compare$Dose)

      # recycling & depletion ratio
      list <- lapply(1:length(split), function (y){
        C <- (combn(split[[y]]$Name,2))
        D <- B[match(C[2, ], B$Name), "LxTx"] / B[match(C[1, ], B$Name), "LxTx"]
        names(D) <-paste(C[2, ], C[1, ], sep = "/")
        E <- sqrt((B[match(C[2, ], B$Name), "LxTx.Error"] / B[match(C[2, ], B$Name), "LxTx"])^2 + (B[match(C[1, ], B$Name), "LxTx.Error"] / B[match(C[1, ], B$Name), "LxTx"])^2) * D
        names(E) <-paste(paste(C[2, ], C[1, ], sep = "/"), ".Error", sep = "")
        return(c(D, E))
      })
      Recyc_Depl <- unlist(list)

      # recuperation ratio
      origin <- B[B$Name %in% c("Natural", "R0"), c("Name", "LxTx", "LxTx.Error")]
      Recup <- c(
        "R0.Natural"        = origin$LxTx[2] / origin$LxTx[1],
        "R0.Natural.Error"  = origin$LxTx[2] / origin$LxTx[1] * sqrt( (origin$LxTx.Error[2] / origin$LxTx[2])^2 + (origin$LxTx.Error[1] / origin$LxTx[1])^2)
      )
      return(c(Recyc_Depl, Recup))
    })
  }
  if(recuprecyc == FALSE) {
    Recyc_Depl_Recup <- list(NA)
  }
  if(is.factor(data$De.Error)) {
    data$De.Error <- as.numeric(as.character(data$De.Error))
  }

  n.sig <- length(signal.integral)
  n.BG <- length(background.integral)

  sigma.Tn            <- sqrt(Tn)
  sigma.Tn.BG         <- sqrt(Tn.BG*n.BG/n.sig) / (n.BG / n.sig)
  sigma.Tn.Corrected  <- sqrt(sigma.Tn^2  +  sigma.Tn.BG^2)
  STn20               <- sigma.Tn.Corrected / (Tn-Tn.BG)

  Recyc_Depl_Recup <- do.call(rbind, Recyc_Depl_Recup)
  RDRnames <- sub("/", ".", colnames(Recyc_Depl_Recup))
  colnames(Recyc_Depl_Recup) <- RDRnames
  info <- data.frame(Disc, Grain, PH, Ln, Ln.BG, Tn, Tn.BG)
  DR <- data$De/drdose
  De.data <- data[, c("De", "De.Error", "D01", "D01.ERROR")]*readerrate
  criteria <- STn20
  ID <- data[, c("signal.range", "background.range", "UID")]

  cat("\n----------------------------------------------------------------------------------\n")
  cat(">> All recycle cycles << ")
  cat("\n----------------------------------------------------------------------------------\n")
  print(colnames(Recyc_Depl_Recup))

  return(list(results = results, data = cbind(info, De.data, DR, STn20, Recyc_Depl_Recup, ID, label)))
}

# crit ####

#' Crit
#'
#' Does things.
#' And returns things.
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
    index <- (criteria + abs(error) >= limits[1]) & (criteria - abs(error) < limits[2])
  }
  if(limittype == ">") {
    index <- (criteria - abs(error) > limits[1])
  }
  if(limittype == "<") {
    index <- (criteria - abs(error) < limits[1])
  }
  index
}


# my criteria function ####

#' My criteria function
#'
#' Does things.
#' And returns things.
#'
#' @import ggplot2
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

# my blanks checker ####

#' My Blanks Checker
#'
#' Does things.
#' And returns things.
#'
#' @export
my_blanks_checker <- function(x, signal.threshold, interval = 6:10) {
  my.data.frame <- data.frame("Signal" = sapply(1:length(x@DATA), function(y) { sum(x@DATA[[y]][interval]) }))
  my.data.frame <- cbind("Disc" = x@METADATA$POSITION, "Grain" = x@METADATA$GRAIN, my.data.frame)
  my.data.frame <- my.data.frame[my.data.frame$"Signal" > signal.threshold,]
  return(my.data.frame)
}

# DR plotter ####

#' My DR Plotter
#'
#' Does things.
#' And returns things.
#'
#' @import ggplot2
#' @export
my_plotter_dr <- function(data.object, by = "mean", plot = TRUE, drdose = 1, width = 3){
  if(ncol(data.object) == 3) data.object$label <- "all"
  if(ncol(data.object) == 4) colnames(data.object) <- c("De", "De.Error", "PH", "label")
  data.object <- data.object[complete.cases(data.object[, c("De", "De.Error", "PH")]), ]
  DR <- aggregate(data.object$De, by = list(data.object$PH, data.object$label), FUN = function(x) x/drdose, simplify = F)
  DR.Error <- aggregate(data.object$De.Error, by = list(data.object$PH, data.object$label), FUN = function(x) x/drdose, simplify = F)
  n <- aggregate(data.object$De, by = list(data.object$PH, data.object$label), FUN = function(x) length(x))

  if(by == "mean") {
    data.temp <- lapply(1:nrow(DR), function(j) {
      mean <- mean(DR$x[[j]], na.rm = T)
      sigma <- sd(DR$x[[j]])/sqrt(length(DR$x[[j]]))
      data.frame("PH" = DR$Group.1[j], "DR" = mean, "DR.Error"= sigma, "n" = n$x[[j]], type = by, "label" = DR$Group.2[j])
    })
    data <- do.call(rbind, data.temp)
  }
  if(by == "cam") {
    data.temp <- lapply(1:nrow(DR), function(j) {
      DR.DR.Error <- get_RLum(calc_CentralDose(cbind.data.frame(DR$x[[j]], DR.Error$x[[j]]), verbose = FALSE, plot = FALSE))[1:2]
      data.frame("PH" = DR$Group.1[j], "DR" = DR.DR.Error[1], "DR.Error"= DR.DR.Error[2], "n" = n$x[[j]], type = by, "label" = DR$Group.2[j])
    })
    data <- do.call(rbind, data.temp)
  }
  colnames(data) <- c("PH", "DR", "DR.Error", "n", "type", "label")
  if(isTRUE(plot)){
    p <- ggplot(data, aes(x=PH, y=DR, fill = label)) +
      geom_errorbar(aes(ymin=DR-DR.Error, ymax=DR+DR.Error), width = width) +
      geom_point(size = 3, shape = 23) +
      geom_hline(yintercept = c(0.9, 1, 1.1), linetype = c("dashed", "solid", "dashed")) +
      theme_classic() +
      scale_x_continuous(breaks = data$PH) +
      ylim(0, max(rowSums(data[, c("DR", "DR.Error")]))) +
      ylab("Dose recovery ratio") +
      xlab(bquote("Preheat temperature ("*degree*"C)"))
  }
  if(length(unique(data.object$label)) == 1) {p <- p + theme(legend.position = "none")}
  print(data)
  return(list("results" = data, "plot" = p))
}

# my signal function ####

#' My Signal Function
#'
#' Does things.
#' And returns things.
#'
#' @export
my_signal_function <- function(Risoe.object, LTYPE = "OSL", ID = NULL, set = NULL, run = NULL, position = NULL, S.interval = 6:10, BG.interval = 346:395) {

  if(is.null(ID)) { ID <- Risoe.object@METADATA$FNAME[1] }

  curveID <- Risoe.object@METADATA$ID[Risoe.object@METADATA$LTYPE == LTYPE & Risoe.object@METADATA$SET %in% set & Risoe.object@METADATA$RUN %in% run & Risoe.object@METADATA$POSITION %in% position]
  CTS <- lapply(curveID, function(x) Risoe.object@DATA[[x]])
  CTS <- do.call(cbind, CTS)
  BG <- colSums(CTS[BG.interval, ]) / length(BG.interval)
  n <- ncol(CTS)

  CTS.S <- CTS[c(1:(min(S.interval)-1), (max(BG.interval)+1):nrow(CTS)), ]

  CTS <- matrix(CTS, ncol = ncol(CTS)) - matrix(rep(BG, each = nrow(CTS)), ncol = ncol(CTS))
  CTS[c(1:(min(S.interval)-1), (max(BG.interval)+1):nrow(CTS)), ] <- CTS.S

  CTS.max <- sapply(1:ncol(CTS), function(x) { max(CTS[, x]) } )
  CTS <- matrix(CTS, ncol = ncol(CTS)) / matrix(rep(CTS.max, each = nrow(CTS)), ncol = ncol(CTS))

  CTS <- rowMeans(CTS)

  HIGH <- Risoe.object@METADATA$HIGH[[curveID[1]]]
  NPOINTS <- Risoe.object@METADATA$NPOINTS[[curveID[1]]]
  t <- seq(HIGH/NPOINTS, HIGH, by = HIGH/NPOINTS)

  data.frame("t" = t, "CTS.norm" = CTS, "n" = n, TYPE = LTYPE, "ID" = ID)
}
# my readerrate retriever ####

#' My readerrate retriever
#'
#' Does things.
#' And returns things.
#'
#' @export
my_readerrate <- function(dr.model, date, method = "average", data.path = "C:/Users/fhaba/OneDrive - Danmarks Tekniske Universitet/OSLdata/Calibration/calibration_info.csv", day.lim = 60){

  calib.file = read.csv(data.path)

  if(!(dr.model %in% calib.file$dr.model)) {
    print("chosen dr.model must be one of:")
    print(unique(calib.file$dr.model))
    stop("operation cancelled")
  }

  modeldata <- calib.file[calib.file$dr.model == dr.model,]

  if(method == "lm") {
    model <- lm(readerrate ~ as.Date(day), modeldata)
    prediction <- predict.lm(model, data.frame("day" = as.Date(date)), se.fit = T)
    output <- round(c(prediction[[1]], prediction[[2]]), 4)
    output = c(output, nrow(modeldata))
  }
  if(method == "exp") {
    model <- lm(log(readerrate) ~ as.Date(day), modeldata)
    prediction <- predict.lm(model, data.frame("day" = as.Date(date)), se.fit = T)
    output <- round(c(exp(prediction[[1]]), exp(prediction[[2]])), 4)
    output = c(output, nrow(modeldata))
  }
  if(method == "average") {
    modeldata = modeldata[modeldata$day <= as.Date(date) + day.lim & modeldata$day >= as.Date(date) - day.lim,]
    mean = mean(modeldata[, "readerrate"])
    if(nrow(modeldata) > 1) se = my_calc_sem(modeldata[, "readerrate"])
    if(nrow(modeldata) == 1) se = modeldata$readerrate.error
    output = round(c(mean, se), 4)
    output = c(output, nrow(modeldata))
  }

  output = data.frame("readerrate" = output[1], "se" = output[2], "n" = output[3], "method" = method)
  output
}

my_excel_copypaste <- function(object) {
  print(object)
  write.table(object, "clipboard", row.names = F, col.names = F, sep = "\t")
}

# my BayLum files creator ####

#'My write BayLum files
#'
#'Does things
#' And returns things
#'
#' @export
my_write_BayLum_files <- function(path,
                                  subsample.folder.names,
                                  Disc = NULL,
                                  DiscPos = NULL,
                                  DoseEnv, DoseEnv.error,
                                  DoseSource, DoseSource.error,
                                  begin.signal.integral, end.signal.integral,
                                  begin.background.integral, end.background.integral,
                                  inflatePercent = 0.025,
                                  nbOfLastCycleToRemove=2) {

  info = data.frame("sample" = subsample.folder.names,
             "DoseEnv" = DoseEnv,
             "DoseEnv.error" = DoseEnv.error,
             "DoseSource" = DoseSource,
             "DoseSource.error" = DoseSource.error,
             "begin.signal.integral" = begin.signal.integral,
             "end.signal.integral" = end.signal.integral,
             "begin.background.integral" = begin.background.integral,
             "end.background.integral" = end.background.integral,
             "inflatePercent" = inflatePercent,
             "nbOfLastCycleToRemove" = nbOfLastCycleToRemove)

  # setup folder paths
  folder = paste(path, subsample.folder.names, sep ="")
  for (j in 1:length(folder)){ dir.create(path = folder[j]) }
  BayLum_files_path = paste(folder,"/",sep="")

  if(!is.null(Disc)) {
    # Disc
    if(length(Disc) == 1) { Disc = rep(Disc, length(BayLum_files_path)) }
    for (j in 1:length(BayLum_files_path)) { write.csv(data.frame("position" = Disc[[j]]), paste(BayLum_files_path[j],"Disc.csv",sep=""), row.names = F) }
  }

  if(!is.null(DiscPos)) {
    # DiscPos
    if(length(DiscPos) == 1) { DiscPos = rep(DiscPos, length(BayLum_files_path)) }
    for (j in 1:length(BayLum_files_path)) { write.csv(data.frame("position" = DiscPos[[j]][,1], "grain" = DiscPos[[j]][,2]), paste(BayLum_files_path[j],"DiscPos.csv",sep=""), row.names = F) }
  }

  # DoseEnv
  for (j in 1:length(BayLum_files_path)) { write.csv(data.frame("obs"=info$DoseEnv[j] , "var" = info$DoseEnv.error[j]^2), paste(BayLum_files_path[j],"DoseEnv.csv",sep=""), row.names = F) }

  # DoseSource
  for (j in 1:length(BayLum_files_path)) { write.csv(data.frame("obs"=info$DoseSource[j] , "var" = info$DoseSource.error[j]^2), paste(BayLum_files_path[j],"DoseSource.csv",sep=""), row.names = F) }

  # rule
  for (j in 1:length(BayLum_files_path)) {
    write.csv(data.frame("[Param]" = c(
      paste("beginSignal=",info$begin.signal.integral[j], sep=" "),
      paste("endSignal=",info$end.signal.integral[j], sep=" "),
      paste("beginBackground=",info$begin.background.integral[j], sep=" "),
      paste("endBackground=",info$end.background.integral[j], sep=" "),
      paste("beginTest=",info$begin.signal.integral[j], sep=" "),
      paste("endTest=",info$end.signal.integral[j], sep=" "),
      paste("beginTestBackground=",info$begin.background.integral[j], sep=" "),
      paste("endTestBackground=",info$end.background.integral[j], sep=" "),
      paste("inflatePercent=",info$inflatePercent[j], sep=" "),
      paste("nbOfLastCycleToRemove=",info$nbOfLastCycleToRemove[j], sep=" ")
    ), check.names = FALSE), paste(BayLum_files_path[j],"rule.csv",sep=""), row.names = F, quote = F)
  }
}
