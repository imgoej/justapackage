# my_SAR_analyzer2 ####

#' my_SAR_analyzer2
#' @import ggplot2
#' @import Luminescence
#' @import stringr
#' @export
my_SAR_analyzer2 <- function(Risoe.object, position = NULL, grain = NULL, run = NULL, set = NULL, ltype = "OSL", dtype = NULL, protocol = "unknown", keep.empty = TRUE, txtProgressBar = FALSE,
                             signal.integral = c(6:10), background.integral = c(346:395), rejection.criteria = list(
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
                             fit.weights = FALSE,
                             fit.includingRepeatedRegPoints = FALSE,
                             plot = FALSE,
                             plot.single = FALSE,
                             onlyLxTxTable = FALSE,
                             readerrate = 1,
                             drdose = 1,
                             label = "sample",
                             recuprecyc = TRUE,
                             plot.drc = FALSE
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

  func.results = list(results = results, data = cbind(info, De.data, DR, STn20, Recyc_Depl_Recup, ID, label))

  if(plot.drc == TRUE) {
    # select drc to plot
    drc.data = func.results$results$LnLxTnTx.table[func.results$results$LnLxTnTx.table$UID == unique(func.results$results$LnLxTnTx.table$UID)[1],]
    # extract formula
    formula = as.formula(paste("y ~",as.character(func.results$results$formula[[1]])))

    # define start points and end points for lines to be drawn (from Natural on y axis to Dose on x-axis)
    drc.natural.lxtx = drc.data$LxTx[drc.data$Name == "Natural"]
    drc.natural.lxtx.error = drc.data$LxTx.Error[drc.data$Name == "Natural"]
    drc.de = func.results$data$De[1]/readerrate
    drc.de.error = func.results$data$De.Error[1]/readerrate

    # design new data table for LnLxTnTx data
    drc.data$Group[drc.data$Repeated == "TRUE"] = "Recycle point(s)"
    drc.data$Group[drc.data$Repeated == "FALSE"] = "Regenerative point(s)"
    drc.data$Group[drc.data$Name == "R0"] = "Recuperation point(s)"
    drc.data$Group[drc.data$Name == "Natural"] = "Natural"
    drc.data$Group = factor(drc.data$Group,
                            levels = unique(drc.data$Group))

    # specify shape and color attributes to categories of data points
    drc.attr = data.frame("Group" = c("Natural","Regenerative point(s)","Recuperation point(s)","Recycle point(s)"),
                          "Shape" = c(10,23,24,23),
                          "Color" = c("black","black","white","white"))

    # select only attributes for points that are actually present
    drc.attr = drc.attr[drc.attr$Group %in% drc.data$Group,]

    # plot drc
    p <- ggplot(drc.data, aes(x=Dose,y=LxTx,fill=Group,shape=Group)) +
      stat_function(fun = function(x) eval(func.results$results$formula[[1]]), inherit.aes = F) +
      geom_errorbar(aes(ymin = LxTx-LxTx.Error, ymax = LxTx+LxTx.Error), width = (max(drc.data$Dose)-min(drc.data$Dose))*0.025) +
      geom_point(size = 3) +
      scale_x_continuous(breaks = drc.data$Dose, labels = round(drc.data$Dose*readerrate, 0)) +
      scale_shape_manual(values = drc.attr$Shape) + scale_fill_manual(values = drc.attr$Color) +
      geom_segment(aes(x = 0, # de plotting x1 to x2
                       xend = drc.de,
                       y = drc.natural.lxtx,
                       yend = drc.natural.lxtx),
                   color = "red", linetype = "dashed", size = 0.4)+
      geom_segment(aes(x = drc.de, # de plotting y1 to y2
                       xend = drc.de,
                       y = 0,
                       yend = drc.natural.lxtx),
                   color = "red", linetype = "dashed", size = 0.4)  +
      geom_segment(aes(x = 0, # de upper error plotting y1 to y2
                       xend = drc.de+drc.de.error,
                       y = drc.natural.lxtx+drc.natural.lxtx.error,
                       yend = drc.natural.lxtx+drc.natural.lxtx.error),
                   color = "red", linetype = "dashed", size = 0.25) +
      geom_segment(aes(x = drc.de+drc.de.error, # de upper error plotting y1 to y2
                       xend = drc.de+drc.de.error,
                       y = 0,
                       yend = drc.natural.lxtx+drc.natural.lxtx.error),
                   color = "red", linetype = "dashed", size = 0.25) +
      geom_segment(aes(x = 0, # de lower error plotting x1 to x2
                       xend = drc.de-drc.de.error,
                       y = drc.natural.lxtx-drc.natural.lxtx.error,
                       yend = drc.natural.lxtx-drc.natural.lxtx.error),
                   color = "red", linetype = "dashed", size = 0.25) +
      geom_segment(aes(x = drc.de-drc.de.error, # de lower error plotting y1 to y2
                       xend = drc.de-drc.de.error,
                       y = 0,
                       yend = drc.natural.lxtx-drc.natural.lxtx.error),
                   color = "red", linetype = "dashed", size = 0.25) +
      xlab("Regenerative Dose (Gy)") + ylab(bquote("L"["x"]*"/T"["x"])) + theme_classic() +
      theme(legend.title = element_blank())


    # normalize
    LxTx.norm.faktor = 1 / as.numeric(str_extract_all(func.results$results$formula[[1]],"\\(?[0-9,.]+\\)?")[[1]][1])
    Dose.norm.faktor = 1 / max(drc.data$Dose)

    drc.data$LxTx = drc.data$LxTx*LxTx.norm.faktor
    drc.data$LxTx.Error = drc.data$LxTx.Error*LxTx.norm.faktor

    # normalize lxtx dataset
    norm.drc = cbind.data.frame(
      "dose"=drc.data$Dose*readerrate,
      "LxTx"=drc.data$LxTx,
      "LxTx.Error"=drc.data$LxTx.Error,
      "Group"=drc.data$Group,
      "UID"=drc.data$UID
    )
    func.results$norm.lxtx = norm.drc
    func.results$formula = formula

    return(list("plot" = p, "results" = func.results)) } else{
      return(func.results)
    }
}
