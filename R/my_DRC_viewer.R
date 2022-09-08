# my_DRC_viewer ####

#' my_DRC_viewer
#' @import ggplot2
#' @import Luminescence
#' @import stringr
#' @import cowplot
#' @export
my_DRC_viewer <- function(file, my_SAR_analyzer_output,index, readerrate){
  # ggplot2 settings
  theme_ggplot2 = theme_classic() +
    theme(text = element_text(family = "serif"),
          legend.title = element_blank(),
          axis.text.y=element_text(size = 12),
          axis.text.x=element_text(size = 12),
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),panel.border = element_rect(fill = NA))
  # open pdf graphic device
  pdf(file = file, onefile = TRUE,height = 4, width = 11)

  # create list of all DRCs
  LxTxTnTx.table.IDs <- unique(my_SAR_analyzer_output$results$LnLxTnTx.table$UID)

  # select DRCs to evaluate
  DRCs.chosen <- as.numeric(rownames(my_SAR_analyzer_output$data[index,]))

  # loop to create all plots
  for (i in DRCs.chosen) {
    # isolate LnLxTnTx-table for chosen DRC
    drc.data = my_SAR_analyzer_output$results$LnLxTnTx.table[my_SAR_analyzer_output$results$LnLxTnTx.table$UID == LxTxTnTx.table.IDs[[i]],]
    # extract formula
    formula = as.formula(paste("y ~",as.character(my_SAR_analyzer_output$results$formula[[i]])))
    # define start points and end points for lines to be drawn (from Natural on y axis to Dose on x-axis)
    drc.natural.lxtx = drc.data$LxTx[drc.data$Name == "Natural"]
    drc.natural.lxtx.error = drc.data$LxTx.Error[drc.data$Name == "Natural"]
    drc.de = my_SAR_analyzer_output$data$De[i] / readerrate
    drc.de.error = my_SAR_analyzer_output$data$De.Error[i] / readerrate
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
    plotDRC <- ggplot(drc.data, aes(x=Dose,y=LxTx,fill=Group,shape=Group)) +
      stat_function(fun = function(x) eval(my_SAR_analyzer_output$results$formula[[i]]), inherit.aes = F) +
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
      xlab("Regenerative Dose (Gy)") + ylab(bquote("L"["x"]*"/T"["x"])) + theme_ggplot2 +
      theme(legend.title = element_blank(),legend.position = c(0.98,0.02), legend.justification = c(0.98,0.02), legend.background = element_rect(fill = "transparent")) +
      ggtitle(label = paste("index = ", i, " |  Disc = ", my_SAR_analyzer_output$data$Disc[i], " |  Grain = ", my_SAR_analyzer_output$data$Grain[i]))

    # PLOT EVOLUTION OF Tx/Tn ratio
    TnTx.ratio <- drc.data$Net_TnTx / drc.data$Net_TnTx[drc.data$Group == "Natural"]

    plotTnTx <- ggplot(drc.data, aes(x = Name, y = TnTx.ratio)) +
      geom_point(size = 2) +
      ylab(bquote("T"["x"] * "/" * "T"["n"])) +
      xlab("Step") +
      ylim(0, 2) +
      theme_ggplot2

    # PLOT EVOULUTIN OF Tx COUNTS PER REG GY GIVEN
    plotTnTxcounts <- ggplot(drc.data, aes(x = Name, y = TnTx)) +
      geom_bar(stat = "identity") +
      ylab(bquote("T"["x"]*"T"["n"]*" (counts)")) +
      xlab("Step") +
      theme_ggplot2

    print(plot_grid(
      plotDRC, plotTnTx, plotTnTxcounts, align = "hv", nrow = 1
    ))

  }
  dev.off()
}
