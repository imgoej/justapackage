# my_D0_filter ####

#' my_D0_filter
#' @import Luminescence
#' @export
my_D0_filter <- function(OSLdata, range,method ="CAM_UL", ...) {
  colnames(OSLdata) <- c("De", "De.Error", "D0")
  OSLdata <- OSLdata[complete.cases(OSLdata[, c("De", "De.Error", "D0")]), ]
  if(method == "CAM") {
    dose.model <- lapply(range, function(x) {
      OSLdata.tmp <- OSLdata[OSLdata[, "D0"] >= x, ]
      if(nrow(OSLdata.tmp) <= 1) CAM.tmp <- NA
      else {
        CAM.tmp <- get_RLum(calc_CentralDose(   # --
          OSLdata.tmp[, c("De", "De.Error")],       #   |
          log = TRUE,                          #   |
          sigmab = 0,                           #   | --- CAM settings
          na.rm = TRUE,                         #   |
          plot = FALSE,                         #   |
          verbose = FALSE))[1:4]                # --
        return(CAM.tmp)
      }})
  }
  if(method == "CAM_UL") {
    dose.model <- lapply(range, function(x) {
      OSLdata.tmp <- OSLdata[OSLdata[, "D0"] >= x, ]
      if(nrow(OSLdata.tmp) <= 1) CAM.tmp <- NA
      else {
        CAM_UL.tmp <- get_RLum(calc_CentralDose(   # --
          OSLdata.tmp[, c("De", "De.Error")],       #   |
          log = FALSE,                          #   |
          sigmab = 0,                           #   | --- CAM settings
          na.rm = TRUE,                         #   |
          plot = FALSE,                         #   |
          verbose = FALSE))[1:4]                # --
        return(CAM_UL.tmp)
      }})
  }
  if(method == "ADM") {
    dose.model <- lapply(range, function(x) {
      OSLdata.tmp <- OSLdata[OSLdata[, "D0"] >= x, ]
      if(nrow(OSLdata.tmp) <= 1) ADM.tmp <- NA
      else {
        ADM.tmp <- get_RLum(calc_AverageDose(   # --
          OSLdata.tmp[, c("De", "De.Error")],       #   |
          ...))[1:4]                # --
        return(ADM.tmp)
      }})  }
  dose.model <- do.call(rbind, dose.model)
  dose.model$method = method
  dose.model <- cbind.data.frame("minD0" = range, dose.model)
  cbind.data.frame(dose.model, "n" = sapply(range, function(x) nrow(OSLdata[OSLdata[, "D0"] >= x, ])))
}
