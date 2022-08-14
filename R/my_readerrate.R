# my_readerrate ####

#' my_readerrate
#'
#' @export
my_readerrate <- function(dr.model, date, method = "average", day.lim = 60, data.path = "C:/Users/fhaba/OneDrive - Danmarks Tekniske Universitet/OSLdata/Calibration/calibration.csv"){

  calib.file = read.csv(data.path)

  if(!(dr.model %in% calib.file$model)) {
    print("chosen dr.model must be one of:")
    print(unique(calib.file$model))
    stop("operation cancelled")
  }

  modeldata <- calib.file[calib.file$model == dr.model,]

  if(method == "lm") {
    model <- lm(readerrate ~ as.Date(day), modeldata)
    prediction <- predict.lm(model, data.frame("day" = as.Date(date)), se.fit = T)
    output <- round(c(prediction[[1]], prediction[[2]]), 4)
    output = c(output, nrow(modeldata))
  }
  if(method == "exp") {
    model <- lm(log(readerrate) ~ as.Date(day), modeldata)
    prediction <- predict.lm(model, data.frame("day" = as.Date(date)), se.fit = T)
    output <- round(c(exp(prediction[[1]]), (prediction[[2]])), 4)
    output = c(output, nrow(modeldata))
  }

  if(method == "average") {
    modeldata = modeldata[modeldata$date <= as.Date(date) + day.lim & modeldata$date >= as.Date(date) - day.lim,]
    mean = mean(modeldata[, "readerrate"])
    if(nrow(modeldata) > 1) se = my_calc_sem(modeldata[, "readerrate"])
    if(nrow(modeldata) == 1) se = modeldata$readerrate.error
    output = round(c(mean, se), 4)
    output = c(output, nrow(modeldata))
  }

  output = data.frame("readerrate" = output[1], "se" = output[2], "n" = output[3], "method" = method)
  output
}
