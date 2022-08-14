# my_signal_function ####

#' my_signal_function
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
