### analyst_to_BayFolders ####

#' analyst_to_BayFolders
#'
#' Writes files necessary to use Generate_DataFile() function of the BayLum R-package.
#'
#' @import BayLum
#' @import Luminescence
#' @import data.table
#'@param data Data.frame containing information about each disc and grain to be used in BayLum. This function is optimized for use with the output of Analyst. Importantly, the data.frame must contain columns: "Labcode", "DoseEnv", "DoseEnv_error", "DoseSource", "DoseSource_error", "NbOfLastCycleToRemove", "Filename", "Disc#", "Grain#", "Sig1", "Sig2", "BG1", "BG2". 
#'@param folder.name Specify a name for a folder to be created in the current working directory which will contain the output of this function (defaults to "BayLum_files").
#'@param inflatePercent Specify measurement uncertainty that should be incorporated (defaults to 0.025). 
#'@return Returns a list of vectors containing the required input for BayLum functions Generate_DataFile() and AgeS_Computation(). 
#'@export
analyst_to_BayFolders <- function(data, 
                                  folder.name = "BayLum_files", 
                                  inflatePercent = 0.025) {
  # create storage object
  WriteBF.output = data.frame()
  Disc.list = list()
  
  # quick check that dataset is proper
  graincheck <- sapply(1:nrow(data), function(x) {
    
    if(is.na(data$`Grain#`[x])) {
      check <- FALSE
    } else {
      check <- data$`Grain#`[x] >= 0
    }
    check
  } )
  
  if(sum(graincheck) < nrow(data)) {
    cat(paste("I detect rows which don't have zero or above zero numeric input for \"Grain#\"\nHere are the infractions:\n\n"))
    faildata_2 <- data[which(!graincheck), c("Labcode","Filename","Disc#","Grain#")]
    faildata_1 <- data.frame("Excel_row" = which(!graincheck))
    print(cbind.data.frame(faildata_1, faildata_2))
    stop("Execution halted")
    
  }
  
  # extract baylum info
  for (x in unique(data$Labcode)) {
    # select Labcode
    data.temp <- data[data$Labcode == x, ]
    
    # make a datatable (from datatable library) - this is to keep order between values and filenames
    data.temp.datatable <- as.data.table(data.temp)
    
    # extract Filenames
    Filenames <- unique(data.temp$Filename)
    
    # extract basic baylum info
    DoseEnv <- unique(data.temp$DoseEnv[data.temp$Labcode == x])
    DoseEnv_error <- unique(data.temp$DoseEnv_error[data.temp$Labcode == x])
    
    # Error check: DoseEnv
    if(length(DoseEnv) > 1) {
      stop(paste("Labcode",x,"has more than one unique entry for \"DoseEnv\""))
    }
    
    if(length(DoseEnv_error) > 1) {
      stop(paste("Labcode",x,"has more than one unique entry for \"DoseEnv_error\""))
    }
    
    # Error check: DoseSource
    # I aggregate a table of DoseSource by filename. Then I check that there is only one entry.  
    if(sum(1 < table(data.temp.datatable[, unique(DoseSource), by=Filename]$Filename)) > 0) {
      cat(paste("Oops. Something wrong. Function will be halted.\nHere is a table of binx-files for Labcode",x,"along with DoseSource input:\n\n"))
      print(data.temp.datatable[, unique(DoseSource), by=Filename])
      stop(paste("Labcode",x,"has at least one binx file with multiple entries for \"DoseSource\""))
    }
    
    # Error check: DoseSource_error
    # same as above, but for DoseSource_error
    if(sum(1 < table(data.temp.datatable[, unique(DoseSource_error), by=Filename]$Filename)) > 0) {
      cat(paste("Oops. Something wrong. Function will be halted.\nHere is a table of binx-files for Labcode",x,"along with DoseSource_error input:\n\n"))
      print(data.temp.datatable[, unique(DoseSource_error), by=Filename])
      stop(paste("Labcode",x,"has at least one binx file with multiple entries for \"DoseSource_error\""))
    }
    
    # If there is no DoseEnv, I put NA. Don't remember why.
    if(length(DoseEnv) == 0) {
      DoseEnv <- NA
      DoseEnv_error <- NA
    }
    
    # Extract reader rates
    DoseSource <- data.temp.datatable[, unique(DoseSource), by=Filename]$V1
    DoseSource_error <- data.temp.datatable[, unique(DoseSource_error), by=Filename]$V1
    
    
    if(length(DoseSource) == 0) {
      DoseSource <- NA
      DoseSource_error <- NA
    }
    
    # Log number of binfiles per labcode
    BinPerLabcode <- 1:length(unique(data.temp$Filename))
    
    # multi-grain
    for (b in unique(data.temp$Filename)) {
      #multi-grain
      if(sum(data.temp$`Grain#`[data.temp$Filename == b]) == 0) {
        NewList <- list(data.frame("position" = data.temp$`Disc#`[data.temp$Filename == b]))
        Disc.list <- c(Disc.list, NewList)
      }
      #single-grain
      if(sum(data.temp$`Grain#`[data.temp$Filename == b]) > 0) {
        NewList <- list(data.frame("position" = data.temp$`Disc#`[data.temp$Filename == b],
                                   "grain" = data.temp$`Grain#`[data.temp$Filename == b]))
        Disc.list <- c(Disc.list, NewList)
      }
    }
    
    # write to WriteBF.output.temp for labcode x
    WriteBF.output.temp <- data.frame("Labcode" = x, DoseEnv, DoseEnv_error, DoseSource, DoseSource_error, "BinPerLabcode" = max(BinPerLabcode), "Bin" = BinPerLabcode, "NbLastCycle" = data.temp.datatable[, unique(NbOfLastCycleToRemove), by=Filename]$V1, "Filename" = Filenames)
    
    # determine if target BayLum folder will have subfolders or not
    if(length(WriteBF.output.temp$Labcode) == 1) {
      WriteBF.output.temp$BayLPath <- x
    } else {
      WriteBF.output.temp$BayLPath <- paste(x,WriteBF.output.temp$Bin,sep="/")
    }
    
    # Write signal summation data
    WriteBF.output.temp$Sig1 <- data.temp.datatable[, unique(Sig1), by=Filename]$V1
    WriteBF.output.temp$Sig2 <- data.temp.datatable[, unique(Sig2), by=Filename]$V1
    WriteBF.output.temp$Bg1 <- data.temp.datatable[, unique(BG1), by=Filename]$V1
    WriteBF.output.temp$Bg2 <- data.temp.datatable[, unique(BG2), by=Filename]$V1
    
    # write to main object
    WriteBF.output <- rbind.data.frame(WriteBF.output, na.omit(WriteBF.output.temp))
  }
  
  # 
  DT <- as.data.table(WriteBF.output)
  
  BinPerLabcode_in_order <- DT[, max(BinPerLabcode), by=Labcode]
  BayLPath_in_order <- WriteBF.output$BayLPath
  SampleNames_in_order <- DT[, unique(Labcode), by=Labcode]
  
  # Write BayLum files
  BayLum::write_BayLumFiles(
    folder = folder.name,
    SampleNames = unique(WriteBF.output$Labcode),
    BinPerSample = as.numeric(BinPerLabcode_in_order$V1),
    DiscPos = Disc.list,
    DRenv = WriteBF.output$DoseEnv,
    DRenv.error = WriteBF.output$DoseEnv_error,
    DRsource = WriteBF.output$DoseSource,
    DRsource.error = WriteBF.output$DoseSource_error,
    signal.integral.min = WriteBF.output$Sig1,
    signal.integral.max = WriteBF.output$Sig2,
    background.integral.min = WriteBF.output$Bg1,
    background.integral.max = WriteBF.output$Bg2,
    inflatePercent = inflatePercent,
    nbOfLastCycleToRemove = WriteBF.output$NbLastCycle
  )
  
  # write text to specify filename for each
  folders <- paste(folder.name,"/",WriteBF.output$BayLPath, sep = "")
  
  for (i in 1:length(folders)){
    write.csv(NA, paste(folders[i],"/","NEEDS BINX-FILE---",WriteBF.output$Filename[i], ".txt",sep = ""))
  }
  
  output <- list(
    "SampleNames" = SampleNames_in_order$V1,
    "BinPerSample" = BinPerLabcode_in_order$V1,
    "BinxPaths" = BayLPath_in_order
  )
  
  saveRDS(output, paste(folder.name,"/","input_for_BayLum.RDS",sep = ""))
  
  output
}
