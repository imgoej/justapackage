### autorun_BayLum ####

#' autorun_BayLum
#'
#' autorun_BayLum() takes the output of "Generate_DataFile()", "Generate_DataFile_MG()" or "Create_DataFile()" - which may have been build on several samples - and runs BayLum on each sample individually. \cr
#' BayLum is run until the Rubin-Gelman statistics are below 1.050 for "A", "D" and "sD" parameters for a sample. If a particular run fails to converge, the function will proceed to double the number of "Iter" used. The doubling proceeds until convergence is reached. \cr
#' The function then stitches together the "A"-parameter MCMC-samples for each sample into one csv-file (this is the input for Age-Depth model of the ArchaeoPhases-package). The same is true for the "D"-parameter.
#'
#' @import BayLum
#' @import Luminescence
#' @import coda
#' @import runjags
#'@param DataFile The output of "Generate_DataFile()", "Generate_DataFile_MG" or "Create_DataFile" of the BayLum R-package.
#'@param SampleNames A vector of names to be attached to attached to each sample.
#'@param BinPerSample A vector of numbers indicating how many binx-files exist for samples (in order of samples structered by the DataFile-object). BinPerSample = c(1,2,1) would indicate that sample 1 and sample 3 have 1 binx-file each. Sample 2 as 2 binx-files.
#'@param csv_name The prefix-label to put on each produced csv-file. csv_name = "example" would create csv-files named: "example_MCMC_A.csv" and "example_MCMC_D.csv".
#'@param Iter (with default = 10000) The number of sampling iterations that the MCMC should run for (an equal number of burn-in iterations is defined).
#'@param create_new_files (with default = TRUE) Whether all new MCMC csv-files should be created or not. If "FALSE", the function will built on presumably already existing csv_files with labels that match csv_name and the suffix. This could be useful secondary runs are needed (for example if new samples are added after having run this function already)
#'@param Origin_fit (with default = TRUE) Should DRC-fits be forced through zero?
#'@param LIN_fit (with default = FALSE) Should DRC-fits be a saturating exponential plus a linear component?
#'@param distribution (with default = gaussian) Which dose-dispersion model within BayLum to use? (possible are "gaussian", "lognormal_A", "lognormal_M" and "cauchy")
#'@param (with default = c(1, 100)) The prior for the age. User can only specify one interval which is forced on all samples. So user must consider the most extreme samples. Since samples are run one at a time, this should not have major ramifications on convergence-time.
#'@return This function returns a summary of samples that have been run along with a status-marker. More importantly, csv-files are written - one for the "A" parameter, one for the "D" parameter - in which converged output MCMC samples are sticthed together in the original sample order dictated by the DataFile-input. Credible intervals can then be computed from them or they can be used in Age-depth models. A csv-file with Rubin-Gelman statistics for all samples is also produced.
#'@export

autorun_BayLum <- function(DataFile, SampleNames, BinPerSample, csv_name, Iter = 10000, create_new_files = TRUE, Origin_fit = TRUE, LIN_fit = FALSE, distribution = "gaussian", PriorAge = c(1, 100)) {
  # Ensure min number of iterations (limitation of Age_Computation()) in order to ensure that chains lengths are the same between samples
  Iter <- max(Iter, 10000)

  if(create_new_files) {
    # create storage files
    write.csv(data.frame("Labcode" = SampleNames, "RG_A" = NA, "RG_D" = NA, "RG_sD" = NA, "nbins" = NA), paste(csv_name, "_RuGel.csv",sep = ""),row.names = FALSE)
    write.csv(data.frame("iteration" = 1:9999), paste(csv_name, "_MCMC_A.csv",sep = ""), row.names = FALSE)
    write.csv(data.frame("iteration" = 1:9999), paste(csv_name, "_MCMC_D.csv",sep = ""), row.names = FALSE)
  }
  # create function to extract DataFile for only one sample
  DataFile_extractor_function <- function(DataFile, whichbins, whichqueue) {
    DataFile.new <- DataFile
    DataFile.new$LT <- DataFile$LT[whichqueue]
    DataFile.new$sLT <- DataFile$sLT[whichqueue]
    DataFile.new$ITimes <- DataFile$ITimes[whichqueue]
    DataFile.new$dLab <- matrix(DataFile$dLab[,whichbins], nrow = 2)
    DataFile.new$ddot_env <- matrix(DataFile$ddot_env[,whichbins], nrow = 2)
    DataFile.new$regDose <- DataFile$regDose[whichqueue]
    DataFile.new$J <- DataFile$J[whichbins]
    DataFile.new$K <- DataFile$K[whichbins]
    DataFile.new$Nb_measurement <- DataFile$Nb_measurement[whichbins]
    DataFile.new
  }

  BinPerSample.df <- data.frame(SampleNames, BinPerSample)

  # add column to keep track of cumulative bin number
  BinPerSample.df$bcumsum <- cumsum(BinPerSample)

  # create list for each sample with index for where sample is located in DataFile
  whichbinnumber.LIST <- lapply(BinPerSample.df$Sample, function(s){
    BinPerSample.df_temp <- BinPerSample.df[BinPerSample.df$Sample == s,]
    if(BinPerSample.df_temp$BinPerSample == 1) {
      BinPerSample.df_temp$bcumsum
    } else {
      (BinPerSample.df_temp$bcumsum-BinPerSample.df_temp$BinPerSample+1):BinPerSample.df_temp$bcumsum
    }
  })

  # add name to list of sample index in DataFile
  names(whichbinnumber.LIST) <- paste("s",BinPerSample.df$Sample,sep="")

  # For each sample, run BayLum till completion
  A <- lapply(1:length(whichbinnumber.LIST), function(x) {
    # read storage files
    csv_validation <- read.csv(paste(csv_name, "_RuGel.csv",sep = ""),header = T)
    MCMC_A <- read.csv(paste(csv_name, "_MCMC_A.csv",sep = ""),header = T)
    MCMC_D <- read.csv(paste(csv_name, "_MCMC_D.csv",sep = ""),header = T)

    # extract DataFile for sample x
    targetbins <- whichbinnumber.LIST[[x]]
    DataFile_extracted <- DataFile_extractor_function(DataFile, whichbins = targetbins, whichqueue = x)

    # create Rubin-Gelman object
    RuGel <- list("psrf" = data.frame(c(10,10,10),c(10,10,10)))

    # initialize number of iterations to run for
    Iter.count <- Iter

    # run BayLum
    while(sum(RuGel$psrf[, 2] < 1.050) != length(RuGel$psrf[,2])) {

      Bay.result <- BayLum::Age_Computation(
        DATA = DataFile_extracted,
        SampleName = names(whichbinnumber.LIST)[x],
        PriorAge = PriorAge,
        BinPerSample = length(DataFile_extracted$J),
        LIN_fit = LIN_fit,
        Origin_fit = Origin_fit,
        distribution = distribution,
        Iter = Iter.count,
        t = 3
      )
      RuGel <- coda::gelman.diag(Bay.result$Sampling, multivariate = FALSE)
      Iter.count = Iter.count*2
    }

    # add Rubin-Gelman results to file
    csv_validation[x, ] <- c(names(whichbinnumber.LIST)[x], RuGel$psrf[1,2],RuGel$psrf[2,2], RuGel$psrf[3,2], length(DataFile_extracted$J))

    # extract marginal posterios
    MCMC.output <- Bay.result$Sampling
    MCMC.output <- runjags::combine.mcmc(MCMC.output)
    MCMC.output_A <- as.data.frame(MCMC.output[,1])
    MCMC.output_D <- as.data.frame(MCMC.output[,2])
    names(MCMC.output_A) <- names(whichbinnumber.LIST)[x]
    names(MCMC.output_D) <- names(whichbinnumber.LIST)[x]

    MCMC_A <- coda::as.mcmc(cbind.data.frame(MCMC_A, MCMC.output_A))
    MCMC_D <- coda::as.mcmc(cbind.data.frame(MCMC_D, MCMC.output_D))


    write.csv(csv_validation, paste(csv_name, "_RuGel.csv",sep = ""),row.names = FALSE)
    write.csv(MCMC_A, paste(csv_name, "_MCMC_A.csv",sep = ""), row.names = FALSE)
    write.csv(MCMC_D, paste(csv_name, "_MCMC_D.csv",sep = ""), row.names = FALSE)

    status = data.frame("labcode" = names(MCMC.output_A), status = "OK")
    return(status)
  })
  print(do.call(rbind.data.frame, A))
}
