#' XCMS peak identification (part1)
#'
#' @param pd Data frame with meta data (inc. a column called "File" with path to file for analysis). Assumes that the data to be analysed are centroided and on disk. 
#' @param peakwidth For CentWaveParam. numeric(2) with the expected approximate peak width in chromatographic space. Given as a range (min, max) in seconds.
#' @param ppm For CentWaveParam. numeric(1) defining the maximal tolerated m/z deviation in consecutive scans in parts per million (ppm) for the initial ROI definition.
#' @param noise For CentWaveParam. numeric(1) allowing to set a minimum intensity required for centroids to be considered in the first analysis step (centroids with intensity < noise are omitted from ROI detection).
#' @param snthresh For CentWaveParam. numeric(1) defining the signal to noise ratio cutoff.
#' @param integrate For CentWaveParam. Integration method. For integrate = 1 peak limits are found through descent on the mexican hat filtered data, for integrate = 2 the descent is done on the real data. The latter method is more accurate but prone to noise, while the former is more robust, but less exact.
#' @param prefilter For CentWaveParam. numeric(2): c(k, I) specifying the prefilter step for the first analysis step (ROI detection). Mass traces are only retained if they contain at least k peaks with intensity >= I.
#' @param mzdiff For CentWaveParam. numeric(1) representing the minimum difference in m/z dimension required for peaks with overlapping retention times; can be negative to allow overlap. During peak post-processing, peaks defined to be overlapping are reduced to the one peak with the largest signal.
#' @param fitgauss For CentWaveParam. logical(1) whether or not a Gaussian should be fitted to each peak. This affects mostly the retention time position of the peak.
#' @param mzCenterFun For CentWaveParam. Name of the function to calculate the m/z center of the chromatographic peak. Allowed are: "wMean": intensity weighted mean of the peak’s m/z values, "mean": mean of the peak’s m/z values, "apex": use the m/z value at the peak apex, "wMeanApex3": intensity weighted mean of the m/z value at the peak apex and the m/z values left and right of it and "meanApex3": mean of the m/z value of the peak apex and the m/z values left and right of it.
#' @param ReFine For FilterIntensityParam. logical(1) whether or not peak intensity filtering should be performed (default = FALSE).
#' @param thr For FilterIntensityParam. numeric(1) defining the minimal required intensity for a peak to be retained. Defaults to threshold = 0.
#' @param Align For PeakDensityParam in alignment. logical(1) whether or not peak alignment should be performed (default = FALSE).
#' @param minFraction For PeakDensityParam in alignment. numeric(1) between 0 and 1 defining the minimum required fraction of samples in which peaks for the peak group were identified. Peak groups passing this criteria will aligned across samples and retention times of individual spectra will be adjusted based on this alignment. 
#' For minFraction = 1 the peak group has to contain peaks in all samples of the experiment. Note that if subset is provided, the specified fraction is relative to the defined subset of samples and not to the total number of samples within the experiment (i.e. a peak has to be present in the specified proportion of subset samples).
#' @param subsetAdjust For PeakDensityParam in alignment. Options: "previous" or "average". How the non-subset samples are adjusted bases also on the parameter subsetAdjust: with subsetAdjust = "previous", each non-subset sample is adjusted based on the closest previous subset sample which results in most cases with adjusted retention times of the non-subset sample being identical to the subset sample on which the adjustment bases. 
#' The second, default, option is to use subsetAdjust = "average" in which case each non subset sample is adjusted based on the average retention time adjustment from the previous and following subset sample. For the average a weighted mean is used with weights being the inverse of the distance of the non-subset sample to the subset samples used for alignment.
#' @param sampleGroupsType For PeakDensityParam in alignment. The column which the alignment models should be estimated. 
#' @param sampleGroupsType_sel For PeakGroupsParam in alignment. Type of samples in specified sampleGroupsType column which the alignment models should be estimated. Samples not part of the subset are adjusted based on the closest subset sample. See description xcms for more details (default = "QC")
#' Samples not part of this subset are left out in the estimation of the alignment models, but their retention times are subsequently adjusted based on the alignment results of the closest sample in subset (close in terms of position within the object). 
#' Alignment could thus be performed on only real samples leaving out e.g. blanks, which are then in turn adjusted based on the closest real sample. Here it is up to the user to ensure that the samples within object are ordered correctly (e.g. by injection index).
#' @param span For PeakDensityParam in alignment. numeric(1) defining the degree of smoothing (if smooth = "loess"). This parameter is passed to the internal call to loess.
#' @param smooth For PeakDensityParam in alignment. ncharacter defining the function to be used, to interpolate corrected retention times for all peak groups. Either "loess" or "linear".
#' @param nCore Number of cores for parallelisation. Numeric(1), default is 1.   
#' 
#' @note The user is able to set most xcms function parameters in this script. However, all xcms function parameters might not be included (i.e., is currently using the default settings). Check the manual (https://bioconductor.org/packages/release/bioc/manuals/xcms/man/xcms.pdf) and adjust this script in case you want to modify additional xcms parameters. 
#' @examples https://gitlab.com/stanstrup-teaching/XCMS-course/-/blob/master/1-XCMS.Rmd
#' @examples https://bioconductor.org/packages/devel/bioc/vignettes/xcms/inst/doc/xcms.html
#' 
XCMS_peak_identification <- function(pd, peakwidth = c(20, 50), ppm = 25, noise = 0, snthresh = 10, integrate = 1, prefilter = c(3, 100), mzdiff = -0.001,  fitgauss = FALSE, mzCenterFun = "wMean", ReFine = FALSE, thr = 0, Align = FALSE, minFraction = 0.85, span = 0.4, subsetAdjust = "average", sampleGroupsType = as.vector(), sampleGroupsType_sel = "QC", smooth = "loess", nCore = 1){
  
  # Load packaged
  suppressPackageStartupMessages({
    require(xcms)
    require(magrittr)
    require(dplyr)
    require(tidyverse)
    require(BiocParallel)
    })

  # Check on mandatory 'file' column 
  if("file" %in% colnames(pd)){
    message("Renaming 'file' column to 'File'")
    pd <- pd %>% rename(., File = "file")
  }
  if(!"File" %in% colnames(pd)){
    stop("'File' column missing")
  }
  
  # Load data 
  message("Load files")
  data <- readMSData(pd$File, 
                     pdata = new("NAnnotatedDataFrame", pd),
                     mode = "onDisk", 
                     centroided. = TRUE)  

  # Set CentWave parameters and find peaks
  cwp <- CentWaveParam()
  cwp <- CentWaveParam(peakwidth = peakwidth, 
                       ppm = ppm, 
                       noise = noise, 
                       snthresh = snthresh, 
                       integrate = integrate, 
                       prefilter = prefilter, 
                       mzdiff = mzdiff,
                       mzCenterFun = mzCenterFun,
                       fitgauss = fitgauss) 
  
  message("XCMS peak identification")
  xdata <- findChromPeaks(data, 
                          param = cwp,
                          BPPARAM = BiocParallel::MulticoreParam(workers = nCore)) 
  
  # Refine peaks
  if(ReFine){
    message("Refine peaks")
    inf <- FilterIntensityParam(threshold = thr)
    xdata_refined <- refineChromPeaks(xdata,
                                      param = inf,
                                      BPPARAM = BiocParallel::MulticoreParam(workers = nCore))
  }
  
  if(Align & ReFine){
      message("Align refined peaks")
      # Initial peak grouping. Use sample_type as grouping variable
      pdp_subs <- PeakDensityParam(sampleGroups = sampleGroupsType, 
                                   minFraction = minFraction)
      xdata_align <- groupChromPeaks(xdata, 
                                     param = pdp_subs)
      xdata_refined_align <- groupChromPeaks(xdata_refined, 
                                        param = pdp_subs)
    
      # Define subset-alignment options and perform the alignment
      pgp_subs <- PeakGroupsParam(minFraction = minFraction, 
                                  subset = which(sampleGroupsType == sampleGroupsType_sel),
                                  subsetAdjust = subsetAdjust, 
                                  span = span, 
                                  smooth = smooth)
      xdata_align <- adjustRtime(xdata_align, 
                                 param = pgp_subs)
      xdata_refined_align <- adjustRtime(xdata_refined_align, 
                                    param = pgp_subs)
      } else if(Align & !ReFine){
        message("Align peaks")
        # Initial peak grouping. Use sample_type as grouping variable
        pgp_subs <- PeakDensityParam(sampleGroups = sampleGroupsType, 
                                     minFraction = minFraction)
        xdata_align <- groupChromPeaks(xdata, 
                                       param = pgp_subs)

        # Define subset-alignment options and perform the alignment
        pgp_subs <- PeakGroupsParam(minFraction = minFraction,
                                    subset = which(sampleGroupsType == sampleGroupsType_sel),
                                    subsetAdjust = subsetAdjust, 
                                    span = span, 
                                    smooth = smooth)
        xdata_align <- adjustRtime(xdata_align, 
                                   param = pgp_subs)
      }

  # Store results 
  if(Align & ReFine){
  list(
    "Raw_data" = data,
    "Peak_data" = xdata,
    "Peak_data_refined" = xdata_refined,
    "Peak_data_align" = xdata_align,
    "Peak_data_refined_align" = xdata_refined_align,
    "PeakGroups_Param_subset" = pdp_subs,
    "PeakGroups_Param_subset_selected" = sampleGroupsType_sel,
    "PeakGroups_Alignment_Param" = pgp_subs, 
    "Meta_info" = pData(xdata),
    "CentWave_Param" = cwp,
    "FilterIntensity_Param" = inf)
  } else if(Align & !ReFine){
    list(
      "Raw_data" = data,
      "Peak_data" = xdata,
      "Peak_data_align" = xdata_align,
      "PeakGroups_Param_subset" = pdp_subs,
      "PeakGroups_Param_subset_selected" = sampleGroupsType_sel,
      "PeakGroups_Alignment_Param" = pgp_subs, 
      "Meta_info" = pData(xdata),
      "CentWave_Param" = cwp)
    } else if(!Align & ReFine){
      list(
        "Raw_data" = data,
        "Peak_data" = xdata,
        "Peak_data_refined" = xdata_refined,
        "Meta_info" = pData(xdata),
        "CentWave_Param" = cwp,
        "FilterIntensity_Param" = inf)
    } else{
    list(
      "Raw_data" = data,
      "Peak_data" = xdata,
      "Meta_info" = pData(xdata),
      "CentWave_Param" = cwp)
  }
}

#' XCMS correspondance (part2)
#'
#' @param peakID_obj For groupChromPeaks. XCMSnExp object containing the results from a previous peak detection analysis (see findChromPeaks()).
#' @param sampleGroups For PeakDensityParam. A vector of the same length than samples defining the sample group assignments (i.e. which samples belong to which sample group). This parameter is mandatory for the PeakDensityParam and has to be provided also if there is no sample grouping in the experiment (in which case all samples should be assigned to the same group).
#' @param minFraction For PeakDensityParam. numeric(1) defining the minimum fraction of samples in at least one sample group in which the peaks have to be present to be considered as a peak group (feature). 
#' @param bw For PeakDensityParam. numeric(1) defining the bandwidth (standard deviation ot the smoothing kernel) to be used. This argument is passed to the density() method.
#' @param binSize For PeakDensityParam. numeric(1) defining the size of the overlapping slices in mz dimension.
#' @param minSamples For PeakDensityParam. numeric(1) with the minimum number of samples in at least one sample group in which the peaks have to be detected to be considered a peak group (feature).
#' @param maxFeatures For PeakDensityParam. numeric(1) with the maximum number of peak groups to be identified in a single mz slice.
#' @param value For quantify. character specifying the name of the column in chromPeaks(object) that should be returned. Defaults to "into" in which case the integrated peak area is returned. To get the index of the peak in the chromPeaks(object) matrix use "index".
#' @param PeakFill For fillChromPeaks. logical(1) specifying whether values for filled-in peaks should be returned or not. If filled = FALSE, an NA is returned in the matrix for the respective peak. See fillChromPeaks for details on peak filling.
#' @param nCore Number of cores for parallelisation. Numeric(1), default is 1.   
#' 
#' @note The user is able to set most xcms function parameters in this script. However, all xcms function parameters might not be included (i.e., is currently using the default settings). Check the manual (https://bioconductor.org/packages/release/bioc/manuals/xcms/man/xcms.pdf) and adjust this script in case you want to modify additional xcms parameters. 
#' @note To extract info from Peak_overview object, use 'SummarizedExperiment::colData(res_cor$Peak_overview)' or 'SummarizedExperiment::rowData(res_cor$Peak_overview)'.
#' @examples https://gitlab.com/stanstrup-teaching/XCMS-course/-/blob/master/1-XCMS.Rmd
#' @examples https://bioconductor.org/packages/devel/bioc/vignettes/xcms/inst/doc/xcms.html
#' 
XCMS_correspondance <- function(peakID_obj, sampleGroups = numeric(), minFraction = 0.5, bw = 30, binSize = 0.25, minSamples = 1, maxFeatures = 50, value = "into", PeakFill = TRUE, nCore = 1){
  
  # Load packages 
  suppressPackageStartupMessages({
    require(xcms)
    require(SummarizedExperiment)
    require(BiocParallel)
  })

  # Check if sampleGroups exists
  if(length(sampleGroups) == 0){
    stop("The variable 'sampleGroups' is missing")
  }
  
  # Set group parameters for correspondence
  pdp <- PeakDensityParam(sampleGroups = sampleGroups,
                          minFraction = minFraction, 
                          bw = bw,
                          binSize = binSize,
                          minSamples = minSamples,
                          maxFeatures = maxFeatures) 
  
  # Perform the correspondence
  message("Correspondence")
  peakID_obj_gr <- groupChromPeaks(peakID_obj,
                             param = pdp)
  res <- quantify(peakID_obj_gr, 
                  value = value)
  
  # Fill peaks 
 if(PeakFill){
    message("Fill chromatogram peaks")
    setMSnbaseFastLoad(FALSE) 
    peakID_obj_gr_filled <- fillChromPeaks(peakID_obj_gr, 
                                     param = ChromPeakAreaParam(),
                                     BPPARAM = BiocParallel::MulticoreParam(workers = nCore))
    
    # Add feature value matrix with the filled-in data for missing peaks to res
    assays(res)$raw_filled <- featureValues(peakID_obj_gr_filled, 
                                            filled = TRUE) 
    } 
  
  # Store results 
  if(PeakFill){
    list(
      "PeakID_obj" = peakID_obj,
      "PeakID_obj_gr" = peakID_obj_gr,
      "PeakID_obj_gr_filled" = peakID_obj_gr_filled,
      "Peak_overview" = res,
      "Meta_info" = pData(peakID_obj),
      "PeakDensity_Param" = pdp,
      "Filled_peak" = PeakFill,
      "Quant_val" = value
    )
  } else{
  list(
    "PeakID_obj" = peakID_obj,
    "PeakID_obj_gr" = peakID_obj_gr,
    "Peak_overview" = res,
    "Meta_info" = pData(peakID_obj),
    "PeakDensity_Param" = pdp,
    "Filled_peak" = PeakFill,
    "Quant_val" = value
  )
  }
}


#' XCMS peak plotting 
#' 
#' @param XCMS_obj A XCMSnExp object with identified chromatographic peaks (see findChromPeaks()).
#' @param trait_df Dataframe with peak into (peak ID, mz, mzmin, mzmax, rt, rtmin, rtmax, etc.) of peaks.
#' @param title For plotting. Default title = NULL. 
#' @param nCore Number of cores for parallelisation. Numeric(1), default is 1.   
#' 
#' @examples https://gitlab.com/stanstrup-teaching/XCMS-course/-/blob/master/1-XCMS.Rmd
#' @examples https://bioconductor.org/packages/devel/bioc/vignettes/xcms/inst/doc/xcms.html
#'
peak_plot <- function(XCMS_obj, trait_df, title = NULL, nCore = 1){
  
  # Load packaged
  suppressPackageStartupMessages({
    require(xcms)
    require(tidyverse)
    require(BiocParallel)
  })
  
  # Extract mz and rt 
  mzr <- c(trait_df$mzmin %>% min(), trait_df$mzmax %>% max())
  rtr <- c(trait_df$rtmin %>% min(), trait_df$rtmax %>% max())
  
  # Create the feature plot 
  my_plot <- XCMS_obj %>% filterFile(trait_df$sample) %>% 
    filterMz(mzr) %>% 
    filterRt(rtr) %>% 
    xcms::chromatogram(BPPARAM = BiocParallel::MulticoreParam(workers = nCore)) %>% 
    plot(., main = title) 
  return(my_plot)
}


#' XCMS ChromPeakDensity plotting
#' 
#' @param XCMS_obj_gr A XCMSnExp object with identified chromatographic peaks.
#' @param trait_df Dataframe with peak into (peak ID, mz, mzmin, mzmax, rt, rtmin, rtmax, etc.) of peaks.
#' @param sample_colors For plotting. A character vector with names. 
#' @param PeakDensity_Param PeakDensity Parameters from XCMS_correspondance run.
#' @param title For plotting. Default title = NULL. 
#' @param nCore Number of cores for parallelisation. Numeric(1), default is 1.   
#' 
#' @examples https://gitlab.com/stanstrup-teaching/XCMS-course/-/blob/master/1-XCMS.Rmd
#' @examples https://bioconductor.org/packages/devel/bioc/vignettes/xcms/inst/doc/xcms.html
#'
peak_ChromPeakDensity_plot <- function(XCMS_obj_gr, trait_df, sample_colors, PeakDensity_Param, title = NULL, nCore = 1){
  
  # Load packages 
  suppressPackageStartupMessages({
    require(xcms)
    require(magrittr)
    require(BiocParallel)
  })
  
  # Extract mz and rt 
  mzr <- c(trait_df$mzmin %>% min(), trait_df$mzmax %>% max())
  rtr <- c(trait_df$rtmin %>% min(), trait_df$rtmax %>% max())
  
  # Create the feature plot 
  chr_cpc <- xcms::chromatogram(XCMS_obj_gr, 
                                mz = mzr, 
                                rt = rtr, 
                                BPPARAM = BiocParallel::MulticoreParam(workers = nCore))
  plotChromPeakDensity(chr_cpc, 
                       col = sample_colors, 
                       peakBg = sample_colors[chromPeaks(chr_cpc)[, "sample"]],
                       peakCol = sample_colors[chromPeaks(chr_cpc)[, "sample"]],
                       param = PeakDensity_Param,
                       peakPch = 16,
                       main = title)
  }

