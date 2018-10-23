library(easypackages)
libraries('VennDiagram', 'reshape2', 'GenomicRanges', 'ggplot2', 'limma', 'tidyverse', 'ggforce', 'eulerr')

# Filtering out the data based on the -Log10(p-val) and M-val cut off
filter_data <- function(x, p_val, M_cutoff){
  Filter_data_x <- x[x$log10_p_value > p_val | x$log10_p_value == 0, ]
  Filtered_data_x <- Filter_data_x[Filter_data_x$M_value_rescaled < M_cutoff,]
  return(Filtered_data_x)
}

# Comparing genomics ranges and peaks
compare_peaks <- function(sample_1, sample_2, percent = FALSE, counts = FALSE){
  gr0 = unique(GRanges(sample_1$chrom, IRanges(start = sample_1$start, end = sample_1$end)))
  gr1 = unique(GRanges(sample_2$chrom, IRanges(start = sample_2$start, end = sample_2$end)))
  denom <- sum(length(gr0), length(gr1))
  
  #counts
  sample_1_area <- length(gr0) - length(findOverlaps(gr1, gr0))
  sample_2_area <- length(gr1) - length(findOverlaps(gr1, gr0))
  overlap_area <- length(findOverlaps(gr1, gr0))
  
  # percentages
  sample_1_per <- (sample_1_area/denom)*100
  sample_2_per <- (sample_2_area/denom)*100
  overlap_per <- (overlap_area/denom)*100
  
  if (counts == TRUE){
    fit <- euler(c(AF4 = sample_1_area, AF9 = sample_2_area, "AF4&AF9" = overlap_area))
    plot(fit, fills = list(fill = c("white", "steelblue4"), alpha = 0.5), quantities = TRUE)
  }
  else if (percent == TRUE){
    fit_per <- euler(c(AF4 = round(sample_1_per), AF9 = round(sample_2_per), "AF4&AF9" = round(overlap_per)))
    plot(fit_per, fills = list(fill = c("white", "steelblue4"), alpha = 0.5), quantities = TRUE)
  } 
}
