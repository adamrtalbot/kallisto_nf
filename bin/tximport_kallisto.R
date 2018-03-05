#!/usr/bin/env Rscript

output_folder <- commandArgs(trailingOnly = TRUE)

## Check for presence of tximport, if so load, if not install and load
check_bioconductor_packages <- function(packages){
  for(package in packages){
    if(!require(package, character.only = TRUE)){
      source("https://bioconductor.org/biocLite.R")
      biocLite(package)
    }
  }
}

check_bioconductor_packages(c("tximport", "rhdf5"))


library(tximport)

load_kallisto <- function(output_folder){

  # Gather list of abundance files
  kallisto_output_files <- file.path(list.files(pattern = "abundance.h5",
                                                recursive = TRUE,
                                                full.names = TRUE))


  # Set names of kallisto_output_files

  names(kallisto_output_files) <- unlist(lapply(strsplit(x = kallisto_output_files, split = '/'), function(x) x[length(x) - 1]))

  # Extract count data from kallisto output folders

  kallist_counts <- tximport(files = kallisto_output_files,
                             type = "kallisto",
                             txIn = TRUE,
                             txOut = TRUE)

  # Return just pseudocount data

  return(kallist_counts$counts)
}


# Get files and write to csv
write.csv(
  x = load_kallisto(output_folder = output_folder),
  file = "count_data.csv"
)
