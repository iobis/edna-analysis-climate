library(dplyr)
library(stringr)
library(purrr)
library(arrow)

#' Download occurrence data
download_occurrence <- function(force = FALSE) {
  if (!dir.exists("data/output") || force) {
    cat("Downloading occurrence data ======\n")
    system("rm -r output output.zip")
    download.file(
        "https://obis-edna-results.s3.amazonaws.com/output.zip",
        "output.zip"
    )
    unzip("output.zip", exdir = "data/")
    file.remove("output.zip")
  } else {
    cat("Occurrence data already downloaded ======\n")
  }
}

#' Reads the full occurrence dataset from TSV.
read_occurrence_tsv <- function() {
  occurrence_files <- list.files("data/output", "*Occurrence*", full.names = TRUE)
  
  dna_files <- list.files("data/output", "*DNADerivedData*", full.names = TRUE)
  
  dna <- map(dna_files, read.table, sep = "\t", quote = "", header = TRUE) %>%
    bind_rows() %>%
    mutate_if(is.character, na_if, "")

  occurrence <- map(occurrence_files, read.table, sep = "\t", quote = "", header = TRUE, comment.char = "") %>%
    bind_rows() %>%
    mutate_if(is.character, na_if, "")%>%
    left_join(dna, by = "occurrenceID")

  return(occurrence)
}

#' Converts the occurrence dataset to parquet.
convert_occurrence_data <- function() {
  message("Converting TSV to parquet...")
  occurrence <- read_occurrence_tsv()
  write_parquet(occurrence, "data/output/occurrence.parquet")
  return(occurrence)
}

#' Reads the full occurrence dataset from parquet, converts TSV to parquet if the parquet file does not exist.
read_occurrence_data <- function() {
  if (!file.exists("data/output/occurrence.parquet")) {
    occurrence <- convert_occurrence_data()  
  } else {
    occurrence <- read_parquet("data/output/occurrence.parquet")
  }
  return(occurrence)
}
