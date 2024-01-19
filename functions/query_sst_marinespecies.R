#################### eDNA expeditions - scientific analysis ####################
########################## Environmental data download #########################
# January of 2024
# Author: Silas C. Principe
# Contact: s.principe@unesco.org
#
################### Query SST values/summaries for species #####################


#' Query temperature data from the OBIS/GBIF database
#'
#' @param species scientific name of the species
#' @param species_gbif equivalent scientific name for GBIF. If null, the
#'   \code{species} will be used
#' @param return_summary if \code{TRUE}, then a summary is returned instead of
#'   the full data. Default is \code{FALSE}
#' @param depth depth for which the variables should be retrieved 
#'   (depthsurf, depthmean or depthmax). If \code{NULL}, all are returned
#' @param variant the variant (i.e., min, max or mean) for which to query. If
#'   \code{NULL}, all are returned
#'
#' @return a data frame with the full data or the summaries
#' @export
#'
#' @examples
#' \dontrun {
#' query_sst("Acanthurus chirurgus")
#' }
query_sst <- function(species,
                      species_gbif = NULL,
                      return_summary = FALSE,
                      depth = NULL,
                      variant = NULL){
  library(DBI)
  
  database <- "~/Downloads/protectedseas/database.sqlite"
  con <- dbConnect(RSQLite::SQLite(), database)
  
  gbif_occurrence_table <- "gbif_occurrence"
  obis_occurrence_table <- "obis_occurrence"
  
  env_vars <- c(
    "thetao_baseline_depthmax_max", "thetao_baseline_depthmax_mean",
    "thetao_baseline_depthmax_min", "thetao_baseline_depthmean_max",
    "thetao_baseline_depthmean_mean", "thetao_baseline_depthmean_min",
    "thetao_baseline_depthsurf_max", "thetao_baseline_depthsurf_mean",
    "thetao_baseline_depthsurf_min"
  )
  
  if (!is.null(depth)) {
    if (!depth %in% c("depthsurf", "depthmean", "depthmax")) {
      stop('Depth should be one of "depthsurf", "depthmean", "depthmax"')
    } else {
      env_vars <- env_vars[grepl(depth, env_vars)]
    }
  }
  
  if (!is.null(variant)) {
    if (!variant %in% c("min", "max", "mean")) {
      stop('Variant should be one of "min", "max", "mean"')
    } else {
      env_vars <- env_vars[grepl(variant, env_vars)]
    }
  }
  
  cli::cat_line("Running species ", cli::col_cyan(species))
  
  cli::cli_progress_step("Querying OBIS")
  
  sel_species <- species
  
  obis_query <- glue::glue("
      with filtered_data as (
        select species, h3
        from {obis_occurrence_table}
        where {obis_occurrence_table}.species = '{sel_species}'
        )
      select species, filtered_data.h3, {paste(env_vars, collapse = ', ')}
      from filtered_data
      inner join env_current on filtered_data.h3 = env_current.h3;
                         ")
  
  obis_res <- dbSendQuery(con, obis_query)
  obis_species <- dbFetch(obis_res)
  
  cli::cli_progress_step("Querying GBIF")
  sel_species <- ifelse(!is.null(species_gbif), species_gbif, species)
  
  gbif_query <- glue::glue("
      with filtered_data as (
        select species, h3
        from {gbif_occurrence_table}
        where {gbif_occurrence_table}.species = '{sel_species}'
        )
      select species, filtered_data.h3, {paste(env_vars, collapse = ', ')}
      from filtered_data
      inner join env_current on filtered_data.h3 = env_current.h3;
                         ")
  
  gbif_res <- dbSendQuery(con, gbif_query)
  gbif_species <- dbFetch(gbif_res)
  
  all_result <- bind_rows(obis_species, gbif_species)
  
  if (nrow(all_result) == 0) {
    final_result <- NULL
    warning("No results found for the species")
  } else {
    if (return_summary) {
      q25_fun <- function(x) quantile(x, .25)
      q75_fun <- function(x) quantile(x, .75)
      final_result <- all_result %>%
        select(-h3, -species) %>%
        summarise(across(starts_with("thetao"),
                         list(max = max, min = min,
                              mean = mean, sd = sd, median = median,
                              q25 = q25_fun, q75 = q75_fun),
                         .names = "{.col}${.fn}")) %>%
        pivot_longer(1:length(.), names_to = "variable", values_to = "value") %>%
        separate_wider_delim(cols = "variable", names = c("variable", "variant"),
                             delim = "$") %>%
        pivot_wider(names_from = variant, values_from = value) %>%
        mutate(scientificName = species)
    } else {
      final_result <- all_result
    }
  }
  
  cli::cli_progress_done()
  
  dbDisconnect(con)
  
  return(final_result)
  
}

# Test
#query_sst("Albula glossodonta")

### END