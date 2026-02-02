salmon_out<- read.table("salmon.merged.gene_counts_length_scaled.tsv", header = T)

library(readr)
library(tidyverse)
samples <- read_csv("samples.csv")

# Extract sample descriptions
sample_names <- samples$sample_description

# Create a function to parse the sample names
parse_sample_names <- function(sample_name) {
  # Remove the date and replicate prefix (e.g., "20241211_R1_")
  cleaned <- sub("^\\d{8}_R\\d+_", "", sample_name)
  return(cleaned)
}

clean_treatments <- sapply(sample_names, parse_sample_names, USE.NAMES = FALSE)

simple_metadata<-data.frame(sample=paste0("S",samples$requests_sample_sample_id), 
                            treatment=clean_treatments,
                            experiment=ifelse(grepl("DMSO", clean_treatments) | grepl("SB50", clean_treatments) 
                                              , 2,1))

# Function to extract concentration and time from treatment names
extract_concentration_time <- function(treatment) {
  # Extract concentration (ng/ml for Activin, uM for SB50)
  concentration <- case_when(
    grepl("^0ngmlActivin", treatment) ~ "0 ng/ml",
    grepl("^5ngmlActivin", treatment) ~ "5 ng/ml",
    grepl("^10ngmlActivin", treatment) ~ "10 ng/ml",
    grepl("^15ngmlActivin", treatment) ~ "15 ng/ml",
    grepl("^50uMSB50", treatment) ~ "15 ng/ml",  # Default to 15 ng/ml for SB50
    grepl("DMSO", treatment) ~ "0 ng/ml",  # DMSO is control, so 0 ng/ml
    TRUE ~ NA_character_
  )
  
  # Extract exposure time (in minutes)
  time <- str_extract(treatment, "\\d+min$")
  time <- str_remove(time, "min")
  
  return(data.frame(concentration = concentration, 
                    exposure_time_min = as.numeric(time)))
}

# Apply the function to each treatment
concentration_time <- do.call(rbind, lapply(simple_metadata$treatment, extract_concentration_time))

# Add the new columns to simple_metadata
simple_metadata <- cbind(simple_metadata, concentration_time)


library(readxl)
nodal_score<- read_excel("docs/nodal-score-genes_complete.xlsx", 
                                         skip = 1)


