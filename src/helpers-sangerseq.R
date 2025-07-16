# functions for viral sanger sequencing cleanup

#### Function to correct consensus sequences by removing non-ATCG bases that have been identified ####
# in either fwd or rev sequences
correct_consensus <- function(forward, reverse, consensus) {
  requireNamespace("Biostrings", quietly = TRUE)
  
  # Convert to character vectors
  fwd <- strsplit(as.character(forward), "")[[1]]
  rev <- strsplit(as.character(reverse), "")[[1]]
  cons <- strsplit(as.character(consensus), "")[[1]]
  
  # Identify ambiguous bases (not A, T, C, G)
  ambiguous_pos <- which(!cons %in% IUPAC_CODE_MAP[1:4])
  
  # Replace ambiguous bases with unambiguous ones from forward or reverse
  for (i in ambiguous_pos) {
    if (fwd[i] %in% IUPAC_CODE_MAP[1:4]) {
      cons[i] <- fwd[i]
    } else if (rev[i] %in% IUPAC_CODE_MAP[1:4]) {
      cons[i] <- rev[i]
    }
  }
  
  return(DNAStringSet(paste0(cons, collapse = "")))
}

#### Function to get consensus sequence from a pair of sanger sequences ####
get_consensus <- function(forward, reverse){
  
  
}

#### Function to update entries on metadata based on the three obbligatory fields ####
## filename, date, folder
update_meta_turtles <- function (df, metadata){
  # df <- df %>% select(filename, date_of_seq, folder)
  new_rows <- anti_join(df, metadata, by = c("filename", "date_of_seq", "folder"))
  
  updated_metadata <- bind_rows(metadata, new_rows)
  
  return (updated_metadata)
}

#### Asign target species based on primer ####
assign_patho_target <- function(primer) {
  primer <- tolower(primer)
  
  dplyr::case_when(
    primer == "mycseq" ~ "Mycoplasma",
    primer %in% c("tgv", "iyg", "hp3") ~ "Herpesvirus",
    primer %in% c("pol-inner", "pol inner", "polr", "polf") ~ "Adenovirus",
    TRUE ~ NA_character_
  )
}

