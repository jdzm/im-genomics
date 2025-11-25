library(tidyverse)
library(OrthoDB) #remotes::install_gitlab('ezlab/orthodb_r',build_vignettes = TRUE, dependencies = TRUE)

sample = "AXEK_9"
busco_path = "C:/Data/cadagno-lake/submission_genomes/busco"
busco_list_file <- list.files(busco_path, 
                              pattern = "missing_busco_list.tsv", 
                              recursive = TRUE, full.names = TRUE) 
chosen_one = grep(sample, busco_list_file, value=T)
busco_set = str_extract(chosen_one, "(?<=busco/)[^/]+")
busco_version <- str_extract(busco_set, "\\d+$") %>% paste0("v", .)

busco_list = read.table(chosen_one)[,1]
description_df = tibble()
for (busco_id in busco_list){
  api <- OrthoDB::OdbAPI$new(busco_version)
  OGS <- api$orthologs(busco_id)
  
  holder = tibble()
  for (i in seq_along(OGS$df$data$genes)){
    holder = bind_rows(holder, OGS$df$data$genes[[i]] %>% mutate(id = i))
  }
  descriptions = holder$description %>% unique %>% paste(collapse = ",")
  description_df = bind_rows(description_df, tibble("version"=busco_set,"busco_id" = busco_id, "description" = descriptions))
}

description_df
