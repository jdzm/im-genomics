library (tidyverse)
library (readxl)

# path = "C:/Data/databases/GlobalFungiDB/"
path = "I:/Common/COLLABORATORI/Calore Red/Digital Freezer/GlobalFungiDB/"
setwd(path)

taxa_path <- "I:/Common/COLLABORATORI/Calore Red/Digital Freezer/Aquatic Hyphomycetes list/"
taxa_names <- read_excel(paste0(taxa_path, "Taxa_AH_Chapter 1_PhD.xlsx"), sheet = "synonyms", na = "NA") 

synonyms <- taxa_names %>% filter (Present_GF) %>% 
  mutate (Species_GF = ifelse (Spp_GF=="Obligate synonym 1", `Obligate synonym 1`, 
                               ifelse(Spp_GF=="Obligate synonym 3", `Obligate synonym 3`, 
                                      ifelse(Spp_GF=="Basionym", `Basionym`, Species)))) %>%
  select (Species, Species_GF)


files_list<- list.files(pattern = "*txt") %>% as_tibble() %>% separate (value, into = c("type","Species_GF"), remove = F, sep = "_") %>%
  mutate (Species_GF = gsub (".txt","", Species_GF)) %>% left_join(synonyms, by = "Species_GF") %>% rename(filename=value)

global_fungi <- tibble()
for (spp in unique (files_list$Species_GF)){
  extract <- files_list %>% filter(grepl("sh",filename), Species_GF==spp)
  sh_info <- paste0((read_tsv(extract$filename, col_types = cols()) %>% select(SH) %>% distinct())$SH, collapse=";") 
  
  extract <- files_list %>% filter(grepl("sample",filename), Species_GF==spp)
  sample_info <- read_tsv(extract$filename, col_types = cols()) %>% mutate (Species_GF = spp, SH = sh_info, unite_DB = "v10.0") %>%
    select(Species_GF, SH, unite_DB, 1:18)
  
  global_fungi <- bind_rows(global_fungi, sample_info)
}

db_summary <- global_fungi %>% group_by (Species_GF, sample_type) %>% count() 

## Save summary and 
write_csv(global_fungi, file = paste0(path, "GlobalFungiDB_", Sys.Date(), ".csv"))
write_csv(db_summary, file = paste0(path, "GlobalFungiDB_summary_", Sys.Date(), ".csv"))
