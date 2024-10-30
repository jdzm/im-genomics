library(tidyverse); library(readxl); library(ggsci)
library(rgbif)

dig_freez = "C:/Data/digital-freezer/databases/" # for local work
dig_freez = "I:/Common/COLLABORATORI/Calore Red/Digital Freezer/" # for remote work

# import all UNITE sequence names from the V10 2024 04 release.
unite = read_delim(paste0(dig_freez, "UNITE/sh_headers_dynamic.txt"),delim = "|", col_types = cols(),
                 col_names = c("species","accession","sh","sh_type","taxonomy")) %>% 
  mutate(species = gsub(">", "",species), species = gsub("_", " ", species, fixed = T))

# Import the global Fungi list
global_fun = read_csv(paste0(dig_freez, "GlobalFungiDB/GlobalFungiDB_latest.csv"), col_types = cols()) %>% 
  select (Species_GF, SH) %>% distinct()

# Import Duarte's list
duarte_list = read_xlsx(paste0(dig_freez, "Aquatic Hyphomycetes list/Lista F_Duarte.xlsx"), sheet = 1, skip = 2) %>% 
  select(species=`Accepted name`, synonym = Synonyms)%>% 
  separate_longer_delim(synonym, delim = " / ") # expand the synonyms to long format

# Import our table
ah_list = read_excel(paste0(dig_freez, "Aquatic Hyphomycetes list/Taxa_AH_Chapter 1_PhD.xlsx"), sheet = "sh_ids") %>% 
  select (species)

# Objective 1. Get all the UNITE SHs related to our species. 
# Label those with no SH for further investigation. We are possibly missing some 
# due to not checking synonyms here. Keep in mind that some spps have multiple SHs
ah_db <- ah_list %>% mutate(in_unite = species %in% unite$species) %>% 
  left_join(unite, by = "species") 


# Articulospora genus has a lot of SHs that are Articulospora sp. Not sure what to do 
# with this. # here is a checker
ah_db %>% group_by(species) %>% select (-accession, -sh_type) %>% 
  mutate (all_shs = str_c(sh, collapse = ";"), n = n()) %>% 
  arrange(desc(n)) %>% ungroup %>% select(-sh) %>% distinct() %>% select (-taxonomy)

ah_db_collapsed <- ah_db %>% group_by(species,in_unite) %>% select (-accession, -sh_type) %>% 
  summarise (all_shs = str_c(sh, collapse = ";")) %>% ungroup

# Objective 2. 
all_sh_gf = global_fun$SH %>% strsplit(split = ";") %>% unlist

unite %>% mutate (in_gf = sh %in% all_sh_gf) %>% 
  ggplot(aes(x=in_gf, fill = sh_type))+
  geom_bar(position = "fill")+
  theme_bw()+xlab(NULL)+
  coord_flip()+scale_fill_d3()+ggtitle("AH SHs in GlobalFungi")

unite %>% mutate (in_gf = sh %in% all_sh_gf) %>% filter (in_gf) %>% 
  ggplot(aes(x=in_gf, fill = sh_type))+
  geom_bar(position = "dodge")+
  theme_bw()+xlab(NULL)+
  coord_flip()+scale_fill_d3()+ggtitle("AH SHs in GlobalFungi")

# let's build the table of interest 
taxa_phd_list ="C:/Data/digital-freezer/Aquatic Hyphomycetes list/Taxa_AH_Chapter 1_PhD.xlsx"
ah_unite = read_excel(taxa_phd_list, sheet = "sh_ids") %>% 
  mutate(in_GF = species %in% global_fun$Species_GF) %>% select (-ref_sh, -other_sh)

mega_merge = ah_unite %>% mutate (type = "AH_set") %>% 
  full_join(unite %>% select(-taxonomy), by = "species") %>% 
  mutate (in_unite = species %in% unite$species)

mega_merge %>% filter (type == "AH_set") %>% left_join()

global_fun %>% separate_longer_delim(SH, delim = ";") %>% 
  mutate (sh_in_unite = SH %in% unite$sh) %>% left_join(unite, by = c("SH"="sh")) %>% 
  rename (sh_GF=SH)



## let's define some lists
# get duarte and expand the synonyms to long format
duarte_list = read_xlsx("C:/Data/digital-freezer/Aquatic Hyphomycetes list/Lista F_Duarte.xlsx", sheet = 1, skip = 2) %>% 
  select(species=`Accepted name`, synonym = Synonyms)%>% 
  separate_longer_delim(synonym, delim = " / ")

## cross Duarte and Calore
cross_duarte = duarte_list %>% 
  mutate(in_AH = species %in% unique(ah_unite$species), 
         syn_in_AH = synonym %in% unique(ah_unite$species))

ah_unite



###########
# get data from GBIF
# Download separately the ah list and the duarte list. Then save backbone checks for name tracing
# # To set credentials
# usethis::edit_r_environ()
library (rgbif)

setwd(paste0(dig_freez, "GBIF"))
# First Red's list 
gbif_backbone_check <- ah_list %>% 
  name_backbone_checklist()  # match to backbone 
  
gbif_taxon_keys <- gbif_backbone_check %>% filter(!matchType == "NONE") %>% # get matched names
  pull(usageKey) 

occ_download(
  pred_in("taxonKey", gbif_taxon_keys), # important to use pred_in
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  format = "SIMPLE_CSV"
)

runkey='0031875-241007104925546' ## Copypasted from the occ_download message
occ_download_wait(runkey)
d <- occ_download_get(runkey) %>%
  occ_download_import()

### and same procedure for Duarte's list
gbif_backbone_check <- duarte_list %>% select(species) %>% distinct() %>% 
  name_backbone_checklist()  # match to backbone 

gbif_taxon_keys <- gbif_backbone_check %>% filter(!matchType == "NONE") %>% # get matched names
  pull(usageKey) ## This aso matcheshio

occ_download(
  pred_in("taxonKey", gbif_taxon_keys), # important to use pred_in
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  format = "SIMPLE_CSV"
)

runkey='0032506-241007104925546'
occ_download_wait(runkey)
d <- occ_download_get(runkey) %>%
  occ_download_import()

##### Downloads completed. Saved with weird names that I change manually. Also keep in mind
# that it saves everything as tsv but it names it csv.
# let's save the backbone key matching of what I checked. Just for traceability

gbif_downloads = duarte_list %>% select(species) %>% 
  left_join(ah_list, by = 'species') %>% distinct() 

gbif_backbone_check <- gbif_downloads %>%  
  name_backbone_checklist() %>% 
  mutate (in_duarte = verbatim_name %in% duarte_list$species, in_ah = verbatim_name %in% ah_list$species) 
# write_tsv(gbif_backbone_check, file = paste0(dig_freez,"GBIF/gbif-backbone-match.tsv"))

################

