library(tidyverse); library(readxl); library(ggsci); library(rgbif)

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
  select(species=`Accepted name`, synonym = Synonyms)
# %>% separate_longer_delim(synonym, delim = " / ") # expand the synonyms to long format

# Import our table
ah_list = read_excel(paste0(dig_freez, "Aquatic Hyphomycetes list/Taxa_AH_Chapter 1_PhD.xlsx"), sheet = "sh_ids") %>% 
  select (species)


# Objective 1. Get all the UNITE SHs related to our species. 
# join Duarte's list with Unite. I want SH, SH type and Synonyms
cross_duarte <- duarte_list %>% left_join(unite, by = "species")

duarte_syns <- duarte_list %>% separate_longer_delim(synonym, delim = " / ") %>% 
  filter(!is.na(synonym), synonym %in% unite$species) %>% rename (species_d = species) %>% 
  left_join(unite,by = c("synonym"="species")) %>% filter (!is.na(sh)) %>% mutate (notes = "annotated synonyms") %>% 
  rename(species=synonym, synonym = species_d)

export <- cross_duarte %>% bind_rows(duarte_syns) %>% 
  arrange (species)

write_csv(export, paste0 (dig_freez, "Aquatic Hyphomycetes list/duarte_unite_annotated.csv"))


# Objective 2. 
all_sh_gf = global_fun$SH %>% strsplit(split = ";") %>% unlist

unite %>% mutate (in_gf = sh %in% all_sh_gf) %>% 
  ggplot(aes(x=in_gf, fill = sh_type))+
  geom_bar(position = "fill")+
  theme_bw(base_size = 16)+xlab(NULL)+
  coord_flip()+scale_fill_d3()+ggtitle("AH SHs in GlobalFungi")

unite %>% mutate (in_gf = sh %in% all_sh_gf) %>% filter (in_gf) %>% 
  ggplot(aes(x=in_gf, fill = sh_type))+
  geom_bar(position = "dodge")+
  theme_bw()+xlab(NULL)+
  coord_flip()+scale_fill_d3()+ggtitle("AH SHs in GlobalFungi")


###


silva_lsu = read_tsv("C:/Data/digital-freezer/databases/SILVA/silva_lsu_species.txt", col_names = c("ids"), col_types = cols()) %>% 
  separate(ids, into =  c("Accession", "Taxonomy"), sep = " ",extra = "merge")

get_spp_silva = function(x){
  my_x = tail(unlist (strsplit(x, split = ";")),1)
  return (my_x)
}

silva_spps = silva_lsu %>% rowwise() %>% mutate (Species=get_spp_silva(Taxonomy), Accession = gsub(">","", Accession)) %>% ungroup 

silva_spps %>% filter (Species %in% duarte_list$species)

silva_spps %>% filter(str_detect(Species, str_c(duarte_list$species, collapse="|")))

silva_spps %>% filter (grepl ("Alatospora", Taxonomy))



#Pleosporales;Clavariopsis aquatica


# let's build the table of interest 
taxa_phd_list="C:/Data/digital-freezer/databases/Aquatic Hyphomycetes list/Taxa_AH_Chapter 1_PhD.xlsx"
ah_unite = read_excel(taxa_phd_list, sheet = "sh_ids") %>% 
  mutate(in_GF = species %in% global_fun$Species_GF, in_unite = species %in% unite$species) %>% 
  select (-ref_sh, -other_sh)

mega_merge = ah_unite %>% mutate (type = "AH_set") %>% 
  full_join(unite %>% select(-taxonomy), by = "species") %>% 
  mutate (in_unite = species %in% unite$species)

mega_merge %>% filter (type == "AH_set") %>% left_join()

global_fun %>% separate_longer_delim(SH, delim = ";") %>% 
  mutate (sh_in_unite = SH %in% unite$sh) %>% left_join(unite, by = c("SH"="sh")) %>% 
  rename (sh_GF=SH)


## let's define some lists
# get duarte and expand the synonyms to long format
duarte_list = read_xlsx("C:/Data/digital-freezer/databases/Aquatic Hyphomycetes list/Lista F_Duarte.xlsx", sheet = 1, skip = 2) %>% 
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

