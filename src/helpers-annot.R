extract_gff_attribute <- function(gtf_attributes, att_of_interest){
  # takes in the 9th column from a GTF with the attributes and resturns a vector 
  # with the selected attibute. The attribute name must match exactly.
  att_flat <- unlist(gtf_attributes)
  att_of_interest <- paste0(att_of_interest,"=")
  holder = vector()
  for (i in 1:length(att_flat)){
    if(grepl(att_of_interest, att_flat[i])){
      dirty = unlist(strsplit(att_flat[i],";"))[grepl(att_of_interest,unlist(strsplit(att_flat[i],";")))]
      holder = c(holder,gsub(att_of_interest,"",dirty))
    }else{
      holder = c(holder,NA)
    }
  }
  return(holder)
}

# # Example usage to load a gff resulting from PGAP as tibble and extract gene IDs
# # and fields related to gene ontology
# library (tidyverse)
# gff_pgap <- read_tsv("annot.gff", col_names=F, comment = "#", col_types = cols()) %>%
#   filter (X2!="Local", X3!="gene")
# 
# 
# att_list = c("ID","gene","Name","go_function","go_process","go_component","Ontology_term",
#              "inference","product")
# # it takes a couple of seconds to run this over a 7.000 features GFF. Think twice 
# # before extracting a lot of fields over a long dataframe
# for (att in att_list){
#   gff_pgap$newcol <- extract_gff_attribute(gff_pgap$X9, att)
#   names(gff_pgap)[names(gff_pgap) == "newcol"] <- att
# } 

## read in gff file from Gruenenberger Scripts ====
read_in_gff <- function(input_file){
  ape::read.gff(input_file) %>%
    filter(!type %in% c("exon", "gene", "region", "origin of replication")) %>%
    as_tibble() %>%
    mutate(start_feature = start, end_feature = end,strand_feature = strand) %>%
    mutate(Parent = str_split_fixed(str_split_fixed(attributes, ";Parent=",2)[,2],";Dbxref",2)[,1],
                  ecogene = str_split_fixed(str_split_fixed(attributes, ",GeneID", 2)[,1], "EcoGene:",2)[,2],
                  short_gene = str_split_fixed(str_split_fixed(attributes, ";locus_tag", 2)[,1], "gene=",2)[,2],
                  id_name = ifelse(type %in% "repeat_region", str_split_fixed(str_split_fixed(attributes, ";Note=", 2)[,1], "ID=", 2)[,2],
                                   ifelse(type %in% "pseudogene", str_split_fixed(str_split_fixed(attributes, ";Dbxref=", 2)[,1], "ID=", 2)[,2],
                                          ifelse(type %in% "sequence_feature", str_split_fixed(str_split_fixed(attributes, ";Dbxref=", 2)[,1], "ID=", 2)[,2],
                                                 ifelse(type %in% "mobile_genetic_element", str_split_fixed(str_split_fixed(attributes, ";gbkey=", 2)[,1], "ID=", 2)[,2],
                                                        str_split_fixed(str_split_fixed(attributes, ";Parent=", 2)[,1], "ID=", 2)[,2])))),
                  locus_name = ifelse(type %in% c("CDS","mobile_genetic_element", "ncRNA", "recombination_feature"), str_split_fixed(str_split_fixed(attributes, ";product=", 2)[,2], ";", 2)[,1],
                                      ifelse(type ==  "pseudogene",  str_split_fixed(str_split_fixed(attributes, ";gene_biotype", 2)[,1], "gene=", 2)[,2],
                                             ifelse(type == "repeat_region", str_split_fixed(str_split_fixed(attributes, ";gbkey", 2)[,1], "Note=", 2)[,2],
                                                    ifelse(type %in% "sequence_feature", str_split_fixed(str_split_fixed(attributes, ";locus_tag=", 2)[,1], "gene=", 2)[,2],
                                                           ifelse(type %in% "mobile_genetic_element", str_split_fixed(attributes, "insertion sequence:", 2)[,2],
                                                                  ifelse(type == "rRNA", str_split_fixed(str_split_fixed(attributes, ";product=", 2)[,2], " ", 2)[,1], 
                                                                         ifelse(type == "tRNA", str_split_fixed(attributes, ";product=", 2)[,2], NA ))))))),
                  width = abs(start_feature - end_feature)) %>%
    select(seqid, id_name, locus_name, start_feature, end_feature, strand_feature, Parent, type, width, ecogene, short_gene) %>%
    mutate(gene = str_split_fixed(Parent,"-",2)[,2])
}

## list to ;-separated values (useful for GO)
list_to_char <- function (x){return (paste0(unlist(x), collapse = ";"))}

