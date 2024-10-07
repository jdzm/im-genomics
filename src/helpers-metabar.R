##### Collection of functions to assist metabarcoding analysis with R ######
#### 1. Load taxonomy from Nanopore data ####
load_tax_nanopore <- function(path, silva = TRUE) {
  # Define the columns to split based on whether SILVA taxonomy is used
  # If using SILVA taxonomy, process taxonomy to the 'Genus' rank (without 'Species')
  tax_cols <- if (silva) {
    c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
  } else {
    c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  }

  # Separate the 'tax' column into the defined taxonomic ranks
  abundances <- read_tsv(path, col_types = cols()) %>%
    separate(tax, into = tax_cols, sep = ";", extra = "merge")

  # Return the data with taxonomic ranks as separate columns
  return(abundances)
}

#### 2. Plot taxonomy data
plot_taxonomy <- function(df, taxa_level = "Genus", n_taxa = 15) {
    require(tidyverse)
    library(RColorBrewer)  # For color palettes
    # Validate taxonomic level and required columns
    valid_taxa_levels <- c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    if (!taxa_level %in% valid_taxa_levels) stop("Invalid Taxonomic category")
    if (!all(c("count", "sample") %in% colnames(df))) stop("Input df needs to contain the columns 'sample' and 'count'")
    
    # Subset and rename columns
    subset_df <- df %>% select(all_of(taxa_level), count, sample)
    colnames(subset_df) <- c("taxa", "count", "sample")
    
    # Identify the top taxa and collapse the rest into "Other"
    top_taxa <- subset_df %>%
        group_by(taxa) %>%
        summarise(count = sum(count), .groups = 'drop') %>%
        arrange(desc(count)) %>%
        slice_head(n = n_taxa) %>%
        pull(taxa)
    
    mylevels <- c(top_taxa, "Other")
    mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(n_taxa + 1)
    
    # Create a summarized composition of taxa
    composition <- subset_df %>%
        mutate(taxa = ifelse(taxa %in% top_taxa, taxa, "Other")) %>%
        group_by(taxa, sample) %>%
        summarise(count = sum(count), .groups = 'drop') %>%
        mutate(taxa = factor(taxa, levels = mylevels))
    
    # Plot the taxonomy composition
    taxaplot <- composition %>%
        ggplot(aes(x = sample, y = count, fill = taxa)) +
        geom_col(position = 'fill') +
        scale_fill_manual(values = mycolors, name = taxa_level) +
        theme_bw() + 
        xlab(NULL) + 
        ylab("Relative abundance")
    
    return(taxaplot)
}

#### 3. ONT abundance to microeco
get_tax_table = function (abundance, digits = 4, silva = T){
  # extract taxonomy from raw abundance table and add unique names
  tax_cols <- if (silva) {
    c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
  } else {
    c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  }
  taxtable <- abundance %>%
    select(all_of(tax_cols)) %>% distinct() %>%
    mutate(taxname = paste0("TAX", sprintf(paste0("%0", digits, "d"), row_number()))) %>%
    select (taxname, all_of(tax_cols))
  return (taxtable)
}

