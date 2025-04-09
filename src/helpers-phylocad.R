#### Phylogeny Helpers ####

library(Biostrings)
library(GenomicFeatures)
library(rtracklayer)
library(stringr)

#### 1. Extract 16S sequences from AXEK experiments ####
# this function is cool but it is a bit too specific to the AXEK project. Maybe I 
# generalize it in the future
extract_16S_rRNA_sequences <- function(fasta_files, gff_files) {
  if (length(fasta_files) != length(gff_files)) {
    stop("Lengths of fasta_files and gff_files must be the same")
  }
  require ("Biostrings"); require (GenomicFeatures); require (rtracklayer); require (stringr)
  sequences <- DNAStringSet()
  for (i in seq_along(fasta_files)) {
    sample_id <- str_extract(fasta_files[i], "AXEK_\\d+")
    fasta <- readDNAStringSet(fasta_files[i])
    
    # Extract just the contig name (everything after "lcl|" and before the first space)
    fasta_names_clean <- sub("^lcl\\|(\\S+).*", "\\1", names(fasta))
    
    gff <- import.gff(gff_files[i])
    rRNA_features <- subset(gff, type == "rRNA" & grepl("16S", gff$product))
    
    for (j in seq_along(rRNA_features)) {
      seq_id <- as.character(seqnames(rRNA_features[j]))
      start <- start(rRNA_features[j])
      end <- end(rRNA_features[j])
      strand <- as.character(strand(rRNA_features[j]))
      
      # Find matching FASTA index
      fasta_index <- match(seq_id, fasta_names_clean)
      if (!is.na(fasta_index)) {
        extracted_seq <- DNAStringSet(subseq(fasta[[fasta_index]], start=start, end=end))
        # Reverse complement if the strand is "-"
        if (strand == "-") {
          extracted_seq <- reverseComplement(extracted_seq)
        }
        rrna_id <- paste0(sample_id, " 16S rRNA [", seq_id, ":", start, "-", end, "]")
        names(extracted_seq) <- rrna_id
        sequences <- c(sequences, extracted_seq)
      }
    }
    # writeXStringSet(sequences, "16S_AXEK2_all.fasta")
  }
  return(sequences)
}

#### 2. Generalization of the function by gene name
extract_gene_sequences <- function(fasta_files, gff_files, sample_ids, gene_name = "16S") {
  if (length(fasta_files) != length(gff_files) || length(fasta_files) != length(sample_ids)) {
    stop("Lengths of fasta_files, gff_files, and sample_ids must be the same")
  }
  require ("Biostrings"); require (GenomicFeatures); require (rtracklayer); require (stringr)
  # # Example usage:
  # 
  # # File paths
  # fasta_files <- list.files("C:/Data/cadagno-lake/reference_genomes/", pattern = "*fna$", recursive = T, full.names = T)[-1]
  # gff_files <- list.files("C:/Data/cadagno-lake/reference_genomes/", pattern = "*gff$", recursive = T, full.names = T)
  # sample_ids <- c("cmok dsm169", "cmok luedin GCA", "cmok luedin GCF", "thsyn GCA", "thsyn GCF")  # Custom sample names
  # 
  # # Extract 16S rRNA sequences
  # sequences_16S <- extract_gene_sequences(fasta_files, gff_files, sample_ids, gene_name = "16S")
  sequences <- DNAStringSet()
  
  for (i in seq_along(fasta_files)) {
    fasta <- readDNAStringSet(fasta_files[i])
    
    # Extract contig names from FASTA (e.g., "lcl|contig_1 ..." → "contig_1")
    # fasta_names_clean <- sub("^lcl\\|(\\S+).*", "\\1", names(fasta))
    fasta_names_clean <- sapply(strsplit(names(fasta), " "), `[`, 1) # for 
    
    gff <- import.gff(gff_files[i])
    sample_id <- sample_ids[i]
    
    # Extract gene features matching the specified gene_name
    gene_features <- subset(gff, type == "rRNA" & grepl(gene_name, gff$product, ignore.case = TRUE))
    
    for (j in seq_along(gene_features)) {
      seq_id <- as.character(seqnames(gene_features[j]))
      start <- start(gene_features[j])
      end <- end(gene_features[j])
      strand <- as.character(strand(gene_features[j]))
      
      # Find matching sequence in FASTA
      fasta_index <- match(seq_id, fasta_names_clean)
      if (!is.na(fasta_index)) {
        extracted_seq <- subseq(fasta[[fasta_index]], start=start, end=end)
        # Reverse complement if the strand is "-"
        if (strand == "-") {
          extracted_seq <- reverseComplement(extracted_seq)
        }
        
        # Convert to DNAStringSet and assign a meaningful name
        extracted_seq <- DNAStringSet(extracted_seq)
        seq_label <- paste0(sample_id, "_", seq_id, "_", start, "_", end, "_", gene_name)
        names(extracted_seq) <- seq_label
        
        sequences <- c(sequences, extracted_seq)
      }
    }
  }
  
  return(sequences)
}
