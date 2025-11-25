### Extract 16S sequences from AXEK experiment ###
source ("C:/Data/im-genomics/src/helpers-phylocad.R")

## Extract from the PGAP annotations of the AXEK project
# need to point to an "annotation" style folder 
path <- "C:/Data/cadagno-lake/submission_genomes/pgap-annot/"
outpath <- "C:/Data/cadagno-lake/submission_genomes/sequences_16S/"
fasta_files <- paste0(path, "AXEK_", c(1,3:9), "/annot.fna")
gff_files <- paste0(path, "AXEK_", c(1,3:9), "/annot.gtf")

gff_files <- "C:/Users/juan.diaz/Downloads/chmok_luedin/ncbi_dataset/data/GCF_002958735.1/GCF_002958735.1_genomic.gtf"
fasta_files <-   "C:/Users/juan.diaz/Downloads/chmok_luedin/ncbi_dataset/data/GCF_002958735.1/GCF_002958735.1_ASM295873v1_genomic.fna"
seqname = "GCF_002958735.1"

# rework of extraction function to be done
require ("Biostrings"); require (GenomicFeatures); require (rtracklayer); require (stringr)
sequences <- DNAStringSet()

i=1
sample_id <- seqname
fasta <- readDNAStringSet(fasta_files[i])

# Extract just the contig name (everything after "lcl|" and before the first space)
fasta_names_clean <- sub("^lcl\\|(\\S+).*", "\\1", names(fasta)) 
extractor = function(str){return (unlist(strsplit(str, split = " "))[1])}
fasta_names_clean <- sapply(fasta_names_clean, extractor) %>% unname()


gff <- import.gff(gff_files[i])
rRNA_features <- subset(gff, gbkey == "rRNA" & grepl("16S", gff$product))

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

writeXStringSet(unique(sequences), paste0(outpath,seqname,".fasta"))
