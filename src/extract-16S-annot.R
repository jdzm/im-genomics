### Extract 16S sequences from AXEK experiment ###
source ("C:/Data/im-genomics/src/helpers-phylocad.R")

## Extract from the PGAP annotations of the AXEK project
# need to point to an "annotation" style folder 
path <- "C:/Data/cadagno-lake/submission_genomes/pgap-annot/"
outpath <- "C:/Data/cadagno-lake/submission_genomes/"
fasta_files <- paste0(path, "AXEK_", c(1,3:9), "/annot.fna")
gff_files <- paste0(path, "AXEK_", c(1,3:9), "/annot.gtf")

# Extract sequences
sequences <- extract_16S_rRNA_sequences(fasta_files, gff_files)

# Write to FASTA
writeXStringSet(sequences, paste0(outpath,"16S_sequences_AXEK.fasta"))
writeXStringSet(unique(sequences), paste0(outpath,"16S_sequences_AXEK_unique.fasta"))


### Extract 16S from the AXEK2 separate samples
# Define paths
path <- "C:/Data/cadagno-lake/investigate/annotation/"
outpath <- "C:/Data/cadagno-lake/investigate/"
fasta_files <- paste0(path, "asm_AXEK_2", c("","_contig_1","_contig_2"), "/annot.fna")
gff_files <- paste0(path, "asm_AXEK_2", c("","_contig_1","_contig_2"), "/annot.gtf")
sequences <- extract_16S_rRNA_sequences(fasta_files, gff_files)

# Write to FASTA

writeXStringSet(sequences, paste0(outpath,"16S_AXEK2_all.fasta"))
writeXStringSet(unique(sequences), paste0(outpath, "16S_AXEK2_unique.fasta"))

