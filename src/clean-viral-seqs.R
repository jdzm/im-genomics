#### Part 0. Set up the stage for a quick consensus sequence and tree adventure. ####
# Load libraries and set seed
my_packages <- c ("tidyverse", "Biostrings", "DECIPHER", "phangorn", "vegan", "msa")
msg <- paste0(suppressPackageStartupMessages(lapply (my_packages, require, character.only = T)))
set.seed(1225)

# Define input/output dirs
seqpath <- "I:/Common/COLLABORATORI/_ALLIEVI LABORATORISTI/2024/Giovannini Sofia/studio Emys/sequenze Adeno/sequenze corrette/"
outpath <- paste0(seqpath, "results/") # specify any outpath that you like. Make sure it ends with "/"
setwd(seqpath)
dir.create(outpath, recursive = T, showWarnings = F)

# list all files and extract sample IDs
my_seqs = list.files(seqpath, pattern = "Adeno") %>% as_tibble() %>% 
  separate (value, into = c("Adeno", "id", "sense", "polType"), sep = "_", remove = F) %>% select (-polType,-Adeno) %>% 
  mutate (id = as.integer(id)) %>% arrange(id) 


### Part 1. Get consensus sequences ####
# initialize vectors for saving the sequences. We will be using the dnastrings type 
# of object mostly because it prints the letters colored by nucleotide and I like it
ids = unique(my_seqs$id)
all_forward=DNAStringSet()
all_reverse=DNAStringSet()
all_consensus=DNAStringSet()
all_sequences=DNAStringSet()
all_alignments=c() # this one is a list. it will contain nested DNAstrings alignment objects for visual review
for (sid in ids){
  # Read forward and reverse sequences from FASTA files
  forward_seq <- readDNAStringSet(paste0(seqpath, "Adeno_",sid,"_F_polF inner AdenoC.fas"))
  reverse_seq <- readDNAStringSet(paste0(seqpath, "Adeno_",sid,"_R_polR inner AdenoC.fas"))
  
  # Reverse complement the reverse sequences
  reverse_seq_rc <- reverseComplement(reverse_seq)
  sequences <- c(forward_seq, reverse_seq_rc)
  
  # Align the forward sequence to the rc of the reverse sequence and output the consensus
  aligned_sequences <- AlignSeqs(sequences)
  consensus_seq <- ConsensusSequence(aligned_sequences)
  names(consensus_seq) <- paste0("Consensus_Adeno_",sid)

  # Append each sequence to their corresponding vector 
  all_sequences=c(all_sequences, sequences)
  all_alignments=c(all_alignments, aligned_sequences)
  all_consensus=c(all_consensus, consensus_seq)
  all_forward=c(all_forward,forward_seq)
  all_reverse=c(all_reverse,reverse_seq)
}

# Write sequences sequence to FASTA files
writeXStringSet(all_forward, paste0(outpath,"forward_sequences_all.fasta"))
writeXStringSet(all_reverse, paste0(outpath,"reverse_sequences_all.fasta"))
writeXStringSet(all_sequences, paste0(outpath,"sequences_all.fasta"))
writeXStringSet(all_consensus, paste0(outpath,"consensus_sequences_all.fasta"))

# # Sequences are ready for BLASTing 

#### Part 2. Align our sequences using msa package. ####
# I am doing here two alignments. One with all the sequences as we got them and a 
# second one with all of our consensus sequences.
msa_all <- AlignSeqs(all_sequences)
msa_consensus <- AlignSeqs(all_consensus)

# Save the alignments for visualization with JalView or MEGA
writeXStringSet(msa_all, paste0(outpath,"alignment_all.fasta"))
writeXStringSet(msa_consensus, paste0(outpath,"alignment_consensus.fasta"))

#### Part 3. Create an ML tree with Vegan and Phangorn
# # You can also load the alignment directly 
# phy_sequences <- read.dna(paste0(outpath,"alignment_consensus.fasta"), format = "fasta")
# phy_data <- phyDat(phy_sequences, type = "DNA")

# The following lines transform our alignment into the phangorn format, create an initial 
# tree with NJ algorithm. Then we perform a Maximum likelihood (ML) estimation and we optimize it. 
phy_data <- phyDat(as(msa_consensus, "matrix"), type = "DNA")
initial_tree <- nj(dist.ml(phy_data))
ml_tree <- pml(initial_tree, data = phy_data)
optimized_ml_tree <- optim.pml(ml_tree)

# Plot the maximum likelihood tree and save it
pdf(file = paste0(outpath,"consensus_tree.pdf"), width = 6, height = 4)
plot(optimized_ml_tree$tree, main = "ML Tree of Consensus Sequences")
dev.off()

# Save the tree to files
write.tree(optimized_ml_tree$tree, paste0(outpath,"consensus_sequences.tree"))

# Repeat with the alingment of all sequences for completeness
phy_data <- phyDat(as(msa_all, "matrix"), type = "DNA")
initial_tree <- nj(dist.ml(phy_data))
ml_tree <- pml(initial_tree, data = phy_data)
optimized_ml_tree <- optim.pml(ml_tree)

# Plot the maximum likelihood tree and save it
pdf(file = paste0(outpath,"sequences_all_tree.pdf"), width = 6, height = 4)
plot(optimized_ml_tree$tree, main = "ML Tree of Corrected Sequences")
dev.off()

# Save the tree to files
write.tree(optimized_ml_tree$tree, paste0(outpath,"sequences_all.tree"))

