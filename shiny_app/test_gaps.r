# Install packages
pacman::p_load("DECIPHER", "plot3D", "plot3Drgl", "ade4", "muscle", "data.table",
               "phangorn", "cowplot", "RColorBrewer", "phylobase", "treeio", "ggtree", 
               "Biostrings", "readxl", "tidyverse", "gtools", "rentrez")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the dataset
sqs <- readAAStringSet("output/20192305_1693_small_substrate_grp_for_rf_duplicates_removed.fasta")

# Extract 34 aa using muscle
source("src/extract_34_aa.r")
query_fils <- sapply(1:length(sqs), function(x) {tempfile(pattern = "", fileext = ".fasta")})
sapply(1:length(sqs), function(x) {writeXStringSet(sqs[x], query_fils[x])})
extract_34_list <- lapply(1:length(sqs), function(x) { extract_34_aa(query_fils[x]) })
extract_34_df <- data.frame(matrix(unlist(extract_34_list), nrow = length(extract_34_list), byrow=T), 
                            stringsAsFactors=FALSE)

# Read in the sequences 
musc_aa <- readAAStringSet("output/20192305_MUSCLE_34aa_extracted.faa")
length(grep("-", musc_aa)) # 752 seqs with at least one gap
table(duplicated(musc_aa))

# Calculate the number of sequences with gaps
hmm_aa <- readAAStringSet("~/Documents/Wageningen_UR/github/amplicon_pred/20192305_1693_small_substrate_grp_for_rf_34extracted_a_dom_hmm.faa")
length(grep("-", hmm_aa)) # 251 seqs with at least one gap
# Compare to the number of sequences with gaps from HMMAlign with NRPS A-domain specific hmm
table(duplicated(hmm_aa))
