install.packages("ape")

library(ape)
seq_length <- function(fasta_string_path){
  
  sequence <- read.dna(fasta_string_path, "fasta", as.character = TRUE, skip = 0)
  number_of_base
  
}


sequence <- read.dna("C:/Users/nicod/Documents/VisualNicolas/RStudioPhylo/NC_004456.txt", "fasta", as.character = TRUE, skip = 0)