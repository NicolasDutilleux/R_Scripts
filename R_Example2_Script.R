install.packages("ape")

library(ape)
seq_length <- function(fasta_file){
  
  sequence <- read.dna(fasta_file, "fasta", as.character = TRUE, skip = 0)
  number_of_base <- length(sequence)
}
sequence_length <- seq_length("C:/Users/nicod/Documents/VisualNicolas/RStudioPhylo/NC_004456.fasta")

?read.dna
