#install.packages("ape")
SARS_2Example_function <- function(){
  
getwd()
setwd("C:/Users/nicod/OneDrive - Haute Ecole en Hainaut/Documents/VisualStudio/R_scripts")
library("ape")
acc_numbers <- scan("AccessionNumbers.txt", what = "character", quiet = FALSE)

?scan
acc_numbers_sorted <- sort(acc_numbers)
#str(acc_numbers)

sequences_dna <- read.GenBank(acc_numbers_sorted)
#str(sequences_dna)

sequences_char <- as.character(sequences_dna)

#str(sequences_char)

is_clean_sequence <- function(seq_vector){
    # Collapse the vector of letters into a single string
    seq_string <- paste(seq_vector, collapse = "")
    return(grepl("^[atcgATCG]+$", seq_string))
}
#str(is_clean_sequence)
clean_indices <- sapply(sequences_char, is_clean_sequence)
#str(clean_indices)
#table(clean_indices)
?sapply
sequences_clean <- sequences_dna[clean_indices]
clean_seq_strings <- sapply(as.character(sequences_clean),function(x) paste(x, collapse = ""))
distinct_indices <- !duplicated(clean_seq_strings)
#str(distinct_indices)
final_sequences <- sequences_clean[distinct_indices]
#str(final_sequences)
cat("Number of valid, unique SARS-CoV2 sequences obtained:", length(final_sequences), "\n")
length(final_sequences[1])
# Calculate the length of each sequence
seq_lengths <- mean(sapply(final_sequences, length))

# View the lengths
print(seq_lengths)

return(final_sequences)
}
