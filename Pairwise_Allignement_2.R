#install.packages("seqinr")
install.packages("rtools")
library(seqinr)
# 1. Install the Bioconductor Manager (if you haven't already)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# 2. Use BiocManager to install Biostrings
BiocManager::install("Biostrings")
BiocManager::install("pwalign")
library(pwalign)
library("Biostrings")
?seqinr

string1 = unlist(strsplit("ATTCCGA", split = ""))
string2 = unlist(strsplit("TATCG", split = ""))
#Scoring matrix : 2 for match, -1 for mismatch. Linear gap penalty d = 4

# Note: We added 'i' and 'j' to the arguments so the function knows 
# which position it is currently calculating.
SARS_2Example_function <- function(){
  
  # This opens a window. Select your "AccessionNumbers.txt" file.
  file_path <- file.choose()
  
  # This automatically sets the working directory to where that file is.
  setwd(dirname(file_path))
  
  # Check if it worked
  getwd()
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
recursive_global_alignment <- function(seq1, seq2, i, j) {
  
  # --- 1. BASE CASES ---
  if (i == 0 && j == 0) return(0)       # Origin
  if (i == 0) return(j * -4)            # Top edge (Gap penalty accumulation)
  if (j == 0) return(i * -4)            # Left edge (Gap penalty accumulation)
  
  # --- 2. RECURSIVE STEPS ---
  
  # A. Diagonal Move (Match/Mismatch)
  # Note: seq1[i] works here because we will convert string to vector below
  if (seq1[i] == seq2[j]) {
    val_diag <- recursive_global_alignment(seq1, seq2, i - 1, j - 1) + 2
  } else {
    val_diag <- recursive_global_alignment(seq1, seq2, i - 1, j - 1) - 1
  }
  
  # B. Up Move (Gap in seq2)
  val_up <- recursive_global_alignment(seq1, seq2, i - 1, j) - 4
  
  # C. Left Move (Gap in seq1)
  val_left <- recursive_global_alignment(seq1, seq2, i, j - 1) - 4
  
  # --- 3. RETURN MAX ---
  return(max(val_diag, val_up, val_left))
}

print(recursive_global_alignment(string1, string2, length(string1), length(string2)))

create_score_matrix_global <- function(seq1_vec, seq2_vec, match = 2, mismatch = -1, gap = -4) {
  
  # 1. Initialize the Matrix
  # Rows = seq1 length + 1 (for the gap row)
  # Cols = seq2 length + 1 (for the gap col)
  rows <- length(seq1_vec) + 1
  cols <- length(seq2_vec) + 1
  
  # Create a matrix filled with NA
  score_mat <- matrix(NA, nrow = rows, ncol = cols)
  
  # Optional: Label rows/cols for easier reading
  rownames(score_mat) <- c("Gap", seq1_vec)
  colnames(score_mat) <- c("Gap", seq2_vec)
  
  # 2. Base Cases (Initialize Edges)
  
  # Top-left corner is 0
  score_mat[1, 1] <- 0
  
  # Fill First Column (Vertical / i axis)
  for (i in 2:rows) {
    score_mat[i, 1] <- score_mat[i-1, 1] + gap
  }
  
  # Fill First Row (Horizontal / j axis)
  for (j in 2:cols) {
    score_mat[1, j] <- score_mat[1, j-1] + gap
  }
  
  # 3. Fill the rest of the matrix
  for (i in 2:rows) {
    for (j in 2:cols) {
      
      # A. Calculate Diagonal Score (Match or Mismatch)
      # Note: We use i-1 and j-1 for vectors because the matrix has an extra "Gap" row/col
      if (seq1_vec[i-1] == seq2_vec[j-1]) {
        score_diag <- score_mat[i-1, j-1] + match
      } else {
        score_diag <- score_mat[i-1, j-1] + mismatch
      }
      
      # B. Calculate Up Score (Gap in seq2)
      score_up <- score_mat[i-1, j] + gap
      
      # C. Calculate Left Score (Gap in seq1)
      score_left <- score_mat[i, j-1] + gap
      
      # D. Fill the cell with the MAX score
      score_mat[i, j] <- max(score_diag, score_up, score_left)
    }
  }
  
  return(score_mat)
}

find_best_route_global <- function(my_matrix, gap_penalty, match_score, match_penalty){
  
  
  number_of_match <- 0
  number_of_gaps <- 0
  number_of_mismatch <- 0
  col <- length(my_matrix[1, ])
  row <- length(my_matrix[ ,1])
  results <- 0
  results[1] <- my_matrix[row, col]
  while(col != 1 || row != 1){
    if(col == 1){
      number_of_gaps <- number_of_gaps + 1
      row <- row - 1
    }
    if(row == 1){
      number_of_gaps <- number_of_gaps + 1
      col <- col - 1
    }
    
    
    
    if(my_matrix[row, col] == (my_matrix[row - 1, col] + gap_penalty)){
      number_of_gaps <- number_of_gaps + 1
      row <- row - 1
    }
    else if(my_matrix[row, col] == (my_matrix[row, col - 1] + gap_penalty)){
      number_of_gaps <- number_of_gaps + 1
      col <- col - 1
    }
    else{
      if(my_matrix[row, col] == (my_matrix[row - 1, col - 1] + match_penalty)){
        number_of_mismatch <- number_of_mismatch + 1
      }
      else{
        number_of_match <- number_of_match + 1
      }
      col <- col - 1
      row <- row - 1
    }
    
    
  }
  
  results[2] <- number_of_match
  results[3] <- number_of_mismatch
  results[4] <- number_of_gaps
  
  return(results)
}

create_score_matrix_semi_global <- function(seq1_vec, seq2_vec, match = 2, mismatch = -1, gap = -4) {
  
  # 1. Initialize the Matrix
  # Rows = seq1 length + 1 (for the gap row)
  # Cols = seq2 length + 1 (for the gap col)
  rows <- length(seq1_vec) + 1
  cols <- length(seq2_vec) + 1
  
  # Create a matrix filled with NA
  score_mat <- matrix(NA, nrow = rows, ncol = cols)
  
  # Optional: Label rows/cols for easier reading
  rownames(score_mat) <- c("Gap", seq1_vec)
  colnames(score_mat) <- c("Gap", seq2_vec)
  
  # 2. Base Cases (Initialize Edges)
  
  # Top-left corner is 0
  score_mat[1, 1] <- 0
  
  # Fill First Column (Vertical / i axis)
  for (i in 2:rows) {
    score_mat[i, 1] <- 0
  }
  
  # Fill First Row (Horizontal / j axis)
  for (j in 2:cols) {
    score_mat[1, j] <- 0
  }
  
  # 3. Fill the rest of the matrix
  for (i in 2:rows) {
    for (j in 2:cols) {
      
      # A. Calculate Diagonal Score (Match or Mismatch)
      # Note: We use i-1 and j-1 for vectors because the matrix has an extra "Gap" row/col
      if (seq1_vec[i-1] == seq2_vec[j-1]) {
        score_diag <- score_mat[i-1, j-1] + match
      } else {
        score_diag <- score_mat[i-1, j-1] + mismatch
      }
      
      # B. Calculate Up Score (Gap in seq2)
      score_up <- score_mat[i-1, j] + gap
      
      # C. Calculate Left Score (Gap in seq1)
      score_left <- score_mat[i, j-1] + gap
      
      # D. Fill the cell with the MAX score
      score_mat[i, j] <- max(score_diag, score_up, score_left)
    }
  }
  
  return(score_mat)
}

# --- TEST THE CODE ---
find_best_route_semi_global <- function(my_matrix, gap_penalty, match_score, match_penalty){
  
  
  number_of_match <- 0
  number_of_gaps <- 0
  number_of_mismatch <- 0
  col <- length(my_matrix[1, ])
  row <- length(my_matrix[ ,1])
  col_best <- 0
  row_best <- 0
  results <- 0
  max <- -1000
  
  for(i in 1:row){
    if(max < my_matrix[i, col]){
      max <- my_matrix[i, col]
      col_best <- col
      row_best <- i
    }
  }
  for(i in 1:col){
    if(max < my_matrix[row, i]){
      max <- my_matrix[row, i]
      col_best <- i
      row_best <- row
    }
  }
  col <- col_best
  row <- row_best
  results[1] <- my_matrix[row, col]
  while(col != 1 && row != 1){
    if(col == 1){
      number_of_gaps <- number_of_gaps + 1
      row <- row - 1
    }
    if(row == 1){
      number_of_gaps <- number_of_gaps + 1
      col <- col - 1
    }
    
    
    
    if(my_matrix[row, col] == (my_matrix[row - 1, col] + gap_penalty)){
      number_of_gaps <- number_of_gaps + 1
      row <- row - 1
    }
    else if(my_matrix[row, col] == (my_matrix[row, col - 1] + gap_penalty)){
      number_of_gaps <- number_of_gaps + 1
      col <- col - 1
    }
    else{
      if(my_matrix[row, col] == (my_matrix[row - 1, col - 1] + match_penalty)){
        number_of_mismatch <- number_of_mismatch + 1
      }
      else{
        number_of_match <- number_of_match + 1
      }
      col <- col - 1
      row <- row - 1
    }
    
    
  }
  results[2] <- number_of_match
  results[3] <- number_of_mismatch
  results[4] <- number_of_gaps
  
  return(results)
}

pairwiseAlignmentResult <- function(pattern, subject, type, substitutionMatrix, gapOpening, gapExtension) {
  
  # 1. Perform the alignment
  # We pass the arguments directly to the Bioconductor function
  # Note: Ensure library(pwalign) or library(Biostrings) is loaded before running
  alignment_object <- pairwiseAlignment(
    pattern = pattern,
    subject = subject,
    type = type,
    substitutionMatrix = substitutionMatrix,
    gapOpening = gapOpening,
    gapExtension = gapExtension
  )
  
  # 2. Extract the Score
  # The score() function extracts the optimal alignment score
  final_score <- score(alignment_object)
  
  # 3. Extract the Alignment Length
  # nchar() on the alignment object returns the total length (matches + mismatches + gaps)
  alignment_len <- nchar(alignment_object)
  
  # 4. Handle Edge Case (Null alignment length)
  # To avoid division by zero
  if (alignment_len == 0) {
    return(0)
  }
  
  # 5. Return Normalized Score
  return(final_score / alignment_len)
}


#########CODES TESTED #########

my_matrix_global <- create_score_matrix_global(string2, string1)
print(my_matrix_global)
results_global <- find_best_route_global(my_matrix_global, -4, 2, -1)
print(results_global)

my_matrix_semi_global <- create_score_matrix_semi_global(string2, string1)
print(my_matrix_semi_global)
results_semi_global <- find_best_route_semi_global(my_matrix_semi_global, -4, 2, -1)
print(results_semi_global)

############ SARS COV 2 ALLIGNEMENT################ 


MW090859 <- read.GenBank("MW090859.1")
ref_string <- toupper(paste(as.character(MW090859)[[1]], collapse = ""))
cleaned_SARS_2 <- SARS_2Example_function()
str(cleaned_SARS_2)

# 3. Print the result
sub_mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
str(cleaned_SARS_2)
# 1. Initialize storage
all_scores <- numeric(length(cleaned_SARS_2))

# 2. Corrected Loop
for(i in 1:length(cleaned_SARS_2)){
  
  # Step A: Get the raw data
  temp_seq <- cleaned_SARS_2[[i]]
  
  # Step B: FORCE the class to be DNAbin (The Magic Fix)
  # This tells R: "Treat these numbers as DNA, not math!"
  class(temp_seq) <- "DNAbin"
  
  # Step C: Now convert to string (It will correctly become "ATCG...")
  current_subject <- toupper(paste(as.character(temp_seq), collapse = ""))
  
  # DEBUG: Uncomment this once to make sure it looks like "ATCG" now
  # if (i == 1) print(substr(current_subject, 1, 50))
  
  # Step D: Run Alignment
  all_scores[i] <- pairwiseAlignmentResult(
    pattern = ref_string, 
    subject = current_subject, 
    type = "global", 
    substitutionMatrix = sub_mat, 
    gapOpening = 10000, 
    gapExtension = 2
  )
}

# 3. Print Final Results
print(head(all_scores))
cat("Average Normalized Score:", mean(all_scores), "\n")


resulted <- pairwiseAlignment(
  pattern = toupper(paste(string1, collapse = '')), 
  subject = toupper(paste(string2, collapse = '')),
  type = "global", 
  substitutionMatrix = sub_mat, 
  gapOpening = 10000, 
  gapExtension = 2
)

str(resulted)
final_score <- score(resulted)
alignment_len <- nchar(resulted)
?pairwiseAlignment()
