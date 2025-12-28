library("ape")
library("Biostrings")
library("pwalign")
sequences_dna <- read.GenBank("MW090859.1")
final_seq_MW090859.1 <- toupper(paste(as.character(sequences_dna)[[1]], collapse = ""))

sequences_dna <- read.GenBank("Y10051.1")
final_seq_Y10051.1 <- toupper(paste(as.character(sequences_dna)[[1]], collapse = ""))

sequences_dna <- read.GenBank("CY125943.1")
final_seq_CY125943.1 <- toupper(paste(as.character(sequences_dna)[[1]], collapse = ""))

sequences_dna <- read.GenBank("A22884.1")
final_seq_A22884.1 <- toupper(paste(as.character(sequences_dna)[[1]], collapse = ""))


pairwiseAlignmentResult = function(pattern, subject, type, substitutionMatrix, gapOpening, gapExtension) {
  align <- pairwiseAlignment(pattern, subject,
                             type = type,
                             substitutionMatrix = substitutionMatrix,
                             gapOpening = gapOpening,
                             gapExtension = gapExtension)
  L <- nchar(alignedPattern(align))
  if (L == 0) return(0)
  score(align) / L
}

pairwiseSignificance <- function(alpha, pattern, subject, type, substitutionMatrix, gapOpening, gapExtension){
  
  result_one <- pairwiseAlignmentResult(pattern, subject, type, substitutionMatrix, gapOpening, gapExtension)
  number_of_samples <- 500
  result_list <- numeric(number_of_samples)
  subject_length <- nchar(subject)
  
  subject_chars <- unlist(strsplit(subject, ""))
  for(i in (1:number_of_samples)){
    
    sample_seq <- paste(sample(subject_chars), collapse = "")
    result_list[i] <- pairwiseAlignmentResult(pattern, sample_seq, type, substitutionMatrix, gapOpening, gapExtension)
    
  }
  
  #relevant if i have a p value of less than 5%. 0.05 of all samples. 
  
  p_value <- (sum((result_list > result_one)) + 1)/(number_of_samples+1)
  print(result_one)
  print(p_value)
  print(result_list)
  return(p_value <= alpha)
}
sub_mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)

number <- 0

for(a in c(final_seq_Y10051.1, final_seq_CY125943.1, final_seq_A22884.1)){
  if(pairwiseSignificance(0.05, final_seq_MW090859.1, a, "global", sub_mat, 5, 2)){
    number <- number + 1
  }
}


