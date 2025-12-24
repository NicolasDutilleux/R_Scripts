nb_sliding_windows <- function(l, size, shift){
  if (l < size) {
    return(0)
  }
  
  # This part of your code correctly adds the extra partial window
  number <- floor((l - size) / shift) + 1
  if((l-size) %% shift > 0){
    number <- number + 1
  }
  return(number)
}

frequency_in_sliding_window <- function(sequence, ntd, size, shift){
  
  l <- length(sequence)
  N <- nb_sliding_windows(l, size, shift)
  result <- numeric(N)
  
  for (i in 1:N){
    number <- 0
    
    # We still attempt to loop 'size' times...
    for (j in 1:size){ 
      
      # Calculate the current index
      current_idx <- ((i - 1) * shift) + j
      
      # SAFETY CHECK: If we went past the end of the sequence, STOP looking!
      if (current_idx > l) {
        break
      }
      
      if (sequence[[current_idx]] == ntd){
        number <- number + 1
      }
    }
    
    # Note: This divides by the full 'size' even for the short window.
    # If the window is only 1 base long, but size is 3, a match counts as 0.33.
    # If you want 100% for that case, change 'size' to 'min(size, l - start_index + 1)'
    percent <- number / (min(size, l - ((i - 1) * shift)))
    
    result[i] <- round(percent, 2)
  }
  
  return(result)
}