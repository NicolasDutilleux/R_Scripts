add_one <- function(x, y = 1){
  
  return (x + y)
  
}

x = 42

saveRDS(x, file = "my_data.rds")

getwd()
setwd(dir = "C:/Users/nicod/Documents/VisualNicolas/RStudioPhylo")
y <- seq(0, 60, 21)
saveRDS(y, file = "seq_result.rds")
print(y)
?seq
