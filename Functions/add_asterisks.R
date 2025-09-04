###This function just adds asterisks of significance based on p-value 

## Inputs:

## p_values=A vector of p-values

##Outputs:

## asterisks=A vector of asterisks corresponding to significance

add_asterisks <- function(p_values) {
  asterisks <- character(length(p_values))  # Initialize an empty character vector
  
  #It goes through the a column of p-values
  for (i in seq_along(p_values)) {
    # If p-value is below 0.001 three stars are provided
    if (p_values[i] < 0.001) {
      asterisks[i] <- "***"
    # If p-value is below not under 0.001 but is under 0.01 two stars are provided
    } else if (p_values[i] < 0.01) {
      asterisks[i] <- "**"
    # If p-value is below not under 0.01 but is under 0.05 one star is provided
    } else if (p_values[i] < 0.05) {
      asterisks[i] <- "*"
    # If p-value is below not under 0.05, "ns" provided to indicate not significant
    } else {
      asterisks[i] <- "ns"
    }
  }
  #Finally the asterisks are returned
  return(asterisks)
} 
