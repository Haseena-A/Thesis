install.packages("bnlearn")
library(bnlearn)

# Load your data from the specified path
data3 <- read.csv("C:/Users/bdeva/Desktop/Haseena/Mesenchymal/Stochastic/final_top_correlation_adipose_stochastic.csv")

# Check the structure of the data
str(data3)

# Convert integer and character columns to factors
data3 <- as.data.frame(lapply(data3, function(x) {
  if (is.integer(x) || is.character(x)) as.factor(x) else x
}))

# Verify the conversion
str(data3)

# Perform bootstrapping to estimate the strength of network arcs
bay_main2_updated_adj_values_11 <- boot.strength(data = data3, R = 5, algorithm = "hc")

# Save the result to a CSV file
write.csv(bay_main2_updated_adj_values_11, file = "C:/Users/bdeva/Desktop/Haseena/Mesenchymal/Stochastic/bn5_final_top_correlation_adipose_stochastic.csv", row.names = FALSE)

