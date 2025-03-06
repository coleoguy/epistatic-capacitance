# Load necessary libraries
library(ggplot2)

# Read data
axa <- read.csv("df_axa.csv")
axd <- read.csv("df_axd (1).csv")
dxd <- read.csv("df_dxd (1).csv")

# Combine data into a single dataframe with an identifier
axa$group <- "axa"
axd$group <- "axd"
dxd$group <- "dxd"

# Merge all datasets
data <- rbind(axa[, c("loci", "mean_epi", "group")], 
              axd[, c("loci", "mean_epi", "group")], 
              dxd[, c("loci", "mean_epi", "group")])

# Plot the data
ggplot(data, aes(x = loci, y = mean_epi, color = group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Mean Epi Values Across Loci", 
       x = "Loci", 
       y = "Mean Epi",
       color = "Group") +
  theme_minimal()
