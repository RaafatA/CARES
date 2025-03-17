#####################################################
#####################################################
##
##
##    Rafat A. Eissa
##    Venn Diagram Plotting 
##    Fish Diet Project 
##    17 March 2025 
##
##
#####################################################
#####################################################


library(VennDiagram)
library(ggvenn)
# Load data
otu_data <- read.csv("Results/denoise-dada2/table-abundace-ASVs.biom/table-asv-abundance.csv")

# Extract OTU presence (non-zero abundances) for each sample
otu_lists <- lapply(otu_data[, -1], function(x) otu_data$OTU_ID[x > 0])

# Name the lists by sample
names(otu_lists) <- colnames(otu_data)[-1]

# Generate Venn diagram
venn.plot <- venn.diagram(
  x = otu_lists,
  filename = NULL,  # Use NULL to plot directly in RStudio
  fill = c("#0073C2FF", "#EFC000FF", "#868686FF", "yellow", "purple"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  main = "Shared OTUs Across Samples"
)

# Display the plot
grid.draw(venn.plot)




# Install and load ggvenn package
if (!requireNamespace("ggvenn", quietly = TRUE)) {
  install.packages("ggvenn")
}
library(ggvenn)

# Load data
otu_data <- read.csv("Results/denoise-dada2/table-abundace-ASVs.biom/table-asv-abundance.csv")

# Extract OTU presence (non-zero abundances) for each sample
otu_lists <- lapply(otu_data[, -1], function(x) otu_data$OTU_ID[x > 0])

# Name the lists by sample
names(otu_lists) <- colnames(otu_data)[-1]

# Generate the Venn diagram
ggvenn(
  data = otu_lists,
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "yellow", "purple"),
  stroke_size = 0.5,
  set_name_size = 5,
  text_size = 3)


# Install and load the venn package
if (!requireNamespace("venn", quietly = TRUE)) {
  install.packages("venn")
}
library(venn)

# Install and load the venn package
if (!requireNamespace("venn", quietly = TRUE)) {
  install.packages("venn")
}
library(venn)

# Load data
otu_data <- read.csv("Results/denoise-dada2/table-abundace-ASVs.biom/table-asv-abundance.csv")

# Prepare lists of OTUs for each sample
# Each sample will have a list of OTUs with non-zero abundance
otu_lists <- list()
for (col in colnames(otu_data)[-1]) {
  otu_lists[[col]] <- otu_data$OTU_ID[otu_data[[col]] > 0]
}
otu_lists$NT10.raw
# Check if the lists are correctly populated
print(otu_lists)

# Generate the Venn diagram
venn_result <- venn(
  otu_lists,
  zcolor = "style", # Use a styled color scheme
  ilabels = TRUE,   # Include intersection labels
  snames = names(otu_lists) # Set names from the list
)

# Display intersection details (useful for debugging)
print(attr(venn_result, "intersections"))




#############################################################################3
library(dplyr)
library(tidyr)

# Load the ASV abundance table
data <- read.csv("Results/denoise-dada2/table-abundace-ASVs.biom/table-asv-abundance.csv")

# Extract relevant columns
taxa_levels <- c("Phylum", "Class", "Order", "Family", "Genus")
sample_cols <- grep(".raw$", names(data), value = TRUE)

# Create a table of shared taxa across samples
shared_taxa_table <- data %>%
  pivot_longer(cols = all_of(sample_cols), names_to = "Sample", values_to = "Abundance") %>%
  filter(Abundance > 0) %>%
  group_by(across(all_of(taxa_levels)), Sample) %>%
  summarise(Present = n(), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Present, values_fill = 0) %>%
  mutate(Shared_Across_All = rowSums(select(., all_of(sample_cols)) > 0) == length(sample_cols))

# Save the result to a CSV file
write.csv(shared_taxa_table, "shared_taxa_table.csv", row.names = FALSE)

cat("Shared taxa table has been saved as 'shared_taxa_table.csv'.")




setwd("/home/raafat/Documents/CARES/Fish_intestine")



# Install and load required packages
if (!requireNamespace("UpSetR", quietly = TRUE)) {
  install.packages("UpSetR")
}
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

library(UpSetR)
library(tidyverse)

# Load ASV table with taxonomic annotations
file_path <- "Results/denoise-dada2/table-asv-abundance with taxa.csv"
otu_data <- read.csv(file_path)

# Set the taxonomic level to analyze (modify as needed)
taxonomy_level <- "Family"  # Options: "Phylum", "Class", "Order", "Family", "Genus"

# Remove missing values in the selected taxonomic level
otu_data <- otu_data %>%
  filter(!is.na(.data[[taxonomy_level]]))

# Aggregate ASV abundances by the selected taxonomy level
taxa_data <- otu_data %>%
  group_by(.data[[taxonomy_level]]) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  ungroup()

# Convert to presence/absence (1 if taxon is present, 0 if absent)
taxa_binary <- taxa_data %>%
  mutate(across(where(is.numeric), ~ ifelse(. > 0, 1, 0)))

# Convert to a format suitable for UpSetR
upset_data <- taxa_binary %>%
  column_to_rownames(var = taxonomy_level)  # Use taxonomic level as row names

# Generate UpSet plot for shared taxa
upset(upset_data, 
      sets = colnames(upset_data),  # Sample names as sets
      sets.bar.color = "#0073C2FF", 
      order.by = "freq",  # Order intersections by frequency
      mainbar.y.label = paste("Number of Shared", taxonomy_level, "Taxa"), 
      sets.x.label = "Samples",
      keep.order = TRUE)
