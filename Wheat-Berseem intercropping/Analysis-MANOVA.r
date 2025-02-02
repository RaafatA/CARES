# Load necessary libraries
library(agricolae)  # For post-hoc tests
library(dplyr)      # For data manipulation


# Group traits into components
# Vegetative Growth Component
vegetative_growth <- data %>% select(Plant_height, Leaf_length, Leaf_width, Leaf_area, Stalk_diameter)

# Tillering and Density Component
tillering_density <- data %>% select(Plant_per_hill, Tillers_per_plant, Tillers_per_hill)

# Yield Component
yield <- data %>% select(Shoot_dry_weight, Spike_length, Spike_width, Root_length, Root_weight)

# Perform MANOVA for Vegetative Growth Component with multiple factors
manova_vegetative <- manova(as.matrix(vegetative_growth) ~ Treatment * Irrigation * Genotype, data = data)
summary(manova_vegetative, test = "Pillai")  # Use Pillai's trace test

# Perform MANOVA for Tillering and Density Component with multiple factors
manova_tillering <- manova(as.matrix(tillering_density) ~ Treatment * Irrigation * Genotype, data = data)
summary(manova_tillering, test = "Pillai")

# Perform MANOVA for Yield Component with multiple factors
manova_yield <- manova(as.matrix(yield) ~ Treatment * Irrigation * Genotype, data = data)
summary(manova_yield, test = "Pillai")

# Post-hoc analysis (if MANOVA is significant)
# Use agricolae for post-hoc tests on individual traits
if (summary(manova_vegetative, test = "Pillai")$stats[1, "Pr(>F)"] < 0.05) {
  print("Significant MANOVA for Vegetative Growth Component")
  # Perform post-hoc tests for each trait in the component
  for (trait in colnames(vegetative_growth)) {
    print(paste("Post-hoc test for:", trait))
    model <- aov(data[[trait]] ~ Treatment * Irrigation * Genotype, data = data)
    print(HSD.test(model, "Treatment", group = TRUE))
    print(HSD.test(model, "Irrigation", group = TRUE))
    print(HSD.test(model, "Genotype", group = TRUE))
  }
}

if (summary(manova_tillering, test = "Pillai")$stats[1, "Pr(>F)"] < 0.05) {
  print("Significant MANOVA for Tillering and Density Component")
  for (trait in colnames(tillering_density)) {
    print(paste("Post-hoc test for:", trait))
    model <- aov(data[[trait]] ~ Treatment * Irrigation * Genotype, data = data)
    print(HSD.test(model, "Treatment", group = TRUE))
    print(HSD.test(model, "Irrigation", group = TRUE))
    print(HSD.test(model, "Genotype", group = TRUE))
  }
}

if (summary(manova_yield, test = "Pillai")$stats[1, "Pr(>F)"] < 0.05) {
  print("Significant MANOVA for Yield Component")
  for (trait in colnames(yield)) {
    print(paste("Post-hoc test for:", trait))
    model <- aov(data[[trait]] ~ Treatment * Irrigation * Genotype, data = data)
    print(HSD.test(model, "Treatment", group = TRUE))
    print(HSD.test(model, "Irrigation", group = TRUE))
    print(HSD.test(model, "Genotype", group = TRUE))
  }
}
