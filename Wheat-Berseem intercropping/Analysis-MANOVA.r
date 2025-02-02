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


#########################################################################################################
#########################################################################################################
## Visualization 
##
#########################################################################################################
#########################################################################################################
library(ggplot2)    
library(dplyr)      
library(ggpubr)     
library(scales)     



# Custom theme for publication-quality plots
theme_publication <- function(base_size = 14, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(color = "black"),
      axis.title = element_text(face = "bold", size = rel(1.2)),
      axis.text = element_text(size = rel(1)),
      axis.line = element_line(color = "black"),
      panel.grid.major = element_line(linetype = "dotted", color = "gray80"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold", size = rel(1.1))
    )
}

# Boxplot: Leaf_width vs. Treatment with statistical annotations
p1 <- ggplot(data, aes(x = Treatment, y = Leaf_width, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
  stat_compare_means(method = "anova", label = "p.format", label.x = 1.5, label.y = max(data$Leaf_width) + 0.5) +
  labs(title = "Effect of Treatment on Leaf Width",
       x = "Treatment",
       y = "Leaf Width (cm)") +
  theme_publication() +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c"))  # Custom colors

# Line plot: Leaf_width over Time, faceted by Treatment and Irrigation
p2 <- data %>%
  group_by(Treatment, Irrigation, Time) %>%
  summarise(Mean_Leaf_width = mean(Leaf_width, na.rm = TRUE),
            SE_Leaf_width = sd(Leaf_width, na.rm = TRUE) / sqrt(n())) %>%
  ggplot(aes(x = Time, y = Mean_Leaf_width, color = Treatment, group = Treatment)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Mean_Leaf_width - SE_Leaf_width, ymax = Mean_Leaf_width + SE_Leaf_width),
                width = 0.2, size = 0.8) +
  facet_grid(~ Irrigation) +
  labs(title = "Leaf Width Over Time by Treatment and Irrigation",
       x = "Time (Days)",
       y = "Mean Leaf Width (cm)") +
  theme_publication() +
  scale_color_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c"))  # Custom colors

# Bar plot: Leaf_width vs. Genotype, grouped by Treatment
p3 <- data %>%
  group_by(Treatment, Genotype) %>%
  summarise(Mean_Leaf_width = mean(Leaf_width, na.rm = TRUE),
            SE_Leaf_width = sd(Leaf_width, na.rm = TRUE) / sqrt(n())) %>%
  ggplot(aes(x = Genotype, y = Mean_Leaf_width, fill = Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.8) +
  geom_errorbar(aes(ymin = Mean_Leaf_width - SE_Leaf_width, ymax = Mean_Leaf_width + SE_Leaf_width),
                position = position_dodge(width = 0.9), width = 0.2, size = 0.8) +
  labs(title = "Leaf Width by Genotype and Treatment",
       x = "Genotype",
       y = "Mean Leaf Width (cm)") +
  theme_publication() +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c"))  # Custom colors

# Arrange plots in a grid
final_plot <- ggarrange(p1, p2, p3, ncol = 1, labels = c("A", "B", "C"))
print(final_plot)

# Save the final plot for publication
ggsave("growth_data_visualization.png", final_plot, width = 10, height = 12, dpi = 300)
