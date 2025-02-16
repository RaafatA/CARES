###########################################################
##
##   Project:  Wheat-Berseem intercropping
##   Date:    02/02/2025
##   Author:  Rafat A. Eissa
##   MANOVA and post-hoc
##
###########################################################











# Load necessary libraries
library(agricolae)  # For post-hoc tests
library(dplyr)      # For data manipulation
setwd("C:/Users/raafat.abdulmajeed/Documents/Wheat - Barseem Interaction")



Anth_Data = read.csv("Anthesis-Data.csv", header =T)
Yield_Data= read.csv("Yield-Data.csv", header =T)

Yield_Data$Ear_Length
# Group traits into components
# Vegetative Growth Component
vegetative_growth <- Yield_Data %>% select(PH, LL, LW, LA,Ear_Length) # nolint

# Tillering and Density Component
tillering_density <- Yield_Data %>% select(Plant_per_hill, Tillers_per_plant, Tillers_per_hill) # nolint # nolint

# Yield Component
yield <- data %>% select(Shoot_dry_weight, Spike_length, Spike_width, Root_length, Root_weight) # nolint

# Perform MANOVA for Vegetative Growth Component with multiple factors
manova_vegetative <- manova(as.matrix(vegetative_growth) ~ Water_quality*Field_capacity*Cropping_system, data = Yield_Data) # nolint
summary(manova_vegetative, test = "Pillai")  # Use Pillai's trace test

# Perform MANOVA for Vegetative Growth Component with multiple factors
manova_vegetative <- manova(as.matrix(vegetative_growth) ~ Water_quality*Field_capacity*Cropping_system, data = Anth_Data) # nolint
summary(manova_vegetative, test = "Pillai")  # Use Pillai's trace test

# Perform ANOVA on each dependent variable to identify significant factors
anova_results <- summary.aov(manova_vegetative)
print(anova_results)

# Post-hoc test (Tukey's HSD) for significant factors
library(multcomp)

tukey_results <- lapply(names(vegetative_growth), function(var) {
  aov_model <- aov(as.formula(paste(var, "~ Water_quality*Field_capacity*Cropping_system")), data = Yield_Data)
  TukeyHSD(aov_model)
})

names(tukey_results) <- names(vegetative_growth)
print(tukey_results)

# Load necessary libraries
library(multcomp)

# Perform MANOVA for Vegetative Growth Component with multiple factors
manova_vegetative <- manova(as.matrix(vegetative_growth) ~ Water_quality*Field_capacity*Cropping_system, data = Yield_Data)
manova_summary <- summary(manova_vegetative, test = "Pillai")

# Save MANOVA summary to a file
sink("MANOVA_results.txt")
cat("### MANOVA Results (Pillai's Trace Test) ###\n")
print(manova_summary)
sink()

# Perform ANOVA on each dependent variable
anova_results <- summary.aov(manova_vegetative)

# Save ANOVA results to a file
sink("ANOVA_results-Harvest-Vegetative.txt")
cat("### ANOVA Results for Each Dependent Variable ###\n")
print(anova_results)
sink()

# Post-hoc test (Tukey's HSD) for significant factors
tukey_results <- lapply(names(vegetative_growth), function(var) {
  aov_model <- aov(as.formula(paste(var, "~ Water_quality*Field_capacity*Cropping_system")), data = Yield_Data)
  TukeyHSD(aov_model)
})

names(tukey_results) <- names(vegetative_growth)

# Save Tukey's HSD results to a file in an organized manner
sink("TukeyHSD_results-harvest-vegetative.txt")
cat("### Tukey's HSD Post-hoc Test Results ###\n\n")
for (var in names(tukey_results)) {
  cat(paste0("### Dependent Variable: ", var, " ###\n"))
  print(tukey_results[[var]])
  cat("\n--------------------------------------------------\n")
}
sink()

# Confirm completion
cat("Results have been saved to 'MANOVA_results.txt', 'ANOVA_results.txt', and 'TukeyHSD_results.txt'.\n")







# Perform MANOVA for Tillering and Density Component with multiple factors
manova_tillering <- manova(as.matrix(tillering_density) ~ Treatment * Irrigation * Genotype, data = data) # nolint
summary(manova_tillering, test = "Pillai")

# Perform MANOVA for Yield Component with multiple factors
manova_yield <- manova(as.matrix(yield) ~ Treatment * Irrigation * Genotype, data = data) # nolint
summary(manova_yield, test = "Pillai")
tukey_results <- TukeyHSD(manova_vegetative)
print(tukey_results)
# Post-hoc analysis (if MANOVA is significant)
# Use agricolae for post-hoc tests on individual traits
if (summary(manova_vegetative, test = "Pillai")$stats[1, "Pr(>F)"] < 0.05) {
  print("Significant MANOVA for Vegetative Growth Component")
  # Perform post-hoc tests for each trait in the component
  for (trait in colnames(vegetative_growth)) {
    print(paste("Post-hoc test for:", trait))
    model <- aov(data[[trait]] ~ Water_quality*Field_capacity*Cropping_system, data = Anth_Data)
    print(HSD.test(model, "Water_quality", group = TRUE))
    print(HSD.test(model, "Field_capacity", group = TRUE))
    print(HSD.test(model, "Cropping_system", group = TRUE))
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






#######################################################################################################
#######################################################################################################


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
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold", size = rel(1.1))
    )
}

# Perform LSD test for Treatment
anova_treatment <- aov(Leaf_width ~ Treatment, data = data)
lsd_treatment <- LSD.test(anova_treatment, "Treatment", alpha = 0.05, group = TRUE)
lsd_treatment_groups <- lsd_treatment$groups
lsd_treatment_groups$Treatment <- rownames(lsd_treatment_groups)

# Merge LSD grouping with the original data
data <- data %>%
  left_join(lsd_treatment_groups, by = "Treatment")

# Boxplot: Leaf_width vs. Treatment with LSD grouping
p1 <- ggplot(data, aes(x = Treatment, y = Leaf_width, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
  geom_text(aes(label = groups, y = max(Leaf_width) + 0.2), 
            position = position_dodge(width = 0.75), vjust = -0.5, size = 5) +
  labs(title = "Effect of Treatment on Leaf Width",
       x = "Treatment",
       y = "Leaf Width (cm)") +
  theme_publication() +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c"))  # Custom colors

# Perform LSD test for Genotype (grouped by Treatment)
anova_genotype <- aov(Leaf_width ~ Treatment + Genotype, data = data)
lsd_genotype <- LSD.test(anova_genotype, "Genotype", alpha = 0.05, group = TRUE)
lsd_genotype_groups <- lsd_genotype$groups
lsd_genotype_groups$Genotype <- rownames(lsd_genotype_groups)

# Merge LSD grouping with the summarized data
summary_data <- data %>%
  group_by(Treatment, Genotype) %>%
  summarise(Mean_Leaf_width = mean(Leaf_width, na.rm = TRUE),
            SE_Leaf_width = sd(Leaf_width, na.rm = TRUE) / sqrt(n())) %>%
  left_join(lsd_genotype_groups, by = "Genotype")

# Bar plot: Leaf_width vs. Genotype, grouped by Treatment with LSD grouping
p2 <- ggplot(summary_data, aes(x = Genotype, y = Mean_Leaf_width, fill = Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.8) +
  geom_errorbar(aes(ymin = Mean_Leaf_width - SE_Leaf_width, ymax = Mean_Leaf_width + SE_Leaf_width),
                position = position_dodge(width = 0.9), width = 0.2, size = 0.8) +
  geom_text(aes(label = groups, y = Mean_Leaf_width + SE_Leaf_width + 0.1), 
            position = position_dodge(width = 0.9), vjust = -0.5, size = 5) +
  labs(title = "Leaf Width by Genotype and Treatment",
       x = "Genotype",
       y = "Mean Leaf Width (cm)") +
  theme_publication() +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c"))  # Custom colors

# Arrange plots in a grid
final_plot <- ggarrange(p1, p2, ncol = 1, labels = c("A", "B"))
print(final_plot)

# Save the final plot for publication
ggsave("growth_data_visualization_with_LSD.png", final_plot, width = 10, height = 10, dpi = 300)
