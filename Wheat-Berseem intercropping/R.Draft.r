###########################################################
##
##   Project:  Wheat-Berseem intercropping
##   Date:    02/02/2025
##   Author:  Rafat A. Eissa
##
###########################################################

# Load necessary package
library(openxlsx)

# Setting the Working Directory 
setwd("/home/raafat/Documents/CARES/Wheat-Berseem intercropping/Wheat - Barseem Interaction")
Anth_Data = read.csv("Anthesis-Data.csv", header =T)
Yield_Data = read.csv("Yield-Data.csv", header =T)

Factor_columns = c("Water_quality", "Field_capacity", "Cropping_system")

Traits = setdiff(names(Anth_Data), Factor_columns)

hist(Anth_Data$StalkDiameter..mm., col="blue",main = "Lead", xlab= "", vlab="Frequency")


attach(Anth_Data)
boxplot(PH..cm.~Water_quality*Field_capacity*Cropping_system, col="green", main="Fruit Weight", xlab="",ylab="Fruit Weight")



anova_results <- list()


# Loop through each trait and perform ANOVA
for (trait in Traits) {
  
  # Formula for ANOVA
  formula <- as.formula(paste(trait, "~ Water_quality*Field_capacity*Cropping_system"))
  
  # Perform ANOVA
  aov_result <- summary(aov(formula, data = Anth_Data))
  
  # Extract ANOVA table
  anova_table <- as.data.frame(aov_result[[1]])
  
  # Rename columns to match your example
  colnames(anova_table) <- c("Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")
  
  # Add a column for significance symbols based on p-values
  anova_table$Significance <- ifelse(anova_table$`Pr(>F)` < 0.001, "***",
                                     ifelse(anova_table$`Pr(>F)` < 0.01, "**",
                                            ifelse(anova_table$`Pr(>F)` < 0.05, "*", "ns")))
  
  # Add the Trait name at the top of the table to identify the trait in the final Excel sheet
  trait_header <- data.frame(Df = NA, Sum.Sq = NA, Mean.Sq = NA, `F value` = NA, `Pr(>F)` = NA, Significance = trait)
  
  # Ensure column names in the header match the ANOVA table
  colnames(trait_header) <- colnames(anova_table)
  
  # Combine the trait header with the ANOVA table for this trait
  trait_anova <- rbind(trait_header, anova_table)
  
  # Append to the results list
  anova_results[[trait]] <- trait_anova
}

# Combine all ANOVA tables into one data frame
combined_anova <- do.call(rbind, anova_results)

# Create a new Excel workbook
wb <- createWorkbook()

# Add a worksheet to store the combined ANOVA results
addWorksheet(wb, "ANOVA Anth_Data")

# Write the combined ANOVA table to the worksheet
writeData(wb, "ANOVA Anth_Data", combined_anova)

# Save the workbook to an Excel file
saveWorkbook(wb, "ANOVA Anth_Data.xlsx", overwrite = TRUE)

# Output message
cat("ANOVA results have been saved to formatted_anova_results.xlsx")






##############################################################################
####    Multivariate Stats
##############################################################################
setwd("/home/raafat/Documents/CARES/Wheat-Berseem intercropping/Wheat - Barseem Interaction")
Anth_Data = read.csv("Anthesis-Data.csv", header =T)
Yield_Data = read.csv("Yield-Data.csv", header =T)
library(dplyr)
library(agricolae)
Anth_Data$StalkDiameter..mm.
res.man = manova(cbind() ~Water_quality*Field_capacity*Cropping_system, data=Anth_Data)

# Group traits into components
# Vegetative Growth Component
vegetative_growth <- Anth_Data %>% select(PH..cm., LL..cm.,LW..cm.,StalkDiameter..mm., LA)

# Tillering and Density Component
tillering_density <- data %>% select(Plant_per_hill, Tillers_per_plant, Tillers_per_hill)

# Yield Component
yield <- data %>% select(Shoot_dry_weight, Spike_length, Spike_width, Root_length, Root_weight)

# Perform MANOVA for Vegetative Growth Component with multiple factors
manova_vegetative <- manova(as.matrix(vegetative_growth) ~ Water_quality*Field_capacity*Cropping_system, data = Anth_Data)
summary(manova_vegetative, test = "Pillai")  # Use Pillai's trace test

cor(vegetative_growth, use = "pairwise.complete.obs")
pca_res <- prcomp(vegetative_growth, scale. = TRUE)
summary(pca_res)

library(car)
vif(lm(as.matrix(vegetative_growth) ~ Water_quality, data = Anth_Data))


# Perform MANOVA for Tillering and Density Component with multiple factors
manova_tillering <- manova(as.matrix(tillering_density) ~ Treatment * Irrigation * Cropping_system, data = data)
summary(manova_tillering, test = "Pillai")

# Perform MANOVA for Yield Component with multiple factors
manova_yield <- manova(as.matrix(yield) ~ Treatment * Irrigation * Cropping_system, data = data)
summary(manova_yield, test = "Pillai")

# Post-hoc analysis (if MANOVA is significant)
# Use agricolae for post-hoc tests on individual traits
if (summary(manova_vegetative, test = "Pillai")$stats[1, "Pr(>F)"] < 0.05) {
  print("Significant MANOVA for Vegetative Growth Component")
  # Perform post-hoc tests for each trait in the component
  for (trait in colnames(vegetative_growth)) {
    print(paste("Post-hoc test for:", trait))
    model <- aov(Anth_Data[[trait]] ~ Water_quality*Field_capacity*Cropping_system, data = Anth_Data)
    print(HSD.test(model, "Water_quality", group = TRUE))
    print(HSD.test(model, "Field_capacity", group = TRUE))
    print(HSD.test(model, "Cropping_system", group = TRUE))
  }
}

if (summary(manova_tillering, test = "Pillai")$stats[1, "Pr(>F)"] < 0.05) {
  print("Significant MANOVA for Tillering and Density Component")
  for (trait in colnames(tillering_density)) {
    print(paste("Post-hoc test for:", trait))
    model <- aov(data[[trait]] ~ Treatment * Irrigation * Cropping_system, data = data)
    print(HSD.test(model, "Treatment", group = TRUE))
    print(HSD.test(model, "Irrigation", group = TRUE))
    print(HSD.test(model, "Genotype", group = TRUE))
  }
}

if (summary(manova_yield, test = "Pillai")$stats[1, "Pr(>F)"] < 0.05) {
  print("Significant MANOVA for Yield Component")
  for (trait in colnames(yield)) {
    print(paste("Post-hoc test for:", trait))
    model <- aov(data[[trait]] ~ Treatment * Irrigation * Cropping_system, data = data)
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
setwd("/home/raafat/Documents/CARES/Wheat-Berseem intercropping/Wheat - Barseem Interaction")
Anth_Data = read.csv("Anthesis-Data.csv", header =T)
Yield_Data = read.csv("Yield-Data.csv", header =T)

library(ggplot2)    
library(dplyr)      
library(ggpubr)     
library(scales)     
library(ggpattern)
install.packages("ggsignif")
library(ggsignif)
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
Anth_Data$PH
# Boxplot: Leaf_width vs. Treatment with statistical annotations
p1 <- ggplot(Anth_Data, aes(x = Water_quality, y = LA, fill = Cropping_system)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
  facet_wrap(~Cropping_system)+
  facet_grid(~ Field_capacity) +
  stat_compare_means(method = "anova", label = "p.format", label.x = 1.5, label.y = max(Anth_Data$LA) + 0.5) +
  labs(title = "Effect of Treatment on Leaf Area",
       x = "Treatment",
       y = "Leaf Area (cm²)") +
  theme_publication() +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c", "#ff7f9e"))  # Custom colors
p1
# 
# p3 <- Anth_Data %>%
#   group_by(Water_quality, Cropping_system, Field_capacity) %>%
#   summarise(Mean_Plant_Height = mean(PH, na.rm = TRUE),
#             SE_Plant_Height = sd(PH, na.rm = TRUE) / sqrt(n())) %>%
#   ggplot(aes(x = Water_quality, y = Mean_Plant_Height, fill = Water_quality)) +
#   geom_bar(stat = "identity", position = position_dodge(), alpha = 0.8) +
#   geom_errorbar(aes(ymin = Mean_Plant_Height - SE_Plant_Height, ymax = Mean_Plant_Height + SE_Plant_Height),
#                 position = position_dodge(width = 0.9), width = 0.2, size = 0.8) +
#   facet_grid(~ Field_capacity*Cropping_system) +
#   labs(title = "Stalk Diameter by Water Quality and Field Capacity an ",
#        x = "Water Quality",
#        y = "Stalk Diameter (mm)") +
#   geom_col_pattern(
#     aes(
#       pattern_fill= Water_quality,
#       fill = Water_quality,
#       pattern= Water_quality
#     ),
#     colour = "black",
#     fill = 'white',
#     pattern_fill = 'black',
#     pattern_density = 0.2
#   )+
#   theme_publication() #+
#   #scale_fill_manual(values = c("black", "black", "black", "black"))  # Custom colors
# p3
# Arrange plots in a grid
# final_plot <- ggarrange(p1, p2, p3, ncol = 1, labels = c("A", "B", "C"))
# print(final_plot)
# 
# # Save the final plot for publication
# ggsave("growth_data_visualization.png", final_plot, width = 10, height = 12, dpi = 300)
# 



#######################################################################################################
#######################################################################################################
setwd("/home/raafat/Documents/CARES/Wheat-Berseem intercropping/Wheat - Barseem Interaction")
Anth_Data = read.csv("Anthesis-Data.csv", header =T)
Yield_Data = read.csv("Yield-Data.csv", header =T)

library(ggplot2)    
library(dplyr)      
library(ggpubr)     
library(scales)     
library(ggpattern)
library(agricolae)

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
Anth_Data$Cropping_system
# Perform LSD test for Treatment
anova_treatment <- aov(PH ~ Water_quality*Field_capacity*Cropping_system, data = Anth_Data)
lsd_treatment <- LSD.test(anova_treatment, "Water_quality", alpha = 0.05, group = TRUE)
lsd_treatment_groups <- lsd_treatment$groups
lsd_treatment_groups$Water_quality <- rownames(lsd_treatment_groups)

# Merge LSD grouping with the original data
data <- Anth_Data %>%
  left_join(lsd_treatment_groups, by = "Water_quality")

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






library(ggplot2)    
library(dplyr)      
library(ggpubr)     
library(scales)     
library(ggpattern)
library(agricolae)
library(tibble)
library(tidyr)
setwd("/home/raafat/Documents/CARES/Wheat-Berseem intercropping/Wheat - Barseem Interaction")
Anth_Data = read.csv("Anthesis-Data.csv", header =T)
Yield_Data = read.csv("Yield-Data.csv", header =T)

anova_model <- aov(PH ~ Water_quality * Field_capacity * Cropping_system, data = Anth_Data)
summary(anova_model)
lsd_results <- LSD.test(anova_model, c("Water_quality","Field_capacity","Cropping_system"), alpha = 0.05, group = TRUE)
lsd_Water_quality_groups <- lsd_results$groups
lsd_Water_quality_groups 

lsd_results$groups <- lsd_results$groups %>%
  rownames_to_column(var = "Water_quality_Field_capacity_Cropping_system") %>%  # Convert row names to a column
  separate(Water_quality_Field_capacity_Cropping_system, into = c("Water_quality", "Field_capacity","Cropping_system"), sep = ":", remove = FALSE)
# Merge LSD grouping with the summarized data
summary_data <- Anth_Data %>%
  group_by(Water_quality, Field_capacity,Cropping_system) %>%
  summarise(Mean_trait = mean(PH, na.rm = TRUE),
            SE_Trait = sd(PH, na.rm = TRUE) / sqrt(n())) 

summary_data <- summary_data %>%
  left_join(lsd_results$groups, by = c("Water_quality","Field_capacity","Cropping_system"))

# Bar plot: Leaf_width vs. Genotype, grouped by Treatment with LSD grouping
p2 <- ggplot(summary_data, aes(x = Water_quality, y = Mean_trait, fill = Water_quality)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.8) +
  geom_errorbar(aes(ymin = Mean_trait - SE_Trait, ymax = Mean_trait + SE_Trait),
                position = position_dodge(width = 0.9), width = 0.2, size = 0.8) +
  geom_text(aes(label = groups, y = Mean_trait + SE_Trait + 0.1), 
            position = position_dodge(width = 0.9), vjust = -0.5, size = 5) +
  facet_grid(~ Field_capacity*Cropping_system) +
  labs(title = "(A)",
       x = "Water Quality",
       y = "Plant Height (cm)") +
  theme_publication() +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c","#1F4529"))  # Custom colors
p2
# Arrange plots in a grid
final_plot <- ggarrange(p1, p2, ncol = 1, labels = c("A", "B"))
print(final_plot)


# Boxplot: Leaf_width vs. Treatment with LSD grouping
p1 <- ggplot(Anth_Data, aes(x = Water_quality, y = PH, fill = Water_quality)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
  geom_text(aes(label = groups, y = max(PH) + 0.2), 
            position = position_dodge(width = 0.75), vjust = -0.5, size = 5) +
  labs(title = "Effect of Treatment on Leaf Width",
       x = "Treatment",
       y = "Leaf Width (cm)") +
  facet_grid(~ Field_capacity*Cropping_system) +
  theme_publication() +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c","#1F4529"))  # Custom colors

p1# Perform LSD test for Genotype (grouped by Treatment)





# Save the final plot for publication
ggsave("growth_data_visualization_with_LSD.png", final_plot, width = 10, height = 10, dpi = 300)

################################################################################################
################################################################################################
################################################################################################
################################################################################################
setwd("/home/raafat/Documents/CARES/Wheat-Berseem intercropping/Wheat - Barseem Interaction")
Anth_Data = read.csv("Anthesis-Data.csv", header =T)
Yield_Data = read.csv("Yield-Data.csv", header =T)

library(ggplot2)    
library(dplyr)      
library(ggpubr)     
library(scales)     
library(ggpattern)
library(agricolae)

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




anova_model <- aov(PH ~ Water_quality*Treatment, data = Anth_Data)
summary(anova_model)
lsd_results <- LSD.test(anova_model, c("Genotype","Treatment"), alpha = 0.05, group = TRUE)
lsd_Water_quality_groups <- lsd_results$groups
lsd_Water_quality_groups 

lsd_results$groups <- lsd_results$groups %>%
  rownames_to_column(var = "Genotype_Treatment") %>%  # Convert row names to a column
  separate(Genotype_Treatment, into = c("Genotype", "Treatment"), sep = ":", remove = FALSE)


# Merge LSD grouping with the summarized data
summary_data <- df %>%
  group_by(Cropping_system, Treatment) %>%
  summarise(Mean_trait = mean(APX, na.rm = TRUE),
            SE_Trait = sd(APX, na.rm = TRUE) / sqrt(n())) 
colnames(lsd_results$groups)

summary_data <- summary_data %>%
  left_join(lsd_results$groups, by = c("Genotype","Treatment"))
# Bar plot: Leaf_width vs. Genotype, grouped by Treatment with LSD grouping
p2 <- ggplot(summary_data, aes(x = Cropping_system, y = Mean_trait, fill = Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.9,color = "black", linewidth = 0.7) +
  geom_errorbar(aes(ymin = Mean_trait - SE_Trait, ymax = Mean_trait + SE_Trait),
                position = position_dodge(width = 0.9), width = 0.2, linewidth = 0.8) +
  geom_text(aes(label = groups, y = Mean_trait + SE_Trait + 0.1), 
            position = position_dodge(width = 0.9), vjust = -0.5, size = 5) +
  labs(title = "Leaf Width by Genotype and Treatment",
       x = "Cultivars",
       y = "Ascorbate Peroxidase (mM / gFW) ") +
  theme_publication() +
  scale_fill_manual(values = custom_colors)   # Custom colors
p2
# Arrange plots in a grid
final_plot <- ggarrange(p1, p2, ncol = 1, labels = c("A", "B"))
print(final_plot)

# Save the final plot for publication
ggsave("growth_data_visualization_with_LSD.png", final_plot, width = 10, height = 10, dpi = 300)

#################################################################################################

# Load necessary libraries
library(ggplot2)
library(factoextra)  # For PCA visualization
library(dplyr)       # For data manipulation

Anth_Data$Cropping_system <- as.factor(Anth_Data$Cropping_system)
Anth_Data$Field_capacity <- as.factor(Anth_Data$Field_capacity)
Anth_Data$Water_quality <- as.factor(Anth_Data$Water_quality)


factor_columns <- c("Water_quality", "Field_capacity", "Cropping_system")
trait_data <- Anth_Data %>% select(-one_of(factor_columns))

pca_result <- prcomp(trait_data, scale. = TRUE)
summary(pca_result)
pca_scores <- as.data.frame(pca_result$x)

# Combine PCA scores with genotype and treatment for plotting
pca_scores$Cropping_system <- Anth_Data$Cropping_system
pca_scores$Field_capacity <- Anth_Data$Field_capacity
pca_scores$Water_quality <- Anth_Data$Water_quality
# Create a PCA biplot using ggplot2 and color by Genotype and shape by Treatment
ggplot(pca_scores, aes(x = PC1, y = PC2, color = Water_quality, shape = Cropping_system)) +
  geom_point(size = 4) +  # Scatter plot with points
  theme_minimal() +       # Clean theme
  labs(title = "PCA Biplot", x = "PC1", y = "PC2") +
  geom_hline(yintercept = 0, linetype = "dashed") +  # Add dashed lines for the origin
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(legend.position = "right") +  # Legend position
  scale_color_brewer(palette = "Set1")  # Color palette for Genotypes

fviz_pca_var(pca_result, col.var = "black", repel = TRUE)  # Arrows showing variable contributions

# Load necessary libraries
library(ggplot2)
library(factoextra)
library(ggrepel)     # For text label repelling (avoids overlap)
library(dplyr)       # For data manipulation
library(gridExtra)   # To combine multiple plots if needed

#Variance explained for the first two PCs
variance_explained <- round(100 * (pca_result$sdev^2 / sum(pca_result$sdev^2)), 1)
xlab_text <- paste0("PC1 (", variance_explained[1], "% Variance)")
ylab_text <- paste0("PC2 (", variance_explained[2], "% Variance)")

# Create a publication-quality PCA biplot
pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2)) +
  
  # Scatter plot with points (color for Genotype, shape for Treatment)
  geom_point(aes(color = Water_quality, shape = Cropping_system), size = 4, alpha = 0.8) +
  
  # Add labels for samples (optional, can be omitted if too cluttered)
  geom_text_repel(aes(label = rownames(pca_scores)), 
                  size = 3, 
                  segment.color = 'grey50', 
                  max.overlaps = 10) +  # Limits overlaps
  
  # Customize the color and shape palettes
  scale_color_manual(values = c("red", "blue", "green", "purple","yellow"))+  # Add more colors if needed +  # Color for Genotype
  scale_shape_manual(values = c(16, 17, 15, 18)) +  # Custom shapes for Treatment
  
  # Add custom labels for axes with variance explained
  labs(title = "PCA Biplot of Genotypes and Treatments",
       x = xlab_text, y = ylab_text) +
  
  # Center the plot at the origin
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  
  # Customize the theme for a clean, publication-quality appearance
  theme_minimal(base_size = 15) +  # Base font size
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "right",  # Adjust legend position
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid = element_blank(),  # Remove grid lines for cleaner look
    axis.line = element_line(color = "black")  # Add axis lines
  )

# Add the biplot arrows for variable contributions (traits)
biplot_arrows <- fviz_pca_var(pca_result, 
                              col.var = "black", 
                              repel = TRUE, 
                              labelsize = 5, 
                              arrowsize = 1.2) + 
  theme_minimal(base_size = 15) +  # Same clean theme
  labs(title = "Contribution of Traits") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold")
  )


#####################################3#------------------------------------------------------------
library(plotly)

# Ensure Genotype and Treatment are factors
pca_scores$Cropping_system <- Anth_Data$Cropping_system
pca_scores$Field_capacity <- Anth_Data$Field_capacity
pca_scores$Water_quality <- Anth_Data$Water_quality

# Create a custom color scale for Genotype
genotype_colors <- c("red", "blue", "green", "purple", "yellow")
genotype_levels <- levels(pca_scores$Water_quality)

# Create a custom symbol mapping for Treatment
treatment_symbols <- c("circle", "square", "diamond", "cross", "x", "star")
treatment_levels <- levels(pca_scores$Field_capacity)

# Check if the number of levels in Treatment matches the number of symbols
if (length(treatment_levels) > length(treatment_symbols)) {
  stop("You need more symbols to match the number of Treatment levels!")
}

# Define color and symbol mapping functions
get_color <- function(genotype) {
  genotype_colors[as.numeric(Water_quality)]
}

get_symbol <- function(treatment) {
  treatment_symbols[as.numeric(Cropping_system)]
}

# Plotly 3D PCA plot with manual color and symbol mapping
pca_3d <- plot_ly(
  data = pca_scores, 
  x = ~PC1, y = ~PC2, z = ~PC3, 
  type = 'scatter3d', 
  mode = 'markers',
  marker = list(
    size = 6, 
    color = ~Water_quality,  # Use Genotype to map color
    colorscale = genotype_colors,
    symbol = ~Cropping_system,  # Use Treatment to map symbol
    symbol = sapply(pca_scores$Cropping_system, get_symbol)
  ),
  text = ~paste("Genotype:", Water_quality, "<br>Treatment:", Cropping_system),  # Hover info
  hoverinfo = "text"
) %>%
  layout(
    title = "3D PCA Biplot of Genotypes and Treatments",
    scene = list(
      xaxis = list(title = xlab_text), 
      yaxis = list(title = ylab_text), 
      zaxis = list(title = "PC3")
    ),
    legend = list(
      title = list(text = "Genotypes"),
      orientation = 'h',  # Horizontal legend
      yanchor = 'bottom',
      y = 1.1,  # Adjust position if needed
      xanchor = 'right',
      x = 1
    )
  ) %>%
  add_trace(
    type = "scatter3d",
    mode = "markers",
    x = pca_scores$PC1,
    y = pca_scores$PC2,
    z = pca_scores$PC3,
    marker = list(
      size = 6,
      color = sapply(pca_scores$Water_quality, get_color),
      symbol = sapply(pca_scores$Cropping_system, get_symbol)
    ),
    text = ~paste("Genotype:", Water_quality, "<br>Treatment:", Cropping_system),
    hoverinfo = "text"
  )

# Display the plot
pca_3d

