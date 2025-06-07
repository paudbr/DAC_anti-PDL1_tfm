# Author: Paula de Blas Rioja
# Description: This script processes single-cell immune metadata and generates various barplots
#              showing cell proportions across treatments and cell lines without downsampling.

library(ggplot2)
library(tidyverse)
library(dplyr)

# Load the dataset
immune_combined <- readRDS("./analysis_augie/20250120_labeled_immune.combined_edit.RDS")

# Replace sample name "KAP25L" with "KPB25L" in the metadata
immune_combined@meta.data$sample <- gsub("KAP25L", "KPB25L", immune_combined@meta.data$sample)

# Create new columns 'cell_line' and 'treatment' by splitting the sample column at the underscore
immune_combined@meta.data$cell_line <- sub("_.*", "", immune_combined@meta.data$sample)
immune_combined@meta.data$treatment <- sub(".*_", "", immune_combined@meta.data$sample)

# Extract metadata for processing
metadata <- immune_combined@meta.data

# Calculate cell counts and percentages excluding 'Cancer cells' and 'Unannotated'
cell_counts_immune <- metadata %>%
  filter(!immuno %in% c("Cancer cells", "Unannotated")) %>%  
  group_by(cell_line, treatment, immuno) %>%
  summarize(n = n(), .groups = "drop") %>%
  group_by(cell_line, treatment) %>%
  mutate(Percentage = 100 * n / sum(n)) %>%
  # Make percentages negative for the KPB25L-UV cell line for mirrored barplot effect
  mutate(
    Percentage = ifelse(cell_line == "KPB25L-UV", -Percentage, Percentage)
  )

# Set factor levels for treatment to control order in plots
cell_counts_immune$treatment <- factor(cell_counts_immune$treatment, levels = c("control", "PD-L1", "DAC", "Combination"))

# Create a barplot with dodge position showing proportions per immune cell type and treatment
p <- ggplot(cell_counts_immune, aes(x = immuno, y = Percentage, fill = treatment)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.3) +
  
  # Manual fill colors for each treatment group
  scale_fill_manual(
    values = c(
      "control" = "#E0E0E0",      # Light gray pastel
      "DAC" = "#FFB6B6",          # Light pink pastel
      "PD-L1" = "#A8C8FF",        # Light blue pastel
      "Combination" = "#C1A9E0"   # Light purple pastel
    ),
    guide = guide_legend(title = "Treatment")  # Show legend title as "Treatment"
  ) +
  
  # Add titles and axis labels
  labs(
    title = "TME Cell Proportion by Treatments in KPB25L and KPB25L-UV",
    x = "Cell Type",
    y = "Cell Proportion (%)",
    fill = "Treatment"
  ) +
  
  # Use minimal theme with custom font and styling
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "italic"),
    legend.text = element_text(size = 10),
    plot.margin = margin(10, 50, 10, 10),
    text = element_text(family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  
  # Set y-axis limits and labels (mirrored negative values for UV cell line)
  scale_y_continuous(
    limits = c(-40, 40),
    breaks = seq(-40, 40, 10), 
    labels = abs(seq(-40, 40, 10))
  ) +
  
  # Add annotations for cell lines on top and bottom of plot
  annotate("text", x = 6.5, y = 40, label = "↑KPB25L", size = 5, color = "black", family = "Times New Roman") +
  annotate("text", x = 6.75, y = -35, label = "↓KPB25L-UV", size = 5, color = "black", family = "Times New Roman")

# Save the barplot as a high-resolution PNG image
ggsave("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/Barplot_immune_cell_prop.png", plot = p, height = 6, width = 8 , dpi = 600)


# ------------------- OTHER TYPES OF PLOTS -------------------

# Re-extract metadata
metadata <- immune_combined@meta.data

# Calculate cell proportions including 'Cancer cells' and 'Unannotated' but rename them as "Tumor"
cell_counts <- metadata %>%
  mutate(
    immuno = case_when(
      immuno == "Cancer cells" ~ "Tumor", 
      immuno == "Unannotated" ~ "Tumor",
      TRUE ~ immuno
    )
  ) %>%
  group_by(sample) %>%
  dplyr::mutate(total = n()) %>%
  group_by(sample, immuno) %>%
  dplyr::mutate(total_type = n()) %>%
  mutate(Percentage = 100 * total_type / total)

# Order immune cell types (immuno) ascending by mean percentage for plotting order
immuno_levels <- cell_counts %>%
  group_by(immuno) %>%
  summarize(mean_percentage = mean(Percentage)) %>%
  arrange(mean_percentage) %>%  # Ascending order
  pull(immuno)

cell_counts <- cell_counts %>%
  mutate(immuno = factor(immuno, levels = immuno_levels))  # Set ordered factor levels

# Stacked barplot of cell proportions per sample, colored by immune cell type
ggplot(cell_counts, aes(x = sample, y = Percentage, fill = immuno)) +
  geom_bar(stat = "identity", position = "stack", color = "black", size = 0.3) +  
  scale_fill_brewer(palette = "Set3") +  
  labs(
    title = "Cell Proportion by Treatments in Different Cell Lines",
    x = "Treatment in Cell Line",
    y = "Cell Proportion (%)",
    fill = "Cell Type"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 10, face = "italic"),
    legend.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# --------------- STACKED BARPLOT EXCLUDING CANCER AND ENDOTHELIAL ---------------

# Filter out cancer and endothelial cells from metadata and calculate proportions
cell_counts_sincancer <- metadata %>%
  filter(!immuno %in% c("Cancer cells", "Unannotated", "Endothelial")) %>%  # Filter out unwanted categories
  group_by(sample, immuno) %>%
  summarize(n = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(Percentage = 100 * n / sum(n)) 

# Order immune cell types by ascending mean percentage
immuno_levels <- cell_counts_sincancer %>%
  group_by(immuno) %>%
  summarize(mean_percentage = mean(Percentage)) %>%
  arrange(mean_percentage) %>%  # Ascending order
  pull(immuno)

cell_counts_sincancer <- cell_counts_sincancer %>%
  mutate(immuno = factor(immuno, levels = immuno_levels))  # Apply ordered factor levels

# Plot stacked barplot excluding cancer and endothelial cells
ggplot(cell_counts_sincancer, aes(x = sample, y = Percentage, fill = immuno)) +
  geom_bar(stat = "identity", position = "stack", color = "black", size = 0.3) +  
  scale_fill_brewer(palette = "Set3") +  
  labs(
    title = "Cell Proportion by Treatments in Different Cell Lines",
    x = "Treatment in Cell Line",
    y = "Cell Proportion (%)",
    fill = "Cell Type"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 10, face = "italic"),
    legend.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# -------------------- DODGE BARPLOT WITH TREATMENT GROUP AND SUBTYPE --------------------


cell_counts_sincancer <- cell_counts %>%
  filter(!immuno %in% c("Tumor", "Unannotated")) %>% 
  mutate(treatment_group = ifelse(startsWith(as.character(sample), "KAP25L-UV"), "KL25UV", "KL25"),
         
         subtipo = case_when(
           grepl("_control", sample) ~ "control",
           grepl("_DAC", sample) ~ "DAC",
           grepl("_PD-L1", sample) ~ "PD-L1",
           grepl("_Combination", sample) ~ "Combination",
           TRUE ~ "Other"
         )
  )
cell_counts$treatment_group <- factor(cell_counts$treatment_group, levels = c("KL25","KL25UV"))
cell_counts$subtipo <- factor(cell_counts$subtipo, levels = c("control", "PD-L1", "DAC", "Combination"))


ggplot(cell_counts, aes(x = immuno, y = Percentage, fill = interaction(treatment_group, subtipo))) +
  geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.3) +
  scale_fill_manual(
    values = c(

      "KL25.control" = "#f6cdfb",  
      "KL25.DAC" = "#f496ec",   
      "KL25.PD-L1" = "#eb84ec",  
      "KL25.Combination" = "#e85eea", 
      
      "KL25UV.control" = "#afffcd",
      "KL25UV.DAC" = "#90e7a3",   
      "KL25UV.PD-L1" = "#7aff6e",  
      "KL25UV.Combination" = "#32ff03"  
    ),
    guide = guide_legend(title = "Treatment Group and Subtype")
  ) +
  labs(
    title = "Cell Proportion by Treatments Comparison in Different Cell Lines",
    x = "Cell Type",
    y = "Cell Proportion",
    fill = "Treatment in Cell Line"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 10, face = "italic"),
    legend.text = element_text(size = 8),
    plot.margin = margin(10, 10, 10, 10),
    text = element_text(family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

cell_counts_sincancer <- metadata %>%
  filter(!immuno %in% c("Cancer cells", "Unannotated")) %>%  
  group_by(sample, immuno) %>%
  summarize(n = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(Percentage = 100 * n / sum(n)) %>%
  mutate(
    treatment_group = ifelse(startsWith(as.character(sample), "KAP25L-UV"), "KL25UV", "KL25"),
    subtipo = case_when(
      grepl("_control", sample) ~ "control",
      grepl("_DAC", sample) ~ "DAC",
      grepl("_PD-L1", sample) ~ "PD-L1",
      grepl("_Combination", sample) ~ "Combination",
      TRUE ~ "Other"
    ))

cell_counts_sincancer$treatment_group <- factor(cell_counts_sincancer$treatment_group, levels = c("KL25", "KL25UV"))
cell_counts_sincancer$subtipo <- factor(cell_counts_sincancer$subtipo, levels = c("control", "PD-L1", "DAC", "Combination"))


cell_counts_sincancer$group_subtype <- factor(
  paste(cell_counts_sincancer$treatment_group, cell_counts_sincancer$subtipo, sep = "."),
  levels = c("KL25.control", "KL25.PD-L1", "KL25.DAC", "KL25.Combination",
             "KL25UV.control", "KL25UV.PD-L1", "KL25UV.DAC", "KL25UV.Combination")
)

ggplot(cell_counts_sincancer, aes(x = immuno, y = Percentage, fill = group_subtype)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.3) +
  scale_fill_manual(
    values = c(
    
      "KL25.control" = "#f6cdfb",  
      "KL25.PD-L1" = "#eb84ec",  
      "KL25.DAC" = "#f496ec",   
      "KL25.Combination" = "#e85eea",  
      
    
      "KL25UV.control" = "#b3e5fc",  
      "KL25UV.PD-L1" = "#4fc3f7",  
      "KL25UV.DAC" = "#0288d1",  
      "KL25UV.Combination" = "#01579b"  
    ),
    guide = guide_legend(title = "Treatment Group and Subtype")
  ) +
  labs(
    title = "Cell Proportion by Treatments Comparison in Different Cell Lines",
    x = "Cell Type",
    y = "Cell Proportion",
    fill = "Treatment in Cell Line"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 10, face = "italic"),
    legend.text = element_text(size = 8),
    plot.margin = margin(10, 10, 10, 10),
    text = element_text(family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

