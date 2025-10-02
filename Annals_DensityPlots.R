# ============================
# Required libraries
# ============================
# Data manipulation (pipes, mutate, case_when, coalesce)
library(dplyr)
# Date/time utilities (time_length)
library(lubridate)
# Plotting (ggplot2), density/arrange helpers (ggpubr), unit() for margins
library(ggplot2)
library(ggpubr)
library(grid)


full <- read.csv("RelBC_data_R2.csv",, sep = ",",dec = "." ) # load UKB data from Annals_Outcome_Definition.R

### Density plots by decades
# Define age groups by decades
full$age_group <- NA
full$age_group[full$Age >= 40 & full$Age <50] <- 1
full$age_group[full$Age >= 50 & full$Age <60] <- 2
full$age_group[full$Age >= 60 & full$Age <70] <- 3
full$age_group[full$Age >= 70] <- 4

# Define the color palettes
color_palettes <- list(
  SAT_BC    = c("#6baed6", "#3182bd"),
  VAT_BC    = c("#fd8d3c", "#e6550d"),
  Muscle_BC = c("#74c476", "#31a354"),
  SMFF      = c("#fb6a4a", "#de2d26")
)

# Define the columns to loop over
columns <- c("SAT_BC", "VAT_BC", "Muscle_BC", "SMFF")

##### Density Plots for Age Groups #####

# Define the age groups
age_groups <- c(1:max(full$age_group))
# define margin
margin = c(-4.8,0,-2,0)
# open list
density_graphs = list()

# Loop over the columns and age groups
for (i in age_groups) {
  for (col in columns) {
    # Subset data for the current age group
    data_subset <- subset(full, full$age_group == i)
    
    # Create the density plot
    plot <- ggdensity(
      data_subset, x = col,
      fill = "Sex", 
      palette = color_palettes[[col]],
      alpha = 0.7
    ) +
      scale_fill_manual(values = color_palettes[[col]], 
                        labels = c("F", "M")) + # Setting legend labels
      theme_minimal() +
      theme(legend.position = "none") +
      labs(y = "", x = "") +
      theme(
        axis.title = element_text(family = "Helvetica", size = 22),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        axis.text.y = element_text(family = "Helvetica", size = 0),
        axis.text.x = element_text(family = "Helvetica", size = 0),
        plot.margin = unit(margin, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank() 
      )
    
    # Add scale limits based on column
    if (col == "SAT_BC") {
      plot <- plot + scale_x_continuous(limits = c(10, 80)) +
        scale_y_continuous(limits = c(0, 0.09))
    } else if (col == "VAT_BC") {
      plot <- plot + scale_x_continuous(limits = c(-0.5, 40)) +
        scale_y_continuous(limits = c(0, 0.14))
    } else if (col == "Muscle_BC") {
      plot <- plot + scale_x_continuous(limits = c(15, 80)) +
        scale_y_continuous(limits = c(0, 0.08))
    } else if (col == "SMFF") {
      plot <- plot + scale_x_continuous(limits = c(5, 30)) +
        scale_y_continuous(limits = c(0, 0.2))
    }
    if (i == max(full$age_group)) {
      plot <- plot + theme(axis.text.x = element_text(family = "Helvetica", size = 22))
    }
    
    # Save the plot in the list
    density_graphs[[length(density_graphs) + 1]] = plot 
  }
}

