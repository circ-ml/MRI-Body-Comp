# ============================
# Required libraries
# ============================
library(survival)     # Cox models (coxph) and Surv()
library(ggplot2)      # Plotting (spline curves and table text)
library(splines)      # Natural cubic splines ns()
library(gridExtra)    # arrangeGrob for assembling multiple plots
library(RColorBrewer) # Color palettes
library(cowplot)      # plot_grid and draw_label to compose figures
library(grid)         # Low-level grid drawing

full <- read.csv("path/RelBC_data_R2.csv", sep = ",",dec = "." ) # load UKB data from Annals_Outcome_Definition.R

# Adjustment models
# Define RHS covariate sets as strings that will be pasted into formulas.
# model0: minimal clinical adjustment; fullmodel: adds WC and BMI_cat.
model_levels <- list(
  model0 = "Age_MRI + factor(Hypertension_history) + factor(Smoking_cats)",
  fullmodel = "Waistcircumference_2014_48 + factor(BMI_cat) + Age_MRI + factor(Hypertension_history) + factor(Smoking_cats)"
)

# Subset by sex
# Pre-split dataset to simplify downstream loops.
male <- subset(full, Sex == "M")
female <- subset(full, Sex == "F")

# Variables
# Predictors to be modeled with splines. pretty_names provides axis labels.
variables <- c("SAT_BC", "VAT_BC", "Muscle_BC", "SMFF")
pretty_names <- c(SAT_BC = "SATrel (%)", VAT_BC = "VATrel (%)", Muscle_BC = "SMrel (%)", SMFF = "SMFF (%)")

# Outcomes
# For each outcome, define time-to-event, event indicator, and sex-specific palettes.
outcomes <- list(
  Diabetes = list(
    follow_up = "FollowUp_Diab_R2",          # time variable for diabetes
    event = "Diab_inc_R2",                   # event indicator (0/1)
    palette_male = RColorBrewer::brewer.pal(9, "Blues")[c(7, 9)],   # CI and mean colors
    palette_female = RColorBrewer::brewer.pal(9, "Blues")[c(6, 8)]
  ),
  MACE = list(
    follow_up = "FollowUp_MACE_R2",          # time variable for MACE
    event = "MACE_R2",                       # event indicator (0/1)
    palette_male = RColorBrewer::brewer.pal(9, "Oranges")[c(7, 9)],
    palette_female = RColorBrewer::brewer.pal(9, "Oranges")[c(6, 8)]
  )
)

# Axis limits
# Sex-specific x-axis limits for each variable (used to crop the spline range consistently).
xlimits <- list(
  male = list(SAT_BC = c(25, 60), VAT_BC = c(5, 25), Muscle_BC = c(25, 65), SMFF = c(10, 25)),
  female = list(SAT_BC = c(40, 70), VAT_BC = c(0, 18), Muscle_BC = c(15, 45), SMFF = c(10, 25))
)

# Max HR
# Compute an upper y-limit for HR across variables to use a shared, comparable y-scale per sex/outcome/model.
# This avoids different y-scales across the 4 panels.
compute_hr_max_sex <- function(data, sex_label, outcome_info, model_text) {
  hr_vals <- c()
  for (var in variables) {
    # Build RHS: spline on var + optional adjustment covariates
    rhs <- if (model_text == "") {
      paste0("ns(", var, ", 2)")              # use 2 df natural spline if no covariates
    } else {
      paste0("ns(", var, ", 2) + ", model_text)
    }
    # Full Cox formula: Surv(time, event) ~ spline(var) + covariates
    formula <- as.formula(paste0("Surv(", outcome_info$follow_up, ", ", outcome_info$event, ") ~ ", rhs))
    fit <- try(coxph(formula, data = data), silent = TRUE)  # robust to failures
    if (inherits(fit, "try-error")) next
    
    # termplot returns partial effect (log-HR) and its SE across x grid for var
    tp <- termplot(fit, se = TRUE, plot = FALSE)[[var]]
    
    # Restrict x-range to predefined axis limits to avoid extreme tails
    lims <- xlimits[[tolower(sex_label)]][[var]]
    tp_lim <- tp[tp$x >= lims[1] & tp$x <= lims[2], ]
    
    # Upper CI on HR scale = exp(effect + SE*1) (here they use +se, not 1.96)
    # Later plots use 1.96 for shaded lines; here we just need an overall max.
    hr_vals <- c(hr_vals, exp(tp_lim$y + tp_lim$se))
  }
  max(hr_vals, na.rm = TRUE)  # global max HR to set y-axis
}

# HR table
# Build a small "table" as a ggplot with text: thresholds, aHR[95%CI], and population %
# Thresholds are chosen at 5th, 20th, 80th, 95th percentiles and summarized with < or >.
create_hr_table <- function(fit, var, data) {
  raw_cutoffs <- quantile(data[[var]], probs = c(0.05, 0.20, 0.80, 0.95), na.rm = TRUE)
  signs <- c("<", "<", ">", ">")
  rounded_cutoffs <- round(raw_cutoffs, 1)
  
  # Get smooth effect and SE over x for the fitted spline term of var
  tp <- termplot(fit, se = TRUE, plot = FALSE)[[var]]
  
  # For each cutoff, find nearest x in termplot grid and compute HR and 95% CI
  values <- sapply(rounded_cutoffs, function(cut) {
    idx <- which.min(abs(tp$x - cut))
    est <- tp$y[idx]              # log-HR (centered around reference)
    se <- tp$se[idx]              # SE of log-HR
    hr <- exp(est)                # HR at cutoff
    lower <- exp(est - 1.96 * se) # 95% CI lower
    upper <- exp(est + 1.96 * se) # 95% CI upper
    c(hr = hr, lower = lower, upper = upper)
  })
  
  # Cohort fraction below/above each cutoff (matching signs vector)
  cohort_pct <- mapply(function(cut, sign) {
    if (sign == "<") mean(data[[var]] < cut, na.rm = TRUE) * 100
    else mean(data[[var]] > cut, na.rm = TRUE) * 100
  }, rounded_cutoffs, signs)
  
  # Assemble display strings for the table
  df <- data.frame(
    Threshold = paste0(signs, " ", format(rounded_cutoffs, nsmall = 1)),
    HR_CI = paste0(format(round(values["hr", ], 2), nsmall = 2), " [",
                   format(round(values["lower", ], 2), nsmall = 2), "–",
                   format(round(values["upper", ], 2), nsmall = 2), "]"),
    Population = paste0(round(cohort_pct), "%")
  )
  
  # Header row + bold styling and manual y spacing for a neat, compact table
  df_display <- rbind(c("Threshold", "aHR [95% CI]", "% cohort"), df)
  font_bold <- c("bold", rep("plain", nrow(df_display) - 1))
  spacing_y <- rev(seq(1, by = 0.3, length.out = nrow(df_display)))
  
  # Three text columns placed by x coordinates; first row shows the variable label
  ggplot() +
    theme_void() +
    geom_text(data = data.frame(y = spacing_y, df_display, face = font_bold),
              aes(x = 1, y = y, label = ifelse(y == max(y), pretty_names[[var]], Threshold), fontface = face),
              size = 5.5, hjust = 0.5, lineheight = 0.3) +
    geom_text(data = data.frame(y = spacing_y, df_display, face = font_bold),
              aes(x = 1.9, y = y, label = HR_CI, fontface = face), size = 5.5, hjust = 0.3, lineheight = 0.5) +
    geom_text(data = data.frame(y = spacing_y, df_display, face = font_bold),
              aes(x = 3.4, y = y, label = Population, fontface = face), size = 5.5, hjust = 0.5, lineheight = 0.3) +
    xlim(0.5, 3.7) +
    ylim(0.5, max(spacing_y) + 0.5)
}

# Plot + table per variable
# Fits Cox model with spline(var) + covariates, constrains x to sex-specific limits,
# draws HR curve with 95% CI (dashed), adds HR=1 reference, and stacks with its table.
make_spline_plot_with_table <- function(data, var, outcome_info, palette, sex_label, hr_max, model_text) {
  limits <- xlimits[[tolower(sex_label)]][[var]]
  
  # Build RHS of the formula
  rhs <- if (model_text == "") {
    paste0("ns(", var, ", 2)")
  } else {
    paste0("ns(", var, ", 2) + ", model_text)
  }
  
  # Cox model
  formula <- as.formula(paste0("Surv(", outcome_info$follow_up, ", ", outcome_info$event, ") ~ ", rhs))
  fit <- try(coxph(formula, data = data), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  
  # Extract partial effects for plotting, then crop to limits
  tp <- termplot(fit, se = TRUE, plot = FALSE)[[var]]
  tp_limited <- tp[tp$x >= limits[1] & tp$x <= limits[2], ]
  
  # Spline HR plot with CI and consistent y-limit (0..hr_max)
  spline_plot <- ggplot(tp_limited, aes(x = x)) +
    geom_line(aes(y = exp(y)), color = palette[2], linewidth = 1.2) +
    geom_line(aes(y = exp(y - 1.96 * se)), linetype = "dashed", color = palette[1], linewidth = 0.8) +
    geom_line(aes(y = exp(y + 1.96 * se)), linetype = "dashed", color = palette[1], linewidth = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
    labs(x = pretty_names[[var]], y = "Hazard Ratio") +
    coord_cartesian(xlim = limits, ylim = c(0, hr_max)) +
    theme_minimal(base_size = 13) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16),
      plot.title = element_blank()
    )
  
  # Companion "table" plot built from the same fit
  table_plot <- create_hr_table(fit, var, data)
  
  # Stack spline and table vertically with custom height ratio
  plot_grid(spline_plot, table_plot, ncol = 1, rel_heights = c(2.5, 1.7))
}

# Create plots per sex/outcome/model
# Returns a named list of 4 grobs (one per variable), each with a title strip above the plot+table.
create_plots_for_group <- function(data, sex_label, outcome_name, outcome_info, model_text) {
  # Compute a common y-maximum for this sex/outcome/model across all variables
  hr_max <- compute_hr_max_sex(data, sex_label, outcome_info, model_text)
  plots <- list()
  for (var in variables) {
    # Choose palette by sex (male vs female colors for each outcome family)
    palette <- if (sex_label == "male") outcome_info$palette_male else outcome_info$palette_female
    
    # Title strip uses pretty_names; plot itself omits title to keep layout clean
    title <- ggdraw() + draw_label(pretty_names[[var]], fontface = "bold", hjust = 0.5, size = 18)
    
    # Build the spline+table grob for this variable
    spline_table <- make_spline_plot_with_table(data, var, outcome_info, palette, sex_label, hr_max, model_text)
    
    # Title stacked above plot
    plots[[var]] <- plot_grid(title, spline_table, ncol = 1, rel_heights = c(0.1, 1))
  }
  return(plots)
}

# Loop over models to create Figures 2 and 3
# For each adjustment model, sex, and outcome, create a 1x4 grid of variable panels and save as JPG.
for (model_name in names(model_levels)) {
  model_text <- model_levels[[model_name]]
  for (sex in c("male", "female")) {
    data <- if (sex == "male") male else female
    for (outcome_name in names(outcomes)) {
      outcome_info <- outcomes[[outcome_name]]
      
      # Build list of 4 panels (SAT_BC, VAT_BC, Muscle_BC, SMFF)
      plots <- create_plots_for_group(data, sex, outcome_name, outcome_info, model_text)
      
      # File base name like: Male_Diabetes_Spline_Grid_model0
      file_base <- paste0(toupper(substr(sex, 1, 1)), substr(sex, 2, nchar(sex)), "_", outcome_name,
                          "_Spline_Grid_test", model_name)
      
      # Arrange 4 panels in one row
      combined_plot <- arrangeGrob(grobs = plots, nrow = 1, ncol = 4)
      
      # Optional PDF export (currently commented out)
      #pdf(paste0(file_base, ".pdf"), width = 16, height = 5)
      #grid.draw(combined_plot)
      #dev.off()
      
      # Save high-resolution JPEG
      ggsave(paste0(file_base, ".jpg"), combined_plot, width = 16, height = 5, dpi = 300, units = "in")
    }
  }
}
