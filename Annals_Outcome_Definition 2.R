# ============================
# Required libraries
# ============================
# Data manipulation (pipes, mutate, case_when, coalesce)
library(dplyr)
# Date/time utilities (time_length)
library(lubridate)

#prepare data for UKBB: body comp and relevant chronic disease

ukbb <- read.csv("ukbb_ICD.csv", sep = ",",dec = "." ) # UKB data
icd <- read.csv("ICD_pull_zScore_EARLIEST.csv", sep = ",",dec = "." ) # from Annals_ICD_extraction.R
smoke <- read.csv("ukbb_smoke.csv", sep = ",",dec = "." ) # UKB smoking data

### Merge data and keep distinct rows only
full <- merge(ukbb, icd, by='PSID')
full <- merge(full, smoke, by='PSID')

# ============================
# Date parsing: convert selected columns to Date class (YYYY-MM-DD)
# ============================
# Define as Date
full$Death_Date          <- as.Date(as.character(full$Death_Date),          format = "%Y-%m-%d") # data field 40000
full$Date_Lost_FU_191    <- as.Date(as.character(full$Date_Lost_FU_191),    format = "%Y-%m-%d") # data field 191
full$diabetes_ICD        <- as.Date(as.character(full$diabetes_ICD),        format = "%Y-%m-%d") # see Annals_ICD_extraction.R
full$StudyDate           <- as.Date(full$StudyDate,                          format = "%Y-%m-%d")# see Annals_ICD_extraction.R
full$hypertension_ICD    <- as.Date(as.character(full$hypertension_ICD),    format = "%Y-%m-%d") # see Annals_ICD_extraction.R
full$stroke_ICD.y        <- as.Date(as.character(full$stroke_ICD.y),        format = "%Y-%m-%d") # see Annals_ICD_extraction.R
full$myo_infarct_ICD.y   <- as.Date(as.character(full$myo_infarct_ICD.y),   format = "%Y-%m-%d") # see Annals_ICD_extraction.R

# ============================
# Censoring rule for major adverse cardiovascular events (MACE) analysis
# - If atherosclerotic cardiovascular disease (ASCVD) death occurs after 
# 2022-10-31 (cencoring date for ICD codes), set ASCVD death to 0 and Death_Date to NA
# ============================
# Update rows where ASCVD == 1 and Death_Date > 2022-10-31
to_update                 <- full$ASCVD == 1 & full$Death_Date > "2022-10-31" # see Annals_ICD_extraction.R 
full$ASCVD[to_update]     <- 0
full$Death_Date[to_update] <- NA

# ============================
# Follow-up (FU) time calculations (years) from StudyDate (date of MRI) to event dates
# ============================
full$FU_Death        <- time_length(full$Death_Date        - full$StudyDate, unit = "years")
full$FU_LostFU       <- time_length(full$Date_Lost_FU_191  - full$StudyDate, unit = "years")
full$FU_Diab         <- time_length(full$diabetes_ICD      - full$StudyDate, unit = "years")
full$FU_MyoInfarct   <- time_length(full$myo_infarct_ICD.y - full$StudyDate, unit = "years")
full$FU_Stroke       <- time_length(full$stroke_ICD.y      - full$StudyDate, unit = "years")
full$FU_Hypertension <- time_length(full$hypertension_ICD  - full$StudyDate, unit = "years")

# ============================
# Remove positive FU times for history variables
# (Note: as written, positive FU times are set to NA.)
# ============================
full$FU_Diab[full$FU_Diab > 0]               <- NA
full$FU_MyoInfarct[full$FU_MyoInfarct > 0]   <- NA
full$FU_Stroke[full$FU_Stroke > 0]           <- NA
full$FU_Hypertension[full$FU_Hypertension > 0] <- NA

# ============================
# Disease history flags based on FU availability
# 1 if non-missing FU indicates disease before MRI, else 0
# ============================
full$Diab_history        <- ifelse(!is.na(full$FU_Diab),         1, 0)

full$MACE_history        <- NA
full$MACE_history        <- ifelse(!is.na(full$FU_MyoInfarct) | !is.na(full$FU_Stroke), 1, 0)

full$Stroke_history      <- ifelse(!is.na(full$FU_Stroke),       1, 0)
full$MyoInfarct_history  <- ifelse(!is.na(full$FU_MyoInfarct),   1, 0)
full$Hypertension_history <- ifelse(!is.na(full$FU_Hypertension), 1, 0)

# ============================
# Recompute FU times (years) for incidence analyses
# ============================
full$FU_Diab         <- time_length(full$diabetes_ICD      - full$StudyDate, unit = "years")
full$FU_MyoInfarct   <- time_length(full$myo_infarct_ICD.y - full$StudyDate, unit = "years")
full$FU_Stroke       <- time_length(full$stroke_ICD.y      - full$StudyDate, unit = "years")
full$FU_Hypertension <- time_length(full$hypertension_ICD  - full$StudyDate, unit = "years")

# ============================
# Remove negative FU times (i.e., events occurring before MRI)
# ============================
full$FU_Diab[full$FU_Diab < 0]                 <- NA
full$FU_MyoInfarct[full$FU_MyoInfarct < 0]     <- NA
full$FU_Stroke[full$FU_Stroke < 0]             <- NA
full$FU_Hypertension[full$FU_Hypertension < 0] <- NA

# ============================
# Incident event flags (1 if non-missing FU indicates disease after MRI, else 0)
# ============================
full$Diab_inc_R2       <- ifelse(!is.na(full$FU_Diab),         1, 0)
full$Hypertension_inc  <- ifelse(!is.na(full$FU_Hypertension), 1, 0)

# ============================
# MACE composite (ASCVD death, myocardial infarction (MI), or stroke)
# ============================
full$MACE_R2 <- NA
full$MACE_R2 <- ifelse(
  full$ASCVD == 1 | !is.na(full$FU_MyoInfarct) | !is.na(full$FU_Stroke),
  1, 0
)

# ============================
# ICD-based outcomes: define "Last_Alive_ICD" as ICD censoring date
# ============================
full$Last_Alive_ICD <- NA
full <- full %>%
  mutate(
    Last_Alive_ICD = case_when(
      is.na(FU_Death) & is.na(FU_LostFU) ~ "2022-10-31",  # ICD censoring date
      TRUE ~ NA
    )
  )

# Parse Last_Alive_ICD and compute alive FU for ICD analyses
full$Last_Alive_ICD <- as.Date(as.character(full$Last_Alive_ICD), format = "%Y-%m-%d")
full$FU_alive_ICD   <- time_length(full$Last_Alive_ICD - full$StudyDate, unit = "years")

# ============================
# Incident diabetes analysis: define follow-up time
# If no death, no diabetes event, and no loss to follow-up: use FU_alive_ICD
# Otherwise: pmin(...) returns earliest of diabetes, lost to FU, or death
# ============================
full$FollowUp_Diab_R2 <- NA
full <- full %>%
  mutate(
    FollowUp_Diab_R2 = case_when(
      is.na(FU_Death) & is.na(FU_Diab) & is.na(FU_LostFU) ~ FU_alive_ICD,
      TRUE ~ pmin(FU_Diab, FU_LostFU, FU_Death, na.rm = TRUE)
    )
  )

# ============================
# MACE analysis: define follow-up time
# If no death, no MI, no stroke, and no loss to follow-up: use FU_alive_ICD
# Otherwise: pmin(...) returns earliest of death, MI, stroke, or loss to FU
# ============================
full$FollowUp_MACE_R2 <- NA
full <- full %>%
  mutate(
    FollowUp_MACE_R2 = case_when(
      is.na(FU_Death) & is.na(FU_MyoInfarct) & is.na(FU_Stroke) & is.na(FU_LostFU) ~ FU_alive_ICD,
      TRUE ~ pmin(FU_Death, FU_MyoInfarct, FU_Stroke, FU_LostFU, na.rm = TRUE)
    )
  )


full <- full %>%
  filter(is.na(Date_Lost_FU_191) | Date_Lost_FU_191 >= StudyDate)


full <- subset(full, Race != 2) # only keep non-Whites
full <- subset(full, MyoInfarct_history != 1) # exclude those with history of MI
full <- subset(full, Stroke_history != 1) # exclude those with history of stroke
full <- subset(full, Diab_history != 1) # exclude those with history of diabetes

#ascvd_date <- subset(full, ASCVD == 1, select = Death_Date)
ascvd_date_new <- subset(full, ASCVD == 1, select = Death_Date)


## index body comp by body comp
full$SAT_BC <- full$SAT/(full$SAT+full$VAT+full$Muscle)*100
full$VAT_BC<- full$VAT/(full$SAT+full$VAT+full$Muscle)*100
full$Muscle_BC <- full$Muscle/(full$SAT+full$VAT+full$Muscle)*100


## define BMI categories
full$BMI_cat <- NA
full$BMI_cat[full$BMI < 25] <- 0
full$BMI_cat[full$BMI >= 25 & full$BMI < 30] <- 1
full$BMI_cat[full$BMI >= 30] <- 2




write.csv(full, "RelBC_data_R2.csv", row.names=FALSE)

