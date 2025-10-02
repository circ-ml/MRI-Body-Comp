# ============================
# Required libraries
# ============================
# Data manipulation (pipes, mutate, case_when, coalesce)
library(dplyr)
# Date/time utilities (time_length)
library(lubridate)

# ============================
# Diabetes (ICD-10 E10–E14, ICD-9 250*)
# ============================

# Initialize target column as Date
ICD$diabetes_ICD <- as.Date(NA)

# Column name sets
cols_41270  <- sprintf("X41270.0.%d", 0:258)  # ICD-10 codes (UKB)
dates_41280 <- sprintf("X41280.0.%d", 0:258)  # dates for 41270

cols_41271  <- sprintf("X41271.0.%d", 0:46)   # ICD-9 codes
dates_41281 <- sprintf("X41281.0.%d", 0:46)   # dates for 41271

# Patterns to match (anchored)
pat_icd10 <- paste0("^E(", paste(10:14, collapse="|"), ")")  # E10–E14
pat_icd9  <- "^250"                                          # 250*

for (i in seq_len(nrow(ICD))) {
  # --- ICD-10 (E10–E14) ---
  matches10 <- vapply(cols_41270, function(cn) {
    x <- ICD[i, cn]
    !is.na(x) && grepl(pat_icd10, as.character(x))
  }, logical(1))
  
  dates10 <- if (any(matches10)) {
    as.Date(as.character(unlist(ICD[i, dates_41280[matches10]])), format = "%Y-%m-%d")
  } else as.Date(character())
  
  # --- ICD-9 (250*) ---
  matches9 <- vapply(cols_41271, function(cn) {
    x <- ICD[i, cn]
    !is.na(x) && grepl(pat_icd9, as.character(x))
  }, logical(1))
  
  dates9 <- if (any(matches9)) {
    as.Date(as.character(unlist(ICD[i, dates_41281[matches9]])), format = "%Y-%m-%d")
  } else as.Date(character())
  
  # Earliest date across both systems
  all_dates <- c(dates10, dates9)
  all_dates <- all_dates[!is.na(all_dates)]
  
  ICD$diabetes_ICD[i] <- if (length(all_dates)) min(all_dates) else as.Date(NA)
}

# ============================
# Myocardial Infarction (ICD-10 I21–I22; ICD-9 410–411)
# Stroke (ICD-10 I63; ICD-9 433–434)
# ============================

# Initialize target column as Date
ICD$myo_infarct_ICD <- as.Date(NA)
ICD$stroke_ICD      <- as.Date(NA)

cols_41270  <- sprintf("X41270.0.%d", 0:258)  # ICD-10 codes
dates_41280 <- sprintf("X41280.0.%d", 0:258)  # dates for 41270

cols_41271  <- sprintf("X41271.0.%d", 0:46)   # ICD-9 codes
dates_41281 <- sprintf("X41281.0.%d", 0:46)   # dates for 41271

## Patterns (anchored)
pat_mi_icd10   <- "^I(21|22)"     # I21–I22
pat_mi_icd9    <- "^(410|411)"    # 410–411
pat_str_icd10  <- "^I63"          # I63
pat_str_icd9   <- "^(433|434)"    # 433–434

for (i in seq_len(nrow(ICD))) {
  ## ---- Myocardial infarction ----
  mi10_matches <- vapply(cols_41270, function(cn) {
    x <- ICD[i, cn]
    !is.na(x) && grepl(pat_mi_icd10, as.character(x))
  }, logical(1))
  
  mi10_dates <- if (any(mi10_matches)) {
    as.Date(as.character(unlist(ICD[i, dates_41280[mi10_matches]])), format = "%Y-%m-%d")
  } else as.Date(character())
  
  mi9_matches <- vapply(cols_41271, function(cn) {
    x <- ICD[i, cn]
    !is.na(x) && grepl(pat_mi_icd9, as.character(x))
  }, logical(1))
  
  mi9_dates <- if (any(mi9_matches)) {
    as.Date(as.character(unlist(ICD[i, dates_41281[mi9_matches]])), format = "%Y-%m-%d")
  } else as.Date(character())
  
  mi_dates <- c(mi10_dates, mi9_dates)
  mi_dates <- mi_dates[!is.na(mi_dates)]
  ICD$myo_infarct_ICD[i] <- if (length(mi_dates)) min(mi_dates) else as.Date(NA)
  
  ## ---- Stroke ----
  st10_matches <- vapply(cols_41270, function(cn) {
    x <- ICD[i, cn]
    !is.na(x) && grepl(pat_str_icd10, as.character(x))
  }, logical(1))
  
  st10_dates <- if (any(st10_matches)) {
    as.Date(as.character(unlist(ICD[i, dates_41280[st10_matches]])), format = "%Y-%m-%d")
  } else as.Date(character())
  
  st9_matches <- vapply(cols_41271, function(cn) {
    x <- ICD[i, cn]
    !is.na(x) && grepl(pat_str_icd9, as.character(x))
  }, logical(1))
  
  st9_dates <- if (any(st9_matches)) {
    as.Date(as.character(unlist(ICD[i, dates_41281[st9_matches]])), format = "%Y-%m-%d")
  } else as.Date(character())
  
  st_dates <- c(st10_dates, st9_dates)
  st_dates <- st_dates[!is.na(st_dates)]
  ICD$stroke_ICD[i] <- if (length(st_dates)) min(st_dates) else as.Date(NA)
}

# ============================
# Hypertension (ICD-10 I10–I15; ICD-9 401–405)
# ============================

# Earliest hypertension diagnosis date from all ICD columns (ICD-10 I10–I15, ICD-9 401–405)

# Initialize target column as Date
ICD$hypertension_ICD <- as.Date(NA)

# Column name sets
cols_41270  <- sprintf("X41270.0.%d", 0:258)  # ICD-10 codes
dates_41280 <- sprintf("X41280.0.%d", 0:258)  # corresponding dates

cols_41271  <- sprintf("X41271.0.%d", 0:46)   # ICD-9 codes
dates_41281 <- sprintf("X41281.0.%d", 0:46)   # corresponding dates

# Patterns (anchored)
pat_htn_icd10 <- "^I1[0-5]"   # I10–I15 (e.g., I10, I11.0, I15.9)
pat_htn_icd9  <- "^40[1-5]"   # 401–405 (e.g., 4019)

for (i in seq_len(nrow(ICD))) {
  ## --- ICD-10 (I10–I15) ---
  m10 <- vapply(cols_41270, function(cn) {
    x <- ICD[i, cn]
    !is.na(x) && grepl(pat_htn_icd10, as.character(x))
  }, logical(1))
  
  d10 <- if (any(m10)) {
    as.Date(as.character(unlist(ICD[i, dates_41280[m10]])), format = "%Y-%m-%d")
  } else as.Date(character())
  
  ## --- ICD-9 (401–405) ---
  m9 <- vapply(cols_41271, function(cn) {
    x <- ICD[i, cn]
    !is.na(x) && grepl(pat_htn_icd9, as.character(x))
  }, logical(1))
  
  d9 <- if (any(m9)) {
    as.Date(as.character(unlist(ICD[i, dates_41281[m9]])), format = "%Y-%m-%d")
  } else as.Date(character())
  
  ## --- Assign earliest non-NA date ---
  all_dates <- c(d10, d9)
  all_dates <- all_dates[!is.na(all_dates)]
  ICD$hypertension_ICD[i] <- if (length(all_dates)) min(all_dates) else as.Date(NA)
}

# ============================
# ASCVD death definition (ICD-10 I0–I6, I70–I78)
# ============================

# Add target column next to CauseOfDeath_ICD
ICD <- add_column(ICD, ASCVD = NA, .after = "CauseOfDeath_ICD") # data field 40001

# Flag ASCVD death (1) if cause of death ICD starts with I0–I6 or any of I70–I78; else 0
ICD$ASCVD <- ifelse(
  grepl("I0", full$CauseOfDeath_ICD) |
    grepl("I1", full$CauseOfDeath_ICD) |
    grepl("I2", full$CauseOfDeath_ICD) |
    grepl("I3", full$CauseOfDeath_ICD) |
    grepl("I4", full$CauseOfDeath_ICD) |
    grepl("I5", full$CauseOfDeath_ICD) |
    grepl("I6", full$CauseOfDeath_ICD) |
    grepl("I70", full$CauseOfDeath_ICD) |
    grepl("I71", full$CauseOfDeath_ICD) |
    grepl("I72", full$CauseOfDeath_ICD) |
    grepl("I73", full$CauseOfDeath_ICD) |
    grepl("I74", full$CauseOfDeath_ICD) |
    grepl("I75", full$CauseOfDeath_ICD) |
    grepl("I76", full$CauseOfDeath_ICD) |
    grepl("I77", full$CauseOfDeath_ICD) |
    grepl("I78", full$CauseOfDeath_ICD),
  1, 0
)



ICD_pull_z <- select(ICD, PSID, myo_infarct_ICD, stroke_ICD, 
                     diabetes_ICD, hypertension_ICD, ASCVD)

write.csv(ICD_pull_z, "ICD_pull_zScore_EARLIEST.csv", row.names=FALSE)

