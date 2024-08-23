# creates table of CKD (stage III-V) EIDs along with source of earliest
# indication of CKD (stage III-V)

# data inputs:
# this script uses CKD diagnoses from primary (Read codes) and secondary care (ICD10)
# as well as eGFR values from UK Biobank assessment centre + primary care to identify cases of
# CKD (stage III-V) and their earliest indication of disease.

# libraries

library(tidyverse)
library(data.table)
library(magrittr)
library(lubridate)
library(rstatix)
library(car)

work_dir = ""

read_dict_file = paste0(work_dir, "/", "inputs/all_lkps_maps_v3.txt") # dictionary of read codes
icd10_codes_file = paste0(work_dir, "/", "inputs/icd10_codes.tsv") # dictionary of ICD10 codes
icd10_blocks_file = paste0(work_dir, "/", "inputs/icd10_blocks.tsv") # dictionary of ICD10 blocks
read_2_icd10_file = paste0(work_dir, "/", "inputs/read_2_icd10.tsv") # mapping of Read 2 to ICD10
read_3_icd10_file = paste0(work_dir, "/", "inputs/read_3_icd10.tsv") # mapping of Read 3 to ICD10
ukb_data_file = "" # UKB file
gp_egfr_file = "" # UKB GP eGFR data
gp_clinic_file = "" # UKB GP clinical data
ckd_pheno_file = "" # CKD phenotypic codes for primary + secondary care

ckd_eids_file = ""

# which field codes to read in from raw UKB file, string matches for ones with multiple instances
ukb_cols_use <- c("f.eid", "f.41270", "f.41280", "f.53.0.0",
                  "f.31.0.0", "f.21003.0.0", "f.30700.0.0")

# functions

source(paste0(work_dir, "/", "functions/t2dm_time_to_ckd_func.R"))

# read in column names of UKB

ukb_cols <- fread(ukb_data_file, nrows = 0)
ukb_cols %<>%
  select(contains(ukb_cols_use)) %>%
  colnames()

# Load UKB data with columns to get

ukb_data <- fread(ukb_data_file, select = ukb_cols)

# read in ICD10 codes dictionary

icd10_codes <- fread(icd10_codes_file)

# filter for 3-letter ICD10 codes

icd10_codes %<>%
  filter(str_length(coding) == 3) %>%
  mutate(meaning = str_sub(meaning, start = 5)) # remove code at start of meaning

# read in ICD10 blocks dictionary

icd10_blocks_df <- fread(icd10_blocks_file)

# get block codes

icd10_block_codes <- icd10_blocks_df$block

# read in Read 2/3 to ICD10 mappings

read_2_icd10 <- fread(read_2_icd10_file)
read_3_icd10 <- fread(read_3_icd10_file)

# read in UKB GP clinical data

gp_clinic <- fread(gp_clinic_file)

# read in UKB GP eGFR data

gp_egfr <- fread(gp_egfr_file)

# read in read dictionary

read_dict <- fread(read_dict_file)

# read in CKD phenotype codes for primary/secondary care

ckd_pheno <- fread(ckd_pheno_file)

# filter for CKD stage 3-5 and specific diagnostic codes

ckd_pheno %<>%
  filter(!str_detect(description, "stage 2|stage 1|G1|G2|Predicted|disease monitoring|plan|annual review|H/O|History"))

# get CKD primary care codes

ckd_primary <- ckd_pheno %>%
  filter(code_type == "read")

# get CKD secondary care ICD10 codes

ckd_secondary <- ckd_pheno %>%
  filter(code_type == "icd10")

# filter GP clinical data for instances of CKD stage III-V primary care codes

gp_clinic_ckd <- gp_clinic %>%
  filter(read_2 %in% ckd_primary$code | read_3 %in% ckd_primary$code)

# format event dates

gp_clinic_ckd %<>%
  mutate(event_dt = as.Date(event_dt, format = "%d/%m/%Y"))

# get earliest indication of CKD primary care diagnosis for each person

ckd_prim_indic <- gp_clinic_ckd %>%
  # group by individual
  group_by(eid) %>%
  # filter for earliest date
  filter(event_dt == min(event_dt)) %>%
  # assign event type as 'primary_care_diagnosis'
  mutate(event_type = "primary_care_diagnosis") %>%
  select(eid, event_type, event_dt)

# check distribution of dates for weirdness

boxplot(ckd_prim_indic$event_dt)
# there is an extreme outlier

# filter for the extreme values to identify EID
ckd_prim_indic %>% filter(event_dt < "1950-08-01")
# EID 4701144 has CKD event as "1902-02-02" which is UKB code for diagnosis was same as date of birth

# what is the diagnosis?
gp_clinic_ckd %>% filter(eid == "4701144")
ckd_primary %>% filter(code == "K05..")
# chronic renal failure

# remove individual with chronic renal failure at birth; want individuals who follow T2DM -> CKD

# filter out this individual from primary care CKD indications dataframe
ckd_prim_indic %<>%
  filter(eid != "4701144")

# remove duplicate values from primary care CKD indications

ckd_prim_indic <- ckd_prim_indic[-which(duplicated(ckd_prim_indic)),]

# get secondary care ICD10 data from UKB file

icd10 <- ukb_data %>%
  select(contains("f.41270"))

# get secondary care dates of diagnoses

icd10_dates <- ukb_data %>%
  select(contains("f.41280"))

# get dates of CKD diagnoses in secondary care

icd10_l <- lapply(1:nrow(icd10), function(x) vector()) # create list of diagnoses for each individual
for(i in 1:length(icd10_l)) { # for each individual
  icd10_l[[i]] <- unlist(icd10[i,]) # add diagnoses to list
}

icd10_dates_l <- lapply(1:nrow(icd10_dates), function(x) vector()) # create list of diagnosis dates for each individual
for(i in 1:length(icd10_dates_l)) { # for each individual
  icd10_dates_l[[i]] <- unlist(icd10_dates[i,]) # add diagnosis dates to list
}

# iterate through lists of diagnoses and dates simultaneously

ckd_t <- map2(icd10_l, icd10_dates_l, function(d, t) {
  ckd_ind <- which(d %in% ckd_secondary$code) # identify indices of CKD diagnoses
  if(is_empty(ckd_ind)) {
    return(NA)
  } else {
    return(t[ckd_ind])
  } # if no CKD diagnoses, return NA, else, return date of CKD diagnosis
})

# check for multiple CKD diagnoses
sum(lengths(ckd_t) > 1) # 4329 have more than one CKD diagnosis

# get earliest time of CKD

ckd_t <- lapply(ckd_t, min)

# convert to vector

ckd_t <- unlist(ckd_t)

# convert to date format, currently in number of days from 1970-01-01

ckd_t <- as.Date(ckd_t, origin = "1970-01-01")

# create secondary care CKD diagnosis indication dataframe

ckd_second_indic <- data.frame(eid = ukb_data$f.eid, event_type = "secondary_care_diagnosis", event_dt = ckd_t)

# remove those with no secondary care CKD indication

ckd_second_indic %<>%
  filter(!is.na(event_dt))

# check distributions for weirdness
boxplot(ckd_second_indic$event_dt)
boxplot(ckd_prim_indic$event_dt)
# we would expect secondary care diagnoses of CKD to tend to occur at later points in time in the UKB population
# looks fine

# get UKB assessment centre eGFR values

# get date of baseline assessment centre attendance 

bl_dates <- ukb_data %>%
  select(f.eid, "f.53.0.0")

# get variables for calculating eGFR values

# serum creatinine, sex, age
egfr_vars <- ukb_data %>%
  select(f.eid, "f.31.0.0", "f.21003.0.0", "f.30700.0.0")

# calculate eGFR using eGFR function implementing CKD-EPI equation

egfr_vars$egfr <- egfr_calc(sc_var = "f.30700.0.0", age_var = "f.21003.0.0", sex_var = "f.31.0.0", data = egfr_vars)

# join assessment centre eGFR values with date of assessment centre attendance

egfr_vars %<>%
  full_join(bl_dates, by = "f.eid") %>%
  mutate(trait = "Glomerular_filtration_rate") %>%
  select(eid = f.eid,
         trait,
         value = egfr,
         date = f.53.0.0)

# join assessment centre eGFR values with GP eGFR values

gp_egfr <- rbind(gp_egfr, egfr_vars)

# remove NAs, those without eGFR value at baseline assessment centre

gp_egfr %<>%
  filter(!is.na(value))

# identify those individuals who have 2 sequential eGFR measurements <60 and separated by >90 days (3 months)

# get EIDs of individuals with >1 eGFR value
gt_2_val <- gp_egfr %>%
  group_by(eid) %>%
  count(trait) %>%
  filter(n > 1) %>%
  pull(eid) %>%
  unique()
  
ckd_prim_egfr_indic <- lapply(gt_2_val, function(x) {
  return(egfr_persist(gp_egfr, x))
})

# join

ckd_prim_egfr_indic <- lapply(ckd_prim_egfr_indic, function(x) {
  y <- x
  y[,3] <- as.Date(y[,3])
  return(y)
})

ckd_prim_egfr_indic <- rbindlist(ckd_prim_egfr_indic)

# remove NAs, those without 2 such eGFR values

ckd_prim_egfr_indic %<>%
  filter(!is.na(event_dt))

# combine CKD primary diagnosis, primary eGFR and secondary diagnosis indicators

ckd_indic <- rbind(as.data.frame(ckd_prim_indic), ckd_prim_egfr_indic, ckd_second_indic)

# how many cases of CKD do we have?

ckd_indic %>%
  group_by(eid) %>%
  filter(event_dt == min(event_dt)) %>%
  pull(eid) %>%
  length()
# 30715 cases

# write out table

ckd_indic %<>%
  transmute(eid = eid, date = event_dt, source = event_type)

write_tsv(ckd_indic, ckd_eids_file)

