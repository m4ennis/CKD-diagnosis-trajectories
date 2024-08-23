# Modeling of effect of comorbidities between T2DM -> CKD on progression to CKD

# libraries

library(tidyverse)
library(data.table)
library(magrittr)
library(lubridate)
library(UpSetR)
library(survival)
library(survminer)
library(GGally)
library(ComplexHeatmap)
library(FactoMineR)
library(factoextra)
library(rstatix)
library(car)

work_dir = ""

read_dict_file = paste0(work_dir, "/", "inputs/all_lkps_maps_v3.txt") # dictionary of read codes
phecode_icd10_file = paste0(work_dir, "/", "inputs/phecode_icd10.csv") # mapping of ICD10 bodes to phecodes
read_2_icd10_file = paste0(work_dir, "/", "inputs/read_2_icd10.tsv") # mapping of Read 2 to ICD10
read_3_icd10_file = paste0(work_dir, "/", "inputs/read_3_icd10.tsv") # mapping of Read 3 to ICD10
ukb_data_file = "" # UKB file
gp_egfr_file = "" # UKB GP eGFR data
gp_clinic_file = "" # UKB GP clinical data
ckd_pheno_file = "" # CKD phenotypic codes for primary + secondary care
t2dm_pheno_file = "" # T2DM phenotypic codes for primary + secondary care
primary_hba1c_file = "" # UKB GP HbA1c data
t2dm_eids_file = paste0(work_dir, "/", "t2dm_criteria_eids.tsv") # UKB EIDs of individuals passing T2DM inclusion/exclusion criteria
sig_dis_file = paste0(work_dir, "/", "intermediates/t2dm_enrich_diseases.tsv") # Significantly enriched diseases in T2DM
plot_dir = paste0(work_dir, "/", "plots/") # set plot output directory

# which field codes to read in from raw UKB file, string matches for ones with multiple instances
ukb_cols_use <- c("f.eid", "f.20002", "f.20008", "f.41270", "f.41280", "f.53.0.0",
                  "f.31.0.0", "f.21003.0.0", "f.30700.0.0", "f.30750.0.0",
                  "f.40000.0.0", "f.21001")

# functions

source(paste0(work_dir, "/", "functions/t2dm_time_to_ckd_func.R"))

# read in column names of UKB

ukb_cols <- fread(ukb_data_file, nrows = 0)
ukb_cols %<>%
  select(contains(ukb_cols_use)) %>%
  colnames()

# Load UKB data with columns to get

ukb_data <- fread(ukb_data_file, select = ukb_cols)

# read in Read 2/3 to ICD10 mappings

read_2_icd10 <- fread(read_2_icd10_file)
read_3_icd10 <- fread(read_3_icd10_file)

# read in ICD10-phecode mapping

phecode_icd10 <- fread(phecode_icd10_file)

# get 3 letter ICD10 codes

phecode_icd10 %<>%
  filter(str_length(ICD10) == 3)

# read in UKB GP clinical data

gp_clinic <- fread(gp_clinic_file)

# read in UKB GP eGFR data

gp_egfr <- fread(gp_egfr_file)

# read in read dictionary

read_dict <- fread(read_dict_file)

# read in T2DM EID inclusion/exclusion criteria file

t2dm_eids <- fread(t2dm_eids_file)

# filter for those who pass inclusion/exclusion criteria

t2dm_eids %<>%
  filter(t2dm_final == 1)

# filter UKB data for T2DM EIDs

ukb_data %<>%
  filter(f.eid %in% t2dm_eids$eid)

# filter UKB GP clinical for T2DM EIDs

gp_clinic %<>%
  filter(eid %in% t2dm_eids$eid)

# filter UKB GP clinical for T2DM EIDs

gp_egfr %<>%
  filter(eid %in% t2dm_eids$eid)

# get self-reported non-cancer diagnoses from UKB file

self_diag <- ukb_data %>%
  select(contains("f.20002"))

# get self-reported time of non-cancer diagnosis from UKB file

self_dates <- ukb_data %>%
  select(contains("f.20008"))

# read in CKD phenotype codes for primary/secondary care

ckd_pheno <- fread(ckd_pheno_file)

# read in T2DM phenotype codes for primary/secondary care

t2dm_pheno <- fread(t2dm_pheno_file)

# get T2DM primary care codes

t2dm_primary <- t2dm_pheno %>%
  filter(code_type %in% c("read_2", "read_3"))

# get T2DM secondary care ICD10 codes

t2dm_secondary <- t2dm_pheno %>%
  filter(code_type == "ICD10")

# read in primary care HbA1c data

primary_hba1c <- fread(primary_hba1c_file)

# filter primary care HbA1c data for T2DM EIDs

primary_hba1c %<>%
  filter(eid %in% t2dm_eids$eid)

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
sum(lengths(ckd_t) > 1) # 1336 have more than one CKD diagnosis

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

ckd_prim_egfr_indic <- lapply(unique(gp_egfr$eid), function(x) {
  return(egfr_persist(gp_egfr, x))
})

# join

ckd_prim_egfr_indic <- reduce(ckd_prim_egfr_indic, rbind)

# remove NAs, those without 2 such eGFR values

ckd_prim_egfr_indic %<>%
  filter(!is.na(event_dt))

# convert to date

ckd_prim_egfr_indic %<>%
  mutate(event_dt = as.Date(event_dt, "1970-01-01"))

# combine CKD primary diagnosis, primary eGFR and secondary diagnosis indicators

ckd_indic <- rbind(as.data.frame(ckd_prim_indic), ckd_prim_egfr_indic, ckd_second_indic)

# how many cases of CKD do we have?

ckd_indic %>%
  group_by(eid) %>%
  filter(event_dt == min(event_dt)) %>%
  pull(eid) %>%
  length()
# 6949 cases

# get the earliest indicator for each individual and plot the counts of each indicator type

ckd_indic %>%
  # group by individual
  group_by(eid) %>%
  # filter for earliest CKD event for each individual
  filter(event_dt == min(event_dt)) %>%
  # plot
  ggplot(aes(x = event_type)) +
  geom_bar(col = "black", fill = "darkgrey", width = 0.5) +
  theme_bw() +
  xlab("First indicator of CKD stage 3-5") +
  ylab("Number of individuals") +
  scale_x_discrete(labels = c(primary_care_diagnosis = "Primary care diagnosis\nof CKD stage 3-5",
                              primary_egfr_persist_lt_60 = "Two primary care/assessment centre\neGFR values <60 separated by >90 days",
                              secondary_care_diagnosis = "Secondary care diagnosis\nof CKD stage 3-5")) +
  scale_y_continuous(breaks = seq(0, 3000, 500))
# largest is secondary care, makes sense considering
# that we don't have access to ~50% of GP data in UKB


# upset plot of overlap between CKD indicator types

# create sets of EIDs within each event type
ckd_event_map <- data.frame(id = c("primary_care_diagnosis", "primary_egfr_persist_lt_60", "secondary_care_diagnosis"),
                            name = c("Primary care CKD\nstage 3-5","Persistent eGFR <60", "Secondary care CKD\nstage 3-5 diagnosis"))
ckd_upset_in <- lapply(unique(ckd_indic$event_type), function(x) ckd_indic %>% filter(event_type == x) %>% pull(eid))
names(ckd_upset_in) <- ckd_event_map$name[match(unique(ckd_indic$event_type), ckd_event_map$id)]
ckd_upset_in <- fromList(ckd_upset_in)

# plot
upset(ckd_upset_in)
# largest group is those with only secondary care diagnosis
# again, makes sense given the missing ~50% GP data in entire UKB


# filter UKB GP clinical for T2DM primary care codes 

gp_clinic_t2dm <- gp_clinic %>%
  filter(read_2 %in% t2dm_primary$read_code | read_3 %in% t2dm_primary$read_code)

# format dates 

gp_clinic_t2dm %<>%
  mutate(event_dt = as.Date(event_dt, format = "%d/%m/%Y"))

# get earliest primary care code for T2DM for each individual

t2dm_prim_indic <- gp_clinic_t2dm %>%
  # group by individual
  group_by(eid) %>%
  # filter for earliest T2DM diagnosis for each individual
  filter(event_dt == min(event_dt)) %>%
  # add event type
  mutate(event_type = "primary_care_diagnosis") %>%
  select(eid, event_type, event_dt)

# check T2DM primary care dates distribution for weirdness

boxplot(t2dm_prim_indic$event_dt)
# outliers present

# check outliers
t2dm_prim_indic %>%
  filter(event_dt < "1950-01-01")
# 7 individuals with T2DM date as 1902-02-02, which is UKB code for date same as person's date of birth
# don't want "T2DM" individuals diagnosed at birth


# remove those 7 individuals with 1902-02-02 as date of T2DM diagnoses

t2dm_prim_indic %<>%
  filter(event_dt != "1902-02-02")

# remove duplicates

t2dm_prim_indic <- t2dm_prim_indic[-which(duplicated(t2dm_prim_indic)),]

# iterate through lists of secondary care diagnoses and dates simultaneously

t2dm_t <- map2(icd10_l, icd10_dates_l, function(d, t) {
  t2dm_ind <- which(d %in% t2dm_secondary$read_code) # get indices of T2DM ICD10 codes
  if(is_empty(t2dm_ind)) {
    return(NA)
  } else {
    return(t[t2dm_ind])
  } # if no diagnosis of T2DM, return NA, else return date of T2DM ICD10 diagnosis
})

# check for multiple diagnoses of T2DM
sum(lengths(t2dm_t) > 1) # 3028 with multiple T2DM diagnoses

# get earliest T2DM diagnosis

t2dm_t <- lapply(t2dm_t, min)

# convert to vector
t2dm_t <- unlist(t2dm_t)

# convert to date format, currently as number of days since 1970-01-01

t2dm_t <- as.Date(t2dm_t, origin = "1970-01-01")

# create T2DM secondary care indication dataframe

t2dm_second_indic <- data.frame(eid = ukb_data$f.eid, event_type = "secondary_care_diagnosis", event_dt = t2dm_t)

# remove NAs, those without secondary care diagnosis
t2dm_second_indic %<>%
  filter(!is.na(event_dt))

# self-reported T2DM indication

# identify instances of self-reported T2DM, coded as 1223 in UKB field 6
self_t2dm <- t(apply(self_diag, 1, function(x) ifelse(x == "1223", T, F)))
# identify instances of self-reported Diabetes, coded as 1220 in UKB field 6
self_diab <- t(apply(self_diag, 1, function(x) ifelse(x == "1220", T, F)))

# get date of self-reported T2DM

self_t2dm_l <- lapply(1:nrow(self_t2dm), function(x) vector()) # create lists of booleans for self-reported T2DM for each individual
for(i in 1:length(self_t2dm_l)) {
  self_t2dm_l[[i]] <- unlist(self_t2dm[i,])
}

self_t2dm_dates_l <- lapply(1:nrow(self_dates), function(x) vector()) # create lists of self-reported dates for each individual
for(i in 1:length(self_t2dm_dates_l)) {
  self_t2dm_dates_l[[i]] <- unlist(self_dates[i,])
}

# iterate through boolean lists and self-reported dates simultaneously

self_t2dm_t <- map2(self_t2dm_l, self_t2dm_dates_l, function(d, t) {
  self_ind <- which(d) # indices of T2DM
  if(is_empty(self_ind)) {
    return(NA)
  } else {
    return(t[self_ind])
  } # if no self-reported T2DM, return NA, else return self-reported date of T2DM
})

# check for multiple self-reported T2DM

sum(lengths(self_t2dm_t) > 1) # 89 with multiple self-reported T2Dm

# get earliest self-reported T2DM

self_t2dm_t <- lapply(self_t2dm_t, min)

# convert to vector
self_t2dm_t <- unlist(self_t2dm_t)

# convert self-reported dates to date format
# eg. 2007.0 to 2007-01-01
# uses float_date_conv() function

self_t2dm_t <- as.Date(sapply(self_t2dm_t, float_date_conv), origin = "1970-01-01")

# create self-reported T2DM indication data.frame

self_t2dm_indic <- data.frame(eid = ukb_data$f.eid, event_type = "self_reported_t2dm", event_dt = self_t2dm_t)

# remove NAs, those without self-reported T2DM
self_t2dm_indic %<>%
  filter(!is.na(event_dt))

# get date of self-reported diabetes

self_diab_l <- lapply(1:nrow(self_diab), function(x) vector()) # create list of booleans indicating self-reported diabetes
for(i in 1:length(self_diab_l)) {
  self_diab_l[[i]] <- unlist(self_diab[i,])
}

self_diab_dates_l <- lapply(1:nrow(self_dates), function(x) vector()) # # create lists of self-reported dates for each individual
for(i in 1:length(self_diab_dates_l)) {
  self_diab_dates_l[[i]] <- unlist(self_dates[i,])
}

# iterate through boolean and dates simultaneously

self_diab_t <- map2(self_diab_l, self_diab_dates_l, function(d, t) {
  self_ind <- which(d) # get indices of self-reported diabetes 
  if(is_empty(self_ind)) {
    return(NA)
  } else {
    return(t[self_ind])
  } # if no self-reported diabetes, return NA, else return self-reported date of diabetes
})

# check for multiple self-reported diabetes

self_diab_t <- lapply(self_diab_t, min)

# convert to vector
self_diab_t <- unlist(self_diab_t)

# convert self-reported dates to date format
# eg. 2007.0 to 2007-01-01
# uses float_date_conv() function
self_diab_t <- as.Date(sapply(self_diab_t, float_date_conv), origin = "1970-01-01")

# create self-reported diabetes indication dataframe

self_diab_indic <- data.frame(eid = ukb_data$f.eid, event_type = "self_reported_diabetes", event_dt = self_diab_t)

# remove NAs, those with no self-reported diabetes
self_diab_indic %<>%
  filter(!is.na(event_dt))


# get baseline assessment centre HbA1cs values from UKB data

ukb_hba1c <- ukb_data %>%
  select(f.eid, contains("f.30750.0.0")) %>%
  full_join(bl_dates, by = "f.eid")

# format column names

ukb_hba1c %<>%
  mutate(trait = "HbA1c") %>%
  select(eid = f.eid,
         trait,
         value = f.30750.0.0,
         date = f.53.0.0)

# combine baseline asssessment centre HbA1c values with primary care HbA1c values

primary_hba1c <- rbind(primary_hba1c, ukb_hba1c)

# do any not have a date?
primary_hba1c %>%
  filter(is.na(date))
# no

# filter for individuals with at least one HbA1c >= 48 and
# get earliest instance of HbA1c value >= 48 for each individual

t2dm_prim_hba1c_indic <- primary_hba1c %>%
  # filter for HbA1c >= 48
  filter(value >= 48) %>%
  # group by individual
  group_by(eid) %>%
  # filter for earliest date of HbA1c >= 48
  filter(date == min(date, na.rm = T)) %>%
  select(eid,
         event_type = trait,
         event_dt = date)

# combine T2DM primary care diagnosis, primary care HbA1c, secondary care diagnosis,
# self-reported T2DM diagnosis, and self-reported diabetes diagnosis indications

t2dm_indic <- rbind(t2dm_prim_indic, t2dm_prim_hba1c_indic %>% mutate(event_dt = as.Date(event_dt)),
                    t2dm_second_indic, self_t2dm_indic, self_diab_indic)

# get number of T2DM cases

t2dm_indic %>%
  group_by(eid) %>%
  filter(event_dt == min(event_dt)) %>%
  pull(eid) %>%
  length()
# 43224

# identify earliest indicator of T2DM for each individual and 
# plot counts of the T2DM indicator type

t2dm_indic %>%
  # group by individual
  group_by(eid) %>%
  # filter for earliest indicator of T2DM
  filter(event_dt == min(event_dt)) %>%
  # plot
  ggplot(aes(x = event_type)) +
  geom_bar(col = "black", fill = "darkgrey", width = 0.5) +
  theme_bw() +
  xlab("First indicator of T2DM") +
  ylab("Number of individuals") +
  scale_x_discrete(labels = c(primary_care_diagnosis = "Primary care\ndiagnosis of T2DM",
                              primary_egfr_persist_lt_60 = "Primary care or\nassessment centre HbA1c\nmeasurement > 48 mmol/mol",
                              secondary_care_diagnosis = "Secondary care\ndiagnosis of T2DM",
                              self_reported_t2dm = "Self-reported T2DM\nat assessment centre",
                              self_reported_diabetes = "Self-reported Diabetes\nat assessment centre")) +
  scale_y_continuous(breaks = seq(0, 14000, 2000))
# relatively even spread, top 2 are self-reported diabetes and secondary care diagnosis of T2DM


# upset plot of overlap of T2DM indicators

# create T2DM event ID to event visualisation name map

t2dm_event_map <- data.frame(id = c("self_reported_t2dm", "primary_care_diagnosis", "self_reported_diabetes", 
                                    "HbA1c", "secondary_care_diagnosis"),
                             name = c("Self-reported T2DM", "Primary care T2DM diagnosis", "Self-reported Diabetes",
                                      "HbA1c >= 48 mmol/mol", "Secondary care T2DM diagnosis"))


# create list of EID sets for each T2DM indicator
t2dm_upset_in <- lapply(unique(t2dm_indic$event_type), function(x) t2dm_indic %>% filter(event_type == x) %>% pull(eid))
names(t2dm_upset_in) <- t2dm_event_map$name[match(unique(t2dm_indic$event_type), t2dm_event_map$id)]

t2dm_upset_in <- fromList(t2dm_upset_in)

# plot upset
upset(t2dm_upset_in)
# quite complicated, but I think reliance
# on T2DM secondary care diagnosis would be pretty
# unreliable in the subset without GP data


# combine T2DM and CKD indications data.frames

# re-set column names
colnames(t2dm_indic) <- c("eid", "t2dm_event_type", "t2dm_event_dt")
colnames(ckd_indic) <- c("eid", "ckd_event_type", "ckd_event_dt")

# join dataframes
t2dm_ckd_indic <- t2dm_indic %>%
  full_join(ckd_indic, by = "eid")

# get earliest indication of T2DM and CKD for those with both

t2dm_ckd_indic %<>%
  # remove those without CKD event
  filter(!is.na(ckd_event_type)) %>%
  # group by individual
  group_by(eid) %>%
  # filter for earliest indication of T2DM and CKD
  filter(t2dm_event_dt == min(t2dm_event_dt),
         ckd_event_dt == min(ckd_event_dt))

# plot the distributions of time to CKD from T2DM for different T2DM/CKD earliest indications
t2dm_ckd_indic %>%
  mutate(event_diff = ckd_event_dt - t2dm_event_dt) %>%
  ggplot(aes(x = event_diff/365)) +
  geom_histogram(col = "black", fill ="darkgrey") +
  theme_bw() +
  facet_grid(ckd_event_type ~ t2dm_event_type) +
  xlab("Time of CKD indication in years since T2DM indication") +
  ylab("Number of individuals") +
  scale_x_continuous(breaks = seq(-40, 40, 20)) +
  scale_y_continuous(breaks = seq(0, 600, 100)) +
  geom_segment(xend = 0, x = 0, y = 0, yend = 250, lty = 2, inherit.aes = F) +
  annotate(x = 0, y = 350, geom = "text", label = "Time of T2DM\ndiagnosis")
# individuals where secondary care diagnoses of 
# T2DM are their earliest indicators of T2DM 
# look problematic; too many CKD diagnoses before T2DM
# seems pretty good evidence should use the subset of
# UKB who has GP data

# map read codes to ICD10 codes

read_icd10 <- read_2_icd10 %>%
  transmute(read_2_code = coding, meaning = meaning) %>%
  full_join(read_3_icd10, by = "meaning") %>%
  transmute(read_2_code = read_2_code, read_3_code = coding, icd10_code = meaning)

# selecting for diseases with >= 10% prevalence within T2DM + CKD

# convert read codes to icd10

gp_clinic_block <- gp_clinic %>%
  # filter for T2DM + CKD subgroup
  filter(eid %in% ckd_indic$eid) %>%
  # convert read diagnostic codes to ICD10 diagnostic codes
  mutate(icd10_2 = read_icd10$icd10_code[match(read_2, read_icd10$read_2_code)],
         icd10_3 = read_icd10$icd10_code[match(read_3, read_icd10$read_3_code)]) %>%
  # select for physiological ICD10 blocks
  filter(str_detect(icd10_2, "A|B|C|D|E|F|G|H|I|J|K|L|M|N") | str_detect(icd10_3, "A|B|C|D|E|F|G|H|I|J|K|L|M|N")) %>%
  select(eid, event_dt, icd10_2, icd10_3)

# convert icd10 to phecode

gp_clinic_block %<>%
  mutate(phecode_2 = phecode_icd10$PheCode[match(icd10_2, phecode_icd10$ICD10)],
         phecode_3 = phecode_icd10$PheCode[match(icd10_3, phecode_icd10$ICD10)])

# get prevalence of ICD10 blocks within GP data

gp_block_prev <- gp_clinic_block %>%
  select(-event_dt, -icd10_2, -icd10_3) %>%
  pivot_longer(-eid, names_to = "read_type", values_to = "block") %>%
  group_by(eid, block) %>%
  select(-read_type) %>%
  unique.data.frame() %>%
  filter(!is.na(block)) %>%
  pull(block) %>%
  table() %>%
  as.data.frame()

# get instances of ICD10 blocks for each individual in GP records

gp_block_filt <- gp_clinic_block %>%
  select(-event_dt, -icd10_2, -icd10_3) %>%
  pivot_longer(-eid, names_to = "read_type", values_to = "block") %>%
  group_by(eid, block) %>%
  select(-read_type) %>%
  unique.data.frame() %>%
  filter(!is.na(block))


# most prevalent ICD10 secondary care codes

# filter ICD10 data for T2DM + CKD subgroup
icd10_filt <- as.matrix(icd10[ukb_data$f.eid %in% t2dm_ckd_indic$eid,])

# convert to 3 letter code

icd10_filt <- str_trunc(icd10_filt, 3, ellipsis = "")

# get instances of ICD10 blocks in secondary care for each individual
icd10_filt %<>%
  as.data.frame() %>%
  # filter for T2DM + CKD subset
  mutate(eid = ukb_data$f.eid[ukb_data$f.eid %in% t2dm_ckd_indic$eid]) %>%
  # filter for individuals who have GP data
  filter(eid %in% unique(gp_clinic$eid)[unique(gp_clinic$eid) %in% gp_clinic_block$eid]) %>%
  pivot_longer(-eid, names_to = "diag", values_to = "block") %>%
  select(eid, block) %>%
  unique.data.frame()

# filter for physiological ICD10

icd10_filt %<>%
  filter(str_detect(block, "A|B|C|D|E|F|G|H|I|J|K|L|M|N"))

# convert ICD10 codes to phecodes

icd10_filt %<>%
  mutate(block = phecode_icd10$PheCode[match(block, phecode_icd10$ICD10)]) %>%
  select(eid, block)

# combine primary + secondary care ICD10 block instances

diag_prev <- rbind(icd10_filt, gp_block_filt) %>%
  # remove duplicate phecodes
  unique.data.frame()

# get counts
diag_prev %<>%
  pull(block) %>%
  table() %>%
  as.data.frame()

colnames(diag_prev) <- c("block", "n")

# calculate frequencies
diag_prev %<>%
  arrange(desc(n)) %>%
  mutate(freq = n/length(unique(icd10_filt$eid)))

# filter for those diseases which are > 10% prevalence rate within T2DM + CKD subset

diag_prev_filt <- diag_prev %>%
  arrange(desc(freq)) %>%
  filter(freq > 0.10)

# plot distribution of CKD times relative to T2DM faceted by source of earliest T2DM/CKD indication

# get those with GP data
eid_use <- (unique(gp_clinic$eid)[unique(gp_clinic$eid) %in% gp_clinic_block$eid])

plot <- t2dm_ckd_indic %>%
  # filter for those with GP data
  filter(eid %in% eid_use) %>%
  # get time difference between CKD event date and T2DM event date
  mutate(event_diff = ckd_event_dt - t2dm_event_dt) %>%
  # convert to years
  ggplot(aes(x = event_diff/365)) +
  geom_histogram(col = "black", fill ="darkgrey") +
  theme_bw() +
  facet_grid(ckd_event_type ~ t2dm_event_type, labeller = labeller(ckd_event_type = c(primary_care_diagnosis = "Primary care\ndiagnosis of CKD", 
                                                                                      primary_egfr_persist_lt_60 = "2 eGFR measurements\n<60 separated by\n>90 days",
                                                                                      secondary_care_diagnosis = "Secondary care\ndiagnosis of CKD"),
                                                                   t2dm_event_type = c(HbA1c = "HbA1c measurement\n> 48 mmol/mol",
                                                                                       primary_care_diagnosis = "Primary care\ndiagnosis of T2DM",
                                                                                       secondary_care_diagnosis = "Secondary care\ndiagnosis of T2DM",
                                                                                       self_reported_diabetes = "Self-reported\ndiagnosis of diabetes",
                                                                                       self_reported_t2dm = "Self-reported\ndiagnosis of T2DM"))) +
  xlab("Time of CKD indication in years since T2DM indication") +
  ylab("Number of individuals") +
  scale_x_continuous(breaks = seq(-40, 60, 20)) +
  scale_y_continuous(breaks = seq(0, 300, 50)) +
  geom_segment(xend = 0, x = 0, y = 0, yend = 225, lty = 2, inherit.aes = F) +
  annotate(x = 0, y = 300, geom = "text", label = "Time of T2DM\nindication") +
  coord_cartesian(ylim = c(0, 350))
ggsave(filename = paste0(plot_dir, "t2dm_ckd_prim_care_filt_timeline_indicators.png"), plot = plot,
       units = "px", width = 2400, height = 1600)
# secondary care issue with T2DM occurring after primary care sources of CKD is significantly 
# smaller in this subset, so focusing on GP subset has alleviated this issue a good deal
# still the issue of a number of individuals having CKD before indications before their T2DM indication
# in the subset whose earliest T2DM indication is an HbA1c >48; could be non-T2DM dependent
# CKD that is developing and clinicians measuring HbA1c as a check on cause of CKD

# get diseases which we will look at based on enrichment in T2DM + >10% prevalence in T2DM + CKD

dis_use <- diag_prev_filt %>%
  pull(block)

# create T2DM block data

# convert read codes to blocks

gp_clinic_block <- gp_clinic %>%
  # filter for T2DM individuals
  filter(eid %in% t2dm_indic$eid) %>%
  # convert read codes to ICD10 codes
  mutate(icd10_2 = read_icd10$icd10_code[match(read_2, read_icd10$read_2_code)],
         icd10_3 = read_icd10$icd10_code[match(read_3, read_icd10$read_3_code)]) %>%
  # select for physiological codes
  filter(str_detect(icd10_2, "A|B|C|D|E|F|G|H|I|J|K|L|M|N") | str_detect(icd10_3, "A|B|C|D|E|F|G|H|I|J|K|L|M|N")) %>%
  select(eid, event_dt, icd10_2, icd10_3)

# convert ICD10 to phecodes

gp_clinic_block %<>%
  mutate(phecode_2 = phecode_icd10$PheCode[match(icd10_2, phecode_icd10$ICD10)],
         phecode_3 = phecode_icd10$PheCode[match(icd10_3, phecode_icd10$ICD10)])

# get instances of blocks for each individual

gp_block_filt <- gp_clinic_block %>%
  select(-event_dt, -icd10_2, -icd10_3) %>%
  pivot_longer(-eid, names_to = "read_type", values_to = "block") %>%
  group_by(eid, block) %>%
  select(-read_type) %>%
  unique.data.frame() %>%
  filter(!is.na(block))

# filter ICD10 data for T2DM individuals
# NOTE: this excludes individuals who meet our T2DM inclusion criteria but don't have any temporal data on
# their T2DM; eg. someone who reported being on T2DM-specific medication but has no diagnosis or self-report

icd10_filt <- as.matrix(icd10[ukb_data$f.eid %in% t2dm_indic$eid,])

# convert to 3 letter code

icd10_filt <- str_trunc(icd10_filt, 3, ellipsis = "")

# create df

icd10_filt %<>%
  as.data.frame() %>%
  # add eids
  mutate(eid = ukb_data$f.eid[ukb_data$f.eid %in% t2dm_indic$eid]) %>%
  # select for those with GP data
  filter(eid %in% unique(gp_clinic$eid)[unique(gp_clinic$eid) %in% gp_clinic_block$eid]) %>%
  pivot_longer(-eid, names_to = "diag", values_to = "block")

# get diagnosis order
icd10_filt$diag <- sapply(str_split(icd10_filt$diag, "\\."), function(x) tail(x, 1))

# get dates of diagnoses

icd10_t_filt <- as.matrix(icd10_dates[ukb_data$f.eid %in% t2dm_indic$eid,])

# create df

icd10_t_filt %<>%
  as.data.frame() %>%
  # add EIDs
  mutate(eid = ukb_data$f.eid[ukb_data$f.eid %in% t2dm_indic$eid]) %>%
  # filter for those with GP data
  filter(eid %in% unique(gp_clinic$eid)[unique(gp_clinic$eid) %in% gp_clinic_block$eid]) %>%
  pivot_longer(-eid, names_to = "diag", values_to = "event_dt")

# get diagnosis order

icd10_t_filt$diag <- sapply(str_split(icd10_t_filt$diag, "\\."), function(x) tail(x, 1))

# combine diagnoses and diagnosis dates

icd10_filt <- data.frame(icd10_filt, event_dt = icd10_t_filt$event_dt)

# filter out NAs
icd10_filt %<>%
  filter(!is.na(block) & !is.na(event_dt))

# filter for physiological codes

icd10_filt %<>%
  filter(str_detect(block, "A|B|C|D|E|F|G|H|I|J|K|L|M|N"))

# convert ICD10 codes to phecodes

icd10_filt %<>%
  mutate(block = phecode_icd10$PheCode[match(block, phecode_icd10$ICD10)]) %>%
  select(eid, event_dt, block)

# get gp blocks

gp_clinic_block %<>%
  select(-icd10_2, -icd10_3) %>%
  pivot_longer(-c(eid, event_dt), names_to = "read_type", values_to = "block") %>%
  group_by(eid, block) %>%
  # get earliest date of block occurrence for each individual
  filter(event_dt == min(event_dt)) %>%
  select(eid, event_dt, block) %>%
  unique.data.frame()

# format dates 

gp_clinic_block %<>%
  mutate(event_dt = as.Date(event_dt, format = "%d/%m/%Y"))

icd10_filt %<>%
  mutate(event_dt = as.Date(event_dt))

# join primary care and secondary care blocks

diags <- rbind(icd10_filt, gp_clinic_block)

# filter for physiological blocks

diags %<>% filter(block %in% dis_use) %>%
  mutate(eid = as.character(eid)) %>%
  # remove redundant diabetes + CKD codes and unknown/uncertain/same as date of birth diagnoses
  filter(block != "250.2" & block != "585.3" & event_dt != "1901-01-01" & event_dt != "1902-02-02" & event_dt != "1903-03-03" & !is.na(event_dt)) %>%
  group_by(eid, block) %>%
  # filter for earliest date of block across primary + secondary care diagnoses for each individual
  filter(event_dt == min(event_dt, na.rm = T))

# remove redundant rows

diags <- unique.data.frame(diags)

# get earliest indication times of T2DM and CKD

t2dm_indic_min <- t2dm_indic %>%
  mutate(eid = as.character(eid)) %>%
  group_by(eid) %>%
  # get earliest T2DM event date
  filter(t2dm_event_dt == min(t2dm_event_dt, na.rm = T)) %>%
  select(eid, t2dm_event_dt) %>%
  unique.data.frame()

ckd_indic_min <- ckd_indic %>%
  mutate(eid = as.character(eid)) %>%
  group_by(eid) %>%
  # get earliest CKD event date
  filter(ckd_event_dt == min(ckd_event_dt, na.rm = T)) %>%
  select(eid, ckd_event_dt) %>%
  unique.data.frame()

# join with T2DM + CKD earliest times with diagnosis data

diags %<>%
  left_join(t2dm_indic_min, by = "eid") %>%
  left_join(ckd_indic_min, by = "eid")

# create end point df
# if individual develops CKD, this is their end point
# if they do not, label them as "End"

ckd_event <- diags %>%
  ungroup() %>%
  select(eid, t2dm_event_dt, ckd_event_dt) %>%
  unique.data.frame() %>%
  mutate(block = ifelse(is.na(ckd_event_dt), "End", "CKD"), event_dt = ckd_event_dt) %>%
  select(eid, event_dt, block, t2dm_event_dt, ckd_event_dt)

# join

diags %<>%
  mutate(block = as.character(block))

diags <- rbind(diags, ckd_event)

# filter out diseases <= 10% prevalence

diags %<>%
  filter(block %in% c(as.character(dis_use), "CKD", "End"))

# filter out diagnoses which occur after CKD diagnosis for those that develop CKD

diags %<>%
  filter(!event_dt > ckd_event_dt | is.na(ckd_event_dt))

# convert dates of diagnoses to relative time from T2DM indication

diags %<>%
  ungroup() %>%
  mutate(event_dt = as.Date(event_dt) - as.Date(t2dm_event_dt))

# filter out those individuals whose earliest CKD indication is before their earliest T2DM indication

ckd_bef <- diags %>%
  filter(block == "CKD" & event_dt <= 0) %>%
  pull(eid)
# how many?
length(ckd_bef)
# 1458 individuals with CKD -> T2DM

# filter out those individuals
diags %<>%
  filter(!eid %in% ckd_bef)

# get age at which T2DM indication occurred

# get age at baseline UK attendance
age <- ukb_data %>%
  transmute(eid = as.character(f.eid), age = f.21003.0.0)

# get sex

sex <- ukb_data %>%
  transmute(eid = as.character(f.eid), sex = f.31.0.0)

# join with diagnoses df

bl_dates %<>%
  transmute(eid = as.character(f.eid), bl_date = f.53.0.0)

diags %<>%
  left_join(bl_dates, by = "eid")

# get time difference between date of attendance of assessment centre and date of T2DM indication

diags %<>%
  # get difference
  mutate(t2dm_diff = as.Date(t2dm_event_dt) - as.Date(bl_date)) %>%
  # convert to nearest year
  mutate(t2dm_diff = round(as.numeric(t2dm_diff/365)))

# join diagnoses with sex and age

diags %<>%
  left_join(age, by = "eid") %>%
  left_join(sex, by = "eid")

# get age at earliest T2DM indication

diags %<>%
  mutate(t2dm_age = age + t2dm_diff)

# get death variable from UKB data

death <- ukb_data %>%
  transmute(eid = as.character(f.eid), death = as.Date(f.40000.0.0))

# join diagnoses df with death

diags %<>%
  left_join(death, by = "eid")

# NA those deaths after CKD

diags$death[diags$death > diags$ckd_event_dt] <- NA

# get relative time of death from T2DM indication

diags %<>%
  mutate(death = death - t2dm_event_dt)

# add death event

tmp <- diags %>%
  # filter for those who died
  filter(!is.na(death)) %>%
  # add death event
  mutate(event_dt = death, block = "Death") %>%
  unique.data.frame()

# add in death events

diags <- rbind(diags, tmp) %>%
  arrange(eid, event_dt)

# remove End

diags %<>%
  filter(block != "End")

# check again and remove extremely early age of T2DM diagnoses (same age cutoff as our exclusion criteria)

early_t2dm <- diags %>%
  filter(t2dm_age < 36) %>%
  pull(eid) %>%
  unique()

diags %<>%
  filter(!eid %in% early_t2dm)

# plot distribution time to CKD from T2DM

diags %>%
  filter(block == "CKD") %>%
  ggplot(aes(x = event_dt/365)) +
  stat_ecdf() +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 35, 5)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1))

# ~95% of cases of CKD stage 1-3 happen within 20 years of T2DM diagnosis




# Cox PH modelling of comorbidities




# model the effect of history of comorbidities

# create model df

# create model input
model_df <- diags %>%
  select(eid, event_dt, block, age, sex) %>%
  mutate(event_dt = as.numeric(event_dt))

# create coefficient dataframe

coef_df <- data.frame(matrix(nrow = 0, ncol = 8))

# get comorbidities
comorbs <- unique(model_df$block)
# remove CKD and Death
comorbs <- comorbs[!comorbs %in% c("Death", "CKD")]

# create follow-up times
futimes <- model_df %>%
  group_by(eid) %>%
  # select for either time of CKD or time of last diagnosis or time of death, these are the follow-up times
  filter(event_dt == max(event_dt)) %>%
  transmute(eid, futime = event_dt) %>%
  unique.data.frame()

# for each comorbidity history
for(i in 1:length(comorbs)) {
  
  # select for comorbidity i
  cmb <- comorbs[i]
  
  # create df for testing comorbidity i
  test <- model_df %>%
    # convert comorbidity i name to "comorb"
    mutate(block = ifelse(block == cmb, "comorb", block)) %>%
    # pivot wider
    pivot_wider(names_from = block, values_from = event_dt) %>%
    # filter out other comorbidities than comorbidity i
    select(eid, sex, age, comorb, CKD, Death) %>%
    # join with follow-up times
    full_join(futimes, by = "eid") %>%
    # convert follow-up times of t = 0 to t = 1
    mutate(futime = ifelse(futime == 0, futime + 1, futime))
  
  # remove those individuals where T2DM is their last diagnosis
  # as there is no comorbidities after T2DM to model and convert
  # comorbidity time to whether it was before T2DM or not/never occurred
  
  test %<>%
    mutate(comorb = ifelse(comorb >= 0 | is.na(comorb), 0, 1)) %>%
    filter(!futime < 0)
  
  # create survival table input for Cox PH
  
  test_new <- tmerge(data1 = test[, 1:4], data2 = test, id = eid, tstop = futime)
  # create CKD event
  test_new <- tmerge(test_new, test, id = eid, CKD = event(CKD))
  
  # format variables for input to model
  
  test_new %<>%
    mutate(age = as.numeric(age), sex = as.factor(ifelse(sex == 1, "male", "female")), comorb = factor(ifelse(comorb == 1, "cmb", "no_cmb"), levels = c("no_cmb", "cmb")))
  
  # run Cox PH model, modelling the effect of history of comorbidity i before T2DM on risk of CKD after T2DM,
  # controlling for age at T2DM indication and sex
  
  cmb_model <- coxph(Surv(tstart, tstop, CKD) ~ age + sex + comorb,
                     data = test_new, cluster = eid)
  
  # create model summary df
  model_sum <- summary(cmb_model)
  model_sum <- data.frame(model_sum$coefficients)
  
  # get confidence intervals of hazard ratios
  cmb_confint <- confint(cmb_model)
  colnames(cmb_confint) <- c("ci_2_5_pct", "ci_97_5_pct")
  # create output df
  model_sum <- data.frame(comorb = cmb, var = rownames(model_sum), model_sum, cmb_confint)
  
  # add coefficient output to coefficient df
  coef_df <- rbind(coef_df, model_sum)
  print(i)
}

# adjust p-values using holm

coef_df %<>%
  mutate(Pr...z.. = p.adjust(Pr...z.., method = "holm"))

# get significant comorbidity histories

sig_cmb_bef <- coef_df %>%
  filter(var == "comorbcmb") %>%
  filter(Pr...z.. < 0.05) %>%
  pull(comorb)

# get hazard ratios of significant comorbidity histories
sig_bef_hr <- coef_df %>%
  filter(var == "comorbcmb") %>%
  filter(Pr...z.. < 0.05) %>%
  select(comorb, exp.coef.)

# get hazard ratios of non-significant comorbidity histories

non_sig_bef_hr <- coef_df %>%
  filter(var == "comorbcmb") %>%
  filter(Pr...z.. >= 0.05) %>%
  select(comorb, exp.coef.)

# get non-significant comorbidity histories

non_sig_bef <- coef_df %>%
  filter(var == "comorbcmb") %>%
  filter(Pr...z.. >= 0.05) %>%
  pull(comorb)

# visualise hazard ratios of significant comorbidity histories

plot <- coef_df %>%
  filter(var == "comorbcmb") %>%
  filter(Pr...z.. < 0.05) %>%
  mutate(name = phecode_icd10$Phenotype[match(comorb, phecode_icd10$PheCode)]) %>%
  mutate(name = str_trunc(name, 40)) %>%
  mutate(name = factor(name, levels = name[order(exp.coef.)])) %>%
  ggplot(aes(x = exp.coef., y = name, xmax = exp(ci_97_5_pct), xmin = exp(ci_2_5_pct))) +
  geom_point() +
  geom_errorbar() +
  xlab("Hazard Ratio") +
  ylab("History of comorbidity before T2DM") +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  geom_vline(xintercept = 1)
ggsave(filename = paste0(plot_dir, "t2dm_to_ckd_all_cox_ph_phecode_cmb_bef_hr_plot.png"), units = "px",
       width = 1800, height = 1600)

# visualise hazard ratios of non-significant comorbidity histories

plot <- coef_df %>%
  filter(var == "comorbcmb") %>%
  filter(Pr...z.. >= 0.05) %>%
  mutate(name = phecode_icd10$Phenotype[match(comorb, phecode_icd10$PheCode)]) %>%
  mutate(name = str_trunc(name, 40)) %>%
  mutate(name = factor(name, levels = name[order(exp.coef.)])) %>%
  ggplot(aes(x = exp.coef., y = name, xmax = exp(ci_97_5_pct), xmin = exp(ci_2_5_pct))) +
  geom_point() +
  geom_errorbar() +
  xlab("Hazard Ratio") +
  ylab("History of comorbidity before T2DM") +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 2.5, 0.1)) +
  geom_vline(xintercept = 1)
ggsave(filename = paste0(plot_dir, "t2dm_to_ckd_all_cox_ph_phecode_cmb_bef_non_sig_hr_plot.png"), units = "px",
       width = 1800, height = 1600)



# re-create model df

model_df <- diags %>%
  select(eid, event_dt, block, age, sex) %>%
  mutate(event_dt = as.numeric(event_dt))

# get size of T2DM
t2dm_size <- model_df %>%
  pull(eid) %>%
  unique() %>%
  length()

# plot prevalence rates of signficant comorbidity histories before T2DM
# ordered by their risk of CKD after T2DM

plot <- model_df %>%
  filter(event_dt < 0 & block %in% sig_cmb_bef) %>%
  ungroup() %>%
  count(block) %>%
  mutate(name = str_trunc(phecode_icd10$Phenotype[match(block, phecode_icd10$PheCode)], 40)) %>%
  full_join(sig_bef_hr, by = c("block" = "comorb")) %>%
  mutate(name = factor(name, levels = name[order(exp.coef.)]),
         prev = n/t2dm_size) %>%
  ggplot(aes(y = name, x = prev)) +
  geom_col() +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 0.3, 0.05)) +
  xlab("Frequency of T2DM cases with\nhistory of disease before T2DM") +
  ylab("Disease history ordered by risk of CKD after T2DM")
ggsave(filename = paste0(plot_dir, "t2dm_to_ckd_all_sig_phecode_cmb_bef_prev_rate_at_t2dm.png"), plot = plot,
       units = "px", width = 1800, height = 1600)


# plot prevalence rates of non-signficant comorbidity histories before T2DM
# ordered by their risk of CKD after T2DM
plot <- model_df %>%
  filter(event_dt < 0 & block %in% non_sig_bef) %>%
  ungroup() %>%
  count(block) %>%
  mutate(name = str_trunc(phecode_icd10$Phenotype[match(block, phecode_icd10$PheCode)], 40)) %>%
  full_join(non_sig_bef_hr, by = c("block" = "comorb")) %>%
  mutate(name = factor(name, levels = name[order(exp.coef.)]),
         prev = n/t2dm_size) %>%
  ggplot(aes(y = name, x = prev)) +
  geom_col() +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 0.3, 0.05)) +
  xlab("Frequency of T2DM cases with\nhistory of disease before T2DM") +
  ylab("Non-significant disease history ordered by risk of CKD after T2DM")
ggsave(filename = paste0(plot_dir, "t2dm_to_ckd_all_non_sig_phecode_cmb_bef_prev_rate_at_t2dm.png"), plot = plot,
       units = "px", width = 1800, height = 1600)


# get which comorbidity histories we can ignore based on 
# low evidence of effect on CKD progression after T2DM
# i.e do not show significant effect on progression to CKD
# in single comorbidity history models and do not retain
# significance when controlling for each of the other
# comorbidity histories

# get the comorbidity histories which do have evidence of
# significant effect on CKD progression after T2DM and
# which we should control for in forward models

cmb_bef_cntrl <- sig_cmb_bef
cmb_bef_ignor <- non_sig_bef

# forward modelling of time-dependent comorbidity events
# after T2DM on progression to CKD
# for this, use MCA for dimensionality reduction of the
# important comorbidity histories to attain a set of
# reduce features for controlling for comorbidity histories
# in the foward model

# recreate model df

model_df <- diags %>%
  select(eid, event_dt, block, age, sex) %>%
  mutate(event_dt = as.numeric(event_dt))

# only consider those comorbidity histories which have evidence of
# effect on CKD progression after T2DM 

model_df <- model_df[-which(model_df$event_dt < 0 & model_df$block %in% cmb_bef_ignor),]

# create comorbidity history presence/absence matrix
# split this into year blocks; <2.5 years before T2DM, <5 years before T2DM, < 10 years, etc.. 

cmb_mm_1 <- model_df %>%
  pivot_wider(names_from = block, values_from = event_dt) %>%
  ungroup() %>%
  # convert comorbidity times after T2DM or which are NA to 0 and
  # those comorbidity times before T2DM to 1
  mutate(across(-all_of(c("eid", "age", "sex")), .fns = function(x) ifelse(is.na(x) | x >= 0, 0, 1))) %>%
  # select for important histories
  select(all_of(cmb_bef_cntrl)) %>%
  as.matrix()
colnames(cmb_mm_1) <- paste0(colnames(cmb_mm_1), "_2_5_yr")

cmb_mm_2 <- model_df %>%
  pivot_wider(names_from = block, values_from = event_dt) %>%
  ungroup() %>%
  # convert comorbidity times after T2DM or which are NA to 0 and
  # those comorbidity times before T2DM to 1
  mutate(across(-all_of(c("eid", "age", "sex")), .fns = function(x) ifelse(is.na(x) | x >= -2.5*(365) | x < -5*(365), 0, 1))) %>%
  # select for important histories
  select(all_of(cmb_bef_cntrl)) %>%
  as.matrix()
colnames(cmb_mm_2) <- paste0(colnames(cmb_mm_2), "_5_yr")

cmb_mm_3 <- model_df %>%
  pivot_wider(names_from = block, values_from = event_dt) %>%
  ungroup() %>%
  # convert comorbidity times after T2DM or which are NA to 0 and
  # those comorbidity times before T2DM to 1
  mutate(across(-all_of(c("eid", "age", "sex")), .fns = function(x) ifelse(is.na(x) | x >= -5*(365), 0, 1))) %>%
  # select for important histories
  select(all_of(cmb_bef_cntrl)) %>%
  as.matrix()
colnames(cmb_mm_3) <- paste0(colnames(cmb_mm_3), "_gt_5_yr")

# combine

cmb_mm <- cbind(cmb_mm_1, cmb_mm_2, cmb_mm_3)

# convert to logical

cmb_mm <- apply(cmb_mm, 2, as.logical)

# run MCA

cmb_mca <- MCA(cmb_mm)

# plot scree

plot(cmb_mca$eig[,1])

# select number of components using eigenvalues of components > 1/nvar rule

ncp_use <- sum(cmb_mca$eig[,1] > 1/ncol(cmb_mm))

cmb_mca$eig[ncp_use, 3]
# 50.33329% of variance with 31
cmb_mca$eig[7, 3]
# 17.154% with 7 based on scree plot

# re-reun MCA with number of components to retain

cmb_mca <- MCA(cmb_mm, ncp = ncp_use)

# get coordinates of individuals individual

indiv_coor <- as.data.frame(cmb_mca$ind$coord)

# run Cox PH forward following previous processing; controlling for MCA coordinates

# add MCA coordinates to model

model_df %<>% 
  pivot_wider(names_from = block, values_from = event_dt) %>%
  ungroup()
colnames(indiv_coor) <- paste0("mca_", 1:ncol(indiv_coor))
model_df <- cbind(model_df, indiv_coor)
model_df %<>%
  pivot_longer(-contains(c("eid", "age", "sex", "mca")),
               names_to = "block", values_to = "event_dt") %>%
  filter(!is.na(event_dt))

# create coefficient dataframe for forward modelling of comorbidities after T2DM

coef_df <- data.frame(matrix(nrow = 0, ncol = 8))

comorbs <- unique(model_df$block)
comorbs <- comorbs[!comorbs %in% c("Death", "CKD")]

# for each comorbidity event after T2DM
for(i in 1:length(comorbs)) {
  
  cmb <- comorbs[i]
  
  # create test data.frame for comorbidity i
  
  test <- model_df %>%
    # convert comorbidity i to  "comorb"
    mutate(block = ifelse(block == cmb, "comorb", block)) %>%
    # pivot wider
    pivot_wider(names_from = block, values_from = event_dt) %>%
    # select for relevant variables
    select(eid, sex, age, contains("mca"), comorb, CKD, Death) %>%
    # join with follow-up times
    full_join(futimes, by = "eid") %>%
    # convert follow-up times which = 0 to = 1
    mutate(futime = ifelse(futime == 0, futime + 1, futime))
  
  # remove those where T2DM is last disease as there are no
  # comorbidities to model after T2DM
  # set comorbidities before T2DM to NA as we are controlling for
  # the ones affecting CKD progression with the MCA coordinates
  
  test %<>%
    filter(futime > 0) %>%
    mutate(comorb = ifelse(comorb < 0, NA, comorb))
  
  # create survival input table
  
  test_new <- tmerge(data1 = select(test, contains(c("eid", "sex", "age", "mca"))), data2 = test, id = eid, tstop = futime)
  # create CKD event
  test_new <- tmerge(test_new, test, id = eid, CKD = event(CKD))
  # create time-dependent CKD event
  test_new <- tmerge(test_new, test, id = eid, comorb = tdc(comorb))
  
  # format variables
  
  test_new %<>%
    mutate(age = as.numeric(age), sex = as.factor(ifelse(sex == 1, "male", "female")), comorb = factor(ifelse(comorb == 1, "cmb", "no_cmb"), levels = c("no_cmb", "cmb")))
  
  # run Cox PH modelling of time-dependent comorbidity development after T2DM on effect on CKD progression
  # controlling for age at T2DM indication, sex, and previous comorbidity history 
  # affecting CKD progression represented in MCA coordinates
  
  # create model formula
  model_form <- test_new %>%
    select(contains(c("age", "sex", "comorb", "mca"))) %>%
    colnames()
  model_form <- as.formula(paste0("Surv(tstart, tstop, CKD) ~ ", paste(model_form, collapse = "+")))
  
  # run model
  
  cmb_model <- coxph(model_form,
                     data = test_new, cluster = eid)
  
  # create model summary
  model_sum <- summary(cmb_model)
  # get coefficients
  model_sum <- data.frame(model_sum$coefficients)
  # get confidence intervals of hazard ratios
  cmb_confint <- confint(cmb_model)
  colnames(cmb_confint) <- c("ci_2_5_pct", "ci_97_5_pct")
  # create model summary output
  model_sum <- data.frame(comorb = cmb, var = rownames(model_sum), model_sum, cmb_confint)
  
  # add model summary to coefficientdf
  coef_df <- rbind(coef_df, model_sum)
  print(i)
}

# adjust p-values

coef_df %<>%
  mutate(Pr...z.. = p.adjust(Pr...z.., method = "holm"))

# get significant comorbidity which occur after T2DM and 
# are affecting CKD progression after T2DM independent of
# age, sex, and comorbidity history before T2DM

fin_sig_hr <- coef_df %>%
  filter(var %in% c("comorbcmb")) %>%
  filter(Pr...z.. < 0.05) %>%
  select(comorb, exp.coef.)

control_sig_hr <- coef_df %>%
  filter(var %in% c("comorbcmb")) %>%
  select(comorb, exp.coef., Pr...z..)

sig_cmb <- fin_sig_hr %>%
  pull(comorb)

hr_order <- fin_sig_hr$comorb[order(fin_sig_hr$exp.coef.)]
hr_order <- phecode_icd10$Phenotype[match(hr_order, phecode_icd10$PheCode)]

# visualise hazard ratios of comorbidity events on CKD progression

plot <- coef_df %>%
  filter(var == "comorbcmb") %>%
  filter(Pr...z.. < 0.05) %>%
  mutate(name = phecode_icd10$Phenotype[match(comorb, phecode_icd10$PheCode)]) %>%
  mutate(name = str_trunc(name, 40)) %>%
  mutate(name = factor(name, levels = name[order(exp.coef.)])) %>%
  ggplot(aes(x = exp.coef., y = name, xmax = exp(ci_97_5_pct), xmin = exp(ci_2_5_pct))) +
  geom_point() +
  geom_errorbar() +
  xlab("Hazard Ratio") +
  ylab("Comorbidity") +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  geom_vline(xintercept = 1)
ggsave(filename = paste0(plot_dir, "t2dm_to_ckd_final_processed_cox_ph_phecode_cmb_hr_plot.png"), units = "px",
       width = 1800, height = 1600)

# incidence rates of diseases between T2DM -> CKD

# get those with T2DM -> CKD

ckd_eid <- model_df %>%
  filter(block == "CKD") %>%
  pull(eid)

cmb_mm <- model_df %>%
  filter(eid %in% ckd_eid) %>%
  pivot_wider(names_from = block, values_from = event_dt) %>%
  ungroup() %>%
  select(-age, -sex, -eid) %>%
  select(-contains("mca")) %>%
  mutate(across(everything(), .fns = function(x) ifelse(is.na(x) | x < 0, 0, 1))) %>%
  select(all_of(sig_cmb)) %>%
  as.matrix()

# get incidence rates between T2DM -> CKD
cmb_prev <- data.frame(comorb = colnames(cmb_mm), n = colSums(cmb_mm), prev = colSums(cmb_mm)/nrow(cmb_mm))
cmb_prev %<>%
  arrange(desc(prev))

plot <- cmb_prev %>%
  mutate(comorb = phecode_icd10$Phenotype[match(comorb, phecode_icd10$PheCode)]) %>%
  mutate(comorb = factor(comorb, levels = hr_order)) %>%
  ggplot(aes(y = comorb, x = n)) +
  geom_col(fill = "grey", col = "black") +
  scale_y_discrete(labels = function(x) str_trunc(x, 40, "center")) +
  scale_x_continuous(breaks = seq(0, 1000, 100)) +
  xlab("Number of events between T2DM -> CKD") +
  ylab("Comorbidity ordered by risk of CKD\n----------------->") +
  theme_bw()
ggsave(filename = paste0(plot_dir, "t2dm_to_ckd_final_processed_phecode_prev_between_t2dm_ckd_plot.png"), units = "px",
       width = 1900, height = 1600)


# what does this look like when we don't control for the MCA coordinates?

# create coefficient dataframe for forward modelling of comorbidities after T2DM

coef_df <- data.frame(matrix(nrow = 0, ncol = 8))

comorbs <- unique(model_df$block)
comorbs <- comorbs[!comorbs %in% c("Death", "CKD")]

# for each comorbidity event after T2DM
for(i in 1:length(comorbs)) {
  
  cmb <- comorbs[i]
  
  # create test data.frame for comorbidity i
  
  test <- model_df %>%
    # convert comorbidity i to  "comorb"
    mutate(block = ifelse(block == cmb, "comorb", block)) %>%
    # pivot wider
    pivot_wider(names_from = block, values_from = event_dt) %>%
    # select for relevant variables
    select(eid, sex, age, contains("mca"), comorb, CKD, Death) %>%
    # join with follow-up times
    full_join(futimes, by = "eid") %>%
    # convert follow-up times which = 0 to = 1
    mutate(futime = ifelse(futime == 0, futime + 1, futime))
  
  # remove those where T2DM is last disease as there are no
  # comorbidities to model after T2DM
  # set comorbidities before T2DM to NA as we are controlling for
  # the ones affecting CKD progression with the MCA coordinates
  
  test %<>%
    filter(futime > 0) %>%
    mutate(comorb = ifelse(comorb < 0, NA, comorb))
  
  # create survival input table
  
  test_new <- tmerge(data1 = select(test, contains(c("eid", "sex", "age", "mca"))), data2 = test, id = eid, tstop = futime)
  # create CKD event
  test_new <- tmerge(test_new, test, id = eid, CKD = event(CKD))
  # create time-dependent CKD event
  test_new <- tmerge(test_new, test, id = eid, comorb = tdc(comorb))
  
  # format variables
  
  test_new %<>%
    mutate(age = as.numeric(age), sex = as.factor(ifelse(sex == 1, "male", "female")), comorb = factor(ifelse(comorb == 1, "cmb", "no_cmb"), levels = c("no_cmb", "cmb")))
  
  # run Cox PH modelling of time-dependent comorbidity development after T2DM on effect on CKD progression
  # controlling for age at T2DM indication, sex, and previous comorbidity history 
  # affecting CKD progression represented in MCA coordinates
  
  # create model formula
  model_form <- test_new %>%
    select(contains(c("age", "sex", "comorb"))) %>%
    colnames()
  model_form <- as.formula(paste0("Surv(tstart, tstop, CKD) ~ ", paste(model_form, collapse = "+")))
  
  # run model
  
  cmb_model <- coxph(model_form,
                     data = test_new, cluster = eid)
  
  # create model summary
  model_sum <- summary(cmb_model)
  # get coefficients
  model_sum <- data.frame(model_sum$coefficients)
  # get confidence intervals of hazard ratios
  cmb_confint <- confint(cmb_model)
  colnames(cmb_confint) <- c("ci_2_5_pct", "ci_97_5_pct")
  # create model summary output
  model_sum <- data.frame(comorb = cmb, var = rownames(model_sum), model_sum, cmb_confint)
  
  # add model summary to coefficient df
  coef_df <- rbind(coef_df, model_sum)
  print(i)
}

# adjust p-values

coef_df %<>%
  mutate(Pr...z.. = p.adjust(Pr...z.., method = "holm"))

# get significant comorbidity which occur after T2DM and 
# are affecting CKD progression after T2DM independent of
# age, sex, and comorbidity history before T2DM

uncontrol_sig_hr <- coef_df %>%
  filter(var %in% c("comorbcmb")) %>%
  select(comorb, exp.coef., Pr...z..)

# visualise difference in hazard ratios between MCA-controlled and non-controlled results

plot <- control_sig_hr %>%
  full_join(uncontrol_sig_hr, "comorb", suffix = c("control", "uncontrol")) %>%
  mutate(diff = exp.coef.control/exp.coef.uncontrol) %>%
  mutate(sig_change = ifelse(Pr...z..uncontrol < 0.05 & Pr...z..control >= 0.05, "change_to_non_sig", "no_change")) %>%
  mutate(sig_change = ifelse(Pr...z..uncontrol >= 0.05 & Pr...z..control < 0.05, "change_to_sig", sig_change)) %>%
  mutate(comorb = str_trunc(phecode_icd10$Phenotype[match(comorb, phecode_icd10$PheCode)], 30)) %>%
  mutate(comorb = factor(comorb, levels = comorb[order(exp.coef.control)])) %>%
  ggplot(aes(x = diff, y = comorb, col = sig_change)) +
  geom_point(pch = 3) +
  theme_bw() +
  geom_vline(xintercept = 1, lty = 2) +
  xlab("Relative change in hazard\nratio after controlling for\nMCA coordinates") +
  ylab("Comorbidity event after T2DM") +
  scale_x_continuous(breaks = seq(0.85, 1.075, 0.05)) +
  labs(col = "Change in significance\nof hazard ratio") +
  coord_cartesian(xlim = c(0.85, 1.075)) +
  annotate(x = 0.925, y = 5, xend = 0.86, yend = 5, geom = "segment", arrow = arrow(length = unit(0.1, "inches"))) +
  annotate(x = 1.0125, y = 5, xend = 1.075, yend = 5, geom = "segment", arrow = arrow(length = unit(0.1, "inches"))) +
  annotate(x = 0.89, y = 9, geom = "text", label = "Decrease\nin HR", size = 4) +
  annotate(x = 1.045, y = 9, geom = "text", label = "Increase\nin HR", size = 4) +
  scale_color_discrete(labels = c(change_to_non_sig = "Change to\nnon-significant",
                                  no_change = "No change in\nsignificance"))
ggsave(filename = paste0(plot_dir, "t2dm_to_ckd_final_processed_phecode_cox_ph_change_in_hr_control.png"), plot = plot,
       units = "px", width = 2000, height = 2000)


# filter out events after T2DM in model_df that are not significant

df_out <- model_df %>%
  filter(event_dt >= 0 & block %in% c(sig_cmb, "CKD", "Death") | event_dt < 0 & block %in% cmb_bef_cntrl)

# save environment

save.image(paste0(work_dir, "/intermediates/t2dm_time_to_ckd_envir.RData"))

load(paste0(work_dir, "/intermediates/t2dm_time_to_ckd_envir.RData"))
non_t2dm <- readRDS(paste0(work_dir, "/intermediates/non_t2dm_to_ckd_output.rds"))
# 
# # investigation of modifiers of comorbidities and associated risk of CKD
# 
# # create coefficient dataframe for forward modelling of comorbidities after T2DM
# 
# mod_coef_df <- data.frame(matrix(nrow = 0, ncol = 8))
# 
# comorbs <- unique(model_df$block)
# comorbs <- comorbs[!comorbs %in% c("Death", "CKD")]
# 
# # for each comorbidity event after T2DM
# for(i in 1:length(comorbs)) {
#   
#   cmb <- comorbs[i]
#   
#   # create test data.frame for comorbidity i
#   
#   test <- model_df %>%
#     # convert comorbidity i to  "comorb"
#     mutate(block = ifelse(block == cmb, "comorb", block)) %>%
#     # pivot wider
#     pivot_wider(names_from = block, values_from = event_dt) %>%
#     # select for relevant variables
#     select(eid, sex, age, contains("mca"), comorb, CKD, Death) %>%
#     # join with follow-up times
#     full_join(futimes, by = "eid") %>%
#     # convert follow-up times which = 0 to = 1
#     mutate(futime = ifelse(futime == 0, futime + 1, futime))
#   
#   # remove those where T2DM is last disease as there are no
#   # comorbidities to model after T2DM
#   # set comorbidities before T2DM to NA as we are controlling for
#   # the ones affecting CKD progression with the MCA coordinates
#   
#   test %<>%
#     filter(futime > 0) %>%
#     mutate(comorb = ifelse(comorb < 0, NA, comorb))
#   
#   # create survival input table
#   
#   test_new <- tmerge(data1 = select(test, contains(c("eid", "sex", "age", "mca"))), data2 = test, id = eid, tstop = futime)
#   # create CKD event
#   test_new <- tmerge(test_new, test, id = eid, CKD = event(CKD))
#   # create time-dependent CKD event
#   test_new <- tmerge(test_new, test, id = eid, comorb = tdc(comorb))
#   
#   # format variables
#   
#   test_new %<>%
#     mutate(age = as.numeric(age), sex = as.factor(ifelse(sex == 1, "male", "female")), comorb = factor(ifelse(comorb == 1, "cmb", "no_cmb"), levels = c("no_cmb", "cmb")))
#   
#   # run Cox PH modelling of time-dependent comorbidity development after T2DM on effect on CKD progression
#   # controlling for age at T2DM indication, sex, and previous comorbidity history 
#   # affecting CKD progression represented in MCA coordinates
#   
#   # create model formula
#   model_form <- test_new %>%
#     select(contains(c("age", "sex", "comorb"))) %>%
#     colnames()
#   model_form <- paste0("Surv(tstart, tstop, CKD) ~ ", paste(model_form, collapse = " + "))
#   model_form <- as.formula(paste0(model_form, " + age:comorb")) # add age interaction term
#   
#   # run model
#   
#   cmb_model <- coxph(model_form,
#                      data = test_new, cluster = eid)
#   
#   # create model summary
#   model_sum <- summary(cmb_model)
#   # get coefficients
#   model_sum <- data.frame(model_sum$coefficients)
#   # get confidence intervals of hazard ratios
#   cmb_confint <- confint(cmb_model)
#   colnames(cmb_confint) <- c("ci_2_5_pct", "ci_97_5_pct")
#   # create model summary output
#   model_sum <- data.frame(comorb = cmb, var = rownames(model_sum), model_sum, cmb_confint)
#   
#   # add model summary to coefficient df
#   mod_coef_df <- rbind(mod_coef_df, model_sum)
#   print(i)
# }
# 
# # adjust p-values
# 
# mod_coef_df %<>%
#   mutate(Pr...z.. = p.adjust(Pr...z.., method = "holm"))
# 
# mod_coef_df %>%
#   filter(var == "age:comorbcmb") %>%
#   filter(Pr...z.. < 0.05) %>%
#   arrange(desc(exp.coef.))
# # no significant interactions between age and comorbidity risk
# 
# 
# # create coefficient dataframe for forward modelling of comorbidities after T2DM
# 
# mod_coef_df <- data.frame(matrix(nrow = 0, ncol = 8))
# 
# comorbs <- unique(model_df$block)
# comorbs <- comorbs[!comorbs %in% c("Death", "CKD")]
# 
# # for each comorbidity event after T2DM
# for(i in 1:length(comorbs)) {
#   
#   cmb <- comorbs[i]
#   
#   # create test data.frame for comorbidity i
#   
#   test <- model_df %>%
#     # convert comorbidity i to  "comorb"
#     mutate(block = ifelse(block == cmb, "comorb", block)) %>%
#     # pivot wider
#     pivot_wider(names_from = block, values_from = event_dt) %>%
#     # select for relevant variables
#     select(eid, sex, age, contains("mca"), comorb, CKD, Death) %>%
#     # join with follow-up times
#     full_join(futimes, by = "eid") %>%
#     # convert follow-up times which = 0 to = 1
#     mutate(futime = ifelse(futime == 0, futime + 1, futime))
#   
#   # remove those where T2DM is last disease as there are no
#   # comorbidities to model after T2DM
#   # set comorbidities before T2DM to NA as we are controlling for
#   # the ones affecting CKD progression with the MCA coordinates
#   
#   test %<>%
#     filter(futime > 0) %>%
#     mutate(comorb = ifelse(comorb < 0, NA, comorb))
#   
#   # create survival input table
#   
#   test_new <- tmerge(data1 = select(test, contains(c("eid", "sex", "age", "mca"))), data2 = test, id = eid, tstop = futime)
#   # create CKD event
#   test_new <- tmerge(test_new, test, id = eid, CKD = event(CKD))
#   # create time-dependent CKD event
#   test_new <- tmerge(test_new, test, id = eid, comorb = tdc(comorb))
#   
#   # format variables
#   
#   test_new %<>%
#     mutate(age = as.numeric(age), sex = as.factor(ifelse(sex == 1, "male", "female")), comorb = factor(ifelse(comorb == 1, "cmb", "no_cmb"), levels = c("no_cmb", "cmb")))
#   
#   # run Cox PH modelling of time-dependent comorbidity development after T2DM on effect on CKD progression
#   # controlling for age at T2DM indication, sex, and previous comorbidity history 
#   # affecting CKD progression represented in MCA coordinates
#   
#   # create model formula
#   model_form <- test_new %>%
#     select(contains(c("age", "sex", "comorb", "mca"))) %>%
#     colnames()
#   model_form <- as.formula(paste0("Surv(tstart, tstop, CKD) ~ ", paste(model_form, collapse = " + ")))
#   # model_form <- as.formula(paste0(model_form, " + sex:comorb")) # add sex interaction term
#   
#   # run model
#   
#   cmb_model <- coxph(model_form,
#                      data = test_new, cluster = eid)
#   
#   # create model summary
#   model_sum <- summary(cmb_model)
#   # get coefficients
#   model_sum <- data.frame(model_sum$coefficients)
#   # get confidence intervals of hazard ratios
#   cmb_confint <- confint(cmb_model)
#   colnames(cmb_confint) <- c("ci_2_5_pct", "ci_97_5_pct")
#   # create model summary output
#   model_sum <- data.frame(comorb = cmb, var = rownames(model_sum), model_sum, cmb_confint)
#   
#   # add model summary to coefficient df
#   mod_coef_df <- rbind(mod_coef_df, model_sum)
#   print(i)
# }
# 
# # adjust p-values
# 
# mod_coef_df %<>%
#   mutate(Pr...z.. = p.adjust(Pr...z.., method = "holm"))
# 
# mod_coef_df %>%
#   filter(var == "comorbcmb") %>%
#   filter(Pr...z.. < 0.05) %>%
#   arrange(desc(exp.coef.))
# # no significant sex interactions
# 
# 
# 
# 
# # interaction between history of diseases vs not
# 
# cmb_mm <- model_df %>%
#   pivot_wider(names_from = block, values_from = event_dt) %>%
#   ungroup() %>%
#   # convert comorbidity times after T2DM or which are NA to 0 and
#   # those comorbidity times before T2DM to 1
#   mutate(across(-all_of(c("eid", "age", "sex")), .fns = function(x) ifelse(is.na(x) | x >= 0, 0, 1))) %>%
#   # select for important histories
#   select(all_of(cmb_bef_cntrl)) %>%
#   as.matrix()
# 
# dis_bef <- apply(cmb_mm, 1, function(x) all(x == 1))
# 
# model_df %<>% 
#   pivot_wider(names_from = block, values_from = event_dt) %>%
#   ungroup()
# 
# model_df %<>%
#   mutate(dis_bef = dis_bef)
# model_df %<>%
#   pivot_longer(-contains(c("eid", "age", "sex", "mca", "dis_bef")),
#                names_to = "block", values_to = "event_dt") %>%
#   filter(!is.na(event_dt))
# 
# 
# # create coefficient dataframe for forward modelling of comorbidities after T2DM
# 
# mod_coef_df <- data.frame(matrix(nrow = 0, ncol = 8))
# 
# comorbs <- unique(model_df$block)
# comorbs <- comorbs[!comorbs %in% c("Death", "CKD")]
# 
# # for each comorbidity event after T2DM
# for(i in 1:length(comorbs)) {
#   
#   cmb <- comorbs[i]
#   
#   # create test data.frame for comorbidity i
#   
#   test <- model_df %>%
#     # convert comorbidity i to  "comorb"
#     mutate(block = ifelse(block == cmb, "comorb", block)) %>%
#     # pivot wider
#     pivot_wider(names_from = block, values_from = event_dt) %>%
#     # select for relevant variables
#     select(eid, sex, age, dis_bef, comorb, CKD, Death) %>%
#     # join with follow-up times
#     full_join(futimes, by = "eid") %>%
#     # convert follow-up times which = 0 to = 1
#     mutate(futime = ifelse(futime == 0, futime + 1, futime))
#   
#   # remove those where T2DM is last disease as there are no
#   # comorbidities to model after T2DM
#   # set comorbidities before T2DM to NA as we are controlling for
#   # the ones affecting CKD progression with the MCA coordinates
#   
#   test %<>%
#     filter(futime > 0) %>%
#     mutate(comorb = ifelse(comorb < 0, NA, comorb))
#   
#   # create survival input table
#   
#   test_new <- tmerge(data1 = select(test, contains(c("eid", "sex", "age", "dis_bef"))), data2 = test, id = eid, tstop = futime)
#   # create CKD event
#   test_new <- tmerge(test_new, test, id = eid, CKD = event(CKD))
#   # create time-dependent CKD event
#   test_new <- tmerge(test_new, test, id = eid, comorb = tdc(comorb))
#   
#   # format variables
#   
#   test_new %<>%
#     mutate(age = as.numeric(age), sex = as.factor(ifelse(sex == 1, "male", "female")), comorb = factor(ifelse(comorb == 1, "cmb", "no_cmb"), levels = c("no_cmb", "cmb")), dis_bef = factor(ifelse(dis_bef == 0, "no_hist", "hist")))
#   
#   # run Cox PH modelling of time-dependent comorbidity development after T2DM on effect on CKD progression
#   # controlling for age at T2DM indication, sex, and previous comorbidity history 
#   # affecting CKD progression represented in MCA coordinates
#   
#   # create model formula
#   model_form <- test_new %>%
#     select(contains(c("age", "sex", "dis_bef", "comorb"))) %>%
#     colnames()
#   model_form <- paste0("Surv(tstart, tstop, CKD) ~ ", paste(model_form, collapse = " + "))
#   model_form <- as.formula(paste0(model_form, " + dis_bef:comorb")) # add sex interaction term
#   
#   # run model
#   
#   cmb_model <- coxph(model_form,
#                      data = test_new, cluster = eid)
#   
#   # create model summary
#   model_sum <- summary(cmb_model)
#   # get coefficients
#   model_sum <- data.frame(model_sum$coefficients)
#   # get confidence intervals of hazard ratios
#   cmb_confint <- confint(cmb_model)
#   colnames(cmb_confint) <- c("ci_2_5_pct", "ci_97_5_pct")
#   # create model summary output
#   model_sum <- data.frame(comorb = cmb, var = rownames(model_sum), model_sum, cmb_confint)
#   
#   # add model summary to coefficient df
#   mod_coef_df <- rbind(mod_coef_df, model_sum)
#   print(i)
# }
# 
# # adjust p-values
# 
# mod_coef_df %<>%
#   mutate(Pr...z.. = p.adjust(Pr...z.., method = "holm"))
# 
# mod_coef_df %>%
#   filter(var == "sex:comorbcmb") %>%
#   filter(Pr...z.. < 0.05) %>%
#   arrange(desc(exp.coef.))
# 
# 
# 
# # stratified analysis approach
# 
# # create coefficient dataframe for forward modelling of comorbidities after T2DM
# 
# mod_coef_df <- data.frame(matrix(nrow = 0, ncol = 8))
# 
# comorbs <- unique(model_df$block)
# comorbs <- comorbs[!comorbs %in% c("Death", "CKD")]
# 
# # for each comorbidity event after T2DM
# for(i in 1:length(comorbs)) {
#   
#   cmb <- comorbs[i]
#   
#   # create test data.frame for comorbidity i
#   
#   test <- model_df %>%
#     # convert comorbidity i to  "comorb"
#     mutate(block = ifelse(block == cmb, "comorb", block)) %>%
#     # pivot wider
#     pivot_wider(names_from = block, values_from = event_dt) %>%
#     # select for relevant variables
#     select(eid, sex, age, dis_bef, comorb, CKD, Death) %>%
#     # join with follow-up times
#     full_join(futimes, by = "eid") %>%
#     # convert follow-up times which = 0 to = 1
#     mutate(futime = ifelse(futime == 0, futime + 1, futime)) %>%
#     filter(dis_bef == 0)
#   
#   # remove those where T2DM is last disease as there are no
#   # comorbidities to model after T2DM
#   # set comorbidities before T2DM to NA as we are controlling for
#   # the ones affecting CKD progression with the MCA coordinates
#   
#   test %<>%
#     filter(futime > 0) %>%
#     mutate(comorb = ifelse(comorb < 0, NA, comorb))
#   
#   # create survival input table
#   
#   test_new <- tmerge(data1 = select(test, contains(c("eid", "sex", "age", "dis_bef"))), data2 = test, id = eid, tstop = futime)
#   # create CKD event
#   test_new <- tmerge(test_new, test, id = eid, CKD = event(CKD))
#   # create time-dependent CKD event
#   test_new <- tmerge(test_new, test, id = eid, comorb = tdc(comorb))
#   
#   # format variables
#   
#   test_new %<>%
#     mutate(age = as.numeric(age), sex = as.factor(ifelse(sex == 1, "male", "female")), comorb = factor(ifelse(comorb == 1, "cmb", "no_cmb"), levels = c("no_cmb", "cmb")), dis_bef = factor(ifelse(dis_bef == 0, "no_hist", "hist")))
#   
#   # run Cox PH modelling of time-dependent comorbidity development after T2DM on effect on CKD progression
#   # controlling for age at T2DM indication, sex, and previous comorbidity history 
#   # affecting CKD progression represented in MCA coordinates
#   
#   # create model formula
#   model_form <- test_new %>%
#     select(contains(c("age", "sex", "comorb"))) %>%
#     colnames()
#   model_form <- as.formula(paste0("Surv(tstart, tstop, CKD) ~ ", paste(model_form, collapse = " + ")))
#   # model_form <- as.formula(paste0(model_form, " + dis_bef:comorb")) # add sex interaction term
#   
#   # run model
#   
#   cmb_model <- coxph(model_form,
#                      data = test_new, cluster = eid)
#   
#   # create model summary
#   model_sum <- summary(cmb_model)
#   # get coefficients
#   model_sum <- data.frame(model_sum$coefficients)
#   # get confidence intervals of hazard ratios
#   cmb_confint <- confint(cmb_model)
#   colnames(cmb_confint) <- c("ci_2_5_pct", "ci_97_5_pct")
#   # create model summary output
#   model_sum <- data.frame(comorb = cmb, var = rownames(model_sum), model_sum, cmb_confint)
#   
#   # add model summary to coefficient df
#   mod_coef_df <- rbind(mod_coef_df, model_sum)
#   print(i)
# }
# 
# # adjust p-values
# 
# mod_coef_df %<>%
#   mutate(Pr...z.. = p.adjust(Pr...z.., method = "holm"))
# 
# mod_coef_df %>%
#   filter(var == "comorbcmb") %>%
#   filter(Pr...z.. < 0.05) %>%
#   arrange(desc(exp.coef.))
# 
# 
# 
# 
# # create coefficient dataframe for forward modelling of comorbidities after T2DM
# 
# mod_coef_df <- data.frame(matrix(nrow = 0, ncol = 8))
# 
# comorbs <- unique(model_df$block)
# comorbs <- comorbs[!comorbs %in% c("Death", "CKD")]
# 
# # for each comorbidity event after T2DM
# for(i in 1:length(comorbs)) {
#   
#   cmb <- comorbs[i]
#   
#   # create test data.frame for comorbidity i
#   
#   test <- model_df %>%
#     # convert comorbidity i to  "comorb"
#     mutate(block = ifelse(block == cmb, "comorb", block)) %>%
#     # pivot wider
#     pivot_wider(names_from = block, values_from = event_dt) %>%
#     # select for relevant variables
#     select(eid, sex, age, dis_bef, comorb, CKD, Death) %>%
#     # join with follow-up times
#     full_join(futimes, by = "eid") %>%
#     # convert follow-up times which = 0 to = 1
#     mutate(futime = ifelse(futime == 0, futime + 1, futime)) %>%
#     filter(dis_bef == 1)
#   
#   # remove those where T2DM is last disease as there are no
#   # comorbidities to model after T2DM
#   # set comorbidities before T2DM to NA as we are controlling for
#   # the ones affecting CKD progression with the MCA coordinates
#   
#   test %<>%
#     filter(futime > 0) %>%
#     mutate(comorb = ifelse(comorb < 0, NA, comorb))
#   
#   # create survival input table
#   
#   test_new <- tmerge(data1 = select(test, contains(c("eid", "sex", "age", "dis_bef"))), data2 = test, id = eid, tstop = futime)
#   # create CKD event
#   test_new <- tmerge(test_new, test, id = eid, CKD = event(CKD))
#   # create time-dependent CKD event
#   test_new <- tmerge(test_new, test, id = eid, comorb = tdc(comorb))
#   
#   # format variables
#   
#   test_new %<>%
#     mutate(age = as.numeric(age), sex = as.factor(ifelse(sex == 1, "male", "female")), comorb = factor(ifelse(comorb == 1, "cmb", "no_cmb"), levels = c("no_cmb", "cmb")), dis_bef = factor(ifelse(dis_bef == 0, "no_hist", "hist")))
#   
#   # run Cox PH modelling of time-dependent comorbidity development after T2DM on effect on CKD progression
#   # controlling for age at T2DM indication, sex, and previous comorbidity history 
#   # affecting CKD progression represented in MCA coordinates
#   
#   # create model formula
#   model_form <- test_new %>%
#     select(contains(c("age", "sex", "comorb"))) %>%
#     colnames()
#   model_form <- as.formula(paste0("Surv(tstart, tstop, CKD) ~ ", paste(model_form, collapse = " + ")))
#   # model_form <- as.formula(paste0(model_form, " + dis_bef:comorb")) # add sex interaction term
#   
#   # run model
#   
#   cmb_model <- coxph(model_form,
#                      data = test_new, cluster = eid)
#   
#   # create model summary
#   model_sum <- summary(cmb_model)
#   # get coefficients
#   model_sum <- data.frame(model_sum$coefficients)
#   # get confidence intervals of hazard ratios
#   cmb_confint <- confint(cmb_model)
#   colnames(cmb_confint) <- c("ci_2_5_pct", "ci_97_5_pct")
#   # create model summary output
#   model_sum <- data.frame(comorb = cmb, var = rownames(model_sum), model_sum, cmb_confint)
#   
#   # add model summary to coefficient df
#   mod_coef_df <- rbind(mod_coef_df, model_sum)
#   print(i)
# }
# 
# # adjust p-values
# 
# mod_coef_df %<>%
#   mutate(Pr...z.. = p.adjust(Pr...z.., method = "holm"))
# 
# mod_coef_df %>%
#   filter(var == "comorbcmb") %>%
#   filter(Pr...z.. < 0.05) %>%
#   arrange(desc(exp.coef.))
# 
# 
# 
# # age stratification
# 
# t2dm_age <-  model_df %>%
#   select(eid, age) %>%
#   unique.data.frame()
# 
# t2dm_age <- t2dm_age %>%
#   mutate(age_grp = ifelse(age < median(age), 0, 1)) %>%
#   select(eid, age_grp)
# 
# model_df %<>%
#   full_join(t2dm_age, "eid")
# 
# # create coefficient dataframe for forward modelling of comorbidities after T2DM
# 
# mod_coef_df <- data.frame(matrix(nrow = 0, ncol = 8))
# 
# comorbs <- unique(model_df$block)
# comorbs <- comorbs[!comorbs %in% c("Death", "CKD")]
# 
# # for each comorbidity event after T2DM
# for(i in 1:length(comorbs)) {
#   
#   cmb <- comorbs[i]
#   
#   # create test data.frame for comorbidity i
#   
#   test <- model_df %>%
#     # convert comorbidity i to  "comorb"
#     mutate(block = ifelse(block == cmb, "comorb", block)) %>%
#     # pivot wider
#     pivot_wider(names_from = block, values_from = event_dt) %>%
#     # select for relevant variables
#     select(eid, sex, age, age_grp, contains("mca"), comorb, CKD, Death) %>%
#     # join with follow-up times
#     full_join(futimes, by = "eid") %>%
#     # convert follow-up times which = 0 to = 1
#     mutate(futime = ifelse(futime == 0, futime + 1, futime)) %>%
#     filter(age_grp == 1) %>%
#     select(-age_grp)
#   
#   # remove those where T2DM is last disease as there are no
#   # comorbidities to model after T2DM
#   # set comorbidities before T2DM to NA as we are controlling for
#   # the ones affecting CKD progression with the MCA coordinates
#   
#   test %<>%
#     filter(futime > 0) %>%
#     mutate(comorb = ifelse(comorb < 0, NA, comorb))
#   
#   # create survival input table
#   
#   test_new <- tmerge(data1 = select(test, contains(c("eid", "sex", "age", "mca"))), data2 = test, id = eid, tstop = futime)
#   # create CKD event
#   test_new <- tmerge(test_new, test, id = eid, CKD = event(CKD))
#   # create time-dependent CKD event
#   test_new <- tmerge(test_new, test, id = eid, comorb = tdc(comorb))
#   
#   # format variables
#   
#   test_new %<>%
#     mutate(age = as.numeric(age), sex = as.factor(ifelse(sex == 1, "male", "female")), comorb = factor(ifelse(comorb == 1, "cmb", "no_cmb"), levels = c("no_cmb", "cmb")))
#   
#   # run Cox PH modelling of time-dependent comorbidity development after T2DM on effect on CKD progression
#   # controlling for age at T2DM indication, sex, and previous comorbidity history 
#   # affecting CKD progression represented in MCA coordinates
#   
#   # create model formula
#   model_form <- test_new %>%
#     select(contains(c("age", "sex", "comorb", "mca"))) %>%
#     colnames()
#   model_form <- as.formula(paste0("Surv(tstart, tstop, CKD) ~ ", paste(model_form, collapse = " + ")))
#   # model_form <- as.formula(paste0(model_form, " + dis_bef:comorb")) # add sex interaction term
#   
#   # run model
#   
#   cmb_model <- coxph(model_form,
#                      data = test_new, cluster = eid)
#   
#   # create model summary
#   model_sum <- summary(cmb_model)
#   # get coefficients
#   model_sum <- data.frame(model_sum$coefficients)
#   # get confidence intervals of hazard ratios
#   cmb_confint <- confint(cmb_model)
#   colnames(cmb_confint) <- c("ci_2_5_pct", "ci_97_5_pct")
#   # create model summary output
#   model_sum <- data.frame(comorb = cmb, var = rownames(model_sum), model_sum, cmb_confint)
#   
#   # add model summary to coefficient df
#   mod_coef_df <- rbind(mod_coef_df, model_sum)
#   print(i)
# }
# 
# # adjust p-values
# 
# mod_coef_df %<>%
#   mutate(Pr...z.. = p.adjust(Pr...z.., method = "holm"))
# 
# mod_coef_df %>%
#   filter(var == "comorbcmb") %>%
#   filter(Pr...z.. < 0.05) %>%
#   arrange(desc(exp.coef.))
# 
# 
# 
# 
# 
# 
# # sex stratification
# 
# # create coefficient dataframe for forward modelling of comorbidities after T2DM
# 
# mod_coef_df <- data.frame(matrix(nrow = 0, ncol = 8))
# 
# comorbs <- unique(model_df$block)
# comorbs <- comorbs[!comorbs %in% c("Death", "CKD")]
# 
# # for each comorbidity event after T2DM
# for(i in 1:length(comorbs)) {
#   
#   cmb <- comorbs[i]
#   
#   # create test data.frame for comorbidity i
#   
#   test <- model_df %>%
#     # convert comorbidity i to  "comorb"
#     mutate(block = ifelse(block == cmb, "comorb", block)) %>%
#     # pivot wider
#     pivot_wider(names_from = block, values_from = event_dt) %>%
#     # select for relevant variables
#     select(eid, sex, age, contains("mca"), comorb, CKD, Death) %>%
#     # join with follow-up times
#     full_join(futimes, by = "eid") %>%
#     # convert follow-up times which = 0 to = 1
#     mutate(futime = ifelse(futime == 0, futime + 1, futime)) %>%
#     filter(sex == 1) %>%
#     select(-sex)
#   
#   if(sum(!is.na(test$comorb)) == 0) {
#     next
#   }
#   
#   # remove those where T2DM is last disease as there are no
#   # comorbidities to model after T2DM
#   # set comorbidities before T2DM to NA as we are controlling for
#   # the ones affecting CKD progression with the MCA coordinates
#   
#   test %<>%
#     filter(futime > 0) %>%
#     mutate(comorb = ifelse(comorb < 0, NA, comorb))
#   
#   # create survival input table
#   
#   test_new <- tmerge(data1 = select(test, contains(c("eid", "age", "mca"))), data2 = test, id = eid, tstop = futime)
#   # create CKD event
#   test_new <- tmerge(test_new, test, id = eid, CKD = event(CKD))
#   # create time-dependent CKD event
#   test_new <- tmerge(test_new, test, id = eid, comorb = tdc(comorb))
#   
#   # format variables
#   
#   test_new %<>%
#     mutate(age = as.numeric(age), comorb = factor(ifelse(comorb == 1, "cmb", "no_cmb"), levels = c("no_cmb", "cmb")))
#   
#   # run Cox PH modelling of time-dependent comorbidity development after T2DM on effect on CKD progression
#   # controlling for age at T2DM indication, sex, and previous comorbidity history 
#   # affecting CKD progression represented in MCA coordinates
#   
#   # create model formula
#   model_form <- test_new %>%
#     select(contains(c("age", "comorb", "mca"))) %>%
#     colnames()
#   model_form <- as.formula(paste0("Surv(tstart, tstop, CKD) ~ ", paste(model_form, collapse = " + ")))
#   # model_form <- as.formula(paste0(model_form, " + dis_bef:comorb")) # add sex interaction term
#   
#   # run model
#   
#   cmb_model <- coxph(model_form,
#                      data = test_new, cluster = eid)
#   
#   # create model summary
#   model_sum <- summary(cmb_model)
#   # get coefficients
#   model_sum <- data.frame(model_sum$coefficients)
#   # get confidence intervals of hazard ratios
#   cmb_confint <- confint(cmb_model)
#   colnames(cmb_confint) <- c("ci_2_5_pct", "ci_97_5_pct")
#   # create model summary output
#   model_sum <- data.frame(comorb = cmb, var = rownames(model_sum), model_sum, cmb_confint)
#   
#   # add model summary to coefficient df
#   mod_coef_df <- rbind(mod_coef_df, model_sum)
#   print(i)
# }
# 
# # adjust p-values
# 
# mod_coef_df %<>%
#   mutate(Pr...z.. = p.adjust(Pr...z.., method = "holm"))
# 
# mod_coef_df %>%
#   filter(var == "comorbcmb") %>%
#   filter(Pr...z.. < 0.05) %>%
#   arrange(desc(exp.coef.))
# 
# 
# 
# 
# 
# 
# 
# # multimorbidity at T2DM diagnosis
# 
# library(cluster)
# library(fpc)
# 
# # join primary care and secondary care blocks
# 
# diags <- rbind(icd10_filt, gp_clinic_block)
# 
# # filter for physiological blocks
# 
# diags %<>% filter(block %in% dis_use) %>%
#   mutate(eid = as.character(eid)) %>%
#   # remove redundant diabetes + CKD codes and unknown/uncertain/same as date of birth diagnoses
#   filter(block != "250.2" & block != "585.3" & event_dt != "1901-01-01" & event_dt != "1902-02-02" & event_dt != "1903-03-03" & !is.na(event_dt)) %>%
#   group_by(eid, block) %>%
#   # filter for earliest date of block across primary + secondary care diagnoses for each individual
#   filter(event_dt == min(event_dt, na.rm = T))
# 
# # remove redundant rows
# 
# diags <- unique.data.frame(diags)
# 
# # get earliest diagnosis times of T2DM and CKD
# 
# t2dm_indic_min <- t2dm_indic %>%
#   filter(t2dm_event_type %in% c("primary_care_diagnosis", "secondary_care_diagnosis")) %>%
#   mutate(eid = as.character(eid)) %>%
#   group_by(eid) %>%
#   # get earliest T2DM event date
#   filter(t2dm_event_dt == min(t2dm_event_dt, na.rm = T)) %>%
#   select(eid, t2dm_event_dt) %>%
#   unique.data.frame()
# 
# ckd_indic_min <- ckd_indic %>%
#   filter(ckd_event_type != "primary_egfr_persist_lt_60") %>%
#   mutate(eid = as.character(eid)) %>%
#   group_by(eid) %>%
#   # get earliest CKD event date
#   filter(ckd_event_dt == min(ckd_event_dt, na.rm = T)) %>%
#   select(eid, ckd_event_dt) %>%
#   unique.data.frame()
# 
# # join with T2DM + CKD earliest times with diagnosis data
# 
# diags %<>%
#   left_join(t2dm_indic_min, by = "eid") %>%
#   left_join(ckd_indic_min, by = "eid")
# 
# # create end point df
# # if individual develops CKD, this is their end point
# # if they do not, label them as "End"
# 
# ckd_event <- diags %>%
#   ungroup() %>%
#   select(eid, t2dm_event_dt, ckd_event_dt) %>%
#   unique.data.frame() %>%
#   mutate(block = ifelse(is.na(ckd_event_dt), "End", "CKD"), event_dt = ckd_event_dt) %>%
#   select(eid, event_dt, block, t2dm_event_dt, ckd_event_dt)
# 
# # join
# 
# diags %<>%
#   mutate(block = as.character(block))
# 
# diags <- rbind(diags, ckd_event)
# 
# # filter out diseases <= 10% prevalence
# 
# diags %<>%
#   filter(block %in% c(as.character(dis_use), "CKD", "End"))
# 
# # convert dates of diagnoses to relative time from T2DM indication
# 
# diags %<>%
#   ungroup() %>%
#   mutate(event_dt = as.Date(event_dt) - as.Date(t2dm_event_dt))
# 
# # get age at which T2DM indication occurred
# 
# # get age at baseline UK attendance
# age <- ukb_data %>%
#   transmute(eid = as.character(f.eid), age = f.21003.0.0)
# 
# # get sex
# 
# sex <- ukb_data %>%
#   transmute(eid = as.character(f.eid), sex = f.31.0.0)
# 
# # join with diagnoses df
# 
# # bl_dates %<>%
# #   transmute(eid = as.character(f.eid), bl_date = f.53.0.0)
# 
# diags %<>%
#   left_join(bl_dates, by = "eid")
# 
# # get time difference between date of attendance of assessment centre and date of T2DM indication
# 
# diags %<>%
#   # get difference
#   mutate(t2dm_diff = as.Date(t2dm_event_dt) - as.Date(bl_date)) %>%
#   # convert to nearest year
#   mutate(t2dm_diff = round(as.numeric(t2dm_diff/365)))
# 
# # join diagnoses with sex and age
# 
# diags %<>%
#   left_join(age, by = "eid") %>%
#   left_join(sex, by = "eid")
# 
# # get age at earliest T2DM indication
# 
# diags %<>%
#   mutate(t2dm_age = age + t2dm_diff)
# 
# # get death variable from UKB data
# 
# death <- ukb_data %>%
#   transmute(eid = as.character(f.eid), death = as.Date(f.40000.0.0))
# 
# # join diagnoses df with death
# 
# diags %<>%
#   left_join(death, by = "eid")
# 
# # NA those deaths after CKD
# 
# diags$death[diags$death > diags$ckd_event_dt] <- NA
# 
# # get relative time of death from T2DM indication
# 
# diags %<>%
#   mutate(death = death - t2dm_event_dt)
# 
# # add death event
# 
# tmp <- diags %>%
#   # filter for those who died
#   filter(!is.na(death)) %>%
#   # add death event
#   mutate(event_dt = death, block = "Death") %>%
#   unique.data.frame()
# 
# # add in death events
# 
# diags <- rbind(diags, tmp) %>%
#   arrange(eid, event_dt)
# 
# # remove End
# 
# diags %<>%
#   filter(block != "End")
# 
# # check again and remove extremely early age of T2DM diagnoses (same age cutoff as our exclusion criteria)
# 
# early_t2dm <- diags %>%
#   filter(t2dm_age < 36) %>%
#   pull(eid) %>%
#   unique()
# 
# diags %<>%
#   filter(!eid %in% early_t2dm)
# 
# # plot distribution time to CKD from T2DM
# 
# plot <- diags %>%
#   filter(block == "CKD") %>%
#   ggplot(aes(x = event_dt/365)) +
#   stat_ecdf() +
#   theme_bw() +
#   scale_x_continuous(breaks = seq(-35, 35, 5)) +
#   scale_y_continuous(breaks = seq(0, 1, 0.1)) +
#   xlab("Time of CKD (stage 3-5) diagnosis in years from T2DM diagnosis") +
#   ylab("Cumulative proportion")
# ggsave(filename = paste0(plot_dir, "t2dm_time_to_ckd.png"), plot = plot,
#        units = "px", width = 1800, height = 1600)
# 
# # ~95% of cases of CKD stage 1-3 happen within 15 years of T2DM diagnosis
# 
# # create multimorbidity matrix at T2DM diagnosis
# 
# # recreate model df
# 
# model_df <- diags %>%
#   select(eid, event_dt, block, age, sex) %>%
#   mutate(event_dt = as.numeric(event_dt))
# 
# 
# cmb_mm <- model_df %>%
#   pivot_wider(names_from = block, values_from = event_dt) %>%
#   ungroup() %>%
#   # convert comorbidity times after T2DM or which are NA to 0 and
#   # those comorbidity times before T2DM to 1
#   select(-all_of(c("eid", "age", "sex", "Death"))) %>%
#   mutate(across(everything(), .fns = function(x) ifelse(is.na(x) | x > 0, 0, 1))) %>%
#   as.matrix()
# 
# # convert to logical
# 
# cmb_mm <- apply(cmb_mm, 2, as.logical)
# 
# # run MCA
# 
# cmb_mca <- MCA(cmb_mm)
# 
# # plot scree
# 
# plot(cmb_mca$eig[,1])
# 
# # select number of components using eigenvalues of components > 1/nvar rule
# 
# ncp_use <- sum(cmb_mca$eig[,1] > 1/ncol(cmb_mm))
# 
# cmb_mca$eig[ncp_use, 3]
# # 39.65% of variance with 14
# 
# # re-reun MCA with number of components to retain
# 
# cmb_mca <- MCA(cmb_mm, ncp = ncp_use)
# 
# # get coordinates of individuals individual
# 
# indiv_coor <- as.data.frame(cmb_mca$ind$coord)
# 
# # plot
# 
# indiv_coor %>%
#   ggplot(aes(x = `Dim 1`, y = `Dim 2`, col = factor(cmb_mm[, which(colnames(cmb_mm) == "401.1")]))) +
#   geom_point(pch = ".")
# 
# # get coordinates of diseases
# 
# dis_coor <-as.data.frame(cmb_mca$var$coord)
# 
# # remove FALSE
# 
# dis_coor <- dis_coor[-which(str_detect(rownames(dis_coor), "FALSE")),]
# rownames(dis_coor) <- str_remove(rownames(dis_coor), "_TRUE")
# 
# # convert to names
# 
# tmp_dis <- phecode_icd10$Phenotype[match(rownames(dis_coor), phecode_icd10$PheCode)]
# tmp_dis[is.na(tmp_dis)] <- "Chronic renal failure [CKD]"
# 
# rownames(dis_coor) <- tmp_dis
# 
# # cluster diseases
# 
# hclust_model <- hclust(dist(dis_coor), method = "ward.D2")
# 
# fviz_dend(hclust_model, cex = 0.7, labels_track_height = 12)
# abline(h = 20)
# 
# # cluster individuals
# 
# hclust_model <- hclust(dist(indiv_coor), method = "ward.D2")
# 
# ch_k <- sapply(2:15, function(k) calinhara(indiv_coor, cutree(hclust_model, k)))
# 
# plot(ch_k, type = "b")
# 
# # cl <- cutree(hclust_model, h = 20)
# # 
# # # umap
# # 
# # library(umap)
# # 
# # umap_coor <- umap(indiv_coor, n_neighbors = 200)
# # 
# # # plot
# # 
# # umap_coor$layout %>%
# #   as.data.frame() %>%
# #   mutate(cl = cl) %>%
# #   ggplot(aes(x = V1, y = V2, col = factor(cl))) +
# #   geom_point() +
# #   theme_bw()
# # 
# # 
# # umap_coor$layout %>%
# #   as.data.frame() %>%
# #   mutate(cl = cl) %>%
# #   filter(V2 > -12.5 & V1 < 25 & V1 > -10 & V2 < 10) %>%
# #   ggplot(aes(x = V1, y = V2, col = factor(cl))) +
# #   geom_point(pch = ".") +
# #   theme_bw()
# #   
# 
# # multimorbidity at CKD diagnosis before T2DM
# 
# 
# # create multimorbidity matrix at CKD diagnosis
# 
# # recreate model df
# 
# model_df <- diags %>%
#   select(eid, event_dt, block, age, sex) %>%
#   mutate(event_dt = as.numeric(event_dt))
# 
# # get those with CKD before T2DM
# 
# ckd_bef <- model_df %>%
#   filter(block == "CKD") %>%
#   filter(event_dt < 0) %>%
#   pull(eid)
# 
# model_df %<>%
#   filter(eid %in% ckd_bef)
# 
# # get CKD time
# 
# ckd_bef_t <- model_df %>%
#   filter(block == "CKD") %>%
#   transmute(eid = eid, CKD_t = event_dt)
# 
# # convert diagnoses to relative to CKD diagnosis before T2DM
# 
# model_df %<>%
#   left_join(ckd_bef_t, by = "eid") %>%
#   mutate(event_dt = event_dt - CKD_t)
# 
# cmb_mm <- model_df %>%
#   pivot_wider(names_from = block, values_from = event_dt) %>%
#   ungroup() %>%
#   select(-CKD_t) %>%
#   # convert comorbidity times after T2DM or which are NA to 0 and
#   # those comorbidity times before T2DM to 1
#   select(-all_of(c("eid", "age", "sex"))) %>%
#   mutate(across(everything(), .fns = function(x) ifelse(is.na(x) | x > 0, 0, 1))) %>%
#   as.matrix()
# 
# # convert to logical
# 
# cmb_mm <- apply(cmb_mm, 2, as.logical)
# 
# # run MCA
# 
# cmb_mca <- MCA(cmb_mm)
# 
# # plot scree
# 
# plot(cmb_mca$eig[,1])
# 
# # select number of components using eigenvalues of components > 1/nvar rule
# 
# ncp_use <- sum(cmb_mca$eig[,1] > 1/ncol(cmb_mm))
# 
# cmb_mca$eig[ncp_use, 3]
# # 56.61% of variance with 21
# cmb_mca$eig[4, 3]
# # 17.9% with 4 based on scree plot
# 
# # re-reun MCA with number of components to retain
# 
# cmb_mca <- MCA(cmb_mm, ncp = ncp_use)
# 
# # get coordinates of individuals individual
# 
# indiv_coor <- as.data.frame(cmb_mca$ind$coord)
# 
# # plot
# 
# indiv_coor %>%
#   ggplot(aes(x = `Dim 1`, y = `Dim 2`)) +
#   geom_point(pch = "o")
# 
# # get coordinates of diseases
# 
# dis_coor <-as.data.frame(cmb_mca$var$coord)
# 
# # remove FALSE
# 
# dis_coor <- dis_coor[-which(str_detect(rownames(dis_coor), "FALSE")),]
# rownames(dis_coor) <- str_remove(rownames(dis_coor), "_TRUE")
# 
# # convert to names
# 
# tmp_dis <- phecode_icd10$Phenotype[match(rownames(dis_coor), phecode_icd10$PheCode)]
# tmp_dis[is.na(tmp_dis)] <- "Chronic renal failure [CKD]"
# 
# rownames(dis_coor) <- tmp_dis
# 
# # cluster diseases
# 
# hclust_model <- hclust(dist(dis_coor), method = "ward.D2")
# 
# fviz_dend(hclust_model, cex = 0.7, labels_track_height = 12)
# abline(h = 20)
# 
# # cluster individuals
# 
# hclust_model <- hclust(dist(indiv_coor), method = "ward.D2")
# 
# ch_k <- sapply(2:15, function(k) calinhara(indiv_coor, cutree(hclust_model, k)))
# 
# plot(ch_k, type = "b")
# 
# # disease prevalences at CKD before T2DM
# 
# prev_ckd_bef <- colSums(cmb_mm)/length(ckd_bef)
# 
# prev_ckd_bef <- data.frame(dis = names(prev_ckd_bef), prev = prev_ckd_bef)
# 
# prev_ckd_bef %<>%
#   mutate(name = phecode_icd10$Phenotype[match(dis, phecode_icd10$PheCode)])
# 
# prev_ckd_bef %>%
#   arrange(desc(prev))
# 
# prev_ckd_bef %>%
#   filter(!is.na(name)) %>%
#   mutate(name = factor(name, levels = name[order(prev)])) %>%
#   ggplot(aes(y = name, x = prev)) +
#   geom_col(fill = "lightgrey", col = "black") +
#   theme_bw()
# 
# # get age at CKD (rough)
# 
# model_df %<>%
#   mutate(ckd_age = age + round(CKD_t/365))
# 
# # get age at CKD diagnosis in non-T2DM
# 
# non_t2dm_df <- non_t2dm$model_df
# 
# non_t2dm_ckd_age <- non_t2dm_df %>% 
#   filter(block == "CKD") %>%
#   ungroup() %>%
#   transmute(eid = eid, CKD_age = round(event_dt/365))
# 
# # match up every pre-T2DM CKD with a non-T2DM CKD diagnosed at similar age
# 
# tmp <- non_t2dm_ckd_age
# 
# matched <- character()
# 
# set.seed(35)
# 
# eids <- unique(model_df$eid)
# for(i in 1:length(eids)) {
#   eid_i <- eids[i]
#   age_i <- model_df %>% filter(eid == eid_i) %>%
#     pull(ckd_age) %>%
#     unique()
#   match_i <- which(abs(tmp$CKD_age - age_i) < 5)
#   samp_i <- sample(match_i, 3)
#   eid_match <- tmp$eid[samp_i]
#   matched <- c(matched, eid_match)
#   tmp <- tmp[-samp_i, ]
# }
# 
# # filtered non-T2DM
# 
# non_t2dm_filt <- non_t2dm_df %>%
#   filter(eid %in% matched)
# 
# # sex differences
# 
# # non-T2DM
# non_t2dm_filt %>%
#   ungroup() %>%
#   select(eid, sex) %>%
#   unique.data.frame() %>%
#   count(sex)
# 1747/(1747+1235)
# 
# # pre-T2DM CKD
# model_df %>%
#   select(eid, sex) %>%
#   unique.data.frame() %>%
#   count(sex)
# 492/(492+502)
# 
# # T2DM all CKD
# diags %>% 
#   filter(block == "CKD") %>% 
#   select(eid, sex) %>% 
#   unique.data.frame() %>% 
#   count(sex)
# 1417/(1417+1903)
# 
# # sex differences; more females in non-T2DM group
# # pre-T2DM CKD is intermediate sex composition between post-T2DM CKD and non-T2DM CKD
# 
# # we should probably match on sex as well in order to compare multimorbidity composition
# 
# # get sex of non-T2DM
# 
# non_t2dm_ckd_sex <- non_t2dm_df %>% 
#   filter(block == "CKD") %>%
#   ungroup() %>%
#   select(eid, sex)
# 
# tmp <- non_t2dm_ckd_age
# tmp2 <- non_t2dm_ckd_sex
# 
# matched <- character()
# 
# set.seed(35)
# 
# eids <- unique(model_df$eid)
# for(i in 1:length(eids)) {
#   eid_i <- eids[i]
#   age_i <- model_df %>% filter(eid == eid_i) %>%
#     pull(ckd_age) %>%
#     unique()
#   sex_i <- model_df %>% filter(eid == eid_i) %>%
#     pull(sex) %>%
#     unique()
#   match_i <- which(abs(tmp$CKD_age - age_i) < 5 & tmp2$sex == sex_i)
#   samp_i <- sample(match_i, 3)
#   eid_match <- tmp$eid[samp_i]
#   matched <- c(matched, eid_match)
#   tmp <- tmp[-samp_i, ]
#   tmp2 <- tmp2[-samp_i,]
# }
# 
# # filtered non-T2DM
# 
# non_t2dm_filt <- non_t2dm_df %>%
#   filter(eid %in% matched)
# 
# # compare multimorbidity at CKD diagnosis
# 
# # convert non-T2DM times to time since CKD diagnosis
# 
# non_t2dm_ckd_t <- non_t2dm_df %>% 
#   filter(eid %in% matched) %>%
#   filter(block == "CKD") %>%
#   ungroup() %>%
#   transmute(eid = eid, CKD_age = event_dt)
# 
# non_t2dm_filt %<>%
#   left_join(non_t2dm_ckd_t, "eid")
# 
# non_t2dm_filt %<>%
#   mutate(event_dt = event_dt - CKD_age) %>%
#   select(-CKD_age)
# 
# # create non-T2DM mm
# 
# non_t2dm_mm <- non_t2dm_filt %>%
#   filter(block %in% unique(diags$block)) %>%
#   filter(block != "Death") %>%
#   pivot_wider(names_from = block, values_from = event_dt) %>%
#   ungroup() %>%
#   # convert comorbidity times after T2DM or which are NA to 0 and
#   # those comorbidity times before T2DM to 1
#   select(-all_of(c("eid", "age", "sex"))) %>%
#   mutate(across(everything(), .fns = function(x) ifelse(is.na(x) | x > 0, 0, 1))) %>%
#   as.matrix()
# 
# # convert to logical
# 
# non_t2dm_mm <- apply(non_t2dm_mm, 2, as.logical)
# 
# # run MCA
# 
# non_t2dm_mca <- MCA(non_t2dm_mm)
# 
# # plot scree
# 
# plot(non_t2dm_mca$eig[,1])
# 
# # select number of components using eigenvalues of components > 1/nvar rule
# 
# ncp_use <- sum(non_t2dm_mca$eig[,1] > 1/ncol(non_t2dm_mm))
# 
# non_t2dm_mca$eig[ncp_use, 3]
# # 50.19568% of variance with 20
# 
# # re-reun MCA with number of components to retain
# 
# non_t2dm_mca <- MCA(non_t2dm_mm, ncp = ncp_use)
# 
# # get coordinates of individuals individual
# 
# non_indiv_coor <- as.data.frame(non_t2dm_mca$ind$coord)
# 
# # plot
# 
# non_indiv_coor %>%
#   ggplot(aes(x = `Dim 1`, y = `Dim 2`)) +
#   geom_point(pch = "o")
# 
# # get coordinates of diseases
# 
# non_dis_coor <-as.data.frame(non_t2dm_mca$var$coord)
# 
# # remove FALSE
# 
# non_dis_coor <- non_dis_coor[-which(str_detect(rownames(non_dis_coor), "FALSE")),]
# rownames(non_dis_coor) <- str_remove(rownames(non_dis_coor), "_TRUE")
# 
# # convert to names
# 
# tmp_non_dis <- phecode_icd10$Phenotype[match(rownames(non_dis_coor), phecode_icd10$PheCode)]
# tmp_non_dis[is.na(tmp_non_dis)] <- "Chronic renal failure [CKD]"
# 
# rownames(non_dis_coor) <- tmp_non_dis
# 
# # cluster non_diseases
# 
# hclust_model <- hclust(dist(non_dis_coor), method = "ward.D2")
# 
# fviz_dend(hclust_model, cex = 0.7, labels_track_height = 12)
# 
# # get disease prevalences
# 
# # disease prevalences at CKD in non-T2DM at comparable age
# 
# prev_non_t2dm_ckd <- colSums(non_t2dm_mm)/nrow(non_t2dm_mm)
# 
# prev_non_t2dm_ckd <- data.frame(dis = names(prev_non_t2dm_ckd), prev = prev_non_t2dm_ckd)
# 
# prev_non_t2dm_ckd %<>%
#   mutate(name = phecode_icd10$Phenotype[match(dis, phecode_icd10$PheCode)])
# 
# prev_non_t2dm_ckd %>%
#   arrange(desc(prev))
# 
# prev_non_t2dm_ckd %>%
#   filter(!is.na(name)) %>%
#   mutate(name = factor(name, levels = name[order(prev)])) %>%
#   ggplot(aes(y = name, x = prev)) +
#   geom_col(fill = "lightgrey", col = "black") +
#   theme_bw()
# 
# # compare
# 
# prev_non_t2dm_ckd %>%
#   filter(!is.na(name)) %>%
#   full_join(prev_ckd_bef %>% filter(!is.na(name)), by = c("dis", "name"),
#             suffix = c("_non_t2dm", "_t2dm")) %>%
#   mutate(diff = prev_t2dm - prev_non_t2dm) %>%
#   mutate(name = factor(name, levels = name[order(diff)])) %>%
#   ggplot(aes(x = diff, y = name)) +
#   geom_col(fill = "lightgrey", col = "black") +
#   theme_bw()
# 
# # compare multimorbidity
# 
# ckd_bef_morb <- rowSums(cmb_mm)
# non_t2dm_morb <- rowSums(non_t2dm_mm)
# 
# data.frame(group = c(rep("pre_t2dm_ckd", length(ckd_bef_morb)), rep("non_t2dm_ckd", length(non_t2dm_morb))),
#            morb = c(ckd_bef_morb, non_t2dm_morb)) %>%
#   ggplot(aes(x = group, y = morb)) +
#   geom_boxplot() +
#   theme_bw() +
#   xlab("") +
#   ylab("Morbidity at CKD diagnosis") +
#   scale_x_discrete(labels = c(non_t2dm_ckd = "Age + sex matched\nnon-T2DM CKD", pre_t2dm_ckd = "Pre-T2DM CKD")) +
#   scale_y_continuous(breaks = seq(0, 25, 2))
# 
# 
# # multimorbidity at T2DM diagnosis of pre-T2DM CKD patients
# 
# # create multimorbidity matrix at CKD diagnosis
# 
# # recreate model df
# 
# model_df <- diags %>%
#   select(eid, event_dt, block, age, sex) %>%
#   mutate(event_dt = as.numeric(event_dt))
# 
# # get those with CKD before T2DM
# 
# ckd_bef <- model_df %>%
#   filter(block == "CKD") %>%
#   filter(event_dt < 0) %>%
#   pull(eid)
# 
# # create variable indicating pre-T2DM CKD
# 
# model_df$ckd_bef <- "no_ckd_bef"
# model_df$ckd_bef[model_df$eid %in% ckd_bef] <- "ckd_bef"
# 
# tmp <- model_df %>%
#   pivot_wider(names_from = block, values_from = event_dt) %>%
#   select(eid, ckd_bef)
# 
# # create comorbidity at T2DM matrix
# 
# cmb_mm <- model_df %>%
#   pivot_wider(names_from = block, values_from = event_dt) %>%
#   ungroup() %>%
#   # convert comorbidity times after T2DM or which are NA to 0 and
#   # those comorbidity times before T2DM to 1
#   select(-all_of(c("eid", "age", "sex", "ckd_bef"))) %>%
#   mutate(across(everything(), .fns = function(x) ifelse(is.na(x) | x > 0, 0, 1))) %>%
#   as.matrix()
# 
# # convert to logical
# 
# cmb_mm <- apply(cmb_mm, 2, as.logical)
# 
# morb_at_t2dm <- rowSums(cmb_mm)
# 
# data.frame(group = tmp$ckd_bef, morb = morb_at_t2dm) %>%
#   ggplot(aes(x = group, y = morb)) +
#   geom_boxplot() +
#   theme_bw() +
#   xlab("") +
#   ylab("Morbidity at T2DM diagnosis") +
#   scale_x_discrete(labels = c(no_ckd_bef = "No CKD before T2DM", "1" = "Pre-T2DM CKD")) +
#   scale_y_continuous(breaks = seq(0, 34, 2))
# 
# data.frame(group = tmp$ckd_bef, morb = morb_at_t2dm) %>%
#   group_by(group) %>%
#   summarise(sum(morb > 3))
