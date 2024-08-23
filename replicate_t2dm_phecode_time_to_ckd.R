# create replication cohort for GWAS on comorbidity + CKD subgroups

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

# just using secondary care diagnosis indicators

ckd_indic <- ckd_second_indic

# how many cases of CKD from just secondary care do we have?

ckd_indic %>%
  group_by(eid) %>%
  filter(event_dt == min(event_dt)) %>%
  pull(eid) %>%
  length()
# 4470 cases

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

sum(lengths(self_t2dm_t) > 1) # 89 with multiple self-reported T2DM

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

# use just baseline asssessment centre HbA1c values

primary_hba1c <- ukb_hba1c

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

# combine assessment centre HbA1c, secondary care diagnosis,
# self-reported T2DM diagnosis, and self-reported diabetes diagnosis indications

t2dm_indic <- rbind(t2dm_prim_hba1c_indic %>% mutate(event_dt = as.Date(event_dt)),
                    t2dm_second_indic, self_t2dm_indic, self_diab_indic)

# get number of T2DM cases

t2dm_indic %>%
  group_by(eid) %>%
  filter(event_dt == min(event_dt)) %>%
  pull(eid) %>%
  length()
# 38665

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
  scale_x_discrete(labels = c(primary_egfr_persist_lt_60 = "Primary care or\nassessment centre HbA1c\nmeasurement > 48 mmol/mol",
                              secondary_care_diagnosis = "Secondary care\ndiagnosis of T2DM",
                              self_reported_t2dm = "Self-reported T2DM\nat assessment centre",
                              self_reported_diabetes = "Self-reported Diabetes\nat assessment centre")) +
  scale_y_continuous(breaks = seq(0, 14000, 2000))
# relatively even spread, top 2 are self-reported diabetes and secondary care diagnosis of T2DM


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

# use just secondary care ICD10 block instances

diag_prev <- rbind(icd10_filt) %>%
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

# filter for those diseases which were run in the GWAS

diag_prev_filt <- diag_prev %>%
  arrange(desc(freq)) %>%
  filter(block %in% c("274.1", "278.1", "280.1", "285",
                      "327", "366", "401.1", "411.2",
                      "411.3", "411.8", "427.2", "428.2",
                      "443.8", "458.9", "480", "507",
                      "535", "569", "585.1"))

# get diseases used in GWAS

dis_use <- diag_prev_filt %>%
  pull(block)

# create T2DM block data

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
  filter(!eid %in% unique(gp_clinic$eid)) %>%
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
  filter(!eid %in% unique(gp_clinic$eid)) %>%
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

# format dates

icd10_filt %<>%
  mutate(event_dt = as.Date(event_dt))

# use just secondary care blocks

diags <- rbind(icd10_filt)

# filter for physiological blocks

diags %<>% filter(block %in% dis_use) %>%
  mutate(eid = as.character(eid)) %>%
  # remove redundant diabetes + CKD codes and unknown/uncertain/same as date of birth diagnoses
  filter(block != "250.2" & block != "585.3" & event_dt != "1901-01-01" & event_dt != "1902-02-02" & event_dt != "1903-03-03" & !is.na(event_dt)) %>%
  group_by(eid, block) %>%
  # filter for earliest date of block across secondary care diagnoses for each individual
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

# filter out diseases not used in GWAS

diags %<>%
  filter(block %in% c(as.character(dis_use), "CKD", "End"))

# filter out diagnoses which occur after CKD diagnosis for those that develop CKD

# diags %<>%
#   filter(!event_dt > ckd_event_dt | is.na(ckd_event_dt))

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
# 521 individuals with CKD -> T2DM

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
  scale_x_continuous(breaks = seq(-20, 35, 5)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1))

# write out df

df_out <- diags %>%
  filter(event_dt >= 0)

df_out <- diags

write_tsv(df_out, paste0(work_dir, "/intermediates/replicate_t2dm_time_to_ckd_out.txt"))
