# create replication cohort for comorbidity + CKD GWAS in non-T2DM population

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
library(parallel)
library(foreach)
library(doParallel)

work_dir = ""

read_dict_file = paste0(work_dir, "/", "inputs/all_lkps_maps_v3.txt") # dictionary of read codes
phecode_icd10_file = paste0(work_dir, "/", "inputs/phecode_icd10.csv") # mapping of ICD10 bodes to phecodes
read_2_icd10_file = paste0(work_dir, "/", "inputs/read_2_icd10.tsv") # mapping of Read 2 to ICD10
read_3_icd10_file = paste0(work_dir, "/", "inputs/read_3_icd10.tsv") # mapping of Read 3 to ICD10
ukb_data_file = "" # UKB file
gp_egfr_file = "" # UKB GP eGFR data
# gp_bmi_file = "/nfs_home/projects/departments/nnrco/comp_bio/ml_obesity/data_extractions/clinical/biomarkers/primary_care/combined_codes_bmi.txt"
gp_clinic_file = "" # UKB GP clinical data
ckd_pheno_file = "" # CKD phenotypic codes for primary + secondary care
t2dm_pheno_file = "" # T2DM phenotypic codes for primary + secondary care
primary_hba1c_file = "" # UKB GP HbA1c data
# egfr_persist_file = paste0(work_dir, "/", "ckd_prim_egfr_indic.rds")
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

# read in phecode/ICD10 mapping

phecode_icd10 <- fread(phecode_icd10_file)

# filter for 3 letter code

phecode_icd10 %<>%
  filter(str_count(ICD10) == 3)

# read in UKB GP clinical data

gp_clinic <- fread(gp_clinic_file)

# read in UKB GP eGFR data

gp_egfr <- fread(gp_egfr_file)

# read in UKB GP BMI data

gp_bmi <- fread(gp_bmi_file)

# read in read dictionary

read_dict <- fread(read_dict_file)

# read in T2DM EID inclusion/exclusion criteria file

t2dm_eids <- fread(t2dm_eids_file)
# needs redone

# filter for those who do not have evidence of any diabetes

non_t2dm_eids <- t2dm_eids %>%
  filter(!self_t1dm & !self_diab & !diab_medic & !self_6153_insulin & !hba1c & !primary_care_t1dm &!diab_age_lt_50 &
           !diab_age_lteq_30 & !diab_age_lteq_36 & !icd10_t1dm & !icd10_unspec_dm & !self_t2dm & !self_gest & !diab_diag &
           !diab_atc & !self_6177_insulin & !primary_care_hba1c & !primary_care_t2dm & !insulin_lt1yr & !icd10_t2dm &
           !icd10_other_spec_dm)

# filter UKB data for non-T2DM EIDs

ukb_data %<>%
  filter(f.eid %in% non_t2dm_eids$eid)

# filter UKB GP clinical for non-T2DM EIDs

gp_clinic %<>%
  filter(eid %in% non_t2dm_eids$eid)

# filter UKB GP clinical for non-T2DM EIDs

gp_egfr %<>%
  filter(eid %in% non_t2dm_eids$eid)

# get self-reported non-cancer diagnoses from UKB file

self_diag <- ukb_data %>%
  select(contains("f.20002"))

# get self-reported time of non-cancer diagnosis from UKB file

self_dates <- ukb_data %>%
  select(contains("f.20008"))

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
# there is extreme outliers
# one with CKD at birth

ckd_prim_indic %<>%
  filter(event_dt > "1920-01-01")

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
sum(lengths(ckd_t) > 1) # 2329 have more than one CKD diagnosis

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
  transmute(eid = f.eid, bl_date = `f.53.0.0`)


# use just secondary diagnosis indicators

ckd_indic <- rbind(ckd_second_indic)

# how many cases of CKD do we have with just secondary care?

ckd_indic %>%
  group_by(eid) %>%
  filter(event_dt == min(event_dt)) %>%
  pull(eid) %>%
  length()
# 11037 cases

colnames(ckd_indic) <- c("eid", "ckd_event_type", "ckd_event_dt")

# get earliest indication of CKD for those with both

ckd_indic %<>%
  # remove those without CKD event
  filter(!is.na(ckd_event_type)) %>%
  # group by individual
  group_by(eid) %>%
  # filter for earliest indication of CKD
  filter(ckd_event_dt == min(ckd_event_dt))

# map read codes to ICD10 blocks

read_icd10 <- read_2_icd10 %>%
  transmute(read_2_code = coding, meaning = meaning) %>%
  full_join(read_3_icd10, by = "meaning") %>%
  transmute(read_2_code = read_2_code, read_3_code = coding, icd10_code = meaning)

# selecting for diseases with >= 10% prevalence within CKD

# most prevalent ICD10 secondary care codes

# filter ICD10 data for CKD subgroup
icd10_filt <- as.matrix(icd10[ukb_data$f.eid %in% ckd_indic$eid,])

# convert to 3 letter code

icd10_filt <- str_trunc(icd10_filt, 3, ellipsis = "")

# get instances of ICD10 blocks in secondary care for each individual
icd10_filt %<>%
  as.data.frame() %>%
  # filter for CKD subset
  mutate(eid = ukb_data$f.eid[ukb_data$f.eid %in% ckd_indic$eid]) %>%
  # filter for individuals who do not have GP data
  filter(!eid %in% unique(gp_clinic$eid)) %>%
  pivot_longer(-eid, names_to = "diag", values_to = "block") %>%
  select(eid, block) %>%
  unique.data.frame()

# convert ICD10 codes to phecodes

icd10_filt %<>%
  filter(str_detect(block, "A|B|C|D|E|F|G|H|I|J|K|L|M|N")) %>%
  mutate(block = phecode_icd10$PheCode[match(block, phecode_icd10$ICD10)]) %>%
  select(eid, block)

# combine primary + secondary care ICD10 block instances

diag_prev <- rbind(icd10_filt) %>%
  # remove duplicate ICD10 blocks
  unique.data.frame()

# get counts
diag_prev %<>%
  pull(block) %>%
  table() %>%
  as.data.frame()

colnames(diag_prev) <- c("block", "n")

# calculate frequencies
diag_prev %<>%
  mutate(freq = n/length(unique(icd10_filt$eid)),
         name = phecode_icd10$Phenotype[match(block, phecode_icd10$PheCode)])

# filter for those diseases which were used in GWAS

diag_prev_filt <- diag_prev %>%
  arrange(desc(freq)) %>%
  filter(block %in% c("274.1", "278.1", "280.1", "285",
                      "327", "366", "401.1", "411.2",
                      "411.3", "411.8", "427.2", "428.2",
                      "443.8", "458.9", "480", "507",
                      "535", "569", "585.1"))
# get diseases which were used in GWAS

dis_use <- diag_prev_filt %>%
  pull(block)

# create block data

icd10_filt <- as.matrix(icd10)

# convert to 3 letter code

icd10_filt <- str_trunc(icd10_filt, 3, ellipsis = "")

# create df

icd10_filt %<>%
  as.data.frame() %>%
  # add eids
  mutate(eid = ukb_data$f.eid) %>%
  # select for those without GP data
  filter(!eid %in% unique(gp_clinic$eid)) %>%
  pivot_longer(-eid, names_to = "diag", values_to = "block")

# get diagnosis order
icd10_filt$diag <- sapply(str_split(icd10_filt$diag, "\\."), function(x) tail(x, 1))

# get dates of diagnoses

icd10_t_filt <- as.matrix(icd10_dates)

# create df

icd10_t_filt %<>%
  as.data.frame() %>%
  # add EIDs
  mutate(eid = ukb_data$f.eid) %>%
  # filter for those without GP data
  filter(!eid %in% unique(gp_clinic$eid)) %>%
  pivot_longer(-eid, names_to = "diag", values_to = "event_dt")

# get diagnosis order

icd10_t_filt$diag <- sapply(str_split(icd10_t_filt$diag, "\\."), function(x) tail(x, 1))

# combine diagnoses and diagnosis dates

icd10_filt <- data.frame(icd10_filt, event_dt = icd10_t_filt$event_dt)

# filter out NAs
icd10_filt %<>%
  filter(!is.na(block) & !is.na(event_dt))

# convert ICD10 codes to phecodes

icd10_filt %<>%
  mutate(block = phecode_icd10$PheCode[match(block, phecode_icd10$ICD10)]) %>%
  select(eid, event_dt, block)

# format dates 

icd10_filt %<>%
  mutate(event_dt = as.Date(event_dt))

# use just secondary care blocks

diags <- rbind(icd10_filt)

# filter out weird times

diags %<>%
  mutate(eid = as.character(eid)) %>%
  # remove unknown/uncertain/same as date of birth diagnoses
  filter(event_dt != "1901-01-01" & event_dt != "1902-02-02" & event_dt != "1903-03-03" & !is.na(event_dt)) %>%
  group_by(eid, block) %>%
  # filter for earliest date of block across primary + secondary care diagnoses for each individual
  filter(event_dt == min(event_dt, na.rm = T))

# remove redundant rows

diags <- unique.data.frame(diags)

# get earliest indication times of CKD

ckd_indic_min <- ckd_indic %>%
  mutate(eid = as.character(eid)) %>%
  group_by(eid) %>%
  # get earliest CKD event date
  filter(ckd_event_dt == min(ckd_event_dt, na.rm = T)) %>%
  select(eid, ckd_event_dt) %>%
  unique.data.frame()

# join with CKD earliest times with diagnosis data

diags %<>%
  left_join(ckd_indic_min, by = "eid")

# create end point df
# if individual develops CKD, this is their end point
# if they do not, label them as "End"

ckd_event <- diags %>%
  ungroup() %>%
  select(eid, ckd_event_dt) %>%
  unique.data.frame() %>%
  mutate(block = ifelse(is.na(ckd_event_dt), "End", "CKD"), event_dt = ckd_event_dt) %>%
  select(eid, event_dt, block, ckd_event_dt)

# join

diags <- rbind(diags %>% mutate(block = as.character(block)), ckd_event)

# filter for diseases used in GWAS

diags %<>%
  filter(block %in% c(as.character(dis_use), "CKD", "End"))

# filter out diagnoses which occur after CKD diagnosis for those that develop CKD

# diags %<>%
#   filter(!event_dt > ckd_event_dt | is.na(ckd_event_dt))

# get age at baseline UK attendance

age <- ukb_data %>%
  transmute(eid = as.character(f.eid), age = f.21003.0.0)

# get sex

sex <- ukb_data %>%
  transmute(eid = as.character(f.eid), sex = f.31.0.0)

# join with diagnoses df

bl_dates %<>%
  transmute(eid = as.character(eid), bl_date = bl_date)

diags %<>%
  left_join(bl_dates, by = "eid")

# get time difference between date of attendance of assessment centre and diagnoses

diags <- diags %>%
  # get difference
  mutate(diag_diff = as.numeric(event_dt) - as.numeric(bl_date))

# join diagnoses with sex and age

diags %<>%
  left_join(age, by = "eid") %>%
  left_join(sex, by = "eid")

# get death variable from UKB data

death <- ukb_data %>%
  transmute(eid = as.character(f.eid), death = as.Date(f.40000.0.0))

# join diagnoses df with death

diags %<>%
  left_join(death, by = "eid")

# NA those deaths after CKD

diags$death[diags$death > diags$ckd_event_dt] <- NA

# convert age to days since birth (roughly)

diags %<>%
  mutate(age = round((age-1)*((365.5))+ 182.75)) %>%
  mutate(event_dt = age + diag_diff)

# get relative time of death from birth

diags %<>%
  mutate(death = as.numeric(death) - as.numeric(bl_date)) %>%
  mutate(death = death + age)

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

# plot distribution time to CKD from birth

diags %>%
  filter(block == "CKD") %>%
  ggplot(aes(x = event_dt/365)) +
  stat_ecdf() +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 90, 5)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1))

# set those events with negative time to t = 0

diags %<>%
  mutate(event_dt = ifelse(event_dt < 0, 0, event_dt))

# write out

df_out <- diags %>%
  select(eid, event_dt, block, age, sex) %>%
  mutate(event_dt = as.numeric(event_dt)) %>%
  unique.data.frame()

write_tsv(df_out, paste0(work_dir, "/intermediates/replicate_non_t2dm_time_to_ckd_out.txt"))
