
# Generates T2DM criteria eIDs files

# Load libraries...

library(readr)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(data.table)

ukb_data_file = ""# ukb data file
ukb_insulin_lt_yr_file = "ukb_insulin_lt_1yr.tsv"
gp_clinic_file = "" # GP clinical data file
t2dm_pheno_file = "" # Type 2 diabetes phenotypic codes
t1dm_pheno_file = "" # Type 1 diabetes phenotypic codes
read_2_icd10_file = "read_2_icd10.tsv" # mapping of Read 2 to ICD10
read_3_icd10_file = "read_3_icd10.tsv" # mapping of Read 3 to ICD10
primary_hba1c_file = "combined_codes_HbA1c.txt"# GP primary care HbA1c values
medication_codes_file = "" # top most prevalent medication codes indicative of T2DM
ukb_medic_atc_map_file = "atc_all_matches.csv" # ATC mapping of UKB medication codes
atc_onto_file = "ATC.csv" # ATC ontology table
t2dm_crit_eids_df_file = "t2dm_criteria_eids.tsv" # output criteria EIDs

ukb_cols_use <- c("f.eid", "f.20002", "f.2443", "f.20003", "f.6153", "f.6177", "f.30750", "f.2976", "f.21000", "f.4041",
                  "f.41270")

# read in column names of UKB

ukb_cols <- fread(ukb_data_file, nrows = 0)
ukb_cols %<>%
  select(contains(ukb_cols_use)) %>%
  colnames()

# Load data...

ukb_data <- fread(ukb_data_file, select = ukb_cols)

# read in UKB data containing field code f.2986 insulin started within less than 1 year of diagnosis

insulin_lt_yr_data <- fread(ukb_insulin_lt_yr_file)

# left join insulin data

ukb_data <- ukb_data %>%
  left_join(insulin_lt_yr_data, by = "f.eid")

# get eIDs

eids <- ukb_data$f.eid

# read in GP data

gp_clinic <- fread(gp_clinic_file)

# read in read 2/3 icd10 mappings

read_2_icd10 <- fread(read_2_icd10_file)
read_3_icd10 <- fread(read_3_icd10_file)

#
##
### GP T2DM diagnoses workflow
##
#

# get T2DM diagnoses phenotype file

t2dm_pheno <- fread(t2dm_pheno_file)

t2dm_primary <- t2dm_pheno %>%
  filter(str_detect(source, "Primary"))

# get those with T2DM primary care diagnosis

t2dm_primary_eid <- t2dm_primary %>%
  pull(eid) %>%
  unique()

# use read diagnostic codes from Read/ICD10 mapping to catch any missing

read_2_t2dm <- read_2_icd10 %>%
  filter(meaning == "E11") %>%
  pull(coding)

read_3_t2dm <- read_3_icd10 %>%
  filter(meaning == "E11") %>%
  pull(coding)

t2dm_primary_eid2 <- gp_clinic %>%
  filter(read_2 %in% read_2_t2dm | read_3 %in% read_3_t2dm) %>%
  pull(eid) %>%
  unique()

# how many more?
sum(!t2dm_primary_eid2 %in% t2dm_primary_eid)
# 6 additional

# create unique set

t2dm_primary_eid <- unique(t2dm_primary_eid, t2dm_primary_eid2)

#
##
### GP T1DM diagnoses workflow
##
#

# get T2DM diagnoses phenotype file

t1dm_pheno <- fread(t1dm_pheno_file)

t1dm_primary <- t1dm_pheno %>%
  filter(code_type == "read")

# get those with t1dm primary care diagnosis

t1dm_primary_eid <- gp_clinic %>%
  filter(read_2 %in% t1dm_primary$code | read_3 %in% t1dm_primary$code) %>%
  pull(eid) %>%
  unique()

read_2_t1dm <- read_2_icd10 %>%
  filter(meaning == "E10") %>%
  pull(coding)

read_3_t1dm <- read_3_icd10 %>%
  filter(meaning == "E10") %>%
  pull(coding)

t1dm_primary_eid2 <- gp_clinic %>%
  filter(read_2 %in% read_2_t1dm | read_3 %in% read_3_t1dm) %>%
  pull(eid) %>%
  unique()

# how many more?
sum(!t1dm_primary_eid2 %in% t1dm_primary_eid)
# 33 additional

# create unique set

t1dm_primary_eid <- unique(t1dm_primary_eid, t1dm_primary_eid2)



#
##
### GP HbA1c value workflow
##
#

# read in primary care HbA1c values

primary_hba1c <- fread(primary_hba1c_file)

# filter for those with HbA1c >= 48 mmol/mol

primary_hba1c_eid <- primary_hba1c %>%
  filter(value >= 48) %>%
  pull(eid) %>%
  unique()


#
##
### Self-reported diabetes workflow 
##
#

# Filter for self-reported T2DM diabetes column instance 0

ukb_data_self <- ukb_data %>%
  select(contains("f.20002."))
                                            

# Get the self-reported T2DM diabetes participants indices

self_t2dm_indices <- which(apply(ukb_data_self, 1, function(x) 1223 %in% x))

self_t2dm_eids <- eids[self_t2dm_indices]

# get the self-reported T1DM diabetes participants indices

self_t1dm_indices <- which(apply(ukb_data_self, 1, function(x) 1222 %in% x))

self_t1dm_eids <- eids[self_t1dm_indices]

# get the self-reported gestational diabetes participants indices

self_gest_indices <- which(apply(ukb_data_self, 1, function(x) 1221 %in% x))

self_gest_eids <- eids[self_gest_indices]

# get the self-reported diabetes participants indices

self_diab_indices <- which(apply(ukb_data_self, 1, function(x) 1220 %in% x))

self_diab_eids <- eids[self_diab_indices]



#
##
### Diabetes diagnosed by doctor (interview) workflow 
##
#


# Filter for diabetes diagnosed by doctor column instance 0

ukb_data_diag <- ukb_data %>%
  select(contains("f.2443"))

# Get the diabetes diagnosed by doctor T2DM participant indices

diabetes_diag_indices <- which(apply(ukb_data_diag, 1, function(x) 1 %in% x))

diabetes_diag_eids <- eids[diabetes_diag_indices]


#
##
### Diabetes from Treatment/medication code workflow 
##
#

# read in medication codes mapping

medic_codes <- read_tsv(medication_codes_file)

# read in UKB ATC mapping

ukb_atc <- read_csv(ukb_medic_atc_map_file, col_names = F)

colnames(ukb_atc) <- c("ukb_code", "ukb_name", "atc_code")

# read in atc ontology

atc_onto <- read_csv(atc_onto_file)

# parse class ID

atc_onto <- atc_onto %>%
  mutate(`Class ID` = str_remove(`Class ID`, "..*ATC/"))

atc_onto <- atc_onto %>%
  select(`Class ID`, `Preferred Label`)

colnames(atc_onto) <- c("atc_code", "atc_name")

# get ukb codes to use

# top most prevalent (>1% prevalence) UKB medication codes related to T2DM
medic_use <- c("1140883066", "1140874744", "1141171646", "1141177600", "1141152590",
               "1141189090", "1141189094", "1141168660", "1140874646", "1140910566", "1141171652", "1140874746", "1141177606",
               "1140884600", "1140874686", "1140868902", "1141156984", "1140874674", "1141173882")

medic_use_names <- c("insulin product", "gliclazide", "pioglitazone", "rosiglitazone", "glimepiride", "rosiglitazone 1mg / metformin 500mg tablet",
                     "avandamet 1mg / 500mg tablet", "repaglinide", "glipizide", "glyclizide", "actos 15mg tablet", "diamicron 80mg tablet", "avandia 4mg tablet",
                     "metformin", "glucophage 500mg tablet", "acarbose", "amaryl 1mg tablet", "tolbutamide", "nateglinide")

medic_use_df <- data.frame(medic_use, medic_use_names)

# Filter for treatment/medication column instance 0

ukb_data_medi <- ukb_data %>%
  select(contains("f.20003."))

# Get the diabetes-related medication recipient indices

medi_indices <- which(apply(ukb_data_medi, 1, function(x) any(medic_use %in% x)))

diabetes_medi_eids <- eids[medi_indices]

# ATC selection

# Filter for treatment/medication column instance 0

ukb_data_medi <- as.matrix(ukb_data_medi)

# convert to atc

ukb_data_atc <- apply(ukb_data_medi, 2, function(x) ukb_atc$atc_code[match(x, ukb_atc$ukb_code)])

diabetes_atc <- atc_onto %>% filter(str_starts(atc_code, "A10")) %>% pull(atc_code)

ukb_data_atc <- as.data.frame(ukb_data_atc)

# Get the diabetes-related medication recipient indices

atc_indices <- which(apply(ukb_data_atc, 1, function(x) any(diabetes_atc %in% x)))

diabetes_atc_eids <- eids[atc_indices]


#
##
### Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones workflow  
##
#


# Filter for field 6153 (medication for chol...) column instance 0

ukb_data_6153 <- ukb_data %>%
  select(contains("f.6153"))

# Get the field 6153 insulin therapy recipient indices

ukb_6153_indices <- which(apply(ukb_data_6153, 1, function(x) 3 %in% x))

self_6153_insulin_eids <- eids[ukb_6153_indices]


#
##
### Medication for cholesterol, blood pressure, diabetes workflow  
##
#


# Filter for field 6177 (medication for chol...) column instance 0

ukb_data_6177 <- ukb_data %>%
  select(contains("f.6177"))

# Get the field 6177 insulin therapy recipient indices

ukb_6177_indices <- which(apply(ukb_data_6177, 1, function(x) 3 %in% x))

self_6177_insulin_eids <- eids[ukb_6177_indices]


#
##
### Glycated haemoglobin workflow 
##
#

# Filter for field 30750 (HbA1c) column instance 0

ukb_data_hba1c <- ukb_data %>%
  select(contains("f.30750"))


# Get indices of participants with HbA1c >= 48 mmol/mol

# ukb_data_hba1c_indices_vector <- which(ukb_data_hba1c >= 0.065)

hba1c_indices <- apply(ukb_data_hba1c, 1, function(x) any(x >= 48, na.rm = T))

hba1c_eids <- eids[hba1c_indices]

#
##
### Age diabetes diagnosed workflow 
##
#

# Filter for field 2976 (age diabetes diagnosed) column instance 0

ukb_data_2976 <- ukb_data %>%
  select(contains("f.2976"))

# get indices for different age brackets

f2976_unknown <- t(apply(ukb_data_2976, 1, function(x) x < 0))

f2976_unknown_any <- apply(f2976_unknown, 1, any)

f2976_age <- rep(NA, length = nrow(f2976_unknown))
unknown_any <- f2976_unknown[!is.na(f2976_unknown_any),]
for(i in 1:nrow(unknown_any)) {
  unknown_i <- unknown_any[i,]
  if(all(unknown_i[!is.na(unknown_i)])) {
    f2976_age[!is.na(f2976_unknown_any)][i] <- NA
  }
}

f2976_known <- t(apply(ukb_data_2976, 1, function(x) x >= 0))
known_rs <- rowSums(f2976_known, na.rm = T)
for(i in 1:length(f2976_age)) {
  if(known_rs[i] >= 1) {
    f2976_age[i] <- min(ukb_data_2976[i,], na.rm = T)
  }
}

ukb_2976_age_lt_50_indices <- which(f2976_age < 50)
diabetes_age_lt_50_eids <- eids[ukb_2976_age_lt_50_indices]

ukb_2976_age_lteq_36_indices <- which(f2976_age <= 36)
diabetes_age_lteq_36_eids <- eids[ukb_2976_age_lteq_36_indices]

ukb_2976_age_lteq_30_indices <- which(f2976_age <= 30)
diabetes_age_lteq_30_eids <- eids[ukb_2976_age_lteq_30_indices]

#
##
### Ethnic background workflow 
##
#

eur<-c(1001,1002,1003)
mix<-c(2001,2002,2003)
south_asian<-c(3001,3002,3003)
african_caribbean<-c(4001,4002,4003)

# Filter for ethnic background column instance 0

ukb_data_ethnic <- ukb_data %>%
  select(contains("f.21000"))

# Get indices of participants who are European

ethnic_euro_indices <- which(apply(ukb_data_ethnic, 1, function(x) any(eur %in% x)))

euro_eids <- eids[ethnic_euro_indices]

# Get indices of participants who are mixed ethnicity

ethnic_mix_indices <- which(apply(ukb_data_ethnic, 1, function(x) any(mix %in% x)))

mixed_eids <- eids[ethnic_mix_indices]

# Get indices of participants who are South Asian

ethnic_sa_indices <- which(apply(ukb_data_ethnic, 1, function(x) any(south_asian %in% x)))

sa_eids <- eids[ethnic_sa_indices]

# Get indices of participants who are African-Caribbean

ethnic_ac_indices <- which(apply(ukb_data_ethnic, 1, function(x) any(african_caribbean %in% x)))

ac_eids <- eids[ethnic_ac_indices]



#
##
### Started insulin within one year diagnosis of diabetes workflow 
##
#

# Filter for field 2986 (started insulin within...) column instance 0

ukb_2986_insulin_lt1yr <- ukb_data %>%
  select(contains("f.2986"))

ukb_2986_insulin_lt1yr_indices <- which(apply(ukb_2986_insulin_lt1yr, 1, function(x) 1 %in% x))

insulin_lt1yr_eids <- eids[ukb_2986_insulin_lt1yr_indices]


#
##
### Gestational diabetes only workflow 
##
#

# Filter for field gestational diabetes only column instance 0

ukb_data_gest <- ukb_data %>%
  select(contains("f.4041"))

# Get the gestational diabetes only participant indices

gest_indices <- which(apply(ukb_data_gest, 1, function(x) 1 %in% x))

gest_eids <- eids[gest_indices]


#
##
### ICD10 diabetes workflow 
##
#


###### For ICD10 E11 T2DM

# Filter for field 41270 diagnoses ICD10 column instance 0

ukb_data_icd10 <- ukb_data %>%
  select(contains("f.41270"))

# Get the E11 (Non-insulin dependent diabetes) participant indices

icd10_t2dm_indices <- apply(as.matrix(ukb_data_icd10), 1, function(x) any(str_detect(x, "E11")))

icd10_t2dm_eids <- eids[which(icd10_t2dm_indices)]


###### For ICD10 E10 T1DM

icd10_t1dm_indices <- apply(as.matrix(ukb_data_icd10), 1, function(x) any(str_detect(x, "E10")))

icd10_t1dm_eids <- eids[which(icd10_t1dm_indices)]


###### For ICD10 E12 malnutrition dm

icd10_maln_dm_indices <- apply(as.matrix(ukb_data_icd10), 1, function(x) any(str_detect(x, "E12")))

icd10_maln_dm_eids <- eids[which(icd10_maln_dm_indices)]


###### For ICD10 E13 other specified dm

icd10_other_spec_dm_indices <- apply(as.matrix(ukb_data_icd10), 1, function(x) any(str_detect(x, "E13")))

icd10_other_spec_dm_eids <- eids[which(icd10_other_spec_dm_indices)]


###### For ICD10 E14 unspecified dm

icd10_unspec_dm_indices <- apply(as.matrix(ukb_data_icd10), 1, function(x) any(str_detect(x, "E14")))

icd10_unspec_dm_eids <- eids[which(icd10_unspec_dm_indices)]

# make eids data.frame

eids_list <- list(self_t2dm_eids,
                  self_t1dm_eids,
                  self_gest_eids,
                  self_diab_eids,
                  diabetes_diag_eids,
                  diabetes_medi_eids,
                  diabetes_atc_eids,
                  self_6153_insulin_eids,
                  self_6177_insulin_eids,
                  hba1c_eids,
                  primary_hba1c_eid,
                  t1dm_primary_eid,
                  t2dm_primary_eid,
                  diabetes_age_lt_50_eids,
                  diabetes_age_lteq_36_eids,
                  diabetes_age_lteq_30_eids,
                  euro_eids,
                  mixed_eids,
                  sa_eids,
                  ac_eids,
                  insulin_lt1yr_eids,
                  gest_eids,
                  icd10_t2dm_eids,
                  icd10_t1dm_eids,
                  icd10_other_spec_dm_eids,
                  icd10_unspec_dm_eids)

names(eids_list) <- c("self_t2dm", "self_t1dm", "self_gest", "self_diab", "diab_diag", "diab_medic", "diab_atc", "self_6153_insulin",
                      "self_6177_insulin", "hba1c", "primary_care_hba1c", "primary_care_t1dm", "primary_care_t2dm", "diab_age_lt_50", "diab_age_lteq_36", "diab_age_lteq_30",
                      "euro", "mixed", "south_asian", "afro_caribbean", "insulin_lt1yr", "gest_diab", "icd10_t2dm", "icd10_t1dm",
                      "icd10_other_spec_dm", "icd10_unspec_dm")

# Add eIDs

eids_df <- as.data.frame(matrix(nrow = length(eids), ncol = length(eids_list)+1))
colnames(eids_df) <- c("eid", names(eids_list))
eids_df$eid <- eids

for(i in 1:length(eids_list)) {
  crit_i <- eids_list[[i]]
  eid_in <- eids_df$eid %in% crit_i
  eids_df[, i+1] <- eid_in
}

# write out eIDs T2DM criteria table

write_tsv(eids_df, t2dm_crit_eids_df_file)

