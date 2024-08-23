# generate T2DM eids

print("Loading libraries...")

library(readr)
library(dplyr)
library(stringr)
library(purrr)

# Variables

t2dm_crit_eids_df_file = "t2dm_criteria_eids.tsv"

# Read in criteria eids

t2dm_crit_eids <- read_tsv(t2dm_crit_eids_df_file)

# inclusion criteria

t2dm_crit_incl <- which(t2dm_crit_eids$self_t2dm | t2dm_crit_eids$self_diab | t2dm_crit_eids$diab_diag |
                        t2dm_crit_eids$diab_medic | t2dm_crit_eids$diab_atc | t2dm_crit_eids$hba1c | t2dm_crit_eids$icd10_t2dm |
                        t2dm_crit_eids$primary_care_hba1c | t2dm_crit_eids$primary_care_t2dm)

t2dm_incl_eids <- t2dm_crit_eids$eid[t2dm_crit_incl]

# exclusion

t2dm_crit_excl <- which(t2dm_crit_eids$self_t1dm | t2dm_crit_eids$self_gest | t2dm_crit_eids$diab_age_lteq_36 |
                        t2dm_crit_eids$insulin_lt1yr | t2dm_crit_eids$gest_diab | t2dm_crit_eids$icd10_t1dm |
                        t2dm_crit_eids$primary_care_t1dm)

t2dm_excl_eids <- t2dm_crit_eids$eid[t2dm_crit_excl]

# get final eIDs

t2dm_eids_final <- t2dm_incl_eids[-which(t2dm_incl_eids %in% t2dm_excl_eids)]

t2dm_crit_eids$t2dm_final <- 0
t2dm_crit_eids$t2dm_final[which(t2dm_crit_eids$eid %in% t2dm_eids_final)] <- 1

# Upset plot of final eids criterias

t2dm_filt <- t2dm_crit_eids %>%
  filter(t2dm_final == 1)

t2dm_filt <- t2dm_filt %>%
  select(eid, self_t2dm, diab_medic, diab_atc, hba1c, icd10_t2dm, primary_care_hba1c, primary_care_t2dm)

crit_list <- list(self_t2dm = t2dm_filt$eid[t2dm_filt$self_t2dm],
                  diab_medic = t2dm_filt$eid[t2dm_filt$diab_medic],
                  diab_atc = t2dm_filt$eid[t2dm_filt$diab_atc],
                  hba1c = t2dm_filt$eid[t2dm_filt$hba1c],
                  icd10_t2dm = t2dm_filt$eid[t2dm_filt$icd10_t2dm],
                  primary_hba1c = t2dm_filt$eid[t2dm_filt$primary_care_hba1c],
                  primary_t2dm = t2dm_filt$eid[t2dm_filt$primary_care_t2dm])

names(crit_list) <- c("Self-reported T2DM",
                      "Self-reported T2DM medication",
                      "Self-reported T2DM medication (ATC mapping)",
                      "HbA1C > 6.5%",
                      "HES ICD10 T2DM",
                      "Primary care HbA1c > 6.5%",
                      "Primary care Read T2DM")

upset_in <- UpSetR::fromList(crit_list)

UpSetR::upset(upset_in, nsets = 7)

UpSetR::upset(upset_in, sets = c("HES ICD10 T2DM", "Primary care HbA1c > 6.5%", "Primary care Read T2DM"))


# write out final inds

write_tsv(as.data.frame(t2dm_crit_eids), t2dm_crit_eids_df_file)



