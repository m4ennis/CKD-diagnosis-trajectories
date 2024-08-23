# generate replication GWAS input

# generate GWAS inputs

# libraries

library(tidyverse)
library(data.table)
library(magrittr)
library(forcats)

# variables

work_dir = ""

intermediates_file = paste0(work_dir, "intermediates/replicate_t2dm_time_to_ckd_out.txt")

phecode_icd10_file = paste0(work_dir, "phecode_icd10.csv")

non_t2dm_output = paste0(work_dir, "intermediates/replicate_non_t2dm_time_to_ckd_out.txt")

# load intermediates

df_out <- fread(intermediates_file)

# read in non-T2DM outputs

non_t2dm <- fread(non_t2dm_output)

# filter for relevant variables

df_out <- df_out %>%
  select(eid, sex, age, block, event_dt) %>%
  filter(event_dt >= 0)

df_out %<>%
  mutate(cmb_age = round(age + event_dt/365))

non_t2dm_df <- non_t2dm %>%
  mutate(cmb_age = round(event_dt/365)) %>%
  select(eid, sex, age, block, event_dt, cmb_age)

# add T2DM variable

df_out %<>%
  mutate(t2dm = 1)

non_t2dm_df %<>%
  mutate(t2dm = 0)

# combine

df_out <- rbind(df_out,
                non_t2dm_df)

ckd_eids <- df_out %>%
  filter(block == "CKD") %>%
  pull(eid)

ckd_df <- df_out %>%
  filter(eid %in% ckd_eids)

# format comorbidities

df_out %<>%
  mutate(block = str_replace(block, "\\.", "_"))

cmbs <- unique(df_out$block)
cmbs <- cmbs[!cmbs %in% c("CKD", "Death")]

# plot case-control sizes

phecode_icd10 <- fread(phecode_icd10_file)

plot <- df_out %>%
  select(eid, sex, cmb_age, block) %>%
  pivot_wider(all_of(c("eid", "sex")), names_from = "block", values_from = "cmb_age") %>%
  mutate(across(-all_of(c("eid", "sex")), .fns = function(x) ifelse(is.na(x), 0, 1))) %>%
  pivot_longer(-all_of(c("eid", "sex", "CKD")), values_to = "pres", names_to = "cmb") %>%
  group_by(CKD, cmb) %>%
  summarise(freq = sum(pres == 1)) %>%
  arrange(cmb) %>%
  filter(cmb !="Death") %>%
  mutate(cmb = str_replace(cmb, "_", ".")) %>%
  mutate(cmb = phecode_icd10$Phenotype[match(cmb, phecode_icd10$PheCode)]) %>%
  mutate(CKD = factor(CKD, levels = c("1", "0"))) %>%
  mutate(cmb = factor(cmb, levels = c("Acute renal failure", "Heart failure NOS", "Atrial fibrillation and flutter",
                                      "Pleurisy; pleural effusion", "Myocardial infarction", "Other anemias",
                                      "Essential hypertension", "Other chronic ischemic heart disease, unspecified",
                                      "Gout", "Other specified peripheral vascular diseases", "Hypotension NOS",
                                      "Angina pectoris", "Pneumonia", "Iron deficiency anemias, unspecified or not due to blood loss",
                                      "Obesity", "Sleep disorders", "Gastritis and duodenitis", "Other disorders of intestine", "Cataract"))) %>%
  ggplot(aes(y = fct_rev(cmb), x = freq, fill = factor(CKD))) +
  geom_col(colour = "black") +
  scale_x_continuous(breaks = seq(0, 80000, 5000)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 7)) +
  scale_y_discrete(label = function(x) str_trunc(x, 40)) +
  scale_fill_viridis_d(option = "G", label = c("0" = "Control",
                                               "1" = "Case"), direction = -1) +
  xlab("Number of individuals") +
  ylab("High risk CKD disease background") +
  labs(fill = "Case-control")
ggsave("replication_comorbidity_ckd_case_control_size.png", plot = plot,
       width = 3000, height = 1200, units = "px")

case_control_df <- df_out %>%
  select(eid, sex, cmb_age, block) %>%
  pivot_wider(all_of(c("eid", "sex")), names_from = "block", values_from = "cmb_age") %>%
  mutate(across(-all_of(c("eid", "sex")), .fns = function(x) ifelse(is.na(x), 0, 1))) %>%
  pivot_longer(-all_of(c("eid", "sex", "CKD")), values_to = "pres", names_to = "cmb") %>%
  group_by(CKD, cmb) %>%
  summarise(freq = sum(pres == 1)) %>%
  arrange(cmb) %>%
  filter(cmb !="Death") %>%
  mutate(cmb = str_replace(cmb, "_", ".")) %>%
  mutate(cmb = phecode_icd10$Phenotype[match(cmb, phecode_icd10$PheCode)]) %>%
  mutate(CKD = factor(CKD, levels = c("1", "0"))) %>%
  mutate(cmb = factor(cmb, levels = c("Acute renal failure", "Heart failure NOS", "Atrial fibrillation and flutter",
                                      "Pleurisy; pleural effusion", "Myocardial infarction", "Other anemias",
                                      "Essential hypertension", "Other chronic ischemic heart disease, unspecified",
                                      "Gout", "Other specified peripheral vascular diseases", "Hypotension NOS",
                                      "Angina pectoris", "Pneumonia", "Iron deficiency anemias, unspecified or not due to blood loss",
                                      "Obesity", "Sleep disorders", "Gastritis and duodenitis", "Other disorders of intestine", "Cataract"))) %>%
  mutate(CKD = ifelse(CKD == 1, "cases", "controls")) %>%
  pivot_wider(cmb, names_from = CKD, values_from = freq) %>%
  mutate(cc_ratio = cases/controls)
write_tsv(case_control_df, paste0(work_dir, "intermediates/replication_case_control_sizes.txt"))


# for each comorbidity event

for(i in 1:length(cmbs)) {
  cmb_i <- cmbs[i]
  # eid_test <- df_out %>%
  #   filter(event_dt >= 0) %>%
  #   filter(block == cmb_i) %>%
  #   pull(eid)
  
  # create model df
  
  model_df <- df_out %>%
    select(eid, sex, cmb_age, block) %>%
    filter(block %in% c("CKD", cmb_i)) %>%
    pivot_wider(all_of(c("eid", "sex")), names_from = "block", values_from = "cmb_age") %>%
    filter(!is.na(get(cmb_i))) %>%
    rename(age = cmb_i) %>%
    mutate(CKD = ifelse(is.na(CKD), 0, 1),
           comorb_event = paste0("replicate_", cmb_i)) %>%
    mutate(background = "all")
  
  write_tsv(model_df, paste0("t2dm_time_to_ckd_replicate_", cmb_i, "_gwas_in.tsv"))
}

