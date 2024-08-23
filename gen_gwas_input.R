# generate GWAS inputs

# libraries

library(tidyverse)
library(data.table)
library(magrittr)
library(ggrepel)

# variables

work_dir = ""

intermediates_file = paste0(work_dir, "intermediates/t2dm_time_to_ckd_envir.RData")

phecode_icd10_file = paste0(work_dir, "phecode_icd10.csv")

non_t2dm_output = paste0(work_dir, "intermediates/non_t2dm_to_ckd_output.rds")

# load intermediates

load(intermediates_file)

# read in non-T2DM outputs

non_t2dm <- readRDS(non_t2dm_output)

# filter for relevant variables

df_out <- df_out %>%
  select(eid, sex, age, block, event_dt) %>%
  filter(event_dt >= 0)

df_out %<>%
  mutate(cmb_age = round(age + event_dt/365))

non_t2dm_df <- non_t2dm$model_df %>%
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

# filter for blocks which are prevalent + significant in T2DM

df_out %<>%
  filter(block %in% c("CKD", "Death", fin_sig_hr$comorb))

# format comorbidities

df_out %<>%
  mutate(block = str_replace(block, "\\.", "_"))

# plot case-control sizes

phecode_icd10 <- fread(phecode_icd10_file)

t2dm_hr <- fin_sig_hr %>%
  mutate(bg = "t2dm")

non_t2dm_hr <- non_t2dm$fin_sig_hr %>%
  mutate(bg = "non_t2dm")

avg_hrs <- rbind(t2dm_hr, non_t2dm_hr %>% filter(comorb %in% t2dm_hr$comorb)) %>%
  mutate(comorb = phecode_icd10$Phenotype[match(comorb, phecode_icd10$PheCode)]) %>%
  group_by(comorb) %>%
  mutate(avg_hr = mean(exp.coef.))

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
  mutate(cmb = factor(cmb, levels = unique(avg_hrs$comorb[order(avg_hrs$avg_hr)]))) %>%
  ggplot(aes(y = cmb, x = freq, fill = factor(CKD))) +
  geom_col(colour = "black") +
  scale_x_continuous(breaks = seq(0, 60000, 4000)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 7)) +
  scale_y_discrete(label = function(x) str_trunc(x, 40)) +
  scale_fill_viridis_d(option = "G", label = c("0" = "Control",
                                               "1" = "Case"), direction = -1) +
  xlab("Number of individuals") +
  ylab("High risk CKD disease background") +
  labs(fill = "Case-control")
ggsave("comorbidity_ckd_case_control_size.png", plot = plot,
         width = 3000, height = 1200, units = "px")

# cases-controls size df

case_contr_df <- df_out %>%
  select(eid, sex, t2dm, cmb_age, block) %>%
  pivot_wider(all_of(c("eid", "sex", "t2dm")), names_from = "block", values_from = "cmb_age")

tot_size <- length(unique(case_contr_df$eid))
t2dm_size <- sum(case_contr_df$t2dm)
non_t2dm_size <- sum(case_contr_df$t2dm == 0)
t2dm_ckd_size <- sum(case_contr_df$t2dm == 1 & !is.na(case_contr_df$CKD))
non_t2dm_ckd_size <- sum(case_contr_df$t2dm == 0 & !is.na(case_contr_df$CKD))

case_contr_df <- data.frame(tot_size,
                            t2dm_size,
                            non_t2dm_size,
                            t2dm_ckd_size,
                            non_t2dm_ckd_size)

write_tsv(case_contr_df, "ckd_t2dm_cases_controls_sizes.txt")

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
write_tsv(case_control_df, paste0(work_dir, "/intermediates/discovery_case_control_sizes.txt"))


# plot hazard ratios in T2DM + non-T2DM

plot <- rbind(t2dm_hr, non_t2dm_hr %>% filter(comorb %in% t2dm_hr$comorb)) %>%
  mutate(comorb = phecode_icd10$Phenotype[match(comorb, phecode_icd10$PheCode)]) %>%
  group_by(comorb) %>%
  mutate(avg_hr = mean(exp.coef.)) %>%
  ungroup() %>%
  mutate(comorb = factor(comorb, levels = unique(comorb[order(avg_hr)]))) %>%
  ggplot(aes(x = avg_hr, y = comorb)) +
  geom_point() +
  theme_bw() +
  scale_y_discrete(label = function(x) str_trunc(x, 40)) +
  xlab("Average significant hazard ratio of  \ndisease background for CKD risk in  \nT2D and non-diabetic patients") +
  ylab("High risk CKD disease background") +
  scale_x_continuous(breaks = seq(1, 7, 1)) +
  coord_cartesian(xlim = c(1, 7))
ggsave("t2dm_non_t2dm_ckd_risk_sig_hr.png", plot = plot,
       width = 1600, height = 1200, units = "px")

plot <- rbind(t2dm_hr, non_t2dm_hr %>% filter(comorb %in% t2dm_hr$comorb)) %>%
  mutate(comorb = phecode_icd10$Phenotype[match(comorb, phecode_icd10$PheCode)]) %>%
  group_by(comorb) %>%
  mutate(avg_hr = mean(exp.coef.)) %>%
  mutate(bg = ifelse(bg == "t2dm", "T2D", "Non-diabetic")) %>%
  ungroup() %>%
  mutate(comorb = factor(comorb, levels = unique(comorb[order(avg_hr)]))) %>%
  ggplot(aes(x = exp.coef., y = comorb, col = factor(bg))) +
  geom_point() +
  theme_bw() +
  scale_y_discrete(label = function(x) str_trunc(x, 40)) +
  xlab("Significant hazard ratio of  \ndisease background for CKD risk in  \nT2D and non-diabetic patients") +
  ylab("High risk CKD disease background") +
  scale_x_continuous(breaks = seq(1, 9, 1)) +
  coord_cartesian(xlim = c(1, 9)) +
  labs(col = "Patient population")
ggsave("t2dm_non_t2dm_ckd_risk_ungrouped_sig_hr.png", plot = plot,
       width = 2400, height = 1200, units = "px")

# plot disease burden against CKD risk

risk <- rbind(t2dm_hr, non_t2dm_hr %>% filter(comorb %in% t2dm_hr$comorb)) %>%
  mutate(comorb = phecode_icd10$Phenotype[match(comorb, phecode_icd10$PheCode)])

burden <- df_out %>%
  select(eid, sex, t2dm, cmb_age, block) %>%
  pivot_wider(all_of(c("eid", "sex", "t2dm")), names_from = "block", values_from = "cmb_age") %>%
  mutate(across(-all_of(c("eid", "sex", "t2dm")), .fns = function(x) ifelse(is.na(x), 0, 1))) %>%
  pivot_longer(-all_of(c("eid", "sex", "t2dm", "CKD")), values_to = "pres", names_to = "cmb") %>%
  group_by(t2dm, cmb) %>%
  summarise(freq = sum(pres == 1))%>%
  mutate(cmb = str_replace(cmb, "_", ".")) %>%
  mutate(comorb = phecode_icd10$Phenotype[match(cmb, phecode_icd10$PheCode)]) %>%
  filter(comorb %in% risk$comorb)

burden %<>%
  mutate(freq = ifelse(t2dm == 0, freq/case_contr_df$non_t2dm_size, freq)) %>%
  mutate(freq = ifelse(t2dm == 1, freq/case_contr_df$t2dm_size, freq))

burden$t2dm <- ifelse(burden$t2dm, "non_t2dm", "t2dm")
burden$bg <- burden$t2dm
burden %<>%
  ungroup() %>%
  select(-t2dm)

risk %<>%
  full_join(burden, by = c("bg", "comorb"))

cols <- c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(8, "Set1"), RColorBrewer::brewer.pal(3, "Set2"))

lvls <- rbind(t2dm_hr, non_t2dm_hr %>% filter(comorb %in% t2dm_hr$comorb)) %>%
  mutate(comorb = phecode_icd10$Phenotype[match(comorb, phecode_icd10$PheCode)]) %>%
  group_by(comorb) %>%
  mutate(avg_hr = mean(exp.coef.)) %>%
  mutate(bg = ifelse(bg == "t2dm", "T2D", "Non-diabetic")) %>%
  ungroup() %>%
  mutate(comorb = factor(comorb, levels = unique(comorb[order(avg_hr)]))) %>%
  pull(comorb) %>% levels()


plot <- risk %>%
  mutate(contrib = exp.coef.*freq) %>%
  mutate(comorb = factor(comorb, levels = lvls)) %>%
  mutate(bg = ifelse(bg == "t2dm", "T2D", "Non-diabetic")) %>%
  ggplot(aes(y = comorb, x = contrib, col = factor(bg), label = comorb, group = comorb)) +
  geom_point() +
  theme_bw() +
  # scale_x_continuous(labels = scales::percent) +
  ylab("High risk CKD disease background") +
  xlab("CKD risk burden") +
  labs(color = "Patient population") +
  scale_y_discrete(label = function(x) str_trunc(x, 40))
ggsave("t2dm_non_t2dm_ckd_risk_burden.png", plot = plot,
       width = 2400, height = 1200, units = "px")


# for each comorbidity event

for(i in 1:length(fin_sig_hr$comorb)) {
  cmb_i <- str_replace(fin_sig_hr$comorb[i], "\\.", "_")
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
           comorb_event = cmb_i) %>%
    mutate(background = "all")
  
  write_tsv(model_df, paste0("t2dm_time_to_ckd_", cmb_i, "_gwas_in.tsv"))
}
