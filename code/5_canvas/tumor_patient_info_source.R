setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('../source/source_code.R')


df_canvas <- rio::import('../5_canvas/20251031TPHP_canvas_plot.xlsx')
df_canvas %<>% mutate(FileName = str_replace(FileName, '\\.', '-'))

meta <- rio::import('../6_3paired_tumor_nontumor/output/T_NT_info_check.xlsx')

info5 <- rio::import('../0_process_DIA-NN_data/input/20251009_PUH_sample_information_3005files_V6.xlsx')
info.sample <- info5 %>%
  filter(sample_type != 'p') %>% 
  select(-cancer, -cancer_abbr) %>% 
  left_join(meta %>% distinct(tissue_name, cancer, cancer_abbr, cancer_subtype))

info.ga <- rio::import('../0_process_DIA-NN_data/input/patient_info/patient_age_gender.xlsx')



# table 1B ------
# patient information with tumor
tbl1b <- info.sample %>%
  filter(sample_type %in% c('T', 'NT')) %>% 
  select(FileName, Age, Gender, tissue_name, cancer,
         anatomical_classification, patient_ID) %>% 
  inner_join(df_canvas %>% select(FileName, DDA_lib_type))
# tbl1b %>% count(Age)
# tbl1b %>% count(Gender)
pid.missed <- tbl1b %>% filter(is.na(Age)) %>% pull(patient_ID)

info.ga.rescue <- info.ga %>% filter(`Admission ID` %in% pid.missed) %>% 
  rename(patient_ID = `Admission ID`)

tbl1b.rescue <- tbl1b %>%
  filter(is.na(Age)) %>% select(-Age, -Gender) %>%
  inner_join(info.ga.rescue)
tbl1b.remain <- tbl1b %>% anti_join(tbl1b.rescue, by = 'patient_ID')
tbl1b.new <- rbind(tbl1b.remain, tbl1b.rescue) %>% select(-FileName) %>% distinct()

c2_pid <- c('PUH_CA 511', 'PUH_CA 512', 'PUH_CA 513', 'PUH_CA 514', 'PUH_CA 515', 'PUH_CA 516', 'PUH_CA 517', 'PUH_CA 518', 'PUH_CA 519')
tbl1b.new %<>% mutate(Center = ifelse(patient_ID %in% c2_pid, 'C2-HS', 'C1-HM'))

tbl1b.new.rename <- tbl1b.new %>% filter(str_detect(patient_ID, '^\\d+$'))
tbl1b.new.named <- tbl1b.new %>% filter(!str_detect(patient_ID, '^\\d+$'))

tbl1b.new.named %<>% mutate(
  nameid = as.numeric(str_extract(patient_ID, ' (\\d+)$', group = 1)),
  PatientID = patient_ID
)
former.nameid <- max(tbl1b.new.named$nameid, na.rm = T) # 519

tbl1b.new.rename %<>% mutate(
  nameid = seq(former.nameid+1, length.out = nrow(.)),
  PatientID = str_c('PUH_CA ', nameid)
)

tbl1b.new.final <- rbind(tbl1b.new.named, tbl1b.new.rename) %>% 
  arrange(desc(is.na(nameid)), nameid, PatientID) %>% 
  select(PatientID, Age, Gender, tissue_name, cancer, anatomical_classification, DDA_lib_type, Center, everything())
.write_excel(tbl1b.new.final, 'tumor_patient_info_rename.xlsx')


t1b <- tbl1b.new.final %>% 
  filter(!is.na(nameid)) %>%
  select(PatientID:Center) %>%
  rename(TissueName = tissue_name, Cancer = cancer,
         AnatomicalClassification = anatomical_classification,
         DDALibraryType = DDA_lib_type)


# table 1C ------
tbl1c <- t1b %>% select(Cancer, Age, Gender)

# Compute overall summary
overall <- tbl1c %>%
  summarise(
    Cancer = "Overall",
    n = n(),
    male_count_pct = paste0(sum(Gender == "M"), " (", round(100 * mean(Gender == "M"), 1), ")"),
    age_mean_sd = paste0(round(mean(Age), 2), " (", round(sd(Age), 2), ")")
  )

# Compute per-cancer summary, sorted alphabetically by Cancer
by_cancer <- tbl1c %>%
  group_by(Cancer) %>%
  summarise(
    n = n(),
    male_count_pct = paste0(sum(Gender == "M"), " (", round(100 * mean(Gender == "M"), 1), ")"),
    age_mean_sd = paste0(round(mean(Age), 2), " (", round(sd(Age), 2), ")")
  ) %>%
  arrange(Cancer)

# Combine overall and per-cancer rows
summary_tbl <- bind_rows(overall, by_cancer)

# Compute p-value for age (ANOVA for mean differences across cancer types)
age_aov <- aov(Age ~ Cancer, data = tbl1c)
p_age <- summary(age_aov)[[1]][["Pr(>F)"]][1]
p_age_str <- ifelse(p_age < 0.001, "<0.001", format(round(p_age, 3), nsmall = 3))

# Compute p-value for gender (Chi-square test for proportion differences across cancer types)
gender_table <- table(tbl1c$Cancer, tbl1c$Gender)
chi_gender <- chisq.test(gender_table)
p_gender <- chi_gender$p.value
p_gender_str <- ifelse(p_gender < 0.001, "<0.001", format(round(p_gender, 3), nsmall = 3))

# Create p-value row
p_row <- data.frame(
  Cancer = "p",
  n = NA_integer_,
  male_count_pct = p_gender_str,
  age_mean_sd = p_age_str
)

# Combine into final data.frame
t1c <- bind_rows(summary_tbl, p_row) %>% 
  setNames(c("Cancer", "n", "Gender = male (%)", "Age (mean (SD))"))


# source data output -----
list(Cancer.patient.information = t1b,
     Cancer.patient.stat = t1c) %>% 
  .write_excel('source_data_table1BC.xlsx')


