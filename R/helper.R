### 01 read original data sets-----
lung_gene <- read.table("data/lung_gene.txt",
                        header = TRUE, sep = "\t"
) %>%
  data.frame()

lung_phe <- read.table("data/phenotype.txt",
                       header = TRUE, sep = "\t"
) %>%
  data.frame()

clinical_dat <- read_delim("data/lusc/clinical.tsv",
                           delim = "\t", escape_double = FALSE,
                           trim_ws = TRUE
) %>%
  bind_rows(
    read_delim("data/luad/clinical.tsv",
               delim = "\t", escape_double = FALSE,
               trim_ws = TRUE
    )
  ) %>%
  select(
    "case_id", "case_submitter_id",
    "project_id", "age_at_index", "age_at_diagnosis", "gender",
    "race", "year_of_birth",
    "vital_status",
    "ajcc_pathologic_m", "ajcc_pathologic_n", "ajcc_pathologic_stage",
    "ajcc_pathologic_t", "ajcc_staging_system_edition"
  ) %>%
  mutate(
    gender = ifelse(gender == "'--", NA, gender),
    ajcc_pathologic_stage = ifelse(ajcc_pathologic_stage == "'--", NA, ajcc_pathologic_stage),
    age_at_diagnosis = ifelse(age_at_diagnosis == "'--", NA, age_at_diagnosis)
  ) %>%
  unique()


####### 02 manipulate data

library(rlang) # set_names()
lung_gene_tmp <- lung_gene %>%
  t() %>%
  data.frame() %>%
  mutate(
    sample = rownames(.),
    sample = str_replace_all(sample, "\\.", "-")
  ) %>%
  select(sample, everything(.))

lung_gene <- lung_gene_tmp %>%
  set_names(lung_gene_tmp[1, ]) %>%
  filter(sample != "sample")

rm(lung_gene_tmp)


phe_gen_dat <- inner_join(lung_phe, lung_gene, by = "sample")



phe_gen_dat <- cbind(str_split(phe_gen_dat$sample, "-") %>% data.frame() %>%
                       t() %>% data.frame() %>% select(X4), phe_gen_dat) %>%
  filter(X4 == "01")

work_dat1 <- inner_join(clinical_dat, phe_gen_dat,
                        by = c("case_submitter_id" = "X_PATIENT")
) %>%
  filter(!is.na(OS.time)) %>%
  filter(!(OS == 1 & OS.time == 0)) %>%
  filter(!(OS == 0 & OS.time == 0)) %>%
  mutate(OS.time_yr = OS.time / 365)

work_dat1 <- work_dat1 %>%
  mutate(OS.time_month = OS.time/30,
         OS.time_yr = OS.time/365) %>%
  mutate(age_at_diagnosis = as.numeric(age_at_diagnosis)/356,
         cancer_stage = case_when(ajcc_pathologic_stage=='Stage IA'| ajcc_pathologic_stage=='Stage IB' ~'Stage I',
                                  ajcc_pathologic_stage=='Stage IIA'| ajcc_pathologic_stage=='Stage IIB' ~'Stage II',
                                  ajcc_pathologic_stage=='Stage IIIA'| ajcc_pathologic_stage=='Stage IIIB' ~'Stage III',
                                  TRUE ~ajcc_pathologic_stage
         ),
         cancer_type = ifelse(project_id =='TCGA-LUAD','LUAD','LUSC'),
         race = ifelse(race =='not reported',NA,race)
  ) %>%
  mutate(race_cat2 = ifelse(race=='white' |is.na(race),race,'other'),
         survival_status = ifelse(OS=='1',"Dead","Alive"),
         age_at_diagnosis1 = round(age_at_diagnosis,0),
         age_group  = case_when(age_at_diagnosis1< 65 ~ 1,
                                age_at_diagnosis1>=65 &age_at_diagnosis1 <75 ~ 2,
                                age_at_diagnosis1>=75 & age_at_diagnosis1<85 ~ 3,
                                age_at_diagnosis1>=85 ~ 4,
                                TRUE ~ age_at_diagnosis1
         ) %>% as.character()) %>%
  mutate(age_group1 = case_when(age_group==1~'Under 65',
                                age_group==2~'65-74',
                                age_group==3~'75-84',
                                age_group==4~'85 and over',
                                TRUE ~age_group)) %>%
  mutate(cancer_stage1 = case_when(cancer_stage=='Stage I'~'1',
                                   cancer_stage=='Stage II'~'2',
                                   cancer_stage=='Stage III'~'3',
                                   cancer_stage=='Stage IV'~'4',
  )) %>%
  mutate(age_group2 = ifelse((age_group==1 |age_group==2),1,2) %>% as.factor()) %>% mutate(cancer_stage_cat2 = case_when(cancer_stage=='Stage I'| cancer_stage=='Stage II' ~'Early stage',
                                                                                                                         cancer_stage=='Stage III'| cancer_stage=='Stage IV' ~'Late stage',
                                                                                                                         TRUE ~cancer_stage
  ))
uni_data <- work_dat1 %>%
  select(OS,OS.time,project_id,age_at_diagnosis,gender,race,ajcc_pathologic_stage)

gene_dat <- work_dat1[, 26:ncol(work_dat1)] %>%
  apply(.,2,as.numeric)
colnames(gene_dat) <- rlang::set_names(paste0('X_',
                                              make.names(colnames(gene_dat))))# Assuming your independent variables start from the third column


outcome_dat <- work_dat1 %>% select(OS.time_month,OS) %>%
  rlang::set_names("time","status")

cov_dat <- work_dat1 %>%
  select(age_at_diagnosis,cancer_stage_cat2) %>%
  rlang::set_names('age_at_diagnosis','cancer_stage')

cox_dat <- cbind(outcome_dat,cov_dat, gene_dat) %>%
  filter(!is.na(age_at_diagnosis)&!is.na(cancer_stage))
save(work_dat1,file="data/work_dat1.rdata")
save(uni_data,file="data/uni_data.rdata")
save(cox_dat,file='data/cox_dat.rdata')


####
