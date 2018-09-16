#### Preprocessing - part II ####
# - determine medication and diagnosis  
# - frequency filter
# - other fiters

rm(list=ls()); gc()
source("./R/util.R")
require_libraries(c("dplyr",
                    "tidyr",
                    "plyr",
                    "magrittr"
                    ))
pat_tbl<-readRDS("./data/DKD_heron_pats_prep.rda")
other_fact_stack<-readRDS("./data/DKD_heron_other_facts_prep.rda")
med_stack<-readRDS("./data/DKD_heron_med_facts_prep.rda")
dx_stack<-readRDS("./data/DKD_heron_dx_facts_prep.rda")

#RX at SCDF/SCBF(level4)
#DX at concept(level5)
fact_stack <- other_fact_stack %>%
  bind_rows(med_stack %>%
              dplyr::filter(CONCEPT_LEVEL==4) %>%
              mutate(CONCEPT_CD=paste0(CONCEPT_CD,"@",RX_MODIFIER)) %>%
              dplyr::select(PATIENT_NUM,
                            VARIABLE_CATEG,
                            C_VISUAL_PATH,
                            CONCEPT_CD,
                            C_NAME,
                            NVAL_NUM,
                            START_DATE)) %>%
  bind_rows(dx_stack %>%
              dplyr::filter(CONCEPT_LEVEL==5) %>%
              mutate(CONCEPT_CD=paste0(CONCEPT_CD,"@",DX_MODIFIER)) %>%
              dplyr::select(PATIENT_NUM,
                            VARIABLE_CATEG,
                            C_VISUAL_PATH,
                            CONCEPT_CD,
                            C_NAME,
                            NVAL_NUM,
                            START_DATE))

# View(fact_stack %>% group_by(VARIABLE_CATEG) %>%
#        dplyr::summarize(fact_cnt=n(),
#                         pat_cnt=length(unique(PATIENT_NUM)),
#                         cd_cnt=length(unique(CONCEPT_CD))))


## Preprocess 7 -- remove features with only non-DKD cases
overall_summary<-fact_stack %>%
  left_join((pat_tbl %>% dplyr::select(PATIENT_NUM,DKD_IND)),
            by="PATIENT_NUM") %>%
  group_by(VARIABLE_CATEG,C_VISUAL_PATH,CONCEPT_CD, C_NAME) %>%
  dplyr::summarize(pat_cnt = length(unique(PATIENT_NUM)),
                   dkd_cnt = length(unique(PATIENT_NUM * DKD_IND))-1) %>%
  mutate(nondkd_cnt = pat_cnt,
         dkd_rt = dkd_cnt/pat_cnt,
         nondkd_rt = 1-dkd_cnt/pat_cnt) %>%
  mutate(odd_ratio_emp=dkd_rt/nondkd_rt) %>%
  arrange(desc(odd_ratio_emp))

fcnt_before<-length(unique(fact_stack$CONCEPT_CD))  
cnt_before<-nrow(fact_stack)
fact_stack %<>% 
  semi_join((overall_summary %>% dplyr::filter(dkd_cnt>0)), by="CONCEPT_CD")
cnt_after<-nrow(fact_stack)
fcnt_after<-length(unique(fact_stack$CONCEPT_CD)) 
#report 
cat("Identified and removed",cnt_before - cnt_after,"facts and",
    fcnt_before - fcnt_after,"features for not containint DKD cases.\n")
# Identified and removed 82731 facts and 27928 features for not containint DKD cases.
gc()

# View(fact_stack %>% group_by(VARIABLE_CATEG) %>%
#        dplyr::summarize(fact_cnt=n(),
#                         pat_cnt=length(unique(PATIENT_NUM)),
#                         cd_cnt=length(unique(CONCEPT_CD))))

## Preprocess 8 -- frequency filtering
pass_ratio<-0.01
pass_filter<-fact_stack %>% 
  group_by(VARIABLE_CATEG) %>%
  do(mutate(.,pat_cnt_categ=length(unique(PATIENT_NUM)))) %>%
  group_by(VARIABLE_CATEG,C_VISUAL_PATH,CONCEPT_CD,C_NAME,pat_cnt_categ) %>%
  dplyr::summarize(pat_cnt = length(unique(PATIENT_NUM))) %>%
  dplyr::filter(pat_cnt >= pass_ratio*pat_cnt_categ)

fcnt_before<-length(unique(fact_stack$CONCEPT_CD))  
cnt_before<-nrow(fact_stack)
fact_stack %<>% 
  semi_join(pass_filter, by=c("CONCEPT_CD"))
cnt_after<-nrow(fact_stack)
fcnt_after<-length(unique(fact_stack$CONCEPT_CD)) 
#report 
cat("Filtered out",cnt_before - cnt_after,"facts and",
    fcnt_before - fcnt_after,"due to low frequency.\n")
# Filtered out 1233254 facts and 43487 due to low frequency.

View(fact_stack %>% group_by(VARIABLE_CATEG) %>%
       dplyr::summarize(fact_cnt=n(),
                        pat_cnt=length(unique(PATIENT_NUM)),
                        cd_cnt=length(unique(CONCEPT_CD))))

saveRDS(fact_stack, file="./data/DKD_heron_facts_prep.rda")
