# collect 440 features being selected by current method #
rm(list=ls()); gc()
source("./R/util.R")
require_libraries(c("Matrix",
                    "dplyr",
                    "tidyr",
                    "plyr",
                    "magrittr",
                    "httr",
                    "XML",
                    "jsonlite"
))


##=============== collect facts ======================
#load data
fact_stack<-readRDS("./data/DKD_heron_facts_prep.rda")
pat_tbl<-readRDS("./data/DKD_heron_pats_prep.rda")
opt_set<-readRDS("./result/feature_set_opt.rda")

fact_stack_sel<-fact_stack %>%
  dplyr::filter(!is.na(START_DATE)) %>%
  semi_join(opt_set %>% select(Feature) %>% unique,
            by=c("CONCEPT_CD"="Feature")) %>%
  group_by(PATIENT_NUM,CONCEPT_CD) %>%
  dplyr::filter(START_DATE==max(START_DATE)) %>%
  dplyr::select(-START_DATE) %>% unique %>% ungroup %>% 
  dplyr::mutate(VARIABLE_CATEG = ifelse(VARIABLE_CATEG %in% c("CARDIOLABTESTS",
                                                              "LABTESTS",
                                                              "MICROBIOLOGY"),
                                        "LABORATORY TESTS",VARIABLE_CATEG)) %>%
  bind_rows(pat_tbl %>%
              gather(CONCEPT_CD,NVAL_NUM,-PATIENT_NUM,-DKD_IND,-year) %>%
              dplyr::filter(NVAL_NUM!=0) %>%
              dplyr::mutate(VARIABLE_CATEG="DEMOGRAPHICS",
                            C_VISUAL_PATH="patient_dimension",
                            C_NAME=CONCEPT_CD) %>%
              select(PATIENT_NUM,VARIABLE_CATEG,C_VISUAL_PATH,
                     CONCEPT_CD,C_NAME,NVAL_NUM) %>%
              semi_join(opt_set %>% select(Feature) %>% unique,
                        by=c("CONCEPT_CD"="Feature"))) %>%
  arrange(PATIENT_NUM) %>%
  dplyr::mutate(PID = as.numeric(as.factor(PATIENT_NUM))) %>%
  left_join(pat_tbl %>% select(PATIENT_NUM,DKD_IND),
            by="PATIENT_NUM")

patient_cw<-fact_stack_sel %>% select(PATIENT_NUM, PID) %>% unique
#save patient_num,pid crosswalk
saveRDS(patient_cw,file="./data/selDKD_fs440_patient_cw.rda")

#discretization
fact_stack_sel %<>%
  left_join(readRDS("./data/feature_dict.rda") %>% 
              ungroup %>% select(CONCEPT_CD,distinct_val,q_5,q_15) %>% unique,
            by="CONCEPT_CD") %>%
  dplyr::mutate(NVAL_GRP=case_when(is.na(q_5) ~ NVAL_NUM,
                                   NVAL_NUM < q_5 ~ 1,
                                   NVAL_NUM >= q_5 & NVAL_NUM <= q_mid ~ 2,
                                   NVAL_NUM > q_mid & NVAL_NUM <= q_15 ~ 3,
                                   NVAL_NUM > q_15 ~ 4)) %>%
  dplyr::mutate(NVAL_GRP=ifelse(CONCEPT_CD!="AGE",NVAL_GRP,
                                case_when(NVAL_NUM >=18 & NVAL_NUM <=25 ~ 1,
                                          NVAL_NUM >=26 & NVAL_NUM <=35 ~ 2,
                                          NVAL_NUM >=36 & NVAL_NUM <=45 ~ 3,
                                          NVAL_NUM >=46 & NVAL_NUM <=55 ~ 4,
                                          NVAL_NUM >=56 & NVAL_NUM <=65 ~ 5,
                                          NVAL_NUM >65 ~ 6)))
#point check
#age
age<-fact_stack_sel %>% filter(CONCEPT_CD=="AGE")
table(age$NVAL_GRP)

#numerical var (e.g. SBP)
sbp<-fact_stack_sel %>% filter(CONCEPT_CD=="KUH|PAT_ENC:BP_SYSTOLIC")
table(sbp$NVAL_GRP)

fact_stack_out<-fact_stack_sel %>%
  select(PID,DKD_IND,VARIABLE_CATEG,C_VISUAL_PATH,
         CONCEPT_CD,C_NAME,NVAL_NUM) %>% unique

fact_dict<-fact_stack_sel %>%
  select(VARIABLE_CATEG,CONCEPT_CD,C_NAME,q_5,q_15) %>% unique %>%
  dplyr::rename(DESCRIPTION=C_NAME,Q1=q_5,Q3=q_15) %>%
  dplyr::mutate(DATA_TYPE=ifelse(!is.na(Q1),"NUM","CAT")) %>%
  select(VARIABLE_CATEG,CONCEPT_CD,DESCRIPTION,DATA_TYPE,Q1,Q3) %>%
  arrange(DATA_TYPE,VARIABLE_CATEG)

#pivot by original values
dkd_ind<-fact_stack_sel %>% select(PID,DKD_IND) %>% unique

dkd_pivot_nval<-dkd_ind %>%
  left_join(fact_stack_sel %>% filter(distinct_val==1) %>% 
              select(PID,CONCEPT_CD,NVAL_NUM) %>% unique %>%
              spread(CONCEPT_CD,NVAL_NUM,fill=0),by="PID") %>%
  left_join(fact_stack_sel %>% filter(distinct_val>1) %>%
              select(PID,CONCEPT_CD,NVAL_NUM) %>% unique %>%
              spread(CONCEPT_CD,NVAL_NUM),by="PID") %>%
  #attach calendar year of onsets
  left_join(readRDS("./data_sel/DKD_fs440_patient_cw.rda"),
            by="PID") %>%
  inner_join(readRDS("./data/DKD_heron_pats_prep.rda") %>%
               dplyr::select(PATIENT_NUM,year),
             by="PATIENT_NUM") %>%
  dplyr::select(-PATIENT_NUM)


#pivot by discritized values
dkd_pivot_nval_grp<-dkd_ind %>%
  left_join(fact_stack_sel %>% filter(distinct_val==1) %>% 
              select(PID,CONCEPT_CD,NVAL_GRP) %>% unique %>%
              spread(CONCEPT_CD,NVAL_GRP,fill=0),by="PID") %>%
  left_join(fact_stack_sel %>% filter(distinct_val>1) %>%
              select(PID,CONCEPT_CD,NVAL_GRP) %>% unique %>%
              spread(CONCEPT_CD,NVAL_GRP),by="PID") %>%
  #attach calendar year of onsets
  left_join(readRDS("./data_sel/DKD_fs440_patient_cw.rda"),
            by="PID") %>%
  inner_join(readRDS("./data/DKD_heron_pats_prep.rda") %>%
               dplyr::select(PATIENT_NUM,year),
             by="PATIENT_NUM") %>%
  dplyr::select(-PATIENT_NUM)

write.csv(fact_dict,file="./data_sel/DKD_fs440_dict.csv",row.names = F)
write.csv(fact_stack_out,file="./data_sel/DKD_fs440.csv",row.names = F)
write.csv(dkd_pivot_nval,file="./data_sel/DKD_fs440_orig.csv",row.names = F)
write.csv(dkd_pivot_nval_grp,file="./data_sel/DKD_fs440_discr.csv",row.names = F)

#=========re-pivot; not map to loinc=============
fact_stack_out<-read.csv("./data_sel/DKD_fs440.csv",stringsAsFactors = F)
fact_dict<-read.csv("./data_sel/DKD_fs440_dict.csv",stringsAsFactors = F)

dkd_ind<-fact_stack_out %>% select(PID,DKD_IND) %>% unique %>%
  #attach calendar year of onsets
  left_join(readRDS("./data_sel/DKD_fs440_patient_cw.rda"),
            by="PID") %>%
  inner_join(readRDS("./data/DKD_heron_pats_prep.rda") %>%
               dplyr::select(PATIENT_NUM,year),
             by="PATIENT_NUM") %>%
  dplyr::select(-PATIENT_NUM)

dkd_pivot_nval<-dkd_ind %>%
  left_join(fact_stack_out %>% 
              semi_join(fact_dict %>% filter(DATA_TYPE=="CAT"),by="CONCEPT_CD") %>% 
              select(PID,CONCEPT_CD,NVAL_NUM) %>% unique %>%
              spread(CONCEPT_CD,NVAL_NUM,fill=0),by="PID") %>%
  left_join(fact_stack_out %>% 
              semi_join(fact_dict %>% filter(DATA_TYPE=="NUM"),by="CONCEPT_CD") %>% 
              select(PID,CONCEPT_CD,NVAL_NUM) %>% unique %>%
              group_by(PID,CONCEPT_CD) %>%
              dplyr::summarise(NVAL_NUM=mean(NVAL_NUM,na.rm=T)) %>%
              ungroup %>%
              spread(CONCEPT_CD,NVAL_NUM),by="PID")


write.csv(dkd_pivot_nval,file="./data_sel/DKD_fs440_orig_v2.csv",row.names = F)
#=======================================================================================


##=============== map lab to loinc =====================
#collect component_to_loinc mapping
# config<-read.csv('./config.csv')
# require_libraries(c("DBI",
#                     "ROracle"))
# c_connect<-dbConnect(Oracle(),config$username,config$password,config$access)
# component_to_loinc<-dbGetQuery(c_connect,"select * from component_to_loinc")
# saveRDS(component_to_loinc,file="./data/component_to_loinc.rda")
component_to_loinc<-readRDS("./data/component_to_loinc.rda") %>%
  filter(!is.na(LOINCCODE))
fact_stack_out<-read.csv("./data_sel/DKD_fs440.csv",stringsAsFactors = F)
fact_dict<-read.csv("./data_sel/DKD_fs440_dict.csv",stringsAsFactors = F)

#get lab and vital concepts out
# lab_vital_cd<-fact_dict %>%
#   filter((VARIABLE_CATEG == "LABORATORY TESTS" | grepl("(KUH\\|PAT_ENC)+",CONCEPT_CD)) & DATA_TYPE=="NUM") %>%
#   left_join(readRDS("./data/feature_dict.rda") %>% ungroup %>%
#               dplyr::select(CONCEPT_CD,C_VISUAL_PATH),
#             by = "CONCEPT_CD") %>%
#   unique
  
fact_stack_out_cat<-fact_stack_out %>%
  semi_join(fact_dict %>% filter(DATA_TYPE=="CAT"),by="CONCEPT_CD")

fact_stack_num_manual<-fact_stack_out %>%
  filter(CONCEPT_CD == "AGE") %>%
  mutate(NVAL_GRP=case_when(NVAL_NUM>=18 & NVAL_NUM<=25 ~ 1,
                            NVAL_NUM>=26 & NVAL_NUM<=35 ~ 2,
                            NVAL_NUM>=36 & NVAL_NUM<=45 ~ 3,
                            NVAL_NUM>=46 & NVAL_NUM<=55 ~ 4,
                            NVAL_NUM>=56 & NVAL_NUM<=65 ~ 5,
                            NVAL_NUM>65 ~ 6),
         TVAL_CHAR=case_when(NVAL_NUM>=18 & NVAL_NUM<=25 ~ "18-25",
                             NVAL_NUM>=26 & NVAL_NUM<=35 ~ "26-35",
                             NVAL_NUM>=36 & NVAL_NUM<=45 ~ "36-45",
                             NVAL_NUM>=46 & NVAL_NUM<=55 ~ "46-55",
                             NVAL_NUM>=56 & NVAL_NUM<=65 ~ "56-65",
                             NVAL_NUM>65 ~ "66=<")) %>%
  bind_rows(fact_stack_out %>%
              filter(CONCEPT_CD == "KUMC|PACK_PER_DAY") %>%
              mutate(NVAL_GRP=case_when(NVAL_NUM==0 ~ 1,
                                        NVAL_NUM>0 & NVAL_NUM<=1 ~ 2,
                                        NVAL_NUM>1 & NVAL_NUM<=2 ~ 3,
                                        NVAL_NUM>2 & NVAL_NUM<=3 ~ 4,
                                        NVAL_NUM>=3 ~ 5),
                     TVAL_CHAR=case_when(NVAL_NUM==0 ~ "0 pack",
                                         NVAL_NUM>0 & NVAL_NUM<=1 ~ "less than 1 pack",
                                         NVAL_NUM>1 & NVAL_NUM<=2 ~ "1 to 2 packs",
                                         NVAL_NUM>2 & NVAL_NUM<=3 ~ "2 to 3 packs",
                                         NVAL_NUM>=3 ~ "at least 3 packs"))) %>%
  bind_rows(fact_stack_out %>%
              filter(CONCEPT_CD == "KUMC|TOBACCO_USED_YEARS") %>%
              mutate(NVAL_GRP=case_when(NVAL_NUM<5 ~ 1,
                                        NVAL_NUM>=5 & NVAL_NUM<=10 ~ 2,
                                        NVAL_NUM>10 & NVAL_NUM<=20 ~ 3,
                                        NVAL_NUM>20 & NVAL_NUM<=40 ~ 4,
                                        NVAL_NUM>40 ~ 5),
                     TVAL_CHAR=case_when(NVAL_NUM<5 ~ "less than 5 years",
                                        NVAL_NUM>=5 & NVAL_NUM<=10 ~ "5 to 10 years",
                                        NVAL_NUM>10 & NVAL_NUM<=20 ~ "10 to 20 years",
                                        NVAL_NUM>20 & NVAL_NUM<=40 ~ "20 to 40 years",
                                        NVAL_NUM>40 ~ "more than 40 years"))) %>%
  bind_rows(fact_stack_out %>%
              filter(CONCEPT_CD == "KUH|PAT_ENC:BMI") %>%
              mutate(NVAL_GRP=case_when(NVAL_NUM<18.5 ~ 1,
                                        NVAL_NUM>=18.5 & NVAL_NUM<25 ~ 2,
                                        NVAL_NUM>=25 & NVAL_NUM<30 ~ 3,
                                        NVAL_NUM>=30 ~ 4),
                     TVAL_CHAR=case_when(NVAL_NUM<18.5 ~ "underweight",
                                         NVAL_NUM>=18.5 & NVAL_NUM<25 ~ "normal",
                                         NVAL_NUM>=25 & NVAL_NUM<30 ~ "overweight",
                                         NVAL_NUM>=30 ~ "obese"))) %>%
  bind_rows(fact_stack_out %>%
              filter(CONCEPT_CD == "KUH|PAT_ENC:BP_DIASTOLIC") %>%
              mutate(NVAL_GRP=case_when(NVAL_NUM<80 ~ 1,
                                        NVAL_NUM>=80 & NVAL_NUM<90 ~ 2,
                                        NVAL_NUM>=90 & NVAL_NUM<100 ~ 3,
                                        NVAL_NUM>=100 ~ 4),
                     TVAL_CHAR=case_when(NVAL_NUM<80 ~ "normal",
                                        NVAL_NUM>=80 & NVAL_NUM<90 ~ "pre-hypertension",
                                        NVAL_NUM>=90 & NVAL_NUM<100 ~ "stage1-hypertension",
                                        NVAL_NUM>=100 ~ "stage2-hypertension"))) %>%
  bind_rows(fact_stack_out %>%
              filter(CONCEPT_CD == "KUH|PAT_ENC:BP_SYSTOLIC") %>%
              mutate(NVAL_GRP=case_when(NVAL_NUM<120 ~ 1,
                                        NVAL_NUM>=120 & NVAL_NUM<140 ~ 2,
                                        NVAL_NUM>=140 & NVAL_NUM<160 ~ 3,
                                        NVAL_NUM>=160 ~ 4),
                     TVAL_CHAR=case_when(NVAL_NUM<120 ~ "normal",
                                        NVAL_NUM>=120 & NVAL_NUM<140 ~ "pre-hypertension",
                                        NVAL_NUM>=140 & NVAL_NUM<160 ~ "stage1-hypertension",
                                        NVAL_NUM>=160 ~ "stage2-hypertension"))) %>%
  bind_rows(fact_stack_out %>%
              filter(CONCEPT_CD == "KUH|PAT_ENC:TEMPERATURE") %>%
              mutate(NVAL_GRP=case_when(NVAL_NUM<97.8 ~ 1,
                                        NVAL_NUM>=97.8 & NVAL_NUM<99.1 ~ 2,
                                        NVAL_NUM>=99.1 ~ 3),
                     TVAL_CHAR=case_when(NVAL_NUM<97.8 ~ "low",
                                        NVAL_NUM>=97.8 & NVAL_NUM<99.1 ~ "normal",
                                        NVAL_NUM>=99.1 ~ "high"))) %>%
  bind_rows(fact_stack_out %>%
              filter(CONCEPT_CD == "KUH|PAT_ENC:RESPIRATIONS") %>%
              mutate(NVAL_GRP=case_when(NVAL_NUM<12 ~ 1,
                                        NVAL_NUM>=12 & NVAL_NUM<=20 ~ 2,
                                        NVAL_NUM>20 ~ 3),
                     TVAL_CHAR=case_when(NVAL_NUM<12 ~ "low",
                                        NVAL_NUM>=12 & NVAL_NUM<=20 ~ "normal",
                                        NVAL_NUM>20 ~ "high"))) %>%
  mutate(NVAL_NUM = NVAL_GRP) %>% dplyr::select(-NVAL_GRP)

#attach LOINC reference range
lab_loinc_cd<-fact_dict %>%
  filter(DATA_TYPE=="NUM" & VARIABLE_CATEG=="LABORATORY TESTS") %>%
  dplyr::select(CONCEPT_CD) %>%
  filter(grepl("(COMPONENT_ID)+",CONCEPT_CD)) %>%
  dplyr::mutate(component_id = as.numeric(gsub(".*:","",gsub("#.*","",CONCEPT_CD))),
                suffix=gsub(".*#","",CONCEPT_CD)) %>%
  unique %>%
  inner_join(component_to_loinc,by=c("component_id"="COMPONENT_ID")) %>%
  dplyr::mutate(CONCEPT_CD2=paste0("LOINC:",LOINCCODE,"#",suffix)) %>%
  dplyr::select(CONCEPT_CD,CONCEPT_CD2,LOINCCODE) %>%
  filter(!is.na(CONCEPT_CD2)) %>%
  mutate(CONCEPT_CD2=paste0(CONCEPT_CD2,"_",
                            gsub("KUH\\|COMPONENT_ID\\:","",
                                 gsub("#.*","",CONCEPT_CD))))
# write.csv(lab_loinc_cd %>% dplyr::select(LOINCCODE),
#           file="./data_sel/loinc_sel.csv",
#           row.names=F)
lab_loinc_cd %<>%
  inner_join(read.csv("./data_sel/loinc_sel_ref_rg.csv",stringsAsFactors=F),
             by="LOINCCODE")

fact_stack_loinc<-fact_stack_out %>% 
  inner_join(lab_loinc_cd,by="CONCEPT_CD") %>%
  mutate(CONCEPT_CD=CONCEPT_CD2) %>% 
  mutate(NVAL_GRP=case_when(include_lb==1&NVAL_NUM < lb ~ 1,
                            include_lb==0&NVAL_NUM<=lb ~ 1,
                            include_lb==1&NVAL_NUM>=lb&include_ub==1&NVAL_NUM<=ub ~ 2,
                            include_lb==1&NVAL_NUM>=lb&include_ub==0&NVAL_NUM < ub ~ 2,
                            include_lb==0&NVAL_NUM > lb&include_ub==1&NVAL_NUM<=ub ~ 2,
                            include_lb==0&NVAL_NUM > lb&include_ub==1&NVAL_NUM < ub ~ 2,
                            include_ub==1&NVAL_NUM > ub ~ 3,
                            include_ub==0&NVAL_NUM>=ub ~ 3),
         TVAL_CHAR=case_when(include_lb==1&NVAL_NUM < lb ~ "low",
                             include_lb==0&NVAL_NUM<=lb ~ "low",
                             include_lb==1&NVAL_NUM>=lb&include_ub==1&NVAL_NUM<=ub ~ "normal",
                             include_lb==1&NVAL_NUM>=lb&include_ub==0&NVAL_NUM < ub ~ "normal",
                             include_lb==0&NVAL_NUM > lb&include_ub==1&NVAL_NUM<=ub ~ "normal",
                             include_lb==0&NVAL_NUM > lb&include_ub==1&NVAL_NUM < ub ~ "normal",
                             include_ub==1&NVAL_NUM > ub ~ "high",
                             include_ub==0&NVAL_NUM>=ub ~ "high")) %>%
  mutate(NVAL_NUM = NVAL_GRP) %>% dplyr::select(-NVAL_GRP) %>%
  dplyr::select(-CONCEPT_CD2,-LOINCCODE)

fact_stack_num_quant<-fact_stack_out %>%
  semi_join(fact_dict %>% filter(DATA_TYPE=="NUM"),by="CONCEPT_CD") %>%
  anti_join(fact_stack_num_manual,by="CONCEPT_CD") %>%
  anti_join(lab_loinc_cd,by="CONCEPT_CD") %>%
  left_join(readRDS("./data/feature_dict.rda") %>% 
              ungroup %>% select(CONCEPT_CD,distinct_val,q_5,q_mid,q_15) %>% unique,
            by="CONCEPT_CD") %>%
  dplyr::mutate(NVAL_GRP=case_when(NVAL_NUM < q_5 ~ 1,
                                   NVAL_NUM >= q_5 & NVAL_NUM < q_mid ~ 2,
                                   NVAL_NUM >= q_mid & NVAL_NUM < q_15 ~ 3,
                                   NVAL_NUM >= q_15 ~ 4),
                TVAL_CHAR=case_when(NVAL_NUM < q_5 ~ "low (below 25th)",
                                    NVAL_NUM >= q_5 & NVAL_NUM < q_mid ~ "mid low (25th-50th)",
                                    NVAL_NUM >= q_mid & NVAL_NUM < q_15 ~ "mid high (50th-75th)",
                                    NVAL_NUM >= q_15 ~ "high (above 75th)")) %>%
  mutate(NVAL_NUM = NVAL_GRP) %>% dplyr::select(-NVAL_GRP)

fact_stack_out2<-fact_stack_out_cat %>% 
  bind_rows(fact_stack_num_manual) %>%
  bind_rows(fact_stack_loinc) %>%
  bind_rows(fact_stack_num_quant)

fact_dict %<>% dplyr::select(-Q1,-Q3)
fact_dict2<-fact_dict %>% 
  semi_join(fact_stack_out_cat,by="CONCEPT_CD") %>%
  dplyr::mutate(VALUESET=1,VALUESET_DESCRIPTION="presence") %>%
  bind_rows(fact_dict %>%
              inner_join(fact_stack_num_manual %>% 
                          dplyr::select(CONCEPT_CD,NVAL_NUM,TVAL_CHAR) %>% unique,
                        by="CONCEPT_CD") %>%
              dplyr::rename(VALUESET=NVAL_NUM,VALUESET_DESCRIPTION=TVAL_CHAR)) %>%
  bind_rows(fact_dict %>%
              inner_join(lab_loinc_cd %>% dplyr::select(CONCEPT_CD,CONCEPT_CD2),
                         by="CONCEPT_CD") %>%
              inner_join(fact_stack_loinc %>% 
                           dplyr::select(CONCEPT_CD,NVAL_NUM,TVAL_CHAR) %>% unique,
                        by=c("CONCEPT_CD2"="CONCEPT_CD")) %>%
              dplyr::select(-CONCEPT_CD) %>%
              dplyr::rename(CONCEPT_CD=CONCEPT_CD2,VALUESET=NVAL_NUM,VALUESET_DESCRIPTION=TVAL_CHAR)) %>%
  bind_rows(fact_dict %>%
              inner_join(fact_stack_num_quant %>% 
                          dplyr::select(CONCEPT_CD,NVAL_NUM, TVAL_CHAR) %>% unique,
                        by="CONCEPT_CD") %>%
              dplyr::rename(VALUESET=NVAL_NUM,VALUESET_DESCRIPTION=TVAL_CHAR)) %>%
  unique

fact_dict2 %<>%
  bind_rows(fact_dict2 %>% dplyr::select(-VALUESET,-VALUESET_DESCRIPTION) %>% 
              unique %>% mutate(VALUESET=0,VALUESET_DESCRIPTION="absence")) %>%
  unique %>% arrange(VARIABLE_CATEG,CONCEPT_CD)

saveRDS(fact_stack_out2,file="./data_sel/DKD_fs440_v2.rda")
saveRDS(fact_dict2,file="./data_sel/DKD_fs440_dict_v2.rda")


#pivot by original values
dkd_ind<-fact_stack_out2 %>% select(PID,DKD_IND) %>% unique %>%
  #attach calendar year of onsets
  left_join(readRDS("./data_sel/DKD_fs440_patient_cw.rda"),
            by="PID") %>%
  inner_join(readRDS("./data/DKD_heron_pats_prep.rda") %>%
               dplyr::select(PATIENT_NUM,year),
             by="PATIENT_NUM") %>%
  dplyr::select(-PATIENT_NUM)

dkd_pivot_nval2<-dkd_ind %>%
  left_join(fact_stack_out2 %>% 
              select(PID,CONCEPT_CD,NVAL_NUM) %>% unique %>%
              group_by(PID,CONCEPT_CD) %>%
              dplyr::summarise(NVAL_NUM=mean(NVAL_NUM,na.rm=T)) %>%
              ungroup %>%
              spread(CONCEPT_CD,NVAL_NUM,fill=0),by="PID")

write.csv(dkd_pivot_nval2,file="./data_sel/DKD_fs440_discrt_v2.csv",row.names = F)
write.csv(fact_dict2,file="./data_sel/DKD_fs440_dict_v2.csv",row.names = F)


