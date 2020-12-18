#### Preprocessing - part I ####
# - pre-select features
# - impute nval_num based on whether the feature is numerical or binary
# - separate medication and diagnosis from other facts
# - collect different representations for medication and diagnosis

rm(list=ls()); gc()
source("./R/util.R")
require_libraries(c("Matrix",
                    "dplyr",
                    "tidyr",
                    "plyr",
                    "magrittr",
                    "stringr",
                    "moments"
))
## Load in fact_stack and pat_tbl
load("./data/DKD_heron_pats.Rdata")
load("./data/DKD_heron_facts.Rdata")

View(fact_stack %>% group_by(VARIABLE_CATEG) %>%
       dplyr::summarize(fact_cnt=n(),
                        pat_cnt=length(unique(PATIENT_NUM)),
                        cd_cnt=length(unique(CONCEPT_CD))))

## Preprocess 1 --- remove predictors collinearly related to outcome
# - creatinine
# - eGFR
# - ACR, PCR
filter_cd<-c(
   "COMPONENT_ID:2009#","COMPONENT_ID:3730#","COMPONENT_ID:7072#" #creatinine
  ,"COMPONENT_ID:7165#","LOINC:14959-1#" #ACR
  ,"COMPONENT_ID:7085#","LOINC:2889-4#" # PCR
  ,"COMPONENT_ID:191#","COMPONENT_ID:200#" #eGFR
  ,"DEM\\|AGEATV" #age at visit in VISIT_DETAIL
  ,"UHC\\|Gender","UHC\\|Race" #demo in UHC
  ,"NCDR\\|CURRENT_AGE","NCDR:207" #demo in NCDR
  ,"NTDS\\|Age","NTDS\\|AgeUnits","NTDS\\|Race","NTDS\\|Sex","NTDS\\|Ethnicity","NTDS\\|Height","NTDS\\|Weight" #demo in NTDS
)

cnt_before<-nrow(fact_stack)
fact_stack %<>% 
  dplyr::filter(!grepl(paste0("(",paste(filter_cd,collapse=")|("),")"),CONCEPT_CD))
cnt_after<-nrow(fact_stack)
#report 
cat("Removed",cnt_before - cnt_after,"facts as collearly relating to outcome and clearly redundant.\n")
# Removed 105930 facts as collearly relating to outcome.


##Preprocess 1B --- remove predictors from ABRIDGED, NAACCR and Demographic(except for Ethinicity and Place)
cnt_before<-nrow(fact_stack)
fact_stack %<>% 
  dplyr::filter(!(VARIABLE_CATEG %in% c("NAACCR"))) %>%
  dplyr::filter(!(VARIABLE_CATEG == "DEMOGRAPHICS" & !grepl("(ETHNICITY)|(GEO)",CONCEPT_CD)))
cnt_after<-nrow(fact_stack)
#report 
cat("Removed",cnt_before - cnt_after,"facts from ABRIDGED,NAACCR and DEMOGRAPHICS(except Ethnicity and Place).\n")
# Removed 715166 facts as collearly relating to outcome.


##Preprocess 2 --- Impute NULL
nval_summary<-fact_stack %>%
  group_by(VARIABLE_CATEG,CONCEPT_CD) %>%
  dplyr::summarize(pat_cnt = length(unique(PATIENT_NUM)),
                   distinct_val_narm=length(unique(NVAL_NUM))-ifelse(sum(is.na(NVAL_NUM))>0,1,0),
                   na_prop=round(sum(is.na(NVAL_NUM))/n(),4),
                   neg_prop=round(sum((NVAL_NUM<0))/n(),4),
                   zero_prop=round(sum((NVAL_NUM==0))/n(),4),
                   min_val=ifelse(sum(is.na(NVAL_NUM))/n()==1,NA,min(NVAL_NUM,na.rm=T)),
                   max_val=ifelse(sum(is.na(NVAL_NUM))/n()==1,NA,max(NVAL_NUM,na.rm=T))) %>%
  ungroup

nval_curation<-list(neg_val=nval_summary %>% filter(neg_prop>0),
                    num_val_w_na=nval_summary %>% filter(distinct_val_narm>0 & na_prop>0))
saveRDS(nval_curation,file="./data/nval_num_abnormality.rda")

nval_summary %<>% select(CONCEPT_CD,distinct_val_narm)
na_cnt_before<-nrow(fact_stack %>% filter(is.na(NVAL_NUM)))
fact_stack %<>%
  left_join(nval_summary,by="CONCEPT_CD") %>%
  mutate(NVAL_NUM=ifelse(distinct_val_narm==0 & is.na(NVAL_NUM),1,
                         ifelse(is.na(NVAL_NUM),0,NVAL_NUM)))
na_cnt_after<-nrow(fact_stack %>% filter(is.na(NVAL_NUM)))
#report 
cat("Impute",na_cnt_before - na_cnt_after,"values for nval_num. \n")
# Impute 6077522 values for nval_num. 

# View(fact_stack %>% group_by(VARIABLE_CATEG) %>%
#        dplyr::summarize(fact_cnt=n(),
#                         pat_cnt=length(unique(PATIENT_NUM)),
#                         cd_cnt=length(unique(CONCEPT_CD))))

## Preprocess 3 --- Correct possible typos
# - temperature conversion
fact_stack %<>%
  mutate(NVAL_NUM = ifelse(CONCEPT_CD=="KUH|PAT_ENC:TEMPERATURE"&NVAL_NUM<47,NVAL_NUM*1.8+32,NVAL_NUM)) 
gc()

## Preprocess 4 --- Identify and remove possible outliers
# Step 1 - get overall summary (/data dictionary/identify outliers)
overall_summary<-fact_stack %>% 
  group_by(VARIABLE_CATEG,C_VISUAL_PATH,CONCEPT_CD, C_NAME) %>%
  dplyr::summarize(pat_cnt = length(unique(PATIENT_NUM)),
                   distinct_val = length(unique(NVAL_NUM)),
                   overall_mean = mean(NVAL_NUM, na.rm=T),
                   overall_median = median(NVAL_NUM, na.rm=T),
                   overall_mode = get_mode(NVAL_NUM),
                   overall_min = min(NVAL_NUM, na.rm=T),
                   overall_perc5 = quantile(NVAL_NUM,probs=0.05, na.rm=T)[1],
                   overall_max = max(NVAL_NUM, na.rm=T),
                   overall_perc95 = quantile(NVAL_NUM,probs=0.95, na.rm=T)[1],
                   overall_sd = sd(NVAL_NUM, na.rm=T),
                   overall_skew = skewness(NVAL_NUM, na.rm=T),
                   overall_kurt = kurtosis(NVAL_NUM, na.rm=T)) %>%
  ungroup

# Step 2 - Identify outliers
outliers<-overall_summary %>%
  dplyr::select(VARIABLE_CATEG,C_VISUAL_PATH,CONCEPT_CD, 
                distinct_val, overall_perc5, overall_perc95, overall_skew, overall_kurt) %>% 
  dplyr::filter (overall_kurt > 10 & distinct_val > 10) %>%
  mutate(upper = overall_perc95+3*(overall_perc95-overall_perc5),
         lower = pmax(0,overall_perc5-3*(overall_perc95-overall_perc5))) %>%
  dplyr::select(VARIABLE_CATEG,C_VISUAL_PATH,CONCEPT_CD, lower, upper) %>%
  dplyr::arrange(VARIABLE_CATEG)

# Step 3 - Remove outliers
cnt_before<-nrow(fact_stack)
fact_stack %<>% 
  left_join(outliers, by=c("VARIABLE_CATEG","C_VISUAL_PATH","CONCEPT_CD")) %>%
  dplyr::filter(is.na(upper)|(NVAL_NUM<=upper & NVAL_NUM>=lower)) %>%
  dplyr::select(-lower,-upper)
cnt_after<-nrow(fact_stack)
#report 
cat("Identified and removed",cnt_before - cnt_after,"facts as outliers.\n")
# Identified and removed 1366 facts as outliers.
gc()


## Preprocess 5 -- take out medication and diagnoses facts
med_stack<-fact_stack %>%
  dplyr::filter(VARIABLE_CATEG == "MEDICATIONS")

dx_stack<-fact_stack %>%
  dplyr::filter(VARIABLE_CATEG == "DIAGNOSES")

other_fact_stack<-fact_stack %>%
  dplyr::filter(!VARIABLE_CATEG %in% c("MEDICATIONS","DIAGNOSES"))


## Preprocess 6 -- feature engineering
pred_pt<-30
#demographic - cold-coding
pat_tbl %<>%
  semi_join(fact_stack,by="PATIENT_NUM") %>%
  mutate(year = as.numeric(format(as.Date(END_DATE,"%Y-%m-%d %H:%M:%S"),"%Y"))) %>% #for identifying hold-out set
  mutate(AGE = round((as.numeric(difftime(END_DATE,BIRTH_DATE,units="days"))-pred_pt)/365.25)) %>%
  dplyr::select(PATIENT_NUM,DKD_IND,year,AGE,GENDER_MALE,RACE,RELIGION) %>%
  mutate(RACE=paste0("race_",gsub(" ","_",RACE)),
         RELIGION=paste0("religion_",gsub(" ","_",RELIGION))) %>%
  mutate(race_ref=1, relig_ref=1) %>%
  spread(RACE, race_ref, fill=0) %>%
  spread(RELIGION, relig_ref, fill=0)
# hist(pat_tbl$YR_SINCE_DM)


#medication - along the HERON-ontology tree
#-- raw value
med_stack %<>%
  mutate(C_VISUAL_PATH=gsub("\\\\i2b2\\\\","",C_VISUAL_PATH)) %>%
  mutate(NLEVEL=stringr::str_count(C_VISUAL_PATH,"\\\\"),
         VA_NLEVEL=stringr::str_count(C_VISUAL_PATH,"\\[")) %>%
  mutate(RX_MODIFIER=gsub("^.*[@]","",CONCEPT_CD),
         CONCEPT_LEVEL=6)

med_stack2<-med_stack %>%
  dplyr::select(PATIENT_NUM, VARIABLE_CATEG, C_VISUAL_PATH, 
                CONCEPT_CD, C_NAME, NVAL_NUM, START_DATE, 
                RX_MODIFIER, CONCEPT_LEVEL)

#-- concept level
med_stack %<>%
  mutate(CONCEPT_CD=gsub("[#@].*$","",CONCEPT_CD)) %>%
  mutate(NVAL_NUM=1, CONCEPT_LEVEL=5)

med_stack_new<-med_stack %>%
  dplyr::select(PATIENT_NUM, VARIABLE_CATEG, C_VISUAL_PATH, 
                CONCEPT_CD, C_NAME, NVAL_NUM, START_DATE, 
                RX_MODIFIER, CONCEPT_LEVEL)
med_stack2 %<>%
  bind_rows(med_stack_new)

#--SCDF or SCBF
med_stack_new<-med_stack %>%
  mutate(C_VISUAL_PATH = ifelse(NLEVEL==(VA_NLEVEL+1),
                                C_VISUAL_PATH,
                                ifelse(VA_NLEVEL==0,
                                       stringr::str_extract(C_VISUAL_PATH,
                                                            paste0("^([^\\\\]*\\\\){",pmin(3,NLEVEL),"}")),
                                       stringr::str_extract(C_VISUAL_PATH,
                                                            paste0("^([^\\\\]*\\\\){",(VA_NLEVEL+2),"}"))))) %>%
  mutate(C_NAME=gsub("[^\\\\]*\\\\","",gsub("\\\\$","",C_VISUAL_PATH))) %>%
  mutate(CONCEPT_CD=C_NAME, NVAL_NUM=1, CONCEPT_LEVEL=4) %>%
  dplyr::select(PATIENT_NUM, VARIABLE_CATEG, C_VISUAL_PATH, 
                CONCEPT_CD, C_NAME, NVAL_NUM, START_DATE, 
                RX_MODIFIER, CONCEPT_LEVEL)
med_stack2 %<>%
  bind_rows(med_stack_new)

# # #eyeball random sample
# View(med_stack_new[sample(seq_len(nrow(med_stack_new)),50),])

#--VA class levels
va_max<-max(med_stack$VA_NLEVEL)
for(lev in 1:va_max){
  med_stack_new<-med_stack %>%
    mutate(C_VISUAL_PATH=ifelse(VA_NLEVEL==0, 
                                stringr::str_extract(C_VISUAL_PATH,
                                                     paste0("^([^\\\\]*\\\\){2}")),
                                stringr::str_extract(C_VISUAL_PATH,
                                                     paste0("^([^\\\\]*\\\\){",
                                                            pmin((VA_NLEVEL+1),(lev+1)),
                                                            "}")))) %>%
    # mutate(RX_CLASS_UD=stringr::str_extract(C_VISUAL_PATH,
    #                                         paste0("^([^\\\\]*\\\\){",
    #                                                NLEVEL-pmin(lev,(VA_NLEVEL-1))+1,
    #                                                "}"))) %>%
    mutate(C_NAME=gsub("[^\\\\]*\\\\","",gsub("\\\\$","",C_VISUAL_PATH))) %>%
    mutate(CONCEPT_CD=coalesce(gsub("\\[|\\]","",stringr::str_extract(C_NAME,"\\[.*?\\]")),C_NAME)) %>%
    mutate(NVAL_NUM=1, CONCEPT_LEVEL=lev) %>%
    dplyr::select(PATIENT_NUM, VARIABLE_CATEG, C_VISUAL_PATH, 
                  CONCEPT_CD, C_NAME, NVAL_NUM, START_DATE, 
                  RX_MODIFIER, CONCEPT_LEVEL)
  
  med_stack2 %<>%
    bind_rows(med_stack_new)
}

#eyeball random sample
View(med_stack2[sample(seq_len(nrow(med_stack2)),100),])


#diagnoses - along the HERON-ontology tree
#-- at raw level
dx_stack_cp<-dx_stack

dx_stack<-dx_stack_cp
dx_stack %<>%
  mutate(C_VISUAL_PATH=str_replace(C_VISUAL_PATH,"\\\\i2b2\\\\","")) %>%
  dplyr::filter(!grepl("^(Diagnoses\\\\A)+", C_VISUAL_PATH)) %>% #remove
  mutate(C_VISUAL_PATH=paste0(str_replace(C_VISUAL_PATH,"Diagnoses\\\\",""),"\\")) %>%
  mutate(DX_type = str_replace(str_extract(C_VISUAL_PATH,"^([^\\\\]*\\\\){1}"),"\\\\","")) %>%
  mutate(NLEVEL=str_count(C_VISUAL_PATH,"\\\\")) %>%
  mutate(DX_MODIFIER=gsub("^.*[@]","",CONCEPT_CD),
         CONCEPT_LEVEL=7)

dx_stack2<-dx_stack %>%
  dplyr::select(PATIENT_NUM, VARIABLE_CATEG, C_VISUAL_PATH, 
                CONCEPT_CD, C_NAME, NVAL_NUM, START_DATE, 
                DX_MODIFIER, CONCEPT_LEVEL)

#-- at concept level
dx_stack %<>% 
  mutate(CONCEPT_CD=gsub("[@].*$","",CONCEPT_CD),
         CONCEPT_LEVEL=6) %>%
  mutate(DECIMAL=ifelse(grepl("^(Other Diagnoses Concepts)+",C_VISUAL_PATH), 1,
                        ifelse(grepl("DX_ID",CONCEPT_CD),
                               pmax(0,
                                    nchar(str_replace(paste0(str_replace(str_replace(C_VISUAL_PATH,
                                                                                     paste0("^([^\\\\]*\\\\){",(NLEVEL-2),"}"),""),
                                                                         "\\s.*",""),
                                                             "."),
                                          "[^\\.]*(\\.)",""))-1) + 1,
                               pmax(0,
                                    nchar(str_replace(paste0(CONCEPT_CD,"."),"[^\\.]*(\\.)","")) - 1))))

dx_stack_new<-dx_stack %>%
  dplyr::select(PATIENT_NUM, VARIABLE_CATEG, C_VISUAL_PATH, 
                CONCEPT_CD, C_NAME, NVAL_NUM, START_DATE, 
                DX_MODIFIER, CONCEPT_LEVEL)
dx_stack2 %<>%
  bind_rows(dx_stack_new)


#-- below DX interger level
# dx_stack %>% dplyr::filter(DECIMAL>3 & grepl("DX_ID",CONCEPT_CD)) 

for(lev in 4:5){
  dx_stack_new<-dx_stack %>%
    mutate(C_VISUAL_PATH=ifelse(grepl("DX_ID",CONCEPT_CD),
                                stringr::str_extract(C_VISUAL_PATH,paste0("^([^\\\\]*\\\\){",
                                                                          pmin((NLEVEL-1),((NLEVEL-DECIMAL)+(lev-3))),"}")),
                                stringr::str_extract(C_VISUAL_PATH,paste0("^([^\\\\]*\\\\){",
                                                                          pmin(NLEVEL,((NLEVEL-DECIMAL)+(lev-3))),"}")))) %>%
    # mutate(C_VISUAL_PATH = stringr::str_extract(C_VISUAL_PATH,paste0("^([^\\\\]*\\\\){",NLEVEL-pmin(lev,(NLEVEL-2))+1,"}"))) %>%
    mutate(C_NAME=str_replace_all(gsub("\\\\$","",C_VISUAL_PATH),"[^\\\\]*\\\\","")) %>%
    mutate(CONCEPT_CD=paste0(DX_type,":",
                             gsub("\\(|\\)","",str_replace(C_NAME,"\\s.*","")))) %>%
    mutate(CONCEPT_LEVEL=lev) %>%
    dplyr::select(PATIENT_NUM, VARIABLE_CATEG, C_VISUAL_PATH, 
                  CONCEPT_CD, C_NAME, NVAL_NUM, START_DATE, 
                  DX_MODIFIER, CONCEPT_LEVEL)
  
  dx_stack2 %<>%
    bind_rows(dx_stack_new)
}


#--at DX integer level
dx_stack_new<-dx_stack %>%
  mutate(C_VISUAL_PATH=stringr::str_extract(C_VISUAL_PATH,
                                            paste0("^([^\\\\]*\\\\){",
                                                   (NLEVEL-DECIMAL),"}"))) %>%
  # mutate(C_VISUAL_PATH = stringr::str_extract(C_VISUAL_PATH,paste0("^([^\\\\]*\\\\){",(NLEVEL-DECIMAL),"}"))) %>%
  mutate(C_NAME=str_replace_all(gsub("\\\\$","",C_VISUAL_PATH),"[^\\\\]*\\\\","")) %>%
  mutate(CONCEPT_CD=paste0(DX_type,":",
                           gsub("\\(|\\)","",str_replace(C_NAME,"\\s.*","")))) %>%
  mutate(CONCEPT_LEVEL=3) %>%
  dplyr::select(PATIENT_NUM, VARIABLE_CATEG, C_VISUAL_PATH, 
                CONCEPT_CD, C_NAME, NVAL_NUM, START_DATE, 
                DX_MODIFIER, CONCEPT_LEVEL)

dx_stack2 %<>%
  bind_rows(dx_stack_new)


#--above DX interger level
# dx_stack %>% dplyr::filter(NLEVEL-DECIMAL>4) %>% head

for(lev in 1:2){
  dx_stack_new<-dx_stack %>%
    mutate(C_VISUAL_PATH=str_extract(C_VISUAL_PATH,
                                     paste0("^([^\\\\]*\\\\){",pmin((NLEVEL-DECIMAL-1),(lev+1)),"}"))) %>%
    # mutate(C_VISUAL_PATH = str_extract(C_VISUAL_PATH,paste0("^([^\\\\]*\\\\){",NLEVEL-pmin(lev,(NLEVEL-2))+1,"}"))) %>%
    mutate(C_NAME=str_replace(gsub("\\\\$","",C_VISUAL_PATH),"[^\\\\]*\\\\","")) %>%
    mutate(CONCEPT_CD=paste0(DX_type,":",
                             gsub("\\(|\\)","",str_replace(C_NAME,"\\s.*","")))) %>%
    mutate(CONCEPT_LEVEL=lev) %>%
    dplyr::select(PATIENT_NUM, VARIABLE_CATEG, C_VISUAL_PATH, 
                  CONCEPT_CD, C_NAME, NVAL_NUM, START_DATE, 
                  DX_MODIFIER, CONCEPT_LEVEL)
  
  dx_stack2 %<>%
    bind_rows(dx_stack_new)
}

#eyeball random sample
# View(dx_stack2[sample(seq_len(nrow(dx_stack2)),100),])

saveRDS(pat_tbl, file="./data/DKD_heron_pats_prep.rda")
saveRDS(other_fact_stack, file="./data/DKD_heron_other_facts_prep.rda")
saveRDS(med_stack2, file="./data/DKD_heron_med_facts_prep.rda")
saveRDS(dx_stack2, file="./data/DKD_heron_dx_facts_prep.rda")
