#### Get Heron facts ####
## require: eligible patient_num (in line with current heron release),
##          eligible encounter_num (in line with current heron release)
rm(list=ls()); gc()
source("./R/util.R")
source("./R/sql_util.R")
require_libraries(c( "DBI"
                    ,"ROracle"
                    ,"dplyr"
                    ,"tidyr"
                    ,"magrittr"))

#set up connection with Oracle db
heronb2_config<-read.csv('../heronb2_config.csv')
c_connect<-dbConnect(Oracle(),heronb2_config$username,heronb2_config$password,heronb2_config$access)


dbGetQuery(c_connect,"select * from t")

# herona1_config<-read.csv('../herona1_config.csv')
# c_connect<-dbConnect(Oracle(),herona1_config$username,herona1_config$password,herona1_config$access)

#heron facts are organized in the following root folders
var_categ<-c(
  "ABRIDGED", #linked on patient_num?
  "ALERTS",
  "ALLERGY", #linked on patient_num?
  "CARDIOLABTESTS",
  # "CF",
  "DEMOGRAPHICS",
  "DIAGNOSES",
  # "FLOWSHEET",
  "HISTORY",
  "LABORATORY TESTS",
  "LABTESTS",
  # "MICROPOSITIVE",
  # "MICRONEGATIVE",
  "MICROBIOLOGY",
  "MEDICATIONS",
  "NAACCR", #linked on patient_num?
  "NCDR", #linked on patient_num?
  "NTDS", #linked on patient_num?
  "PROCEDURES",
  "PROCORDERS",
  "REPORTS",
  # "RESEARCHENROLLMENT",
  "SPECIMENS",
  "UHC",
  "VISIT DETAILS"
) #investigate a sub-list of variable categories

fact_stack<-c() #stack facts
time_bm<-c() #benchmark process time
size_bm<-c() #benchmark object size


#remotely create patient and encounter views/tables
send_enc_tbl(c_connect, 
             pat_tbl="ELIGIBLE_PATIENTS", #patient table
             enc_tbl="ENCOUNTERS_DKD", #encounter table
             pred_pt=30,
             pred_period=0)


#fetch and stack facts relative to encounters of interest
for(j in seq_along(var_categ)){
  cat('Start collect facts from HERON folder:',var_categ[j],'.\n')
  start_j<-Sys.time()
  
  if(var_categ[j] %in% c("ALERTS","DIAGNOSES",'MEDICATIONS',"PROCEDURES","PROCORDERS")){
    pat_fact<-get_recent_Heron_facts(c_connect,
                                     root_folder=var_categ[j],
                                     modifier=T,
                                     freq_bd=0,
                                     restr="ELIGIBLE_ENCOUNTERS") %>% 
      dplyr::select(PATIENT_NUM,VARIABLE_CATEG,C_VISUAL_PATH,
                    CONCEPT_CD,C_NAME,NVAL_NUM,START_DATE)
    
    fact_stack<-rbind(fact_stack,pat_fact) #stack facts
    
  }else if(var_categ[j] %in% c("CARDIOLABTESTS","LABORATORY TESTS","LABTESTS","MICROBIOLOGY","SPECIMENS")){
    pat_fact<-get_recent_LAB_facts(c_connect,
                                   root_folder=var_categ[j],
                                   freq_bd=0,
                                   restr="ELIGIBLE_ENCOUNTERS") %>% 
      dplyr::select(PATIENT_NUM,VARIABLE_CATEG,C_VISUAL_PATH,
                    CONCEPT_CD,C_NAME,NVAL_NUM,START_DATE)
    
    fact_stack<-rbind(fact_stack,pat_fact) #stack facts
    
  # }else if(var_categ[j]=='MEDICATIONS'){
  #   pat_fact<-get_recent_MED_facts(c_connect,
  #                                  ip_id_mod=c("MAR","Inpatient"), 
  #                                  ip_mod=c("MedObs:MAR_Dose"),
  #                                  op_mod=c("MedObs:Dispensed","Surescripts:Amount","MedObs:Dose","RX_REFILLS"),
  #                                  freq_bd=0,
  #                                  restr="ELIGIBLE_ENCOUNTERS") 
  #   
  #   fact_stack<-rbind(fact_stack,
  #                       (pat_fact %>% 
  #                          dplyr::select(PATIENT_NUM,VARIABLE_CATEG,C_VISUAL_PATH,
  #                                        CONCEPT_CD,C_NAME,NVAL_NUM,START_DATE))) #stack facts
  #   
  #   fact_stack<-rbind(fact_stack,
  #                       (pat_fact %>% 
  #                          dplyr::select(PATIENT_NUM,VARIABLE_CATEG,C_VISUAL_PATH,
  #                                        CONCEPT_CD,C_NAME,MED_CNT,START_DATE) %>%
  #                          mutate(CONCEPT_CD = paste0("freq_",CONCEPT_CD)) %>% 
  #                          dplyr::rename(NVAL_NUM = MED_CNT))) #stack facts
    
  }else if(var_categ[j] == "DEMOGRAPHICS"){
    pat_fact<-get_pat_demo_facts(c_connect,
                                 root_folder=var_categ[j],
                                 freq_bd=0,
                                 restr="ELIGIBLE_PATIENTS") %>% 
      dplyr::select(PATIENT_NUM,VARIABLE_CATEG,C_VISUAL_PATH,
                    CONCEPT_CD,C_NAME,NVAL_NUM,START_DATE)
    
    fact_stack<-rbind(fact_stack,pat_fact) #stack facts
    
  }else if(var_categ[j] %in% c("ABRIDGED","ALLERGY","NCDR","NTDS","NAACCR")){
    pat_fact<-get_recent_Heron_facts(c_connect,
                                     root_folder=var_categ[j],
                                     modifier=F,
                                     freq_bd=0,
                                     restr="ELIGIBLE_PATIENTS") %>% 
      dplyr::select(PATIENT_NUM,VARIABLE_CATEG,C_VISUAL_PATH,
                    CONCEPT_CD,C_NAME,NVAL_NUM,START_DATE)
    
    fact_stack<-rbind(fact_stack,pat_fact) #stack facts
    
  }else{
    pat_fact<-get_recent_Heron_facts(c_connect,
                                     root_folder=var_categ[j],
                                     modifier=F,
                                     freq_bd=0,
                                     restr="ELIGIBLE_ENCOUNTERS") %>% 
      dplyr::select(PATIENT_NUM,VARIABLE_CATEG,C_VISUAL_PATH,
                    CONCEPT_CD,C_NAME,NVAL_NUM,START_DATE)
    
    fact_stack<-rbind(fact_stack,pat_fact) #stack facts
  }
  
  #benchmarking
  lapse_j<-Sys.time()-start_j
  time_bm<-c(time_bm,paste0(round(lapse_j,2),units(lapse_j)))
  size_bm<-c(size_bm,nrow(pat_fact))
  
  #report progress
  cat('Finish stack facts of',var_categ[j],'in',lapse_j, units(lapse_j), '.\n')
  gc()
}

bm<-data.frame(time_bm=time_bm, size_bm=size_bm)
row.names(bm)<-var_categ


#save the tables and performance matrics
save(bm, file="DKD_heron_facts_bm.Rdata")
save(fact_stack,file="DKD_heron_facts.Rdata")
