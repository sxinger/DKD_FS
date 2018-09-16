#### Random Sampling - partitioning and boostrapping ####
rm(list=ls()); gc()
source("./R/util.R")
require_libraries(c("dplyr",
                    "tidyr",
                    "magrittr",
                    "plyr"
))


## Load in fact_stack and pat_tbl collected 
pat_tbl<-readRDS("./data/DKD_heron_pats_prep.rda")
fact_stack<-readRDS("./data/DKD_heron_facts_prep.rda")

resamples<-10 #random sampling param
boots<-20 #boostrapping
dat_resample_rand<-list()
dat_resample_rand_boot<-list()

pat_num<-pat_tbl %>%
  semi_join(fact_stack, by=c("PATIENT_NUM")) # drop new patients within 30d

pat_num2<-pat_num %>% dplyr::filter(year < 2017) # hold-out

for (r in 1:resamples){
  start_r<-Sys.time()
  
  #random sampling
  pat_num3<-pat_num2 %>% dplyr::select(PATIENT_NUM) %>%
    mutate(part73=base::sample(c("T","V"),nrow(.),c(0.7,0.3),replace=T)) %>%
    bind_rows(pat_num %>% dplyr::filter(year >= 2017) %>% 
                dplyr::select(PATIENT_NUM) %>% mutate(part73="H")) %>%
    arrange(PATIENT_NUM)
  
  dat_resample_rand[[paste0("resample",r)]]<-pat_num3
  
  pat_num3_T<-pat_num3 %>%
    dplyr::filter(part73=="T")
  
  i<-1
  resample_i_boot<-c()
  while(i<=boots || unique_pat < nrow(pat_num3_T)){ #ending criteria
    pat_boot<-pat_num3_T %>%
      left_join(data.frame(PATIENT_NUM=sample(pat_num3_T$PATIENT_NUM,nrow(pat_num3_T),replace=T), #sample with replacement
                           type=1),by = "PATIENT_NUM") %>%
      mutate(part73=ifelse(is.na(type),"OOB","InB")) %>% dplyr::select(-type) %>%
      arrange(PATIENT_NUM) %>% mutate(PAT_IDX=1:nrow(.)) %>% 
      mutate(boot=i)
    
    resample_i_boot %<>%
      bind_rows(pat_boot)
    
    i<-i+1
    unique_pat<-resample_i_boot %>%
      dplyr::filter(part73=="InB") %>% dplyr::select(PATIENT_NUM) %>% 
      unique %>% unlist %>% length
  }
  
  lapse_r<-Sys.time()-start_r
  cat('finish resample',r,'in',lapse_r,units(lapse_r),'.\n')
  
  dat_resample_rand_boot[[paste0("resample",r)]]<-resample_i_boot
}

save(dat_resample_rand,file=paste0("./data/random_sample",resamples,".Rdata"))
save(dat_resample_rand_boot,file=paste0("./data/random_sample",resamples,"_boots",boots,".Rdata"))
