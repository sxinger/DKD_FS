##################################################################
# Step 1 -- resample and bootstrap within each training          #
#        -- collect modeling results and feature importance rank # 
#        -- 10 resamples, 20 boots                               #
##################################################################

#set up
rm(list=ls()); gc()
source("./helper_functions.R")
require_libraries(c( "Matrix"
                     ,"pROC"
                     ,"xgboost"
                     ,"dplyr"
                     ,"tidyr"
                     ,"magrittr"
                     ,"pdp"
))


## Load in cleaned-up data sets
resamples<-10 #random sampling param
boots<-20 #boostrapping
load("DKD_heron_pats_prep.Rdata")
load("DKD_heron_facts_prep.Rdata")
load(paste0("random_sample",resamples,".Rdata"))
load(paste0("random_sample",resamples,"_boots",boots,".Rdata"))

##global values
#hyper-parameter grid for xgboost
eval_metric<-"auc"
objective<-"binary:logistic"
grid_params<-expand.grid(max.depth=c(2,4,6,8,10,12,14,16,18,20),
                         eta=c(0.5,0.3,0.1,0.05,0.02,0.01,0.005),
                         min_child_weight=1,
                         subsample=0.8,
                         colsample_bytree=0.8, 
                         gamma=1, 
                         nrounds=600)
#variable dictionary
pat_cnt_all<-nrow(pat_tbl)
feat_dict<-fact_stack %>%
  mutate(NVAL_NUM=ifelse(is.na(NVAL_NUM),0,NVAL_NUM)) %>%
  group_by(VARIABLE_CATEG,C_VISUAL_PATH,CONCEPT_CD, C_NAME) %>%
  dplyr::summarize(pat_cnt = length(unique(PATIENT_NUM)),
                   distinct_val=length(unique(NVAL_NUM)),
                   low=ifelse(length(unique(NVAL_NUM))==1,0,quantile(NVAL_NUM,probs=0.10)[1]),
                   mid=ifelse(length(unique(NVAL_NUM))==1,
                              round(length(unique(PATIENT_NUM))/pat_cnt_all,2),
                              median(NVAL_NUM, na.rm=T)),
                   high=ifelse(length(unique(NVAL_NUM))==1,1,quantile(NVAL_NUM,probs=0.90)[1])) %>% 
  bind_rows(pat_tbl %>%
              dplyr::select(-DKD_IND) %>%
              gather(CONCEPT_CD,NVAL_NUM,-PATIENT_NUM) %>%
              dplyr::filter(NVAL_NUM!=0) %>%
              group_by(CONCEPT_CD) %>%
              dplyr::summarize(pat_cnt = length(unique(PATIENT_NUM)),
                               distinct_val=length(unique(NVAL_NUM)),
                               low=ifelse(length(unique(NVAL_NUM))==1,0,quantile(NVAL_NUM,probs=0.10)[1]),
                               mid=ifelse(length(unique(NVAL_NUM))==1,
                                          round(length(unique(PATIENT_NUM))/pat_cnt_all,2),
                                          median(NVAL_NUM, na.rm=T)),
                               high=ifelse(length(unique(NVAL_NUM))==1,1,quantile(NVAL_NUM,probs=0.90)[1])) %>%
              mutate(VARIABLE_CATEG="DEMOGRAPHICS",
                     C_VISUAL_PATH="patient_dimension",
                     C_NAME=CONCEPT_CD) %>%
              dplyr::select(VARIABLE_CATEG,C_VISUAL_PATH,CONCEPT_CD, C_NAME,
                            pat_cnt,distinct_val,low,mid,high))

#eyeball random dample
# View(feat_dict[sample(seq_len(nrow(feat_dict)),50),])

#tune, train, predict and evaluate stability
#initialization
feature_lst<-list() #resample * nfold
pred_real<-list() #resample
time_perf<-c() #
metric_name<-paste0("test_", eval_metric,"_mean")

for(i in 1:10){
  start_i<-Sys.time()
  cat("start resample:",i,"\n")
  time_perf_i<-c()
  
  dat_sample<-dat_resample_rand_boot[[paste0("resample",i)]]
  pred_real_i<-dat_sample %>% mutate(pred=NA,real=NA)
  
  for(j in 1:boots){
    start_j<-Sys.time()
    cat("...start resample:",i,", boots:",j,"\n")
    
    dat_sample_j<-dat_sample %>% 
      dplyr::filter(boot==j) %>% 
      dplyr::select(PATIENT_NUM,PAT_IDX,part73)
    
    #####covert long df to wide sparse matrix (facts)
    x_sparse_val<-fact_stack %>% 
      inner_join(dat_sample_j,by=c("PATIENT_NUM")) %>% 
      dplyr::select(-PATIENT_NUM,-part73) %>%
      group_by(PAT_IDX) %>%
      long_to_sparse_matrix(.,id="PAT_IDX",variable="CONCEPT_CD",val="NVAL_NUM")
    
    x_sparse_pat<-pat_tbl %>%
      inner_join(dat_sample_j,by=c("PATIENT_NUM")) %>% 
      dplyr::select(-PATIENT_NUM,-part73,-year)  %>%
      arrange(PAT_IDX)
    
    y<-x_sparse_pat %>% dplyr::select(DKD_IND)
    x_sparse_pat %<>%
      dplyr::select(-DKD_IND) %>% 
      gather(key,value,-PAT_IDX) %>%
      long_to_sparse_matrix(.,id="PAT_IDX",variable="key",val="value")
    
    ######separate training and testing sets and convert to xgb.DMatrix
    trainx<-cbind(x_sparse_pat[(dat_sample_j$part73=="InB"),],
                  x_sparse_val[(dat_sample_j$part73=="InB"),])
    trainy<-as.vector(y[(dat_sample_j$part73=="InB"),])
    pred_real_i[((pred_real_i$boot==j)&(pred_real_i$part73=="InB")),"real"]<-trainy
    
    testx<-cbind(x_sparse_pat[(dat_sample_j$part73=="OOB"),],
                 x_sparse_val[(dat_sample_j$part73=="OOB"),])
    testy<-as.vector(y[(dat_sample_j$part73=="OOB"),])
    pred_real_i[((pred_real_i$boot==j)&(pred_real_i$part73=="OOB")),"real"]<-testy
    
    dtrain<-xgb.DMatrix(data=trainx,label=trainy)
    dtest<-xgb.DMatrix(data=testx,label=testy)
    
  
    ######train and predict
    start_k<-Sys.time()
    cat("......train and predict \n")
    
    watchlist<-list(train=dtrain, test=dtest)
    xgb_boot<-xgb.train(data=dtrain,
                        max_depth=grid_params$max.depth,
                        eta=grid_params$eta,
                        min_child_weight=grid_params$min_child_weight,
                        subsample=grid_params$subsample,
                        colsample_bytree=grid_params$colsample_bytree,
                        gamma=grid_params$gamma,
                        nrounds=grid_params$nrounds,
                        watchlist=watchlist,
                        maximize = TRUE,
                        eval_metric=eval_metric,
                        objective=objective,
                        print_every_n=100)
    
    pred_real_i[((pred_real_i$boot==j)&(pred_real_i$part73=="InB")),"pred"]<-as.numeric(predict(xgb_boot,dtrain))
    pred_real_i[((pred_real_i$boot==j)&(pred_real_i$part73=="OOB")),"pred"]<-as.numeric(predict(xgb_boot,dtest))
                                                           
    time_perf_i<-c(time_perf_i,paste0(Sys.time()-start_k,units(Sys.time()-start_k)))
    cat("......finish training and predicting in",time_perf_i[length(time_perf_i)],"\n")
                   
                          
    #####feature selections
    start_k<-Sys.time()
    cat("......collect important features \n")
    
    feat_ij<-xgb.importance(colnames(trainx),model=xgb_boot) %>%
      dplyr::filter(Gain > 0) %>%
      left_join(feat_dict,by=c("Feature"="CONCEPT_CD")) %>% unique #decode
    feature_lst[[paste0("resample",i,"boot",j)]]<-feat_ij
    
    time_perf_i<-c(time_perf_i,paste0(Sys.time()-start_k,units(Sys.time()-start_k)))
    cat("......finish collecting important features in",time_perf_i[length(time_perf_i)],"\n")
    
    #end inner loop
    time_perf_i<-c(time_perf_i,paste0(Sys.time()-start_i,units(Sys.time()-start_i)))
    cat("...finish modeling resample",i,",bootstrapped set",j,"in",time_perf_i[length(time_perf_i)],"\n")
  }
  pred_real[[paste0("resample",i)]]<-pred_real_i
  
  
  ####perform a single run on training sample
  start_k<-Sys.time()
  cat("...perform a single run just for feature selection \n")
  
  dat_sample_i<-dat_resample_rand[[paste0("resample",i)]] %>%
    dplyr::select(PATIENT_NUM, part73)
  
  x_sparse_val<-fact_stack %>% 
    inner_join(dat_sample_i,by=c("PATIENT_NUM")) %>% 
    dplyr::select(-part73) %>%
    group_by(PATIENT_NUM) %>%
    long_to_sparse_matrix(.,id="PATIENT_NUM",variable="CONCEPT_CD",val="NVAL_NUM")
  
  x_sparse_pat<-pat_tbl %>%
    inner_join(dat_sample_i,by=c("PATIENT_NUM")) %>% 
    dplyr::select(-part73)  %>%
    arrange(PATIENT_NUM)
  
  y<-x_sparse_pat %>% dplyr::select(DKD_IND)
  x_sparse_pat %<>%
    dplyr::select(-DKD_IND) %>% 
    gather(key,value,-PATIENT_NUM) %>%
    long_to_sparse_matrix(.,id="PATIENT_NUM",variable="key",val="value")
  
  dtrain<-xgb.DMatrix(data=cbind(x_sparse_pat[(dat_sample_i$part73=="T"),],
                                 x_sparse_val[(dat_sample_i$part73=="T"),]),
                      label=as.vector(y[(dat_sample_i$part73=="T"),]))
  
  xgb_single<-xgb.train(data=dtrain,
                        max_depth=grid_params$max.depth,
                        eta=grid_params$eta,
                        min_child_weight=grid_params$min_child_weight,
                        subsample=grid_params$subsample,
                        colsample_bytree=grid_params$colsample_bytree,
                        gamma=grid_params$gamma,
                        nrounds=grid_params$nrounds,
                        maximize = TRUE,
                        eval_metric=eval_metric,
                        objective=objective,
                        print_every_n=100)
  
  feat_ij<-xgb.importance(colnames(trainx),model=xgb_single) %>%
    dplyr::filter(Gain > 0) %>%
    left_join(feat_dict,by=c("Feature"="CONCEPT_CD")) %>% unique #decode
  feature_lst[[paste0("resample",i,"single")]]<-feat_ij
  
  time_perf_i<-c(time_perf_i,paste0(Sys.time()-start_i,units(Sys.time()-start_i)))
  cat("...finish the singel run for feature selection in",time_perf_i[length(time_perf_i)],"\n")
  
  #end outer loop
  time_perf_i<-c(time_perf_i,paste0(Sys.time()-start_i,units(Sys.time()-start_i)))
  cat("finish modeling resample",i,"in",time_perf_i[length(time_perf_i)],"\n")
  
  time_perf<-rbind(time_perf,time_perf_i)
}

colnames(time_perf)<-c(paste0(c("train_and_predict","feature_collect","boots_all"),rep(1:boots,each=3)),
                       c("single_run","resample_all"))
rownames(time_perf)<-paste0("resample",1:resamples)

save(feature_lst,file="xgb_feature_boot.Rdata")
save(pred_real,file="xgb_prediction_boot.Rdata")
save(time_perf,file="xgb_performance_boot.Rdata")

