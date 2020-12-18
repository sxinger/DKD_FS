##################################################################
# Step 1 -- resample and bootstrap within each training          #
#        -- collect modeling results and feature importance rank # 
#        -- 10 resamples, 20 boots                               #
##################################################################

#set up
setwd("~/dkd")
rm(list=ls()); gc()
source("./helper_functions.R")
require_libraries(c( "Matrix"
                     ,"pROC"
                     ,"h2o"
                     ,"dplyr"
                     ,"tidyr"
                     ,"magrittr"
))


## Load in cleaned-up data sets
resamples<-10 #random sampling param
boots<-20 #boostrapping
load("DKD_heron_pats_prep.Rdata")
load("DKD_heron_facts_prep.Rdata")
load(paste0("random_sample",resamples,".Rdata"))
load(paste0("random_sample",resamples,"_boots",boots,".Rdata"))

##global values
# hyper_params<-list(alpha=seq(1,0,by=-0.1))
hyper_params<-list(alpha=c(1,0.5))

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
feature_lst_elastnet<-list() #resample * nfold
feature_lst_lasso<-list() #resample * nfold
pred_real<-list() #resample
time_perf<-c() #

#initialize h2o
h2o.init(nthreads=-1)

for(i in 1:10){
  start_i<-Sys.time()
  cat("start resample:",i,"\n")
  time_perf_i<-c()
  
  dat_sample<-dat_resample_rand_boot[[paste0("resample",i)]]
  pred_real_i<-dat_sample %>% mutate(alpha_elastnet=NA,pred_elastnet=NA,pred_lasso=NA,real=NA)
  
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
      arrange(PAT_IDX) %>%
      gather(key,value,-PAT_IDX) %>%
      long_to_sparse_matrix(.,id="PAT_IDX",variable="key",val="value")

    ######separate training and testing sets and convert to h2o.frame
    train_mt<-cbind(x_sparse_pat[(dat_sample_j$part73=="InB"),],
                    x_sparse_val[(dat_sample_j$part73=="InB"),])
    
    test_mt<-cbind(x_sparse_pat[(dat_sample_j$part73=="OOB"),],
                   x_sparse_val[(dat_sample_j$part73=="OOB"),])
    
    pred_real_i[(pred_real_i$boot==j),"real"]<-x_sparse_pat[,"DKD_IND"]
    

    ######train and predict
    start_k<-Sys.time()
    cat("......inner partition of training set and convert to h2o frame \n")
    
    train_h2o<-as.h2o(train_mt)
    train_h2o_splits<-h2o.splitFrame(data=train_h2o,ratios=0.8)
    train_h2o_tr<-train_h2o_splits[[1]]
    train_h2o_ts<-train_h2o_splits[[2]]
    
    time_perf_i<-c(time_perf_i,paste0(Sys.time()-start_k,units(Sys.time()-start_k)))
    cat("......finish converting to h2o frame in",time_perf_i[length(time_perf_i)],"\n")
    
    start_k<-Sys.time()
    cat("......tune alpha for elastic net model (include lasso) and train the model \n")
    predictors<-colnames(train_h2o_tr)[(colnames(train_mt)!="DKD_IND")]
    response<-colnames(train_h2o_tr)[(colnames(train_mt)=="DKD_IND")]
    alpha_grid<-h2o.grid(x=predictors,
                         y=response,  
                         training_frame=train_h2o_tr,
                         validation_frame=train_h2o_ts,
                         algorithm = "glm",
                         grid_id = paste0("elastnet_alpha_grid",j),
                         family="binomial",
                         solver="COORDINATE_DESCENT",   #same optimization method as glmnet
                         lambda_search=TRUE,
                         early_stopping = TRUE,
                         missing_values_handling="Skip",
                         remove_collinear_columns=TRUE,
                         hyper_params = hyper_params,
                         search_criteria = list(strategy = "Cartesian")
                         )
    ## alpha search for elastic net---if time allowed
    # sortedGrid<-h2o.getGrid(paste0("elastnet_alpha_grid",j),sort_by = "auc",decreasing = TRUE)
    # alpha_opt_model1<-h2o.getModel(sortedGrid@model_ids[[1]])
    # alpha_opt<-alpha_opt_model1@parameters$alpha
    # if(alpha_opt < 1){
    #   alpha_opt_model<-alpha_opt_model1
    #   lasso_model<-h2o.getModel(alpha_grid@model_ids[[length(hyper_params$alpha)]]) #models are queued
    # }else{
    #   alpha_opt_model<-h2o.getModel(sortedGrid@model_ids[[2]])
    #   lasso_model<-alpha_opt_model1
    # }
    
    ## alpha=0.5 for elastic net
    alpha_opt<-0.5
    alpha_grid<-h2o.getGrid(paste0("elastnet_alpha_grid",j))
    alpha_opt_model<-h2o.getModel(alpha_grid@model_ids[[1]])
    lasso_model<-h2o.getModel(alpha_grid@model_ids[[2]])
    rm(alpha_grid); gc()
    
    time_perf_i<-c(time_perf_i,paste0(Sys.time()-start_k,units(Sys.time()-start_k)))
    cat("......finish tuning and traning in",time_perf_i[length(time_perf_i)],"\n")
    
    
    start_k<-Sys.time()
    cat("......predict using elastic net/lasso model \n")
    
    pred_real_i[(pred_real_i$boot==j),"alpha_elastnet"]<-alpha_opt
    
    test_h2o<-as.h2o(test_mt)
    pred_real_i[((pred_real_i$boot==j)&(pred_real_i$part73=="InB")),"pred_elastnet"]<-as.data.frame(predict(alpha_opt_model,train_h2o))$p1
    pred_real_i[((pred_real_i$boot==j)&(pred_real_i$part73=="OOB")),"pred_elastnet"]<-as.data.frame(predict(alpha_opt_model,test_h2o))$p1
    
    pred_real_i[((pred_real_i$boot==j)&(pred_real_i$part73=="InB")),"pred_lasso"]<-as.data.frame(predict(lasso_model,train_h2o))$p1    
    pred_real_i[((pred_real_i$boot==j)&(pred_real_i$part73=="OOB")),"pred_lasso"]<-as.data.frame(predict(lasso_model,test_h2o))$p1
    
    time_perf_i<-c(time_perf_i,paste0(Sys.time()-start_k,units(Sys.time()-start_k)))
    cat("......finish predicting in",time_perf_i[length(time_perf_i)],"\n")
                   
                          
    #####feature selections
    start_k<-Sys.time()
    cat("......collect important features \n")
    
    decode_colnames<-data.frame(Feature=colnames(train_mt),
                                col_code=paste0("C",seq_along(colnames(train_mt))),
                                stringsAsFactors = FALSE)
    
    #important features from lasso
    feat_lasso_ij<-h2o.varimp(lasso_model) %>%
      left_join(decode_colnames,by=c("names"="col_code")) %>%
      dplyr::filter(coefficients != 0) %>%
      dplyr::select(Feature,coefficients,sign) %>%
      left_join(feat_dict,by=c("Feature"="CONCEPT_CD")) %>% unique #decode
    feature_lst_lasso[[paste0("resample",i,"boot",j)]]<-feat_lasso_ij
    
    #important features from elastic net
    feat_elast_ij<-h2o.varimp(alpha_opt_model) %>%
      left_join(decode_colnames,by=c("names"="col_code")) %>%
      dplyr::filter(coefficients != 0) %>%
      dplyr::select(Feature,coefficients,sign) %>%
      left_join(feat_dict,by=c("Feature"="CONCEPT_CD")) %>% unique #decode
    feature_lst_elastnet[[paste0("resample",i,"boot",j)]]<-feat_elast_ij
    
    time_perf_i<-c(time_perf_i,paste0(Sys.time()-start_k,units(Sys.time()-start_k)))
    cat("......finish collecting important features in",time_perf_i[length(time_perf_i)],"\n")
    
    #end inner loop
    time_perf_i<-c(time_perf_i,paste0(Sys.time()-start_i,units(Sys.time()-start_i)))
    cat("...finish modeling resample",i,",bootstrapped set",j,"in",time_perf_i[length(time_perf_i)],"\n")
    
    h2o.removeAll()
  }
  pred_real[[paste0("resample",i)]]<-pred_real_i
  
  ####perform a single run on training sample
  start_k<-Sys.time()
  cat("...perform a single run just for feature selection \n")
  
  #convert and partition
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
    arrange(PATIENT_NUM) %>%
    gather(key,value,-PATIENT_NUM) %>%
    long_to_sparse_matrix(.,id="PATIENT_NUM",variable="key",val="value")
  
  train_mt<-cbind(x_sparse_pat[(dat_sample_i$part73=="T"),],
                  x_sparse_val[(dat_sample_i$part73=="T"),])
  

  train_h2o<-as.h2o(train_mt)
  train_h2o_splits<-h2o.splitFrame(data=train_h2o,ratios=0.8)
  train_h2o_tr<-train_h2o_splits[[1]]
  train_h2o_ts<-train_h2o_splits[[2]]
  gc()
  
  #tune and train
  alpha_grid_single<-h2o.grid(x=colnames(train_h2o_tr)[-2],
                              y=colnames(train_h2o_tr)[2],   #colnames are automatically converted to "C##", must use column index
                              training_frame=train_h2o_tr,
                              validation_frame=train_h2o_ts,
                              algorithm = "glm",
                              grid_id = "elastnet_alpha_grid",
                              family="binomial",
                              solver="COORDINATE_DESCENT",
                              lambda_search=TRUE,
                              early_stopping = TRUE,
                              missing_values_handling="Skip",
                              remove_collinear_columns=TRUE,
                              hyper_params = hyper_params,
                              search_criteria = list(strategy = "Cartesian")
                              )
  
  ## alpha search for elastic net---if time allowed
  # sortedGrid<-h2o.getGrid("elastnet_alpha_grid",sort_by = "auc",decreasing = FALSE)
  # alpha_opt_model1<-h2o.getModel(sortedGrid@model_ids[[1]])
  # alpha_opt<-alpha_opt_model1@alpha
  # if(alpha_opt < 1){
  #   alpha_opt_model<-alpha_opt_model1
  #   lasso_model<-h2o.getModel(alpha_grid@model_ids[[1]])
  # }else{
  #   alpha_opt_model<-h2o.getModel(sortedGrid@model_ids[[2]])
  #   lasso_model<-h2o.getModel(sortedGrid@model_ids[[1]])
  # }
  
  ## alpha=0.5 for elastic net
  alpha_grid_single<-h2o.getGrid("elastnet_alpha_grid")
  alpha_opt_model<-h2o.getModel(alpha_grid_single@model_ids[[1]])
  lasso_model<-h2o.getModel(alpha_grid_single@model_ids[[2]])
  rm(alpha_grid_single); gc()
  
  #important features from lasso
  feat_lasso_ij<-h2o.varimp(lasso_model) %>%
    left_join(decode_colnames,by=c("names"="col_code")) %>%
    dplyr::filter(coefficients != 0) %>%
    dplyr::select(Feature,coefficients,sign) %>%
    left_join(feat_dict,by=c("Feature"="CONCEPT_CD")) %>% unique #decode
  feature_lst_lasso[[paste0("resample",i,"single")]]<-feat_lasso_ij
  
  #important features from elastic net
  feat_elast_ij<-h2o.varimp(alpha_opt_model) %>%
    left_join(decode_colnames,by=c("names"="col_code")) %>%
    dplyr::filter(coefficients != 0) %>%
    dplyr::select(Feature,coefficients,sign) %>%
    left_join(feat_dict,by=c("Feature"="CONCEPT_CD")) %>% unique #decode
  feature_lst_elastnet[[paste0("resample",i,"single")]]<-feat_elast_ij
  
  time_perf_i<-c(time_perf_i,paste0(Sys.time()-start_i,units(Sys.time()-start_i)))
  cat("...finish the single run for feature selection in",time_perf_i[length(time_perf_i)],"\n")
  
  #end outer loop
  time_perf_i<-c(time_perf_i,paste0(Sys.time()-start_i,units(Sys.time()-start_i)))
  cat("finish modeling resample",i,"in",time_perf_i[length(time_perf_i)],"\n")
  
  time_perf<-rbind(time_perf,time_perf_i)
  
  h2o.removeAll()
}

colnames(time_perf)<-c(paste0(c("partition","tune_and_train","predict","feature_collect","boots_all"),rep(1:boots,each=5)),
                       c("single_run","resample_all"))
rownames(time_perf)<-paste0("resample",1:resamples)

save(feature_lst_lasso,file="lasso_feature_boot.Rdata")
save(feature_lst_elastnet,file="elastnet05_feature_boot.Rdata")
save(pred_real,file="elastnet_prediction_boot.Rdata")
save(time_perf,file="elastnet_performance_boot.Rdata")

h2o.shutdown(prompt = FALSE)
