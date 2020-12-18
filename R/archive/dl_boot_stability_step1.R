##################################################################
# Step 1 -- resample and bootstrap within each training          #
#        -- collect modeling results and feature importance rank # 
#        -- 10 resamples, 20 boots                               #
##################################################################

#set up
setwd("~/dkd")
rm(list=ls()); gc()
source("./helper_functions.R")
require_libraries("devtools")
slam_url<-"https://cran.r-project.org/src/contrib/Archive/slam/slam_0.1-37.tar.gz"
install_url(slam_url)
require_libraries(c( "Matrix"
                    ,"pROC"
                    ,"h2o"
                    ,"dplyr"
                    ,"tidyr"
                    ,"magrittr"
                    ,"data.table"
                    ,"slam"
))
options("h2o.use.data.table"=TRUE)

## Load in cleaned-up data sets
resamples<-10 #random sampling param
boots<-20 #boostrapping
load("DKD_heron_pats_prep.Rdata")
load("DKD_heron_facts_prep.Rdata")
load("feature_dict.Rdata")
load(paste0("random_sample",resamples,".Rdata"))
load(paste0("random_sample",resamples,"_boots",boots,".Rdata"))
load("feature_dict.Rdata") #variable dictionary

##global values
hyper_params<-list(
  activation=c("Rectifier"), #default
  hidden=list(rep(100,2),rep(200,2),rep(300,2),rep(400,2),rep(500,2),
              rep(100,3),rep(200,3),rep(300,3),rep(400,3),rep(500,3)),
  input_dropout_ratio=0.1,   #common choice: 0.1,0.2
  l1=1e-5,  #default
  l2=1e-5   #default
)

search_criteria<-list(
                      strategy = "RandomDiscrete",
                      max_runtime_secs = 360,
                      max_models = 100,
                      stopping_rounds=5,
                      stopping_tolerance=1e-2
                      )


#tune, train, predict and evaluate stability
#initialization
feature_lst<-list() #resample * nfold
pred_real<-list() #resample
time_perf<-c()

#initialize h2o
h2o.init(nthreads=-1)

for(i in 1:10){
  start_i<-Sys.time()
  cat("start resample:",i,"\n")
  time_perf_i<-c()
  
  hyper_param_lst<-c()
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
    pred_idx<-which(colnames(train_mt)!="DKD_IND")
    target_idx<-which(colnames(train_mt)=="DKD_IND")
    train_h2o[,target_idx]<-as.factor(train_h2o[,target_idx]) 
    
    train_h2o_splits<-h2o.splitFrame(data=train_h2o,ratios=0.8)
    train_h2o_tr<-train_h2o_splits[[1]]
    train_h2o_ts<-train_h2o_splits[[2]]
    
    time_perf_i<-c(time_perf_i,paste0(Sys.time()-start_k,units(Sys.time()-start_k)))
    cat("......finish converting to h2o frame in",time_perf_i[length(time_perf_i)],"\n")
    
    start_k<-Sys.time()
    cat("......tune for deep learning and train the model \n")
    
    dl_random_grid <- h2o.grid(
      algorithm="deeplearning",
      grid_id = paste0("dl_grid_random",j),
      training_frame=train_h2o_tr,   ## make it faster
      validation_frame=train_h2o_ts, ## make it faster
      # training_frame=train_h2o,       ## cv
      # nfolds=5,                       ## cv
      x=pred_idx,
      y=target_idx,
      distribution="bernoulli",
      standardize=T,
      epochs=100,                     ## make it fast
      stopping_metric="logloss",
      stopping_tolerance=1e-2,        ## stop when logloss does not improve by >=1% for 2 scoring events
      stopping_rounds=2,
      score_duty_cycle=0.025,         ## don't score more than 2.5% of the wall time
      max_w2=10,                      ## can help improve stability for Rectifier
      sparse=T,
      hyper_params = hyper_params,
      search_criteria = search_criteria,
      variable_importances=T          ## track variable importance
    )
    
    ## pick out the optimal model
    dlr_grid<-h2o.getGrid(paste0("dl_grid_random",j),sort_by="auc",decreasing=T)
    flr_opt_model<-h2o.getModel(dlr_grid@model_ids[[1]])
    flr_hyper_params<-flr_opt_model@parameters
    hyper_param_lst<-rbind(hyper_param_lst,
                           cbind(hidden=paste(flr_hyper_params$hidden,collapse = "_"),
                                 input_dropout_ratio=flr_hyper_params$input_dropout_ratio,
                                 l1=flr_hyper_params$l1,
                                 l2=flr_hyper_params$l2))
    rm(dlr_grid); gc()
    
    time_perf_i<-c(time_perf_i,paste0(Sys.time()-start_k,units(Sys.time()-start_k)))
    cat("......finish tuning and traning in",time_perf_i[length(time_perf_i)],"\n")
    
    start_k<-Sys.time()
    cat("......make prediction \n")
    
    test_h2o<-as.h2o(test_mt)
    pred_real_i[(pred_real_i$boot==j),"hidden"]<-paste(flr_hyper_params$hidden,collapse = "_")
    pred_real_i[((pred_real_i$boot==j)&(pred_real_i$part73=="InB")),"pred"]<-as.data.frame(predict(flr_opt_model,train_h2o))$p1
    pred_real_i[((pred_real_i$boot==j)&(pred_real_i$part73=="OOB")),"pred"]<-as.data.frame(predict(flr_opt_model,test_h2o))$p1

    time_perf_i<-c(time_perf_i,paste0(Sys.time()-start_k,units(Sys.time()-start_k)))
    cat("......finish predicting in",time_perf_i[length(time_perf_i)],"\n")
                   
                          
    #####feature selections
    start_k<-Sys.time()
    cat("......collect important features \n")

    decode_colnames<-data.frame(Feature=colnames(train_mt),
                                col_code=paste0("C",seq_along(colnames(train_mt))),
                                stringsAsFactors = FALSE)
    
    #important features 
    feat_ij<-h2o.varimp(flr_opt_model) %>%
      left_join(decode_colnames,by=c("variable"="col_code")) %>%
      dplyr::select(Feature,relative_importance,scaled_importance,percentage) %>%
      left_join(feat_dict,by=c("Feature"="CONCEPT_CD")) %>% unique #decode
    feature_lst[[paste0("resample",i,"boot",j)]]<-feat_ij
    
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
  
  #used hyperparameters from the previous bootstrapped rounds and train
  hyper_param_opt<-as.data.frame(hyper_param_lst)
  hyper_param_opt %<>%
    group_by_(.dots=names(hyper_param_opt)) %>%
    dplyr::summarise(sel_n=n()) %>%
    dplyr::arrange(desc(sel_n)) %>%
    ungroup %>% dplyr::slice(1)
  
  hidden_opt<-as.numeric(strsplit(as.character(hyper_param_opt$hidden),"_")[[1]])
  io_ratio<-as.numeric(as.character(hyper_param_opt$input_dropout_ratio))
  
  pred_idx<-which(colnames(train_mt)!="DKD_IND")
  target_idx<-which(colnames(train_mt)=="DKD_IND")
  train_h2o_tr[,target_idx]<-as.factor(train_h2o_tr[,target_idx]) 
  dl_random_single<-h2o.deeplearning(
    training_frame=train_h2o_tr,
    validation_frame=train_h2o_ts,
    x=pred_idx,
    y=target_idx,
    distribution="bernoulli",
    standardize=T,
    hidden=hidden_opt,              ## tuned 
    input_dropout_ratio=io_ratio,   ## tuned
    epochs=10,                      ## make it fast
    stopping_metric="logloss",
    stopping_tolerance=1e-2,        ## stop when logloss does not improve by >=1% for 2 scoring events
    stopping_rounds=2,
    score_duty_cycle=0.025,         ## don't score more than 2.5% of the wall time
    max_w2=10,                      ## can help improve stability for Rectifier
    sparse=T,
    variable_importances=T          ## track variable importance
    )

  #important features 
  feat_ij<-h2o.varimp(dl_random_single)%>%
    left_join(decode_colnames,by=c("variable"="col_code")) %>%
    dplyr::select(Feature,relative_importance,scaled_importance,percentage) %>%
    left_join(feat_dict,by=c("Feature"="CONCEPT_CD")) %>% unique #decode
  feature_lst[[paste0("resample",i,"single")]]<-feat_ij
  
  time_perf_i<-c(time_perf_i,paste0(Sys.time()-start_i,units(Sys.time()-start_i)))
  cat("...finish the single run for feature selection in",time_perf_i[length(time_perf_i)],"\n")
  
  #end outer loop
  time_perf_i<-c(time_perf_i,paste0(Sys.time()-start_i,units(Sys.time()-start_i)))
  cat("finish modeling resample",i,"in",time_perf_i[length(time_perf_i)],"\n")
  
  time_perf<-rbind(time_perf,time_perf_i)
  
  h2o.removeAll()
}

colnames(time_perf)<-c(paste0(c("partition",
                                "tune_and_train",
                                "predict",
                                "feature_collect",
                                "boots_all"),
                              rep(1:boots,each=5)),
                       c("single_run","resample_all"))
rownames(time_perf)<-paste0("resample",1:resamples)

#save results
save(feature_lst,file="dl_feature_boot.Rdata")
save(pred_real,file="dl_prediction_boot.Rdata")
save(time_perf,file="dl_performance_boot.Rdata")

#shut down h2o cluster
h2o.shutdown(prompt = FALSE)

