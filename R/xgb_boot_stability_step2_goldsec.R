##################################################################
# Step 2 -- estimate percentage of featrue selected              #
#        -- retrain with different subset of features            # 
#        -- test on both hold-in and hold-out set                #
#        -- record out-of-sample AUC and Kuncheva index          #
##################################################################

##set up
rm(list=ls()); gc()
source("./helper_functions.R")
require_libraries(c( "Matrix"
                    ,"pROC"
                    ,"xgboost"
                    ,"dplyr"
                    ,"tidyr"
                    ,"magrittr"
))

##load data
resamples<-10 #random sampling param
boots<-20 #boostrapping
load("DKD_heron_pats_prep.Rdata")
load("DKD_heron_facts_prep.Rdata")

load(paste0("random_sample",resamples,".Rdata"))
load("xgb_feature_boot.Rdata")
load("xgb_prediction_boot.Rdata")
load("feature_dict.Rdata")

##global values
#overall target outcome
overall_y_sort<-pat_tbl %>% semi_join(fact_stack,by="PATIENT_NUM") %>% 
  arrange(PATIENT_NUM) %>% dplyr::select(DKD_IND) %>% unlist

#hyper-parameter grid for xgboost
eval_metric<-"auc"
objective<-"binary:logistic"
nrounds<-1000

#can be set to range of values
grid_params<-expand.grid(max.depth=c(6,8,10),    
                         eta=c(0.05,0.02,0.01),
                         min_child_weight=1,
                         subsample=0.8,
                         colsample_bytree=0.8, 
                         gamma=1)

#optimal feature size search stopping criteria
feat_sel_k_low<-2
feat_sel_k_up<-600
inc_tol_p<-0.01
s<-c(50,100) # number of top rank of interests
phi<-(1+sqrt(5))/2

# ensemble feature, tune, train, predict and evaluate stability
feature_rk<-list()
fs_summary<-list() 
time_perf<-list() 

metric_name<-paste0("test_", eval_metric,"_mean")
metric_sd_name<-paste0("test_", eval_metric,"_std")

#start experiment
for(i in seq_len(resamples)){
  start_i<-Sys.time()
  cat("start resample:",i,"\n")
  time_perf_i_nm<-c() #track task
  time_perf_i<-c()    #track time
  
  #########################################################
  #####ensemble and select features
  start_k<-Sys.time()
  cat("...ensemble features \n")
  
  feature_i<-c()
  pred_real_i<-pred_real[[paste0("resample",i)]]
  oob_auc<-pred_real_i %>% ungroup %>%
    dplyr::filter(part73 == "OOB") %>%
    dplyr::select(boot,pred,real) %>%
    group_by(boot) %>%
    dplyr::summarize(oob_weight=pROC::auc(real,pred))
  
  for(b in 1:boots){
    feature_i %<>% 
      bind_rows(feature_lst[[paste0("resample",i,"boot",b)]] %>%
                  dplyr::select(Feature,Gain)%>%
                  mutate(rank=rank(-Gain),boot=b))
  }
  
  feature_i  %<>%
    left_join(oob_auc,by="boot") %>%
    mutate(wt_rank=rank*oob_weight,
           top_50_f=(rank<=s[1])*1,
           top_100_f=(rank<=s[2])*1,
           wt_top_50_f=(rank<=s[1])*oob_weight,
           wt_top_100_f=(rank<=s[2])*oob_weight,
           top_50_exp_f=exp(-rank/s[1]),
           top_100_exp_f=exp(-rank/s[2]),
           wt_top_50_exp_f=(exp(-rank/s[1]))*oob_weight,
           wt_top_100_exp_f=(exp(-rank/s[2]))*oob_weight) %>%
    group_by(Feature) %>%
    dplyr::summarize(sel_cnt=length(unique(boot)),
                     best_rank=min(rank,na.rm=T),
                     worst_rank=max(rank,na.rm=T),
                     mean_rank=mean(rank,na.rm=T),   
                     wt_mean_rank=sum(wt_rank)/sum(oob_weight),   
                     top_50=-mean(top_50_f),
                     wt_top_50=-sum(wt_top_50_f)/sum(oob_weight),
                     top_100=-mean(top_100_f),
                     wt_top_100=-sum(wt_top_100_f)/sum(oob_weight),
                     top_50_exp=-mean(top_50_exp_f),
                     wt_top_50_exp=-sum(wt_top_50_exp_f)/sum(oob_weight),
                     top_100_exp=-mean(top_100_exp_f),
                     wt_top_100_exp=-sum(wt_top_100_exp_f)/sum(oob_weight)) %>%
    full_join(feature_lst[[paste0("resample",i,"single")]] %>%
                dplyr::select(Feature,Gain) %>%
                mutate(single_rank=rank(-Gain)) %>%
                dplyr::select(Feature,single_rank), by="Feature")
  
  feature_rk[[paste0("resample",i)]]<-feature_i
  
  time_perf_i_nm<-c(time_perf_i_nm,"ensemble_feature")
  time_perf_i<-c(time_perf_i,paste0(Sys.time()-start_k,units(Sys.time()-start_k)))
  cat("...finish ensembling features in",time_perf_i[length(time_perf_i)],"\n")
  
  ################################################################
  ####experiment on different selection ratio
  fs_mth<-colnames(feature_i)[-c(1:4,7:8,11:12)]
  for(fs in seq_along(fs_mth)){
    start_fs<-Sys.time()
    cat("...rank feature importance based on",fs_mth[fs],"\n")
    
    #####train a classifier with selected features
    dat_sample_i<-dat_resample_rand[[paste0("resample",i)]] %>%
      arrange(PATIENT_NUM) #sort by patient_num
    
    #####adaptive feature inclusion
    feat_num<-list()
    
    #initialization
    a<-feat_sel_k_low
    d<-feat_sel_k_up
    b<-a+0.1
    c<-d-0.1
    auc_inc<--Inf
    inc_p<-1
    ROC_obj_new<-pROC::roc(overall_y_sort[(dat_sample_i$part73!="T")],
                           sample(c(0,1),nrow(dat_sample_i[(dat_sample_i$part73!="T"),]),replace=T),
                           direction="<") ## random chance
    ROC_obj_opt<-ROC_obj_new
    opt_size<-1
    global_min<-d
    ROC_opt_update<-0
    track_path<-c()   
    
    while(a < (d-1) && b < (c-1)){
      start_g<-Sys.time()
      cat("...start golden-section search \n")
      
      #track the approximation path
      track_path<-rbind(track_path,cbind(a=a,d=d))
      
      #update golden-section points: b,c
      b<-floor(d-(d-a)/phi)
      c<-floor(a+(d-a)/phi)
      bounds<-c(a,b,c,d)
      
      for(k in seq_along(bounds)){
        feat_sel_k<-bounds[k]
        
        if(!is.null(feat_num[[paste0("size_",feat_sel_k)]])){
          cat("...the case of keeping",feat_sel_k,"features has already been saved. \n")
        }else{
          start_j<-Sys.time()
          cat("...keep",feat_sel_k,"features \n")
          
          #####################################################################################
          start_k<-Sys.time()
          cat("......subset features and transform to wide sparse matrix \n")
          
          fs_sel<-feature_i %>%
            mutate(rk=get(fs_mth[fs])) %>% 
            arrange(rk) %>% dplyr::select(Feature,rk) %>% 
            dplyr::slice(1:feat_sel_k)
          
          x_sparse_pat<-pat_tbl %>%
            semi_join(dat_sample_i,by="PATIENT_NUM") %>%
            dplyr::select(-year) %>%
            arrange(PATIENT_NUM) #sort by patient_num
          y<-x_sparse_pat %>% dplyr::select(DKD_IND) #sort by patient_num
          
          x_sparse<-fact_stack %>% 
            semi_join(fs_sel, by=c("CONCEPT_CD"="Feature")) %>%    #subset features
            inner_join(dat_sample_i,by="PATIENT_NUM") %>% 
            dplyr::select(-part73) %>%
            bind_rows(x_sparse_pat %>%
                        dplyr::select(-DKD_IND) %>% 
                        gather(CONCEPT_CD,NVAL_NUM,-PATIENT_NUM) %>%
                        semi_join(fs_sel, by=c("CONCEPT_CD"="Feature"))) %>%
            group_by(PATIENT_NUM) %>%
            long_to_sparse_matrix(.,id="PATIENT_NUM",
                                  variable="CONCEPT_CD",
                                  val="NVAL_NUM") #sort by patient_num
          
          if(nrow(x_sparse)<dim(dat_sample_i)[1]){
            dat_sample_i %<>%
              semi_join(data.frame(PATIENT_NUM=as.numeric(rownames(x_sparse))),
                        by="PATIENT_NUM") #shrink feature space will reduce training data size as well
          }
          
          #record real y values
          dat_sample_i[,"real"]<-y
          
          time_perf_i_nm<-c(time_perf_i_nm,paste0("subset_feature_transform@",fs,"@",feat_sel_k))
          time_perf_i<-c(time_perf_i,paste0(Sys.time()-start_k,units(Sys.time()-start_k)))
          cat("......finish subsetting and transforming in",time_perf_i[length(time_perf_i)],"\n")
          
          #######################################################################################
          start_k<-Sys.time()
          cat("......separate training and testing sets \n")
          
          trainx<-x_sparse[(dat_sample_i$part73=="T"),]
          # colnames(trainx)<-c(colnames(x_sparse_pat),colnames(x_sparse_val)) #colname may be dropped when only one column is selected
          trainy<-as.vector(y[(dat_sample_i$part73=="T"),])
          
          testx<-x_sparse[(dat_sample_i$part73!="T"),]
          # colnames(testx)<-c(colnames(x_sparse_pat),colnames(x_sparse_val)) #colname may be dropped when only one column is selected
          testy<-as.vector(y[(dat_sample_i$part73!="T"),])
          
          dtrain<-xgb.DMatrix(data=trainx,label=trainy)
          dtest<-xgb.DMatrix(data=testx,label=testy)
          
          time_perf_i_nm<-c(time_perf_i_nm,paste0("partition@",fs,"@",feat_sel_k))
          time_perf_i<-c(time_perf_i,paste0(Sys.time()-start_k,units(Sys.time()-start_k)))
          cat("......finish partitioning in",time_perf_i[length(time_perf_i)],"\n") 
          
          ######################################################################################
          start_k<-Sys.time()
          cat("......tune, train and predict \n")
          
          hyper_perf<-c()
          roc_lst<-c()
          for(params in seq_len(dim(grid_params)[1])){
            bst_tune<-xgb.cv(
              grid_params[params,],
              dtrain, 
              nfold = 3, 
              objective = objective,
              metrics = eval_metric,
              maximize = TRUE, 
              nrounds=nrounds,
              early_stopping_rounds = 50,
              print_every_n = 200,
              prediction=TRUE)
            
            hyper_perf<-rbind(hyper_perf,
                              cbind(grid_params[params,],
                                    metric_avg=max(bst_tune$evaluation_log[[metric_name]]),
                                    metric_sd=max(bst_tune$evaluation_log[[metric_sd_name]]),
                                    steps=bst_tune$best_iteration))
            
            roc_lst[[params]]<-pROC::roc(trainy,bst_tune$pred)
          }
          
          hyper_opt<-hyper_perf[which.max(hyper_perf$metric_avg),] #take the best hyperparameter set
          watchlist<-list(train=dtrain, test=dtest)
          xgb_tune<-xgb.train(
            data=dtrain,
            max_depth=hyper_opt$max.depth,
            eta=hyper_opt$eta,
            min_child_weight=hyper_opt$min_child_weight,
            ubsample=hyper_opt$subsample,
            colsample_bytree=hyper_opt$colsample_bytree,
            gamma=hyper_opt$gamma,
            nrounds=hyper_opt$steps,
            watchlist=watchlist,
            eval_metric=eval_metric,
            objective=objective,
            print_every_n = 200)
          
          dat_sample_i[,"max_depth"]<-hyper_opt$max.depth
          dat_sample_i[,"eta"]<-hyper_opt$eta
          dat_sample_i[,"ntree"]<-hyper_opt$steps
          dat_sample_i[,"fs_num"]<-feat_sel_k
          dat_sample_i[,"fs"]<-fs_mth[fs]
          
          #only record validation results for final comparison, NO intermediate decisions are made based on them
          dat_sample_i[(dat_sample_i$part73=="T"),"pred"]<-as.numeric(predict(xgb_tune,dtrain))
          dat_sample_i[(dat_sample_i$part73!="T"),"pred"]<-as.numeric(predict(xgb_tune,dtest))
          
          # evaluate auc improvement
          ROC_obj_new<-roc_lst[[which.max(hyper_perf$metric_avg)]] #take the ROC curve based on the optimal model
          
          #need to update opt?
          if(ROC_obj_new$auc > ROC_obj_opt$auc){
            ROC_obj_opt<-ROC_obj_new
            opt_size<-feat_sel_k
            ROC_opt_update<-ROC_opt_update+1
          }
          
          # save everything about this senario, in case being called later
          feat_num[[paste0("size_",feat_sel_k)]]<-list(model_summary=dat_sample_i,
                                                       roc_obj=ROC_obj_new,
                                                       ROC_opt_update=ROC_opt_update)
          
          time_perf_i_nm<-c(time_perf_i_nm,paste0("tune_train_predict@",fs,"@",feat_sel_k))
          time_perf_i<-c(time_perf_i,paste0(Sys.time()-start_k,units(Sys.time()-start_k)))
          cat("......finish tuning, training and predicting in",time_perf_i[length(time_perf_i)],"\n")
          ##########################################################################################################
          
          time_perf_i_nm<-c(time_perf_i_nm,paste0("finish_eval_feature_num@",fs_mth[fs],"@",feat_sel_k))
          time_perf_i<-c(time_perf_i,paste0(Sys.time()-start_j,units(Sys.time()-start_j)))
          cat("...finish modeling on",feat_sel_k,"features for",fs_mth[fs],"in",time_perf_i[length(time_perf_i)],"\n")
        }
      }
      
      #compare b(2),c(3)
      auc_b_c<-feat_num[[paste0("size_",bounds[2])]]$roc_obj$auc-feat_num[[paste0("size_",bounds[3])]]$roc_obj$auc
      aucp_b_c<-pROC::roc.test(feat_num[[paste0("size_",bounds[2])]]$roc_obj,feat_num[[paste0("size_",bounds[3])]]$roc_obj,method='delong')$p.value
      
      #compare b(2) with opt
      auc_b_opt<-feat_num[[paste0("size_",bounds[2])]]$roc_obj$auc-ROC_obj_opt$auc
      aucp_b_opt<-pROC::roc.test(feat_num[[paste0("size_",bounds[2])]]$roc_obj,ROC_obj_opt,method='delong')$p.value
      
      #compare c with opt
      auc_c_opt<-feat_num[[paste0("size_",bounds[3])]]$roc_obj$auc-ROC_obj_opt$auc
      aucp_c_opt<-pROC::roc.test(feat_num[[paste0("size_",bounds[3])]]$roc_obj,ROC_obj_opt,method='delong')$p.value
      
      #update a,b,c,d
      if((max(aucp_b_opt,aucp_c_opt) <= inc_tol_p &&
          max(auc_b_opt,auc_c_opt) < 0)){
        a<-b
        d<-opt_size
        local_min<<-opt_size
      }else if(aucp_b_c > inc_tol_p ||
               auc_b_c > 0){
        d<-c
        local_min<<-b
      }else if(aucp_c_opt > inc_tol_p ||
               auc_c_opt >= 0){
        a<-b
        local_min<<-c
      }else{
        stop("conditions are not exhaustive!")
      }
      
      #update global_min?
      if(local_min < global_min){
        global_min <- local_min
        min_model <<- xgb_tune
      }
      
      #report progress
      time_perf_i_nm<-c(time_perf_i_nm,paste0("shrink_search_interval_to_",a,"_",d))
      time_perf_i<-c(time_perf_i,paste0(Sys.time()-start_g,units(Sys.time()-start_g)))
      cat(fs_mth[fs],":shrink interval to[",a,",",d,"] \n")
    }
    
    #end inner loop
    feat_num$track_path<-track_path #record the search track
    feat_num$opt_model<-xgb_tune #only record the model with optimal feature size
    
    fs_summary[[paste0("resample",i,"@",fs_mth[fs])]]<-feat_num
    
    time_perf_i_nm<-c(time_perf_i_nm,paste0("completion_at_resample",i,"@",fs_mth[fs]))
    time_perf_i<-c(time_perf_i,paste0(Sys.time()-start_fs,units(Sys.time()-start_fs)))
    cat("...finish evaluating feature ensemble method:",fs_mth[fs],"in",time_perf_i[length(time_perf_i)],"\n")
  }
  
  #end outer loop
  time_perf_i_nm<-c(time_perf_i_nm,paste0("completion_at_resample",i))
  time_perf_i<-c(time_perf_i,paste0(Sys.time()-start_i,units(Sys.time()-start_i)))
  cat("finish evaluating resample",i,"in",time_perf_i[length(time_perf_i)],"\n")
  
  time_perf[[i]]<-data.frame(task=time_perf_i_nm,
                             time=time_perf_i)
}

# save results
save(fs_summary,file="xgb_fs_summary.Rdata")
save(time_perf,file="xgb_performance2_boot.Rdata")
