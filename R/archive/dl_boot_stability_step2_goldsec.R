##################################################################
# Step 2 -- estimate percentage of featrue selected              #
#        -- retrain with different subset of features            # 
#        -- test on both hold-in and hold-out set                #
#        -- record out-of-sample AUC and Kuncheva index          #
##################################################################

##set up
rm(list=ls()); gc()
source("./helper_functions.R")
# require_libraries("devtools")
# slam_url<-"https://cran.r-project.org/src/contrib/Archive/slam/slam_0.1-37.tar.gz"
# install_url(slam_url)
require_libraries(c( "Matrix"
                    ,"pROC"
                    ,"h2o"
                    ,"dplyr"
                    ,"tidyr"
                    ,"magrittr"
                    ,"data.table"
                    ,"slam"
))

##load data
resamples<-10 #random sampling param
boots<-20 #boostrapping
load("DKD_heron_pats_prep.Rdata")
load("DKD_heron_facts_prep.Rdata")

load(paste0("random_sample",resamples,".Rdata"))
load("dl_feature_boot.Rdata")
load("dl_prediction_boot.Rdata")
load("feature_dict.Rdata")

##global values
#can be set to range of values
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

# feat_sel_rt<-c(0.25,0.1,0.05,0.02,0.01,0.005) #feature selection ratio
# feat_sel_k<-length(unique(fact_stack$CONCEPT_CD))*feat_sel_rt
feat_sel_k_low<-3
feat_sel_k_up<-600
inc_tol_p<-0.01
s<-c(50,100) # number of top rank of interests
phi<-(1+sqrt(5))/2
overall_y_sort<-pat_tbl %>% semi_join(fact_stack,by="PATIENT_NUM") %>% 
  arrange(PATIENT_NUM) %>% dplyr::select(DKD_IND) %>% unlist

# ensemble feature, tune, train, predict and evaluate stability
#initialization
feature_rk<-list()
fs_summary<-list() 
time_perf<-list()

#initialize h2o
h2o.init(nthreads=-1)

#start experiment
for(i in 1:resamples){
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
                  dplyr::select(Feature,scaled_importance)%>%
                  mutate(rank=rank(-scaled_importance),boot=b))
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
                dplyr::select(Feature,scaled_importance) %>%
                mutate(single_rank=rank(-scaled_importance)) %>%
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
    ROC_obj_new<-pROC::roc(overall_y_sort[(dat_sample_i$part73=="T")],
                           sample(c(0,1),nrow(dat_sample_i[(dat_sample_i$part73=="T"),]),replace=T),
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
            dplyr::slice(1:feat_sel_k)  %>%
            bind_rows(list(Feature="DKD_IND",rk=0))
          
          x_sparse<-fact_stack %>% 
            semi_join(fs_sel, by=c("CONCEPT_CD"="Feature")) %>%    #subset features
            inner_join(dat_sample_i,by="PATIENT_NUM") %>% 
            dplyr::select(-part73) %>%
            bind_rows(pat_tbl %>%
                        semi_join(dat_sample_i,by="PATIENT_NUM") %>%
                        dplyr::select(-year) %>%
                        gather(CONCEPT_CD,NVAL_NUM,-PATIENT_NUM) %>%
                        semi_join(fs_sel, by=c("CONCEPT_CD"="Feature"))) %>%
            group_by(PATIENT_NUM) %>%
            long_to_sparse_matrix(.,id="PATIENT_NUM",
                                  variable="CONCEPT_CD",
                                  val="NVAL_NUM") #sort by patient_num
          
          if(nrow(x_sparse)<dim(dat_sample_i)[1]){
            dat_sample_i %<>%
              semi_join(data.frame(PATIENT_NUM=as.numeric(rownames(x_sparse_val))),
                        by="PATIENT_NUM") #shrink feature space will reduce training data size as well
          }
          
          #record real y values
          dat_sample_i[,"real"]<-x_sparse[,"DKD_IND"]
          
          time_perf_i_nm<-c(time_perf_i_nm,paste0("subset_feature_transform@",fs,"@",feat_sel_k))
          time_perf_i<-c(time_perf_i,paste0(Sys.time()-start_k,units(Sys.time()-start_k)))
          cat("......finish subsetting and transforming in",time_perf_i[length(time_perf_i)],"\n")
    
          #######################################################################################
          start_k<-Sys.time()
          cat("......separate training and testing sets \n")
    
          train_mt<-x_sparse[(dat_sample_i$part73=="T"),]
          # colnames(train_mt)<-colnames(x_sparse) #colname may be dropped when only one column is selected
      
          test_mt<-x_sparse[(dat_sample_i$part73!="T"),]
          colnames(test_mt)<-colnames(x_sparse) #colname may be dropped when only one column is selected
      
          time_perf_i_nm<-c(time_perf_i_nm,paste0("partition@",fs,"@",feat_sel_k))
          time_perf_i<-c(time_perf_i,paste0(Sys.time()-start_k,units(Sys.time()-start_k)))
          cat("......finish partitioning in",time_perf_i[length(time_perf_i)],"\n") 
      
          ######################################################################################
          start_k<-Sys.time()
          cat("......tune, train and predict \n")
          
          h2o.removeAll()
          
          #h2o bug: special characters in colum names will result in an additional row
          train_mt_cp<-train_mt
          colnames(train_mt_cp)<-paste0("C",1:ncol(train_mt))
          train_h2o<-as.h2o(as.matrix(train_mt_cp)) #column order doesn't change
          #h2o truncate column space at 100 when apply as.h2o
          #h2o.splitFrame would force index and print out all columns
          train_h2o2<-h2o.splitFrame(data=train_h2o,ratios=0.9999)
          train_h2o_tr<-train_h2o2[[1]]
          
          #h2o bug: special characters in colum names will result in an additional row
          test_mt_cp<-test_mt
          colnames(test_mt_cp)<-paste0("C",1:ncol(test_mt))
          test_h2o<-as.h2o(as.matrix(test_mt_cp)) #column order doesn't change
          pred_idx<-which(colnames(train_mt)!="DKD_IND")
          target_idx<-which(colnames(train_mt)=="DKD_IND")
      
          train_h2o_tr[,target_idx]<-as.factor(train_h2o_tr[,target_idx])
          dl_random_grid <- h2o.grid(
            algorithm="deeplearning",
            grid_id = "dl_grid_random",
            training_frame=train_h2o_tr,
            nfolds=3,
            keep_cross_validation_predictions =T,
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
          dlr_grid<-h2o.getGrid("dl_grid_random",sort_by="auc",decreasing=T)
          flr_opt_model<-h2o.getModel(dlr_grid@model_ids[[1]])
          flr_hyper_params<-flr_opt_model@parameters 
          
          dat_sample_i[,"hidden"]<-paste(flr_opt_model@parameters$hidden,collapse = "_")
          dat_sample_i[,"input_dropout_ratio"]<-flr_opt_model@parameters$input_dropout_ratio
          dat_sample_i[,"fs_num"]<-feat_sel_k
          dat_sample_i[,"fs"]<-fs_mth[fs]
          
          #only record validation results for final comparison, NO intermediate decisions are made based on them
          dat_sample_i[(dat_sample_i$part73=="T"),"pred"]<-as.data.frame(predict(flr_opt_model,train_h2o))$p1
          dat_sample_i[(dat_sample_i$part73!="T"),"pred"]<-as.data.frame(predict(flr_opt_model,test_h2o))$p1
    
          # evaluate auc improvement
          cvpreds<-h2o.getFrame(flr_opt_model@model[["cross_validation_holdout_predictions_frame_id"]][["name"]]) # row sorted, matched with h2o.auc(flr_opt_model,xval=T)
          ROC_obj_new<-pROC::roc(as.vector(train_h2o_tr[,target_idx]),
                                 as.data.frame(cvpreds)$p1)

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
        min_model <<- flr_opt_model
        min_decode<<-data.frame(feat_code=colnames(train_h2o_tr)[-target_idx],
                                feat_name=colnames(train_mt)[-target_idx],stringsAsFactors=F)
      }
      
      #report progress
      time_perf_i_nm<-c(time_perf_i_nm,paste0("shrink_search_interval_to_",a,"_",d))
      time_perf_i<-c(time_perf_i,paste0(Sys.time()-start_g,units(Sys.time()-start_g)))
      cat(fs_mth[fs],":shrink interval to[",a,",",d,"] \n")
    }
    
    #end inner loop
    feat_num$track_path<-track_path #record the search track
    feat_num$opt_model<-list(opt_model=min_model, #only record the model with optimal feature size
                             opt_model_decode=min_decode)
    
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
save(feature_rk,file="dl_feature_rk.Rdata")
save(fs_summary,file="dl_fs_summary.Rdata")
save(time_perf,file="dl_performance2_boot.Rdata")

h2o.shutdown(prompt = FALSE)
