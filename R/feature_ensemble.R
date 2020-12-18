################################################################
# Feature Ensemble Module                                      #
# Parameters:                                                  #
# -- rank_lst:feature ranking list                             #
#    - list: each list can be a data.frame or matrix           #
#       with two columns (var_colnm, rank_colnm)               #
# -- var_colnm: column name for variable in each list          #
# -- rank_colnm: column name for rank in each list             #
# -- ensemble_mth: feature ranking ensemble methods            #
#    - mean: unweighted average of ranks                       #
#    - top_k: hard membership of a feature ranking among top k #                                 
#    - exp_k: soft membership of a  feature ranking among top k#                 
#    - wt_mean: mean adjusted by weights                       #
#    - wt_top_k: top_k adjusted by weights                     #
#    - wt_exp_k: exp_k adjusted by weights                     #
# -- wt: weights associated with each feature set '            #
#    - each row of rank_mt                                     #
# -- k: top k                                                  #
################################################################

feature_ensemble<-function(rank_lst,
                           var_colnm="var",
                           rank_colnm="rk",
                           ensemble_mth=c("mean",
                                          "top_k",
                                          "exp_k",
                                          "wt_mean",
                                          "wt_top_k",
                                          "wt_exp_k"),
                           wt=rep(1,ncol(rank_mt)),
                           k=NULL){
  
  p<-length(rank_lst) #number of different rankings
  varlst<-c()
  rk_stack<-c()
  for(i in 1:p){
    varlst_i<-setnames(rank_lst[[i]],
                       old=c(var_colnm,rank_colnm),
                       new=c("var","rk"))
    varlst<-unique(c(varlst_i$var))
    rk_stack %<>% 
      bind_rows(varlst_i %>% mutate(mod=paste0("v",i)))
  }
  n<-length(varlst) #number of distinct features ever got selected
  
  #===
  rank_df<-rk_stack %>%
    spread(mod,rk) %>%
    replace(is.na(.),Inf) %>%
    gather(mod,rk) %>%
    group_by(var) %>%
    mutate(rk=dense_rank(rk))
    
  #====
  if(grepl("_k",ensemble_mth)){
    if(is.null(k)){
      stop("need to specify k, i.e. number of top features of interest!")
    }
    
    if(grepl("top_k",ensemble_mth)){
      rank_df %<>%
        gather(var,rk) %>%
        mutate(rk=(rk<=k)) %>%
        spread(var,rk)
    }else if(grepl("exp_k",ensemble_mth)){
      rank_df %<>%
        gather(var,rk) %>%
        mutate(rk=(exp(-rank/k))) %>%
        spread(var,rk)
    }else{
      stop("current ensemble method is not supported!")
    }
  }
  
  #====
  if(grepl("^wt",ensemble_mth)){
    #weighted ensemble schemes
    if(p!=length(wt)){
      stop("length of wt should be the same as nrow(rank_mt)!")
    }
    
    #broadcast the weight vector
    wt_mt<-(wt/sum(wt)) %*% rep(1,n)
    rank_df<-rank_df %*% wt_mt
  }
  
  rank_df_pivot<-rank_df %>%
    gather(var,rk) %>%
    group_by(var) %>%
    summarise(agg_rk=mean(rk)) %>%
    ungroup %>%
    mutate(agg_rk=dense_rank(agg_rk))
  
  return(rank_df_pivot)
}