####SQL utility functions####
source("~/dkd//helper_functions.R")
require_libraries(c("DBI",
                    "ROracle",
                    "dplyr"))

##Get cohort summaries over inclusion period##
#pred_period only supports "days"
summarize_target_shift<-function(conn,pat_tbl,col_id,col_time,col_target,pred_pt,pred_period){
  if(!inherits(pred_pt,'Date') & !inherits(pred_pt,'numeric')){
    warning("pred_dt has to be either an enrollment date or number of years since DM onset")
  }else if(inherits(pred_pt,'Date')){
    sql_query<-paste0("select ",col_id,", target_adj ",col_target,", ",col_time,"
                       from (select DM_date, ",col_id,", ",col_target,", ",col_time,", 
                             max(",col_target,") over (partition by ",col_id,") target_adj, ","
                             row_number() over (partition by ",col_id,", ",col_target," order by ",col_time," desc) rn 
                             from ",pat_tbl," 
                             where DM_date < to_date('",pred_pt,"','YYYY-MM-DD') and ",
                                   col_time," < to_date('", pred_pt+pred_period, "','YYYY-MM-DD'))
                       where ",col_target,"= target_adj and rn = 1 and ",
                               col_time," >= to_date('", pred_pt,"','YYYY-MM-DD')")
  }else{
    sql_query<-paste0("select ",col_id,", target_adj ",col_target,", ",col_time,"
                       from (select DM_date, ",col_id,", ",col_target,", ",col_time,", 
                             max(",col_target,") over (partition by ",col_id,") target_adj, ","
                             row_number() over (partition by ",col_id,", ",col_target," order by ",col_time," desc) rn
                             from ",pat_tbl,"
                             where ",col_time," < DM_date + to_number(", pred_pt+pred_period, "))
                       where ",col_target,"= target_adj and rn = 1 and ",
                               col_time," >= DM_date + to_number(", pred_pt,")")
  }

  summ_pat<-dbGetQuery(conn, sql_query) %>%
    summarise(patient_cnt=n(),
              label_rate=round(mean(get(col_target)),2))
  return(summ_pat)
}

##write tables for patients of interest and encounters of interest
#to include the whole pat_tbl and all encounters between DMonset and endpoint, set:
# - pred_pt = 365*20;
# - pred_period = 0
send_pat_tbl<-function(conn,pat_tbl,pred_pt,pred_period){
  if(!inherits(pred_pt,'Date') & !inherits(pred_pt,'numeric')){
    warning("pred_dt has to be either an enrollment date or number of years since DM2 onset")
  }else if(inherits(pred_pt,'Date')){
    sql_pat<-paste0("create table ELIGIBLE_PATIENTS as
                     select x.patient_num, x.DM_date, x.DKD_ind, x.end_date
                     from (select patient_num, DM_date, DKD_ind, end_date,
                                  max(DKD_ind) over (partition by patient_num) DKD_ind_adj,
                                  row_number() over (partition by patient_num, DKD_ind order by end_date desc) rn
                           from ",pat_tbl,"
                           where DM_date < to_date('",pred_pt,"','YYYY-MM-DD') and
                                 end_date < to_date('", pred_pt+pred_period, "','YYYY-MM-DD')) x
                           join blueherondata.patient_dimension pat
                           on pat.patient_num = x.patient_num
                           where x.DKD_ind_adj = x.DKD_ind and x.rn = 1 and 
                                 x.end_date >= to_date('", pred_pt,"','YYYY-MM-DD')")
  }else{
    sql_pat<-paste0("create table ELIGIBLE_PATIENTS as
                     select x.patient_num, x.DM_date, x.DKD_ind, x.end_date
                     from (select patient_num, DM_date, DKD_ind, end_date,
                                  max(DKD_ind) over (partition by patient_num) DKD_ind_adj,
                                  row_number() over (partition by patient_num, DKD_ind order by end_date desc) rn
                           from ",pat_tbl,"
                           where end_date < DM_date + to_number(",pred_pt+pred_period, ")) x 
                     join blueherondata.patient_dimension pat
                     on pat.patient_num = x.patient_num
                     where x.DKD_ind_adj = x.DKD_ind and x.rn = 1 and
                           x.end_date >= x.DM_date + to_number(", pred_pt,")")
  }
  
  if(dbExistsTable(conn,"ELIGIBLE_PATIENTS")){
    dbSendQuery(conn,"drop table ELIGIBLE_PATIENTS PURGE")
  }
  dbSendQuery(conn,sql_pat)
}

send_enc_tbl<-function(conn,pat_tbl,enc_tbl,pred_pt,pred_period){
  if(!inherits(pred_pt,'Date') & !inherits(pred_pt,'numeric')){
    warning("pred_dt has to be either an enrollment date or number of years since DM onset")
  }else{
    sql_enc<-paste0("create table ELIGIBLE_ENCOUNTERS as
                     select distinct pat.patient_num, pat.DKD_ind, pe.encounter_num, pe.enc_type, 
                            (pat.end_date-pe.start_date) days_to_end, 
                            (pe.start_date-pat.DM_date) days_since_DM
                     from ",pat_tbl," pat
                     left join ",enc_tbl," pe
                     on pe.patient_num = pat.patient_num
                     where pe.start_date < pat.end_date - to_number(",pred_pt,")")
  }
  if(dbExistsTable(conn,"ELIGIBLE_ENCOUNTERS")){
    dbSendQuery(conn,"drop table ELIGIBLE_ENCOUNTERS PURGE")
  }
  dbSendQuery(conn,sql_enc)
}


##Get the most recent facts from Heron root folders (except for medications, labs, flowsheet and demographics)
get_pat_demo_facts<-function(conn,root_folder,recency=1,freq_bd=50,
                             restr=c("ELIGIBLE_PATIENTS","ELIGIBLE_ENCOUNTERS",NULL),...){
  if(!is.null(restr)){
    restr_arg<-paste0("and exists (select 1 from ",restr," res 
                                   where obs.patient_num = res.patient_num)")
  }else{
    restr_arg<-""
  }
  sql_facts<-paste0("select patient_num, concept_cd, c_name, variable_categ, c_visual_path, c_hlevel, nval_num, start_date,
                     count(distinct patient_num) over (partition by concept_cd) cd_total from
                     (select /*+leading(res,obs,ht) parallel(8)*/
                              obs.patient_num, obs.concept_cd, obs.nval_num,
                              obs.units_cd, obs.modifier_cd, obs.start_date,
                              row_number() over (partition by obs.patient_num, obs.concept_cd order by obs.start_date desc) recency, 
                              ht.variable_categ, ht.c_name, ht.c_visual_path, ht.c_hlevel
                       from blueherondata.observation_fact obs
                       join heron_terms ht
                       on obs.concept_cd = ht.c_basecode
                       where ht.variable_categ = '",root_folder, "' and
                             obs.concept_cd not like 'DEM|AGE%' and
                             obs.concept_cd not like 'DEM|SEX%' and
                             obs.concept_cd not like 'DEM|RACE%' and
                             obs.concept_cd not like 'HICTR_PARTICIPANT%' ",
                             restr_arg,")
                      where recency<=",recency)
  heron_facts<-dbGetQuery(conn,sql_facts) %>% dplyr::filter(CD_TOTAL > freq_bd)
  return(heron_facts)
}


get_recent_Heron_facts<-function(conn,root_folder,modifier=F,recency=1,freq_bd=50,
                                 restr=c("ELIGIBLE_PATIENTS","ELIGIBLE_ENCOUNTERS",NULL),pred_pt=NULL,...){
  if(modifier){
    concept_cd<-"(obs.concept_cd || '@' || obs.modifier_cd) as concept_cd"
    partition_by<-"obs.patient_num, obs.concept_cd, obs.modifier_cd"
  }else{
    concept_cd<-"obs.concept_cd"
    partition_by<-"obs.patient_num, obs.concept_cd"
  }
  
  if(restr=="ELIGIBLE_ENCOUNTERS"){
    restr_arg<-paste0("and exists (select 1 from ",restr," res 
                                   where obs.encounter_num = res.encounter_num)")
  }else if(restr=="ELIGIBLE_PATIENTS"){
    restr_arg<-paste0("and exists (select 1 from ",restr," res 
                                   where obs.patient_num = res.patient_num and 
                                         trunc(obs.start_date) < res.end_date-to_number(",pred_pt,"))")
  }else{
    restr_arg<-""
  }
  
  sql_facts<-paste0("select patient_num, concept_cd, c_name, variable_categ, c_visual_path, c_hlevel, nval_num, start_date,
                     count(distinct patient_num) over (partition by concept_cd) cd_total from
                     (select /*+leading(res,obs,ht) parallel(8)*/
                              obs.patient_num, ",concept_cd,", obs.nval_num,
                              obs.units_cd, obs.modifier_cd, obs.start_date,
                              row_number() over (partition by ",partition_by," order by obs.start_date desc) recency, 
                              ht.variable_categ, ht.c_name, ht.c_visual_path, ht.c_hlevel
                       from blueherondata.observation_fact obs
                       join heron_terms ht
                       on obs.concept_cd = ht.c_basecode
                       where ht.variable_categ = '",root_folder,"' ",restr_arg,")
                      where recency<=",recency)
  heron_facts<-dbGetQuery(conn,sql_facts) %>% dplyr::filter(CD_TOTAL > freq_bd)
  return(heron_facts)
}

##Get the most recent Lab facts from Heron
get_recent_LAB_facts<-function(conn,root_folder,recency=1,freq_bd=50,
                               restr=c("ELIGIBLE_PATIENTS","ELIGIBLE_ENCOUNTERS",NULL),pred_pt=NULL,...){
  if(restr=="ELIGIBLE_ENCOUNTERS"){
    restr_arg<-paste0("and exists (select 1 from ",restr," res 
                                   where obs.encounter_num = res.encounter_num)")
  }else if(restr=="ELIGIBLE_PATIENTS"){
    restr_arg<-paste0("and exists (select 1 from ",restr," res 
                                   where obs.patient_num = res.patient_num and 
                                         trunc(obs.start_date) < res.end_date-to_number(",pred_pt,"))")
  }else{
    restr_arg<-""
  }
  
  sql_facts<-paste0("select patient_num, concept_cd, c_name, variable_categ, c_visual_path, c_hlevel, nval_num, start_date,
                     count(distinct patient_num) over (partition by concept_cd) cd_total from
                     (select /*+leading(enc,obs,ht) parallel(8)*/ obs.patient_num, 
                             (ht.c_basecode || '#' || upper(obs.units_cd) || '@' || obs.modifier_cd) as concept_cd,
                             obs.nval_num, obs.units_cd, obs.modifier_cd, obs.start_date,
                             row_number() over (partition by obs.patient_num, obs.concept_cd, obs.units_cd, obs.modifier_cd order by obs.start_date desc) recency, 
                             ht.variable_categ, ht.c_name, ht.c_visual_path, ht.c_hlevel
                       from blueherondata.observation_fact obs
                       join heron_terms ht
                       on obs.concept_cd = ht.c_basecode
                       where ht.variable_categ = '",root_folder,"' ",restr_arg," and
                             obs.modifier_cd <> '@')
                      where recency<=",recency)
  heron_facts<-dbGetQuery(conn,sql_facts) %>% dplyr::filter(CD_TOTAL > freq_bd)
  return(heron_facts)
}

##Get the most recent Flowsheet facts from Heron
get_recent_FLOSHEET_facts<-function(conn,root_folder,recency=1,freq_bd=50,
                                    restr=c("ELIGIBLE_PATIENTS","ELIGIBLE_ENCOUNTERS",NULL),pred_pt=NULL,...){
  if(restr=="ELIGIBLE_ENCOUNTERS"){
    restr_arg<-paste0("and exists (select 1 from ",restr," res 
                                   where obs.encounter_num = res.encounter_num)")
  }else if(restr=="ELIGIBLE_PATIENTS"){
    restr_arg<-paste0("and exists (select 1 from ",restr," res 
                                   where obs.patient_num = res.patient_num and 
                                         trunc(obs.start_date) < res.end_date-to_number(",pred_pt,"))")
  }else{
    restr_arg<-""
  }
  
  sql_facts<-paste0("select patient_num, concept_cd, c_name, variable_categ, c_visual_path, c_hlevel, nval_num, start_date,
                     count(distinct patient_num) over (partition by concept_cd) cd_total from
                     (select /*+leading(enc,obs,ht) parallel(8)*/ obs.patient_num, 
                             (ht.c_basecode || '#' || upper(obs.units_cd)) as concept_cd,
                             obs.nval_num, obs.units_cd, obs.start_date,
                             row_number() over (partition by obs.patient_num, obs.concept_cd, obs.units_cd order by obs.start_date desc) recency, 
                             ht.variable_categ, ht.c_name, ht.c_visual_path, ht.c_hlevel
                       from blueherondata.observation_fact obs
                       join heron_terms ht
                       on obs.concept_cd = ht.c_basecode
                       where ht.variable_categ = '",root_folder,"' ",restr_arg,")
                      where recency<=",recency)
  heron_facts<-dbGetQuery(conn,sql_facts) %>% dplyr::filter(CD_TOTAL > freq_bd)
  return(heron_facts)
}


##Get the most recent MEDICATIONS facts from Heron (with grouping of modifiers)
get_recent_MED_facts<-function(conn,ip_id_mod,ip_mod,op_mod,recency=1,freq_bd=50,
                               restr=c("ELIGIBLE_PATIENTS","ELIGIBLE_ENCOUNTERS",NULL),pred_pt=NULL,...){
  if(restr=="ELIGIBLE_ENCOUNTERS"){
    restr_arg<-paste0("and exists (select 1 from ",restr," res 
                                   where obs.encounter_num = res.encounter_num)")
  }else if(restr=="ELIGIBLE_PATIENTS"){
    restr_arg<-paste0("and exists (select 1 from ",restr," res 
                                   where obs.patient_num = res.patient_num and 
                                         trunc(obs.start_date) < res.end_date-to_number(",pred_pt,"))")
  }else{
    restr_arg<-""
  }
  
  # create temp table collecting all medication info for all encounters of interest
  sql_medins<-paste0("create table MED_FULL as
                      select /*+leading(res,obs,ht) parallel(8)*/
                             obs.patient_num, obs.encounter_num, obs.instance_num, obs.concept_cd, 
                             ht.c_name, ht.variable_categ, ht.c_visual_path, ht.c_hlevel,
                             obs.modifier_cd, obs.nval_num, obs.units_cd, obs.start_date
                      from blueherondata.observation_fact obs
                      join heron_terms ht
                      on obs.concept_cd = ht.c_basecode
                      where ht.variable_categ = 'MEDICATIONS' ",restr_arg)
  
  if(dbExistsTable(conn,"MED_FULL")){
    dbSendQuery(conn,"drop table MED_FULL PURGE")
  }
  dbSendQuery(conn,sql_medins)
  
  sql_query_med<-function(in_med=T,mod){
    paste0("with MED_TERMS as
           (select distinct obs.concept_cd, ht.c_name, ht.variable_categ, ht.c_visual_path
            from MED_FULL obs
            join heron_terms ht
            on obs.concept_cd = ht.c_basecode
            where",ifelse(in_med," "," not "),"regexp_like(obs.modifier_cd,'((",paste(ip_id_mod,collapse = ")|("),"))+')),
         
                 MED_INS as
            (select distinct obs.instance_num, obs.concept_cd
            from MED_FULL obs
            join heron_terms ht
            on obs.concept_cd = ht.c_basecode
            where",ifelse(in_med," "," not "),"regexp_like(obs.modifier_cd,'((",paste(ip_id_mod,collapse = ")|("),"))+'))

            select /*+leading(ins,obs,ht) parallel(6)*/
                   distinct x.patient_num, x.concept_cd, x.c_name, x.variable_categ, x.c_visual_path, 
                            x.nval_num, x.start_date, x.med_cnt,
                            count(distinct x.patient_num) over (partition by x.concept_cd) cd_total 
            from (select obs.patient_num, obs.encounter_num,
                         (ht.concept_cd || '#' || upper(obs.units_cd) || '@' || obs.modifier_cd) as concept_cd, 
                         ht.c_name, ht.variable_categ, ht.c_visual_path, obs.nval_num, obs.start_date,
                         count(distinct obs.instance_num) over (partition by obs.patient_num, ht.concept_cd) med_cnt,
                         row_number() over (partition by obs.patient_num, ht.concept_cd, obs.units_cd, obs.modifier_cd order by obs.start_date desc) recency
                  from MED_FULL obs
                  join MED_TERMS ht
                  on obs.concept_cd = ht.concept_cd
                  where obs.nval_num is not null and
                        regexp_like(obs.modifier_cd, '(",paste(mod,collapse = ")|("),")+') and
                        exists (select 1 from MED_INS ins where ins.instance_num = obs.instance_num)) x
             where x.recency<=",recency, "
             union
             select /*+leading(ins,obs,ht) parallel(6)*/
                     distinct y.patient_num, y.concept_cd, y.c_name, y.variable_categ, y.c_visual_path, 
                              1 as nval_num, y.start_date, y.med_cnt,
                              count(distinct y.patient_num) over (partition by y.concept_cd) cd_total 
             from (select obs.patient_num, obs.encounter_num, 
                          (ht.concept_cd || '@' || '",ifelse(in_med,"In","Out"),"Given') as concept_cd,
                          ht.c_name, ht.variable_categ, ht.c_visual_path, obs.start_date,
                          count(distinct obs.instance_num) over (partition by obs.patient_num, ht.concept_cd) med_cnt,
                          row_number() over (partition by obs.patient_num, ht.concept_cd order by obs.start_date desc) recency
                   from MED_FULL obs
                   join MED_TERMS ht
                   on obs.concept_cd = ht.concept_cd
                   where obs.nval_num is null and
                         not regexp_like(obs.modifier_cd, '(",paste(mod,collapse = ")|("),")+') and
                         exists (select 1 from MED_INS ins where ins.instance_num = obs.instance_num)) y
              where y.recency<=",recency)
  }
  
  heron_facts<-bind_rows(dbGetQuery(conn,sql_query_med(in_med=T,ip_mod)) %>% dplyr::filter(CD_TOTAL > freq_bd),
                         dbGetQuery(conn,sql_query_med(in_med=F,op_mod)) %>% dplyr::filter(CD_TOTAL > freq_bd))
  return(heron_facts)
}
