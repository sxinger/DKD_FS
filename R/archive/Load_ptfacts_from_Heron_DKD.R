#### Load DKD cohort and facts from patient_dimension####
rm(list=ls()); gc()
source("./R/util.R")
source("./R/sql_util.R")
require_libraries(c( "DBI"
                    ,"ROracle"
                    ,"dplyr"
                    ,"tidyr"))

#set up connection with Oracle db
config<-read.csv('./config.csv')
c_connect<-dbConnect(Oracle(),config$username,config$password,config$access)

pat_tbl_Oracle<-"PATIENTS_DKD"

pat_sql<-paste0("select dkd.*
                       ,pat.birth_date
                       ,case when pat.vital_status_cd = 'y' then 1 else 0 end as death_ind
                       ,pat.death_date
                       ,case when pat.sex_cd = 'm' then 1 else 0 end as gender_male
                       ,case when pat.race_cd in ('declined','@','not used') then 'unknown' else pat.race_cd end as race
                       ,case when pat.religion_cd = '@' then 'unknown' else pat.religion_cd end as religion
                 from ",pat_tbl_Oracle," dkd
                 join blueherondata.patient_dimension pat
                 on dkd.patient_num = pat.patient_num")

pat_tbl<-dbGetQuery(c_connect, pat_sql)

save(pat_tbl, file="DKD_heron_pats.Rdata")
