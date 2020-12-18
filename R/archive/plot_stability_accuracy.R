######################################
# Accuracy and Stability Measurement #
######################################
setwd("~/dkd")

rm(list=ls()); gc()
source("./helper_functions.R")
require_libraries(c(
                     "pROC"
                    ,"dplyr"
                    ,"tidyr"
                    ,"magrittr"
                    ,"ggplot2"
                    ))
# models under considerations
model_type<-c("lasso",
              "elastnet",
              "xgb",
              "dl"
              )

##################################################################################################

##########################################
#########plot accuracy vs stability#######
##########################################
x_min<-y_min<-1
x_max<-y_max<-0

stable_accuracy<-c()
for(i in seq_along(model_type)){
  load(paste0(model_type[i],"accuracy_stability_opt.Rdata"))
  accu_stab_opt %<>%
    group_by(fs_mth) %>%
    filter(geo_root==max(geo_root)) %>%
    dplyr::mutate(model=model_type[i]) %>%
    dplyr::mutate(label=paste0(model,"+",fs_mth,"+",fs_num)) %>%
    ungroup 
  
  stable_accuracy %<>%
    bind_rows(accu_stab_opt)
  
  x_min<-min(x_min,max(min(accu_stab_opt$auc)-0.01,0))
  x_max<-max(x_max,min(max(accu_stab_opt$auc)+0.01,1))
  
  y_min<-min(y_min,max(min(accu_stab_opt$cw_rel)-0.01,0))
  y_max<-max(y_max,min(max(accu_stab_opt$cw_rel)+0.01,1))
}

size_val<-c(100,200,300,400,500,600)
shape_val<-1:nlevels(as.factor(stable_accuracy$fs_mth))
plot_st_ac_v<-ggplot(stable_accuracy,aes(x=auc,y=cw_rel,
                                         color=fs_mth,shape=fs_mth))+
  geom_point(aes(size=fs_num),stroke=1)+
  scale_color_discrete("feature ensemble")+
  scale_shape_manual("feature ensemble",values=shape_val)+
  scale_size_continuous("feature space",breaks=size_val,labels=size_val)+
  scale_x_continuous("AUC",
                     limits=c(x_min,x_max),breaks=seq(x_min,x_max,by=0.02))+
  scale_y_continuous("Relative Weighted Consistency Index",
                     limits=c(y_min,y_max),breaks=seq(y_min,y_max,by=0.1))+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=90,hjust=1))+
  labs(title="Stability vs Accuracy")+
  annotate("rect",alpha=0.2,
           xmin=auc-0.004,xmax=auc+0.005,
           ymin=cw_rel-0.04,ymax=cw_rel+0.04)+
  annotate("text",vjust=-1.4,hjust=0.8,
           x=auc, y=cw_rel,  label=label)+
  geom_hline(yintercept = cw_rel,linetype="dashed",color="black")+
  geom_vline(xintercept = auc,linetype="dashed",color="black")
#############################################################################################################



###########################################################
############### plot for individual method ################
###########################################################
#xgb
i<-4
load(paste0(model_type[i],"_accuracy_summary.Rdata"))


x_scale<-unique(accuracy_stack$fs_num)
x_scale<-x_scale[order(x_scale)]
shape_val<-1:nlevels(as.factor(accuracy_stack$fs))
y_min<-max(round(min(accuracy_stack[(accuracy_stack$part73!="T"),]$auc_90low),1)-0.05,0.5)
y_max<-min(round(max(accuracy_stack[(accuracy_stack$part73!="T"),]$auc_90up),1)+0.05,1)

plot_auc_v<-ggplot(accuracy_stack %>% dplyr::filter(part73=="V"),
       aes(x=fs_num_grp,y=auc_med,colour=fs,shape=fs,linetype=fs))+
  geom_line()+geom_point(size=3)+
  geom_errorbar(aes(ymin=auc_90low,ymax=auc_90up),width=5)+
  scale_color_discrete("feature ensemble")+
  scale_shape_manual("feature ensemble",values=shape_val)+
  scale_linetype_manual("feature ensemble",values=shape_val)+
  scale_x_continuous(labels=x_scale,breaks=x_scale)+
  scale_y_continuous(limits=c(y_min,y_max),breaks=seq(y_min,y_max,by=0.05))+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=90,hjust=1))+
  labs(x="number of features (k*1.5)",y="AUC",
       title=paste0(model_type[i],":Internal Validation"))

plot_auc_h<-ggplot(accuracy_stack %>% dplyr::filter(part73=="H"),
                   aes(x=fs_num_grp,y=auc_med,colour=fs,shape=fs,linetype=fs))+
  geom_line()+ geom_point(size=3)+
  geom_errorbar(aes(ymin=auc_90low,ymax=auc_90up),width=5)+
  scale_color_discrete("feature ensemble")+
  scale_shape_manual("feature ensemble",values=shape_val)+
  scale_linetype_manual("feature ensemble",values=shape_val)+
  scale_x_continuous(labels=x_scale,breaks=x_scale)+
  scale_y_continuous(limits=c(y_min,y_max),breaks=seq(y_min,y_max,by=0.05))+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=90,hjust=1))+
  labs(x="number of features (k*1.5)",y="AUC",
       title=paste0(model_type[i],":Temporal External Validation"))

grid_arrange_shared_legend(list(auc_v=plot_auc_v,
                                auc_h=plot_auc_h),
                           ncol=2,
                           nrow=1)

#Stability
load(paste0(model_type[i],"_stability_summary.Rdata"))
cw_idx<-stable_stack$cw_idx %>%
  arrange(fs,fs_num) %>%
  group_by(fs) %>%
  mutate(fs_num_grp = cut(fs_num,breaks=brks,include.lowest=T,labels=brks))

gk_kch_idx<-stable_stack$gk_kch_idx %>%
  arrange(fs,fs_num) %>%
  group_by(fs) %>%
  mutate(fs_num_grp = cut(fs_num,breaks=brks,include.lowest=T,labels=brks))

x_scale<-unique(cw_idx$fs_num)
x_scale<-x_scale[order(x_scale)]
brks<-c(2,80,100,120,140,160,180,200,300,400,500,600)
linetype_val<-1:nlevels(as.factor(cw_idx$fs_mth))
y_min<-max(min(round(min(cw_idx$cw),1)-0.1,
               round(min(cw_idx$cw_rel),1)-0.1,
               round(min(gk_kch_idx$ki),1)-0.1,
               round(min(gk_kch_idx$gk_idx),1)-0.1),0)

y_max<-max(round(max(cw_idx$cw),1)+0.1,
           round(max(cw_idx$cw_rel),1)+0.1,
           round(max(gk_kch_idx$ki),1)+0.1,
           round(min(gk_kch_idx$gk_idx),1)+0.1)

plot_cw<-ggplot(cw_idx,aes(x=fs_num,y=cw,color=fs_mth,shape=fs_mth,linetype=fs_mth))+
  geom_line()+geom_point(size=3)+
  scale_color_discrete("feature ensemble")+
  scale_shape_manual("feature ensemble",values=linetype_val) +
  scale_linetype_manual("feature ensemble",values=linetype_val) +
  scale_x_continuous(labels=x_scale,breaks=x_scale)+
  scale_y_continuous(limits=c(y_min,y_max),breaks=seq(y_min,y_max,by=0.1))+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        title =element_text(size=8, face='bold'),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=90,hjust=1))+
  labs(x="number of features (k*1.5)",y="CW Index",
       title=paste0(model_type[i],":Weighted Consistency Index"))
plot_stability[[paste0(model_type[i],"_cw")]]<-plot_cw

plot_cwr<-ggplot(cw_idx,aes(x=fs_num,y=cw_rel,color=fs_mth,shape=fs_mth, linetype=fs_mth))+
  geom_line()+geom_point(size=3)+
  scale_color_discrete("feature ensemble")+
  scale_shape_manual("feature ensemble",values=linetype_val) +
  scale_linetype_manual("feature ensemble",values=linetype_val) +
  scale_x_continuous(labels=x_scale,breaks=x_scale)+
  scale_y_continuous(limits=c(y_min,y_max),breaks=seq(y_min,y_max,by=0.1))+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        title =element_text(size=8, face='bold'),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=90,hjust=1))+
  labs(x="number of features (k*1.5)",y="Relative CW Index",
       title=paste0(model_type[i],":Relative Weighted Consistency Index"))
plot_stability[[paste0(model_type[i],"_cwr")]]<-plot_cwr

plot_ki<-ggplot(gk_kch_idx,aes(x=fs_num,y=ki,color=fs_mth,shape=fs_mth,linetype=fs_mth))+
  geom_line()+geom_point(size=3)+
  scale_color_discrete("feature ensemble")+
  scale_shape_manual("feature ensemble",values=linetype_val) +
  scale_linetype_manual("feature ensemble",values=linetype_val) +
  scale_x_continuous(labels=x_scale,breaks=x_scale)+
  scale_y_continuous(limits=c(y_min,y_max),breaks=seq(y_min,y_max,by=0.1))+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        title =element_text(size=8, face='bold'),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=90,hjust=1))+
  labs(x="number of features (k*1.5)",y="Kuncheva Index",
       title=paste0(model_type[i],":Kuncheva Index"))
plot_stability[[paste0(model_type[i],"_ki")]]<-plot_ki

plot_gk<-ggplot(gk_kch_idx,aes(x=fs_num,y=gk_idx,color=fs_mth,shape=fs_mth, linetype=fs_mth))+
  geom_line()+geom_point(size=3)+
  scale_color_discrete("feature ensemble")+
  scale_shape_manual("feature ensemble",values=linetype_val) +
  scale_linetype_manual("feature ensemble",values=linetype_val) +
  scale_x_continuous(labels=x_scale,breaks=x_scale)+
  scale_y_continuous(limits=c(y_min,y_max),breaks=seq(y_min,y_max,by=0.1))+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        title =element_text(size=8, face='bold'),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=90,hjust=1))+
  labs(x="number of features (k*1.5)",y="GK Index",
       title=paste0(model_type[i],":General Kalousis Index"))

grid_arrange_shared_legend(list(cw=plot_cw,
                                cwr=plot_cwr,
                                kch=plot_ki,
                                gk=plot_gk),
                           ncol=4,nrow=1)

#AUC vs Stability
paste0(model_type[i],"accuracy_stability_opt.Rdata")

opt_pair<-c()
opt_pair %<>%
  bind_rows(accu_stab_opt %>%
              mutate(penalty="euclid",bal=euclid) %>%
              filter(euclid==max(euclid)) %>%
              select(fs_mth,fs_num,auc_med,cw_rel,penalty,bal) %>%
              mutate(label=paste0(fs_mth,"_",fs_num))) %>%
  bind_rows(accu_stab_opt %>%
              mutate(penalty="geo_root",bal=geo_root)  %>% 
              filter(geo_root==max(geo_root)) %>%
              select(fs_mth,fs_num,auc_med,cw_rel,penalty,bal) %>%
              mutate(label=paste0(fs_mth,"_",fs_num)))

shape_val<-1:nlevels(as.factor(stable_accuracy$fs_mth))
y_min<-max(round(min(accu_stab_opt$cw_rel),1),0.1)
y_max<-round(max(accu_stab_opt$cw_rel),1)
x_min<-max(round(min(accu_stab_opt$auc),2)-0.01,0)
x_max<-min(round(max(accu_stab_opt$auc),2)+0.01,1)

plot_st_ac_v<-ggplot(accu_stab_opt,
                     aes(x=auc,y=cw_rel,color=fs_mth,shape=fs_mth))+
  geom_point(aes(size=fs_num),stroke=1)+
  scale_color_discrete("feature ensemble")+
  scale_shape_manual("feature ensemble",values=shape_val)+
  scale_size_continuous("feature space",breaks=size_val,labels=size_val)+
  scale_x_continuous("AUC",
                     limits=c(x_min,x_max),breaks=seq(x_min,x_max,by=0.02))+
  scale_y_continuous("Relative Weighted Consistency Index",
                     limits=c(y_min,y_max),breaks=seq(y_min,y_max,by=0.1))+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=90,hjust=1))+
  labs(title=paste0(model_type[i],":Stability vs Accuracy"))+
  annotate("rect",alpha=0.2,
           xmin=opt_pair[2,]$auc_med-0.004,xmax=opt_pair[2,]$auc_med+0.005,
           ymin=opt_pair[2,]$cw_rel-0.04,ymax=opt_pair[2,]$cw_rel+0.04)+
  annotate("text",vjust=-1.4,hjust=0.8,
           x=opt_pair[2,]$auc_med,
           y=opt_pair[2,]$cw_rel,
           label=opt_pair[2,]$label)+
  geom_hline(yintercept = opt_pair[2,]$cw_rel,linetype="dashed",color="black")+
  geom_vline(xintercept = opt_pair[2,]$auc_med,linetype="dashed",color="black")

plot_st_ac_v

#################################################
############### plot AUC-overall ################
#################################################
#adapt y-scale to data
y_min<-1
y_max<-0
for(i in 1:4){
  load(paste0(model_type[i],"_accuracy_summary.Rdata"))
  accuracy_stack %<>% dplyr::filter(fs_num<=300)
  
  y_min<-min(y_min,max(round(min(accuracy_stack[(accuracy_stack$part73!="T"),]$auc_90low),1)-0.05,0.5))
  y_max<-max(y_max,min(round(max(accuracy_stack[(accuracy_stack$part73!="T"),]$auc_90up),1)+0.05,1))
}

#collect accuracy plots
plot_accuracy<-list()

for(i in 1:4){
  load(paste0(model_type[i],"_accuracy_summary.Rdata"))
  accuracy_stack %<>%
    dplyr::filter(fs_num<=300)
  
  x_scale<-unique(accuracy_stack$fs_num)
  shape_val<-1:nlevels(as.factor(accuracy_stack$fs))
  
  plot_auc_v<-ggplot(accuracy_stack %>% dplyr::filter(part73=="V"),
                     aes(x=fs_num,y=auc_med,colour=fs,shape=fs,linetype=fs))+
    geom_line()+geom_point(size=3)+
    geom_errorbar(aes(ymin=auc_90low,ymax=auc_90up),width=5)+
    scale_color_discrete("feature ensemble")+
    scale_shape_manual("feature ensemble",values=shape_val)+
    scale_linetype_manual("feature ensemble",values=shape_val)+
    scale_x_continuous(labels=x_scale,breaks=x_scale)+
    scale_y_continuous(limits=c(y_min,y_max),breaks=seq(y_min,y_max,by=0.05))+
    theme_bw()+
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x=element_text(angle=90,hjust=1))+
    labs(x="number of features (k*1.5)",y="AUC",
         title=paste0(model_type[i],":Internal Validation"))
  
  plot_accuracy[[paste0(model_type[i],"_auc_v")]]<-plot_auc_v
  
  
  plot_auc_h<-ggplot(accuracy_stack %>% dplyr::filter(part73=="H"),
                     aes(x=fs_num,y=auc_med,colour=fs,shape=fs,linetype=fs))+
    geom_line()+ geom_point(size=3)+
    geom_errorbar(aes(ymin=auc_90low,ymax=auc_90up),width=5)+
    scale_color_discrete("feature ensemble")+
    scale_shape_manual("feature ensemble",values=shape_val)+
    scale_linetype_manual("feature ensemble",values=shape_val)+
    scale_x_continuous(labels=x_scale,breaks=x_scale)+
    scale_y_continuous(limits=c(y_min,y_max),breaks=seq(y_min,y_max,by=0.05))+
    theme_bw()+
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x=element_text(angle=90,hjust=1))+
    labs(x="number of features (k*1.5)",y="AUC",
         title=paste0(model_type[i],":Temporal External Validation"))
  
  plot_accuracy[[paste0(model_type[i],"_auc_h")]]<-plot_auc_h
  
}

# plot all
grid_arrange_shared_legend(plot_accuracy,
                           ncol=2,
                           nrow=4)
#################################################################################################


####################################################################
############### plot Sensitivity vs Specificity (V) ################
####################################################################
#adapt y-scale to data
y_min<-1
y_max<-0
for(i in 1:4){
  load(paste0(model_type[i],"_accuracy_summary.Rdata"))
  accuracy_stack %<>% dplyr::filter(fs_num<=300)
  
  y_min<-min(y_min,max(round(min(c(accuracy_stack[(accuracy_stack$part73=="V"),]$sens_90low,
                                   accuracy_stack[(accuracy_stack$part73=="V"),]$spec_90low)),1)-0.05,0))
  y_max<-max(y_max,min(round(max(c(accuracy_stack[(accuracy_stack$part73=="V"),]$sens_90up,
                                   accuracy_stack[(accuracy_stack$part73=="V"),]$spec_90up)),1)+0.05,1))
}

#collect accuracy plots
plot_sens_spec_v<-list()

for(i in 1:4){
  load(paste0(model_type[i],"_accuracy_summary.Rdata"))
  accuracy_stack %<>%
    dplyr::filter(fs_num<=300)
  
  x_scale<-unique(accuracy_stack$fs_num)
  linetype_val<-1:nlevels(as.factor(accuracy_stack$fs))
  
  plot_sens_v<-ggplot(accuracy_stack %>% dplyr::filter(part73=="V"),
                     aes(x=fs_num,y=auc_med,colour=fs,shape=fs,linetype=fs))+
    geom_line()+geom_point(size=3)+
    geom_errorbar(aes(ymin=auc_90low,ymax=auc_90up))+
    scale_color_discrete("feature ensemble")+
    scale_shape_manual("feature ensemble",values=linetype_val) +
    scale_linetype_manual("feature ensemble",values=linetype_val) +
    scale_x_continuous(labels=x_scale,breaks=x_scale)+
    scale_y_continuous(limits=c(y_min,y_max),breaks=seq(y_min,y_max,by=0.05))+
    theme_bw()+
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x=element_text(angle=90,hjust=1))+
    labs(x="number of features (k*1.5)",y="Sensitivity",
         title=paste0(model_type[i],":Internal Validation"))
  
  plot_sens_spec_v[[paste0(model_type[i],"_sens_v")]]<-plot_sens_v
  
  
  plot_spec_v<-ggplot(accuracy_stack %>% dplyr::filter(part73=="V"),
                     aes(x=fs_num,y=auc_med,colour=fs,shape=fs,linetype=fs))+
    geom_line()+ geom_point(size=3)+
    geom_errorbar(aes(ymin=auc_90low,ymax=auc_90up))+
    scale_color_discrete("feature ensemble")+
    scale_shape_manual("feature ensemble",values=linetype_val) +
    scale_linetype_manual("feature ensemble",values=linetype_val) +
    scale_x_continuous(labels=x_scale,breaks=x_scale)+
    scale_y_continuous(limits=c(y_min,y_max),breaks=seq(y_min,y_max,by=0.05))+
    theme_bw()+
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x=element_text(angle=90,hjust=1))+
    labs(x="number of features (k*1.5)",y="Specificity",
         title=paste0(model_type[i],":Internal Validation"))
  
  plot_sens_spec_v[[paste0(model_type[i],"_spec_v")]]<-plot_spec_v
  
}

# plot all
grid_arrange_shared_legend(plot_sens_spec_v,
                           ncol=2,
                           nrow=4)
##############################################################################################

####################################################################
############### plot Sensitivity vs Specificity (H) ################
####################################################################
#adapt y-scale to data
y_min<-1
y_max<-0
for(i in 1:4){
  load(paste0(model_type[i],"_accuracy_summary.Rdata"))
  accuracy_stack %<>% dplyr::filter(fs_num<=500)
  
  y_min<-min(y_min,max(round(min(c(accuracy_stack[(accuracy_stack$part73=="H"),]$sens_90low,
                                   accuracy_stack[(accuracy_stack$part73=="H"),]$spec_90low)),1)-0.05,0))
  y_max<-max(y_max,min(round(max(c(accuracy_stack[(accuracy_stack$part73=="H"),]$sens_90up,
                                   accuracy_stack[(accuracy_stack$part73=="H"),]$spec_90up)),1)+0.05,1))
}

#collect accuracy plots
plot_sens_spec_h<-list()

for(i in 1:4){
  load(paste0(model_type[i],"_accuracy_summary.Rdata"))
  accuracy_stack %<>%
    dplyr::filter(fs_num<=500)
  
  x_scale<-unique(accuracy_stack$fs_num)
  linetype_val<-1:nlevels(as.factor(accuracy_stack$fs))
  
  plot_sens_h<-ggplot(accuracy_stack %>% dplyr::filter(part73=="H"),
                      aes(x=fs_num,y=auc_med,colour=fs,shape=fs,linetype=fs))+
    geom_line()+geom_point(size=3)+
    geom_errorbar(aes(ymin=auc_90low,ymax=auc_90up))+
    scale_color_discrete("feature ensemble")+
    scale_shape_manual("feature ensemble",values=linetype_val) +
    scale_linetype_manual("feature ensemble",values=linetype_val) +
    scale_x_continuous(labels=x_scale,breaks=x_scale)+
    scale_y_continuous(limits=c(y_min,y_max),breaks=seq(y_min,y_max,by=0.05))+
    theme_bw()+
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x=element_text(angle=90,hjust=1))+
    labs(x="number of features (k*1.5)",y="Sensitivity",
         title=paste0(model_type[i],":Temporal External Validation"))
  
  plot_sens_spec_h[[paste0(model_type[i],"_sens_h")]]<-plot_sens_h
  
  
  plot_spec_h<-ggplot(accuracy_stack %>% dplyr::filter(part73=="H"),
                      aes(x=fs_num,y=auc_med,colour=fs,shape=fs,linetype=fs))+
    geom_line()+ geom_point(size=3)+
    geom_errorbar(aes(ymin=auc_90low,ymax=auc_90up))+
    scale_color_discrete("feature ensemble")+
    scale_shape_manual("feature ensemble",values=linetype_val) +
    scale_linetype_manual("feature ensemble",values=linetype_val) +
    scale_x_continuous(labels=x_scale,breaks=x_scale)+
    scale_y_continuous(limits=c(y_min,y_max),breaks=seq(y_min,y_max,by=0.05))+
    theme_bw()+
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x=element_text(angle=90,hjust=1))+
    labs(x="number of features (k*1.5)",y="Specificity",
         title=paste0(model_type[i],":Temporal External Validation"))
  
  plot_sens_spec_h[[paste0(model_type[i],"_spec_h")]]<-plot_spec_h
  
}

# plot all
grid_arrange_shared_legend(plot_sens_spec_h,
                           ncol=2,
                           nrow=4)
################################################################################################

########################################################
############### plot stability measures ################
########################################################
y_min<-1
y_max<-0
for(i in 1:4){
  load(paste0(model_type[i],"_stability_summary.Rdata"))
  cw_idx<-stable_stack$cw_idx %<>% dplyr::filter(fs_num<=300)
  gk_kch_idx<-stable_stack$gk_kch_idx %<>% dplyr::filter(fs_num<=300)
  
  y_min<-min(y_min,max(min(round(min(cw_idx$cw),1)-0.1,
                           round(min(cw_idx$cw_rel),1)-0.1,
                           round(min(gk_kch_idx$ki),1)-0.1,
                           round(min(gk_kch_idx$gk_idx),1)-0.1),0))
  
  y_max<-max(y_max,max(round(max(cw_idx$cw),1)+0.1,
                       round(max(cw_idx$cw_rel),1)+0.1,
                       round(max(gk_kch_idx$ki),1)+0.1,
                       round(min(gk_kch_idx$gk_idx),1)+0.1))
}



plot_stability<-list()
for(i in 1:4){
  load(paste0(model_type[i],"_stability_summary.Rdata"))
  cw_idx<-stable_stack$cw_idx %<>% dplyr::filter(fs_num<=300)
  gk_kch_idx<-stable_stack$gk_kch_idx %<>% dplyr::filter(fs_num<=300)
  
  x_scale<-unique(cw_idx$fs_num)
  linetype_val<-1:nlevels(as.factor(cw_idx$fs_mth))
  
  plot_cw<-ggplot(cw_idx,aes(x=fs_num,y=cw,color=fs_mth,shape=fs_mth,linetype=fs_mth))+
    geom_line()+geom_point(size=3)+
    scale_color_discrete("feature ensemble")+
    scale_shape_manual("feature ensemble",values=linetype_val) +
    scale_linetype_manual("feature ensemble",values=linetype_val) +
    scale_x_continuous(labels=x_scale,breaks=x_scale)+
    scale_y_continuous(limits=c(y_min,y_max),breaks=seq(y_min,y_max,by=0.1))+
    theme_bw()+
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          title =element_text(size=8, face='bold'),
          axis.line = element_line(colour = "black"),
          axis.text.x=element_text(angle=90,hjust=1))+
    labs(x="number of features (k*1.5)",y="CW Index",
         title=paste0(model_type[i],":Weighted Consistency Index"))
  plot_stability[[paste0(model_type[i],"_cw")]]<-plot_cw
  
  plot_cwr<-ggplot(cw_idx,aes(x=fs_num,y=cw_rel,color=fs_mth,shape=fs_mth, linetype=fs_mth))+
    geom_line()+geom_point(size=3)+
    scale_color_discrete("feature ensemble")+
    scale_shape_manual("feature ensemble",values=linetype_val) +
    scale_linetype_manual("feature ensemble",values=linetype_val) +
    scale_x_continuous(labels=x_scale,breaks=x_scale)+
    scale_y_continuous(limits=c(y_min,y_max),breaks=seq(y_min,y_max,by=0.1))+
    theme_bw()+
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          title =element_text(size=8, face='bold'),
          axis.line = element_line(colour = "black"),
          axis.text.x=element_text(angle=90,hjust=1))+
    labs(x="number of features (k*1.5)",y="Relative CW Index",
         title=paste0(model_type[i],":Relative Weighted Consistency Index"))
  plot_stability[[paste0(model_type[i],"_cwr")]]<-plot_cwr
  
  plot_ki<-ggplot(gk_kch_idx,aes(x=fs_num,y=ki,color=fs_mth,shape=fs_mth,linetype=fs_mth))+
    geom_line()+geom_point(size=3)+
    scale_color_discrete("feature ensemble")+
    scale_shape_manual("feature ensemble",values=linetype_val) +
    scale_linetype_manual("feature ensemble",values=linetype_val) +
    scale_x_continuous(labels=x_scale,breaks=x_scale)+
    scale_y_continuous(limits=c(y_min,y_max),breaks=seq(y_min,y_max,by=0.1))+
    theme_bw()+
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          title =element_text(size=8, face='bold'),
          axis.line = element_line(colour = "black"),
          axis.text.x=element_text(angle=90,hjust=1))+
    labs(x="number of features (k*1.5)",y="Kuncheva Index",
         title=paste0(model_type[i],":Kuncheva Index"))
  plot_stability[[paste0(model_type[i],"_ki")]]<-plot_ki
  
  plot_gk<-ggplot(gk_kch_idx,aes(x=fs_num,y=gk_idx,color=fs_mth,shape=fs_mth, linetype=fs_mth))+
    geom_line()+geom_point(size=3)+
    scale_color_discrete("feature ensemble")+
    scale_shape_manual("feature ensemble",values=linetype_val) +
    scale_linetype_manual("feature ensemble",values=linetype_val) +
    scale_x_continuous(labels=x_scale,breaks=x_scale)+
    scale_y_continuous(limits=c(y_min,y_max),breaks=seq(y_min,y_max,by=0.1))+
    theme_bw()+
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          title =element_text(size=8, face='bold'),
          axis.line = element_line(colour = "black"),
          axis.text.x=element_text(angle=90,hjust=1))+
    labs(x="number of features (k*1.5)",y="GK Index",
         title=paste0(model_type[i],":General Kalousis Index"))
  plot_stability[[paste0(model_type[i],"_gk")]]<-plot_gk
}

# plot all
grid_arrange_shared_legend(plot_stability,
                           ncol=4,nrow=4)
