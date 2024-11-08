# Code accompanying Benedict et al., 2024
# Written by Jessie Benedict
# Date reviewed: 2/8/24
# Spatial histology random intercepts model and plotting for Figure 7

#optional: to clear workspace 
rm(list=ls())

#libraries needed
#install.packages("blank") for any uninstalled packages, library to call for use in this script
library("plyr")
library("tidyverse")
library("reshape2")
library("ggplot2")
library("ggpubr")
library("ggprism")
library("ggsignif")
library("readr")
library("pkgconfig")
library("tibble")
library("stringr")
library("dplyr")
library("knitr")
library("lattice")
library("scales")
library("ggsci")
library("RColorBrewer")
library("lme4")
library("rstatix")
library("cowplot")
library("lmerTest")

#set working directory to wherever Manuscript_Code folder is located
setwd("/Volumes/AIM1")

#read in spatial histology and behavior csv
KAMBT_AP<-read.csv2("Manuscript_Code/KAMBT_SPATIAL-HISTOLOGY_BEHAVIOR-DATA.csv",check.names=FALSE, sep=",")

# excluding Animal X306 since we are testing maternal behavior outcomes and spatial histology, and her maternal
#behavior outcomes are comparable to saline group but her spatial histology may or may not be comparable to 
#saline group since we can't confirm whether her lesion failed or was more anterior than the other females in 
#the kainic acid group, so it would be inappropriate to group her with the salines

# renaming to KAMBT
KAMBT_AP2<-data.frame(KAMBT_AP[!(KAMBT_AP[,1]=="306"),], check.names=FALSE)

# make vector with animal IDs that start with "x" so that R doesn't try to treat them as integers
ID<-as.vector(rep(c("X101","X102", "X104", "X107", "X301", "X302", 
                    "X303", "X304", "X305", "X308"),each=10))

# add vector to dataframe
KAMBT_AP2$ID<-ID
#define data type as factor
KAMBT_AP2$ID=factor(KAMBT_AP2$ID,levels=c("X101","X102", "X104", "X107", "X301", "X302", 
                                  "X303", "X304", "X305", "X308"))
#verify the structure of the column
str(KAMBT_AP2$ID)

#define any other variable types that need to be defined

KAMBT_AP2$GROUP=factor(KAMBT_AP2$GROUP,levels=c("SAL","KA"))
KAMBT_AP2$TREATMENTGROUP<-ifelse(KAMBT_AP2$GROUP=="SAL",paste(0),paste(1)) 
KAMBT_AP2$TREATMENTGROUP<-as.double(KAMBT_AP2$TREATMENTGROUP)
KAMBT_AP2$AVERAGE<-as.double(KAMBT_AP2$AVERAGE)
KAMBT_AP2$`HEMI-1`<-as.double(KAMBT_AP2$`HEMI-1`)
KAMBT_AP2$`HEMI-2`<-as.double(KAMBT_AP2$`HEMI-2`)

print(KAMBT_AP2$GROUP)
print(KAMBT_AP2$TREATMENTGROUP[1:100])

#rename for ease
KAMBT<-KAMBT_AP2

#creating AP axis variable with 10th percentile as y-intercept
KAMBT$AP10Y<-KAMBT$AP-10
KAMBT$AP100Y<-abs(KAMBT$AP-100)

#verifying
print(KAMBT$AP10Y)
print(KAMBT$AP100Y)

#LMER FOR ANT-POST
require(lmerTest)
# define the terms: 
# LMER.. = named object with the results of the model; 
# lmer = the code used to run a linear mixed model; 
# HEMI-1/2 = the hemisphere the cell counts are from; 
# outcome variable ~ response variable
# * = interaction between two variables
# AP10Y/AP100Y = the AP axis scale either starting anterior - posterior or posterior - anterior; 
# GROUP = treatment group, which is the same now as maternal behavior outcomes; 
# (1| `ID`) = a random effects clause saying the data may be subject to random variation of intercept at the level of the individual animal)
# since lmer only allows continuous data as the outcome variable, the statistician consulted recommended running HEMI-1/HEMI-2 as the outcome
# variable and treatment group as the response variable. Just as it doesn't affect the model, it also logically makes no difference whether 
# maternal behavior outcomes predict the cell counts or cell counts predict the maternal behavior outcome to the interpretation of these data

str(KAMBT$`HEMI-1`)

LMER_HEMI1_AP10Y = lmer(`HEMI-1` ~ `AP10Y`*`TREATMENTGROUP` + (1|`ID`), data = KAMBT)
LMER_HEMI1_AP100Y = lmer(`HEMI-1` ~ `AP100Y`*`TREATMENTGROUP` + (1|`ID`), data = KAMBT)
LMER_HEMI2_AP10Y = lmer(`HEMI-2` ~ `AP10Y`*`TREATMENTGROUP` + (1|`ID`), data = KAMBT)
LMER_HEMI2_AP100Y = lmer(`HEMI-2` ~ `AP100Y`*`TREATMENTGROUP` + (1|`ID`), data = KAMBT)

# LMER summary function provides contrast values. In the summary function under fixed effects, the first fixed 
# effect "(Intercept)" corresponds to the intercept of the reference group (saline), and the third fixed effect 
# "GROUPKA" is the contrast value for the intercept of the non-reference group (kainic acid). A contrast value
# is the difference between the two groups' values, so the contrast value given in "GROUPKA" must be added to
# the value of "(Intercept)" to obtain the kainic acid group's intercept. The reference group (saline) slope is
#the 2nd fixed effect value, "AP10Y"/"AP100Y", and finally, "AP10Y:GROUPKA"/"AP100Y:GROUPKA" is the non-reference 
#group contrast value for the slope 

# Note: Running the model with from anterior-posterior and from posterior-anterior is actually the same model 
# (doing the math from AP10Y of adding the slope to the intercept 90 times (to get from 10 to 100) results in the
# exact intercept value given in AP100Y, and the slopes are identical). The only reason we run them both is to 
# collect significance testing on the model's confidence about group differences at the intercept, which differ 
# significantly (aka the model is not sure there's a group difference at the anterior aspect, but is sure there's
# a group difference at the posterior aspect.

#all the data needed for Extended Figure 7-1, adding the contrast outputs to the saline group term values was done in excel. 
summary(LMER_HEMI1_AP10Y)
summary(LMER_HEMI1_AP100Y)
summary(LMER_HEMI2_AP10Y)
summary(LMER_HEMI2_AP100Y)

confint(LMER_HEMI1_AP10Y)
confint(LMER_HEMI1_AP100Y)
confint(LMER_HEMI2_AP10Y)
confint(LMER_HEMI2_AP100Y)

# validating the use of a mixed linear model over a linear model by looking at t-scores and AIC values. 

linear_model_HEMI1_AP10Y<-lm(`HEMI-1`~`AP10Y`*`GROUP`, data=KAMBT)
linear_model_HEMI1_AP100Y<-lm(`HEMI-1`~`AP100Y`*`GROUP`, data=KAMBT)
linear_model_HEMI2_AP10Y<-lm(`HEMI-2`~`AP10Y`*`GROUP`, data=KAMBT)
linear_model_HEMI2_AP100Y<-lm(`HEMI-2`~`AP100Y`*`GROUP`, data=KAMBT)

summary(linear_model_HEMI1_AP10Y)
AIC(linear_model_HEMI1_AP10Y)
AIC(LMER_HEMI1_AP10Y)

summary(linear_model_HEMI1_AP100Y)
AIC(linear_model_HEMI1_AP100Y)
AIC(LMER_HEMI1_AP100Y)

summary(linear_model_HEMI2_AP10Y)
AIC(linear_model_HEMI2_AP10Y)
AIC(LMER_HEMI2_AP10Y)

summary(linear_model_HEMI2_AP100Y)
AIC(linear_model_HEMI2_AP100Y)
AIC(LMER_HEMI2_AP100Y)

# writing values to csv for excel table
dfKAMBTlmer<-data.frame(rbind(coef(summary(LMER_HEMI1_AP10Y)),coef(summary(LMER_HEMI1_AP100Y)),coef(summary(LMER_HEMI2_AP10Y)),coef(summary(LMER_HEMI2_AP100Y))))
Model_ID<-rep(c("HEMI1_AP10Y","HEMI1_AP100Y","HEMI2_AP10Y", "HEMI2_AP100Y"),each=4)
KAMBT_LMER_TABLE<-cbind(dfKAMBTlmer,Model_ID)

#set path for saving csv ouput
write.csv(KAMBT_LMER_TABLE,file="/Volumes/AIM1/Manuscript_Code/KAMBT_LMER_TABLE.csv",row.names=TRUE)

# note: adding the contrast values given by R for kainic-acid group values to the saline-group values
# to get the raw kainic acid group values was conducted in excel as the table was made

# extracting the data for Figure 7 plots
# Anterior intercept; hemi 1
summary(LMER_HEMI1_AP10Y)
coef(summary(LMER_HEMI1_AP10Y))
fixef(LMER_HEMI1_AP10Y)
confint(LMER_HEMI1_AP10Y)
random_effects<-unlist(ranef(LMER_HEMI1_AP10Y))
ranef(LMER_HEMI1_AP10Y)

# Anterior intercept; hemi 2
summary(LMER_HEMI2_AP10Y)
coef(summary(LMER_HEMI2_AP10Y))
fixef(LMER_HEMI2_AP10Y)
confint(LMER_HEMI2_AP10Y)
random_effects<-unlist(ranef(LMER_HEMI2_AP10Y))
ranef(LMER_HEMI2_AP10Y)

#extracting the components for the plotted group lines
# hemi 1
fixef_hemi1<-unlist(fixef(LMER_HEMI1_AP10Y))

saline_intercept_hemi1<-as.numeric(paste(fixef_hemi1[1]))
saline_slope_hemi1<-as.numeric(paste(fixef_hemi1[2]))

contrast_intercept_hemi1<-as.numeric(fixef_hemi1[3])
contrast_slope_hemi1<-as.numeric(fixef_hemi1[4])

kainicacid_intercept_hemi1<- saline_intercept_hemi1+contrast_intercept_hemi1
kainicacid_slope_hemi1<-saline_slope_hemi1+contrast_slope_hemi1

#hemi 2
fixef_hemi2<-unlist(fixef(LMER_HEMI2_AP10Y))

saline_intercept_hemi2<-as.numeric(paste(fixef_hemi2[1]))
saline_slope_hemi2<-as.numeric(paste(fixef_hemi2[2]))

contrast_intercept_hemi2<-as.numeric(fixef_hemi2[3])
contrast_slope_hemi2<-as.numeric(fixef_hemi2[4])

kainicacid_intercept_hemi2<- saline_intercept_hemi2+contrast_intercept_hemi2
kainicacid_slope_hemi2<-saline_slope_hemi2+contrast_slope_hemi2

#pulling out their IDs for the dataframe
KAMBT_ANIMAL_IDS_ONCE=as.character(KAMBT$ID[seq(from=1,to=100,by=10)],levels=c("X101","X102", "X104", "X107", "X301", "X302", 
                                                                               "X303", "X304", "X305", "X308"), order=FALSE)
#pulling out the treatment group IDs by extracting 1 out of 10 rows in the large KAMBT$GROUP vector
KAMBT_TREATMENT_GROUP_ONCE<-KAMBT$GROUP[seq(from = 1, to = 100, by= 10)]
KAMBT_TREATMENT_GROUP_ONCE=as.character(KAMBT_TREATMENT_GROUP_ONCE,levels=c("SAL","KA"),order=FALSE)

#hemi 1 group line coordinates:
saline_group_startpoint1<-saline_intercept_hemi1
saline_group_endpoint1<-saline_intercept_hemi1+saline_slope_hemi1*90
kainicacid_group_startpoint1<-kainicacid_intercept_hemi1
kainicacid_group_endpoint1<-kainicacid_intercept_hemi1+kainicacid_slope_hemi1*90   

#hemi 2 group line coordinates:
saline_group_startpoint2<-saline_intercept_hemi2
saline_group_endpoint2<-saline_intercept_hemi2+saline_slope_hemi2*90
kainicacid_group_startpoint2<-kainicacid_intercept_hemi2
kainicacid_group_endpoint2<-kainicacid_intercept_hemi2+kainicacid_slope_hemi2*90      

## code for each individual animal's random effects
#looking at each animals random intercept

#hemi 1
random_intercepts_raw1<-print(coef(LMER_HEMI1_AP10Y)$ID[,1])
#hemi 2
random_intercepts_raw2<-print(coef(LMER_HEMI2_AP10Y)$ID[,1])

#making individual animal vector:
random_intercepts_raw1<-as.vector(c(random_intercepts_raw1))
random_intercepts_raw2<-as.vector(c(random_intercepts_raw2))
#the raw random intercepts are given as they diverge from the reference group (saline) intercept:
mean(random_intercepts_raw1)
mean(random_intercepts_raw2)
# subtracting out the saline intercept from each random intercept
#first making matrix
raw_intercepts_matrix1<-matrix(random_intercepts_raw1,10,1)
raw_intercepts_matrix2<-matrix(random_intercepts_raw2,10,1)
#converting all to mean of 0
meanzerointercepts1<-raw_intercepts_matrix1[, 1] - mean(random_intercepts_raw1)
meanzerointerceptsmatrix1<-matrix(meanzerointercepts1,10,1)
meanzerointercepts2<-raw_intercepts_matrix2[, 1] - mean(random_intercepts_raw2)
meanzerointerceptsmatrix2<-matrix(meanzerointercepts2,10,1)

meanzerointercepts1
meanzerointercepts2

#making dataframe
KAMBT_MODEL_DATA<-data.frame(cbind(KAMBT_ANIMAL_IDS_ONCE,KAMBT_TREATMENT_GROUP_ONCE))
colnames(KAMBT_MODEL_DATA)<-c("ID","GROUP")

KAMBT_MODEL_DATA$GROUP_INTERCEPT1<- ifelse(KAMBT_MODEL_DATA$GROUP=="SAL",paste(saline_intercept_hemi1),paste(kainicacid_intercept_hemi1))
KAMBT_MODEL_DATA$GROUP_INTERCEPT2<- ifelse(KAMBT_MODEL_DATA$GROUP=="SAL",paste(saline_intercept_hemi2),paste(kainicacid_intercept_hemi2))

KAMBT_MODEL_DATA$GROUP_SLOPE1<- ifelse(KAMBT_MODEL_DATA$GROUP=="SAL",paste(saline_slope_hemi1),paste(kainicacid_slope_hemi1))
KAMBT_MODEL_DATA$GROUP_SLOPE2<- ifelse(KAMBT_MODEL_DATA$GROUP=="SAL",paste(saline_slope_hemi2),paste(kainicacid_slope_hemi2))

#adding the correct group intercepts to the mean_zero random effects to get each animals' correct random intercept
group_intercepts1<-matrix(as.double(KAMBT_MODEL_DATA$GROUP_INTERCEPT1,10,1))
individual_final_random_intercepts1<-meanzerointerceptsmatrix1[,1]+group_intercepts1[,1]
group_intercepts2<-matrix(as.double(KAMBT_MODEL_DATA$GROUP_INTERCEPT2,10,1))
individual_final_random_intercepts2<-meanzerointerceptsmatrix2[,1]+group_intercepts2[,1]

#checking work
print(individual_final_random_intercepts1)
print(meanzerointercepts1)

print(individual_final_random_intercepts2)
print(meanzerointercepts2)

#adding the complete random intercepts to the df
KAMBT_MODEL_DATA$INDIVIDUAL_RANDOM_INTERCEPTS1<-individual_final_random_intercepts1

KAMBT_MODEL_DATA$INDIVIDUAL_RANDOM_INTERCEPTS2<-individual_final_random_intercepts2

#doubling its length to make x and y individual columns
KAMBT_MODEL_DATA_LONG<-rbind(KAMBT_MODEL_DATA,KAMBT_MODEL_DATA)

#making X column
KAMBT_MODEL_DATA_LONG$X<-c(rep(10,times=10),rep(100,times=10))

#making the numbers numeric to multiply and add them
KAMBT_MODEL_DATA_LONG$INDIVIDUAL_RANDOM_INTERCEPTS1<-as.numeric(KAMBT_MODEL_DATA_LONG$INDIVIDUAL_RANDOM_INTERCEPTS1)
KAMBT_MODEL_DATA_LONG$GROUP_SLOPE1<-as.numeric(KAMBT_MODEL_DATA_LONG$GROUP_SLOPE1)

KAMBT_MODEL_DATA_LONG$INDIVIDUAL_RANDOM_INTERCEPTS2<-as.numeric(KAMBT_MODEL_DATA_LONG$INDIVIDUAL_RANDOM_INTERCEPTS2)
KAMBT_MODEL_DATA_LONG$GROUP_SLOPE2<-as.numeric(KAMBT_MODEL_DATA_LONG$GROUP_SLOPE2)

# y-values are, at X=10 equal to the individual random intercepts; and then where X=100, equal to the individual final 
# random intercepts + the group slope *90 (to get to the Y value at X=100)

KAMBT_MODEL_DATA_LONG$Y1<-c(KAMBT_MODEL_DATA_LONG$INDIVIDUAL_RANDOM_INTERCEPTS1[1:10],individual_final_random_intercepts1+(KAMBT_MODEL_DATA_LONG$GROUP_SLOPE1[11:20]*90))
KAMBT_MODEL_DATA_LONG$Y2<-c(KAMBT_MODEL_DATA_LONG$INDIVIDUAL_RANDOM_INTERCEPTS2[1:10],individual_final_random_intercepts2+(KAMBT_MODEL_DATA_LONG$GROUP_SLOPE2[11:20]*90))

#adding interaction term for the legend
KAMBT_MODEL_DATA_LONG$interaction<-interaction(KAMBT_MODEL_DATA_LONG$ID,KAMBT_MODEL_DATA_LONG$GROUP, sep=" : ")
KAMBT$interaction <- interaction(KAMBT$ID, KAMBT$GROUP, sep=" : ")
#subsetting by group for plotting:

KAMBT_MODEL_DATA_SAL<-data.frame(KAMBT_MODEL_DATA_LONG[(KAMBT_MODEL_DATA_LONG$GROUP =="SAL"),],check.names=FALSE)
KAMBT_MODEL_DATA_KA<-data.frame(KAMBT_MODEL_DATA_LONG[(KAMBT_MODEL_DATA_LONG$GROUP =="KA"),],check.names=FALSE)

KAMBT_SALINE<-data.frame(KAMBT[(KAMBT$GROUP=="SAL"),],check.names=FALSE)
KAMBT_KAINICACID<-data.frame(KAMBT[(KAMBT$GROUP=="KA"),],check.names=FALSE)
str(KAMBT_KAINICACID)

#plot objects
pd<-position_dodge(width=8)

#these shouldn't be factors, otherwise it combines the blues into one thing
saline_colors<-c("cornflowerblue","cornflowerblue","cornflowerblue","cornflowerblue","cornflowerblue","cornflowerblue")
saline_shapes<-c(24,22,7,23,25,21)
saline_ids<-as.factor(c("X101","X102","X301","X302","X304","X305","X306"))

kainicacid_shapes<-c(24,7,22,21)
kainicacid_ids<-c("X104","X107","X303","X308")
kainicacid_colors<-c("darkgoldenrod1","darkgoldenrod1","darkgoldenrod1","darkgoldenrod1")

# Code for Figure 7A

FIG_7A_HEMI1_SALINE <- ggplot(KAMBT_SALINE, aes(x=`AP`, y=`HEMI-1`,color=`interaction`,group=`interaction`,fill=`interaction`))+
  #x-scale
scale_shape_manual(name="Animal ID",
                     labels=saline_ids,
                     values=saline_shapes)+
scale_colour_manual(name="Animal ID",
                      labels=saline_ids,
                      values=saline_colors)+
scale_fill_manual(name="Animal ID",
                    labels=saline_ids,
                    values=saline_colors)+
  #scale x
  scale_x_continuous(limits= c(2,110), breaks = seq(10,100, by = 10), name="",expand=c(0, 0.2))+
  #y-scale
  scale_y_continuous(guide="prism_minor", limits= c(-50,1505), breaks = seq(0,1500, by = 300), minor_breaks = seq(0,1500, by = 100), name="Number of neurons surviving",expand=c(0.00, 0.4))+
  #points
  ggtitle("Saline")+
  geom_line(aes(x=`AP`, y=`HEMI-1`),linewidth=3,alpha=0.3, position=pd, show.legend=FALSE)+
  geom_jitter(aes(shape=`ID`),size=10,alpha=0.4, stroke=5, position=pd)+
  coord_cartesian(
    xlim = NULL,
    ylim = NULL,
    expand = FALSE,
    default = FALSE,
    clip = "off")+
  theme_prism(base_family="Helvetica")+
  #all the rest of the theme stuff for axes, titles, legends etc
  theme(
    plot.title=element_text
    (size=75,hjust=0.5,margin = margin(t = 20, r = 0, b = 20, l = 0),face="bold"),
    axis.title.y=element_text(size=70, margin=margin(r=15,l=15),hjust=0.3),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size=60,margin=margin(t = 10, r = 0, b = 0, l = 0)),
    axis.text.y=element_text(size=60,margin=margin(t = 0, r = 5, b = 0, l = 15)),
    panel.grid.major = element_blank(), 
    panel.grid.minor= element_blank(), 
    axis.line.y= element_line(linewidth=0, color="black"),
    axis.line.x= element_line(linewidth=0, color="black"), 
    prism.ticks.length=unit(.4,"cm"),
    axis.ticks.length = unit(.6, "cm"),
    axis.ticks.y = element_line(lineend="round",linewidth = 2.5),
    axis.ticks.x= element_line(lineend="round",linewidth = 2.5),
    panel.background=element_rect(fill=FALSE,linetype="solid",colour="black",linewidth=0),
    legend.title=element_text(size=45,face="bold",margin = margin(t = 15, r =0, b = 0, l = -15)),
    legend.text=element_text(size=45,face="bold",margin = margin(t = 2, r =0, b = 15, l = 10)),
    legend.spacing.x = unit(1.0, 'cm'),
    legend.spacing.y = unit(1.0, 'cm'),
    legend.title.align=0.3,
    legend.background=element_rect(fill=NA,
                                   size=0, linetype="solid", 
                                   colour ="black"),
    panel.border=element_blank(), 
    plot.margin = margin(4,0.5,0,0, "cm"))+
  guides(color = guide_legend("Animal ID"), 
         shape = guide_legend("Animal ID"),
         fill=guide_legend("Animal ID"))+
  
  # Y-axis line segment drawn in to have it start at y=0 and stop (yend) at y=4. The x value is 0.61 because default padding expansion for discrete variables is 0.6 units and 5% for continuous variables: 
  annotate("segment", x = 2, xend = 2, y = 0, yend = 1500, colour = "black", linewidth=2.5, lineend="round")+
  annotate("segment", x = 10, xend = 100, y = -50, yend = -50, colour = "black", linewidth=2.5, lineend="round")+border(color="white",size=10)

ggsave("Manuscript_Code/Plots/FIG_7A_HEMI1_SALINE.png", width=20.5, height=13)


#7A KAINIC ACID
set.seed(2)
#all the main data stuff
FIG_7A_HEMI1_KAINICACID <- ggplot(KAMBT_KAINICACID, aes(x=`AP`, y=`HEMI-1`,color=`interaction`,shape=`interaction`,fill=`interaction`)) +
  #x-scale
  scale_colour_manual(name="Animal ID",
                      labels=kainicacid_ids,
                      values=kainicacid_colors)+
  scale_shape_manual(name="Animal ID",
                     labels=kainicacid_ids,
                     values=c(24,7,22,21))+
  scale_fill_manual(name="Animal ID",
                    labels=kainicacid_ids,
                    values=kainicacid_colors)+
  #scale x
  scale_x_continuous(limits= c(2,110), breaks = seq(10,100, by = 10), name="",expand=c(0, 0.2))+
  #y-scale
  scale_y_continuous(guide="prism_minor", limits= c(-50,1505), breaks = seq(0,1500, by = 300), minor_breaks = seq(0,1500, by = 100), name=" ",expand=c(0.00, 0.4))+
  #points
  ggtitle("Kainic acid")+
  geom_line(linewidth=3,alpha=0.3, position=pd,show.legend=FALSE)+
  geom_jitter(size=10,alpha=0.4, stroke=5, position=pd)+
  coord_cartesian(
    xlim = NULL,
    ylim = NULL,
    expand = FALSE,
    default = FALSE,
    clip = "off")+
  theme_prism(base_family="Helvetica")+
  #all the rest of the theme stuff for axes, titles, legends etc
  theme(
    plot.title=element_text
    (size=75,hjust=0.5,margin = margin(t = 20, r = 0, b = 20, l = 0),face="bold"),
    axis.title.y=element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size=60,margin=margin(t = 10, r = 0, b = 0, l = 0)),
    axis.text.y=element_text(size=60,margin=margin(t = 0, r = 5, b = 0, l = 15)),
    panel.grid.major = element_blank(), 
    panel.grid.minor= element_blank(), 
    axis.line.y= element_line(linewidth=0, color="black"),
    axis.line.x= element_line(linewidth=0, color="black"), 
    prism.ticks.length=unit(.4,"cm"),
    axis.ticks.length = unit(.6, "cm"),
    axis.ticks.y = element_line(lineend="round",linewidth = 2.5),
    axis.ticks.x= element_line(lineend="round",linewidth = 2.5),
    panel.background=element_rect(fill=FALSE,linetype="solid",colour="black",linewidth=0),
    legend.title=element_text(size=45,face="bold",margin = margin(t = 15, r =0, b = 0, l = -15)),
    legend.text=element_text(size=45,face="bold",margin = margin(t = 2, r =0, b = 15, l = 10)),
    legend.spacing.x = unit(1.0, 'cm'),
    legend.spacing.y = unit(1.0, 'cm'),
    legend.title.align=0.3,
    legend.background=element_rect(fill=NA,
                                   size=0, linetype="solid", 
                                   colour ="black"),
    panel.border=element_blank(), 
    plot.margin = margin(4,0.5,0,0, "cm"))+
  guides(color = guide_legend("Animal ID"), 
         shape = guide_legend("Animal ID"),
         fill=guide_legend("Animal ID"))+
  
  # Y-axis line segment drawn in to have it start at y=0 and stop (yend) at y=4. The x value is 0.61 because default padding expansion for discrete variables is 0.6 units and 5% for continuous variables: 
  annotate("segment", x = 2, xend = 2, y = 0, yend = 1500, colour = "black", linewidth=2.5, lineend="round")+
  annotate("segment", x = 10, xend = 100, y = -50, yend = -50, colour = "black", linewidth=2.5, lineend="round")+border(color="white",size=10)

ggsave("Manuscript_Code/Plots/FIG_7A_HEMI1_KAINICACID.png", width=20.5, height=13)


#plot code
set.seed(6)

FIG_7B_HEMI1_SALINE<-ggplot(data=KAMBT_SALINE,aes(x=`AP`,y=`HEMI-1`,
                                            # color=`interaction`,
                                            # group=`interaction`
)) +
  #scales
  scale_shape_manual(name="Animal ID",
                     labels=saline_ids,
                     values=saline_shapes)+
  scale_colour_manual(name="Animal ID",
                      labels=saline_ids,
                      values=saline_colors)+
  scale_fill_manual(name="Animal ID",
                    labels=saline_ids,
                    values=saline_colors)+
  #scale x
  scale_x_continuous(limits= c(2,110), breaks = seq(10,100, by = 10), name="",expand=c(0, 0.2))+
  #y-scale
  scale_y_continuous(guide="prism_minor", limits= c(-50,1500), breaks = seq(0,1500, by = 300), minor_breaks = seq(0,1500, by = 100), name="Number of neurons surviving",expand=c(0.00, 0.4))+
  #points
  ggtitle("Saline")+
  #MODELED DATA: 
  #lines:
  geom_line(data=KAMBT_MODEL_DATA_SAL, aes(x=X,y=Y1,group=`interaction`, color=`interaction`),position=pd,alpha=0.4,linewidth=2.5, lineend="round",show.legend = FALSE)+
  #points:
  geom_jitter(data=KAMBT_MODEL_DATA_SAL,aes(x=X, y=Y1, shape=`interaction`, color=`interaction`, group=`interaction`,fill=`interaction`),size=10,alpha=0.4, stroke=5, position=pd, show.legend=TRUE)+
  #trying group line up here
  geom_segment(aes(x=10, xend=100, y = saline_group_startpoint1, yend = saline_group_endpoint1),color="black",lineend="round",linetype="dashed",linewidth=4,show.legend=FALSE)+
  #formatting stuff
  coord_cartesian(
    xlim = NULL,
    ylim = NULL,
    expand = FALSE,
    default = FALSE,
    clip = "off")+
  theme_prism(base_family="Helvetica")+
  #all the rest of the theme stuff for axes, titles, legends etc
  theme(
    plot.title=element_text
    (size=75,hjust=0.5,margin = margin(t = 20, r = 0, b = 20, l = 0),face="bold"),
    axis.title.y=element_text(size=70, margin=margin(r=15,l=15),hjust=0.3),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size=60,margin=margin(t = 10, r = 0, b = 0, l = 0)),
    axis.text.y=element_text(size=60,margin=margin(t = 0, r = 5, b = 0, l = 15)),
    panel.grid.major = element_blank(), 
    panel.grid.minor= element_blank(), 
    axis.line.y= element_line(linewidth=0, color="black"),
    axis.line.x= element_line(linewidth=0, color="black"), 
    prism.ticks.length=unit(.4,"cm"),
    axis.ticks.length = unit(.6, "cm"),
    axis.ticks.y = element_line(lineend="round",linewidth = 2.5),
    axis.ticks.x= element_line(lineend="round",linewidth = 2.5),
    panel.background=element_rect(fill=FALSE,linetype="solid",colour="black",linewidth=0),
    legend.title=element_text(size=45,face="bold",margin = margin(t = 15, r =0, b = 0, l = -15)),
    legend.text=element_text(size=45,face="bold",margin = margin(t = 2, r =0, b = 15, l = 10)),
    legend.spacing.x = unit(1.0, 'cm'),
    legend.spacing.y = unit(1.0, 'cm'),
    legend.title.align=0.3,
    legend.background=element_rect(fill=NA,
                                   size=0, linetype="solid", 
                                   colour ="black"),
    panel.border=element_blank(), 
    plot.margin = margin(4,0.5,0,0, "cm"))+
  guides(color = guide_legend("Animal ID"), 
         shape = guide_legend("Animal ID"),
         fill=guide_legend("Animal ID"),
         linetype = guide_legend(override.aes = list(linetype="dashed")))+
  # Y-axis line segment drawn in to have it start at y=0 and stop (yend) at y=4. The x value is 0.61 because default padding expansion for discrete variables is 0.6 units and 5% for continuous variables: 
  annotate("segment", x = 2, xend = 2, y = 0, yend = 1500, colour = "black", linewidth=2.5, lineend="round")+
  annotate("segment", x = 10, xend = 100, y = -50, yend = -50, colour = "black", linewidth=2.5, lineend="round")+
  border(color="white",size=10)


ggsave("Manuscript_Code/Plots/FIG_7B_SALINE_HEMI1.png", width=20.5, height=13) 

# plot 7B kainic acid

set.seed(7)
FIG_7B_HEMI1_KAINICACID<-ggplot(data=KAMBT_KAINICACID,aes(x=`AP`,y=`HEMI-1`,
                                              # color=`interaction`,
                                              # group=`interaction`
)) +
  #scales
  scale_shape_manual(name="Animal ID",
                     labels=kainicacid_ids,
                     values=kainicacid_shapes)+
  scale_colour_manual(name="Animal ID",
                      labels=kainicacid_ids,
                      values=kainicacid_colors)+
  scale_fill_manual(name="Animal ID",
                    labels=kainicacid_ids,
                    values=kainicacid_colors)+
  #scale x
  scale_x_continuous(limits= c(2,110), breaks = seq(10,100, by = 10), name="",expand=c(0, 0.2))+
  #y-scale
  scale_y_continuous(guide="prism_minor", limits= c(-50,1500), breaks = seq(0,1500, by = 300), minor_breaks = seq(0,1500, by = 100), name="",expand=c(0.00, 0.4))+
  #points
  ggtitle("Kainic acid")+
  #MODELED DATA: 
  #lines:
  geom_line(data=KAMBT_MODEL_DATA_KA, aes(x=X, y=Y1, group=`interaction`, color=`interaction`),alpha=0.4,linewidth=2.5, lineend="round",show.legend = FALSE)+ 
  #points:
  geom_jitter(data=KAMBT_MODEL_DATA_KA, aes(x=X, y=Y1, shape=`interaction`, color=`interaction`, group=`interaction`,fill=`interaction`),size=10,alpha=0.4, stroke=5, position=pd, show.legend=TRUE)+
  #trying group line up here
  geom_segment(aes(x=10, xend=100, y = kainicacid_group_startpoint1, yend = kainicacid_group_endpoint1),color="black",lineend="round",linetype="dashed",linewidth=4,show.legend=FALSE)+
  #formatting stuff
  coord_cartesian(
    xlim = NULL,
    ylim = NULL,
    expand = FALSE,
    default = FALSE,
    clip = "off")+
  theme_prism(base_family="Helvetica")+
  #all the rest of the theme stuff for axes, titles, legends etc
  theme(
    plot.title=element_text
    (size=75,hjust=0.5,margin = margin(t = 20, r = 0, b = 20, l = 0),face="bold"),
    axis.title.y=element_blank(), 
    axis.title.x=element_blank(),
    axis.text.x=element_text(size=60,margin=margin(t = 10, r = 0, b = 0, l = 0)),
    axis.text.y=element_text(size=60,margin=margin(t = 0, r = 5, b = 0, l = 15)),
    panel.grid.major = element_blank(), 
    panel.grid.minor= element_blank(), 
    axis.line.y= element_line(linewidth=0, color="black"),
    axis.line.x= element_line(linewidth=0, color="black"), 
    prism.ticks.length=unit(.4,"cm"),
    axis.ticks.length = unit(.6, "cm"),
    axis.ticks.y = element_line(lineend="round",linewidth = 2.5),
    axis.ticks.x= element_line(lineend="round",linewidth = 2.5),
    panel.background=element_rect(fill=FALSE,linetype="solid",colour="black",linewidth=0),
    legend.title=element_text(size=45,face="bold",margin = margin(t = 15, r =0, b = 0, l = -15)),
    legend.text=element_text(size=45,face="bold",margin = margin(t = 2, r =0, b = 15, l = 10)),
    legend.spacing.x = unit(1.0, 'cm'),
    legend.spacing.y = unit(1.0, 'cm'),
    legend.title.align=0.3,
    legend.background=element_rect(fill=NA,
                                   size=0, linetype="solid", 
                                   colour ="black"),
    panel.border=element_blank(), 
    plot.margin = margin(4,0.5,0,0, "cm"))+
  guides(color = guide_legend("Animal ID"), 
         shape = guide_legend("Animal ID"),
         fill=guide_legend("Animal ID"),
         #linetype = guide_legend(override.aes = list(linetype="dashed"))
  )+
  # Y-axis line segment drawn in to have it start at y=0 and stop (yend) at y=4. The x value is 0.61 because default padding expansion for discrete variables is 0.6 units and 5% for continuous variables: 
  annotate("segment", x = 2, xend = 2, y = 0, yend = 1500, colour = "black", linewidth=2.5, lineend="round")+
  annotate("segment", x = 10, xend = 100, y = -50, yend = -50, colour = "black", linewidth=2.5, lineend="round")+
  border(color="white",size=10)
#group slope line, just easier to do it this way


ggsave("Manuscript_Code/Plots/FIG_7B_HEMI1_KAINICACID.png", width=20.5, height=13) 

# 7C

set.seed(5)
FIG_7C_HEMI2_SALINE <- ggplot(KAMBT_SALINE, aes(x=`AP`, y=`HEMI-2`,color=`interaction`,group=`interaction`,fill=`interaction`))+
  #x-scale
  scale_shape_manual(name="Animal ID",
                     labels=saline_ids,
                     values=saline_shapes)+
  scale_colour_manual(name="Animal ID",
                      labels=saline_ids,
                      values=saline_colors)+
  scale_fill_manual(name="Animal ID",
                    labels=saline_ids,
                    values=saline_colors)+
  #scale x
  scale_x_continuous(limits= c(2,110), breaks = seq(10,100, by = 10), name="",expand=c(0, 0.2))+
  #y-scale
  scale_y_continuous(guide="prism_minor", limits= c(-50,1505), breaks = seq(0,1500, by = 300), minor_breaks = seq(0,1500, by = 100), name="Number of neurons surviving",expand=c(0.00, 0.4))+
  #points
  ggtitle("Saline")+
  geom_line(aes(x=`AP`, y=`HEMI-2`),linewidth=3,alpha=0.3, position=pd, show.legend=FALSE)+
  geom_jitter(aes(shape=`ID`),size=10,alpha=0.4, stroke=5, position=pd)+
  coord_cartesian(
    xlim = NULL,
    ylim = NULL,
    expand = FALSE,
    default = FALSE,
    clip = "off")+
  theme_prism(base_family="Helvetica")+
  #all the rest of the theme stuff for axes, titles, legends etc
  theme(
    plot.title=element_text
    (size=75,hjust=0.5,margin = margin(t = 20, r = 0, b = 20, l = 0),face="bold"),
    axis.title.y=element_text(size=70, margin=margin(r=15,l=15),hjust=0.3),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size=60,margin=margin(t = 10, r = 0, b = 0, l = 0)),
    axis.text.y=element_text(size=60,margin=margin(t = 0, r = 5, b = 0, l = 15)),
    panel.grid.major = element_blank(), 
    panel.grid.minor= element_blank(), 
    axis.line.y= element_line(linewidth=0, color="black"),
    axis.line.x= element_line(linewidth=0, color="black"), 
    prism.ticks.length=unit(.4,"cm"),
    axis.ticks.length = unit(.6, "cm"),
    axis.ticks.y = element_line(lineend="round",linewidth = 2.5),
    axis.ticks.x= element_line(lineend="round",linewidth = 2.5),
    panel.background=element_rect(fill=FALSE,linetype="solid",colour="black",linewidth=0),
    legend.title=element_text(size=45,face="bold",margin = margin(t = 15, r =0, b = 0, l = -15)),
    legend.text=element_text(size=45,face="bold",margin = margin(t = 2, r =0, b = 15, l = 10)),
    legend.spacing.x = unit(1.0, 'cm'),
    legend.spacing.y = unit(1.0, 'cm'),
    legend.title.align=0.3,
    legend.background=element_rect(fill=NA,
                                   size=0, linetype="solid", 
                                   colour ="black"),
    panel.border=element_blank(), 
    plot.margin = margin(4,0.5,0,0, "cm"))+
  guides(color = guide_legend("Animal ID"), 
         shape = guide_legend("Animal ID"),
         fill=guide_legend("Animal ID"))+
  
  # Y-axis line segment drawn in to have it start at y=0 and stop (yend) at y=4. The x value is 0.61 because default padding expansion for discrete variables is 0.6 units and 5% for continuous variables: 
  annotate("segment", x = 2, xend = 2, y = 0, yend = 1500, colour = "black", linewidth=2.5, lineend="round")+
  annotate("segment", x = 10, xend = 100, y = -50, yend = -50, colour = "black", linewidth=2.5, lineend="round")+border(color="white",size=10)

ggsave("Manuscript_Code/Plots/FIG_7C_HEMI2_SALINE.png", width=20.5, height=13)


#7C KAINIC ACID
set.seed(6)
#all the main data stuff
FIG_7C_HEMI2_KAINICACID <- ggplot(KAMBT_KAINICACID, aes(x=`AP`, y=`HEMI-2`,color=`interaction`,shape=`interaction`,fill=`interaction`)) +
  #x-scale
  scale_colour_manual(name="Animal ID",
                      labels=kainicacid_ids,
                      values=kainicacid_colors)+
  scale_shape_manual(name="Animal ID",
                     labels=kainicacid_ids,
                     values=c(24,7,22,21))+
  scale_fill_manual(name="Animal ID",
                    labels=kainicacid_ids,
                    values=kainicacid_colors)+
  #scale x
  scale_x_continuous(limits= c(2,110), breaks = seq(10,100, by = 10), name="",expand=c(0, 0.2))+
  #y-scale
  scale_y_continuous(guide="prism_minor", limits= c(-50,1505), breaks = seq(0,1500, by = 300), minor_breaks = seq(0,1500, by = 100), name=" ",expand=c(0.00, 0.4))+
  #points
  ggtitle("Kainic acid")+
  geom_line(linewidth=3,alpha=0.3, position=pd,show.legend=FALSE)+
  geom_jitter(size=10,alpha=0.4, stroke=5, position=pd)+
  coord_cartesian(
    xlim = NULL,
    ylim = NULL,
    expand = FALSE,
    default = FALSE,
    clip = "off")+
  theme_prism(base_family="Helvetica")+
  #all the rest of the theme stuff for axes, titles, legends etc
  theme(
    plot.title=element_text
    (size=75,hjust=0.5,margin = margin(t = 20, r = 0, b = 20, l = 0),face="bold"),
    axis.title.y=element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size=60,margin=margin(t = 10, r = 0, b = 0, l = 0)),
    axis.text.y=element_text(size=60,margin=margin(t = 0, r = 5, b = 0, l = 15)),
    panel.grid.major = element_blank(), 
    panel.grid.minor= element_blank(), 
    axis.line.y= element_line(linewidth=0, color="black"),
    axis.line.x= element_line(linewidth=0, color="black"), 
    prism.ticks.length=unit(.4,"cm"),
    axis.ticks.length = unit(.6, "cm"),
    axis.ticks.y = element_line(lineend="round",linewidth = 2.5),
    axis.ticks.x= element_line(lineend="round",linewidth = 2.5),
    panel.background=element_rect(fill=FALSE,linetype="solid",colour="black",linewidth=0),
    legend.title=element_text(size=45,face="bold",margin = margin(t = 15, r =0, b = 0, l = -15)),
    legend.text=element_text(size=45,face="bold",margin = margin(t = 2, r =0, b = 15, l = 10)),
    legend.spacing.x = unit(1.0, 'cm'),
    legend.spacing.y = unit(1.0, 'cm'),
    legend.title.align=0.3,
    legend.background=element_rect(fill=NA,
                                   size=0, linetype="solid", 
                                   colour ="black"),
    panel.border=element_blank(), 
    plot.margin = margin(4,0.5,0,0, "cm"))+
  guides(color = guide_legend("Animal ID"), 
         shape = guide_legend("Animal ID"),
         fill=guide_legend("Animal ID"))+
  
  # Y-axis line segment drawn in to have it start at y=0 and stop (yend) at y=4. The x value is 0.61 because default padding expansion for discrete variables is 0.6 units and 5% for continuous variables: 
  annotate("segment", x = 2, xend = 2, y = 0, yend = 1500, colour = "black", linewidth=2.5, lineend="round")+
  annotate("segment", x = 10, xend = 100, y = -50, yend = -50, colour = "black", linewidth=2.5, lineend="round")+border(color="white",size=10)

ggsave("Manuscript_Code/Plots/FIG_7C_HEMI2_KAINICACID.png", width=20.5, height=13)


#plot code
set.seed(7)

FIG_7D_HEMI2_SALINE<-ggplot(data=KAMBT_SALINE,aes(x=`AP`,y=`HEMI-2`,
                                                  # color=`interaction`,
                                                  # group=`interaction`
)) +
  #scales
  scale_shape_manual(name="Animal ID",
                     labels=saline_ids,
                     values=saline_shapes)+
  scale_colour_manual(name="Animal ID",
                      labels=saline_ids,
                      values=saline_colors)+
  scale_fill_manual(name="Animal ID",
                    labels=saline_ids,
                    values=saline_colors)+
  #scale x
  scale_x_continuous(limits= c(2,110), breaks = seq(10,100, by = 10), name="Anterior-Posterior axis",expand=c(0, 0.2))+
  #y-scale
  scale_y_continuous(guide="prism_minor", limits= c(-50,1500), breaks = seq(0,1500, by = 300), minor_breaks = seq(0,1500, by = 100), name="Number of neurons surviving",expand=c(0.00, 0.4))+
  #points
  ggtitle("Saline")+
  #MODELED DATA: 
  #lines:
  geom_line(data=KAMBT_MODEL_DATA_SAL, aes(x=X,y=Y2,group=`interaction`, color=`interaction`),position=pd,alpha=0.4,linewidth=2.5, lineend="round",show.legend = FALSE)+
  #points:
  geom_jitter(data=KAMBT_MODEL_DATA_SAL,aes(x=X, y=Y2, shape=`interaction`, color=`interaction`, group=`interaction`,fill=`interaction`),size=10,alpha=0.4, stroke=5, position=pd, show.legend=TRUE)+
  #trying group line up here
  geom_segment(aes(x=10, xend=100, y = saline_group_startpoint2, yend = saline_group_endpoint2),color="black",lineend="round",linetype="dashed",linewidth=4,show.legend=FALSE)+
  #formatting stuff
  coord_cartesian(
    xlim = NULL,
    ylim = NULL,
    expand = FALSE,
    default = FALSE,
    clip = "off")+
  theme_prism(base_family="Helvetica")+
  #all the rest of the theme stuff for axes, titles, legends etc
  theme(
    plot.title=element_text
    (size=75,hjust=0.5,margin = margin(t = 20, r = 0, b = 20, l = 0),face="bold"),
    axis.title.y=element_text(size=70, margin=margin(r=15,l=15),hjust=0.5),
    axis.title.x=element_text(size=70, margin=margin(r=0,t=25,b=15),hjust=0.5),
    axis.text.x=element_text(size=60,margin=margin(t = 10, r = 0, b = 0, l = 0)),
    axis.text.y=element_text(size=60,margin=margin(t = 0, r = 5, b = 0, l = 15)),
    panel.grid.major = element_blank(), 
    panel.grid.minor= element_blank(), 
    axis.line.y= element_line(linewidth =0, color="black"),
    axis.line.x= element_line(linewidth =0, color="black"), 
    prism.ticks.length=unit(.4,"cm"),
    axis.ticks.length = unit(.6, "cm"),
    axis.ticks.y = element_line(lineend="round",linewidth = 2.5),
    axis.ticks.x= element_line(lineend="round",linewidth = 2.5),
    panel.background=element_rect(fill=FALSE,linetype="solid",colour="black",linewidth=0),
    legend.title=element_text(size=45,face="bold",margin = margin(t = 15, r =0, b = 0, l = -15)),
    legend.text=element_text(size=45,face="bold",margin = margin(t = 2, r =0, b = 15, l = 10)),
    legend.spacing.x = unit(1.0, 'cm'),
    legend.spacing.y = unit(1.0, 'cm'),
    legend.title.align=0.3,
    legend.background=element_rect(fill=NA,
                                   size=0, linetype="solid", 
                                   colour ="black"),
    panel.border=element_blank(), 
    plot.margin = margin(4,0.5,0,0, "cm"))+
  guides(color = guide_legend("Animal ID"), 
         shape = guide_legend("Animal ID"),
         fill=guide_legend("Animal ID"),
         linetype = guide_legend(override.aes = list(linetype="dashed")))+
  # Y-axis line segment drawn in to have it start at y=0 and stop (yend) at y=4. The x value is 0.61 because default padding expansion for discrete variables is 0.6 units and 5% for continuous variables: 
  annotate("segment", x = 2, xend = 2, y = 0, yend = 1500, colour = "black", linewidth=2.5, lineend="round")+
  annotate("segment", x = 10, xend = 100, y = -50, yend = -50, colour = "black", linewidth=2.5, lineend="round")+
  border(color="white",size=10)



ggsave("Manuscript_Code/Plots/FIG_7D_SALINE_HEMI2.png", width=20.5, height=13.5) 

# plot 7D kainic acid

set.seed(8)
FIG_7D_HEMI2_KAINICACID<-ggplot(data=KAMBT_KAINICACID,aes(x=`AP`,y=`HEMI-2`,
                                                          # color=`interaction`,
                                                          # group=`interaction`
)) +
  #scales
  scale_shape_manual(name="Animal ID",
                     labels=kainicacid_ids,
                     values=kainicacid_shapes)+
  scale_colour_manual(name="Animal ID",
                      labels=kainicacid_ids,
                      values=kainicacid_colors)+
  scale_fill_manual(name="Animal ID",
                    labels=kainicacid_ids,
                    values=kainicacid_colors)+
  #scale x
  scale_x_continuous(limits= c(2,110), breaks = seq(10,100, by = 10), name="Anterior-Posterior axis",expand=c(0, 0.2))+
  #y-scale
  scale_y_continuous(guide="prism_minor", limits= c(-50,1500), breaks = seq(0,1500, by = 300), minor_breaks = seq(0,1500, by = 100), name="",expand=c(0.00, 0.4))+
  #points
  ggtitle("Kainic acid")+
  #MODELED DATA: 
  #lines:
  geom_line(data=KAMBT_MODEL_DATA_KA, aes(x=X, y=Y2, group=`interaction`, color=`interaction`),alpha=0.4,linewidth=2.5, lineend="round",show.legend = FALSE)+ 
  #points:
  geom_jitter(data=KAMBT_MODEL_DATA_KA, aes(x=X, y=Y2, shape=`interaction`, color=`interaction`, group=`interaction`,fill=`interaction`),size=10,alpha=0.4, stroke=5, position=pd, show.legend=TRUE)+
  #trying group line up here
  geom_segment(aes(x=10, xend=100, y = kainicacid_group_startpoint2, yend = kainicacid_group_endpoint2),color="black",lineend="round",linetype="dashed",linewidth=4,show.legend=FALSE)+
  #formatting stuff
  coord_cartesian(
    xlim = NULL,
    ylim = NULL,
    expand = FALSE,
    default = FALSE,
    clip = "off")+
  theme_prism(base_family="Helvetica")+
  #all the rest of the theme stuff for axes, titles, legends etc
  theme(
    plot.title=element_text
    (size=75,hjust=0.5,margin = margin(t = 20, r = 0, b = 20, l = 0),face="bold"),
    axis.title.y=element_blank(), 
    axis.title.x=element_text(size=70, margin=margin(r=0,t=25,b=15),hjust=0.5),
    axis.text.x=element_text(size=60,margin=margin(t = 10, r = 0, b = 0, l = 0)),
    axis.text.y=element_text(size=60,margin=margin(t = 0, r = 5, b = 0, l = 15)),
    panel.grid.major = element_blank(), 
    panel.grid.minor= element_blank(), 
    axis.line.y= element_line(linewidth=0, color="black"),
    axis.line.x= element_line(linewidth=0, color="black"), 
    prism.ticks.length=unit(.4,"cm"),
    axis.ticks.length = unit(.6, "cm"),
    axis.ticks.y = element_line(lineend="round",linewidth = 2.5),
    axis.ticks.x= element_line(lineend="round",linewidth = 2.5),
    panel.background=element_rect(fill=FALSE,linetype="solid",colour="black",linewidth=0),
    legend.title=element_text(size=45,face="bold",margin = margin(t = 15, r =0, b = 0, l = -15)),
    legend.text=element_text(size=45,face="bold",margin = margin(t = 2, r =0, b = 15, l = 10)),
    legend.spacing.x = unit(1.0, 'cm'),
    legend.spacing.y = unit(1.0, 'cm'),
    legend.title.align=0.3,
    legend.background=element_rect(fill=NA,
                                   size=0, linetype="solid", 
                                   colour ="black"),
    panel.border=element_blank(), 
    plot.margin = margin(4,0.5,0,0, "cm"))+
  guides(color = guide_legend("Animal ID"), 
         shape = guide_legend("Animal ID"),
         fill=guide_legend("Animal ID"),
         linetype = guide_legend(override.aes = list(linetype="dashed")))+
  # Y-axis line segment drawn in to have it start at y=0 and stop (yend) at y=4. The x value is 0.61 because default padding expansion for discrete variables is 0.6 units and 5% for continuous variables: 
  annotate("segment", x = 2, xend = 2, y = 0, yend = 1500, colour = "black", linewidth=2.5, lineend="round")+
  annotate("segment", x = 10, xend = 100, y = -50, yend = -50, colour = "black", linewidth=2.5, lineend="round")+
  border(color="white",size=10)

ggsave("Manuscript_Code/Plots/FIG_7D_HEMI2_KAINICACID.png", width=20.5, height=13.5)


##ARRANGING THE FIGURE

# adding left and right panels together for each row of plots

FIG_7A <- ggarrange(FIG_7A_HEMI1_SALINE, FIG_7A_HEMI1_KAINICACID,widths=c(23.75,22.25),heights=c(17.35,17.35),
                                      ncol = 2, nrow = 1,font.label = list(size = 100, color = "black", face = "bold", family = "Helvetica"))

FIG_7B <- ggarrange(FIG_7B_HEMI1_SALINE, FIG_7B_HEMI1_KAINICACID,widths=c(23.75,22.25),heights=c(17.35,17.35),
                    ncol = 2, nrow = 1,font.label = list(size = 100, color = "black", face = "bold", family = "Helvetica"))

FIG_7C <- ggarrange(FIG_7C_HEMI2_SALINE, FIG_7C_HEMI2_KAINICACID,widths=c(23.75,22.25),heights=c(17.35,17.35),
                    ncol = 2, nrow = 1,font.label = list(size = 100, color = "black", face = "bold", family = "Helvetica"))

FIG_7D <- ggarrange(FIG_7D_HEMI2_SALINE, FIG_7D_HEMI2_KAINICACID,widths=c(23.75,22.25),heights=c(18.85,18.85),
                    ncol = 2, nrow = 1,font.label = list(size = 100, color = "black", face = "bold", family = "Helvetica"))

# adding communal titles to each row, ignore warning about converting grob

FIG_7A_LABELED <-ggdraw(aes(family="Helvetica",face="bold"))+
  draw_plot(FIG_7A)+
  draw_label("Hemisphere 1 Data", colour = "black", size = 85,vjust=-8.2, hjust=0.5, fontfamily = "Helvetica",
             fontface = "bold")

FIG_7B_LABELED <-ggdraw(aes(family="Helvetica",face="bold"))+
  draw_plot(FIG_7B)+
  draw_label("Hemisphere 1 Model", colour = "black", size = 85,vjust=-8.2, hjust=0.5, fontfamily = "Helvetica",
             fontface = "bold")

FIG_7C_LABELED <-ggdraw(aes(family="Helvetica",face="bold"))+
  draw_plot(FIG_7C)+
  draw_label("Hemisphere 2 Data", colour = "black", size = 85,vjust=-8.2, hjust=0.5, fontfamily = "Helvetica",
             fontface = "bold")

FIG_7D_LABELED <-ggdraw(aes(family="Helvetica",face="bold"))+
  draw_plot(FIG_7D)+
  draw_label("Hemisphere 2 Model", colour = "black", size = 85,vjust=-9.2, hjust=0.5, fontfamily = "Helvetica",
             fontface = "bold")
# arranging Hemisphere 1 and Hemisphere 2 separately
FIG_7AB <- ggarrange(FIG_7A_LABELED,FIG_7B_LABELED,
                     ncol = 1, nrow = 2,labels = c("A","B"),font.label = list(size = 100, color = "black", face = "bold", family = "Helvetica"))

FIG_7CD <- ggarrange(FIG_7C_LABELED,FIG_7D_LABELED,widths=c(23.75,22.25),heights=c(17.35,18.85),
                     ncol = 1, nrow = 2,labels = c("C","D"),font.label = list(size = 100, color = "black", face = "bold", family = "Helvetica"))

#merging all into one

FIG_7 <- ggarrange(FIG_7AB, FIG_7CD, widths=c(48,48),heights=c(34.7,36.2),
                   ncol = 1, nrow = 2)
#note shift in width in order to accommodate printing margins at 6.93" wide
ggsave("Manuscript_Code/Plots/FIG_7.png", width=55, height=70.9, limitsize=FALSE)

