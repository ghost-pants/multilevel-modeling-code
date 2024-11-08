# Code accompanying Benedict et al., 2024
# Written by Jessie Benedict
# Date reviewed: 2/8/24
# R: Version 4.3.2
# R-Studio: Version 2023.12.1+402
# Code for chronic chemogenetic inactivation behavior data, random intercepts model and plotting for Figure 5. 

#optional clear workspace

rm(list=ls())

#install.packages("blank") for any uninstalled packages, library to call for use in this script
#LIBRARIES

library("emmeans")
library("lmerTest")
library("lme4")
library("ggprism")
library("ggpubr")
library("ggplot2")
library("patchwork")
library("dplyr")
library("magrittr")
library("glue")
library("plyr")
library("tidyverse")
library("reshape2")
library("ggplot2")
library("ggsignif")
library("tibble")
library("ggsci")
library("readr")
library("EnvStats")
library("vtable")
library("rstatix")
library("knitr")
library("lattice")
library("scales")
library("RColorBrewer")
#library("lmerModLmerTest")

##PALETTES
my_color_palette_controls<-c("#EF3E36","#fddbc7","#f06719","#f8b620","#12664F","#A6D3A0","#56B4E9","#114CE8","#558EA6","#2B4570","#7873c0","#ce69be")
control_shapes<-c(24,22,23,25,21,24,22,23,25,21,24,22)

my_color_palette_dreadds<-c("#ce69be","#7873c0","#114CE8","#56B4E9","#A6D3A0","#12664F","#f8b620","#f06719","#fddbc7","#EF3E36")
dreadd_shapes<-c(24,23,22,24,21,25,23,22,23,21)

#READ IN CDMBT_MASTER
#set working directory to wherever Manuscript_Code folder is located
setwd("/Volumes/AIM3-ANALYS")

CDMBT <- data.frame(read.csv2("Manuscript_Code/CDMBT_MASTER_DATA_FRAME.csv", 
                                  sep=","), check.names=FALSE)
CDMBT$GROUP<- factor(CDMBT$GROUP, levels = c("CONTROL","DREADD"), order=FALSE)
CDMBT$LD<-c(rep(1:4,each=22))

#verifying the structure of:
str(CDMBT$`LACTATION.DAY`)
str(CDMBT$`GROUP`)
str(CDMBT$`NEST.SCORE`)
str(CDMBT)

#subsetting data frames
CDMBT_DREADDS<-data.frame(CDMBT[(CDMBT$`GROUP` =="DREADD"),], check.names=FALSE)
#setting factor order
CDMBT_DREADDS$ANIMAL<-factor(CDMBT_DREADDS$ANIMAL,levels=CDMBT_DREADDS$ANIMAL[1:10])
#arranging by animal order
CDMBT_DREADDS_REORDERED<-arrange(CDMBT_DREADDS,CDMBT_DREADDS$ANIMAL)
#arranging new df by ld order
CDMBT_DREADDS_REORDERED<-arrange(CDMBT_DREADDS_REORDERED,CDMBT_DREADDS_REORDERED$LD)
#setting as data frame
CDMBT_DREADDS_REORDERED<-data.frame(CDMBT_DREADDS_REORDERED)
CDMBT_CONTROLS$ANIMAL[1:12]
#controls
CDMBT_CONTROLS<-data.frame(CDMBT[(CDMBT$`GROUP` =="CONTROL"),], check.names=FALSE)
#setting the new order
CDMBT_CONTROLS$ANIMAL<-factor(CDMBT_CONTROLS$ANIMAL,levels=CDMBT_CONTROLS$ANIMAL[1:12])
#reordering by animal order
CDMBT_CONTROLS_REORDERED<-arrange(CDMBT_CONTROLS,CDMBT_CONTROLS$ANIMAL)
#reordering by LD
CDMBT_CONTROLS_REORDERED<-arrange(CDMBT_CONTROLS_REORDERED,CDMBT_CONTROLS_REORDERED$LD)
CDMBT_CONTROLS_REORDERED<-data.frame(CDMBT_CONTROLS_REORDERED)

#ordered according to y axis from bottom to top
DREADD_ANIMALS<-paste(CDMBT_DREADDS_REORDERED$ANIMAL[1:10])
CONTROL_ANIMALS<-paste(CDMBT_CONTROLS_REORDERED$ANIMAL[1:12])
DREADD_GROUP<-paste(CDMBT_DREADDS_REORDERED$GROUP[1:10])
CONTROL_GROUP<-paste(CDMBT_CONTROLS_REORDERED$GROUP[1:12])

ANIMALS<-c(CONTROL_ANIMALS,DREADD_ANIMALS)
GROUPS<-c(CONTROL_GROUP,DREADD_GROUP)

CDMBT_REORDERED<-rbind(CDMBT_CONTROLS_REORDERED,CDMBT_DREADDS_REORDERED)

# define the terms: 
# LMER_GATHER = named object with the results of the model; 
# lmer = the code used to run a linear mixed model; 
# ALL.PUPS.GATHERED = the scored windows out of 30 when the dam had all pups retrieved; 
# outcome variable ~ response variable
# * = interaction between two variables
# LD = lactation day 1-4; 
# GROUP = treatment group, controls and DREADDs; 
# (1| `ANIMAL`) = a random effects clause saying the data may be subject to random variation of intercept at the level of the individual animal)

#GETTING P VALUES USING THE SATTERTHWAITE APPROXIMATION
require(lmerTest)

#ALL PUPS GATHERED
LMER_GATHER = lmerTest::lmer(`ALL.PUPS.GATHERED`~ `LD` * `GROUP` + (1|`ANIMAL`), data=CDMBT_REORDERED)
summary(LMER_GATHER)
confint(LMER_GATHER)

#NEST BUILDING
LMER_NESTS = lmerTest::lmer(`NEST.BUILDING`~ `LD` * `GROUP` + (1|`ANIMAL`), data=CDMBT_REORDERED)
summary(LMER_NESTS)
confint(LMER_NESTS)

#linear models
linear_model_gathers<-lm(`ALL.PUPS.GATHERED`~ `LD` * `GROUP`,data=CDMBT_REORDERED)
linear_model_nests<-lm(`NEST.BUILDING`~ `LD` * `GROUP`,data=CDMBT_REORDERED)

# LOOKING FOR LOWEST AIC
AIC(LMER_GATHER)
AIC(linear_model_gathers)

AIC(LMER_NESTS)
AIC(linear_model_nests)


gatherdf<-data.frame(coef(summary(LMER_GATHER)))
nestdf<-data.frame(coef(summary(LMER_NESTS)))
lmer_cdmbt<-rbind(gatherdf,nestdf)
lmer_cdmbt$TESTID<-c("Gathers","Gathers","Gathers","Gathers","Nests","Nests","Nests","Nests")

write.csv(lmer_cdmbt,file="/Volumes/AIM3-ANALYS/Manuscript_Code/lmer_cdmbt_r_table_output.csv",row.names=FALSE)

#extracting model terms for plotting
#gather
summary(LMER_GATHER)
fixef_gather<-unlist(fixef(LMER_GATHER))

control_intercept_gather<-as.numeric(paste(fixef_gather[1]))
control_slope_gather<-as.numeric(paste(fixef_gather[2]))

contrast_intercept_gather<-as.numeric(fixef_gather[3])
contrast_slope_gather<-as.numeric(fixef_gather[4])

dreadd_intercept_gather<- control_intercept_gather+contrast_intercept_gather
dreadd_slope_gather<-control_slope_gather+contrast_slope_gather

#nests
summary(LMER_NESTS)
fixef_nests<-unlist(fixef(LMER_NESTS))

control_intercept_nests<-as.numeric(paste(fixef_nests[1]))
control_slope_nests<-as.numeric(paste(fixef_nests[2]))

contrast_intercept_nests<-as.numeric(fixef_nests[3])
contrast_slope_nests<-as.numeric(fixef_nests[4])

dreadd_intercept_nests<-control_intercept_nests+contrast_intercept_nests
dreadd_slope_nests<-control_slope_nests+contrast_slope_nests
#pulling out info for model dataframe 

# gather group line coordinates (startpoint is LD1, and intercept is LD0 
# (which is unmeasured) so we add intercept + slope for LD1 startpoint, 
# and endpoint is intercept + slope*4 for LD4)
gather_controls_startpoint<-control_intercept_gather+control_slope_gather
gather_controls_endpoint<-control_intercept_gather+control_slope_gather*4

gather_dreadds_startpoint<-dreadd_intercept_gather+dreadd_slope_gather
gather_dreadds_endpoint<-dreadd_intercept_gather+dreadd_intercept_gather*4 

#nests group line coordinates:
nests_controls_startpoint<-control_intercept_nests+control_slope_nests
nests_controls_endpoint<-control_intercept_nests+control_slope_nests*4

nests_dreadds_startpoint<-dreadd_intercept_nests+dreadd_slope_nests
nests_dreadds_endpoint<-dreadd_intercept_nests+dreadd_slope_nests*4  

## code for each individual animal's random effects
#looking at each animals random intercept
#gather
random_intercepts_raw_gather<-print(coef(LMER_GATHER)$ANIMAL[,1])
#nests
random_intercepts_raw_nests<-print(coef(LMER_NESTS)$ANIMAL[,1])

#making individual animal vector:
random_intercepts_raw_gather<-as.vector(c(random_intercepts_raw_gather))
random_intercepts_raw_nests<-as.vector(c(random_intercepts_raw_nests))
#the raw random intercepts are given as they diverge from the reference group (saline) intercept:
mean(random_intercepts_raw_gather)
mean(random_intercepts_raw_nests)
# subtracting out the control intercept from each random intercept
#first making matrix
raw_intercepts_matrix_gather<-matrix(random_intercepts_raw_gather,22,1)
raw_intercepts_matrix_nests<-matrix(random_intercepts_raw_nests,22,1)
#converting all to mean of 0
meanzerointercepts_gather<-raw_intercepts_matrix_gather[, 1] - mean(random_intercepts_raw_gather)
meanzerointerceptsmatrix_gather<-matrix(meanzerointercepts_gather,22,1)
meanzerointercepts_nests<-raw_intercepts_matrix_nests[, 1] - mean(random_intercepts_raw_nests)
meanzerointerceptsmatrix_nests<-matrix(meanzerointercepts_nests,22,1)

#making dataframe
CDMBT_MODEL_DATA<-data.frame(cbind(ANIMALS,GROUPS))
colnames(CDMBT_MODEL_DATA)<-c("ANIMAL","GROUP")

#verifying group ids and animal ids were transferred correctly
print(CDMBT_MODEL_DATA)
print(CDMBT[1:22,1:2])

#filling in the group terms
CDMBT_MODEL_DATA$GROUP_INTERCEPT_GATHER<- ifelse(CDMBT_MODEL_DATA$GROUP=="CONTROL",paste(control_intercept_gather),paste(dreadd_intercept_gather))
CDMBT_MODEL_DATA$GROUP_INTERCEPT_NESTS<- ifelse(CDMBT_MODEL_DATA$GROUP=="CONTROL",paste(control_intercept_nests),paste(dreadd_intercept_nests))

CDMBT_MODEL_DATA$GROUP_SLOPE_GATHER<- ifelse(CDMBT_MODEL_DATA$GROUP=="CONTROL",paste(control_slope_gather),paste(dreadd_slope_gather))
CDMBT_MODEL_DATA$GROUP_SLOPE_NESTS<- ifelse(CDMBT_MODEL_DATA$GROUP=="CONTROL",paste(control_slope_nests),paste(dreadd_slope_nests))

#adding the correct group intercepts to the mean_zero random effects to get each animals' correct random intercept
group_intercepts_gather<-matrix(as.double(CDMBT_MODEL_DATA$GROUP_INTERCEPT_GATHER,22,1))
group_slopes_gather<-matrix(as.double(CDMBT_MODEL_DATA$GROUP_SLOPE_GATHER,22,1))
individual_final_random_intercepts_gather<-meanzerointerceptsmatrix_gather[,1]+group_intercepts_gather[,1]
group_intercepts_nests<-matrix(as.double(CDMBT_MODEL_DATA$GROUP_INTERCEPT_NESTS,22,1))
group_slopes_nests<-matrix(as.double(CDMBT_MODEL_DATA$GROUP_SLOPE_NESTS,22,1))
individual_final_random_intercepts_nests<-meanzerointerceptsmatrix_nests[,1]+group_intercepts_nests[,1]

#startpoints
group_startpoints_gather<- group_intercepts_gather+group_slopes_gather
individual_startpoints_gather<-individual_final_random_intercepts_gather+group_slopes_gather
individual_endpoints_gather<-individual_final_random_intercepts_gather+group_slopes_gather*4
group_startpoints_nests<-group_intercepts_nests+group_slopes_nests
individual_startpoints_nests<-individual_final_random_intercepts_nests+group_slopes_nests
individual_endpoints_nests<-individual_final_random_intercepts_nests+group_slopes_nests*4

#checking work
print(individual_final_random_intercepts_gather)
print(meanzerointercepts_gather)

print(individual_final_random_intercepts_nests)
print(meanzerointercepts_nests)

#adding the complete random intercepts to the df
CDMBT_MODEL_DATA$INDIVIDUAL_RANDOM_INTERCEPTS_GATHER<-individual_final_random_intercepts_gather
CDMBT_MODEL_DATA$INDIVIDUAL_RANDOM_INTERCEPTS_NESTS<-individual_final_random_intercepts_nests
#adding the complete startpoints to the df
CDMBT_MODEL_DATA$INDIVIDUAL_RANDOM_INTERCEPT_BASED_STARTPOINTS_GATHER<-individual_startpoints_gather
CDMBT_MODEL_DATA$INDIVIDUAL_RANDOM_INTERCEPT_BASED_STARTPOINTS_NESTS<-individual_startpoints_nests

#doubling its length to make x and y individual columns
CDMBT_MODEL_DATA_LONG<-rbind(CDMBT_MODEL_DATA,CDMBT_MODEL_DATA)

#making X column
CDMBT_MODEL_DATA_LONG$X<-c(rep(1,times=22),rep(4,times=22))

#making the numbers numeric to multiply and add them
CDMBT_MODEL_DATA_LONG$INDIVIDUAL_RANDOM_INTERCEPT_BASED_STARTPOINTS_GATHER<-as.numeric(CDMBT_MODEL_DATA_LONG$INDIVIDUAL_RANDOM_INTERCEPT_BASED_STARTPOINTS_GATHER)
CDMBT_MODEL_DATA_LONG$GROUP_SLOPE_GATHER<-as.numeric(CDMBT_MODEL_DATA_LONG$GROUP_SLOPE_GATHER)

CDMBT_MODEL_DATA_LONG$INDIVIDUAL_RANDOM_INTERCEPT_BASED_STARTPOINTS_NESTS<-as.numeric(CDMBT_MODEL_DATA_LONG$INDIVIDUAL_RANDOM_INTERCEPT_BASED_STARTPOINTS_NESTS)
CDMBT_MODEL_DATA_LONG$GROUP_SLOPE_NESTS<-as.numeric(CDMBT_MODEL_DATA_LONG$GROUP_SLOPE_NESTS)

# y-values are, at X=1 equal to the individual random intercepts+1 slope; and then where X=4, equal to the individual final 
# random intercepts + the group slope *4

CDMBT_MODEL_DATA_LONG$Y_GATHER<-c(individual_startpoints_gather[1:22],individual_endpoints_gather)
CDMBT_MODEL_DATA_LONG$Y_NESTS<-c(individual_startpoints_nests[1:22],individual_endpoints_nests)

#rearranging by order the animals appear in gathers plots so that legend is easier to parse
CDMBT_MODEL_DATA_LONG<-arrange(CDMBT_MODEL_DATA_LONG,CDMBT_MODEL_DATA_LONG$INDIVIDUAL_RANDOM_INTERCEPTS_GATHER)
CDMBT_MODEL_DATA_LONG<-arrange(CDMBT_MODEL_DATA_LONG,CDMBT_MODEL_DATA_LONG$X)
CDMBT_MODEL_DATA_CONTROLS<-data.frame(CDMBT_MODEL_DATA_LONG[(CDMBT_MODEL_DATA_LONG$GROUP=="CONTROL"),], check.names=FALSE)
CDMBT_MODEL_DATA_DREADDS<-data.frame(CDMBT_MODEL_DATA_LONG[(CDMBT_MODEL_DATA_LONG$GROUP=="DREADD"),], check.names=FALSE)

#making vectors for labeling the scales in the plots
control_animal_ids<-CDMBT_MODEL_DATA_CONTROLS$ANIMAL[1:12]
dreadd_animal_ids<-CDMBT_MODEL_DATA_DREADDS$ANIMAL[1:10]

#setting order for animals to be plotted
CDMBT_CONTROLS_REORDERED$ANIMAL=factor(CDMBT_CONTROLS_REORDERED$ANIMAL,levels=c(control_animal_ids))
CDMBT_DREADDS_REORDERED$ANIMAL=factor(CDMBT_DREADDS_REORDERED$ANIMAL,levels=c(dreadd_animal_ids))

CDMBT_MODEL_DATA_CONTROLS$ANIMAL=factor(CDMBT_MODEL_DATA_CONTROLS$ANIMAL,levels=c(control_animal_ids))
CDMBT_MODEL_DATA_DREADDS$ANIMAL=factor(CDMBT_MODEL_DATA_DREADDS$ANIMAL,levels=c(dreadd_animal_ids))

#for jitter
pd<-position_dodge(0.3)

set.seed(1)
#Fig 5A_LEFT
#CONTROL: point and lines for actual data
FIG_5A_LEFT<-ggplot(data=CDMBT_CONTROLS_REORDERED,aes(x=LD,y=ALL.PUPS.GATHERED, color=`ANIMAL`, shape=`ANIMAL`, fill=`ANIMAL`))+
  #scales
  scale_colour_manual(name="Animal ID",
                      labels=control_animal_ids,
                      values=my_color_palette_controls)+
  scale_fill_manual(name="Animal ID",
                    labels=control_animal_ids,
                    values=my_color_palette_controls)+
  scale_shape_manual(name="Animal ID",
                     labels=control_animal_ids,
                     values=control_shapes)+
  #axes
  scale_x_continuous(name="", limits=c(0.50,4.6))+
  scale_y_continuous(guide="prism_minor",limits= c(-5,33),breaks = seq(-5,30, by = 5), minor_breaks = seq(-5,30, by = 1), name="All pups gathered")+
  #plot data
  geom_line(linewidth=3,alpha=0.9, position=pd)+
  geom_jitter(aes(color=`ANIMAL`, shape=`ANIMAL`, fill=`ANIMAL`), size=10, alpha=0.5, stroke=5, position=pd)+
  
  coord_cartesian(
    xlim = NULL,
    ylim = NULL,
    expand = FALSE,
    default = FALSE,
    clip = "off")+
  #title
  ggtitle("Control Data")+
  theme_prism(base_family="Helvetica")+
  theme(
    plot.title=element_text
    (size=75,hjust=0.5,margin = margin(t = 20, r = 0, b = 0, l = 0),face="bold"),
    axis.title.y=element_text(size=70, margin=margin(r=15,l=15),hjust=0.5),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size=60,margin=margin(t = 10, r = 0, b = 0, l = 0)),
    axis.text.y=element_text(size=60,margin=margin(t = 0, r = 5, b = 0, l = 0)),
    panel.grid.major = element_blank(), 
    panel.grid.minor= element_blank(), 
    axis.line.y= element_line(linewidth =0, color="black"),
    axis.line.x= element_line(linewidth=0, color="black"), 
    prism.ticks.length=unit(.4,"cm"),
    axis.ticks.length = unit(.6, "cm"),
    axis.ticks.y = element_line(lineend="round",linewidth = 2.5),
    axis.ticks.x= element_line(lineend="round",linewidth = 2.5),
    panel.background=element_rect(fill=FALSE,linetype="solid",colour="black",linewidth=0),
    legend.title=element_text(size=40,face="bold",margin = margin(t = 15, r = 0, b = 0, l = 10)), 
    legend.text=element_text(size=40,face="bold",margin = margin(t = 2, r = 0, b = 15, l = 10)),
    legend.spacing.x = unit(1.0, 'cm'),
    legend.spacing.y = unit(1.0, 'cm'),
    legend.background=element_rect(fill=NA,
                                   size=0, linetype="solid", 
                                   colour ="black"),
    panel.border=element_blank(), 
    plot.margin = margin(t=1.5, l=1, r=0, b=0.5, unit = "cm"))+
  guides(fill=guide_legend(reverse=TRUE), 
         shape = guide_legend(reverse=TRUE),
         color = guide_legend(reverse=TRUE))+
  
  # Y-axis line segment drawn in to have it start at y=0 and stop (yend) at y=4. The x value is 0.61 because default padding expansion for discrete variables is 0.6 units and 5% for continuous variables: 
  annotate("segment", x = 0.5, xend = 0.5, y = 0, yend = 30, colour = "black", linewidth=2.5, lineend="round")+
  annotate("segment", x = 0.5, xend = 0.5, y = -5, yend = 0, colour = "grey", linewidth=2.5, lineend="round")+
  annotate("segment", x = 1, xend = 4.0, y = -5, yend = -5, colour = "black", linewidth=2.5, lineend="round")

ggsave("Manuscript_Code/Plots/FIG_5A_LEFT.png", width=20.5, height=13)

#Fig 5A_RIGHT
#DREADDS
set.seed(2)
FIG_5A_RIGHT <-ggplot(data=na.pass(CDMBT_DREADDS_REORDERED),aes(x=LD,y=ALL.PUPS.GATHERED,color=`ANIMAL`,shape=`ANIMAL`,fill=`ANIMAL`))+
  #x-scale
  scale_colour_manual(name="Animal ID",
                      labels=dreadd_animal_ids,
                      values=c(my_color_palette_dreadds))+
  scale_fill_manual(name="Animal ID",
                    labels=dreadd_animal_ids,
                    values=c(my_color_palette_dreadds))+
  scale_shape_manual(name="Animal ID",
                     labels=dreadd_animal_ids,
                     values=dreadd_shapes)+
  scale_x_continuous(name="", limits=c(0.50,4.6))+
  #y-scale
  scale_y_continuous(guide="prism_minor",limits= c(-5,33),breaks = seq(-5,30, by = 5), minor_breaks = seq(-5,30, by = 1), name="")+
  #points
  ggtitle("DREADD Data")+
  geom_line(aes(color=`ANIMAL`),linewidth=3,alpha=1.0, position=pd)+
  geom_jitter(aes(color=`ANIMAL`,shape=`ANIMAL`),size=10,alpha=0.5,stroke=5, position=pd)+
  coord_cartesian(
    xlim = NULL,
    ylim = NULL,
    expand = FALSE,
    default = FALSE,
    clip = "off")+
  theme_prism(base_family="Helvetica")+
  theme(
    plot.title=element_text
    (size=75,hjust=0.5,margin = margin(t = 20, r = 0, b = 0, l = 0),face="bold"),
    axis.title.y=element_text(size=70, margin=margin(r=15,l=15),hjust=0.5),
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
    legend.title=element_text(size=40,face="bold",margin = margin(t = 15, r = 0, b = 0, l = 10)), 
    legend.text=element_text(size=40,face="bold",margin = margin(t = 2, r = 0, b = 15, l = 10)),
    legend.spacing.x = unit(1.0, 'cm'),
    legend.spacing.y = unit(1.0, 'cm'),
    legend.background=element_rect(fill=NA,
                                   size=0, linetype="solid", 
                                   colour ="black"),
    panel.border=element_blank(), 
    plot.margin = margin(t=1.5, l=1, r=0, b=0.5, unit = "cm"))+
  guides(fill=guide_legend(reverse=TRUE), 
         shape = guide_legend(reverse=TRUE),
         color = guide_legend(reverse=TRUE))+
  
  # Y-axis line segment drawn in to have it start at y=0 and stop (yend) at y=4. The x value is 0.61 because default padding expansion for discrete variables is 0.6 units and 5% for continuous variables: 
  annotate("segment", x = 0.5, xend = 0.5, y = 0, yend = 30, colour = "black", linewidth=2.5, lineend="round")+
  annotate("segment", x = 0.5, xend = 0.5, y = -5, yend = 0, colour = "grey", linewidth=2.5, lineend="round")+
  annotate("segment", x = 1, xend = 4.0, y = -5, yend = -5, colour = "black", linewidth=2.5, lineend="round")

ggsave("Manuscript_Code/Plots/FIG_5A_RIGHT.png", width=20.5, height=13)

#FIG 5B_LEFT
set.seed(3)
FIG_5B_LEFT<-ggplot(data=CDMBT_CONTROLS_REORDERED,aes(x=`LD`,y=`ALL.PUPS.GATHERED`))+
  
  scale_colour_manual(name="Animal ID",
                      labels=c(control_animal_ids),
                      values=my_color_palette_controls)+
  scale_shape_manual(name="Animal ID",
                     labels=c(control_animal_ids),
                     values=c(control_shapes))+  
  scale_fill_manual(name="Animal ID",
                    labels=control_animal_ids,
                    values=my_color_palette_controls)+
  scale_x_continuous(name="", limits=c(0.50,4.6))+
  scale_y_continuous(guide="prism_minor",limits= c(-5,33),breaks = seq(-5,30, by = 5), minor_breaks = seq(-5,30, by = 1), name="All pups gathered")+
  ggtitle("Control Model")+
  
  #plotted data
  geom_line(data=CDMBT_MODEL_DATA_CONTROLS, aes(x=X,y=Y_GATHER, group=`ANIMAL` ,color=`ANIMAL`), linewidth=2.5, lineend="round",alpha=0.8,show.legend=TRUE,position=pd)+
  geom_jitter(data=CDMBT_MODEL_DATA_CONTROLS,aes(x= X, y= Y_GATHER, color=`ANIMAL`,shape=`ANIMAL`,fill=`ANIMAL`),size=10,alpha=0.4, stroke=5, position=pd, show.legend=TRUE)+
  #groupline
  geom_segment(data=CDMBT_MODEL_DATA_CONTROLS, aes(x=1, xend=4, y = control_intercept_gather+control_slope_gather, yend = control_intercept_gather+control_slope_gather*4),color="black",lineend="round",linetype="dashed",linewidth=4,show.legend=FALSE)+
  #formatting stuff
  coord_cartesian(
    xlim = NULL,
    ylim = NULL,
    expand = FALSE,
    default = FALSE,
    clip = "off")+
  theme_prism(base_family="Helvetica")+
  theme(
    plot.title=element_text
    (size=75,hjust=0.5,margin = margin(t = 20, r = 0, b = 0, l = 0),face="bold"),
    axis.title.y=element_text(size=70, margin=margin(r=15,l=15),hjust=0.5),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size=60,margin=margin(t = 10, r = 0, b = 0, l = 0)),
    axis.text.y=element_text(size=60,margin=margin(t = 0, r = 5, b = 0, l = 0)),
    panel.grid.major = element_blank(), 
    panel.grid.minor= element_blank(), 
    axis.line.y= element_line(linewidth =0, color="black"),
    axis.line.x= element_line(linewidth=0, color="black"), 
    prism.ticks.length=unit(.4,"cm"),
    axis.ticks.length = unit(.6, "cm"),
    axis.ticks.y = element_line(lineend="round",linewidth = 2.5),
    axis.ticks.x= element_line(lineend="round",linewidth = 2.5),
    panel.background=element_rect(fill=FALSE,linetype="solid",colour="black",linewidth=0),
    legend.title=element_text(size=40,face="bold",margin = margin(t = 15, r = 0, b = 0, l = 10)), 
    legend.text=element_text(size=40,face="bold",margin = margin(t = 2, r = 0, b = 15, l = 10)),
    legend.spacing.x = unit(1.0, 'cm'),
    legend.spacing.y = unit(1.0, 'cm'),
    legend.background=element_rect(fill=NA,
                                   size=0, linetype="solid", 
                                   colour ="black"),
    panel.border=element_blank(), 
    plot.margin = margin(t=1.5, l=1, r=0, b=0.5, unit = "cm"))+
  guides(color = guide_legend("Animal ID"), 
         shape = guide_legend("Animal ID"),
         fill=guide_legend("Animal ID"),
         linetype = guide_legend(override.aes = list(linetype="dashed")))+
  guides(fill=guide_legend(reverse=TRUE), 
         shape = guide_legend(reverse=TRUE),
         color = guide_legend(reverse=TRUE))+
  # Y-axis line segment drawn in to have it start at y=0 and stop (yend) at y=4. The x value is 0.61 because default padding expansion for discrete variables is 0.6 units and 5% for continuous variables: 
  annotate("segment", x = 0.5, xend = 0.5, y = -5, yend = 30, colour = "black", linewidth=2.5, lineend="round")+
  annotate("segment", x = 1, xend = 4.0, y = -5, yend = -5, colour = "black", linewidth=2.5, lineend="round")

ggsave("Manuscript_Code/Plots/FIG_5B_LEFT.png",width=20.5, height=13)

#FIG 5B_RIGHT
set.seed(4)
FIG_5B_RIGHT<-ggplot(data=CDMBT_DREADDS_REORDERED,aes(x=`LD`,y=`ALL.PUPS.GATHERED`))+
  scale_colour_manual(name="Animal ID",
                      labels=dreadd_animal_ids,
                      values=my_color_palette_dreadds)+
  scale_shape_manual(name="Animal ID",
                     labels=dreadd_animal_ids,
                     values=dreadd_shapes)+  
  scale_fill_manual(name="Animal ID",
                    labels=dreadd_animal_ids,
                    values=my_color_palette_dreadds)+
  scale_x_continuous(name="", limits=c(0.50,4.6))+
  scale_y_continuous(guide="prism_minor",limits= c(-5,33),breaks = seq(-5,30, by = 5), minor_breaks = seq(-5,30, by = 1), name="")+
  ggtitle("DREADD Model")+
  
  #plotted data
  geom_line(data=CDMBT_MODEL_DATA_DREADDS, aes(x=X,y=Y_GATHER, group=`ANIMAL` ,color=`ANIMAL`), linewidth=2.5, lineend="round",alpha=0.8,show.legend=TRUE,position=pd)+
  geom_jitter(data=CDMBT_MODEL_DATA_DREADDS,aes(x= X, y= Y_GATHER, color=`ANIMAL`,shape=`ANIMAL`,fill=`ANIMAL`),size=10,alpha=0.4, stroke=5, position=pd, show.legend=TRUE)+

  #group line
  geom_segment(data=CDMBT_MODEL_DATA_DREADDS, aes(x=1, xend=4, y = dreadd_intercept_gather+dreadd_slope_gather, yend = dreadd_intercept_gather+dreadd_slope_gather*4),color="black",lineend="round",linetype="dashed",linewidth=4,show.legend=FALSE)+
  #formatting stuff
  coord_cartesian(
    xlim = NULL,
    ylim = NULL,
    expand = FALSE,
    default = FALSE,
    clip = "off")+
  theme_prism(base_family="Helvetica")+
  theme(
    plot.title=element_text(size=75,hjust=0.5,margin = margin(t = 20, r = 0, b = 0, l = 0),face="bold"),
    axis.title.y=element_text(size=70, margin=margin(r=15,l=15),hjust=0.5),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size=60,margin=margin(t = 10, r = 0, b = 0, l = 0)),
    axis.text.y=element_text(size=60,margin=margin(t = 0, r = 5, b = 0, l = 0)),
    panel.grid.major = element_blank(), 
    panel.grid.minor= element_blank(), 
    axis.line.y= element_line(linewidth =0, color="black"),
    axis.line.x= element_line(linewidth=0, color="black"), 
    prism.ticks.length=unit(.4,"cm"),
    axis.ticks.length = unit(.6, "cm"),
    axis.ticks.y = element_line(lineend="round",linewidth = 2.5),
    axis.ticks.x= element_line(lineend="round",linewidth = 2.5),
    panel.background=element_rect(fill=FALSE,linetype="solid",colour="black",linewidth=0),
    legend.title=element_text(size=40,face="bold",margin = margin(t = 15, r = 0, b = 0, l = 10)), 
    legend.text=element_text(size=40,face="bold",margin = margin(t = 2, r = 0, b = 15, l = 10)),
    legend.spacing.x = unit(1.0, 'cm'),
    legend.spacing.y = unit(1.0, 'cm'),
    legend.background=element_rect(fill=NA,
                                   size=0, linetype="solid", 
                                   colour ="black"),
    panel.border=element_blank(), 
    plot.margin = margin(t=1.5, l=1, r=0, b=0.5, unit = "cm"))+
  guides(color = guide_legend("Animal ID"), 
         shape = guide_legend("Animal ID"),
         fill=guide_legend("Animal ID"),
         linetype = guide_legend(override.aes = list(linetype="dashed")))+
  guides(fill=guide_legend(reverse=TRUE), 
         shape = guide_legend(reverse=TRUE),
         color = guide_legend(reverse=TRUE))+
  # Y-axis line segment drawn in to have it start at y=0 and stop (yend) at y=4. The x value is 0.61 because default padding expansion for discrete variables is 0.6 units and 5% for continuous variables: 
  annotate("segment", x = 0.5, xend = 0.5, y = -5, yend = 30, colour = "black", linewidth=2.5, lineend="round")+
  annotate("segment", x = 1, xend = 4.0, y = -5, yend = -5, colour = "black", linewidth=2.5, lineend="round")

ggsave("Manuscript_Code/Plots/FIG_5B_RIGHT.png",width=20.5, height=13)

#5C_LEFT
#CONTROL: point and lines for actual data
set.seed(5)
FIG_5C_LEFT<-ggplot(data=CDMBT_CONTROLS_REORDERED,aes(x=LD,y=NEST.BUILDING, color=`ANIMAL`,shape=`ANIMAL`,fill=`ANIMAL`))+
  #scales
  scale_colour_manual(name="Animal ID",
                      labels=control_animal_ids,
                      values=my_color_palette_controls)+
  scale_fill_manual(name="Animal ID",
                    labels=control_animal_ids,
                    values=my_color_palette_controls)+
  scale_shape_manual(name="Animal ID",
                     labels=control_animal_ids,
                     values=control_shapes)+
  #axes
  scale_x_continuous(name="", limits=c(0.50,4.6))+
  scale_y_continuous(guide="prism_minor",limits= c(-5,33),breaks = seq(-5,30, by = 5), minor_breaks = seq(-5,30, by = 1), name="Nest building")+
  #plot data
  geom_line(linewidth=3,alpha=0.9, position=pd)+
  geom_jitter(aes(color=`ANIMAL`,shape=`ANIMAL`,fill=`ANIMAL`),size=10,alpha=0.5,stroke=5, position=pd)+
  
  coord_cartesian(
    xlim = NULL,
    ylim = NULL,
    expand = FALSE,
    default = FALSE,
    clip = "off")+
  #title
  ggtitle("Control Data")+
  theme_prism(base_family="Helvetica")+
  theme(
    plot.title=element_text
    (size=75,hjust=0.5,margin = margin(t = 20, r = 0, b = 0, l = 0),face="bold"),
    axis.title.y=element_text(size=70, margin=margin(r=15,l=15),hjust=0.5),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size=60,margin=margin(t = 10, r = 0, b = 0, l = 0)),
    axis.text.y=element_text(size=60,margin=margin(t = 0, r = 5, b = 0, l = 0)),
    panel.grid.major = element_blank(), 
    panel.grid.minor= element_blank(), 
    axis.line.y= element_line(linewidth =0, color="black"),
    axis.line.x= element_line(linewidth=0, color="black"), 
    prism.ticks.length=unit(.4,"cm"),
    axis.ticks.length = unit(.6, "cm"),
    axis.ticks.y = element_line(lineend="round",linewidth = 2.5),
    axis.ticks.x= element_line(lineend="round",linewidth = 2.5),
    panel.background=element_rect(fill=FALSE,linetype="solid",colour="black",linewidth=0),
    legend.title=element_text(size=40,face="bold",margin = margin(t = 15, r = 0, b = 0, l = 10)), 
    legend.text=element_text(size=40,face="bold",margin = margin(t = 2, r = 0, b = 15, l = 10)),
    legend.spacing.x = unit(1.0, 'cm'),
    legend.spacing.y = unit(1.0, 'cm'),
    legend.background=element_rect(fill=NA,
                                   size=0, linetype="solid", 
                                   colour ="black"),
    panel.border=element_blank(), 
    plot.margin = margin(t=1.5, l=1, r=0, b=0.5, unit = "cm"))+
  guides(fill=guide_legend(reverse=TRUE), 
         shape = guide_legend(reverse=TRUE),
         color = guide_legend(reverse=TRUE))+
  
  # Y-axis line segment drawn in to have it start at y=0 and stop (yend) at y=4. The x value is 0.61 because default padding expansion for discrete variables is 0.6 units and 5% for continuous variables: 
  annotate("segment", x = 0.5, xend = 0.5, y = 0, yend = 30, colour = "black", linewidth=2.5, lineend="round")+
  annotate("segment", x = 0.5, xend = 0.5, y = -5, yend = 0, colour = "grey", linewidth=2.5, lineend="round")+
  annotate("segment", x = 1, xend = 4.0, y = -5, yend = -5, colour = "black", linewidth=2.5, lineend="round")

ggsave("Manuscript_Code/Plots/FIG_5C_LEFT.png", width=20.5, height=13)

#Fig 5C_RIGHT
#DREADDS
set.seed(6)
FIG_5C_RIGHT<-ggplot(data=na.pass(CDMBT_DREADDS_REORDERED),aes(x=LD,y=NEST.BUILDING,color=`ANIMAL`,shape=`ANIMAL`,fill=`ANIMAL`))+
  #x-scale
  scale_colour_manual(name="Animal ID",
                      labels=dreadd_animal_ids,
                      values=my_color_palette_dreadds)+
  scale_fill_manual(name="Animal ID",
                    labels=dreadd_animal_ids,
                    values=my_color_palette_dreadds)+
  scale_shape_manual(name="Animal ID",
                     labels=dreadd_animal_ids,
                     values=dreadd_shapes)+
  scale_x_continuous(name="", limits=c(0.50,4.6))+
  #y-scale
  scale_y_continuous(guide="prism_minor",limits= c(-5,33),breaks = seq(-5,30, by = 5), minor_breaks = seq(-5,30, by = 1), name="")+
  #points
  ggtitle("DREADD Data")+
  geom_line(aes(color=`ANIMAL`),linewidth=3,alpha=1.0, position=pd)+
  geom_jitter(aes(color=`ANIMAL`,shape=`ANIMAL`),size=10,alpha=0.5,stroke=5, position=pd)+
  coord_cartesian(
    xlim = NULL,
    ylim = NULL,
    expand = FALSE,
    default = FALSE,
    clip = "off")+
  theme_prism(base_family="Helvetica")+
  theme(
    plot.title=element_text
    (size=75,hjust=0.5,margin = margin(t = 20, r = 0, b = 0, l = 0),face="bold"),
    axis.title.y=element_text(size=70, margin=margin(r=15,l=15),hjust=0.5),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size=60,margin=margin(t = 10, r = 0, b = 0, l = 0)),
    axis.text.y=element_text(size=60,margin=margin(t = 0, r = 5, b = 0, l = 15)),
    panel.grid.major = element_blank(), 
    panel.grid.minor= element_blank(), 
    axis.line.y= element_line(linewidth =0, color="black"),
    axis.line.x= element_line(linewidth=0, color="black"), 
    prism.ticks.length=unit(.4,"cm"),
    axis.ticks.length = unit(.6, "cm"),
    axis.ticks.y = element_line(lineend="round",linewidth = 2.5),
    axis.ticks.x= element_line(lineend="round",linewidth = 2.5),
    panel.background=element_rect(fill=FALSE,linetype="solid",colour="black",linewidth=0),
    legend.title=element_text(size=40,face="bold",margin = margin(t = 15, r = 0, b = 0, l = 10)), 
    legend.text=element_text(size=40,face="bold",margin = margin(t = 2, r = 0, b = 15, l = 10)),
    legend.spacing.x = unit(1.0, 'cm'),
    legend.spacing.y = unit(1.0, 'cm'),
    legend.background=element_rect(fill=NA,
                                   size=0, linetype="solid", 
                                   colour ="black"),
    panel.border=element_blank(), 
    plot.margin = margin(t=1.5, l=1, r=0, b=0.5, unit = "cm"))+
  guides(fill=guide_legend(reverse=TRUE), 
         shape = guide_legend(reverse=TRUE),
         color = guide_legend(reverse=TRUE))+
  
  # Y-axis line segment drawn in to have it start at y=0 and stop (yend) at y=4. The x value is 0.61 because default padding expansion for discrete variables is 0.6 units and 5% for continuous variables: 
  annotate("segment", x = 0.5, xend = 0.5, y = 0, yend = 30, colour = "black", linewidth=2.5, lineend="round")+
  annotate("segment", x = 0.5, xend = 0.5, y = -5, yend = 0, colour = "grey", linewidth=2.5, lineend="round")+
  annotate("segment", x = 1, xend = 4.0, y = -5, yend = -5, colour = "black", linewidth=2.5, lineend="round")

ggsave("Manuscript_Code/Plots/FIG_5C_RIGHT.png", width=20.5, height=13)

#FIG 5D_LEFT
set.seed(7)
FIG_5D_LEFT<-ggplot(data=CDMBT_CONTROLS_REORDERED,aes(x=`LD`,y=`NEST.BUILDING`))+
  
  scale_colour_manual(name="Animal ID",
                      labels=control_animal_ids,
                      values=my_color_palette_controls)+
  scale_shape_manual(name="Animal ID",
                     labels=control_animal_ids,
                     values=control_shapes)+  
  scale_fill_manual(name="Animal ID",
                    labels=control_animal_ids,
                    values=my_color_palette_controls)+
  scale_x_continuous(name="Lactation day", limits=c(0.50,4.6))+
  scale_y_continuous(guide="prism_minor",limits= c(-5,33),breaks = seq(-5,30, by = 5), minor_breaks = seq(-5,30, by = 1), name="Nest building")+
  ggtitle("Control Model")+

  #plotted data
  geom_line(data=CDMBT_MODEL_DATA_CONTROLS, aes(x=X,y=Y_NESTS, group=`ANIMAL` ,color=`ANIMAL`), linewidth=2.5, lineend="round",alpha=0.8,show.legend=TRUE,position=pd)+
  geom_jitter(data=CDMBT_MODEL_DATA_CONTROLS,aes(x= X, y= Y_NESTS, color=`ANIMAL`,shape=`ANIMAL`,fill=`ANIMAL`),size=10,alpha=0.4, stroke=5, position=pd, show.legend=TRUE)+
  
  #group line
  geom_segment(data=CDMBT_MODEL_DATA_CONTROLS, aes(x=1, xend=4, y = control_intercept_nests+control_slope_nests, yend = control_intercept_nests+control_slope_nests*4),color="black",lineend="round",linetype="dashed",linewidth=4,show.legend=FALSE)+
  #formatting stuff
  
  coord_cartesian(
    xlim = NULL,
    ylim = NULL,
    expand = FALSE,
    default = FALSE,
    clip = "off")+
  theme_prism(base_family="Helvetica")+
  theme(
    plot.title=element_text
    (size=75,hjust=0.5,margin = margin(t = 20, r = 0, b = 0, l = 0),face="bold"),
    axis.title.y=element_text(size=70, margin=margin(r=15,l=15),hjust=0.5),
    axis.title.x=element_text(size=70, margin=margin(r=0,t=15,b=15),hjust=0.5),
    axis.text.x=element_text(size=60,margin=margin(t = 10, r = 0, b = 0, l = 0)),
    axis.text.y=element_text(size=60,margin=margin(t = 0, r = 5, b = 0, l = 0)),
    panel.grid.major = element_blank(), 
    panel.grid.minor= element_blank(), 
    axis.line.y= element_line(linewidth = 0, color="black"),
    axis.line.x= element_line(linewidth = 0, color="black"), 
    prism.ticks.length=unit(.4,"cm"),
    axis.ticks.length = unit(.6, "cm"),
    axis.ticks.y = element_line(lineend="round",linewidth = 2.5),
    axis.ticks.x= element_line(lineend="round",linewidth = 2.5),
    panel.background=element_rect(fill=FALSE,linetype="solid",colour="black",linewidth=0),
    legend.title=element_text(size=40,face="bold",margin = margin(t = 15, r = 0, b = 0, l = 10)), 
    legend.text=element_text(size=40,face="bold",margin = margin(t = 2, r = 0, b = 15, l = 10)),
    legend.spacing.x = unit(1.0, 'cm'),
    legend.spacing.y = unit(1.0, 'cm'),
    legend.background=element_rect(fill=NA,
                                   size=0, linetype="solid", 
                                   colour ="black"),
    panel.border=element_blank(), 
    plot.margin = margin(t=1.5, l=1, r=0, b=0.5, unit = "cm"))+
  guides(color = guide_legend("Animal ID"), 
         shape = guide_legend("Animal ID"),
         fill=guide_legend("Animal ID"),
         linetype = guide_legend(override.aes = list(linetype="dashed")))+
  guides(fill=guide_legend(reverse=TRUE), 
         shape = guide_legend(reverse=TRUE),
         color = guide_legend(reverse=TRUE))+
  # Y-axis line segment drawn in to have it start at y=0 and stop (yend) at y=4. The x value is 0.61 because default padding expansion for discrete variables is 0.6 units and 5% for continuous variables: 
  annotate("segment", x = 0.5, xend = 0.5, y = -5, yend = 30, colour = "black", linewidth=2.5, lineend="round")+
  annotate("segment", x = 1, xend = 4.0, y = -5, yend = -5, colour = "black", linewidth=2.5, lineend="round")

ggsave("Manuscript_Code/Plots/FIG_5D_LEFT.png",width=20.5, height=13.5)

#FIG_5D_RIGHT
FIG_5D_RIGHT<-ggplot(data=CDMBT_DREADDS_REORDERED,aes(x=`LD`,y=`NEST.BUILDING`))+
  
  scale_colour_manual(name="Animal ID",
                      labels=dreadd_animal_ids,
                      values=my_color_palette_dreadds)+
  scale_shape_manual(name="Animal ID",
                     labels=dreadd_animal_ids,
                     values=dreadd_shapes)+  
  scale_fill_manual(name="Animal ID",
                    labels=dreadd_animal_ids,
                    values=my_color_palette_dreadds)+
  scale_x_continuous(name="Lactation day", limits=c(0.50,4.6))+
  scale_y_continuous(guide="prism_minor",limits= c(-5,33),breaks = seq(-5,30, by = 5), minor_breaks = seq(-5,30, by = 1), name="")+
  ggtitle("DREADD Model")+

  #plotted data
  geom_line(data=CDMBT_MODEL_DATA_DREADDS, aes(x=X,y=Y_NESTS, group=`ANIMAL` ,color=`ANIMAL`), linewidth=2.5, lineend="round",alpha=0.8,show.legend=TRUE,position=pd)+
  geom_jitter(data=CDMBT_MODEL_DATA_DREADDS,aes(x= X, y= Y_NESTS, color=`ANIMAL`,shape=`ANIMAL`,fill=`ANIMAL`),size=10,alpha=0.4, stroke=5, position=pd, show.legend=TRUE)+
  
  #group line
  geom_segment(data=CDMBT_MODEL_DATA_DREADDS, aes(x=1, xend=4, y = dreadd_intercept_nests+dreadd_slope_nests, yend = dreadd_intercept_nests+dreadd_slope_nests*4),color="black",lineend="round",linetype="dashed",linewidth=4,show.legend=FALSE)+
  #formatting stuff
  
  coord_cartesian(
    xlim = NULL,
    ylim = NULL,
    expand = FALSE,
    default = FALSE,
    clip = "off")+
  theme_prism(base_family="Helvetica")+
  theme(
    plot.title=element_text
    (size=75,hjust=0.5,margin = margin(t = 20, r = 0, b = 0, l = 0),face="bold"),
    axis.title.y=element_text(size=70, margin=margin(r=15,l=15),hjust=0.5),
    axis.title.x=element_text(size=70, margin=margin(r=0,t=15,b=15),hjust=0.5),
    axis.text.x=element_text(size=60,margin=margin(t = 10, r = 0, b = 0, l = 0)),
    axis.text.y=element_text(size=60,margin=margin(t = 0, r = 5, b = 0, l = 0)),
    panel.grid.major = element_blank(), 
    panel.grid.minor= element_blank(), 
    axis.line.y= element_line(linewidth =0, color="black"),
    axis.line.x= element_line(linewidth =0, color="black"), 
    prism.ticks.length=unit(.4,"cm"),
    axis.ticks.length = unit(.6, "cm"),
    axis.ticks.y = element_line(lineend="round",linewidth = 2.5),
    axis.ticks.x= element_line(lineend="round",linewidth = 2.5),
    panel.background=element_rect(fill=FALSE,linetype="solid",colour="black",linewidth=0),
    legend.title=element_text(size=40,face="bold",margin = margin(t = 15, r = 0, b = 0, l = 10)), 
    legend.text=element_text(size=40,face="bold",margin = margin(t = 2, r = 0, b = 15, l = 10)),
    legend.spacing.x = unit(1.0, 'cm'),
    legend.spacing.y = unit(1.0, 'cm'),
    legend.background=element_rect(fill=NA,
                                   size=0, linetype="solid", 
                                   colour ="black"),
    panel.border=element_blank(), 
    plot.margin = margin(t=1.5, l=1, r=0, b=0.5, unit = "cm"))+
  guides(color = guide_legend("Animal ID"), 
         shape = guide_legend("Animal ID"),
         fill=guide_legend("Animal ID"),
         linetype = guide_legend(override.aes = list(linetype="dashed")))+
  guides(fill=guide_legend(reverse=TRUE), 
         shape = guide_legend(reverse=TRUE),
         color = guide_legend(reverse=TRUE))+
  # Y-axis line segment drawn in to have it start at y=0 and stop (yend) at y=4. The x value is 0.61 because default padding expansion for discrete variables is 0.6 units and 5% for continuous variables: 
  annotate("segment", x = 0.5, xend = 0.5, y = -5, yend = 30, colour = "black", linewidth=2.5, lineend="round")+
  annotate("segment", x = 1, xend = 4.0, y = -5, yend = -5, colour = "black", linewidth=2.5, lineend="round")

ggsave("Manuscript_Code/Plots/FIG_5D_RIGHT.png",width=20.5, height=13.5)


#FIG 5 GATHER
FIG_5_GATHER <- ggarrange(FIG_5A_LEFT, FIG_5A_RIGHT, FIG_5B_LEFT, FIG_5B_RIGHT,
                            labels = c("A", "", "B",""),
                            ncol = 2, nrow = 2,font.label = list(size = 80, color = "black", face = "bold", family = "Helvetica"))

#FIG 5 NESTS
FIG_5_NESTS <- ggarrange(FIG_5C_LEFT, FIG_5C_RIGHT, FIG_5D_LEFT, FIG_5D_RIGHT,
                          labels = c("C", "", "D",""),
                          ncol = 2, nrow = 2,font.label = list(size = 80, color = "black", face = "bold", family = "Helvetica"))

#FIG 5 ARRANGE
FIG_5_COMBO <- ggarrange(FIG_5_GATHER,FIG_5_NESTS,
                         
                         ncol = 1, nrow = 2,font.label = list(size = 80, color = "black", face = "bold", family = "Helvetica"))

ggsave("Manuscript_Code/Plots/FIG_5.png",width=41, height=52.5,limitsize=FALSE)
