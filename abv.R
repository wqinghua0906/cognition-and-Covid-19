rm(list = ls())
setwd("E:/他人数据分析/陈祖齐/20240104/陈祖齐MMSE/统计用数据-整理/20240120/计算各组平均变化率")
library(emmeans)
library(optimx)
library(lme4)
library(lmerTest)
library(tidyverse)
library(lubridate)
library(stringr)
data<-read.csv('mmseresult600.csv',header = T)

data$followupyear<-data$dd_followup/12
group<-read.csv('group.csv',header = T)
data<-merge(data,group,by='RID',all.x = T)
data$birthday<-as.Date(data$birthday)
data$petbldate<-as.Date(data$petbldate)
#计算出生日期与PET阳性时的时间差值，保留zheng位。
data$pet_age<-round(time_length(interval(data$birthday,data$petbldate), "year"))
write.csv(data,"mmseresult600.csv")



data$group<-as.factor(data$group)
data1<-subset(data,data$group==1)
data2<-subset(data,data$group==2)
data3<-subset(data,data$group==3)
data12<-subset(data,data$group==1|data$group==2)
data13<-subset(data,data$group==1|data$group==3)
data23<-subset(data,data$group==2|data$group==3)
names(data)

###############MMSE
###每个组的斜率
##CN
mmse1_0<-lmer(MMSECORE~followupyear+(1+followupyear|RID),data=data1,control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')))
summary(mmse1_0)
confint(mmse1_0,method="Wald")
mmse1_1<-lmer(MMSECORE~followupyear+pet_age+sex+edu+apoe4+(1+followupyear|RID),data=data1,control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')))
summary(mmse1_1)
confint(mmse1_1,method="Wald")


#MCI
mmse2_0<-lmer(MMSECORE~followupyear+(1+followupyear|RID),data=data2,control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')))
summary(mmse2_0)
confint(mmse2_0,method="Wald")
mmse2_1<-lmer(MMSECORE~followupyear+pet_age+sex+edu+apoe4+(1+followupyear|RID),data=data2,control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')))
summary(mmse2_1)
confint(mmse2_1,method="Wald")

#AD
mmse3_0<-lmer(MMSECORE~followupyear+(1+followupyear|RID),data=data3,control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')))
summary(mmse3_0)
confint(mmse3_0,method="Wald")
mmse3_1<-lmer(MMSECORE~followupyear+pet_age+sex+edu+apoe4+(1+followupyear|RID),data=data3,control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')))
summary(mmse3_1)
confint(mmse3_1,method="Wald")




##两两对比
##CN VS MCI
mmse12_0<-lmer(MMSECORE~followupyear*group+(1|RID),data=data12,control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')))
summary(mmse12_0)
confint(mmse12_0,method="Wald")
mmse12_1<-lmer(MMSECORE~followupyear*group+pet_age+sex+edu+apoe4+(1|RID),data=data12,control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')))
summary(mmse12_1)
confint(mmse12_1,method="Wald")


##CN VS AD
mmse13_0<-lmer(MMSCORE~followupyear*group+(1|RID),data=data13,control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')))
summary(mmse13_0)
confint(mmse13_0,method="Wald")
mmse13_1<-lmer(MMSCORE~followupyear*group+pet_age+sex+edu+apoe4+(1|RID),data=data13,control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')))
summary(mmse13_1)
confint(mmse13_1,method="Wald")

##MCI VS AD
mmse23_0<-lmer(MMSCORE~followupyear*group+(1|RID),data=data23,control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')))
summary(mmse23_0)
confint(mmse23_0,method="Wald")
mmse23_1<-lmer(MMSCORE~followupyear*group+pet_age+sex+edu+apoe4+(1|RID),data=data23,control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')))
summary(mmse23_1)
confint(mmse23_1,method="Wald")






######
#####CDRSB

###每个组的斜率
##CN
CDRSB1_0<-lmer(CDRSB~followupyear+(1+followupyear|RID),data=data1,control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')))
summary(CDRSB1_0)
confint(CDRSB1_0,method="Wald")
CDRSB1_1<-lmer(CDRSB~followupyear+pet_age+sex+edu+apoe4+(1+followupyear|RID),data=data1,control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')))
summary(CDRSB1_1)
confint(CDRSB1_1,method="Wald")


#MCI
CDRSB2_0<-lmer(CDRSB~followupyear+(1+followupyear|RID),data=data2,control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')))
summary(CDRSB2_0)
confint(CDRSB2_0,method="Wald")
CDRSB2_1<-lmer(CDRSB~followupyear+pet_age+sex+edu+apoe4+(1+followupyear|RID),data=data2,control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')))
summary(CDRSB2_1)
confint(CDRSB2_1,method="Wald")

#AD
CDRSB3_0<-lmer(CDRSB~followupyear+(1+followupyear|RID),data=data3,control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')))
summary(CDRSB3_0)
confint(CDRSB3_0,method="Wald")
CDRSB3_1<-lmer(CDRSB~followupyear+pet_age+sex+edu+apoe4+(1+followupyear|RID),data=data3,control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')))
summary(CDRSB3_1)
confint(CDRSB3_1,method="Wald")




##两两对比
##CN VS MCI
CDRSB12_0<-lmer(CDRSB~followupyear*group+(1|RID),data=data12,control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')))
summary(CDRSB12_0)
confint(CDRSB12_0,method="Wald")
CDRSB12_1<-lmer(CDRSB~followupyear*group+pet_age+sex+edu+apoe4+(1|RID),data=data12,control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')))
summary(CDRSB12_1)
confint(CDRSB12_1,method="Wald")


##CN VS AD
CDRSB13_0<-lmer(CDRSB~followupyear*group+(1|RID),data=data13,control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')))
summary(CDRSB13_0)
confint(CDRSB13_0,method="Wald")
CDRSB13_1<-lmer(CDRSB~followupyear*group+pet_age+sex+edu+apoe4+(1|RID),data=data13,control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')))
summary(CDRSB13_1)
confint(CDRSB13_1,method="Wald")

##MCI VS AD
CDRSB23_0<-lmer(CDRSB~followupyear*group+(1|RID),data=data23,control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')))
summary(CDRSB23_0)
confint(CDRSB23_0,method="Wald")
CDRSB23_1<-lmer(CDRSB~followupyear*group+pet_age+sex+edu+apoe4+(1|RID),data=data23,control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')))
summary(CDRSB23_1)
confint(CDRSB23_1,method="Wald")

