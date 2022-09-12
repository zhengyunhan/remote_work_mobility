library(miceadds)
library(lmtest)
library(multiwayvcov)
library(modelsummary)
library(fabricatr)
library(mediation)
library(ivpack)
library(tidyr)
library(lubridate)
library(dplyr)
library(stargazer)
library(plm)
library(lfe)

#data pre-processing
df<-read.csv('./data/data_main.csv')
df$date<-as.Date(df$date)
df$trend<-(year(df$date)-2020)*12+month(df$date)-3
df$upt_pct<-(df$ratio_upt_pre-1)*100
df$upt_pct_bus<-(df$ratio_upt_pre_bus-1)*100
df$upt_pct_rail<-(df$ratio_upt_pre_rail-1)*100
df$vrm_pct<-(df$ratio_vrm_pre-1)*100
df$vrm_pct_bus<-(df$ratio_vrm_pre_bus-1)*100
df$vrm_pct_rail<-(df$ratio_vrm_pre_rail-1)*100
df$workplaces_time<-(df$workplaces_time)*100
df$flex_pct<-(df$flex_pct)*100
df$INDIVIDUAL_gap_z<-(df$INDIVIDUAL_gap_z)*1000

#get data for bus and rail
df_bus<-df[df$val_vrm_bus>0,]
df_rail<-df[(df$val_vrm_rail>0)&(df$val_upt_rail>0)&!is.na((df$val_vrm_rail>0)),]

MSA.Code_rail<-df_rail %>%
  group_by(MSA.Code)%>%
  dplyr::summarise(count=n())%>%
  filter(count==21)%>%
  dplyr::select(MSA.Code)
df_rail<-df_rail%>%
  filter(MSA.Code %in% MSA.Code_rail$MSA.Code)

#First stage IV regression
st1_main<-lm(workplaces_time~flex_pct+reopened+INDIVIDUAL_gap_z+cases_p+dose2_cum+vrm_pct+log(pop_month/1000000)+factor(MSA.Code)+factor(month),data=df)
st1_bus<-lm(workplaces_time~flex_pct+reopened+INDIVIDUAL_gap_z+cases_p+dose2_cum+vrm_pct_bus+log(pop_month/1000000)+factor(MSA.Code)+factor(month),data=df_bus)
st1_rail<-lm(workplaces_time~flex_pct+reopened+INDIVIDUAL_gap_z+cases_p+dose2_cum+vrm_pct_rail+log(pop_month/1000000)+factor(MSA.Code)+factor(month),data=df_rail)

rob.st1_main<-coeftest(st1_main,function(x){cluster.vcov(x, df$MSA.Code)})
rob.st1_bus<-coeftest(st1_bus,function(x){cluster.vcov(x, df_bus$MSA.Code)})
rob.st1_rail<-coeftest(st1_rail,function(x){cluster.vcov(x, df_rail$MSA.Code)})

stargazer(st1_main,st1_bus,st1_rail,se = list(rob.st1_main[,"Std. Error"],rob.st1_bus[,"Std. Error"],rob.st1_rail[,"Std. Error"]))

#First stage F stats
fs_main = lm(workplaces_time~flex_pct+reopened+INDIVIDUAL_gap_z+cases_p+dose2_cum+vrm_pct+log(pop_month/1000000)+factor(MSA.Code)+factor(month),data=df)
fn_main=lm(workplaces_time~reopened+INDIVIDUAL_gap_z+cases_p+dose2_cum+vrm_pct+log(pop_month/1000000)+factor(MSA.Code)+factor(month),data=df)

fs_bus = lm(workplaces_time~flex_pct+reopened+INDIVIDUAL_gap_z+cases_p+dose2_cum+vrm_pct_bus+log(pop_month/1000000)+factor(MSA.Code)+factor(month),data=df_bus)
fn_bus=lm(workplaces_time~reopened+INDIVIDUAL_gap_z+cases_p+dose2_cum+vrm_pct_bus+log(pop_month/1000000)+factor(MSA.Code)+factor(month),data=df_bus)

fs_rail = lm(workplaces_time~flex_pct+reopened+INDIVIDUAL_gap_z+cases_p+dose2_cum+vrm_pct_rail+log(pop_month/1000000)+factor(MSA.Code)+factor(month),data=df_rail)
fn_rail=lm(workplaces_time~reopened+INDIVIDUAL_gap_z+cases_p+dose2_cum+vrm_pct_rail+log(pop_month/1000000)+factor(MSA.Code)+factor(month),data=df_rail)

wald_main<-lmtest::waldtest(fs_main, fn_main, vcov = cluster.vcov(fs_main, df$MSA.Code))
wald_bus<-lmtest::waldtest(fs_bus, fn_bus, vcov = cluster.vcov(fs_bus, df_bus$MSA.Code))
wald_rail<-lmtest::waldtest(fs_rail, fn_rail, vcov = cluster.vcov(fs_rail, df_rail$MSA.Code))

print(wald_main)
print(wald_bus)
print(wald_rail)

#### OLS regression
ols_main<-lm(upt_pct~workplaces_time+reopened+INDIVIDUAL_gap_z+cases_p+dose2_cum+vrm_pct+log(pop_month/1000000)+factor(month)+factor(MSA.Code),data=df)
ols_bus<-lm(upt_pct_bus~workplaces_time+reopened+INDIVIDUAL_gap_z+cases_p+dose2_cum+vrm_pct_bus+log(pop_month/1000000)+factor(month)+factor(MSA.Code),data=df_bus)
ols_rail<-lm(upt_pct_rail~workplaces_time+reopened+INDIVIDUAL_gap_z+cases_p+dose2_cum+vrm_pct_rail+log(pop_month/1000000)+factor(month)+factor(MSA.Code),data=df_rail)

rob.ols_main<-coeftest(ols_main,function(x){cluster.vcov(x, df$MSA.Code)})
rob.ols_bus<-coeftest(ols_bus,function(x){cluster.vcov(x, df_bus$MSA.Code)})
rob.ols_rail<-coeftest(ols_rail,function(x){cluster.vcov(x, df_rail$MSA.Code)})

#### IV regression
iv_main<-ivreg(upt_pct~workplaces_time+reopened+INDIVIDUAL_gap_z+cases_p+dose2_cum+vrm_pct+log(pop_month/1000000)+factor(month)+factor(MSA.Code)|flex_pct+reopened+INDIVIDUAL_gap_z+cases_p+dose2_cum+vrm_pct+log(pop_month/1000000)+factor(month)+factor(MSA.Code),data=df)
iv_bus<-ivreg(upt_pct_bus~workplaces_time+reopened+INDIVIDUAL_gap_z+cases_p+dose2_cum+vrm_pct_bus+log(pop_month/1000000)+factor(month)+factor(MSA.Code)|flex_pct+reopened+INDIVIDUAL_gap_z+cases_p+dose2_cum+vrm_pct_bus+log(pop_month/1000000)+factor(month)+factor(MSA.Code),data=df_bus)
iv_rail<-ivreg(upt_pct_rail~workplaces_time+reopened+INDIVIDUAL_gap_z+cases_p+dose2_cum+vrm_pct_rail+log(pop_month/1000000)+factor(month)+factor(MSA.Code)|flex_pct+reopened+INDIVIDUAL_gap_z+cases_p+dose2_cum+vrm_pct_rail+log(pop_month/1000000)+factor(month)+factor(MSA.Code),data=df_rail)

rob.fit_main<-coeftest(iv_main,function(x){cluster.vcov(x, df$MSA.Code)})
rob.fit_bus<-coeftest(iv_bus,function(x){cluster.vcov(x, df_bus$MSA.Code)})
rob.fit_rail<-coeftest(iv_rail,function(x){cluster.vcov(x, df_rail$MSA.Code)})

stargazer(ols_main,iv_main,ols_bus,iv_bus,ols_rail,iv_rail,se = list(rob.ols_main[,"Std. Error"],rob.fit_main[,"Std. Error"],rob.ols_bus[,"Std. Error"],rob.fit_bus[,"Std. Error"],rob.ols_rail[,"Std. Error"],rob.fit_rail[,"Std. Error"]))


