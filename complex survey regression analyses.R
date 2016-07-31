library(foreign)#import from SAS
library(stats)
library(epiDisplay)
library(survival)
library(Hmisc)
library(survey)
library (plyr) #use to revalue variables
library(reshape)#use to rename variables
library (ggplot2)
library (psych)
library(nlme)
library (survey)

setwd ("C:/Users/laujohns/Dropbox/Vitamin D - NHANES/Analyses/R code") #set working directory

#### Load datasets ####

vitd<-read.xport('data0110.xpt') #observations and variables match those in SAS; N=52,195 participants in cycles 2001-2010
names(vitd)
names(vitd)<-tolower(names(vitd))
names(vitd)#check to make sure lower case

save(vitd,file="vitd.rda")#save as permanent dataset
load("vitd.rda")
load("vitd2_ph_comp.rda")

#### Subset for years 2005-2010 ####

vitd2<-subset(vitd, cycle!=1 & cycle!=2)#only work with data from 2005-2010 (i.e., get rid of 2001-2004 data)
describe.vector(factor(vitd2$cycle), exclude.missing=TRUE)#frequency by cycle to make sure subsetting worked
save(vitd2, file="vitd2.rda")#save dataset

#### Subset for adults >=20 years ####

vitd2_adults<-subset(vitd2, ridageyr>=20)#create dataset of adults; N=17132 which matches NHANES files for total population interviewed!
names(vitd)

######################## VARIABLE CREATION/DATA CLEANING #################################

##### rename bpa variable [not in cycle 1] and check distributions for normality ####
vitd2$bpa<-vitd2$urxbph #rename variable

#check to make sure renaming worked by calculating mean by cycle for both variables
tapply (vitd2$bpa, vitd2$cycle, mean, na.rm=T)
tapply (vitd2$urxbph, vitd2$cycle, mean, na.rm=T)

#investigate normality of data 
hist(log(vitd2$bpa), breaks=50)

#### rename vitamin d variables lbdvid and lbxvid to vitamin_d ####
summ(vitd2$lbxvidms)#look at mean and sd
vitd2$vitamin_d[vitd2$cycle==3]<-vitd2$lbdvidms[vitd2$cycle==3] #rename lbxvid in cycle 3 to vitamin_d
vitd2$vitamin_d[vitd2$cycle==4 | vitd2$cycle==5]<-vitd2$lbxvidms[vitd2$cycle==4 | vitd2$cycle==5]#rename lbdvid in cycles 4 and 5 to vitamin_d
summ(vitd2$vitamin_d)#look at mean and sd of new variable

#make sure number of obs per cycle for vitamin_d match old and newly named variable
length(vitd2$vitamin_d[which(vitd2$cycle==4)]) 
mean(vitd2$vitamin_d[which(vitd2$cycle==4)], na.rm=T)
length(vitd2$lbdvidms[which(vitd2$cycle==4)]) 
mean(vitd2$lbxvidms[which(vitd2$cycle==4)], na.rm=T)

#find means by cycle and make sure number of vitamin d observations matches (#nonmissing)
describeBy(vitd2$vitamin_d, vitd2$cycle)#from psych package get the means/medians by group 
tapply (vitd2$lbdvidms, vitd2$cycle, mean, na.rm=T)
tapply (vitd2$lbxvidms, vitd2$cycle, mean, na.rm=T)

#investigate normality of vitamin_d variable
hist(vitd2$vitamin_d, breaks=30)

#### rename ridageyr to age ####
vitd2$age<-vitd2$ridageyr
length(vitd2$age)
length(vitd2$age[which(vitd2$age==0)])#make sure number of missing "age"="ridageyr"
length(vitd2$ridageyr[which(vitd2$ridageyr==0)])#make sure number of missing "age"="ridageyr"

#create factor variable for age
vitd2$age_cat<-cut2(vitd2$age, c(20, 40, 60), minmax=T) #[a,b) no need to specify lower
table(vitd2$age_cat) #[ 0,20) [20,40) [40,60) [60,85]

#### rename and create education variable ####
vitd2$edu_cat[vitd2$dmdeduc2==1|vitd2$dmdeduc2==2]<-1 #<HS
vitd2$edu_cat[vitd2$dmdeduc2==3]<-2 #HS/GED
vitd2$edu_cat[vitd2$dmdeduc2==4]<-3 #some college
vitd2$edu_cat[vitd2$dmdeduc2==5]<-4 #college or above
vitd2$edu_cat[vitd2$dmdeduc2==7|vitd2$dmdeduc2==9]<-NA
vitd2$edu_cat <- factor(vitd2$edu_cat,
                    levels = c(1,2,3,4),
                    labels = c("<HS", "HS/GED", "Some college", "College or above"))
#check to make sure variable revalue worked
table(vitd2$dmdeduc2)
table(vitd2$edu_cat)

#### rename gender variable ####
vitd2$sex<-vitd2$riagendr
vitd2$sex <- factor(vitd2$sex,
                        levels = c(1,2),
                        labels = c("Male", "Female"))
#check to make sure revalue worked
table(vitd2$riagendr)#1=male, 2=female
table(vitd2$sex)
table(vitd2$riagendr, vitd2$cycle)
table(vitd2$sex, vitd2$cycle)

#rename race/ethnicity variable
vitd2$race<-vitd2$ridreth1
vitd2$race_cat[vitd2$ridreth1==1|vitd2$ridreth1==2]<-1#hispanic
vitd2$race_cat[vitd2$ridreth1==3]<-2 #NHW
vitd2$race_cat[vitd2$ridreth1==4]<-3 #NHB
vitd2$race_cat[vitd2$ridreth1==5]<-4 #other
vitd2$race_cat <- factor(vitd2$race_cat,
                    levels = c(1,2,3,4),
                    labels = c("Hispanic", "NH White", "NH Black", "Other"))
#See if revalue worked
table(vitd2$race_cat)
table(vitd2$ridreth1) #1=MA, #2=Other Hisp, #3=NHW, #4=NHB, #5 other
table(vitd2$race)

#### rename PIR variable ####
vitd2$pir<-vitd2$indfmpir
#see if recode worked
summary(vitd2$indfmpir)
summary(vitd2$pir)
#rename and revalue alcohol variable
vitd2$alc<-vitd2$alq130
vitd2$alc[vitd2$alq130==999 | vitd2$alq130==777]<-NA
#see if revalue worked
table<-table(vitd2$alc, vitd2$cycle)
table
table2<-table(vitd2$alq130, vitd2$cycle)
table2

#### create smoking variables #### 
vitd2$smk[vitd2$smq020==1]<-1 #yes smoked 100 cigs in lifetime and current smoke
vitd2$smk[vitd2$smq020==2]<-2 #no did not smoke
vitd2$smk[vitd2$smq020==7|vitd2$smq020==9]<-NA
vitd2$smk <- factor(vitd2$smk,
                         levels = c(1,2),
                         labels = c("Yes", "No"))
class(vitd2$smk)
summ(vitd2$smk)

vitd2$smk_cat[vitd2$smq040==1 | vitd2$smq040==2]<-1 #current smoker
vitd2$smk_cat[vitd2$smq040==3]<-2 #former
vitd2$smk_cat[vitd2$smk==2]<-3 #never
vitd2$smk_cat[vitd2$smq040==7 |vitd2$smq040==9]<-NA
vitd2$smk_cat <- factor(vitd2$smk_cat,
                         levels = c(1,2,3),
                         labels = c("Current smoker", "Former smoker", "Never smoker"))

#see if recodes worked
table(vitd2$smk, vitd2$cycle)
table(vitd2$smk_cat, vitd2$cycle)
summary(vitd2$smk)
summary(vitd2$smq020)
table2<-table(vitd2$smq020, vitd2$cycle)
table2

#### create bmi variables ####
vitd2$bmi<-vitd2$bmxbmi #continuous
vitd2$bmi_cat<-cut2(vitd2$bmi, c(18.5 , 25, 30)) #<18.5, 18.5-24.99, 25-29.99, 30+
summ(vitd2$bmi_cat)

#see if recode worked
class(vitd2$bmi_cat)
table(vitd2$bmi_cat)
summary(vitd2$bmi)
summary(vitd2$bmxbmi)

#### create time period taking exam ###
vitd2$season<-vitd2$ridexmon #1=nov1-april30; 2=may1-oct31
vitd2$season <- factor(vitd2$season,
                       levels = c(1,2),
                       labels = c("November 1-April 30", "May 1 - October 31"))

#see if recode worked
table<-table(vitd2$season, vitd2$cycle)
table
table2<-table(vitd2$ridexmon, vitd2$cycle)
table2

##### create cotinine variable ### 
vitd2$cot<-vitd2$lbxcot

#see if recode worked
tapply(vitd2$cot, vitd2$cycle, mean, na.rm=T)
tapply(vitd2$lbxcot, vitd2$cycle, mean, na.rm=T)
summary(vitd2$cot)

#######create survey design#########
summ(vitd2$weight)
vitd2$weight[is.na(vitd2$weight)]<-0 #assign a weight of 0 for missing "weight" variables so can create design

summ(vitd2$weight[which(vitd2$weight==0)])

vitd2$weight6yr<-vitd2$weight/3 #weight variable created in SAS

tapply(vitd2$weight, vitd2$cycle, mean, na.rm=T)
tapply(vitd2$weight6yr, vitd2$cycle, mean, na.rm=T)

vitd.dsn<-svydesign(id=~sdmvpsu, strata=~sdmvstra, weights=~weight6yr, data=vitd2, nest=T) #need to create complete dataset

###########################TABLE 1: DEMOGRAPHICS##############################
#Create complete datasets for BPA and phthalates so can figure out N's and %

#vitd2_ph_comp<- na.omit(vitd2[,c("urxucr","urxmbp", "urxbph", "vitamin_d", "pir", "race", "age", "age_cat", "bmi", "bmi_cat","edu_cat",
                                 #"season","sex", "vitdsup", "cot", "cycle", "sdmvpsu","weight8yr","sdmvstra","weight")]) #already >=20 year olds because of edu_cat variable

vitd2_ph_comp<-subset(vitd2, !is.na(vitd2$urxmep) & !is.na(vitd2$urxucr) & !is.na(vitd2$urxbph) & !is.na(vitd2$vitamin_d) & vitd2$age>=20)#subset data, #4724

vitd.ph.comp.dsn<-svydesign(id=~sdmvpsu, strata=~sdmvstra, weights=~weight6yr, data=vitd2_ph_comp, nest=T) #need to create complete dataset

##### weighted N and % ######

### determine missing

table(is.na(vitd2_ph_comp$vitamin_d))

###AGE###
vitd2_ph_comp$age_cat<-cut2(vitd2_ph_comp$age, c(20, 40, 60), minmax=T) #[a,b) no need to specify lower

prop.table(table(vitd2_ph_comp$age_cat))#for unweighted %
table(vitd2_ph_comp$age_cat)#frequencies
svytable(~vitd2_ph_comp$age_cat, Ntotal=100, vitd.ph.comp.dsn)
table(vitd2_ph_comp$cycle)
summ(vitd2_ph_comp$age)

###RACE###
prop.table(table(vitd2_ph_comp$race))
table(vitd2_ph_comp$race)
svytable(~vitd2_ph_comp$race,Ntotal=100, vitd.ph.comp.dsn)
summ(vitd2_ph_comp$race)

###SEX###
prop.table(table(vitd2_ph_comp$sex))
table(vitd2_ph_comp$sex)
svytable(~vitd2_ph_comp$sex,Ntotal=100, vitd.ph.comp.dsn)
summ(vitd2_ph_comp$sex)

###BMI###
prop.table(table(vitd2_ph_comp$bmi_cat))
table(vitd2_ph_comp$bmi_cat)
svytable(~vitd2_ph_comp$bmi_cat,Ntotal=100, vitd.ph.comp.dsn)
summ(vitd2_ph_comp$bmi)

###Vitamin D Supp###
prop.table(table(vitd2_ph_comp$vitdsup))
table(vitd2_ph_comp$vitdsup)
svytable(~vitd2_ph_comp$vitdsup,Ntotal=100, vitd.ph.comp.dsn)
summ(vitd2_ph_comp$vitdsup)

###Smoker###
prop.table(table(vitd2_ph_comp$smk_cat))
table(vitd2_ph_comp$smk_cat)
svytable(~vitd2_ph_comp$smk_cat,Ntotal=100, vitd.ph.comp.dsn)
summ(vitd2_ph_comp$cot)

###Season###
prop.table(table(vitd2_ph_comp$season))
table(vitd2_ph_comp$season)
svytable(~vitd2_ph_comp$season,Ntotal=100, vitd.ph.comp.dsn)
summ(vitd2_ph_comp$season)

###Education###
table(vitd2_ph_comp$edu_cat)
prop.table(table(vitd2_ph_comp$edu_cat))
svytable(~vitd2_ph_comp$edu_cat, Ntotal=100, vitd.ph.comp.dsn)
summ(vitd2_ph_comp$edu_cat)

###median and IQR of vitamin D by demographics###

#Age
svyquantile(~vitd2_ph_comp$vitamin_d,vitd.ph.comp.dsn,c(0.5,0.25,0.75),na.rm=TRUE)
svyquantile(~vitd2_ph_comp$vitamin_d[which(vitd2_ph_comp$age_cat==1)],vitd.ph.comp.dsn[which(vitd2_ph_comp$age_cat==1)],c(0.5,0.25,0.75),na.rm=TRUE)
svyquantile(~vitd2_ph_comp$vitamin_d[which(vitd2_ph_comp$age_cat==2)],vitd.ph.comp.dsn[which(vitd2_ph_comp$age_cat==2)],c(0.5,0.25,0.75),na.rm=TRUE)
svyquantile(~vitd2_ph_comp$vitamin_d[which(vitd2_ph_comp$age_cat==3)],vitd.ph.comp.dsn[which(vitd2_ph_comp$age_cat==3)],c(0.5,0.25,0.75),na.rm=TRUE)

svyby(~vitamin_d, by=~age_cat, design=vitd.ph.comp.dsn, svymean)#use the following for confidence interval: vartype="ci"

summary(svyglm(vitamin_d~as.factor(age_cat), vitd.ph.comp.dsn, na.action=na.omit))

#Race
svyquantile(~vitd2_ph_comp$vitamin_d[which(vitd2_ph_comp$race==1)],vitd.ph.comp.dsn[which(vitd2_ph_comp$race==1)],c(0.5,0.25,0.75),na.rm=TRUE)
svyquantile(~vitd2_ph_comp$vitamin_d[which(vitd2_ph_comp$race==2)],vitd.ph.comp.dsn[which(vitd2_ph_comp$race==2)],c(0.5,0.25,0.75),na.rm=TRUE)
svyquantile(~vitd2_ph_comp$vitamin_d[which(vitd2_ph_comp$race==3)],vitd.ph.comp.dsn[which(vitd2_ph_comp$race==3)],c(0.5,0.25,0.75),na.rm=TRUE)
svyquantile(~vitd2_ph_comp$vitamin_d[which(vitd2_ph_comp$race==4)],vitd.ph.comp.dsn[which(vitd2_ph_comp$race==4)],c(0.5,0.25,0.75),na.rm=TRUE)
svyquantile(~vitd2_ph_comp$vitamin_d[which(vitd2_ph_comp$race==5)],vitd.ph.comp.dsn[which(vitd2_ph_comp$race==5)],c(0.5,0.25,0.75),na.rm=TRUE)

svyby(~vitamin_d, by=~factor(race, c("3", "1", "2", "4", "5")), design=vitd.ph.comp.dsn, svymean)

summary(svyglm(vitamin_d~factor(race,  c("3", "1", "2", "4", "5")), vitd.ph.comp.dsn, na.action=na.omit))

vitd2_ph_comp$race<-factor(vitd2_ph_comp$race, c("3", "1", "2", "4", "5"))
class(vitd2_ph_comp$race)
vitd2_ph_comp$race<-as.factor(vitd2_ph_comp$race)
table(vitd2_ph_comp$race)

#Sex
svyquantile(~vitd2_ph_comp$vitamin_d[which(vitd2_ph_comp$sex==1)],vitd.ph.comp.dsn[which(vitd2_ph_comp$sex==1)],c(0.5,0.25,0.75),na.rm=TRUE)
svyquantile(~vitd2_ph_comp$vitamin_d[which(vitd2_ph_comp$sex==2)],vitd.ph.comp.dsn[which(vitd2_ph_comp$sex==2)],c(0.5,0.25,0.75),na.rm=TRUE)
svyttest(vitamin_d~factor(sex), vitd.ph.comp.dsn, na.action=na.omit)

svyby(~vitamin_d, ~sex, design=vitd.ph.comp.dsn, svymean)#use the following for confidence interval: vartype="ci"

svyttest(vitamin_d~as.factor(sex), vitd.ph.comp.dsn, na.action=na.omit)

#BMI
svyquantile(~vitd2_ph_comp$vitamin_d[which(vitd2_ph_comp$bmi_cat==1)],vitd.ph.comp.dsn[which(vitd2_ph_comp$bmi_cat==1)],c(0.5,0.25,0.75),na.rm=TRUE)
svyquantile(~vitd2_ph_comp$vitamin_d[which(vitd2_ph_comp$bmi_cat==2)],vitd.ph.comp.dsn[which(vitd2_ph_comp$bmi_cat==2)],c(0.5,0.25,0.75),na.rm=TRUE)
svyquantile(~vitd2_ph_comp$vitamin_d[which(vitd2_ph_comp$bmi_cat==3)],vitd.ph.comp.dsn[which(vitd2_ph_comp$bmi_cat==3)],c(0.5,0.25,0.75),na.rm=TRUE)
svyquantile(~vitd2_ph_comp$vitamin_d[which(vitd2_ph_comp$bmi_cat==4)],vitd.ph.comp.dsn[which(vitd2_ph_comp$bmi_cat==4)],c(0.5,0.25,0.75),na.rm=TRUE)
summary(svyglm(vitamin_d~factor(bmi_cat), vitd.ph.comp.dsn, na.action=na.omit))
table(factor(vitd2_ph_comp$bmi_cat, c("2","1","3","4")))

svyby(~vitamin_d, by=~factor(bmi_cat, c("[ 18.5, 25.0)","[ 11.7, 18.5)","[ 25.0, 30.0)","[ 30.0,130.2]")), design=vitd.ph.comp.dsn, svymean)#use the following for confidence interval: vartype="ci"
levels(vitd2_ph_comp$bmi_cat)

summary(svyglm(vitamin_d~factor(bmi_cat, c("2","1","3","4")), vitd.ph.comp.dsn, na.action=na.omit))

#Vitamin D
svyquantile(~vitd2_ph_comp$vitamin_d[which(vitd2_ph_comp$vitdsup==1)],vitd.ph.comp.dsn[which(vitd2_ph_comp$vitdsup==1)],c(0.5,0.25,0.75),na.rm=TRUE)
svyquantile(~vitd2_ph_comp$vitamin_d[which(vitd2_ph_comp$vitdsup==0)],vitd.ph.comp.dsn[which(vitd2_ph_comp$vitdsup==0)],c(0.5,0.25,0.75),na.rm=TRUE)

svyby(~vitamin_d, by=~factor(vitdsup), design=vitd.ph.comp.dsn, svymean)#use the following for confidence interval: vartype="ci"

summary(svyttest(vitamin_d~factor(vitdsup), vitd.ph.comp.dsn, na.action=na.omit))

#Smoker
svyquantile(~vitd2_ph_comp$vitamin_d[which(vitd2_ph_comp$smk_cat==1)],vitd.ph.comp.dsn[which(vitd2_ph_comp$smk_cat==1)],c(0.5,0.25,0.75),na.rm=TRUE)
svyquantile(~vitd2_ph_comp$vitamin_d[which(vitd2_ph_comp$smk_cat==2)],vitd.ph.comp.dsn[which(vitd2_ph_comp$smk_cat==2)],c(0.5,0.25,0.75),na.rm=TRUE)
svyquantile(~vitd2_ph_comp$vitamin_d[which(vitd2_ph_comp$smk_cat==3)],vitd.ph.comp.dsn[which(vitd2_ph_comp$smk_cat==3)],c(0.5,0.25,0.75),na.rm=TRUE)
table(vitd2_ph_comp$smk_cat)

svyby(~vitamin_d, by=~factor(smk_cat, c("3","2","1")), design=vitd.ph.comp.dsn, svymean)#use the following for confidence interval: vartype="ci"

summary(svyglm(vitamin_d~factor(smk_cat, c("3","2","1")), vitd.ph.comp.dsn, na.action=na.omit))

#Season
svyquantile(~vitd2_ph_comp$vitamin_d[which(vitd2_ph_comp$season==1)],vitd.ph.comp.dsn[which(vitd2_ph_comp$season==1)],c(0.5,0.25,0.75),na.rm=TRUE)
svyquantile(~vitd2_ph_comp$vitamin_d[which(vitd2_ph_comp$season==2)],vitd.ph.comp.dsn[which(vitd2_ph_comp$season==2)],c(0.5,0.25,0.75),na.rm=TRUE)

svyby(~vitamin_d, by=~factor(season), design=vitd.ph.comp.dsn, svymean)#use the following for confidence interval: vartype="ci"

svyttest(vitamin_d~factor(season), vitd.ph.comp.dsn)

#Education
svyquantile(~vitd2_ph_comp$vitamin_d[which(vitd2_ph_comp$edu_cat==1)],vitd.ph.comp.dsn[which(vitd2_ph_comp$edu_cat==1)],c(0.5,0.25,0.75),na.rm=TRUE)
svyquantile(~vitd2_ph_comp$vitamin_d[which(vitd2_ph_comp$edu_cat==2)],vitd.ph.comp.dsn[which(vitd2_ph_comp$edu_cat==2)],c(0.5,0.25,0.75),na.rm=TRUE)
svyquantile(~vitd2_ph_comp$vitamin_d[which(vitd2_ph_comp$edu_cat==3)],vitd.ph.comp.dsn[which(vitd2_ph_comp$edu_cat==3)],c(0.5,0.25,0.75),na.rm=TRUE)
svyquantile(~vitd2_ph_comp$vitamin_d[which(vitd2_ph_comp$edu_cat==4)],vitd.ph.comp.dsn[which(vitd2_ph_comp$edu_cat==4)],c(0.5,0.25,0.75),na.rm=TRUE)
table(vitd2_ph_comp$edu_cat)

svyby(~vitamin_d, by=~factor(edu_cat, c("4","3","2","1")), design=vitd.ph.comp.dsn, svymean)#use the following for confidence interval: vartype="ci"

summary(svyglm(vitamin_d~factor(edu_cat, c("4","3","2","1")), vitd.ph.comp.dsn, na.action=na.omit))

####################TABLE 2: GM AND QUANTILES OF PHTHALATES#######################

#MEHP
prop.table(table(vitd2_ph_comp$urdmhplc))
exp(svymean(~log(urxmhp_c),design=vitd.ph.comp.dsn, na.omit=T))
svyquantile(~urxmhp_c,vitd.ph.comp.dsn,c(0.25,0.5,0.75, 0.9, .95, 1),na.rm=TRUE)
summ(vitd2_ph_comp$urxmhp)
summ(vitd2_ph_comp$urxmhp_c[which(vitd2_ph_comp$urdmhplc==0)])
quantile(vitd2_ph_comp$urxmhp[which(vitd2_ph_comp$urdmhplc==1)], c(0.25, 0.5, .75))
tapply(vitd2_ph_comp$urxmhp, vitd2_ph_comp$cycle, mean, na.rm=T)
1.2/sqrt(2)
table(vitd2_ph_comp$cycle)

#MEHHP
prop.table(table(vitd2_ph_comp$urdmhhlc))
exp(svymean(~log(urxmhh_c),design=vitd.ph.comp.dsn, na.omit=T))
svyquantile(~urxmhh_c,vitd.ph.comp.dsn,c(0.25,0.5,0.75, 0.9, .95, 1),na.rm=TRUE)
quantile(vitd2_ph_comp$urxmhh_c[which(vitd2_ph_comp$urdmhhlc==1)], c(0.25, 0.5))

#MEOHP
prop.table(table(vitd2_ph_comp$urdmohlc))
exp(svymean(~log(urxmoh_c),design=vitd.ph.comp.dsn, na.omit=T))
svyquantile(~urxmoh_c,vitd.ph.comp.dsn,c(0.25,0.5,0.75, 0.9, .95, 1),na.rm=TRUE)
quantile(vitd2_ph_comp$urxmoh_c[which(vitd2_ph_comp$urdmohlc==1)], c(0.25, 0.5))

#MECPP
prop.table(table(vitd2_ph_comp$urdecplc))
exp(svymean(~log(urxecp_c),design=vitd.ph.comp.dsn, na.omit=T))
svyquantile(~urxecp_c,vitd.ph.comp.dsn,c(0.25,0.5,0.75, 0.9, .95, 1),na.rm=TRUE)
quantile(vitd2_ph_comp$urxecp_c[which(vitd2_ph_comp$urdecplc==1)], c(0.25, 0.5))

#MBP
prop.table(table(vitd2_ph_comp$urdmbplc))
exp(svymean(~log(urxmbp_c),design=vitd.ph.comp.dsn, na.omit=T))
svyquantile(~urxmbp_c,vitd.ph.comp.dsn,c(0.25,0.5,0.75, 0.9, .95, 1),na.rm=TRUE)
quantile(vitd2_ph_comp$urxmbp_c[which(vitd2_ph_comp$urdmbplc==1)], c(0.25, 0.5, 1))
quantile(vitd2_ph_comp$urxmbp_c, c(0.25, 0.5, 1))

#MiBP
prop.table(table(vitd2_ph_comp$urdmiblc))
exp(svymean(~log(urxmib_c),design=vitd.ph.comp.dsn, na.omit=T))
svyquantile(~urxmib_c,vitd.ph.comp.dsn,c(0.25,0.5,0.75, 0.9, .95, 1),na.rm=TRUE)
quantile(vitd2_ph_comp$urxmib_c[which(vitd2_ph_comp$urdmiblc==1)], c(0.25, 0.5, 1))
quantile(vitd2_ph_comp$urxmib_c, c(0.25, 0.5, 1))

#MEP
prop.table(table(vitd2_ph_comp$urdmeplc))
exp(svymean(~log(urxmep_c),design=vitd.ph.comp.dsn, na.omit=T))
svyquantile(~urxmep_c,vitd.ph.comp.dsn,c(0.25,0.5,0.75, 0.9, .95, 1),na.rm=TRUE)
quantile(vitd2_ph_comp$urxmep_c[which(vitd2_ph_comp$urdmeplc==1)], c(0.25, 0.5, 1))
quantile(vitd2_ph_comp$urxmep_c, c(0.25, 0.5, 1))

#MBzP
prop.table(table(vitd2_ph_comp$urdmzplc))
exp(svymean(~log(urxmzp_c),design=vitd.ph.comp.dsn, na.omit=T))
svyquantile(~urxmzp_c,vitd.ph.comp.dsn,c(0.25,0.5,0.75, 0.9, .95, 1),na.rm=TRUE)
quantile(vitd2_ph_comp$urxmzp_c[which(vitd2_ph_comp$urdmzplc==1)], c(0.25, 0.5, 1))
quantile(vitd2_ph_comp$urxmzp_c, c(0.25, 0.5, 1))

#MCPP
prop.table(table(vitd2_ph_comp$urdmc1lc))
exp(svymean(~log(urxmc1_c),design=vitd.ph.comp.dsn, na.omit=T))
svyquantile(~urxmc1_c,vitd.ph.comp.dsn,c(0.25,0.5,0.75, 0.9, .95, 1),na.rm=TRUE)
quantile(vitd2_ph_comp$urxmc1_c[which(vitd2_ph_comp$urdmc1lc==1)], c(0.25, 0.5, 1))
quantile(vitd2_ph_comp$urxmc1_c, c(0.25, 0.5, 1))

#MiNP
prop.table(table(vitd2_ph_comp$urdmnplc))
exp(svymean(~log(urxmnp_c),design=vitd.ph.comp.dsn, na.omit=T))
svyquantile(~urxmnp_c,vitd.ph.comp.dsn,c(0.25,0.5,0.75, 0.9, .95, 1),na.rm=TRUE)
quantile(vitd2_ph_comp$urxmnp[which(vitd2_ph_comp$urdmnplc==1)], c(0.25, 0.5, 0.75,  1))
quantile(vitd2_ph_comp$urxmnp, c(0.25, 0.5, 0.75, 0.90, 1))
table(vitd2_ph_comp$urxmnp[which(vitd2_ph_comp$urdmnplc==1)])

#MnOP
prop.table(table(vitd2_ph_comp$urdmoplc))
exp(svymean(~log(urxmop_c),design=vitd.ph.comp.dsn, na.omit=T))
svyquantile(~urxmop_C,vitd.ph.comp.dsn,c(0.25,0.5,0.75, 0.9, .95, 1),na.rm=TRUE)
quantile(vitd2_ph_comp$urxmop[which(vitd2_ph_comp$urdmoplc==1)], c(0.25, 0.5, 0.75,  1))
quantile(vitd2_ph_comp$urxmop, c(0.25, 0.5, 0.75, 0.90, 1))
table(vitd2_ph_comp$urxmop[which(vitd2_ph_comp$urdmoplc==1)])

#MCHP
prop.table(table(vitd2_ph_comp$urdmcplc))
exp(svymean(~log(urxmcp_c),design=vitd.ph.comp.dsn, na.omit=T))
svyquantile(~urxmcp_c,vitd.ph.comp.dsn,c(0.25,0.5,0.75, 0.9, .95, 1),na.rm=TRUE)
quantile(vitd2_ph_comp$urxmcp[which(vitd2_ph_comp$urdmcplc==1)], c(0.25, 0.5, 0.75,  1))
quantile(vitd2_ph_comp$urxmcp, c(0.25, 0.5, 0.75, 0.90, 1))
table(vitd2_ph_comp$urxmcp[which(vitd2_ph_comp$urdmcplc==1)])

#MnMP
prop.table(table(vitd2_ph_comp$urdmnmlc))
exp(svymean(~log(urxmnm_c),design=vitd.ph.comp.dsn, na.omit=T))
svyquantile(~urxmnm_c,vitd.ph.comp.dsn,c(0.25,0.5,0.75, 0.9, .95, 1),na.rm=TRUE)
quantile(vitd2_ph_comp$urxmnm[which(vitd2_ph_comp$urdmnmlc==1)], c(0.25, 0.5, 0.75,  1))
quantile(vitd2_ph_comp$urxmnm, c(0.25, 0.5, 0.75, 0.90, 1))
table(vitd2_ph_comp$urxmnm[which(vitd2_ph_comp$urdmnmlc==1)])

#MnMP
prop.table(table(vitd2_ph_comp$urdmnmlc))
exp(svymean(~log(urxmnm_c),design=vitd.ph.comp.dsn, na.omit=T))
svyquantile(~urxmnm_c,vitd.ph.comp.dsn,c(0.25,0.5,0.75, 0.9, .95, 1),na.rm=TRUE)
quantile(vitd2_ph_comp$urxmnm[which(vitd2_ph_comp$urdmnmlc==1)], c(0.25, 0.5, 0.75,  1))
quantile(vitd2_ph_comp$urxmnm, c(0.25, 0.5, 0.75, 0.90, 1))
table(vitd2_ph_comp$urxmnm[which(vitd2_ph_comp$urdmnmlc==1)])

#MCNP
prop.table(table(vitd2_ph_comp$urdcnplc))
exp(svymean(~log(urxcnp_c),design=vitd.ph.comp.dsn, na.omit=T))
svyquantile(~urxcnp_c,vitd.ph.comp.dsn,c(0.25,0.5,0.75, 0.9, .95, 1),na.rm=TRUE)
quantile(vitd2_ph_comp$urxcnp[which(vitd2_ph_comp$urdcnplc==1)], c(0.25, 0.5, 0.75,  1))
quantile(vitd2_ph_comp$urxcnp, c(0.25, 0.5, 0.75, 0.90, 1))
table(vitd2_ph_comp$urxcnp[which(vitd2_ph_comp$urdcnplc==1)])

#MCOP
prop.table(table(vitd2_ph_comp$urdcoplc))
exp(svymean(~log(urxcop_c),design=vitd.ph.comp.dsn, na.omit=T))
svyquantile(~urxcop_c,vitd.ph.comp.dsn,c(0.25,0.5,0.75, 0.9, .95, 1),na.rm=TRUE)
quantile(vitd2_ph_comp$urxcop[which(vitd2_ph_comp$urdcoplc==1)], c(0.25, 0.5, 0.75,  1))
quantile(vitd2_ph_comp$urxcop, c(0.25, 0.5, 0.75, 0.90, 1))
table(vitd2_ph_comp$urxcop[which(vitd2_ph_comp$urdcoplc==1)])

#BPA
prop.table(table(vitd2_ph_comp$urdbphlc))
exp(svymean(~log(urxbph_c),design=vitd.ph.comp.dsn, na.omit=T))
svyquantile(~urxbph_c,vitd.ph.comp.dsn,c(0.25,0.5,0.75, 0.9, .95, 1),na.rm=TRUE)
quantile(vitd2_ph_comp$urxbph[which(vitd2_ph_comp$urdbphlc==1)], c(0.25, 0.5, 0.75,  1))
quantile(vitd2_ph_comp$urxbph, c(0.25, 0.5, 0.75, 0.90, 1))
table(vitd2_ph_comp$urxbph[which(vitd2_ph_comp$urdbphlc==1)])

tapply(vitd2_ph_comp$vitamin_d, list(vitd2_ph_comp$sex, vitd2_ph_comp$age_cat), mean)
tapply(vitd2_ph_comp$urxbph, list(vitd2_ph_comp$sex, vitd2_ph_comp$age_cat), mean)

######################## SUBSET DATASETS BY AGE AND SEX #################################
vitd2.female<-subset(vitd2_ph_comp, sex==2)#2414
vitd2.female.dsn<-svydesign(id=~sdmvpsu, strata=~sdmvstra, weights=~weight6yr, data=vitd2.female, nest=T)

vitd2.male<-subset(vitd2_ph_comp, sex==1)#2310
vitd2.male.dsn<-svydesign(id=~sdmvpsu, strata=~sdmvstra, weights=~weight6yr, data=vitd2.male, nest=T)

vitd2.2039.female<-subset(vitd2_ph_comp, age>=20 & age<=39 & sex==2)#20-39 female; 897
vitd2.2039.female.dsn<-svydesign(id=~sdmvpsu, strata=~sdmvstra, weights=~weight6yr, data=vitd2.2039.female, nest=T)
summ(vitd2.2039.female$sdmvpsu)

vitd2.4059.female<-subset(vitd2_ph_comp, age>=40 & age<=59 & sex==2)#40-59 female; 760
vitd2.4059.female.dsn<-svydesign(id=~sdmvpsu, strata=~sdmvstra, weights=~weight6yr, data=vitd2.4059.female, nest=T)

vitd2.60.female<-subset(vitd2_ph_comp, age>=60 & sex==2)#60+ female; 757
vitd2.60.female.dsn<-svydesign(id=~sdmvpsu, strata=~sdmvstra, weights=~weight6yr, data=vitd2.60.female, nest=T)

vitd2.2039.male<-subset(vitd2_ph_comp, age>=20 & age<=39 & sex==1)#20-39 male; 749
vitd2.2039.male.dsn<-svydesign(id=~sdmvpsu, strata=~sdmvstra, weights=~weight6yr, data=vitd2.2039.male, nest=T)

vitd2.4059.male<-subset(vitd2_ph_comp, age>=40 & age<=59 & sex==1)#40-59 male; 773
vitd2.4059.male.dsn<-svydesign(id=~sdmvpsu, strata=~sdmvstra, weights=~weight6yr, data=vitd2.4059.male, nest=T)

vitd2.60.male<-subset(vitd2_ph_comp, age>=60 & sex==1)#60+ male; 788
vitd2.60.male.dsn<-svydesign(id=~sdmvpsu, strata=~sdmvstra, weights=~weight6yr, data=vitd2.60.male, nest=T)

vitd2.cycle45<-subset(vitd2_ph_comp, cycle==4 | cycle==5)
vitd2.cycle45$weight4yr<-vitd2.cycle45$weight/2 #weight variable created in SAS
vitd2.cycle45.dsn<-svydesign(id=~sdmvpsu, strata=~sdmvstra, weights=~weight4yr, data=vitd2.cycle45, nest=T)

vitd2.black<-subset(vitd2_ph_comp, race_cat==3)
vitd2.black.dsn<-svydesign(id=~sdmvpsu, strata=~sdmvstra, weights=~weight6yr, data=vitd2.black, nest=T)

######################## UNIVARIARTE ANALYSES #################################

hist(vitd2$vitamin_d)#normally distributed
hist(vitd2$urxucr)#skewed, log transform
hist(vitd2$urxmep)#skewed, log transform
hist(vitd2$bpa)#skewed, log transform
hist(vitd2$cot)#skewed, log transform

######################## BIVARIARTE ANALYSES #################################
summary(svyglm(vitamin_d~factor(age_cat), vitd.ph.comp.dsn, na.action=na.omit)) #p<0.05
summary(svyglm(vitamin_d~factor(sex), vitd.dsn, na.action=na.omit)) #p=0.09; may need to stratify by sex
summary(svyglm(vitamin_d~factor(race), vitd.dsn, na.action=na.omit)) #p<0.01
summary(svyglm(vitamin_d~pir, vitd.dsn, na.action=na.omit)) #p<0.05
summary(svyglm(vitamin_d~factor(edu_cat), vitd.dsn, na.action=na.omit)) #p<0.05
summary(svyglm(vitamin_d~alc, vitd.dsn, na.action=na.omit))#p>0.05
summary(svyglm(vitamin_d~bmi, vitd.dsn, na.action=na.omit))#p<0.05
summary(svyglm(vitamin_d~factor(smk), vitd.dsn, na.action=na.omit))#p<0.05
summary(svyglm(vitamin_d~factor(milk), vitd.dsn, na.action=na.omit))#p<0.05
summary(svyglm(vitamin_d~factor(sun), vitd.dsn, na.action=na.omit))#p<0.05
summary(svyglm(vitamin_d~factor(fish), vitd.dsn, na.action=na.omit)#p<0.05
summary(svyglm(vitamin_d~factor(season), vitd.dsn, na.action=na.omit))#p<0.05
summary(svyglm(vitamin_d~factor(vitdsup), vitd.dsn, na.action=na.omit))#p<0.05
summary(svyglm(vitamin_d~log(cot), vitd.dsn, na.action=na.omit))#p<0.05
summary(svyglm(vitamin_d~factor(cycle), vitd.dsn, na.action=na.omit))#p<0.05
summary(svyglm(vitamin_d~factor(smk), vitd.ph.comp.dsn, na.action=na.omit))

summary(svyglm(log(urxmbp)~age, vitd.dsn, na.action=na.omit)) #p<0.05
summary(svyglm(log(urxmbp)~factor(sex), vitd.dsn, na.action=na.omit)) #p=0.02; may need to stratify by sex
summary(svyglm(log(urxmbp)~factor(race), vitd.dsn, na.action=na.omit)) #p<0.05
summary(svyglm(log(urxmbp)~pir, vitd.dsn, na.action=na.omit)) #p<0.05
summary(svyglm(log(urxmbp)~factor(edu_cat), vitd.dsn, na.action=na.omit)) #p>0.05
summary(svyglm(log(urxmbp)~alc, vitd.dsn, na.action=na.omit))#p<0.05
summary(svyglm(log(urxmbp)~bmi, vitd.ph.comp.dsn, na.action=na.omit))#p<0.05
summary(svyglm(log(urxmbp)~factor(smk), vitd.dsn, na.action=na.omit))#p<0.05
summary(svyglm(log(urxmbp)~factor(milk), vitd.dsn, na.action=na.omit))#p=0.7
summary(svyglm(log(urxmep)~factor(sun), vitd2_ph_comp, na.action=na.omit))#p>05 except for those who didn't use sunscreen - had higher phthalates?!
summary(svyglm(log(urxmbp)~factor(fish), vitd.dsn, na.action=na.omit))#p=0.05
summary(svyglm(log(urxmbp)~factor(season), vitd.dsn, na.action=na.omit))#p<0.05
summary(svyglm(log(urxmbp)~factor(vitdsup), vitd.dsn, na.action=na.omit))#p<0.05
summary(svyglm(log(urxmbp)~factor(cycle), vitd.dsn, na.action=na.omit))#p<0.05
summary(svyglm(log(urxmbp)~factor(smk), vitd.ph.dsn, na.action=na.omit))

summary(lm(log(urxmep)~factor(sun), data=vitd2_ph_comp))

###check correlation to see which phthalate to use for testing these bivariate associations
cor.test(vitd2.adults$vitamin_d, log(vitd2.adults$urxmcp), use="complete.obs")#p<0.001

######################## REGRESSION ANALYSES #################################

#run crude models
summary(lm(vitamin_d~log(urxmbp)+log(urxucr), data=vitd2.adults))#p=0.5
summary(lm(vitamin_d~log(urxmcp)+log(urxucr), data=vitd2.adults))#p<0.05
summary(lm(vitamin_d~log(urxmhp)+log(urxucr), data=vitd2.adults))#p<0.05
summary(lm(vitamin_d~log(urxmnp)+log(urxucr), data=vitd2.adults))#p=0.16
summary(lm(vitamin_d~log(urxmop)+log(urxucr), data=vitd2.adults))#p<0.05
summary(lm(vitamin_d~log(urxmzp)+log(urxucr), data=vitd2.adults))#p=0.09
summary(lm(vitamin_d~log(urxmnm)+log(urxucr), data=vitd2.adults))#p=0.6
summary(lm(vitamin_d~log(urxmc1)+log(urxucr), data=vitd2.adults))#p<0.05
summary(lm(vitamin_d~log(urxmhh)+log(urxucr), data=vitd2.adults))#p=0.049
summary(lm(vitamin_d~log(urxmoh)+log(urxucr), data=vitd2.adults))#p=0.16
summary(lm(vitamin_d~log(urxmib)+log(urxucr), data=vitd2.adults))#p=0.02
summary(lm(vitamin_d~log(urxecp)+log(urxucr), data=vitd2.adults))#p=0.40

#run models with covariates - step-wise approach
summary(svyglm(vitamin_d~log(vitd2.women$urxmcp)+log(vitd2.women$urxucr)+ vitd2.women$age + factor(vitd2.women$race)+ vitd2.women$bmi+ log(vitd2.women$cot) + factor(vitd2.women$suppl)
           + factor(vitd2.women$season) + factor(cycle), data=vitd2.women))#p<0.05

hist(vitd2_ph_comp$vitamin_d[which(vitd2_ph_comp$vitdsu==1)])
###########calculate IQR for conversion##########

svyquantile(~urxbph,vitd.ph.comp.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxmhp,vitd.ph.comp.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxmhh,vitd.ph.comp.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxmoh,vitd.ph.comp.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxecp,vitd.ph.comp.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~dehpsum,vitd.ph.comp.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxmbp,vitd.ph.comp.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxmib,vitd.ph.comp.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxmep,vitd.ph.comp.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxmzp,vitd.ph.comp.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxmc1,vitd.ph.comp.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxcnp,vitd.ph.comp.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxcop,vitd.ph.comp.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~vitamin_d,vitd.ph.comp.dsn,c(.5),na.rm=TRUE)

svyquantile(~urxbph,vitd2.female.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxmhp,vitd2.female.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxmhh,vitd2.female.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxmoh,vitd2.female.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxecp,vitd2.female.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~dehpsum,vitd2.female.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxmbp,vitd2.female.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxmib,vitd2.female.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxmep,vitd2.female.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxmzp,vitd2.female.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxmc1,vitd2.female.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxcnp,vitd2.female.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxcop,vitd2.female.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~vitamin_d,vitd2.female.dsn,c(.5),na.rm=TRUE)

svyquantile(~urxbph,vitd2.male.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxmhp,vitd2.male.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxmhh,vitd2.male.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxmoh,vitd2.male.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxecp,vitd2.male.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~dehpsum,vitd2.male.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxmbp,vitd2.male.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxmib,vitd2.male.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxmep,vitd2.male.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxmzp,vitd2.male.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxmc1,vitd2.male.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxcnp,vitd2.male.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~urxcop,vitd2.male.dsn,c(0.25,0.75),na.rm=TRUE)
svyquantile(~vitamin_d,vitd2.male.dsn,c(.5),na.rm=TRUE)


########### Weighted Models ###############

#########Overall Population#########
options(survey.lonely.psu = "adjust")#this is the conservative approach for strata with only 1 PSU
#from UCLA's stats department:The most conservative approach would be to center any single-PSU strata around the sample grand mean rather than the stratum mean

overall<-svyglm(vitamin_d~log(urxbph)+log(urxucr)+ age + factor(sex) + factor(race) + bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd.ph.comp.dsn, na.action=na.omit)
summary(overall)
nobs(overall)
confint(overall, level=0.95)

overall<-svyglm(vitamin_d~log(urxmhp)+log(urxucr)+ age + factor(sex)+ factor(race) + bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd.ph.comp.dsn, na.action=na.omit)
summary(overall)
nobs(overall)
confint(overall, level=0.95)

overall<-svyglm(vitamin_d~log(urxmhh)+log(urxucr)+ age + factor(sex)+ factor(race) + bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd.ph.comp.dsn, na.action=na.omit)
summary(overall)
confint(overall, level=0.95)

overall<-svyglm(vitamin_d~log(urxmoh)+log(urxucr)+ age + factor(sex) + factor(race) +  bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd.ph.comp.dsn, na.action=na.omit)
summary(overall)
confint(overall, level=0.95)

overall<-svyglm(vitamin_d~log(urxecp)+log(urxucr)+ age + factor(sex) + factor(race) + bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd.ph.comp.dsn, na.action=na.omit)
summary(overall)
confint(overall, level=0.95)

overall<-svyglm(vitamin_d~log(dehpsum)+log(urxucr)+ age + factor(sex) + factor(race)+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd.ph.comp.dsn, na.action=na.omit)
summary(overall)
confint(overall, level=0.95)

overall<-svyglm(vitamin_d~log(urxmbp)+log(urxucr)+ age + factor(sex) + factor(race) + bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd.ph.comp.dsn, na.action=na.omit)
summary(overall)
confint(overall, level=0.95)

overall<-svyglm(vitamin_d~log(urxmib)+log(urxucr)+ age + factor(sex) + factor(race) + bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd.ph.comp.dsn, na.action=na.omit)
summary(overall)
nobs(overall)
confint(overall, level=0.95)

overall<-svyglm(vitamin_d~log(urxmep)+log(urxucr)+ age + factor(sex) + factor(race) +  bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd.ph.comp.dsn, na.action=na.omit)
summary(overall)
confint(overall, level=0.95)

overall<-svyglm(vitamin_d~log(urxmzp)+log(urxucr)+ age + factor(sex) + factor(race) + bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd.ph.comp.dsn, na.action=na.omit)
summary(overall)
confint(overall, level=0.95)

overall<-svyglm(vitamin_d~log(urxmc1)+log(urxucr)+ age + factor(sex) + factor(race) + bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd.ph.comp.dsn, na.action=na.omit)
summary(overall)
confint(overall, level=0.95)

overall<-svyglm(vitamin_d~log(urxcnp)+log(urxucr)+ age + factor(sex) + factor(race) + bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd.ph.comp.dsn, na.action=na.omit)
summary(overall)
confint(overall, level=0.95)

overall<-svyglm(vitamin_d~log(urxcop)+log(urxucr)+ age + factor(sex) + factor(race) + bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd.ph.comp.dsn, na.action=na.omit)
summary(overall)
nobs(overall)
confint(overall, level=0.95)

#########Overall Population with quartiles of exposure#########

#####Create quartiles of exposure###
vitd2_ph_comp$urxbph_q<-cut2(vitd2_ph_comp$urxbph_c, g=5)
table(vitd2_ph_comp$urxbph_q)
table(as.integer(vitd2_ph_comp$urxbph_q))

vitd2_ph_comp$urxmhp_q<-cut2(vitd2_ph_comp$urxmhp_c, g=5)
table(vitd2_ph_comp$urxmhp_q)

vitd2_ph_comp$urxmhh_q<-cut2(vitd2_ph_comp$urxmhh_c, g=5)
table(vitd2_ph_comp$urxmhh_q)

vitd2_ph_comp$urxmoh_q<-cut2(vitd2_ph_comp$urxmoh_c, g=5)
table(vitd2_ph_comp$urxmoh_q)

vitd2_ph_comp$urxecp_q<-cut2(vitd2_ph_comp$urxecp_c, g=5)
table(vitd2_ph_comp$urxecp_q)

vitd.ph.comp.dsn<-svydesign(id=~sdmvpsu, strata=~sdmvstra, weights=~weight6yr, data=vitd2_ph_comp, nest=T) #need to create complete dataset
save(vitd2_ph_comp, file="vitd2_ph_comp.rda")

###models 
load("vitd2_ph_comp.rda")
overall<-svyglm(vitamin_d~as.factor(urxbph_q)+ age + factor(sex) + factor(race)+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd.ph.comp.dsn, na.action=na.omit)
summary(overall)
confint(overall, level=0.95)

overall1<-svyglm(vitamin_d~as.numeric(urxmhp_q)+ age + factor(sex) + factor(race)+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd.ph.comp.dsn, na.action=na.omit)
summary(overall1)
confint(overall1, level=0.95)

overall2<-svyglm(vitamin_d~as.factor(urxmhh_q)+ age + factor(sex)+ factor(race)+  bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd.ph.comp.dsn, na.action=na.omit)
summary(overall2)
confint(overall2, level=0.95)

overall<-svyglm(vitamin_d~as.factor(urxmoh_q)+age + factor(sex) + factor(race)+  bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd.ph.comp.dsn, na.action=na.omit)
summary(overall)
confint(overall, level=0.95)

overall<-svyglm(vitamin_d~as.numeric(urxecp_q)+ age + factor(sex) + factor(race)+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd.ph.comp.dsn, na.action=na.omit)
summary(overall)
confint(overall, level=0.95)

###############Women only###############
rm(vitd2.male)

female<-svyglm(vitamin_d~log(urxbph)+log(urxucr)+ age + factor(race)+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.female.dsn, na.action=na.omit)
summary(female)
nobs(female)
confint(female, level=0.95)

female<-svyglm(vitamin_d~log(urxmhp)+log(urxucr)+ age + factor(race)+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.female.dsn, na.action=na.omit)
summary(female)
confint(female, level=0.95)

female<-svyglm(vitamin_d~log(urxmhh)+log(urxucr)+ age + factor(race)+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.female.dsn, na.action=na.omit)
summary(female)
confint(female, level=0.95)

female<-svyglm(vitamin_d~log(urxmoh)+log(urxucr)+ age + factor(race)+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.female.dsn, na.action=na.omit)
summary(female)
confint(female, level=0.95)

female<-svyglm(vitamin_d~log(urxecp)+log(urxucr)+ age + factor(race)+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.female.dsn, na.action=na.omit)
summary(female)
confint(female, level=0.95)

female<-svyglm(vitamin_d~log(dehpsum)+log(urxucr)+ age + factor(race)+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.female.dsn, na.action=na.omit)
summary(female)
nobs(female)
confint(female, level=0.95)

female<-svyglm(vitamin_d~log(urxmbp)+log(urxucr)+ age + factor(race)+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.female.dsn, na.action=na.omit)
summary(female)
confint(female, level=0.95)

female<-svyglm(vitamin_d~log(urxmib)+log(urxucr)+ age  + factor(race) + bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.female.dsn, na.action=na.omit)
summary(female)
confint(female, level=0.95)

female<-svyglm(vitamin_d~log(urxmep)+log(urxucr)+ age  + factor(race)+  bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.female.dsn, na.action=na.omit)
summary(female)
confint(female, level=0.95)

female<-svyglm(vitamin_d~log(urxmzp)+log(urxucr)+ age + factor(race)+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.female.dsn, na.action=na.omit)
summary(female)
nobs(female)
confint(female, level=0.95)

female<-svyglm(vitamin_d~log(urxmc1)+log(urxucr)+ age  + factor(race)+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.female.dsn, na.action=na.omit)
summary(female)
confint(female, level=0.95)

female<-svyglm(vitamin_d~log(urxcnp)+log(urxucr)+ age  + factor(race)+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.female.dsn, na.action=na.omit)
summary(female)
confint(female, level=0.95)

female<-svyglm(vitamin_d~log(urxcop)+log(urxucr)+ age  + factor(race)+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.female.dsn, na.action=na.omit)
summary(female)
confint(female, level=0.95)

###############Men only###############

male<-svyglm(vitamin_d~log(urxbph)+log(urxucr)+ age + factor(race)+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.male.dsn, na.action=na.omit)
summary(male)
nobs(male)
confint(male, level=0.95)

male<-svyglm(vitamin_d~log(urxmhp)+log(urxucr)+ age + factor(race)+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.male.dsn, na.action=na.omit)
summary(male)
confint(male, level=0.95)

male<-svyglm(vitamin_d~log(urxmhh)+log(urxucr)+ age + factor(race)+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.male.dsn, na.action=na.omit)
summary(male)
confint(male, level=0.95)

male<-svyglm(vitamin_d~log(urxmoh)+log(urxucr)+ age + factor(race)+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.male.dsn, na.action=na.omit)
summary(male)
confint(male, level=0.95)

male<-svyglm(vitamin_d~log(urxecp)+log(urxucr)+ age + factor(race)+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.male.dsn, na.action=na.omit)
summary(male)
confint(male, level=0.95)

male<-svyglm(vitamin_d~log(dehpsum)+log(urxucr)+ age + factor(race) + bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.male.dsn, na.action=na.omit)
summary(male)
nobs(overall)
confint(male, level=0.95)

male<-svyglm(vitamin_d~log(urxmbp)+log(urxucr)+ age + factor(race)+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.male.dsn, na.action=na.omit)
summary(male)
confint(male, level=0.95)

male<-svyglm(vitamin_d~log(urxmib)+log(urxucr)+ age  + factor(race) + bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.male.dsn, na.action=na.omit)
summary(male)
nobs(male)
confint(male, level=0.95)

male<-svyglm(vitamin_d~log(urxmep)+log(urxucr)+ age  + factor(race)+  bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.male.dsn, na.action=na.omit)
summary(male)
confint(male, level=0.95)

male<-svyglm(vitamin_d~log(urxmzp)+log(urxucr)+ age + factor(race)+  bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.male.dsn, na.action=na.omit)
summary(male)
confint(male, level=0.95)

male<-svyglm(vitamin_d~log(urxmc1)+log(urxucr)+ age  + factor(race)+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.male.dsn, na.action=na.omit)
summary(male)
confint(male, level=0.95)

male<-svyglm(vitamin_d~log(urxcnp)+log(urxucr)+ age  + factor(race)+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.male.dsn, na.action=na.omit)
summary(male)
confint(male, level=0.95)

male<-svyglm(vitamin_d~log(urxcop)+log(urxucr)+ age  + factor(race)+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.male.dsn, na.action=na.omit)
summary(male)
confint(male, level=0.95)

#########20-39 Women#########
options(survey.lonely.psu = "adjust")#this is the conservative approach for strata with only 1 PSU
#from UCLA's stats department:The most conservative approach would be to center any single-PSU strata around the sample grand mean rather than the stratum mean

f.2039<-svyglm(vitamin_d~log(urxbph)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.2039.female.dsn, na.action=na.omit)
summary(f.2039)
confint(f.2039, level=0.95)

f.2039<-svyglm(vitamin_d~log(urxmhp)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.2039.female.dsn, na.action=na.omit)
summary(f.2039)
confint(f.2039, level=0.95)

f.2039<-svyglm(vitamin_d~log(urxmhh)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.2039.female.dsn, na.action=na.omit)
summary(f.2039)
confint(f.2039, level=0.95)

f.2039<-svyglm(vitamin_d~log(urxmoh)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.2039.female.dsn, na.action=na.omit)
summary(f.2039)
confint(f.2039, level=0.95)

f.2039<-svyglm(vitamin_d~log(urxecp)+log(urxucr)+ age + factor(race)+ pir + bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.2039.female.dsn, na.action=na.omit)
summary(f.2039)
confint(f.2039, level=0.95)

f.2039<-svyglm(vitamin_d~log(dehpsum)+log(urxucr)+ age + factor(race)+ pir + bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.2039.female.dsn, na.action=na.omit)
summary(f.2039)
confint(f.2039, level=0.95)

f.2039<-svyglm(vitamin_d~log(urxmbp)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.2039.female.dsn, na.action=na.omit)
summary(f.2039)
confint(f.2039, level=0.95)

f.2039<-svyglm(vitamin_d~log(urxmib)+log(urxucr)+ age + factor(race)+ pir + bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.2039.female.dsn, na.action=na.omit)
summary(f.2039)
confint(f.2039, level=0.95)

f.2039<-svyglm(vitamin_d~log(urxmep)+log(urxucr)+ age + factor(race)+ pir+  bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.2039.female.dsn, na.action=na.omit)
summary(f.2039)
confint(f.2039, level=0.95)

f.2039<-svyglm(vitamin_d~log(urxmzp)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.2039.female.dsn, na.action=na.omit)
summary(f.2039)
confint(f.2039, level=0.95)

f.2039<-svyglm(vitamin_d~log(urxmc1)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.2039.female.dsn, na.action=na.omit)
summary(f.2039)
confint(f.2039, level=0.95)

f.2039<-svyglm(vitamin_d~log(urxcnp)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.2039.female.dsn, na.action=na.omit)
summary(f.2039)
confint(f.2039, level=0.95)

f.2039<-svyglm(vitamin_d~log(urxcop)+log(urxucr)+ age + factor(race)+ pir + bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.2039.female.dsn, na.action=na.omit)
summary(f.2039)
confint(f.2039, level=0.95)

#########40-59 Women#########
f.4059<-svyglm(vitamin_d~log(urxbph)+log(urxucr)+ age + factor(race)+ pir +bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.4059.female.dsn, na.action=na.omit)
summary(f.4059)
confint(f.4059, level=0.95)

f.4059<-svyglm(vitamin_d~log(urxmhp)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.4059.female.dsn, na.action=na.omit)
summary(f.4059)
confint(f.4059, level=0.95)

f.4059<-svyglm(vitamin_d~log(urxmhh)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.4059.female.dsn, na.action=na.omit)
summary(f.4059)
confint(f.4059, level=0.95)

f.4059<-svyglm(vitamin_d~log(urxmoh)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.4059.female.dsn, na.action=na.omit)
summary(f.4059)
confint(f.4059, level=0.95)

f.4059<-svyglm(vitamin_d~log(urxecp)+log(urxucr)+ age + factor(race)+ pir + bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.4059.female.dsn, na.action=na.omit)
summary(f.4059)
confint(f.4059, level=0.95)

f.4059<-svyglm(vitamin_d~log(dehpsum)+log(urxucr)+ age + factor(race)+ pir + bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.4059.female.dsn, na.action=na.omit)
summary(f.4059)
confint(f.4059, level=0.95)

f.4059<-svyglm(vitamin_d~log(urxmbp)+log(urxucr)+ age + factor(race)+ pir + bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.4059.female.dsn, na.action=na.omit)
summary(f.4059)
confint(f.4059, level=0.95)

f.4059<-svyglm(vitamin_d~log(urxmib)+log(urxucr)+ age + factor(race)+ pir + bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.4059.female.dsn, na.action=na.omit)
summary(f.4059)
confint(f.4059, level=0.95)

f.4059<-svyglm(vitamin_d~log(urxmep)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.4059.female.dsn, na.action=na.omit)
summary(f.4059)
confint(f.4059, level=0.95)

f.4059<-svyglm(vitamin_d~log(urxmzp)+log(urxucr)+ age + factor(race)+ pir + bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.4059.female.dsn, na.action=na.omit)
summary(f.4059)
confint(f.4059, level=0.95)

f.4059<-svyglm(vitamin_d~log(urxmc1)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.4059.female.dsn, na.action=na.omit)
summary(f.4059)
confint(f.4059, level=0.95)

f.4059<-svyglm(vitamin_d~log(urxcnp)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.4059.female.dsn, na.action=na.omit)
summary(f.4059)
confint(f.4059, level=0.95)

f.4059<-svyglm(vitamin_d~log(urxcop)+log(urxucr)+ age + factor(race)+ pir + bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.4059.female.dsn, na.action=na.omit)
summary(f.4059)
confint(f.4059, level=0.95)

######### 60+ Women#########
f.60<-svyglm(vitamin_d~log(urxbph)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.60.female.dsn, na.action=na.omit)
summary(f.60)
confint(f.60, level=0.95)

f.60<-svyglm(vitamin_d~log(urxmhp)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.60.female.dsn, na.action=na.omit)
summary(f.60)
confint(f.60, level=0.95)

f.60<-svyglm(vitamin_d~log(urxmhh)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.60.female.dsn, na.action=na.omit)
summary(f.60)
confint(f.60, level=0.95)

f.60<-svyglm(vitamin_d~log(urxmoh)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.60.female.dsn, na.action=na.omit)
summary(f.60)
confint(f.60, level=0.95)

f.60<-svyglm(vitamin_d~log(urxecp)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.60.female.dsn, na.action=na.omit)
summary(f.60)
confint(f.60, level=0.95)

f.60<-svyglm(vitamin_d~log(dehpsum)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.60.female.dsn, na.action=na.omit)
summary(f.60)
confint(f.60, level=0.95)

f.60<-svyglm(vitamin_d~log(urxmbp)+log(urxucr)+ age + factor(race)+ pir+bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.60.female.dsn, na.action=na.omit)
summary(f.60)
confint(f.60, level=0.95)

f.60<-svyglm(vitamin_d~log(urxmib)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.60.female.dsn, na.action=na.omit)
summary(f.60)
confint(f.60, level=0.95)

f.60<-svyglm(vitamin_d~log(urxmep)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.60.female.dsn, na.action=na.omit)
summary(f.60)
confint(f.60, level=0.95)

f.60<-svyglm(vitamin_d~log(urxmzp)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.60.female.dsn, na.action=na.omit)
summary(f.60)
confint(f.60, level=0.95)

f.60<-svyglm(vitamin_d~log(urxmc1)+log(urxucr)+ age + factor(race)+ pir+bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.60.female.dsn, na.action=na.omit)
summary(f.60)
confint(f.60, level=0.95)

f.60<-svyglm(vitamin_d~log(urxcnp)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.60.female.dsn, na.action=na.omit)
summary(f.60)
confint(f.60, level=0.95)

f.60<-svyglm(vitamin_d~log(urxcop)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.60.female.dsn, na.action=na.omit)
summary(f.60)
confint(f.60, level=0.95)

#########20-39 Men#########
m.2039<-svyglm(vitamin_d~log(urxbph)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.2039.male.dsn, na.action=na.omit)
summary(m.2039)
confint(m.2039, level=0.95)

m.2039<-svyglm(vitamin_d~log(urxmhp)+log(urxucr)+ age + factor(race)+ pir + bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.2039.male.dsn, na.action=na.omit)
summary(m.2039)
confint(m.2039, level=0.95)

m.2039<-svyglm(vitamin_d~log(urxmhh)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.2039.male.dsn, na.action=na.omit)
summary(m.2039)
confint(m.2039, level=0.95)

m.2039<-svyglm(vitamin_d~log(urxmoh)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.2039.male.dsn, na.action=na.omit)
summary(m.2039)
confint(m.2039, level=0.95)

m.2039<-svyglm(vitamin_d~log(urxecp)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.2039.male.dsn, na.action=na.omit)
summary(m.2039)
confint(m.2039, level=0.95)

m.2039<-svyglm(vitamin_d~log(dehpsum)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.2039.male.dsn, na.action=na.omit)
summary(m.2039)
confint(m.2039, level=0.95)

m.2039<-svyglm(vitamin_d~log(urxmbp)+log(urxucr)+ age + factor(race)+ pir+bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.2039.male.dsn, na.action=na.omit)
summary(m.2039)
confint(m.2039, level=0.95)

m.2039<-svyglm(vitamin_d~log(urxmib)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.2039.male.dsn, na.action=na.omit)
summary(m.2039)
confint(m.2039, level=0.95)

m.2039<-svyglm(vitamin_d~log(urxmep)+log(urxucr)+ age + factor(race)+ bmi+ pir+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.2039.male.dsn, na.action=na.omit)
summary(m.2039)
confint(m.2039, level=0.95)

m.2039<-svyglm(vitamin_d~log(urxmzp)+log(urxucr)+ age + factor(race)+ bmi+ pir +log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.2039.male.dsn, na.action=na.omit)
summary(m.2039)
confint(m.2039, level=0.95)

m.2039<-svyglm(vitamin_d~log(urxmc1)+log(urxucr)+ age + factor(race)+ bmi+ pir + log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.2039.male.dsn, na.action=na.omit)
summary(m.2039)
confint(m.2039, level=0.95)

m.2039<-svyglm(vitamin_d~log(urxcnp)+log(urxucr)+ age + factor(race)+ bmi+ pir + log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.2039.male.dsn, na.action=na.omit)
summary(m.2039)
confint(m.2039, level=0.95)

m.2039<-svyglm(vitamin_d~log(urxcop)+log(urxucr)+ age + factor(race)+ bmi+ pir + log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.2039.male.dsn, na.action=na.omit)
summary(m.2039)
confint(m.2039, level=0.95)

#########40-59 Men#########
m.4059<-svyglm(vitamin_d~log(urxbph)+log(urxucr)+ age + factor(race)+ bmi+ pir+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.4059.male.dsn, na.action=na.omit)
summary(m.4059)
confint(m.4059, level=0.95)

m.4059<-svyglm(vitamin_d~log(urxmhp)+log(urxucr)+ age + factor(race)+ bmi+ pir+log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.4059.male.dsn, na.action=na.omit)
summary(m.4059)
confint(m.4059, level=0.95)

m.4059<-svyglm(vitamin_d~log(urxmhh)+log(urxucr)+ age + factor(race)+ bmi+ pir+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.4059.male.dsn, na.action=na.omit)
summary(m.4059)
confint(m.4059, level=0.95)

m.4059<-svyglm(vitamin_d~log(urxmoh)+log(urxucr)+ age + factor(race)+ bmi+ pir+log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.4059.male.dsn, na.action=na.omit)
summary(m.4059)
confint(m.4059, level=0.95)

m.4059<-svyglm(vitamin_d~log(urxecp)+log(urxucr)+ age + factor(race)+ bmi+ pir+log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.4059.male.dsn, na.action=na.omit)
summary(m.4059)
confint(m.4059, level=0.95)

m.4059<-svyglm(vitamin_d~log(dehpsum)+log(urxucr)+ age + factor(race)+ bmi+ pir+log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.4059.male.dsn, na.action=na.omit)
summary(m.4059)
confint(m.4059, level=0.95)

m.4059<-svyglm(vitamin_d~log(urxmbp)+log(urxucr)+ age + factor(race)+ bmi+ pir+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.4059.male.dsn, na.action=na.omit)
summary(m.4059)
confint(m.4059, level=0.95)

m.4059<-svyglm(vitamin_d~log(urxmib)+log(urxucr)+ age + factor(race)+ bmi+ pir+log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.4059.male.dsn, na.action=na.omit)
summary(m.4059)
confint(m.4059, level=0.95)

m.4059<-svyglm(vitamin_d~log(urxmep)+log(urxucr)+ age + factor(race)+ bmi+pir+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.4059.male.dsn, na.action=na.omit)
summary(m.4059)
confint(m.4059, level=0.95)

m.4059<-svyglm(vitamin_d~log(urxmzp)+log(urxucr)+ age + factor(race)+ bmi+ pir+log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.4059.male.dsn, na.action=na.omit)
summary(m.4059)
confint(m.4059, level=0.95)

m.4059<-svyglm(vitamin_d~log(urxmc1)+log(urxucr)+ age + factor(race)+ bmi+ pir+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.4059.male.dsn, na.action=na.omit)
summary(m.4059)
confint(m.4059, level=0.95)

m.4059<-svyglm(vitamin_d~log(urxcnp)+log(urxucr)+ age + factor(race)+ bmi+ pir+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.4059.male.dsn, na.action=na.omit)
summary(m.4059)
confint(m.4059, level=0.95)

m.4059<-svyglm(vitamin_d~log(urxcop)+log(urxucr)+ age + factor(race)+ bmi+ pir+log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.4059.male.dsn, na.action=na.omit)
summary(m.4059)
confint(m.4059, level=0.95)

######### 60+ Men#########
m.60<-svyglm(vitamin_d~log(urxbph)+log(urxucr)+ age + factor(race)+ bmi+ pir+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.60.male.dsn, na.action=na.omit)
summary(m.60)
confint(m.60, level=0.95)

m.60<-svyglm(vitamin_d~log(urxmhp)+log(urxucr)+ age + factor(race)+ bmi+pir+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.60.male.dsn, na.action=na.omit)
summary(m.60)
confint(m.60, level=0.95)

m.60<-svyglm(vitamin_d~log(urxmhh)+log(urxucr)+ age + factor(race)+ bmi+ pir+log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.60.male.dsn, na.action=na.omit)
summary(m.60)
confint(m.60, level=0.95)

m.60<-svyglm(vitamin_d~log(urxmoh)+log(urxucr)+ age + factor(race)+ bmi+pir+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.60.male.dsn, na.action=na.omit)
summary(m.60)
confint(m.60, level=0.95)

m.60<-svyglm(vitamin_d~log(urxecp)+log(urxucr)+ age + factor(race)+ bmi+ pir+log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.60.male.dsn, na.action=na.omit)
summary(m.60)
confint(m.60, level=0.95)

m.60<-svyglm(vitamin_d~log(dehpsum)+log(urxucr)+ age + factor(race)+ bmi+pir+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.60.male.dsn, na.action=na.omit)
summary(m.60)
confint(m.60, level=0.95)

m.60<-svyglm(vitamin_d~log(urxmbp)+log(urxucr)+ age + factor(race)+ bmi+ pir+log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.60.male.dsn, na.action=na.omit)
summary(m.60)
confint(m.60, level=0.95)

m.60<-svyglm(vitamin_d~log(urxmib)+log(urxucr)+ age + factor(race)+ bmi+pir+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.60.male.dsn, na.action=na.omit)
summary(m.60)
confint(m.60, level=0.95)

m.60<-svyglm(vitamin_d~log(urxmep)+log(urxucr)+ age + factor(race)+ bmi+pir+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.60.male.dsn, na.action=na.omit)
summary(m.60)
confint(m.60, level=0.95)

m.60<-svyglm(vitamin_d~log(urxmzp)+log(urxucr)+ age + factor(race)+ bmi+ pir+log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.60.male.dsn, na.action=na.omit)
summary(m.60)
confint(m.60, level=0.95)

m.60<-svyglm(vitamin_d~log(urxmc1)+log(urxucr)+ age + factor(race)+ bmi+pir+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.60.male.dsn, na.action=na.omit)
summary(m.60)
confint(m.60, level=0.95)

m.60<-svyglm(vitamin_d~log(urxcnp)+log(urxucr)+ age + factor(race)+ bmi+pir+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.60.male.dsn, na.action=na.omit)
summary(m.60)
confint(m.60, level=0.95)

m.60<-svyglm(vitamin_d~log(urxcop)+log(urxucr)+ age + factor(race)+ bmi+pir+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.60.male.dsn, na.action=na.omit)
summary(m.60)
confint(m.60, level=0.95)

###run multiple models
###females only
female.2039.models=function(x){
  vitd.x.female<-lm(vitamin_d~log(x)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + dsqtvd + factor(cycle), data=vitd2.2039.female, na.action=na.omit)
  b.x<-summary(vitd.x.female)$coef[2,1]
  se.x<-summary(vitd.x.female)$coef[2,2]
  p.x<-summary(vitd.x.female)$coef[2,4]
  lo95ci<- exp(b.x-1.96*se.x)
  hi95ci<- exp(b.x+1.96*se.x)
  #nobs <- (vitd.x.female)
  all.x<-cbind(b.x,se.x,p.x, lo95ci, hi95ci)#, nobs)
  return(all.x)}
female.results<-rbind(female.2039.models(vitd2.2039.female$urxbph),female.2039.models(vitd2.2039.female$urxmhp),female.2039.models(vitd2.2039.female$urxmhh),female.2039.models(vitd2.2039.female$urxmoh),female.2039.models(vitd2.2039.female$urxecp),female.2039.models(vitd2.2039.female$dehpsum),
female.2039.models(vitd2.2039.female$urxmbp),female.2039.models(vitd2.2039.female$urxmib),female.2039.models(vitd2.2039.female$urxmep),female.2039.models(vitd2.2039.female$urxmzp),female.2039.models(vitd2.2039.female$urxmc1),female.2039.models(vitd2.2039.female$urxcnp), 
female.2039.models(vitd2.2039.female$urxcop))
row.names(female.results)<-c("BPA","MEHP","MEHHP","MEOHP","MECPP","DEHPsum","MBP","MiBP","MEP","MBzP","MCPP","MCNP","MCOP")
female.results
write.csv(female.results, file="female_2039.csv")

female.4059.models=function(x){
  vitd.x.female<-lm(vitamin_d~log(x)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), data=vitd2.4059.female, na.action=na.omit)
  b.x<-summary(vitd.x.female)$coef[2,1]
  se.x<-summary(vitd.x.female)$coef[2,2]
  p.x<-summary(vitd.x.female)$coef[2,4]
  lo95ci<- exp(b.x-1.96*se.x)
  hi95ci<- exp(b.x+1.96*se.x)
  #nobs <- (vitd.x.female)
  all.x<-cbind(b.x,se.x,p.x, lo95ci, hi95ci)#, nobs)
  return(all.x)}
female.results<-rbind(female.4059.models(vitd2.4059.female$urxbph),female.4059.models(vitd2.4059.female$urxmhp),female.4059.models(vitd2.4059.female$urxmhh),female.4059.models(vitd2.4059.female$urxmoh),female.4059.models(vitd2.4059.female$urxecp),female.4059.models(vitd2.4059.female$dehpsum),
                      female.4059.models(vitd2.4059.female$urxmbp),female.4059.models(vitd2.4059.female$urxmib),female.4059.models(vitd2.4059.female$urxmep),female.4059.models(vitd2.4059.female$urxmzp),female.4059.models(vitd2.4059.female$urxmc1),female.4059.models(vitd2.4059.female$urxcnp), 
                      female.4059.models(vitd2.4059.female$urxcop))
row.names(female.results)<-c("BPA","MEHP","MEHHP","MEOHP","MECPP","DEHPsum","MBP","MiBP","MEP","MBzP","MCPP","MCNP","MCOP")
female.results
write.csv(female.results, file="female_4059.csv")

female.60.models=function(x){
  vitd.x.female<-lm(vitamin_d~log(x)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), data=vitd2.60.female, na.action=na.omit)
  b.x<-summary(vitd.x.female)$coef[2,1]
  se.x<-summary(vitd.x.female)$coef[2,2]
  p.x<-summary(vitd.x.female)$coef[2,4]
  lo95ci<- exp(b.x-1.96*se.x)
  hi95ci<- exp(b.x+1.96*se.x)
  #nobs <- (vitd.x.female)
  all.x<-cbind(b.x,se.x,p.x, lo95ci, hi95ci)#, nobs)
  return(all.x)}
female.results<-rbind(female.60.models(vitd2.60.female$urxbph),female.60.models(vitd2.60.female$urxmhp),female.60.models(vitd2.60.female$urxmhh),female.60.models(vitd2.60.female$urxmoh),female.60.models(vitd2.60.female$urxecp),female.60.models(vitd2.60.female$dehpsum),
                      female.60.models(vitd2.60.female$urxmbp),female.60.models(vitd2.60.female$urxmib),female.60.models(vitd2.60.female$urxmep),female.60.models(vitd2.60.female$urxmzp),female.60.models(vitd2.60.female$urxmc1),female.60.models(vitd2.60.female$urxcnp), 
                      female.60.models(vitd2.60.female$urxcop))
row.names(female.results)<-c("BPA","MEHP","MEHHP","MEOHP","MECPP","DEHPsum","MBP","MiBP","MEP","MBzP","MCPP","MCNP","MCOP")
female.results
write.csv(female.results, file="female_60.csv")

#### Males Only
male.2039.models=function(x){
  vitd.x.male<-lm(vitamin_d~log(x)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), data=vitd2.2039.male, na.action=na.omit)
  b.x<-summary(vitd.x.male)$coef[2,1]
  se.x<-summary(vitd.x.male)$coef[2,2]
  p.x<-summary(vitd.x.male)$coef[2,4]
  lo95ci<- exp(b.x-1.96*se.x)
  hi95ci<- exp(b.x+1.96*se.x)
  #nobs <- (vitd.x.male)
  all.x<-cbind(b.x,se.x,p.x, lo95ci, hi95ci)#, nobs)
  return(all.x)}
male.results<-rbind(male.2039.models(vitd2.2039.male$urxbph),male.2039.models(vitd2.2039.male$urxmhp),male.2039.models(vitd2.2039.male$urxmhh),male.2039.models(vitd2.2039.male$urxmoh),male.2039.models(vitd2.2039.male$urxecp),male.2039.models(vitd2.2039.male$dehpsum),
                      male.2039.models(vitd2.2039.male$urxmbp),male.2039.models(vitd2.2039.male$urxmib),male.2039.models(vitd2.2039.male$urxmep),male.2039.models(vitd2.2039.male$urxmzp),male.2039.models(vitd2.2039.male$urxmc1),male.2039.models(vitd2.2039.male$urxcnp), 
                      male.2039.models(vitd2.2039.male$urxcop))
row.names(male.results)<-c("BPA","MEHP","MEHHP","MEOHP","MECPP","DEHPsum","MBP","MiBP","MEP","MBzP","MCPP","MCNP","MCOP")
male.results
write.csv(male.results, file="male_2039.csv")

male.4059.models=function(x){
  vitd.x.male<-lm(vitamin_d~log(x)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), data=vitd2.4059.male, na.action=na.omit)
  b.x<-summary(vitd.x.male)$coef[2,1]
  se.x<-summary(vitd.x.male)$coef[2,2]
  p.x<-summary(vitd.x.male)$coef[2,4]
  lo95ci<- exp(b.x-1.96*se.x)
  hi95ci<- exp(b.x+1.96*se.x)
  #nobs <- (vitd.x.male)
  all.x<-cbind(b.x,se.x,p.x, lo95ci, hi95ci)#, nobs)
  return(all.x)}
male.results<-rbind(male.4059.models(vitd2.4059.male$urxbph),male.4059.models(vitd2.4059.male$urxmhp),male.4059.models(vitd2.4059.male$urxmhh),male.4059.models(vitd2.4059.male$urxmoh),male.4059.models(vitd2.4059.male$urxecp),male.4059.models(vitd2.4059.male$dehpsum),
                      male.4059.models(vitd2.4059.male$urxmbp),male.4059.models(vitd2.4059.male$urxmib),male.4059.models(vitd2.4059.male$urxmep),male.4059.models(vitd2.4059.male$urxmzp),male.4059.models(vitd2.4059.male$urxmc1),male.4059.models(vitd2.4059.male$urxcnp), 
                      male.4059.models(vitd2.4059.male$urxcop))
row.names(male.results)<-c("BPA","MEHP","MEHHP","MEOHP","MECPP","DEHPsum","MBP","MiBP","MEP","MBzP","MCPP","MCNP","MCOP")
male.results
write.csv(male.results, file="male_4059.csv")

male.60.models=function(x){
  vitd.x.male<-lm(vitamin_d~log(x)+log(urxucr)+ age + factor(race)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), data=vitd2.60.male, na.action=na.omit)
  b.x<-summary(vitd.x.male)$coef[2,1]
  se.x<-summary(vitd.x.male)$coef[2,2]
  p.x<-summary(vitd.x.male)$coef[2,4]
  lo95ci<- exp(b.x-1.96*se.x)
  hi95ci<- exp(b.x+1.96*se.x)
  #nobs <- (vitd.x.male)
  all.x<-cbind(b.x,se.x,p.x, lo95ci, hi95ci)#, nobs)
  return(all.x)}
male.results<-rbind(male.60.models(vitd2.60.male$urxbph),male.60.models(vitd2.60.male$urxmhp),male.60.models(vitd2.60.male$urxmhh),male.60.models(vitd2.60.male$urxmoh),male.60.models(vitd2.60.male$urxecp),male.60.models(vitd2.60.male$dehpsum),
                      male.60.models(vitd2.60.male$urxmbp),male.60.models(vitd2.60.male$urxmib),male.60.models(vitd2.60.male$urxmep),male.60.models(vitd2.60.male$urxmzp),male.60.models(vitd2.60.male$urxmc1),male.60.models(vitd2.60.male$urxcnp), 
                      male.60.models(vitd2.60.male$urxcop))
row.names(male.results)<-c("BPA","MEHP","MEHHP","MEOHP","MECPP","DEHPsum","MBP","MiBP","MEP","MBzP","MCPP","MCNP","MCOP")
male.results
write.csv(male.results, file="male_60.csv")

######################Subset by Race/Ethnicity###############
vitd2.white<-subset(vitd2_ph_comp, race==3)#white; N=2366
vitd2.white.dsn<-svydesign(id=~sdmvpsu, strata=~sdmvstra, weights=~weight6yr, data=vitd2.white, nest=T)

vitd2.black<-subset(vitd2_ph_comp, race==4)#black; N=887
vitd2.black.dsn<-svydesign(id=~sdmvpsu, strata=~sdmvstra, weights=~weight6yr, data=vitd2.black, nest=T)

vitd2.mexican<-subset(vitd2_ph_comp, race==1)#MA; N=879
vitd2.mexican.dsn<-svydesign(id=~sdmvpsu, strata=~sdmvstra, weights=~weight6yr, data=vitd2.mexican, nest=T)

################## Stratified by Race/Ethnicity#####################

####WHITE

white<-svyglm(vitamin_d~log(urxbph)+log(urxucr)+ age + factor (sex)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.white.dsn, na.action=na.omit)
summary(black)
confint(white, level=0.95)

white<-svyglm(vitamin_d~log(urxmhp)+log(urxucr)+ age + factor(sex)+ pir + bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.white.dsn, na.action=na.omit)
summary(white)
confint(white, level=0.95)

white<-svyglm(vitamin_d~log(urxmhh)+log(urxucr)+ age + factor(sex)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.white.dsn, na.action=na.omit)
summary(white)
confint(white, level=0.95)

white<-svyglm(vitamin_d~log(urxmoh)+log(urxucr)+ age + factor(sex)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.white.dsn, na.action=na.omit)
summary(white)
confint(white, level=0.95)

white<-svyglm(vitamin_d~log(urxecp)+log(urxucr)+ age + factor(sex)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.white.dsn, na.action=na.omit)
summary(white)
confint(white, level=0.95)

white<-svyglm(vitamin_d~log(dehpsum)+log(urxucr)+ age + factor(sex)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.white.dsn, na.action=na.omit)
summary(white)
confint(white, level=0.95)

white<-svyglm(vitamin_d~log(urxmbp)+log(urxucr)+ age + factor(sex)+ pir+bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.white.dsn, na.action=na.omit)
summary(white)
confint(white, level=0.95)

white<-svyglm(vitamin_d~log(urxmib)+log(urxucr)+ age + factor(sex)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.white.dsn, na.action=na.omit)
summary(white)
confint(white, level=0.95)

white<-svyglm(vitamin_d~log(urxmep)+log(urxucr)+ age + factor(sex)+ bmi+ pir+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.white.dsn, na.action=na.omit)
summary(white)
confint(white, level=0.95)

white<-svyglm(vitamin_d~log(urxmzp)+log(urxucr)+ age + factor(sex)+ bmi+ pir +log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.white.dsn, na.action=na.omit)
summary(white)
confint(white, level=0.95)

white<-svyglm(vitamin_d~log(urxmc1)+log(urxucr)+ age + factor(sex)+ bmi+ pir + log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.white.dsn, na.action=na.omit)
summary(white)
confint(white, level=0.95)

white<-svyglm(vitamin_d~log(urxcnp)+log(urxucr)+ age + factor(sex)+ bmi+ pir + log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.white.dsn, na.action=na.omit)
summary(white)
confint(white, level=0.95)

white<-svyglm(vitamin_d~log(urxcop)+log(urxucr)+ age + factor(sex)+ bmi+ pir + log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.white.dsn, na.action=na.omit)
summary(white)
confint(white, level=0.95)

####BLACK

black<-svyglm(vitamin_d~log(urxbph)+log(urxucr)+ age + factor(sex)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.black.dsn, na.action=na.omit)
summary(black)
confint(black, level=0.95)

black<-svyglm(vitamin_d~log(urxmhp)+log(urxucr)+ age + factor(sex)+ pir + bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.black.dsn, na.action=na.omit)
summary(black)
confint(black, level=0.95)

black<-svyglm(vitamin_d~log(urxmhh)+log(urxucr)+ age + factor(sex)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.black.dsn, na.action=na.omit)
summary(black)
confint(black, level=0.95)

black<-svyglm(vitamin_d~log(urxmoh)+log(urxucr)+ age + factor(sex)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.black.dsn, na.action=na.omit)
summary(black)
confint(black, level=0.95)

black<-svyglm(vitamin_d~log(urxecp)+log(urxucr)+ age + factor(sex)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.black.dsn, na.action=na.omit)
summary(black)
confint(black, level=0.95)

black<-svyglm(vitamin_d~log(dehpsum)+log(urxucr)+ age + factor(sex)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.black.dsn, na.action=na.omit)
summary(black)
confint(black, level=0.95)

black<-svyglm(vitamin_d~log(urxmbp)+log(urxucr)+ age + factor(sex)+ pir+bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.black.dsn, na.action=na.omit)
summary(black)
confint(black, level=0.95)

black<-svyglm(vitamin_d~log(urxmib)+log(urxucr)+ age + factor(sex)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.black.dsn, na.action=na.omit)
summary(black)
confint(black, level=0.95)

black<-svyglm(vitamin_d~log(urxmep)+log(urxucr)+ age + factor(sex)+ bmi+ pir+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.black.dsn, na.action=na.omit)
summary(black)
confint(black, level=0.95)

black<-svyglm(vitamin_d~log(urxmzp)+log(urxucr)+ age + factor(sex)+ bmi+ pir +log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.black.dsn, na.action=na.omit)
summary(black)
confint(black, level=0.95)

black<-svyglm(vitamin_d~log(urxmc1)+log(urxucr)+ age + factor(sex)+ bmi+ pir + log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.black.dsn, na.action=na.omit)
summary(black)
confint(black, level=0.95)

black<-svyglm(vitamin_d~log(urxcnp)+log(urxucr)+ age + factor(sex)+ bmi+ pir + log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.black.dsn, na.action=na.omit)
summary(black)
confint(black, level=0.95)

black<-svyglm(vitamin_d~log(urxcop)+log(urxucr)+ age + factor(sex)+ bmi+ pir + log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.black.dsn, na.action=na.omit)
summary(black)
confint(black, level=0.95)

####Mexican

mexican<-svyglm(vitamin_d~log(urxbph)+log(urxucr)+ age + factor(sex)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.mexican.dsn, na.action=na.omit)
summary(mexican)
confint(mexican, level=0.95)

mexican<-svyglm(vitamin_d~log(urxmhp)+log(urxucr)+ age + factor(sex)+ pir + bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.mexican.dsn, na.action=na.omit)
summary(mexican)
confint(mexican, level=0.95)

mexican<-svyglm(vitamin_d~log(urxmhh)+log(urxucr)+ age + factor(sex)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.mexican.dsn, na.action=na.omit)
summary(mexican)
confint(mexican, level=0.95)

mexican<-svyglm(vitamin_d~log(urxmoh)+log(urxucr)+ age + factor(sex)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.mexican.dsn, na.action=na.omit)
summary(mexican)
confint(mexican, level=0.95)

mexican<-svyglm(vitamin_d~log(urxecp)+log(urxucr)+ age + factor(sex)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.mexican.dsn, na.action=na.omit)
summary(mexican)
confint(mexican, level=0.95)

mexican<-svyglm(vitamin_d~log(dehpsum)+log(urxucr)+ age + factor(sex)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.mexican.dsn, na.action=na.omit)
summary(mexican)
confint(mexican, level=0.95)

mexican<-svyglm(vitamin_d~log(urxmbp)+log(urxucr)+ age + factor(sex)+ pir+bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.mexican.dsn, na.action=na.omit)
summary(mexican)
confint(mexican, level=0.95)

mexican<-svyglm(vitamin_d~log(urxmib)+log(urxucr)+ age + factor(sex)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.mexican.dsn, na.action=na.omit)
summary(mexican)
confint(mexican, level=0.95)

mexican<-svyglm(vitamin_d~log(urxmep)+log(urxucr)+ age + factor(sex)+ bmi+ pir+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.mexican.dsn, na.action=na.omit)
summary(mexican)
confint(mexican, level=0.95)

mexican<-svyglm(vitamin_d~log(urxmzp)+log(urxucr)+ age + factor(sex)+ bmi+ pir +log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.mexican.dsn, na.action=na.omit)
summary(mexican)
confint(mexican, level=0.95)

mexican<-svyglm(vitamin_d~log(urxmc1)+log(urxucr)+ age + factor(sex)+ bmi+ pir + log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.mexican.dsn, na.action=na.omit)
summary(mexican)
confint(mexican, level=0.95)

mexican<-svyglm(vitamin_d~log(urxcnp)+log(urxucr)+ age + factor(sex)+ bmi+ pir + log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.mexican.dsn, na.action=na.omit)
summary(mexican)
confint(mexican, level=0.95)

mexican<-svyglm(vitamin_d~log(urxcop)+log(urxucr)+ age + factor(sex)+ bmi+ pir + log(cot) + factor(season) + factor(vitdsup) + factor(cycle), vitd2.mexican.dsn, na.action=na.omit)
summary(mexican)
confint(mexican, level=0.95)

############Unweighted analyses##########
white.models=function(x){
  vitd.x.white<-lm(vitamin_d~log(x)+log(urxucr)+ age + factor(sex)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), data=vitd2.white, na.action=na.omit)
  b.x<-summary(vitd.x.white)$coef[2,1]
  se.x<-summary(vitd.x.white)$coef[2,2]
  p.x<-summary(vitd.x.white)$coef[2,4]
  lo95ci<- exp(b.x-1.96*se.x)
  hi95ci<- exp(b.x+1.96*se.x)
  all.x<-cbind(b.x,se.x,p.x, lo95ci, hi95ci)#, nobs)
  return(all.x)}
white.results<-rbind(white.models(vitd2.white$urxbph),white.models(vitd2.white$urxmhp),white.models(vitd2.white$urxmhh),white.models(vitd2.white$urxmoh),white.models(vitd2.white$urxecp),white.models(vitd2.white$dehpsum),
                    white.models(vitd2.white$urxmbp),white.models(vitd2.white$urxmib),white.models(vitd2.white$urxmep),white.models(vitd2.white$urxmzp),white.models(vitd2.white$urxmc1),white.models(vitd2.white$urxcnp), 
                    white.models(vitd2.white$urxcop))
row.names(white.results)<-c("BPA","MEHP","MEHHP","MEOHP","MECPP","DEHPsum","MBP","MiBP","MEP","MBzP","MCPP","MCNP","MCOP")
white.results
write.csv(white.results, file="white.csv")

black.models=function(x){
  vitd.x.black<-lm(vitamin_d~log(x)+log(urxucr)+ age + factor(sex)+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), data=vitd2.black, na.action=na.omit)
  b.x<-summary(vitd.x.black)$coef[2,1]
  se.x<-summary(vitd.x.black)$coef[2,2]
  p.x<-summary(vitd.x.black)$coef[2,4]
  lo95ci<- exp(b.x-1.96*se.x)
  hi95ci<- exp(b.x+1.96*se.x)
  all.x<-cbind(b.x,se.x,p.x, lo95ci, hi95ci)#, nobs)
  return(all.x)}
black.results<-rbind(black.models(vitd2.black$urxbph),black.models(vitd2.black$urxmhp),black.models(vitd2.black$urxmhh),black.models(vitd2.black$urxmoh),black.models(vitd2.black$urxecp),black.models(vitd2.black$dehpsum),
                     black.models(vitd2.black$urxmbp),black.models(vitd2.black$urxmib),black.models(vitd2.black$urxmep),black.models(vitd2.black$urxmzp),black.models(vitd2.black$urxmc1),black.models(vitd2.black$urxcnp), 
                     black.models(vitd2.black$urxcop))
row.names(black.results)<-c("BPA","MEHP","MEHHP","MEOHP","MECPP","DEHPsum","MBP","MiBP","MEP","MBzP","MCPP","MCNP","MCOP")
black.results
write.csv(black.results, file="black.csv")

mexican.models=function(x){
  vitd.x.mexican<-lm(vitamin_d~log(x)+log(urxucr)+  age + factor(sex) + bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), data=vitd2.mexican, na.action=na.omit)
  b.x<-summary(vitd.x.mexican)$coef[2,1]
  se.x<-summary(vitd.x.mexican)$coef[2,2]
  p.x<-summary(vitd.x.mexican)$coef[2,4]
  lo95ci<- exp(b.x-1.96*se.x)
  hi95ci<- exp(b.x+1.96*se.x)
  all.x<-cbind(b.x,se.x,p.x, lo95ci, hi95ci)#, nobs)
  return(all.x)}
mexican.results<-rbind(mexican.models(vitd2.mexican$urxbph),mexican.models(vitd2.mexican$urxmhp),mexican.models(vitd2.mexican$urxmhh),mexican.models(vitd2.mexican$urxmoh),mexican.models(vitd2.mexican$urxecp),mexican.models(vitd2.mexican$dehpsum),
                     mexican.models(vitd2.mexican$urxmbp),mexican.models(vitd2.mexican$urxmib),mexican.models(vitd2.mexican$urxmep),mexican.models(vitd2.mexican$urxmzp),mexican.models(vitd2.mexican$urxmc1),mexican.models(vitd2.mexican$urxcnp), 
                     mexican.models(vitd2.mexican$urxcop))
row.names(mexican.results)<-c("BPA","MEHP","MEHHP","MEOHP","MECPP","DEHPsum","MBP","MiBP","MEP","MBzP","MCPP","MCNP","MCOP")
mexican.results
write.csv(mexican.results, file="mexican.csv")

mexican.models=function(x){
  vitd.x.mexican<-lm(vitamin_d~log(x)+log(urxucr)+ age + factor(sex)+ pir+ bmi+ log(cot) + factor(season) + factor(vitdsup) + factor(cycle), data=vitd2.mexican, na.action=na.omit)
  b.x<-summary(vitd.x.mexican)$coef[2,1]
  se.x<-summary(vitd.x.mexican)$coef[2,2]
  p.x<-summary(vitd.x.mexican)$coef[2,4]
  lo95ci<- exp(b.x-1.96*se.x)
  hi95ci<- exp(b.x+1.96*se.x)
  all.x<-cbind(b.x,se.x,p.x, lo95ci, hi95ci)#, nobs)
  return(all.x)}
mexican.results<-rbind(mexican.models(vitd2.mexican$urxbph),mexican.models(vitd2.mexican$urxmhp),mexican.models(vitd2.mexican$urxmhh),mexican.models(vitd2.mexican$urxmoh),mexican.models(vitd2.mexican$urxecp),mexican.models(vitd2.mexican$dehpsum),
                       mexican.models(vitd2.mexican$urxmbp),mexican.models(vitd2.mexican$urxmib),mexican.models(vitd2.mexican$urxmep),mexican.models(vitd2.mexican$urxmzp),mexican.models(vitd2.mexican$urxmc1),mexican.models(vitd2.mexican$urxcnp), 
                       mexican.models(vitd2.mexican$urxcop))
row.names(mexican.results)<-c("BPA","MEHP","MEHHP","MEOHP","MECPP","DEHPsum","MBP","MiBP","MEP","MBzP","MCPP","MCNP","MCOP")
mexican.results
write.csv(mexican.results, file="mexican.csv")


###########################

women.bpa<-lm(vitamin_d~log(bpa)+ log(urxucr)+ age + factor(race)+ bmi+ log(cot) + factor(season) + factor(suppl) + factor(cycle), data=vitd2.2039.male, na.action=na.omit)
b.x<-summary(women.bpa)$coef[2,1]
se.x<-summary(women.bpa)$coef[2,2]
p.x<-summary(women.bpa)$coef[2,4]
lo95ci<- (b.x-1.96*se.x)
hi95ci<- (b.x+1.96*se.x)
summary(women.bpa)
b.x
se.x
p.x
lo95ci
hi95ci

summary(lm(vitamin_d~log(bpa)+ log(urxucr)+ age + factor(race)+ bmi+ log(cot) + factor(season) + factor(suppl) + factor(cycle), data=vitd2.4059.female, na.action=na.omit))
summary(lm(vitamin_d~log(bpa)+ log(urxucr)+ age + factor(race)+ bmi+ log(cot) + factor(season) + factor(suppl) + factor(cycle), data=vitd2.60.female, na.action=na.omit))

###race, vit d supplement, pir (continuous), bmi (continuous), milk, gender, age
###check stratified by race, gender, age groups

###distributions
dist=function(x){
  geomean.x<-exp(mean(log(x),na.rm=T))
  geosd.x<-exp(sd(log(x),na.rm=T))
  q1.x<-quantile(x,probs=c(0.25),names=FALSE,na.rm=T)
  q2.x<-quantile(x,probs=c(0.5),names=FALSE,na.rm=T)
  q3.x<-quantile(x,probs=c(0.75),names=FALSE,na.rm=T)
  ninety.x<-quantile(x,probs=c(0.9),names=FALSE,na.rm=T)
  ninetyfive.x<-quantile(x,probs=c(0.95),names=FALSE,na.rm=T)
  max.x<-quantile(x,probs=c(1),names=FALSE,na.rm=T)
  all.x<-cbind(geomean.x,geosd.x,q1.x,q2.x,q3.x,ninety.x,ninetyfive.x,max.x)
  return(all.x)}
phthalate.distributions<-rbind(dist(vitd4$urxmbp),dist(vitd4$urxmhp),dist(vitd4$urxmc1))
phthalate.distributions
write.csv(phthalate.distributions,file="output.csv")

vitd<-read.xport('data0106.xport') #load new R file created from SAS
names(vitd)<-tolower(names(vitd)) #change variables to lower case
save(vitd,file='vitd.rda') #save as permanent dataset
load("vitd.rda") #load dataset

setwd("C:Users/laujohns/Desktop")





















































































































































































































































































































































































































































































































































































































































































































































































































































