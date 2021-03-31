library(tidyverse)
library(Hmisc)
library(magrittr)
require(geepack)
require(xtable)
require(doBy)
require(readr)
library(lubridate)
library(matrixStats)
# library(plotly)
library(icd.data)
library(cowplot)
library(emmeans)
library(splines)
library(mitools)
library(Epi)

## Function for data derivations
data.prepare <- function(data,option = 1,emigration.rate=0,seed = 1234) {
  # This function aims to clean the patients outcomes
  # variables: 
  # # data: is the measles data
  # # option = 0 if  we allow all the discharged with presumably worse condition alive
  # # option = 1 if  we allow all the discharged with presumably worse condition dead
  # # option = 2 if  we generate the emigration patient out of the HCMC city
  
  #### ------------------------------ Rule to derive dataset ------------------------------ ####
  # 1. Only the infection cases will be included in the analysis
  # 2. Impute the birth day for who doesn't have birthday record but have the birth years. 
  #     We immputed it as 30June of the birth years. Who doesn't have birth year 
  #     will be removed from the analysis.
  # 3. We only consider a new admission due to a new episode of infection if the 
  #     new admission date after the last discharge at least 2 days, otherwise we
  #     merge them into one espisode of new infections
  # 4. We assume that for each episode of admission is due to one main infection type
  # 5. All measles second infections which happened after 21 days will be excluded 
  #     from the analysis because later we dont know how identify the before and after events
  # 6. All the second measles admission within 21 days after the first admission
  #     will be considered as the measles complicated
  
  
  ## Define ICD.10 code
  icd.infection.diarrhea    <- paste(rep("A",10),c("00","01","02","03","04","05","06","07","08","09"),sep = "")
  icd.infection.respiratory <- paste(rep("J",23),c("00","01","02","03","04","05","06","07","08","09",10:22),sep = "")
  icd.infection.CNS         <- c(paste(rep("A",10),80:89,sep =""),paste(rep("G",9),c("00","01","02","03","04","05","06","07","09"),sep = ""))
  icd.infection             <- c(c(paste(rep("A",100),c("00","01","02","03","04","05","06","07","08","09",10:99),sep = ""),
                               paste(rep("B",100),c("00","01","02","03","04","05","06","07","08","09",10:99),sep = "")),
                               paste(rep("J",23),c("00","01","02","03","04","05","06","07","08","09",10:22),sep = ""),
                               paste(rep("G",10),c("00","01","02","03","04","05","06","07","08","09"),sep = ""))
          
## Rename the id to patid
  data %<>% rename(patid=iduse)

## Remove all duplicated cases with the same admission and discharge. 
## We will remove the first records. Among them, 
## the second admissions are measles therefore we didn't remove the measles case.
  data %<>% mutate(new.id = paste(patid,date.admit,date.dis,sep = "-"))
  index.new.id <- as.numeric(rownames(data)[duplicated(data$new.id)])
  data %<>% slice(-(index.new.id - 1))
  
## Impute all the missing of reasons of discharge --> "unknown".
  data %<>% mutate(outcome = as.character(outcome)) %>% 
    mutate(outcome = as.factor(ifelse(!is.na(outcome),outcome,"unknown")))

##  Merge all the cases have the same admission but different discharge
  ### Arrange the data and create a new id
  data %<>% arrange(patid,date.admit,date.dis) %>% 
            mutate(id.date = paste(patid,date.admit,date.dis,sep = "-")) 
  
  ### Create a new dataset includes the latest discharge of each individual and then create a new id
  data.new <- data %>% 
    group_by(patid,date.admit) %>%
                    dplyr::summarise(date.dis.max=max(date.dis),.groups="keep") %>% 
                    dplyr::mutate(id.date=paste(patid,date.admit,date.dis.max,sep = "-"),
                           select.id=1) %>%
                    ungroup() %>% 
                    select(id.date,select.id)
  
  ### Merge two datasets: data and data.new and create an indicator for selecting the latest discharge
  
  data %<>% left_join(data.new,by="id.date") %>% 
            arrange(patid,date.admit,date.dis) %>% 
            mutate(select.id=ifelse(is.na(select.id),0,select.id)) 
  
  ### Make sure that for all the patients when we merge two identical admissions, we didn't miss the measles cases i.e. icd.10 code ="B05"
  data %<>% mutate(icd.m3=ifelse((lag(select.id,1) == 0  &
                                         lag(patid) == patid) &
                                        (icd.m3 == "B05" | lag(icd.m3,1) == "B05"),"B05",icd.m3))
  data$icd.m3[1] <- "J02"

  data %<>%  filter(select.id == 1)
  data[,c("id.date","select.id")] <- NULL
  
## Re-defined the outcome presumably.worse. If the patient re-admitted to the hospital after 
## being discharged with presumably.worse, the outcome presumably.worse were considered as discharged.
  
  ### Take a subset of patients with outcome presumably.worse (dat),
  ### and identify the last hospital admission of these patients
  
  dat <- data %>% filter(outcome == "presumably.worse")
  id_worse <- unique(dat$patid)

  data %<>% mutate(flag=1) %>% group_by(patid) %>% 
            dplyr::mutate(visit=cumsum(flag),visit.last = max(visit)) %>% ungroup() 
  dat <- data %>% filter(patid %in% id_worse) %>% arrange(patid, visit)
  
  ### if patient re-admited to the hospial after the last visit with presumably.worse outcome then 
  ### they supposed to be discharged and recovery in the previous admission
  dat <- data %>% filter(patid %in% id_worse) %>% arrange(patid, visit) %>% 
                   mutate(presumably_worse = ifelse(outcome == "presumably.worse" & visit == visit.last,"worse.last","discharged"),
                          index.new = paste(patid,date.admit,sep = ":"))
  data %<>% mutate(index.new = paste(patid,date.admit,sep = ":"))
  ### Merge two datset dat and data.
  data %<>%  left_join(subset(dat,select = c("index.new","presumably_worse")),by="index.new") %>% 
             mutate(outcome = ifelse(outcome == "presumably.worse",presumably_worse,as.character(outcome)))
  data %<>%  mutate(outcome = ifelse(outcome == "worse.last","presumably.worse",as.character(outcome)))
  
  if (option >= 1) {
    ## we set all the discharged with presumably worse condition dead
       data %<>% mutate(outcome=ifelse(outcome == "presumably.worse","dead",outcome))
  }
  
  ## Impute the missing birthday information##
  ## and exclude patients without the birth year
  data.imputed.birth <- data %>% filter((is.na(dob)) & (!is.na(birthyear))) %>% 
                                mutate(dob = as.Date(paste("30June",birthyear,sep = ""),"%d%b%Y",origin = "1970-01-01"))
  data <- rbind(subset(data,!is.na(data$dob)),data.imputed.birth)
  
  
  ### Create the date of death: date_death
  data %<>% mutate(date_death=as.Date(ifelse(outcome == "dead",date.dis,NA),origin = "1970-01-01")) 
  
  # In this analysis, we set the time to follow-up is the earliest of 31 December 2015
  # which is the cut-off date of the dataset and the 15th birthday of the child 
  data %<>% mutate(time.at.15birthday = as.Date(paste(birthyear + 15,month(dob),day(dob),sep = "-"),origin = "1970-01-01"),
                   time.at.lastfup = pmin(time.at.15birthday,as.Date("2015-12-31",origin = "1970-01-01")),
                   time.at.lastfup = as.Date(ifelse(outcome == "dead",data$date_death,time.at.lastfup),origin = "1970-01-01"),
                   time.dis.to.lastfup = as.numeric(difftime(time.at.lastfup,date.dis,units = "days"))) 
  
  ### Remove all patients who admitted after the time follow-up
  ### create a variable measles, and infection
  data  %<>% filter(time.dis.to.lastfup >= 0) %>% 
             arrange(patid,date.admit,date.dis) %>% 
             select("patid","sex","district","dob","birthyear","date.admit","date.dis","icd.m3","outcome","year","diagnosis.en","time.at.lastfup","time.dis.to.lastfup","visit") %>% 
             mutate(measles = ifelse(icd.m3 == "B05","yes","no"),
             infection = ifelse(icd.m3 %in% icd.infection,"yes","no"))
  
  ## Select only records with infection and patient who admited due to MeV infection
  data  %<>%  filter(infection == "yes") %>% 
              mutate(age.at.admit.d = as.numeric(difftime(date.admit,dob,units = "days"))) %>% 
              filter(patid %in% unique(filter(.,measles == "yes")$patid))

  ## Select only records from MeV infection and then determine which is the first measles and  second measles infection.
  dat.tmp  <- data %>%  filter(measles == "yes") %>% 
                        arrange(patid,date.admit,date.dis) %>% 
                        mutate(measles.order = ifelse(duplicated(patid),"second","first"))
  
  ### Select non-measles records
  data.without.measles <- data %>% filter(measles == "no") %>%
                                   mutate(measles.order="no") %>%
                                   filter(patid %in% c(unique(dat.tmp$patid)))
  data <- rbind(data.without.measles,dat.tmp)
  data <- arrange(data,patid,date.admit,date.dis)
  
  ##  Merge all the cases discharge and later admit within 2 days. If one of the episode
  ##  is measles then the icd.m3 of merged case will be measles otherwise it will 
  ##  be the icd of the latter case (which is not important to this analysis). 
  data %<>%   arrange(patid,date.admit,date.dis) %>% 
              mutate(lag.patid  = lag(patid,1),
                     lag.admit  = lag(date.admit,1),
                     lag.icd.m3 = lag(icd.m3,1)) 

  data %<>% mutate(id.measles.bf = ifelse((as.numeric(difftime(date.admit,lag(date.dis,1),unit = "day")) <= 2) &
                                          (lag.patid == patid) & (icd.m3 != "B05" & lag.icd.m3 == "B05"),1,0),
                   id.measles.af = ifelse((as.numeric(difftime(date.admit,lag(date.dis,1),unit = "day")) <= 2) &
                                          (lag.patid == patid) & (icd.m3 == "B05" & lag.icd.m3 != "B05"),1,0),
                   id.measles.both = ifelse((as.numeric(difftime(date.admit,lag(date.dis,1),unit = "day")) <= 2) &
                                          (lag.patid == patid) & (icd.m3 == "B05" & lag.icd.m3 == "B05"),1,0),
                   id.non.measles = ifelse((as.numeric(difftime(date.admit,lag(date.dis,1),unit = "day")) <= 2) &
                                          (lag.patid == patid) & (icd.m3 != "B05" & lag.icd.m3 != "B05"),1,0)) %>%
            arrange(patid,desc(date.admit)) %>% group_by(patid) %>% 
            dplyr::mutate(id.measles.bf = lag(id.measles.bf,1),
                   id.measles.af = lag(id.measles.af,1),
                   id.measles.both = lag(id.measles.both,1),
                   id.non.measles = lag(id.non.measles,1)) %>% 
    dplyr::mutate (id.measles.bf= ifelse(is.na(id.measles.bf),0,id.measles.bf),  ## Impute value due to the lag function -->NA value
                    id.measles.af= ifelse(is.na(id.measles.af),0,id.measles.af),
                    id.measles.both= ifelse(is.na(id.measles.both),0,id.measles.both),
                    id.non.measles= ifelse(is.na(id.non.measles),0,id.non.measles)) %>% 
            arrange(patid,date.admit,date.dis) %>% ungroup()
  
  ### For the special cases with consecutive admissions after the pairs of measles-before and before the pairs of measles-after will not be merged
  data %<>% mutate(id.non.measles=ifelse((id.non.measles==1 & (id.non.measles==lag(id.measles.bf,1)|id.non.measles==lead(id.measles.af,1))),0,id.non.measles))

  ## For MEV occured before other admission, we merged the second admission into the measles admission including merge the date admission, ICD.m3="B05". 
  ## For MEV occured after other admission, we merged the measles admission into the first admission including merge the date admission, ICD.m3="B05".
  ## For both admissions were due to MeV, we do the same thing.
  ## For both cases were not MeV admission as well.
  ## The correct diagnosis of the admission is of the second admission
  data %<>% mutate(date.admit=as.Date(ifelse(lag(id.measles.bf,1)==1 & id.measles.bf==0|
                                               lag(id.measles.af,1)==1 & id.measles.af==0|
                                               lag(id.measles.both,1)==1 & id.measles.both==0|
                                               lag(id.non.measles,1)==1 & id.non.measles==0,lag.admit,date.admit),origin = "1970-01-01"),
                   icd.m3 = ifelse(lag(id.measles.bf,1)==1 & id.measles.bf==0|
                                     lag(id.measles.af,1)==1 & id.measles.af==0|
                                     lag(id.measles.both,1)==1 & id.measles.both==0,"B05",icd.m3),
                   measles.order=ifelse(lag(id.measles.bf,1)==1 & id.measles.bf==0|
                                          lag(id.measles.af,1)==1 & id.measles.af==0|
                                          lag(id.measles.both,1)==1 & id.measles.both==0,"first",measles.order))
  data$date.admit[1] <- as.Date("2011-09-21",origin = "1970-01-01")
  
  ## Select only admisions with the correct diagnosis and re-update the age.at.admit.d after merging some admissions.
  data %<>% filter(id.measles.bf==0 & id.measles.af==0 & id.measles.both ==0 & id.non.measles==0) %>%  mutate(age.at.admit.d = as.numeric(difftime(date.admit,dob,units = "days")))
  
  ## Classify the type of infection. for the second measles admission will be considered as measles complicated.
  data <- data %>%
    mutate(
      infection.type = case_when(
        icd.m3 == 'B05' ~ 'Measles',
        icd.m3%in%icd.infection.diarrhea ~ 'Diarrhea',
        icd.m3%in% icd.infection.respiratory ~ 'Respiratory Infection',
        icd.m3%in%icd.infection.CNS ~"CNS Infection",
        TRUE ~ 'Other infection'
      )
    ) %>% 
    mutate(
      infection.type   = ifelse(measles.order == "second","Measles complicated",infection.type)
    )
  
  ## Define the event time of  measles admission
  data.tmp.measles      <- data %>% filter(measles.order == "first") %>% 
                                    mutate(
                                      date.admit.measles = date.admit,
                                      date.dis.measles = date.dis,
                                      age.at.measles.d = age.at.admit.d
                                    )
  
  data <- merge(data,subset(data.tmp.measles,select = c("patid","date.admit.measles","age.at.measles.d","date.dis.measles")),by = "patid",all.x = T)
  data %<>% mutate(hospitalization = as.numeric(difftime(date.dis,date.admit,units = "days")),
                   time.ad.measles.to.ad = as.numeric(difftime(date.admit,date.admit.measles,units = "days")),
                   period = as.factor(ifelse(time.ad.measles.to.ad > 0,"after",ifelse(time.ad.measles.to.ad < 0,"before","measles"))),
                   age = age.at.admit.d/365,
                   sex = as.factor(sex),
                   outcome = as.factor(data$outcome),
                   infection.type = as.factor(data$infection.type),
                   flag=1
                   ) %>% arrange(patid,date.admit,date.dis)
  
  ## Removed all HIV cases. Becauses the immunity of these childrens were already imunocompromised,
  ## therefore we will compilcated the analysis
  HIV_pts <- subset(data,icd.m3=="B20")
  data <- subset(data,!patid %in% unique(HIV_pts$patid))
  
  ## Simulate the emigration population
  if (option == 2){  
    set.seed(seed)
    ### Split data into before and after MeV. Because we focus on people emigrate after MeV
    data.before <- subset(data,period == "before")
    data.after  <- subset(data,period != "before")
    ### Identify the last year which visit the hospital. Because in the later year patients may emigrate out of the city. 
    dat <- data.after %>% arrange(patid,visit)
    dat.sum <- dat %>% group_by(patid) %>% dplyr::summarise(visit.last = max(visit),last.year.ev = max(year),date.dis.last=as.Date(max(date.dis),origin = "1970-01-01"),.groups="keep") %>% ungroup() 
    dat <- merge(dat,dat.sum,by = "patid",all.x=T)
    dat.tmp <- data.frame()
    ### Given the emigration rate = 0.0083 per year, we simulate the time to emigration out of the city for each individual based on exponential distribution since the date of birth.
    data.cen<-data.frame(patid=unique(dat$patid),random.time.from.date.dis.last.to.emigration=round(rexp(length(unique(dat$patid)),emigration.rate)*365))
    dat <- left_join(dat,data.cen,by = "patid") %>% 
                    mutate(date.at.random.emigration=as.Date(date.dis.last+random.time.from.date.dis.last.to.emigration,origin = "1970-01-01"))
     for(i in 2005:2015) {
      ### We focus on patient who has the last visit is year i
      dat.sub <- subset(dat,last.year.ev == i)
      n.pt <- length(unique(dat.sub$patid))
      if(n.pt > 0) {
        dat.sub %<>% mutate(
                            censored.time = as.Date(pmin(time.at.lastfup,date.at.random.emigration),origin = "1970-01-01"),
                            censored  = ifelse(time.at.lastfup>censored.time & visit==visit.last,"yes","no")
                            ) %>% 
                     arrange(patid,date.admit,date.dis)
        dat.tmp<-rbind(dat.sub,dat.tmp)
      }
     }
    dat.tmp %<>% mutate(
                          time.at.lastfup=as.Date(pmin(time.at.lastfup,censored.time),origin = "1970-01-01"),
                          time.dis.to.lastfup=difftime(time.at.lastfup,date.dis.measles,units = "days")
                        ) 
    dat.tmp[,c("random.time.from.date.dis.last.to.emigration","date.at.random.emigration")]<-NULL
    data.before %<>% mutate(
                              visit.last = NA,
                              last.year.ev = NA,
                              date.dis.last =NA,
                              censored.time = NA,
                              censored = NA
                              
                            ) 
    data <- rbind(data.before,dat.tmp) %>% arrange(patid,date.admit,date.dis)
  }
  data[,c("measles.order","lag.patid","lag.admit","lag.icd.m3","date.dis.last","censored.time")]<-NULL
  return(data)
}

data.transformed <- function(data,option=0,emigration.rate=0,seed=1234) {
  # # This function aims to transform the patient records into a weekly time series
  # of hospital admission
  # # variables
  # # data: is the measles data
  # # option = 0 if  we allow all the discharged with presumably worse condition alive
  # # option = 1 if  we allow all the discharged with presumably worse condition dead
  # # option = 2 if  we generate the emigration patient out of the HCMC city
  # # emigration.rate: is the HCMC emigration rate
  
  
  
  ##########################################################################################
  ## Step1: Do the data derivation  
  ########################################################################################## 
  data <- data.prepare(data,option = option,emigration.rate = emigration.rate,seed=seed)
  #########################################################################################
  ## Step2: Do the data transformation
  #########################################################################################  
  #########################################################################################  
  
  ####------------------------------Transformation for  data after measles event------------------------------####
  years.before   <- 2
  time.period    <- 365*years.before
  ### Include only the admission occured 2 year before measles admission
  ### The duration from measles infection to the other admission
  dat.tmp.before <- data %>% filter(period %in% c("before","measles")) %>%
                             arrange(patid,date.admit,date.dis) %>% 
                             mutate(time.dis.to.ad.measles = as.numeric(difftime(date.admit.measles,date.dis,units = "days")),
                                    time.ad.to.ad.measles = as.numeric(difftime(date.admit.measles,date.admit,units = "days")))
  # test if time to measles have na value
  # sum(is.na(dat.tmp.before$time.dis.to.ad.measles))
  
  ## Create dataset of 2 year before MeV event
  ## Impute the follow-up time between the last event to MeV event
  ## flag.index.impute indicate  measles event and there were event before measles
  dat.tmp.before.years.before <- dat.tmp.before %>% filter(time.ad.to.ad.measles<= time.period) %>% 
                                                    arrange(patid,date.admit,date.dis) %>% group_by(patid) %>% 
                                                    dplyr::mutate(status = ifelse(period == "measles",0,1),
                                                           flag.index.impute = ifelse((period=="measles")&(!is.na(lag(patid,1))),1,0))  %>% ungroup()
  ## For a subset of data don't need impute the flup
  dat.tmp.before.years.before.tmp <- dat.tmp.before.years.before %>% 
                                      filter(flag.index.impute == 0) %>% 
                                      mutate(obs.time = ifelse(period == "measles",ifelse(age.at.measles.d >= time.period,time.period,age.at.measles.d),
                                                               ifelse(age.at.measles.d >= time.period,pmax(time.period - abs(time.ad.measles.to.ad),0),
                                                                      age.at.admit.d))) %>% 
                                      arrange(patid,date.admit,date.dis) %>% group_by(patid) %>% 
                                      dplyr::mutate(obs.time = ifelse(date.admit == dob & age.at.measles.d < time.period,1,obs.time),# for special cases the child were born and admitted to hospital
                                             lag.obs.time = lag(obs.time,1),
                                             lag.patid = lag(patid,1),
                                             time.start = ifelse(!(is.na(lag.patid)),
                                                                 lag.obs.time,0),
                                             time.start = ifelse(time.start > 0,
                                                                 time.start + lag(hospitalization,1),
                                                                 time.start)
                                      ) %>% ungroup()
  
  dat.tmp.before.years.before.tmp$time.start[1] <- 0
  dat.tmp.before.years.before.tmp %<>% mutate(time.stop  = obs.time,
                                              time.stop  = ifelse(time.stop == 0,time.start+1,time.stop)) %>% 
                                       filter(time.stop > 0)

  
  #(n10<-length(unique(dat.tmp.before.years.before$patid)))
  ############################## Impute the following up time up to the last flup ############
  dat.impute <- dat.tmp.before.years.before %>% filter(flag.index.impute == 1 )
  index.impute <- unique(dat.impute$patid)
  dat.impute <- dat.tmp.before.years.before.tmp %>% filter(patid %in% index.impute) %>% 
                                                  mutate(time.start = time.stop + hospitalization) %>% 
                                                  arrange(patid,date.admit,date.dis)
  dat.impute <- as.data.frame(dat.impute)
  ### We chose the last visit before measles infection
  #index <- lastobs(dat.impute$patid)
  index <- lastobs(~patid,dat.impute)
  dat.impute %<>% slice(index) %>% 
                  mutate(time.stop = ifelse(age.at.measles.d >= time.period,time.period,age.at.measles.d),# Update time.stop
                         status = 0)

  ## Merge two dataset dat.impute and dat.tmp.before.years.before
  dat.tmp.before.years.before <- rbind(dat.tmp.before.years.before.tmp,dat.impute)
  dat.tmp.before.years.before %<>% filter(time.stop>time.start) %>% arrange(patid,time.start)
  
  (n10 <- length(unique(dat.tmp.before.years.before$patid)))
  
  ## Transforming the dataset into daily interval observations
  dat.tmp.before.split <- survSplit(Surv(time.start,time.stop, status)~.,dat.tmp.before.years.before,cut=c(seq(0,time.period,by=1)),episode ="days")
  
  dat.tmp.before.split$days <- dat.tmp.before.split$days - 1
  dat.tmp.before.split$age.at.infection <- ifelse(dat.tmp.before.split$age.at.measles.d >= time.period,
                                                  dat.tmp.before.split$days/365 + (dat.tmp.before.split$age.at.measles.d - time.period)/365,
                                                  dat.tmp.before.split$days/365)
  
  (n10 <- length(unique(dat.tmp.before.split$patid)))
  
  ####------------------------------Transformation for  data after measles event------------------------------####
  data.tmp.after <- data %>% filter(period %in% c("after","measles")) %>% arrange(patid,date.admit,date.dis)
  (n10<-length(unique(data.tmp.after$patid)))
  
   ## We remove the washout period of 2 weeks=14 days.
  removed <- TRUE
  #removed <- FALSE
  ## collapse two consecutive episode if it happened within 2 weeks
  if(removed == TRUE) {
    # update the new discharge date of measles
    data.tmp.after %<>% mutate(date.last14.measles = date.admit.measles + 14) %>% 
                        arrange(patid,date.admit,date.dis) %>% group_by(patid) %>% 
                        dplyr::mutate(selected.removed = ifelse(!is.na(lag(patid,1)) & (date.last14.measles > date.admit),1,0),
                        date.dis.measles.updated = pmax(date.last14.measles,date.dis.measles)) %>% ungroup() %>% 
                        arrange(patid,desc(date.admit)) 
    
    ### Since there is no 3 admission within 2 week since measles admission, so we defined the updated date.dis.measles
    ### Remove all case happening within 2 weeks
    data.tmp.after %<>% mutate(date.dis.measles = ifelse(period == "measles" & patid == lag(patid,1) & lag(selected.removed) == 1,
                                              pmax(lag(date.dis,1),date.dis.measles.updated),date.dis.measles.updated)) %>% 
                        filter(selected.removed==0) %>% 
                        arrange(patid,date.admit,date.dis) %>% 
                        mutate(date.dis.measles = as.Date(date.dis.measles,origin = "1970-01-01"))
    data.tmp.after$date.dis.measles[1] <- as.Date(data.tmp.after$date.admit.measles[1]+14,origin = "1970-01-01")

    data.tmp <- data.tmp.after %>% filter(period == "measles")
    data.tmp.after$date.dis.measles <- NULL
    data.tmp.after <- left_join(data.tmp.after,subset(data.tmp,select = c("patid","date.dis.measles")),by = "patid")
    ### update the new time at last follow-up for the data after adding  14 days
    data.tmp.after %<>% mutate(time.at.lastfup = as.Date(pmax(time.at.lastfup,date.dis.measles),origin = "1970-01-01"))
  }
  
  ### Define the flag.index.impute variable is the indication of MeV event without other event after 
  data.tmp.after %<>% arrange(patid,desc(date.admit)) %>% group_by(patid) %>% 
                      dplyr::mutate(flag.index.impute = ifelse(is.na(lag(patid,1)) & (period == "measles"),1,0)) %>% ungroup()
  data.tmp.after$flag.index.impute[1] <- 1
  data.tmp.after %<>% arrange(patid,date.admit,date.dis)
  
  
  ### Compute obs.time.from.MeV and set obs.time.from.MeV==-14 for measles cases
  ### Identify time to follow-up for each patient
  data.tmp.after.tmp <- data.tmp.after %>% filter(!((period == "measles") & (flag.index.impute == 0))) %>% 
    arrange(patid,date.admit,date.dis) %>% 
    mutate(obs.time.from.MeV = ifelse(period != "measles",as.numeric(difftime(date.admit,date.dis.measles,units="days")),0),# the duration from measles infection to the other admission
           time.dis.MeV.to.lastfup = as.numeric(difftime(time.at.lastfup,date.dis.measles,units="days")),
           hospitalization.measles = as.numeric(difftime(date.dis.measles,date.admit.measles,units="days"))) 
  (n10<-length(unique(data.tmp.after.tmp$patid)))
  
  data.tmp.after.tmp %<>%  arrange(patid,date.admit,date.dis) %>% 
    mutate(lag.patid  = lag(data.tmp.after.tmp$patid,1),
           obs.time.from.MeV = ifelse(obs.time.from.MeV == 0 & flag.index.impute == 0,1,obs.time.from.MeV),
           lag.obs.time.from.MeV = lag(data.tmp.after.tmp$obs.time.from.MeV,1),
           time.start = ifelse(lag.patid == patid,lag.obs.time.from.MeV+lag(hospitalization,1),0))
  data.tmp.after.tmp$time.start[1] <- 0
  (n10<-length(unique(data.tmp.after.tmp$patid)))
  data.tmp.after.tmp %<>% mutate(time.stop   = ifelse(flag.index.impute == 1,time.dis.MeV.to.lastfup,obs.time.from.MeV),
                                 status   = ifelse(flag.index.impute == 1,0,1)) 
  (n10<-length(unique(data.tmp.after.tmp$patid)))
  
  ############################## Impute the following up time up to the last flup ############
  ############################## We remove all the cases without event after MeV. 
  data.impute <- data.tmp.after %>% filter(flag.index.impute == 0)
  id.impute<-unique(data.impute$patid)
  dat.impute <- data.tmp.after.tmp %>% filter(patid%in%id.impute) %>% 
                                       mutate(time.start = time.stop + hospitalization) %>% 
                                       arrange(patid,date.admit,date.dis) 

  ### chose the last row
  index <- lastobs(dat.impute$patid)
  dat.impute  %<>%  slice(index) %>% 
                    mutate(time.stop = time.dis.MeV.to.lastfup,
                           status = 0) 
  data.tmp.after  <- rbind(data.tmp.after.tmp,dat.impute)
  ### Remove all the case without event due to dead or migration out of the city i.e. no follow-up.
  data.tmp.after  %<>%  arrange(patid,date.admit) %>% filter(time.start < time.stop)
  data.tmp.after  <- subset(data.tmp.after, select = c("patid","period","status","time.start","time.stop","sex","district","age.at.measles.d","hospitalization.measles","date.dis.measles","infection.type","hospitalization"))
  (n10 <- length(unique(data.tmp.after$patid)))
  data.tmp.after.split <- survSplit(Surv(time.start,time.stop, status)~.,data.tmp.after,cut = c(seq(0,max(data.tmp.after$time.stop),by = 1)),episode ="days")
  data.tmp.after.split %<>% mutate(days=days-1,# day 0 is the measle day and day 1 ? the next day after measles event.
                                   time.after.mv = days/365 + hospitalization.measles/365,
                                   age.at.measles.after.hosp=(age.at.measles.d + hospitalization.measles)/365) %>% 
                                   arrange(patid,time.after.mv)

  ####------------------------------Pooling both data before and after MeV to the data.total------------------------------####
  dat.before <- subset(dat.tmp.before.split,select = c("patid","status","age.at.infection","sex","district","age.at.measles.d","infection.type","hospitalization","days"))
  dat.before %<>% mutate(weeks = cut(days,breaks = seq(0,max(days) + 7,by = 7)),
                         time=1)
  dat.before <- dat.before %>% group_by(patid,sex,district,age.at.measles.d,infection.type,hospitalization,weeks) %>% dplyr::summarise(status = sum(status),age.at.infection = mean(age.at.infection),time.interval = sum(time),.groups="keep") %>% ungroup() 
  dat.before %<>% mutate(time.after.mv = -1,period="before") %>% as.data.frame()


  dat.after <- subset(data.tmp.after.split,select=c("patid","status","age.at.measles.after.hosp","time.after.mv","sex","district","age.at.measles.d","infection.type","hospitalization","days","hospitalization.measles"))
  dat.after %<>% mutate(time.after.mv = days/365 + hospitalization.measles/365,
                        age.at.infection = age.at.measles.d/365 + days/365 + hospitalization.measles/365,
                        weeks = cut(days,breaks=seq(0,max(days)+7,by=7)),
                        time = 1)
  
  dat.after <- dat.after %>% group_by(patid,sex,district,age.at.measles.d,infection.type,hospitalization,weeks) %>% dplyr::summarise(status = sum(status),age.at.infection = mean(age.at.infection),time.after.mv = mean(time.after.mv),time.interval = sum(time),.groups="keep")  %>%
    ungroup() %>% mutate(period="after")  %>% as.data.frame()

  data.total <- rbind(dat.before,dat.after)
  
  data.total$patid<-as.factor(data.total$patid)
  data.total  %<>%  mutate(patid=as.factor(patid),
                           age.at.measles.y = age.at.measles.d/365,
                           time.interval = time.interval/365,
                           offset_time = log(time.interval)) %>% 
                    arrange(patid,age.at.infection,time.after.mv)
  data.total <- subset(data.total,select=c(patid,age.at.infection,time.after.mv,offset_time,status,period,age.at.measles.y,infection.type,hospitalization,sex,district,time.interval)) %>% 
                arrange(patid,age.at.infection)
  return(data.total)
}




func2<-function(x,knots,coef,time.max=6){
  ## This function computes the value of log incidence rate ratio of hospital admission after vs. 2yrs before MeV
  ## x is the time value
  ## knots is the knot that used for the spines function in the Poisson model
  ## coef is the coeficent of the time after measles
  ## last.knots is the right boundary knot 
  
  c<-c(1,splines::ns(pmax(x-14/365,0), knots=knots-14/365, Boundary.knots = c(0,time.max-14/365)))
  val<-c%*%coef
  return(val[1])
}


funct<-function(x,knots,coef,time.max=9.8){
  ## This function computes the value of log incidence rate ratio of hospital admission after vs. 2yrs before MeV
  ## x is the time value
  ## knots is the knot that used for the spines function in the Poisson model
  ## coef is the coeficent of the time after measles
  ## last.knots is the right boundary knot 
  
  c<-c(1,splines::ns(pmax(x-14/365,0), knots=knots-14/365, Boundary.knots = c(0,last.knots-14/365)))
  val<-c%*%coef
  return(val[1])
}
Time_estimate<-function(fit,interval,knots,last.knots){
  ## This function find the zeros of the equation "log incidence rate ratio of hospital admission after vs. 2yrs before MeV =0 "
  ## fit is the fitted object from the Poisson mixed effect model
  ## interval is the anticipated range of the interval that the zeros contained
  ## knots is the knot that used for the spines function in the Poisson model
  ## last.knots is the right boundary knot
  extend<-ifelse(interval[1]>0.5,"upX","downX")
  e <- try( d <-uniroot(funct, c(interval[1], interval[2]), tol = 1e-10,extendInt=extend, knots=knots,coef=unname(coef(fit)[-c(1:10)]),last.knots=last.knots))
  if (class(e) == "try-error") {
    return(Inf)
  } else {
    return(e$root)
  }
}



























mySummary.allvar <- function(blvars,group,pooledGroup=F,contSummary="med.IQR",caption=NULL,filename=NA,pval.comparison=F,cont=NA){
  # contSummary can be median (90% range) "med.90" or median (IQR) "med.IQR" or median (range) "med.range" or "MIC"
  require(xtable); require(Hmisc); library(gdata)
  if (is.factor(group)) {
    levelOrder <- levels(group)
    was.factor <- T
  } else was.factor <- F
  group <- as.character(group)
  if (was.factor) group <- factor(group,levels=levelOrder[!is.na(match(levelOrder,group))])
  else group <- factor(group)
  gr.lev <- levels(group)
  if (pooledGroup){ # Add pooled summaries for all patients
    mylabels <- unlist(lapply(blvars,Hmisc::label)) # save them as rbind destroys some of them
    blvars <- rbind(blvars,blvars)
    for (i in 1:ncol(blvars)) label(blvars[,i]) <- mylabels[i] # add labels again
    group <- c(as.character(group),rep("All patients",length(group)))
    group <- factor(group,levels=c("All patients",gr.lev))
    gr.lev <- levels(group)
  }
  header1 <- c("",c(rbind(rep("",length(gr.lev)),paste(gr.lev," (N=",table(group),")",sep=""))))
  header2 <- c("Characteristic",rep(c("n","Summary statistic"),length(gr.lev)))
  result <-  rbind(header1,header2)
  if (pval.comparison) result <- cbind(result,c("Comparison","(p-value)"))
  if (is.na(cont)[1]) cont <- rep(cont,ncol(blvars))
  for (i in 1:ncol(blvars)){
    result.i <- mySummary.onevar(varname=ifelse(Hmisc::label(blvars[,i])!="",label(blvars[,i]),names(blvars)[i]),
                                 blvars[,i],group,contSummary=contSummary,pval.comparison=pval.comparison,cont=cont[i])
    result <- rbind(result,result.i)
  }
  rownames(result) <- rep("",nrow(result))
  if (!is.na(filename)){  # generate html table
    x <- print(xtable(result,caption=caption),type="html",file=filename,
               caption.placement="top",include.rownames=F,include.colnames=F,
               html.table.attributes="border=1")
    cont.summary.options <- c("med.IQR","median.90","median.range","MIC")      
    cont.summary.text <- c("median (IQR)","median (90% range)","median (range)","median, 90% quantile and range")      
    x <- paste(x,"Summary statistic is absolute count (%) for categorical variables and ",
               cont.summary.text[match(contSummary,cont.summary.options)],"for continuous data.",
               ifelse(pval.comparison,"Comparisons based on Fisher's exact test for categorical data and Wilcoxon test for continuous data.",""))
    write(x,file=filename)
  }
  result
}

mySummary.onevar <- function(varname,variable,group,cont=NA,contSummary="med.90",pval.comparison=F){
  require(gdata)
  if (is.na(cont)) cont <- ifelse(is.factor(variable)|length(unique(na.omit(variable)))<=5,F,T)
  ngroup <- length(levels(group))
  mycont.summary <- function(variable,group) {
    if (is.na(match(contSummary,c("med.90","med.IQR","med.range","MIC"))))
      stop("contSummary=",contSummary," not yet implemented")
    n <- c(by(variable,group,function(x) length(na.omit(x))))
    if (contSummary=="med.90") summarystat <- by(variable,group,function(x) quantile(x,c(0.05,0.5,0.95),na.rm=T))
    if (contSummary=="med.IQR") summarystat <- by(variable,group,function(x) quantile(x,c(0.25,0.5,0.75),na.rm=T))
    if (contSummary=="med.range") summarystat <- by(variable,group,function(x) quantile(x,c(0,0.5,1),na.rm=T))
    if (!is.na(match(contSummary,c("med.90","med.IQR","med.range")))){
      summarystat.nice <- lapply(summarystat,function(x){ x <- formatC(round(x,2),2,format="f");
      paste(x[2],"(",x[1],",",x[3],")",sep="")})
      result <- matrix("",ncol=ngroup*2+1,nrow=1)
      result[1,seq(2,ncol(result),by=2)] <- n
      result[1,seq(3,ncol(result),by=2)] <- unlist(summarystat.nice)
    }
    if (contSummary=="MIC") {
      result <- matrix("",ncol=ngroup*2+1,nrow=4)
      result[,1] <- c("","- MIC 50","- MIC 90","- range")
      result[1,seq(2,ncol(result),by=2)] <- n
      medians <- by(variable,group,function(x) quantile(x,c(0.5),na.rm=T))
      q.90s <- by(variable,group,function(x) quantile(x,c(0.9),na.rm=T))
      ranges <- by(variable,group,function(x) quantile(x,c(0,1),na.rm=T))
      result[2,seq(3,ncol(result),by=2)] <- unlist(lapply(medians,function(x){x <- formatC(round(x,2),2,format="f")}))      
      result[3,seq(3,ncol(result),by=2)] <- unlist(lapply(q.90s,function(x){x <-   formatC(round(x,2),2,format="f")}))
      result[4,seq(3,ncol(result),by=2)] <- unlist(lapply(ranges,function(x) {x <- formatC(round(x,2),2,format="f"); 
      paste("(",x[1],",",x[2],")",sep="")}))
      # # Replace "257" by ">256" (and "33" by ">32" for cotrimoxazole) 
      # result <- gsub("257",">256",result)
      # if (varname=="cotrimoxazole") result <- gsub("33",">32",result)
    }     
    if (pval.comparison) {
      pval <- kruskal.test(variable[group!="All patients"]~group[group!="All patients"])$p.value
      if (pval<0.001) pval <- "<0.001" else pval <- formatC(round(pval,3))
      result <- cbind(result,"")
      result[1,ncol(result)] <- pval
    }
    result
  }
  mycat.summary <- function(variable,group) {
    ta <- table(group,variable)
    ta.n <- apply(ta,1,sum)
    ta.prop <- ta/ta.n
    ta.nice <- matrix(paste(ta,"/",ta.n," (",round(100*ta.prop),"%",")",sep=""),nrow=nrow(ta),ncol=ncol(ta))
    result <- matrix("",ncol=ngroup*2+1,nrow=ncol(ta)+1)
    result[2:nrow(result),1] <- paste("- ",colnames(ta))
    result[2:nrow(result),seq(3,ncol(result),by=2)] <- t(ta.nice)
    result[1,seq(2,ncol(result),by=2)] <- apply(ta,1,sum) # n's
    if (pval.comparison) {
      if (ncol(table(group[group!="All patients"],variable[group!="All patients"]))==1) pval <- "NA"
      else {
        # Use simulated p-value if normal fisher-test doesn't work
        options(show.error.messages=F) 
        ft <- try(fisher.test(table(group[group!="All patients"],variable[group!="All patients"])))
        options(show.error.messages=T) 
        if (class(ft)!="try-error"){
          pval <- ft$p.value
        } else {
          warning("Simulated p-values for Fisher test used with B=10000")
          pval <- fisher.test(table(group[group!="All patients"],variable[group!="All patients"]),
                              simulate.p.value=T,B=10000)$p.value
        }  
        if (pval<0.001) pval <- "<0.001" else pval <- formatC(round(pval,3))
      }
      result <- cbind(result,c(pval,rep("",nrow(result)-1)))
    }
    result
  }
  if (cont) r <- mycont.summary(variable,group)
  else r <- mycat.summary(variable,group)
  r[1,1] <- varname
  r
}
geeglm.summary<- function(obj){
  # Simple summary for a geeglm object
  s<-length(obj$coefficients)
  dat<-data.frame()
  if(s>2){
    for(i in 1:s){
      position<-c(rep(0,(i-1)),1,rep(0,(s-i)))
      tmp<-esticon(obj,position)
      dat.tmp<-data.frame(Estimate=tmp$estimate,Robust.se=tmp$std.error,exp.Estimate=exp(tmp$estimate),Lower.CI=exp(tmp$lwr),Upper.CI=exp(tmp$upr))
      dat<-rbind(dat,dat.tmp)
    }
  }else{
    position<-c(1,0)
    tmp<-esticon(obj,position)
    dat.tmp<-data.frame(Estimate=tmp$estimate,Robust.se=tmp$std.error,exp.Estimate=exp(tmp[,1]),Lower.CI=exp(tmp$lwr),Upper.CI=exp(tmp$upr))
    dat<-rbind(dat,dat.tmp)
    position<-c(0,1)
    tmp<-esticon(obj,position)
    dat.tmp<-data.frame(Estimate=tmp$estimate,Robust.se=tmp$std.error,exp.Estimate=exp(tmp$estimate),Lower.CI=exp(tmp$lwr),Upper.CI=exp(tmp$upr))
    dat<-rbind(dat,dat.tmp)
  }
  dat$covariate<-names(obj$coefficients)
  dat<-dat[,c(6,1:5)]
  dat$p.value<-format.pval(2*(1-pnorm(abs(dat$Estimate/dat$Robust.se),0,1)))
  dat<-dat[-1,]
  rownames(dat)<-NULL
  return(dat)
}

my.ci.lin.mi<-function (obj, ctr.mat = NULL, subset = NULL, subint = NULL,
                        diffs = FALSE, fnam = !diffs, vcov = FALSE, alpha = 0.05,
                        df = Inf, Exp = FALSE, sample = FALSE)
{
  if (sample) require(MASS)
  cf <- coef(obj)
  vcv <- obj$variance
  if (any(is.na(cf))) {
    if (inherits(obj, c("coxph"))) {
      wh <- !is.na(cf)
      cf <- cf[wh]
      vcv <- vcv[wh, wh]
    }
    else if (inherits(obj, c("clogistic"))) {
      cf[is.na(cf)] <- 0
    }
    else {
      vM <- matrix(0, length(cf), length(cf))
      dimnames(vM) <- list(names(cf), names(cf))
      vM[!is.na(cf), !is.na(cf)] <- vcv
      cf[is.na(cf)] <- 0
      vcv <- vM
    }
  }
  all.dif <- function(cf, pre = FALSE) {
    nn <- length(cf)
    nr <- nn * (nn - 1)/2
    nam <- names(cf)
    xx <- numeric(0)
    for (j in 2:nn) xx <- c(xx, j:nn)
    ctr <- cbind(rep(1:(nn - 1), (nn - 1):1), xx)
    i <- 1
    while (all(substr(nam, 1, i) == substr(nam[1], 1, i))) i <- i +
      1
    if (is.character(pre)) {
      prefix <- pre
      pre <- TRUE
    }
    else {
      prefix <- substr(nam[1], 1, i - 1)
    }
    rn <- paste(if (pre)
      prefix
      else "", substring(nam[ctr[, 1]], i), "vs.", substring(nam[ctr[,
                                                                     2]], i))
    cm <- matrix(0, nr, nn)
    cm[cbind(1:nr, ctr[, 1])] <- 1
    cm[cbind(1:nr, ctr[, 2])] <- -1
    rownames(cm) <- rn
    cm
  }
  if (diffs) {
    if (is.character(subset)) {
      if (inherits(obj, "lm") & length(grep(subset, names(obj$xlevels))) >
          0) {
        wf <- grep(subset, af <- names(obj$xlevels))
        fn <- obj$xlevels[[af[wf]]]
        pnam <- paste(af[wf], fn, sep = "")
        wh <- match(pnam, names(coef(obj)))
        cf <- coef(obj)[wh]
        cf[is.na(cf)] <- 0
        vcv <- vcov(obj)[wh, wh]
        vcv[is.na(vcv)] <- 0
        names(cf) <- rownames(vcv) <- colnames(vcv) <- paste(subset,
                                                             ": ", fn, sep = "")
      }
      else {
        subset <- grep(subset, names(cf))
        cf <- cf[subset]
        vcv <- vcv[subset, subset]
      }
    }
    else {
      cf <- cf[subset]
      vcv <- vcv[subset, subset]
    }
    ctr.mat <- all.dif(cf, pre = fnam)
  }
  if (!diffs) {
    if (is.character(subset)) {
      sb <- numeric(0)
      for (i in 1:length(subset)) sb <- c(sb, grep(subset[i],
                                                   names(cf)))
      subset <- sb
    }
    if (is.character(subint)) {
      sb <- 1:length(cf)
      for (i in 1:length(subint)) sb <- intersect(sb, grep(subint[i],
                                                           names(cf)))
      subset <- sb
    }
    if (is.null(subset) & is.null(subint))
      subset <- 1:length(cf)
    cf <- cf[subset]
    vcv <- vcv[subset, subset]
    if (is.null(ctr.mat)) {
      ctr.mat <- diag(length(cf))
      rownames(ctr.mat) <- names(cf)
    }
    if (dim(ctr.mat)[2] != length(cf))
      stop(paste("\n Dimension of ", deparse(substitute(ctr.mat)),
                 ": ", paste(dim(ctr.mat), collapse = "x"), ", not compatible with no of parameters in ",
                 deparse(substitute(obj)), ": ", length(cf), sep = ""))
  }
  ct <- ctr.mat %*% cf
  vc <- ctr.mat %*% vcv %*% t(ctr.mat)
  if (sample)
    res <- t(mvrnorm(sample, ct, vc))
  else {
    se <- sqrt(diag(vc))
    ci <- cbind(ct, se) %*% ci.mat(alpha = alpha, df = df)
    t0 <- cbind(se, ct/se, 2 * (1 - pnorm(abs(ct/se))))
    colnames(t0) <- c("StdErr", "z", "P")
    res <- cbind(ci, t0)[, c(1, 4:6, 2:3), drop = FALSE]
    if (Exp) {
      res <- cbind(res[, 1:4, drop = FALSE], exp(res[,
                                                     c(1, 5, 6), drop = FALSE]))
      colnames(res)[5] <- "exp(Est.)"
    }
  }
  if (sample)
    invisible(res)
  else if (vcov)
    invisible(list(est = ct, vcov = vc))
  else res
}

my.ci.lin<-function (obj, ctr.mat = NULL, subset = NULL, subint = NULL,
          diffs = FALSE, fnam = !diffs, vcov = FALSE, alpha = 0.05,
          df = Inf, Exp = FALSE, sample = FALSE)
{
  if (sample) require(MASS)
  cf <- coef(obj)
  vcv <- obj$geese$vbeta
  if (any(is.na(cf))) {
    if (inherits(obj, c("coxph"))) {
      wh <- !is.na(cf)
      cf <- cf[wh]
      vcv <- vcv[wh, wh]
    }
    else if (inherits(obj, c("clogistic"))) {
      cf[is.na(cf)] <- 0
    }
    else {
      vM <- matrix(0, length(cf), length(cf))
      dimnames(vM) <- list(names(cf), names(cf))
      vM[!is.na(cf), !is.na(cf)] <- vcv
      cf[is.na(cf)] <- 0
      vcv <- vM
    }
  }
  all.dif <- function(cf, pre = FALSE) {
    nn <- length(cf)
    nr <- nn * (nn - 1)/2
    nam <- names(cf)
    xx <- numeric(0)
    for (j in 2:nn) xx <- c(xx, j:nn)
    ctr <- cbind(rep(1:(nn - 1), (nn - 1):1), xx)
    i <- 1
    while (all(substr(nam, 1, i) == substr(nam[1], 1, i))) i <- i +
      1
    if (is.character(pre)) {
      prefix <- pre
      pre <- TRUE
    }
    else {
      prefix <- substr(nam[1], 1, i - 1)
    }
    rn <- paste(if (pre)
      prefix
      else "", substring(nam[ctr[, 1]], i), "vs.", substring(nam[ctr[,
                                                                     2]], i))
    cm <- matrix(0, nr, nn)
    cm[cbind(1:nr, ctr[, 1])] <- 1
    cm[cbind(1:nr, ctr[, 2])] <- -1
    rownames(cm) <- rn
    cm
  }
  if (diffs) {
    if (is.character(subset)) {
      if (inherits(obj, "lm") & length(grep(subset, names(obj$xlevels))) >
          0) {
        wf <- grep(subset, af <- names(obj$xlevels))
        fn <- obj$xlevels[[af[wf]]]
        pnam <- paste(af[wf], fn, sep = "")
        wh <- match(pnam, names(coef(obj)))
        cf <- coef(obj)[wh]
        cf[is.na(cf)] <- 0
        vcv <- vcov(obj)[wh, wh]
        vcv[is.na(vcv)] <- 0
        names(cf) <- rownames(vcv) <- colnames(vcv) <- paste(subset,
                                                             ": ", fn, sep = "")
      }
      else {
        subset <- grep(subset, names(cf))
        cf <- cf[subset]
        vcv <- vcv[subset, subset]
      }
    }
    else {
      cf <- cf[subset]
      vcv <- vcv[subset, subset]
    }
    ctr.mat <- all.dif(cf, pre = fnam)
  }
  if (!diffs) {
    if (is.character(subset)) {
      sb <- numeric(0)
      for (i in 1:length(subset)) sb <- c(sb, grep(subset[i],
                                                   names(cf)))
      subset <- sb
    }
    if (is.character(subint)) {
      sb <- 1:length(cf)
      for (i in 1:length(subint)) sb <- intersect(sb, grep(subint[i],
                                                           names(cf)))
      subset <- sb
    }
    if (is.null(subset) & is.null(subint))
      subset <- 1:length(cf)
    cf <- cf[subset]
    vcv <- vcv[subset, subset]
    if (is.null(ctr.mat)) {
      ctr.mat <- diag(length(cf))
      rownames(ctr.mat) <- names(cf)
    }
    if (dim(ctr.mat)[2] != length(cf))
      stop(paste("\n Dimension of ", deparse(substitute(ctr.mat)),
                 ": ", paste(dim(ctr.mat), collapse = "x"), ", not compatible with no of parameters in ",
                 deparse(substitute(obj)), ": ", length(cf), sep = ""))
  }
  ct <- ctr.mat %*% cf
  vc <- ctr.mat %*% vcv %*% t(ctr.mat)
  if (sample)
    res <- t(mvrnorm(sample, ct, vc))
  else {
    se <- sqrt(diag(vc))
    ci <- cbind(ct, se) %*% ci.mat(alpha = alpha, df = df)
    t0 <- cbind(se, ct/se, 2 * (1 - pnorm(abs(ct/se))))
    colnames(t0) <- c("StdErr", "z", "P")
    res <- cbind(ci, t0)[, c(1, 4:6, 2:3), drop = FALSE]
    if (Exp) {
      res <- cbind(res[, 1:4, drop = FALSE], exp(res[,
                                                     c(1, 5, 6), drop = FALSE]))
      colnames(res)[5] <- "exp(Est.)"
    }
  }
  if (sample)
    invisible(res)
  else if (vcov)
    invisible(list(est = ct, vcov = vc))
  else res
}

pred.gee<-function(obj,newdata){
  dat<-newdata
  formula<-paste("~",as.character(formula(obj)[3],sep=""))
  mm<-model.matrix(as.formula(formula),dat)
  pvar1 <- diag(mm %*% tcrossprod(obj$geese$vbeta,mm))
  dat$pred<-mm%*%obj$geese$beta #predict(m,newdat,re.form=NA) would give the same results
  dat$plo<-dat$pred-qnorm(0.975,0,1)*sqrt(pvar1)
  dat$phi <- dat$pred+qnorm(0.975,0,1)*sqrt(pvar1)
  dat$sd<-sqrt(pvar1)
  return(dat)
}


pred.gee.diff<-function(obj,newdata){
  dat<-newdata
  formula<-paste("~",as.character(formula(obj)[3],sep=""))
  mm<-model.matrix(as.formula(formula),dat)
  
  dat$absolute_diff<-exp(mm%*%obj$geese$beta)- exp(mm[,1:10]%*%obj$geese$beta[1:10])
  n<-length(obj$coefficients)
  grad<-matrix(ncol=n,nrow=nrow(dat))
  for (j in 1:n){
    if(j<11){
      grad[,j]<-mm[,j]*(exp(mm%*%obj$geese$beta)-exp(mm[,1:10]%*%obj$geese$beta[1:10]))
    }else{
      grad[,j]<-mm[,j]*(exp(mm%*%obj$geese$beta)) 
    }
  }
  
  vb <- obj$geese$vbeta
  vG <- (grad) %*% vb %*% t(grad)
  pvar1 <- diag(vG)
  dat$absolute_diff_lo<-dat$absolute_diff-qnorm(0.975,0,1)*sqrt(pvar1)
  dat$absolute_diff_high <- dat$absolute_diff+qnorm(0.975,0,1)*sqrt(pvar1)
  return(dat)
}

