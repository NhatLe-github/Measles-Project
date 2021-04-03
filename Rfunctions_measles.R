library(tidyverse)
library(Hmisc)
library(magrittr)
require(geepack)
require(xtable)
require(doBy)
require(readr)
library(lubridate)
library(matrixStats)
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
  # # option = 0 if  the discharged with presumably worse condition alive
  # # option = 1 if  the discharged with presumably worse condition dead
  # # option = 2 if  the emigration patient out of the HCMC city
  
  #### ------------------------------ Rule to derive dataset ------------------------------ ####
  # 1. Only the admissions due to infection will be included in the analysis
  # 2. Impute the birth day for who doesn't have birthday record but have the birth years. 
  #     We immputed it as 30June of the birth years. Who doesn't have birth year 
  #     will be removed from the analysis.
  # 3. We only considered a new admission due to a new episode of infection if the 
  #     new admission date after the last discharge at least 2 days, otherwise we
  #     merge them into one espisode of new infections
  # 4. We assumed that for each episode of admission was due to one main infection type
  # 5. All measles second infections after the first admission
  #     were considered as the measles complicated
  
  
  ## Classify ICD.10 code into common type of infections
  icd.infection.diarrhea    <- paste(rep("A",10),c("00","01","02","03","04","05","06","07","08","09"),sep = "")
  icd.infection.respiratory <- paste(rep("J",23),c("00","01","02","03","04","05","06","07","08","09",10:22),sep = "")
  icd.infection.CNS         <- c(paste(rep("A",10),80:89,sep =""),paste(rep("G",9),c("00","01","02","03","04","05","06","07","09"),sep = ""))
  icd.infection             <- c(c(paste(rep("A",100),c("00","01","02","03","04","05","06","07","08","09",10:99),sep = ""),
                               paste(rep("B",100),c("00","01","02","03","04","05","06","07","08","09",10:99),sep = "")),
                               paste(rep("J",23),c("00","01","02","03","04","05","06","07","08","09",10:22),sep = ""),
                               paste(rep("G",10),c("00","01","02","03","04","05","06","07","08","09"),sep = ""))
          
  ## Rename the id to patid and select the impotrtant variables from the dataset
  ## Re-order the dataset and create viariable visit and visit last
  data %<>% rename(patid=iduse) 
  data %<>% dplyr:: select(patid,sex,dob,birthyear,date.admit,date.dis,diagnosis.m3,diagnosis.en,icd.main,icd.m3,
                           outcome,ward,district,province,year)
  data %<>% arrange(patid,date.admit,date.dis)
  data  <-  visit_provide(data)
  
  ## Remove all duplicated cases with the same admission and discharge. 
  ## We removed the first records: Two patients were excluded 174644 210565 
  
  data %<>% group_by(patid,date.admit,date.dis) %>% dplyr::mutate(visit.max=max(visit)) %>% ungroup()  
  data %<>% mutate(select_id=ifelse(visit==visit.max,1,0))  
  data %<>% filter(select_id==1)  
  
  ## Impute all the missing of reasons of discharge --> "unknown".
  data %<>% mutate(outcome = as.character(outcome)) %>% 
            mutate(outcome = as.factor(ifelse(!is.na(outcome),outcome,"unknown")))
  

  ## Merge all the cases have the same admission but different discharge 
  ## Select the record with the later hospital discharge (--> re-order the data by date.admit and  date.dis )
  ## Make sure that for all the patients when we merged two identical admissions, we didn't miss the measles cases i.e. icd.10 code ="B05"
  ## 6 patients 5800   8330  65942 284552 410822 582807 had two hopsital admissions record with the same admission
  
  data %<>% arrange(patid,date.admit,date.dis) 
  data %<>% group_by(patid,date.admit) %>%
            dplyr::mutate(date.dis.max=max(date.dis)) %>%
            dplyr::mutate(select.id=ifelse(date.dis==date.dis.max,1,0)) %>%
            ungroup() 
  
  data %<>% mutate(icd.m3=ifelse((lag(select.id,1) == 0  &
                                         lag(patid) == patid) &
                                        (icd.m3 == "B05" | lag(icd.m3,1) == "B05"),"B05",icd.m3))
  data$icd.m3[1] <- "J02"

  data %<>%  filter(select.id == 1)
  data[,c("id.date","select.id")] <- NULL
  
  ## Re-defined the outcome presumably.worse. If the patient re-admitted to the hospital after 
  ## being discharged with presumably.worse, the outcome presumably.worse were considered as discharged.
  ## Take a subset of patients with outcome presumably.worse (dat),
  ## and then determine the last hospital admission of these patients
  
  dat <- data %>% filter(outcome == "presumably.worse")
  id_worse <- unique(dat$patid)
  data %<>% arrange(patid,date.admit,date.dis) 
  data  <-  visit_provide(data)
  dat <- data %>% filter(patid %in% id_worse) %>% arrange(patid, visit)
  
  ### If patient re-admited to the hospial after the visit with outcome presumably.worse (i.e. not the last visit) then 
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
  
  ## Impute the missing birthday information
  ## and exclude patients without the birth year
  data.imputed.birth <- data %>% filter((is.na(dob)) & (!is.na(birthyear))) %>% 
                                mutate(dob = as.Date(paste("30June",birthyear,sep = ""),"%d%b%Y",origin = "1970-01-01"))
  data <- rbind(subset(data,!is.na(data$dob)),data.imputed.birth)
  
  
  ### Create the date of death: date_death
  data %<>% mutate(date_death=as.Date(ifelse(outcome == "dead",date.dis,NA),origin = "1970-01-01")) 
  
  ## In this analysis, we set the time to follow-up is the earliest of 31 December 2015 ( i.e., the cut-off date of the dataset)
  ## and the 15th birthday of the child 
  data %<>% mutate(time.at.15birthday = as.Date(paste(birthyear + 15,month(dob),day(dob),sep = "-"),origin = "1970-01-01"),
                   time.at.lastfup = pmin(time.at.15birthday,as.Date("2015-12-31",origin = "1970-01-01")),
                   time.at.lastfup = as.Date(ifelse(outcome == "dead",data$date_death,time.at.lastfup),origin = "1970-01-01")) 
  
  ## Re-assess the time at last follow-up for all individual. Because some one died then the time at last follow-up was the date of death
  data %<>% group_by(patid) %>% dplyr::mutate(time.at.lastfup=min(time.at.lastfup)) %>% ungroup()
  data %<>% mutate(time.dis.to.lastfup = as.numeric(difftime(time.at.lastfup,date.dis,units = "days"))) 
  
  ### Remove all patients who admitted after the time follow-up 
  ### Patient 422164 were older than 15 years old so he was excluded from the analysis
  ### For records with time.dis.to.lastfup= 0 i.e. they died at that admission
  ### create a variable measles, and infection
  
  data  %<>% filter(time.dis.to.lastfup >= 0) %>% 
             arrange(patid,date.admit,date.dis) %>% 
             select("patid","sex","district","dob","birthyear","date.admit","date.dis","icd.m3","outcome","year","diagnosis.en","time.at.lastfup","time.dis.to.lastfup") %>% 
             mutate(measles = ifelse(icd.m3 == "B05","yes","no"),
             infection = ifelse(icd.m3 %in% icd.infection,"yes","no"))
  
  ## Select only records with infection and patient who admited due to MeV infection
  data  %<>%  filter(infection == "yes") %>% 
              mutate(age.at.admit.d = as.numeric(difftime(date.admit,dob,units = "days"))) %>% 
              filter(patid %in% unique(filter(.,measles == "yes")$patid))

  ## Select only records from MeV infection and then determine which were the first measles and second measles infection.
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
  ##  was measles then the icd.m3 of merged case would be measles otherwise it would
  ##  be the icd of the latter case (which was not important to this analysis).
  ##  Compute the time difference between two consecutive admission and
  ##  Set all the value for the first tim_diff to 10000 (a very large value) for each individual.
  
  data %<>% group_by(patid) %>% dplyr::mutate(tim_diff=as.numeric(difftime(date.admit,lag(date.dis,1),unit = "day"))) %>% ungroup() 
  data %<>% mutate(tim_diff=ifelse(is.na(tim_diff),10000,tim_diff))
  
  ## Define a type of pair consecutive admission within 2 days
  ## For each pair, create variable index_pair to indicate: 
  ## index_pair = MeV_bf if MeV event occured before the other event
  ## index_pair = MeV_af if MeV event occured after the other event
  ## index_pair = MeV_both if both admissions were due to MeV
  ## index_pair = non_MeV if both admission were not due to MeV
  ## index_pair = non pairs if not consecutive admission within 2 days
  
  data %<>% group_by(patid) %>% dplyr::mutate(index_pair = case_when(
    tim_diff<2 & (icd.m3 != "B05" & lag(icd.m3,1) == "B05") ~ 'MeV_bf',
    tim_diff<2 & (icd.m3 == "B05" & lag(icd.m3,1) != "B05") ~ 'MeV_af',
    tim_diff<2 & (icd.m3 == "B05" & lag(icd.m3,1) == "B05") ~ 'MeV_both',
    tim_diff<2 & (icd.m3 != "B05" & lag(icd.m3,1) != "B05") ~ "non_MeV",
                                                TRUE ~ 'non pairs'
                                              )) %>% ungroup()

  data %<>% arrange(patid,date.admit,date.dis) 
  data  <-  visit_provide(data)

  
  ## For all records with the same admissions, 
  ## index_select select the second record of the pair of the same admission
  ## visit_select make the pair with the previous visit which has the same admission by set it to the same visit.
  ## index_select_pair indicate the pair of the same admission
  
  data %<>% group_by(patid) %>% 
    dplyr::mutate(index_select=ifelse(index_pair%in%c("MeV_bf","MeV_af","MeV_both","non_MeV"),1,0),
                  visit_select=ifelse(index_select==1,visit-1,visit)) %>% ungroup()
  
  ## Select the pair which links to the second case
  data %<>% group_by(patid,visit_select) %>% dplyr::mutate(index_select_pair=max(index_select)) %>% ungroup()
  
  ## Re-assess the icd.m3 of the second measles  
  data %<>% mutate(icd.m3 =ifelse(measles.order=="second","B05*",icd.m3)) 
  
  ## Merge information of two admission into one by date.admit = minimum of admision date of two admissions and 
  ## date.dis = maximum of discharged date of two admissions 
  ## If one of two admissions was measles then icd.m3 of both admissions will be merged to MeV otherwise it didn't matter
  ## Since we would chose the second record (index_select==1) of each pair so for the case of MeV_bf and MeV_both the icd was not B05 (i.e. measles) so we set it as MeV

  
  dat.tmp.pair<-data %>% filter(index_select_pair==1)
  dat.tmp.pair %<>% group_by(patid) %>% dplyr::mutate(date.admit=min(date.admit),
                                                 date.dis=max(date.dis),
                                                 icd.m3= ifelse(index_pair%in%c("MeV_bf","MeV_both"),"B05",icd.m3)
                                                 ) %>% ungroup()
  ## Re-assess variables measles, and measles order.
  
  dat.tmp.pair %<>% filter(index_select==1) %>% mutate(
                                                  measles=ifelse(icd.m3=="B05","yes","no"),
                                                  measles.order = ifelse(icd.m3=="B05","first",measles.order))
  
  dat.tmp.non.pair<-data%>% filter(index_select_pair==0)
  dat.tmp.non.pair %<>% mutate(
                              measles=ifelse(icd.m3=="B05","yes","no"),
                              measles.order = ifelse(icd.m3=="B05","first",measles.order))
  
  data <-rbind(dat.tmp.non.pair,dat.tmp.pair)
  data %<>% arrange(patid,date.admit,date.dis) 
  data  <-  visit_provide(data)
  
  ## Classify the type of infections. 
  ## For the second measles admission will be considered as measles complicated.
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
                   outcome = as.factor(outcome),
                   infection.type = as.factor(infection.type),
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
    dat <- data.after %>%  arrange(patid,date.admit)
    dat  <-  visit_provide(dat)
    # dat.sum <- dat %>% group_by(patid) %>% dplyr::summarise(visit.last = max(visit),last.year.ev = max(year),date.dis.last=as.Date(max(date.dis),origin = "1970-01-01"),.groups="keep") %>% ungroup() 
    # dat <- merge(dat,dat.sum,by = "patid",all.x=T)
    
    dat <- dat %>% group_by(patid) %>% dplyr::mutate(last.year.ev = max(year),date.dis.last=as.Date(max(date.dis),origin = "1970-01-01")) %>% ungroup() 

    
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
  
  ## test if time to measles have na value
  ## sum(is.na(dat.tmp.before$time.dis.to.ad.measles))

  ## Create dataset of 2 year before MeV event
  ## Impute the follow-up time between the last event to MeV event
  ## flag.index.impute = 1 i.e. MeV event with a previous event
  ## flag.index.impute = 2 i.e. MeV event without a previous event
  ## flag.index.impute = 3 i.e. other type of event 
  
  dat.tmp.before.years.before <- dat.tmp.before %>% filter(time.ad.to.ad.measles<= time.period) %>% 
    arrange(patid,date.admit,date.dis) %>% group_by(patid) %>% 
    dplyr::mutate(status = ifelse(period == "measles",0,1),
                  flag.index.impute = ifelse(period=="measles",ifelse((!is.na(lag(patid,1))),1,2),3))  %>% ungroup()
  
  ## For flag.index.impute = 2 and = 3, we define the obs.time = the period since 2 year before MeV to the current event
  dat.tmp.before.years.before_2_3 <- dat.tmp.before.years.before %>% 
    filter(flag.index.impute %in%c(2,3)) %>% 
    mutate(obs.time = ifelse(flag.index.impute==2,ifelse(age.at.measles.d >= time.period,time.period,age.at.measles.d),
                             ifelse(age.at.measles.d >= time.period,pmax(time.period - abs(time.ad.measles.to.ad),0),
                                    age.at.admit.d)))
  
  ## during the time.period = 2 years
  ## For the first event of flag.index.impute = 3 or flag.index.impute = 2, the time start is 0.
  ## For subsequent event with flag.index.impute = 3, the time start is the discharge time of the previous event, i.e lag(obs.time,1)+lag(hospialization,1)
  ## And the time stop of the current event is the obs.time
  
  dat.tmp.before.years.before_2_3 %<>% 
    arrange(patid,date.admit,date.dis) %>% group_by(patid) %>% 
    dplyr::mutate(obs.time = ifelse(date.admit == dob & age.at.measles.d < time.period,1,obs.time),# for special cases the child were born and admitted to hospital
                  time.start = ifelse(!(is.na(lag(patid,1))),
                                      lag(obs.time,1)+lag(hospitalization,1),0)
                  ) %>% ungroup()
  

  dat.tmp.before.years.before_2_3 %<>% mutate(time.stop  = obs.time,
                                              time.stop  = ifelse(time.stop == 0,time.start+1,time.stop)) %>% 
                                       filter(time.stop > 0)
  
  ## For flag.index.impute ==1, the MeV event occured after a previous event with flag.index.impute ==3,
  ## That was also the last event non MeV occured  befored MeV during the 2 years period.
  ## The time start was the time stop of the previous event + hopitalisation,
  ## The time stop was the time period for children above 2 years ole and = age.at measles for children under 2 years old.
  
  dat.tmp.before.years.before_2_3<-visit_provide(dat.tmp.before.years.before_2_3)
  dat.tmp.before.years.before_1 <- dat.tmp.before.years.before_2_3 %>% filter(visit==visit.last & flag.index.impute ==3) %>% 
                                      mutate(time.start = time.stop + hospitalization,
                                             time.stop = ifelse(age.at.measles.d >= time.period,time.period,age.at.measles.d),
                                             status = 0) %>% 
                                      arrange(patid,date.admit,date.dis)

  ### Merge two dataset dat.tmp.before.years.before_1 and dat.tmp.before.years.before_2_3
  dat.tmp.before.years.before <- rbind(dat.tmp.before.years.before_2_3,dat.tmp.before.years.before_1)
  dat.tmp.before.years.before[,c("measles","period","icd.m3")]<-NULL
  dat.tmp.before.years.before %<>% filter(time.stop > time.start) %>% arrange(patid,time.start)
  
  (n10 <- length(unique(dat.tmp.before.years.before$patid)))
  
  ## Transforming the dataset into daily interval observations
  dat.tmp.before.split <- survSplit(Surv(time.start,time.stop, status)~.,dat.tmp.before.years.before,cut=c(seq(0,time.period,by=1)),episode ="days")
  
  dat.tmp.before.split$days <- dat.tmp.before.split$days - 1
  dat.tmp.before.split$age.at.infection <- ifelse(dat.tmp.before.split$age.at.measles.d >= time.period,
                                                  dat.tmp.before.split$days/365 + (dat.tmp.before.split$age.at.measles.d - time.period)/365,
                                                  dat.tmp.before.split$days/365)
  
  (n10 <- length(unique(dat.tmp.before.split$patid)))
  
  data.tmp.after <- data %>% filter(period %in% c("after","measles")) %>% arrange(patid,date.admit,date.dis)
  (n10<-length(unique(data.tmp.after$patid)))
  data.tmp.after <- visit_provide(data.tmp.after)

  ## We remove the washout period of 2 weeks=14 days.
  removed <- TRUE
  ## collapse two consecutive episode if it happened within 2 weeks
  if(removed == TRUE) {
    # update the new discharge date of measles
    data.tmp.after %<>% mutate(date.last14.measles = date.admit.measles + 14) %>%
      arrange(patid,date.admit,date.dis) %>% group_by(patid) %>%
      dplyr::mutate(selected.removed = ifelse(!is.na(lag(patid,1)) & (date.admit<= date.last14.measles),1,0),
                    date.dis.measles.updated = pmax(date.last14.measles,date.dis.measles)) %>% ungroup() %>%
      arrange(patid,desc(date.admit))

    ### Since there is no 3 admission within 2 week since measles admission, so we defined the updated date.dis.measles
    ### Remove all cases occured within 2 weeks
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

  ### Define the flag.index.impute variable is the indication of type of event
  ### flag.index.impute =1 if MeV without subsequent other event
  ### flag.index.impute =2 if MeV with subsequent other event
  ### flag.index.impute =3 if  other event before the last non-MeV event
  ### flag.index.impute =4 if  the last non-MeV event
  
  data.tmp.after %<>% arrange(patid,desc(date.admit)) 
  data.tmp.after <- visit_provide(data.tmp.after)
  
  data.tmp.after %<>%   group_by(patid) %>%
                        dplyr::mutate(flag.index.impute = ifelse(period == "after",
                                                                 ifelse(visit==visit.last,4,3),
                                                                 ifelse(is.na(lag(patid,1)),1,2)
                                                                 )
                                      ) %>% ungroup()
  
  data.tmp.after %<>% arrange(patid,date.admit,date.dis)



  ### Compute obs.time.since.MeV and set obs.time.since.MeV==-14 for measles cases
  ### Compute time to follow-up for each patient for the group of non MeV event or MeV event with subsequent other events
  ### For all the case with flag.index.impute = 3 and 4, the time start of the event is the time stop of the previous event + hospitalization.
  ### For all the case with flag.index.impute = 1 and 2, the time start of the event is 0.
  data.tmp.after %<>%
    arrange(patid,date.admit,date.dis) %>%
    mutate(#obs.time.since.MeV = ifelse(period == "after",as.numeric(difftime(date.admit,date.dis.measles,units="days")),0),# the duration from measles infection to the other admission
      obs.time.since.MeV = as.numeric(difftime(date.admit,date.dis.measles,units="days")),     
      time.dis.MeV.to.lastfup = as.numeric(difftime(time.at.lastfup,date.dis.measles,units="days")),
           hospitalization.measles = as.numeric(difftime(date.dis.measles,date.admit.measles,units="days")))
  (n10<-length(unique(data.tmp.after.tmp$patid)))

  data.tmp.after %<>%  arrange(patid,date.admit,date.dis) %>%
    mutate(lag.obs.time.since.MeV = lag(obs.time.since.MeV,1),
           time.start = ifelse(flag.index.impute%in%c(3,4),obs.time.since.MeV+hospitalization,0))

  (n10<-length(unique(data.tmp.after$patid)))

  ## The time.stop of flag.index.impute = 1 is the last follow-up time
  ## The time.stop of flag.index.impute = 2 and 3 is the time of the event of the next event

  data.tmp.after %<>% mutate(time.stop   = ifelse(flag.index.impute %in%c(1,4) ,time.dis.MeV.to.lastfup,lead(obs.time.since.MeV,1)),
                             status   = ifelse(flag.index.impute %in% c(1,4),0,1))
  (n10<-length(unique(data.tmp.after.tmp$patid)))

  
  
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
  ## time.max is the right boundary knot 
  
  c<-c(1,splines::ns(pmax(x-14/365,0), knots=knots-14/365, Boundary.knots = c(0,time.max-14/365)))
  val<-c%*%coef
  return(val[1])
}


funct<-function(x,knots,coef,time.max=9.8){
  ## This function computes the value of log incidence rate ratio of hospital admission after vs. 2yrs before MeV
  ## x is the time value
  ## knots is the knot that used for the spines function in the Poisson model
  ## coef is the coeficent of the time after measles
  ## time.max is the right boundary knot 
  
  c<-c(1,splines::ns(pmax(x-14/365,0), knots=knots-14/365, Boundary.knots = c(0,time.max-14/365)))
  val<-c%*%coef
  return(val[1])
}
Time_estimate<-function(fit,interval,knots,time.max){
  ## This function find the zeros of the equation "log incidence rate ratio of hospital admission after vs. 2yrs before MeV =0 "
  ## fit is the fitted object from the Poisson mixed effect model
  ## interval is the anticipated range of the interval that the zeros contained
  ## knots is the knot that used for the spines function in the Poisson model
  ## time.max is the right boundary knot
  extend<-ifelse(interval[1]>0.5,"upX","downX")
  e <- try( d <-uniroot(funct, c(interval[1], interval[2]), tol = 1e-10,extendInt=extend, knots=knots,coef=unname(coef(fit)[-c(1:10)]),time.max=time.max))
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

visit_provide<-function(data){
  data%<>% mutate(flag=1)
  data %<>% mutate(flag=1) %>% group_by(patid) %>% 
    dplyr::mutate(visit=cumsum(flag),visit.last = max(visit)) %>% ungroup() 
  return(data)
}

