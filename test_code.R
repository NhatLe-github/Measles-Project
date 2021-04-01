# N= choose number of test id
data.raw <- as.data.frame(read_csv("C:/Users/nhatlth/OneDrive/Measles_paper/Final Paper/Project_Measles_18EN/R-codes/Data/measa.csv"))
data.raw %<>% mutate(outcome=ifelse(is.na(outcome),"unknown",outcome)) %<>% rename(patid=iduse) %>% arrange(patid,date.admit,date.dis)

N=10
id_sample<-data.raw$patid[sample(1:nrow(data.raw),N)]

time.period<- 2*365
#We will check the number of follow-up day that these patients were included
dat.test<-subset(data.raw,patid%in%id_sample)

dat.test.tmp<-dat.test %>% mutate(date.admit.measles=as.Date(ifelse(icd.m3=="B05",date.admit,NA),origin = "1970-01-01"),
                                  date.dis.measles=as.Date(ifelse(icd.m3=="B05",date.dis,NA),origin = "1970-01-01"),
                                  age.at.measles.d = as.numeric(difftime(date.admit.measles,dob,units = "days"))
                                  ) %>%
                           filter(!is.na(date.admit.measles)) %>% select(patid,date.admit.measles,date.dis.measles,age.at.measles.d)



dat.test <- left_join(dat.test,dat.test.tmp,by="patid")
dat.test %<>% mutate(time.ad.measles.to.ad=as.numeric(difftime(date.admit,date.admit.measles,units = "days"))) %>% 
              mutate(time.at.15birthday = as.Date(paste(birthyear + 15,month(dob),day(dob),sep = "-"),origin = "1970-01-01"),
                     time.at.lastfup = pmin(time.at.15birthday,as.Date("2015-12-31",origin = "1970-01-01")),
                     time.at.lastfup = as.Date(ifelse(outcome == "dead",data$date_death,time.at.lastfup),origin = "1970-01-01"),
                     time.dis.MeV.to.lastfup = as.numeric(difftime(time.at.lastfup,date.dis.measles,units = "days")),
                     age.at.admit.d = as.numeric(difftime(date.admit,dob,units = "days"))) %>% 
              mutate(hospitalization=as.numeric(difftime(date.admit,date.dis,units = "days")))

data.tmp <-dat.test %>% group_by(patid) %>% summarise(time.at.lastfup=min(time.at.lastfup),.groups="keep") %>% ungroup()
dat.test$time.at.lastfup<-NULL
dat.test %<>%left_join(data.tmp,by="patid") 


dat.test %<>% 
  mutate(period = as.factor(ifelse(time.ad.measles.to.ad > 0,"after",ifelse(time.ad.measles.to.ad < 0,"before","measles")))) %>% 
  filter(time.ad.measles.to.ad>= -2*365) 


dat.test_before<-dat.test%>%filter(period!="after") %>% mutate(obs.time=NA) 

dat.test_before <-dat.test_before %>% mutate(obs.time = ifelse(period == "measles",ifelse(age.at.measles.d >= time.period,time.period,age.at.measles.d),
                           ifelse(age.at.measles.d >= time.period,pmax(time.period - abs(time.ad.measles.to.ad),0),
                                  age.at.admit.d))) 
dat.test_before%<>%  group_by(patid) %>% 
  dplyr::summarise(obs.time = max(obs.time),.groups="keep") %>% ungroup()

dat.test_after<-dat.test%>%filter(period!="before") %>%  group_by(patid) %>% 
  dplyr::summarise(obs.time = max(time.dis.MeV.to.lastfup),.groups="keep") %>% ungroup()

(dat.test.res<-rbind(dat.test_before,dat.test_after) %>% group_by(patid) %>% 
  summarise(obs.time = sum(obs.time),.groups="keep")) %>% ungroup() %>% mutate(patid=as.factor(patid)) %>% as.data.frame()

(dat_check<-data.total.main %>% filter(patid%in%id_sample) %>% 
  mutate(flag=1) %>% group_by(patid) %>% 
    dplyr::summarise(obs.time_produced = sum(flag)*7,.groups="keep")) %>% ungroup() %>% as.data.frame()

(data.check<-cbind(dat_check,dat.test.res))
