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
    dplyr::mutate(selected.removed = ifelse(!is.na(lag(patid,1)) & (date.admit<= date.last14.measles),1,0),
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

### Define the flag.index.impute variable is the indication of type of event 
### flag.index.impute =1 if last event record and =0 otherwise

tmp.1 <-data.tmp.after%>% arrange(patid,date.admit) %>% group_by(patid) %>% summarise(visit=cumsum(flag),.groups="keep") %>% ungroup()
data.tmp.after$visit<-NULL
data.tmp.after <-cbind(data.tmp.after,tmp.1[,"visit"])
tmp.2 <-tmp.1 %>% group_by(patid) %>% summarise(last.visit=max(visit),.groups="keep")  %>% ungroup()
data.tmp.after<-left_join(data.tmp.after,tmp.2,by="patid")
data.tmp.after %<>% mutate(flag.index.impute=ifelse(visit==last.visit,1,0))
### Compute obs.time.since.MeV and set obs.time.since.MeV==-14 for measles cases
### Compute time to follow-up for each patient for the group of non MeV event or MeV event with subsequent other events
data.tmp.after<- data.tmp.after %>% filter(flag.index.impute ==0) %>% 
  arrange(patid,date.admit,date.dis) %>% 
  mutate(obs.time.since.MeV = ifelse(period == "after",as.numeric(difftime(date.admit,date.dis.measles,units="days")),0),# the duration from measles infection to the other admission
         time.dis.MeV.to.lastfup = as.numeric(difftime(time.at.lastfup,date.dis.measles,units="days")),
         hospitalization.measles = as.numeric(difftime(date.dis.measles,date.admit.measles,units="days"))) 
(n10<-length(unique(data.tmp.after.tmp$patid)))

data.tmp.after %<>%  arrange(patid,date.admit,date.dis) %>% 
  mutate(lag.patid  = lag(patid,1),
         # obs.time.since.MeV = ifelse(obs.time.since.MeV == 0 & flag.index.impute == 0,1,obs.time.since.MeV),
         lag.obs.time.since.MeV = lag(obs.time.since.MeV,1),
         time.start = ifelse(lag.patid == patid,lag.obs.time.since.MeV+lag(hospitalization,1),0))
data.tmp.after$time.start[1] <- 0
(n10<-length(unique(data.tmp.after.tmp$patid)))
data.tmp.after %<>% mutate(time.stop   = ifelse(flag.index.impute == 1,time.dis.MeV.to.lastfup,obs.time.since.MeV),
                               status   = ifelse(flag.index.impute == 1 ,0,1),
                           time.start  = ifelse(flag.index.impute==0,time.stop + hospitalization,0)) 
(n10<-length(unique(data.tmp.after.tmp$patid)))



