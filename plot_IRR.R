## Derive the data
data.total.main.1 <- data.transformed(data.raw,option = 1,emigration.rate = 0,seed =1234)
(age.max<-max(data.total.main.1$age.at.infection))
(time.max<-max(data.total.main.1$time.after.mv))


##for each set of knots we compute the datset and then merge them and plot them together.

#set knots 0
knots.time.af.mv<-c(1/12,2/12,4/12,6/12,8/12,1,1.5,3)

fit0<- glm(status~ splines::ns(age.at.infection,df=4)+I(time.after.mv>=14/365)+
             splines::ns(pmax(time.after.mv,14/365),knots=knots.time.af.mv,Boundary.knots = c(14/365,last.knots)),family=poisson(link = "log"),data = data.total.main.1,
           offset=offset_time)
summary(fit0)
dat0<-data.frame(time.after.mv=seq(14/365,10,by=0.01),age.at.infection=0,offset_time=1)
a<-splines::ns(data.total.main.1$age.at.infection,df=4)
(knots.age<-attr(a,"knots"))
mm<-model.matrix(~splines::ns(age.at.infection,knots=knots.age,Boundary.knots = c(0,age.max))+I(time.after.mv>=14/365)+
                   splines::ns(pmax(time.after.mv,14/365),knots=knots.time.af.mv,Boundary.knots = c(14/365,last.knots)),dat0)
pvar1 <- diag(mm[,-c(1:5)] %*% tcrossprod(vcov(fit0)[-c(1:5),-c(1:5)],mm[,-c(1:5)]))
dat0$pred<-mm[,-c(1:5)]%*%coef(fit0)[-c(1:5)] #predict(m,newdat,re.form=NA) would give the same results
dat0$plo<-dat0$pred-qnorm(0.975,0,1)*sqrt(pvar1)
dat0$phi = dat0$pred+qnorm(0.975,0,1)*sqrt(pvar1)
dat0$knots<- "c(1/12,2/12,4/12,6/12,8/12,1,1.5,3)"




#set knots 1
knots.time.af.mv<-c(1/12,2/12,5/12,8/12,1,1.5,3)
knots.age<-c(1/2,1,1.5,2,2.5,3,5)

library(splines)
fit1<- glm(status~ ns(age.at.infection,knots=knots.age,Boundary.knots = c(0,age.max))+I(time.after.mv>=14/365)+
             ns(pmax(time.after.mv-14/365,0),knots=knots.time.af.mv-14/365,Boundary.knots = c(0,last.knots-14/365)),family=poisson(link = "log"),data = data.total.main,
           offset=offset_time)
dat1<-data.frame(time.after.mv=seq(2/52,10,by=0.01),age.at.infection=0,offset_time=1)
dat1<-data.frame(time.after.mv=seq(0,10,by=0.01),age.at.infection=0,offset_time=1)
mm<-model.matrix(~ns(age.at.infection,knots=knots.age,Boundary.knots = c(0,age.max))+I(time.after.mv>=14/365)+
                   ns(pmax(time.after.mv-14/365,0),knots=knots.time.af.mv-14/365,Boundary.knots = c(0,last.knots-14/365)),dat1)
pvar1 <- diag(mm[,-c(1:9)] %*% tcrossprod(vcov(fit1)[-c(1:9),-c(1:9)],mm[,-c(1:9)]))
dat1$pred<-mm[,-c(1:9)]%*%coef(fit1)[-c(1:9)] #predict(m,newdat,re.form=NA) would give the same results
dat1$plo<-dat1$pred-qnorm(0.975,0,1)*sqrt(pvar1)
dat1$phi = dat1$pred+qnorm(0.975,0,1)*sqrt(pvar1)
dat1$knots<- "c(1/12,2/12,5/12,8/12,1,1.5,3)"

a<-predict(fit1,newdata=dat1)-predict(fit1,newdata=data.frame(time.after.mv=-1,age.at.infection=0,offset_time=1))


ggplot(dat1,aes(x=time.after.mv, y=pred)) + geom_line(size=.75)+
  geom_ribbon(data=dat1,aes(ymin=plo, ymax=phi, x=time.after.mv), alpha = 0.1)+
  scale_y_continuous("IRR relative to two-year before MeV",breaks = log(c(1/20,1/8,1/3,1,3,8,50)),labels=c("1/20","1/8","1/3",1,3,8,50))+
  scale_x_continuous("Time after measles infection [years]",limits = c(0,10),breaks=c(0,1:10),labels=c("0",1:10))+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", colour="black", size=12, margin = margin(t = 20,r = 20,b =20, l = 20)),
        axis.title.y = element_text(face="bold", colour="black", size=12, margin = margin(t = 20,r = 20,b = 20, l = 20)),
        axis.text.y = element_text(face="bold", colour="black", size=12),
        axis.text.x = element_text(face="bold", colour="black", size=12),
        axis.ticks.x=element_blank(),
        legend.text = element_text(face="bold",size = 8),
        legend.title = element_blank(),
        legend.key = element_rect(colour = NA),
        #legend.position = "right",
        strip.text=element_text(face = "bold",size=12))+coord_cartesian(ylim=c(-3,4))  


### set knots 2


knots.time.af.mv<-c(2/12,5/12,8/12,1,1.5,3)
fit2<- glm(status~ splines::ns(age.at.infection,df=4)+I(time.after.mv>=14/365)+
             splines::ns(pmax(time.after.mv,14/365),knots=knots.time.af.mv,Boundary.knots = c(14/365,last.knots)),family=poisson(link = "log"),data = data.total.main.1,
           offset=offset_time)
summary(fit2)
dat2<-data.frame(time.after.mv=seq(2/52,10,by=0.01),age.at.infection=0,offset_time=1)
a<-splines::ns(data.total.main.1$age.at.infection,df=4)
(knots.age<-attr(a,"knots"))
mm<-model.matrix(~splines::ns(age.at.infection,knots=knots.age,Boundary.knots = c(0,age.max))+I(time.after.mv>=14/365)+
                   splines::ns(pmax(time.after.mv,14/365),knots=knots.time.af.mv,Boundary.knots = c(14/365,last.knots)),dat2)
pvar1 <- diag(mm[,-c(1:5)] %*% tcrossprod(vcov(fit2)[-c(1:5),-c(1:5)],mm[,-c(1:5)]))
dat2$pred<-mm[,-c(1:5)]%*%coef(fit2)[-c(1:5)] #predict(m,newdat,re.form=NA) would give the same results
dat2$plo<-dat2$pred-qnorm(0.975,0,1)*sqrt(pvar1)
dat2$phi = dat2$pred+qnorm(0.975,0,1)*sqrt(pvar1)
dat2$knots<- "c(2/12,5/12,8/12,1,1.5,3)"



### set knots 3


knots.time.af.mv<-c(2/12,1,1.5,3)
fit3<- glm(status~ splines::ns(age.at.infection,df=4)+I(time.after.mv>=14/365)+
             splines::ns(pmax(time.after.mv,14/365),knots=knots.time.af.mv,Boundary.knots = c(14/365,last.knots)),family=poisson(link = "log"),data = data.total.main.1,
           offset=offset_time)
summary(fit3)
dat3<-data.frame(time.after.mv=seq(2/52,10,by=0.01),age.at.infection=0,offset_time=1)
a<-splines::ns(data.total.main.1$age.at.infection,df=4)
(knots.age<-attr(a,"knots"))
mm<-model.matrix(~splines::ns(age.at.infection,knots=knots.age,Boundary.knots = c(0,age.max))+I(time.after.mv>=14/365)+
                   splines::ns(pmax(time.after.mv,14/365),knots=knots.time.af.mv,Boundary.knots = c(14/365,last.knots)),dat3)
pvar1 <- diag(mm[,-c(1:5)] %*% tcrossprod(vcov(fit3)[-c(1:5),-c(1:5)],mm[,-c(1:5)]))
dat3$pred<-mm[,-c(1:5)]%*%coef(fit3)[-c(1:5)] #predict(m,newdat,re.form=NA) would give the same results
dat3$plo<-dat3$pred-qnorm(0.975,0,1)*sqrt(pvar1)
dat3$phi = dat3$pred+qnorm(0.975,0,1)*sqrt(pvar1)


dat3$knots<-"c(2/12,1,1.5,3)"

### set knots 4


knots.time.af.mv<-c(1,1.5,3)
fit4<- glm(status~ splines::ns(age.at.infection,df=4)+I(time.after.mv>=14/365)+
             splines::ns(pmax(time.after.mv,14/365),knots=knots.time.af.mv,Boundary.knots = c(14/365,last.knots)),family=poisson(link = "log"),data = data.total.main.1,
           offset=offset_time)
summary(fit4)
dat4<-data.frame(time.after.mv=seq(2/52,10,by=0.01),age.at.infection=0,offset_time=1)
a<-splines::ns(data.total.main.1$age.at.infection,df=4)
(knots.age<-attr(a,"knots"))
mm<-model.matrix(~splines::ns(age.at.infection,knots=knots.age,Boundary.knots = c(0,age.max))+I(time.after.mv>=14/365)+
                   splines::ns(pmax(time.after.mv,14/365),knots=knots.time.af.mv,Boundary.knots = c(14/365,last.knots)),dat4)
pvar1 <- diag(mm[,-c(1:5)] %*% tcrossprod(vcov(fit4)[-c(1:5),-c(1:5)],mm[,-c(1:5)]))
dat4$pred<-mm[,-c(1:5)]%*%coef(fit4)[-c(1:5)] #predict(m,newdat,re.form=NA) would give the same results
dat4$plo<-dat4$pred-qnorm(0.975,0,1)*sqrt(pvar1)
dat4$phi = dat4$pred+qnorm(0.975,0,1)*sqrt(pvar1)


dat4$knots<-"c(1,1.5,3)"


dat<-rbind(dat4,dat3,dat2,dat1,dat0)

p1<-ggplot(dat,aes(x=time.after.mv, y=pred,group=knots)) + geom_line(size=.75,aes(colour=knots))+
  geom_ribbon(data=dat,aes(ymin=plo, ymax=phi, x=time.after.mv, fill =knots), alpha = 0.1)+
  scale_y_continuous("IRR relative to two-year before MeV",breaks = log(c(1/20,1/8,1/3,1,3,8,50)),labels=c("1/20","1/8","1/3",1,3,8,50))+
  scale_x_continuous("Time after measles infection [years]",limits = c(0,10),breaks=c(0,1:10),labels=c("0",1:10))+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", colour="black", size=12, margin = margin(t = 20,r = 20,b =20, l = 20)),
        axis.title.y = element_text(face="bold", colour="black", size=12, margin = margin(t = 20,r = 20,b = 20, l = 20)),
        axis.text.y = element_text(face="bold", colour="black", size=12),
        axis.text.x = element_text(face="bold", colour="black", size=12),
        axis.ticks.x=element_blank(),
        legend.text = element_text(face="bold",size = 8),
        legend.title = element_blank(),
        legend.key = element_rect(colour = NA),
        #legend.position = "right",
        strip.text=element_text(face = "bold",size=12))+coord_cartesian(ylim=c(-3,4))  
p2<-ggplot(dat,aes(x=time.after.mv, y=pred,group=knots)) + geom_line(size=.75,aes(colour=knots))+facet_wrap(~knots)+
  geom_ribbon(data=dat,aes(ymin=plo, ymax=phi, x=time.after.mv, fill =knots), alpha = 0.4)+
  scale_y_continuous("IRR relative to two-year before MeV",breaks = log(c(1/20,1/8,1/3,1,3,8,50)),labels=c("1/20","1/8","1/3",1,3,8,50))+
  scale_x_continuous("Time after measles infection [years]",limits = c(0,10),breaks=c(0,1:10),labels=c("0",1:10))+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", colour="black", size=12, margin = margin(t = 20,r = 20,b =20, l = 20)),
        axis.title.y = element_text(face="bold", colour="black", size=12, margin = margin(t = 20,r = 20,b = 20, l = 20)),
        axis.text.y = element_text(face="bold", colour="black", size=12),
        axis.text.x = element_text(face="bold", colour="black", size=12),
        axis.ticks.x=element_blank(),
        legend.text = element_text(face="bold",size = 8),
        legend.title = element_blank(),
        legend.key = element_rect(colour = NA),
        #legend.position = "right",
        strip.text=element_text(face = "bold",size=12))+coord_cartesian(ylim=c(-3,4))  

p1




### Extra analysis
```{r eval=FALSE, include=FALSE}
### The Model 2: we randomly selected children to have emigrated from the group that didn't have later hospital admission observed, using an emigration rate of 8.3 per 1000 people per year (14). We censored them at their imputed time of emigration

###----Run the model----###
## We transform/load the data
# start_time <- Sys.time()
# data.total.sens.model2 <- data.transformed(data.raw,option = 2,emigration.rate = 0.00803,seed = 1234)
# end_time <- Sys.time()
# save(data.total.sens.model2,file = "Data/data_total_sens_model2.Rdata")
load(file = "Data/data_total_sens_model2.Rdata")
# 
## Compute the age max and time max of the data 

# (age.max<-max(data.total.sens.model2$age.at.infection))
# (time.max<-max(data.total.sens.model2$time.after.mv))

## We fit the Poisson regression mixed effect model with patient is a random intercept, the current age (age.at.infection) and time after MeV (time.after.mv) 
## are the main covariates.

start_time <- Sys.time()
fit.glmer<- glmer(status~ splines::ns(age.at.infection,df=4)+I(time.after.mv>=14/365)+
                    splines::ns(pmax(time.after.mv,0),knots=knots.time.af.mv,Boundary.knots = c(2/52,last.knots))+
                    (1|patid),family=poisson(link = "log"),data = data.total.sens.model2,
                  offset = offset_time)
end_time <- Sys.time()
# 
# (time_consumed1<-end_time - start_time)
# save(fit.stan_glmer,file="Results/fit_glmer_sens_model2.Rdata")

###---- Read the results from the models----###

# load(file ="Results/fit_glmer_sens_model2.Rdata")
# T.immuno_suppresion <- MCMC_estimate(fit.stan_glmer,interval = c(0,3),knots = knots.time.af.mv,last.knots = last.knots)
# CI.immuno_suppresion <- paste(round(quantile(T.immuno_suppresion,0.025),2),round(quantile(T.immuno_suppresion,0.975),2),sep=",")
# T.estim.immuno_suppresion_model2 <- paste(paste(round(quantile(T.immuno_suppresion,0.5),2),CI.immuno_suppresion,sep=" ("),")",sep="")
# save(T.estim.immuno_suppresion_model2,file = "Results/T.estim.immuno_suppresion_sens_model2.Rdata")
# 
# T.protective_immunity_value <- MCMC_estimate(fit.stan_glmer,interval = c(3,1e8),knots = knots.time.af.mv,last.knots = last.knots)
# T.protective_immunity <- T.protective_immunity_value-T.immuno_suppresion
# CI.protective_immunity <- paste(round(quantile(T.protective_immunity,0.025),2),round(quantile(T.protective_immunity,0.975),2),sep=",")
# T.estim.protective_immunity_model2 <- paste(paste(round(quantile(T.protective_immunity,0.5),2),CI.protective_immunity,sep=" ("),")",sep="")
# save(T.estim.protective_immunity_model2,file = "Results/T.estim.protective_immunity_sens_model2.Rdata")

```




```{r}

data.total.main.1 <- data.transformed(data.raw,option = 0,emigration.rate = 0,seed =1234)
(age.max<-max(data.total.main.1$age.at.infection))
(time.max<-max(data.total.main.1$time.after.mv))
knots.time.af.mv<-c(2/52,1/12,1.5,3)
fit<- glm(status~ splines::ns(age.at.infection,df=4)+I(time.after.mv>=14/365)+
            splines::ns(pmax(time.after.mv,0),knots=knots.time.af.mv,Boundary.knots = c(0,last.knots)),family=poisson(link = "log"),data = data.total.main.1,
          offset=offset_time)
summary(fit)
dat2<-data.frame(time.after.mv=seq(2/52,10,by=0.01),age.at.infection=0)
a<-splines::ns(data.total.main$age.at.infection,df=4)
(knots.age<-attr(a,"knots"))
mm<-model.matrix(~splines::ns(age.at.infection,knots=knots.age,Boundary.knots = c(0,age.max))+I(time.after.mv>=14/365)+
                   splines::ns(pmax(time.after.mv,0),knots=knots.time.af.mv,Boundary.knots = c(0,last.knots)),dat2)
pvar1 <- diag(mm[,-c(1:5)] %*% tcrossprod(vcov(fit)[-c(1:5),-c(1:5)],mm[,-c(1:5)]))
dat2$pred<-mm[,colnames(mm)[-c(1:5)]]%*%coef(fit)[-c(1:5)] #predict(m,newdat,re.form=NA) would give the same results
dat2$plo<-dat2$pred-qnorm(0.975,0,1)*sqrt(pvar1)
dat2$phi = dat2$pred+qnorm(0.975,0,1)*sqrt(pvar1)
dat2$removed="2 week removed"


data.total.main.2 <- data.transformed(data.raw,option = 0,emigration.rate = 0,seed =1234)
(age.max<-max(data.total.main.2$age.at.infection))
(time.max<-max(data.total.main.2$time.after.mv))
fit<- glm(status~ splines::ns(age.at.infection,df=4)+I(time.after.mv>0)+
            splines::ns(pmax(time.after.mv,0),knots=knots.time.af.mv,Boundary.knots = c(0,last.knots)),family=poisson(link = "log"),data = data.total.main.2,
          offset=offset_time)
summary(fit)
dat1<-data.frame(time.after.mv=seq(0.01,10,by=0.01),age.at.infection=0)
a<-splines::ns(data.total.main.2$age.at.infection,df=4)
(knots.age<-attr(a,"knots"))
mm<-model.matrix(~splines::ns(age.at.infection,knots=knots.age,Boundary.knots = c(0,age.max))+I(time.after.mv>0)+
                   splines::ns(pmax(time.after.mv,0),knots=knots.time.af.mv,Boundary.knots = c(0,last.knots)),dat1)
pvar1 <- diag(mm[,-c(1:5)] %*% tcrossprod(vcov(fit)[-c(1:5),-c(1:5)],mm[,-c(1:5)]))
dat1$pred<-mm[,colnames(mm)[-c(1:5)]]%*%coef(fit)[-c(1:5)] #predict(m,newdat,re.form=NA) would give the same results
dat1$plo<-dat1$pred-qnorm(0.975,0,1)*sqrt(pvar1)
dat1$phi = dat1$pred+qnorm(0.975,0,1)*sqrt(pvar1)
dat1$removed="no removed"


dat<-rbind(dat1,dat2)
dat$removed<-factor(dat$removed)
```



```{r fig1,echo=FALSE, out.width='100%',fig.cap="Incidence rate ratio of hospital admision due to non-measles ID pre  vs. post-measles",fig.height=7, fig.width=9}
#save(dat,file="Results/dat.Rdata")


p<-ggplot(dat,aes(x=time.after.mv, y=pred,group=removed)) + geom_line(size=1,aes(colour=removed))+
  geom_ribbon(data=dat,aes(ymin=plo, ymax=phi, x=time.after.mv, fill =removed), alpha = 0.3)+
  scale_y_continuous("IRR relative to two-year before MeV",breaks = log(c(1/20,1/8,1/3,1,3,8)),labels=c("1/20","1/8","1/3",1,3,8))+
  scale_x_continuous("Time after measles infection [years]",limits = c(0,10),breaks=c(0,1:10),labels=c("0",1:10))+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", colour="black", size=20, margin = margin(t = 20,r = 20,b =20, l = 20)),
        axis.title.y = element_text(face="bold", colour="black", size=20, margin = margin(t = 20,r = 20,b = 20, l = 20)),
        axis.text.y = element_text(face="bold", colour="black", size=20),
        axis.text.x = element_text(face="bold", colour="black", size=20),
        axis.ticks.x=element_blank(),
        legend.text = element_text(face="bold",size = 18),
        legend.title = element_blank(),
        legend.key = element_rect(colour = NA),
        #legend.position = "right",
        strip.text=element_text(face = "bold",size=20))+coord_cartesian(ylim=c(-3,2))  
```




```{r}
dat<-data.transformed(data.raw,option = 1,emigration.rate = 0,seed =1234)#using the non-removed data.
dat.new<-subset(dat,time.after.mv<=2/52)
fit<-glmer(status~I(period=="after")+(1|patid),family=poisson(link ="log"),data=dat.new,offset=offset_time)
tab<-summary(fit)$coefficients
Estimate<-paste(tab[2,1],tab[2,1]-qnorm(0.975,0,1)*tab[2,2],tab[2,1]+qnorm(0.975,0,1)*tab[2,2],sep=";")
```


```{r}
dat<-data.transformed(data.raw,option = 1,emigration.rate = 0,seed =1234)#using the non-removed data.
dat.new<-subset(dat,time.after.mv<=2/52)
fit<-glmer(status~I(period=="after")+(1|patid),family=poisson(link ="log"),data=dat.new,offset=offset_time)
tab<-summary(fit)$coefficients
Estimate<-paste(tab[2,1],tab[2,1]-qnorm(0.975,0,1)*tab[2,2],tab[2,1]+qnorm(0.975,0,1)*tab[2,2],sep=";")
```

### test

```{r}
data.total.main.1 <- data.transformed(data.raw,option = 1,emigration.rate = 0,seed =1234)
(age.max<-max(data.total.main.1$age.at.infection))
(time.max<-max(data.total.main.1$time.after.mv))
```

```{r}
#set knots 0
knots.time.af.mv<-c(1/12,2/12,4/12,6/12,8/12,1,1.5,3)

fit0<- glm(status~ splines::ns(age.at.infection,df=4)+I(time.after.mv>=14/365)+
             splines::ns(pmax(time.after.mv,0),knots=knots.time.af.mv,Boundary.knots = c(2/52,last.knots)),family=poisson(link = "log"),data = data.total.main.1,
           offset=offset_time)
summary(fit0)
dat0<-data.frame(time.after.mv=seq(2/52,10,by=0.01),age.at.infection=0,offset_time=1)
a<-splines::ns(data.total.main.1$age.at.infection,df=4)
(knots.age<-attr(a,"knots"))
mm<-model.matrix(~splines::ns(age.at.infection,knots=knots.age,Boundary.knots = c(0,age.max))+I(time.after.mv>=14/365)+
                   splines::ns(pmax(time.after.mv,0),knots=knots.time.af.mv,Boundary.knots = c(2/52,last.knots)),dat0)
pvar1 <- diag(mm[,-c(1:5)] %*% tcrossprod(vcov(fit0)[-c(1:5),-c(1:5)],mm[,-c(1:5)]))
dat0$pred<-mm[,-c(1:5)]%*%coef(fit0)[-c(1:5)] #predict(m,newdat,re.form=NA) would give the same results
dat0$plo<-dat0$pred-qnorm(0.975,0,1)*sqrt(pvar1)
dat0$phi = dat0$pred+qnorm(0.975,0,1)*sqrt(pvar1)
dat0$knots<- "c(1/12,2/12,4/12,6/12,8/12,1,1.5,3)"


#set knots 1
knots.time.af.mv<-c(1/12,2/12,5/12,8/12,1,1.5,3)

fit1<- glm(status~ splines::ns(age.at.infection,df=4)+I(time.after.mv>=14/365)+
             splines::ns(pmax(time.after.mv,0),knots=knots.time.af.mv,Boundary.knots = c(2/52,last.knots)),family=poisson(link = "log"),data = data.total.main.1,
           offset=offset_time)
summary(fit1)
dat1<-data.frame(time.after.mv=seq(2/52,10,by=0.01),age.at.infection=0,offset_time=1)
a<-splines::ns(data.total.main.1$age.at.infection,df=4)
(knots.age<-attr(a,"knots"))
mm<-model.matrix(~splines::ns(age.at.infection,knots=knots.age,Boundary.knots = c(0,age.max))+I(time.after.mv>=14/365)+
                   splines::ns(pmax(time.after.mv,0),knots=knots.time.af.mv,Boundary.knots = c(2/52,last.knots)),dat1)
pvar1 <- diag(mm[,-c(1:5)] %*% tcrossprod(vcov(fit1)[-c(1:5),-c(1:5)],mm[,-c(1:5)]))
dat1$pred<-mm[,-c(1:5)]%*%coef(fit1)[-c(1:5)] #predict(m,newdat,re.form=NA) would give the same results
dat1$plo<-dat1$pred-qnorm(0.975,0,1)*sqrt(pvar1)
dat1$phi = dat1$pred+qnorm(0.975,0,1)*sqrt(pvar1)
dat1$knots<- "c(1/12,2/12,5/12,8/12,1,1.5,3)"

### set knots 2


knots.time.af.mv<-c(2/12,5/12,8/12,1,1.5,3)
fit2<- glm(status~ splines::ns(age.at.infection,df=4)+I(time.after.mv>=14/365)+
             splines::ns(pmax(time.after.mv,0),knots=knots.time.af.mv,Boundary.knots = c(2/52,last.knots)),family=poisson(link = "log"),data = data.total.main.1,
           offset=offset_time)
summary(fit2)
dat2<-data.frame(time.after.mv=seq(2/52,10,by=0.01),age.at.infection=0,offset_time=1)
a<-splines::ns(data.total.main.1$age.at.infection,df=4)
(knots.age<-attr(a,"knots"))
mm<-model.matrix(~splines::ns(age.at.infection,knots=knots.age,Boundary.knots = c(0,age.max))+I(time.after.mv>=14/365)+
                   splines::ns(pmax(time.after.mv,0),knots=knots.time.af.mv,Boundary.knots = c(2/52,last.knots)),dat2)
pvar1 <- diag(mm[,-c(1:5)] %*% tcrossprod(vcov(fit2)[-c(1:5),-c(1:5)],mm[,-c(1:5)]))
dat2$pred<-mm[,-c(1:5)]%*%coef(fit2)[-c(1:5)] #predict(m,newdat,re.form=NA) would give the same results
dat2$plo<-dat2$pred-qnorm(0.975,0,1)*sqrt(pvar1)
dat2$phi = dat2$pred+qnorm(0.975,0,1)*sqrt(pvar1)
dat2$knots<- "c(2/12,5/12,8/12,1,1.5,3)"



### set knots 3


knots.time.af.mv<-c(2/12,1,1.5,3)
fit3<- glm(status~ splines::ns(age.at.infection,df=4)+I(time.after.mv>=14/365)+
             splines::ns(pmax(time.after.mv,0),knots=knots.time.af.mv,Boundary.knots = c(2/52,last.knots)),family=poisson(link = "log"),data = data.total.main.1,
           offset=offset_time)
summary(fit3)
dat3<-data.frame(time.after.mv=seq(2/52,10,by=0.01),age.at.infection=0,offset_time=1)
a<-splines::ns(data.total.main.1$age.at.infection,df=4)
(knots.age<-attr(a,"knots"))
mm<-model.matrix(~splines::ns(age.at.infection,knots=knots.age,Boundary.knots = c(0,age.max))+I(time.after.mv>=14/365)+
                   splines::ns(pmax(time.after.mv,0),knots=knots.time.af.mv,Boundary.knots = c(2/52,last.knots)),dat3)
pvar1 <- diag(mm[,-c(1:5)] %*% tcrossprod(vcov(fit3)[-c(1:5),-c(1:5)],mm[,-c(1:5)]))
dat3$pred<-mm[,-c(1:5)]%*%coef(fit3)[-c(1:5)] #predict(m,newdat,re.form=NA) would give the same results
dat3$plo<-dat3$pred-qnorm(0.975,0,1)*sqrt(pvar1)
dat3$phi = dat3$pred+qnorm(0.975,0,1)*sqrt(pvar1)


dat3$knots<-"c(2/12,1,1.5,3)"

### set knots 4


knots.time.af.mv<-c(1,1.5,3)
fit4<- glm(status~ splines::ns(age.at.infection,df=4)+I(time.after.mv>=14/365)+
             splines::ns(pmax(time.after.mv,0),knots=knots.time.af.mv,Boundary.knots = c(2/52,last.knots)),family=poisson(link = "log"),data = data.total.main.1,
           offset=offset_time)
summary(fit4)
dat4<-data.frame(time.after.mv=seq(2/52,10,by=0.01),age.at.infection=0,offset_time=1)
a<-splines::ns(data.total.main.1$age.at.infection,df=4)
(knots.age<-attr(a,"knots"))
mm<-model.matrix(~splines::ns(age.at.infection,knots=knots.age,Boundary.knots = c(0,age.max))+I(time.after.mv>=14/365)+
                   splines::ns(pmax(time.after.mv,0),knots=knots.time.af.mv,Boundary.knots = c(2/52,last.knots)),dat4)
pvar1 <- diag(mm[,-c(1:5)] %*% tcrossprod(vcov(fit4)[-c(1:5),-c(1:5)],mm[,-c(1:5)]))
dat4$pred<-mm[,-c(1:5)]%*%coef(fit4)[-c(1:5)] #predict(m,newdat,re.form=NA) would give the same results
dat4$plo<-dat4$pred-qnorm(0.975,0,1)*sqrt(pvar1)
dat4$phi = dat4$pred+qnorm(0.975,0,1)*sqrt(pvar1)


dat4$knots<-"c(1,1.5,3)"


dat<-rbind(dat4,dat3,dat2,dat1,dat0)

p<-ggplot(dat,aes(x=time.after.mv, y=pred,group=knots)) + geom_line(size=1,aes(colour=knots))+
  geom_ribbon(data=dat,aes(ymin=plo, ymax=phi, x=time.after.mv, fill =knots), alpha = 0.3)+
  scale_y_continuous("IRR relative to two-year before MeV",breaks = log(c(1/20,1/8,1/3,1,3,8,50)),labels=c("1/20","1/8","1/3",1,3,8,50))+
  scale_x_continuous("Time after measles infection [years]",limits = c(0,10),breaks=c(0,1:10),labels=c("0",1:10))+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", colour="black", size=20, margin = margin(t = 20,r = 20,b =20, l = 20)),
        axis.title.y = element_text(face="bold", colour="black", size=20, margin = margin(t = 20,r = 20,b = 20, l = 20)),
        axis.text.y = element_text(face="bold", colour="black", size=20),
        axis.text.x = element_text(face="bold", colour="black", size=20),
        axis.ticks.x=element_blank(),
        legend.text = element_text(face="bold",size = 18),
        legend.title = element_blank(),
        legend.key = element_rect(colour = NA),
        #legend.position = "right",
        strip.text=element_text(face = "bold",size=20))+coord_cartesian(ylim=c(-3,4))  
p
```

