library(ggplot2)
library(psych)
library(lsmeans)
library(plyr)
library(aod)
library(multcomp)
library(vcd)
library(ggplot2)
library(stringr)
library(aod)
library(lme4)
library(effects)

#data wrangling
setwd("C:/Users/rober/Dropbox/Grad School Research/Manuscripts/ZIKA EIP")
oct <- read.csv("combined_cq_results_master_sheet_october_2016_all_data_cleaned.csv", header=TRUE)
april <- read.csv("april_2017_cq_mastersheet_cleaned_for_R_all.csv", header=TRUE)
june <- read.csv("june_2017_cq_mastersheet_cleaned_for_R_full.csv", header=TRUE)
str(oct)
str(april)
str(june)
oct$temp <- as.factor(oct$temp)
mean(oct$positive)
oct$positive <- as.factor(oct$positive)
oct$org_num <- as.factor(oct$org_num)
oct$plate <- as.factor(oct$plate)
str(oct)

april$temp <- as.factor(april$temp)
april$positive <- as.factor(april$positive)
april$org_num <- as.factor(april$org_num)
april$plate <- as.factor(april$plate)
str(april)

june$temp <- as.factor(june$temp)
june$positive <- as.factor(june$positive)
june$org_num <- as.factor(june$org_num)
june$plate <- as.factor(june$plate)
str(june)


april$month <- as.factor(c("april"))
oct$month <- as.factor(c("oct"))
june$month <- as.factor(c("june"))
ap.oct.vec<-intersect(names(april), names(oct))
ap.jun.vec<-intersect(names(april), names(june))
ap.oc.jun<-rbind(april[, ap.oct.vec],oct[, ap.oct.vec], june[, ap.jun.vec])
ap.oc.jun$organism_num<-as.factor(ap.oc.jun$organism_num)

#subset by tissues
full.abd<- subset(ap.oc.jun, tissue=="Abd")
full.LW<- subset(ap.oc.jun, tissue=="LW")
full.sal <- subset(ap.oc.jun, tissue=="Sal")

#get counts of n-values per experiment

counts.month<-data.frame(table(ap.oc.jun$day, ap.oc.jun$temp, ap.oc.jun$month))
colnames(counts.month)<-c("day", "temp", "month", "number")
counts.month$temp<-sub("C", "", counts.month$temp)
counts.month$temp<-as.integer(paste(counts.month$temp))
counts.month$day<-as.integer(paste(counts.month$day))
str(counts.month)
#counts of positives per experiment
positive.counts<-data.frame(table(ap.oc.jun$day, ap.oc.jun$temp, ap.oc.jun$positive))
colnames(positive.counts)<-c("day", "temp", "positive", "frequency")


###fixed effects models
##fixed effects models

abd.logit.1<-glm(positive~temp*day, data=full.abd, family="binomial")
summary(abd.logit.1)
abd.predict <- expand.grid(temp=c("27C","30C","33C"), day=seq(3,14,1), positive=0)
pd<-predict(abd.logit.1, newdata=abd.predict, type="response", se.fit=TRUE)
abd.predict$output<-pd$fit
abd.predict$se<-pd$se.fit


LW.logit.1<-glm(positive~temp*day, data=full.LW, family="binomial")
summary(LW.logit.1)
LW.predict <- expand.grid(temp=c("27C","30C","33C"), day=seq(3,14,1), positive=0)
pd.lw<-predict(LW.logit.1, newdata=LW.predict, type="response", se.fit=TRUE)
LW.predict$output<-pd.lw$fit
LW.predict$se<-pd.lw$se.fit


sal.logit.1<-glm(positive~temp*day, data=full.sal, family="binomial")
summary(sal.logit.1)
sal.predict <- expand.grid(temp=c("27C","30C","33C"), day=seq(3,14,1), positive=0)
pd.sal<-predict(sal.logit.1, newdata=sal.predict, type="response", se.fit=TRUE)
sal.predict$output<-pd.sal$fit
sal.predict$se<-pd.sal$se.fit

###random effects, percent positive

abd.rand.1<-glmer(positive~temp*day+(1|month), data=full.abd, family="binomial")
summary(abd.rand.1)
cbind(coef(summary(abd.rand.1)),exp(cbind(OR=fixef(abd.rand.1),confint(abd.rand.1,parm="beta_",method="Wald"))))

AIC(abd.logit.1, abd.rand.1)

abd.effect<-Effect(c("temp", "day"), abd.rand.1)

ggplot(as.data.frame(abd.effect),
       aes(day,fit,colour=temp,fill=temp, linetype=temp))+
  geom_line()+
  ## colour=NA suppresses edges of the ribbon
  geom_ribbon(colour=NA,alpha=0.1,
              aes(ymin=lower,ymax=upper))+
  labs(title="Predicted probability of Infection",
       x ="Day", y = "Predicted probabilities (logit)")+
  theme(plot.title = element_text(hjust = 0.5))


LW.rand.1<-glmer(positive~temp*day+(1|month), data=full.LW, family="binomial")
summary(LW.rand.1)
cbind(coef(summary(LW.rand.1)),exp(cbind(OR=fixef(LW.rand.1),confint(LW.rand.1,parm="beta_",method="Wald"))))

AIC(LW.logit.1, LW.rand.1)

LW.effect<-Effect(c("temp", "day"), LW.rand.1)

ggplot(as.data.frame(LW.effect),
       aes(day,fit,colour=temp,fill=temp, linetype=temp))+
  geom_line()+
  ## colour=NA suppresses edges of the ribbon
  geom_ribbon(colour=NA,alpha=0.1,
              aes(ymin=lower,ymax=upper))+
  labs(title="Predicted probability of Dissemination",
       x ="Day", y = "Predicted probabilities (logit)")+
  scale_fill_grey()+
  scale_color_grey()+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))



sal.rand.1<-glmer(positive~temp*day+(1|month), data=full.sal, family="binomial")
summary(sal.rand.1)
cbind(coef(summary(sal.rand.1)),exp(cbind(OR=fixef(sal.rand.1),confint(sal.rand.1,parm="beta_",method="Wald"))))

AIC(sal.logit.1, sal.rand.1)

sal.effect<-Effect(c("temp", "day"), sal.rand.1)

ggplot(as.data.frame(sal.effect),
       aes(day,fit,colour=temp,fill=temp, linetype=temp))+
  geom_line()+
  geom_ribbon(colour=NA, alpha=0.1,
              aes(ymin=lower,ymax=upper))+
  labs(title="Predicted probability of Transmission",
       x ="Day", y = "Predicted probabilities (logit)")+
  scale_fill_grey()+
  scale_color_grey()+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))



#abdomen
positives.abd<-full.abd %>% group_by(temp, day) %>% 
  count(positive) %>% 
  mutate(total=sum(n)) %>% 
  filter(positive != 0) %>%
  mutate(percent = (prop.test(n, total, conf.level=0.95)$estimate)*100) %>%
  mutate(lowCI = (prop.test(n, total, conf.level=0.95, correct=FALSE)$conf.int[1])*100) %>%
  mutate(upCI = (prop.test(n, total, conf.level=0.95, correct=FALSE)$conf.int[2])*100) 
negatives.abd<-full.abd %>% group_by(temp, day) %>% 
  count(positive) %>% 
  mutate(total=sum(n)) %>% 
  filter(positive != 1) %>%
  mutate(percent = (prop.test(n, total, conf.level=0.95)$estimate)*100) %>%
  mutate(lowCI = (prop.test(n, total, conf.level=0.95, correct=FALSE)$conf.int[1])*100) %>%
  mutate(upCI = (prop.test(n, total, conf.level=0.95, correct=FALSE)$conf.int[2])*100)
percents.abd<-bind_rows(positives.abd, negatives.abd) %>%
  arrange(temp, day, desc(positive)) %>%
  mutate(tissue="infection")

###dissemination  
positives.LW<-full.LW %>% group_by(temp, day) %>% 
  count(positive) %>% 
  mutate(total=sum(n)) %>% 
  filter(positive != 0) %>%
  mutate(percent = (prop.test(n, total, conf.level=0.95)$estimate)*100) %>%
  mutate(lowCI = (prop.test(n, total, conf.level=0.95, correct=FALSE)$conf.int[1])*100) %>%
  mutate(upCI = (prop.test(n, total, conf.level=0.95, correct=FALSE)$conf.int[2])*100) 
negatives.LW<-full.LW %>% group_by(temp, day) %>% 
  count(positive) %>% 
  mutate(total=sum(n)) %>% 
  filter(positive != 1) %>%
  mutate(percent = (prop.test(n, total, conf.level=0.95)$estimate)*100) %>%
  mutate(lowCI = (prop.test(n, total, conf.level=0.95, correct=FALSE)$conf.int[1])*100) %>%
  mutate(upCI = (prop.test(n, total, conf.level=0.95, correct=FALSE)$conf.int[2])*100)
percents.LW<-bind_rows(positives.LW, negatives.LW) %>%
  arrange(temp, day, desc(positive)) %>%
  mutate(tissue= "dissemination")


##transmission
positives.sal<-full.sal %>% group_by(temp, day) %>% 
  count(positive) %>% 
  mutate(total=sum(n)) %>% 
  filter(positive != 0) %>%
  mutate(percent = (prop.test(n, total, conf.level=0.95)$estimate)*100) %>%
  mutate(lowCI = (prop.test(n, total, conf.level=0.95, correct=FALSE)$conf.int[1])*100) %>%
  mutate(upCI = (prop.test(n, total, conf.level=0.95, correct=FALSE)$conf.int[2])*100) 
negatives.sal<-full.sal %>% group_by(temp, day) %>% 
  count(positive) %>% 
  mutate(total=sum(n)) %>% 
  filter(positive != 1) %>%
  mutate(percent = (prop.test(n, total, conf.level=0.95)$estimate)*100) %>%
  mutate(lowCI = (prop.test(n, total, conf.level=0.95, correct=FALSE)$conf.int[1])*100) %>%
  mutate(upCI = (prop.test(n, total, conf.level=0.95, correct=FALSE)$conf.int[2])*100)
percents.sal<-bind_rows(positives.sal, negatives.sal) %>%
  arrange(temp, day, desc(positive)) %>%
  mutate(tissue ="transmission")

all.percents<-bind_rows(percents.abd,percents.LW, percents.sal) 
write.csv(all.percents, file="percent.infection.all.csv")
all.percents$day<-as.factor(all.percents$day)


#### percent graphs

full.abd$positive <-as.factor(full.abd$positive)
ggplot(full.abd, aes(temp)) +
  geom_bar(aes(fill=positive), position="fill")+
  labs(x="Temp", y="Percent Positive")+
  ggtitle("ZIKV Infection")+
  theme(plot.title = element_text(hjust = 0.5))+
  facet_grid(~day)+
  scale_fill_grey(name="PCR detection", breaks=c(0,1),
                  labels=c("Negative", "Positive"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))


full.LW$positive <-as.factor(full.LW$positive)
ggplot(full.LW, aes(temp)) +
  geom_bar(aes(fill=positive), position="fill")+
  labs(x="Temp", y="Percent Positive")+
  ggtitle("ZIKV Dissemination")+
  theme(plot.title = element_text(hjust = 0.5))+
  facet_grid(~day)+
  scale_fill_grey(name="PCR detection", breaks=c(0,1),
                  labels=c("Negative", "Positive"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))


full.sal$positive <-as.factor(full.sal$positive)
ggplot(full.sal, aes(temp)) +
  geom_bar(aes(fill=positive), position="fill")+
  labs(x="Temp", y="Percent Positive")+
  ggtitle("ZIKV Transmission")+
  theme(plot.title = element_text(hjust = 0.5))+
  facet_grid(~day)+
  scale_fill_grey(name="PCR detection", breaks=c(0,1),
                  labels=c("Negative", "Positive"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))


#diagnostic plots for false positives

plot(y = april$logSQ, x =april$cq , pch = 1,cex=0.75, col ="grey", xlim=c(20,45),
     ylab="loqSQ value", xlab="Ct value", main="April positives")
with(subset(april, april$positive==1 ), 
     points(y = logSQ, x =cq,
            pch=20, cex = 0.75))
with(subset(april,april$cq>37.1392 & april$positive==1), 
     points(y = logSQ, x =cq,
            pch=4, cex = 0.75))
abline(v=37.1392, col="black", lty=2, lwd=2)

plot(y = oct$logSQ, x =oct$cq , pch = 1,cex=0.75, col ="grey",
     ylab="loqSQ value", xlab="Ct value", main="October positives")
with(subset(oct, oct$positive==1 ), 
     points(y = logSQ, x =cq,
            pch=20, cex = 0.75))
with(subset(oct,oct$cq>33.92242 & oct$positive==1), 
     points(y = logSQ, x =cq,
            pch=4, cex = 0.75))
abline(v=33.92242, col="black", lty=2, lwd=2)

plot(y = june$logSQ, x =june$cq , pch = 1,cex=0.75, col ="grey",
     ylab="loqSQ value", xlab="Ct value", main="June positives")
with(subset(june, june$positive==1 ), 
     points(y = logSQ, x =cq,
            pch=20, cex = 0.75))
with(subset(june,june$cq>36.57538462 & june$positive==1), 
     points(y = logSQ, x =cq,
            pch=4, cex = 0.75))
abline(v=36.57538462, col="black", lty=2, lwd=2)


##fixed effects models

abd.logit.1<-glm(positive~temp*day, data=full.abd, family="binomial")
summary(abd.logit.1)
abd.predict <- expand.grid(temp=c("27C","30C","33C"), day=seq(3,14,1), positive=0)
pd<-predict(abd.logit.1, newdata=abd.predict, type="response", se.fit=TRUE)
abd.predict$output<-pd$fit
abd.predict$se<-pd$se.fit


LW.logit.1<-glm(positive~temp*day, data=full.LW, family="binomial")
summary(LW.logit.1)
LW.predict <- expand.grid(temp=c("27C","30C","33C"), day=seq(3,14,1), positive=0)
pd.lw<-predict(LW.logit.1, newdata=LW.predict, type="response", se.fit=TRUE)
LW.predict$output<-pd.lw$fit
LW.predict$se<-pd.lw$se.fit


sal.logit.1<-glm(positive~temp*day, data=full.sal, family="binomial")
summary(sal.logit.1)
sal.predict <- expand.grid(temp=c("27C","30C","33C"), day=seq(3,14,1), positive=0)
pd.sal<-predict(sal.logit.1, newdata=sal.predict, type="response", se.fit=TRUE)
sal.predict$output<-pd.sal$fit
sal.predict$se<-pd.sal$se.fit



#########logSQ

abd.lm.1<-lm(logSQ~temp*day, data=full.abd)
summary(abd.lm.1)
abd.predict.lm <- expand.grid(temp=c("27C","30C","33C"), day=seq(3,14,1), positive=0)
pd.abd.lm<-predict(abd.lm.1, newdata=abd.predict.lm, type="response", se.fit=TRUE)
abd.predict.lm$output<-pd.abd.lm$fit
abd.predict.lm$se<-pd.abd.lm$se.fit

ggplot(abd.predict.lm, aes(x=day, y=output, linetype=temp))+
  geom_line(aes(color=temp))+
  geom_ribbon(colour=NA,alpha=0.1,aes(ymin=(output-se), ymax=(output+se), fill=temp, alpha=0.1))+
  guides(alpha=FALSE, size=FALSE)+
  labs(title="Predicted logSQ, Abdomen",
       x ="Day", y = "Predicted logSQ")+
  scale_fill_grey()+
  scale_color_grey()+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))


LW.lm.1<-lm(logSQ~temp*day, data=full.LW)
summary(LW.lm.1)
LW.predict.lm <- expand.grid(temp=c("27C","30C","33C"), day=seq(3,14,1), positive=0)
pd.LW.lm<-predict(LW.lm.1, newdata=LW.predict.lm, type="response", se.fit=TRUE)
LW.predict.lm$output<-pd.LW.lm$fit
LW.predict.lm$se<-pd.LW.lm$se.fit

ggplot(LW.predict.lm, aes(x=day, y=output, linetype=temp))+
  geom_line(aes(color=temp))+
  geom_ribbon(colour=NA,alpha=0.1,aes(ymin=(output-se), ymax=(output+se), fill=temp, alpha=0.1))+
  guides(alpha=FALSE, size=FALSE)+
  labs(title="Predicted logSQ, LW",
       x ="Day", y = "Predicted logSQ")+
  scale_fill_grey()+
  scale_color_grey()+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))


sal.lm.1<-lm(logSQ~temp*day, data=full.sal)
summary(sal.lm.1)
sal.predict.lm <- expand.grid(temp=c("27C","30C","33C"), day=seq(3,14,1), positive=0)
pd.sal.lm<-predict(sal.lm.1, newdata=sal.predict.lm, type="response", se.fit=TRUE)
sal.predict.lm$output<-pd.sal.lm$fit
sal.predict.lm$se<-pd.sal.lm$se.fit

ggplot(sal.predict.lm, aes(x=day, y=output, linetype=temp))+
  geom_line(aes(color=temp))+
  geom_ribbon(colour=NA,alpha=0.1,aes(ymin=(output-se), ymax=(output+se), fill=temp, alpha=0.1))+
  guides(alpha=FALSE, size=FALSE)+
  labs(title="Predicted logSQ, Saliva",
       x ="Day", y = "Predicted logSQ")+
  scale_fill_grey()+
  scale_color_grey()+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))