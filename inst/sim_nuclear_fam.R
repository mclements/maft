#Created by Robert Szulkin, 2015-05-14
#
#This file creates families with random components: g2=additive genetic,c2=common environment,e2=random error
#The simulated outcome will be normally distributed with censoring
#

#Load required libraries
library(mvtnorm) #Random multivariate normally distributed data


#Assume mean follow-up time is five years
logT=log(5)

#Set variance components
g2 <- .3; c2 <- .2; e2 <- .5
#g2 <- 0; c2 <- 0; e2 <- 2


#Size of each simulated family structure
N=5000

############################################################################################
######Simulate families with one father and one son

set.seed(34345)
sigma_fa_son<-matrix(c(g2+c2+e2,0.5*g2,0.5*g2,g2+c2+e2),2,2)
logT_fa_son<-data.frame(rmvnorm(N,c(logT,logT),sigma_fa_son))

#Give family id
logT_fa_son$famid<-1:length(logT_fa_son[,1])

###################Make data in long format

#Create data for fathers
outcome_fa<-logT_fa_son[,c("X1","famid")]
names(outcome_fa)<-c("logT","famid")
outcome_fa$parent<-1
outcome_fa$child<-0

#Create data for sons
outcome_son<-logT_fa_son[,c("X2","famid")]
names(outcome_son)<-c("logT","famid")
outcome_son$parent<-0
outcome_son$child<-1

#Combine data and order
fa_son<-rbind(outcome_fa,outcome_son)
fa_son<-fa_son[order(fa_son$famid,-fa_son$parent),]

#Create family size variable (tnr) and number of children variable (knr2) 
fa_son$tnr<-2
fa_son$knr2<-1



############################################################################################
######Simulate families with two brothers

set.seed(35767)
sigma_brothers<-matrix(c(g2+c2+e2,0.5*g2+c2,0.5*g2+c2,g2+c2+e2),2,2)
logT_brothers<-data.frame(rmvnorm(N,c(logT,logT),sigma_brothers))

#Give family id
logT_brothers$famid<-(N+1):(N+length(logT_brothers[,1]))

###################Make data in long format

#Create data for brother 1
outcome_b1<-logT_brothers[,c("X1","famid")]
names(outcome_b1)<-c("logT","famid")
outcome_b1$parent<-0
outcome_b1$child<-1

#Create data for brother 2
outcome_b2<-logT_brothers[,c("X2","famid")]
names(outcome_b2)<-c("logT","famid")
outcome_b2$parent<-0
outcome_b2$child<-1

#Combine data and order
brothers<-rbind(outcome_b1,outcome_b2)
brothers<-brothers[order(brothers$famid,-brothers$parent),]

#Create family size variable (tnr) and number of children variable (knr2) 
brothers$tnr<-2
brothers$knr2<-2



############################################################################################
######Simulate families with 1 father and two brothers

set.seed(8895)
sigma_fa_brothers<-matrix(c(g2+c2+e2,0.5*g2,0.5*g2,0.5*g2,g2+c2+e2,0.5*g2+0.5*c2,0.5*g2,0.5*g2+0.5*c2,g2+c2+e2),3,3)
logT_fa_brothers<-data.frame(rmvnorm(round(N/3),c(logT,logT,logT),sigma_fa_brothers))

#Give family id
logT_fa_brothers$famid<-(2*N+1):(2*N+length(logT_fa_brothers[,1]))

###################Make data in long format

#Create data for father
outcome_father<-logT_fa_brothers[,c("X1","famid")]
names(outcome_father)<-c("logT","famid")
outcome_father$parent<-1
outcome_father$child<-0


#Create data for brother 1
outcome_b1<-logT_fa_brothers[,c("X2","famid")]
names(outcome_b1)<-c("logT","famid")
outcome_b1$parent<-0
outcome_b1$child<-1

#Create data for brother 2
outcome_b2<-logT_fa_brothers[,c("X3","famid")]
names(outcome_b2)<-c("logT","famid")
outcome_b2$parent<-0
outcome_b2$child<-1

#Combine data and order
fa_brothers<-rbind(outcome_father,outcome_b1,outcome_b2)
fa_brothers<-fa_brothers[order(fa_brothers$famid,-fa_brothers$parent),]

#Create family size variable (tnr) and number of children variable (knr2) 
fa_brothers$tnr<-3
fa_brothers$knr2<-2


############################################################################################
######Simulate families with three brothers

set.seed(694)
sigma_3brothers<-matrix(c(g2+c2+e2,0.5*g2+0.5*c2,0.5*g2+0.5*c2,0.5*g2+0.5*c2,g2+c2+e2,0.5*g2+0.5*c2,0.5*g2+0.5*c2,0.5*g2+0.5*c2,g2+c2+e2),3,3)
logT_3brothers<-data.frame(rmvnorm(round(N/3),c(logT,logT,logT),sigma_3brothers))

#Give family id
logT_3brothers$famid<-(3*N+1):(3*N+length(logT_3brothers[,1]))

###################Make data in long format

#Create data for brother 1
outcome_b1<-logT_3brothers[,c("X1","famid")]
names(outcome_b1)<-c("logT","famid")
outcome_b1$parent<-0
outcome_b1$child<-1

#Create data for brother 2
outcome_b2<-logT_3brothers[,c("X2","famid")]
names(outcome_b2)<-c("logT","famid")
outcome_b2$parent<-0
outcome_b2$child<-1

#Create data for brother 3
outcome_b3<-logT_3brothers[,c("X3","famid")]
names(outcome_b3)<-c("logT","famid")
outcome_b3$parent<-0
outcome_b3$child<-1

#Combine data and order
brothers3<-rbind(outcome_b1,outcome_b2,outcome_b3)
brothers3<-brothers3[order(brothers3$famid,-brothers3$parent),]

#Create family size variable (tnr) and number of children variable (knr2) 
brothers3$tnr<-3
brothers3$knr2<-3

#########################################################################################
#Put data togeter
famdat<-rbind(fa_son,brothers,fa_brothers,brothers3)

#Create event variable with no censoring
famdat$death_nocens<-1

#Simulate independent censoring times
set.seed(6249)
famdat$censtime<-rnorm(length(famdat$logT),mean=log(5),sd=1)

#Take the minimum of the observed times and censoring times
famdat$Y<-pmin(famdat$logT,famdat$censtime)

#Create event variable
famdat$event<-0
famdat[famdat$censtime>famdat$logT,"event"]<-1

#Create left truncation variable, entry time at day 1, i.e. at 1/365 years
famdat$left_entry<-log(1/365)

#Simulate age of diagnosis 
#Assume that those who experience an event have a mean age of diagnosis at 69
#Those who are censored are assumed to have a higher mean age at diagnosis,i.e. 
#die from something else
set.seed(4333)
famdat[famdat$event==1,"age"]<-
round(rnorm(length(which(famdat$event==1)),mean=69,sd=5)) 

set.seed(999)
famdat[famdat$event==0,"age"]<-
round(rnorm(length(which(famdat$event==0)),mean=69,sd=5)) 

#Create intercept variable
famdat$prev<-1

require(eha)
famdat2 <- transform(data.frame(famdat),t=exp(Y),
                     entry_time=exp(left_entry))
aftreg(Surv(entry_time,t,event)~I(age-60),data=famdat2,dist="lognormal")
summary(survreg(Surv(t,event)~I(age-60),data=famdat2,dist="lognormal"))

#####################################################################
#####Try data on Ben's code 
#Set work directory 
setwd("Z:/Study7_Heritability/Programs/Bens code") 
library(fampack)
source('functions-BY-11jan2010.R')

beta <- c(prev=1.6,age60=-0.002)


#No idea what the variance components starting values should be
var_comp<-c(g=1,C=1,e=1)

gaussSeidel <- Gauseidel(comp=var_comp,betas=beta,
                         model=logT~prev+age+r(g)+r(C),
                         fdat=as.matrix(famdat[1:10,]),stratum=c('tnr','knr2'),
                         subjects='famid', entry='left_entry', event='death_nocens', ascerp=1, iter=100)


############################################################################
###############Reduce functions 
source('Z:/Study7_Heritability/Programs/Bens_code_modified/functions_BY_modify2015_05_19.R')

famdat=transform(famdat,age60=age-60)
pc_fam<-famdat[,c("Y","prev","age60","famid","left_entry","event","tnr","knr2")]

gaussSeidel <- Gauseidel(comp=var_comp,betas=beta,
                         model=Y~prev+age60+r(g)+r(C),
                         fdat=as.matrix(pc_fam),stratum=c('tnr','knr2'),
                         subjects='famid', entry='left_entry', event='event', ascerp=1, iter=100)

gaussSeidel$var.obj
t(gaussSeidel$summary)

gaussSeidel_small <- Gauseidel(comp=var_comp,betas=beta,
                         model=logT~prev+age+r(g)+r(C),
                         fdat=as.matrix(famdat[1:10,]),stratum=c('tnr','knr2'),
                         subjects='famid', entry='left_entry', event='death_nocens', ascerp=1, iter=100)


#RS: Step1 is essentially step A1-A4 in Yip et al  

betaNew <- beta
opt1 <- list(par=var_comp)
for (i in 1:8) {
step1 <- UP.beta.nu(Y~prev+age60+r(g)+r(C),survdat=as.matrix(famdat),startval=betaNew,
                    subject='famid',family.type='nuclear',varcomp=opt1$par,
                    entry='left_entry', event='event',rel.err=1e-8,max.iter=100)
meldat <- as.matrix(famdat)
mel.strat <- cbind(meldat, U=step1$U, MU=step1$MU, YC=step1$YC, BC=step1$BC)
colnames(mel.strat) <- c(colnames(meldat),'U','MU','YC','BC')
mel.strat <- Stratify(mel.strat,c('tnr','knr2'))
opt1 <- optim(opt1$par,HP.like.nuclear,xdat=mel.strat, subject='famid',
              entry='left_entry', event='event', hessian=TRUE, control=list(reltol=1e-10,trace=1))
betaNew <- as.vector(step1$est.beta)
names(betaNew) <- names(beta)
}
opt1$par
betaNew

objective <- function(varcomp, xdat, subjects, entry, event, famtype, fix.obj, 
                       ascer.prob = 1, ret = c("adj.hlike"))
                      HP.like.nuclear(c(g=0.2,C=0.3,e=varcomp), xdat, subjects, entry, event, famtype, fix.obj, 
                                                     ascer.prob = ascer.prob, ret = ret)
opt1 <- optim(1,objective,xdat=mel.strat, subject='famid',
              entry='left_entry', event='event', hessian=TRUE, control=list(reltol=1e-10,trace=1),
              method="Brent", lower=0,upper=100)
opt1$par

objective(1,xdat=mel.strat, subject='famid', entry='left_entry', event='event')


VB <- Vbeta.nuclear(mel.strat, beta.obj=step1, optim.obj=opt1,entry='b')


