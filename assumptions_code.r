# Purpose: Supporting R code for "The importance of making assumptions in bias analysis"
# Author: Richard MacLehose
# Date: 8/26/20
# 
  
#start with a clean slate
rm(list=ls())
library(triangle)
  
#User-defined functions necessary for  paper
# Ding and VanderWeele bounding factor
bf.dvw<-function(rr.ud,rr.eu){
  rr.ud*rr.eu/(rr.eu+rr.ud-1)
}
  
#Schlesselman bounding factor
bf.s<-function(rr.ud,rr.eu,p1){
  (rr.eu+(rr.ud-1)*rr.eu*p1)/(rr.eu+(rr.ud-1)*p1)
}
  
#e-value
eval<-function(x){
  x+sqrt(x^2-x)
}
  
# g- value
gval<-function(x,p1){
  ((x+p1+x*p1-1)+sqrt((1-p1-x-x*p1)^2-4*x*p1^2))/(2*p1)
}

#prevalence among the confounders based on prevalence overall (pc), prevalence of the exposure (pe) and effect of e on u (RR.eu)
pconf.est<-function(pc,pe,rr.eu){
  pc/(pe+(1-pe)/(rr.eu))
}
  
#Figure 1
#rr.ud=2; rr.eu=2; compare e-val and g-val
rr=2
p=seq(0,1,length=1000)
plot(p,bf.s(rr.ud=rr,rr.eu=rr,p),xlab="Prevalence of the Confounder in the Exposed", ylab="Bias Factor",ylim=c(1,1.5),type='l')
segments(x0=0,x1=1,y0=bf.dvw(rr,rr),y1=bf.dvw(rr,rr),lty=2)
  
  
  
#section: methods - bias and bounding factors
# example with RR.ud RR.eu =2
round(bf.dvw(rr.ud=2,rr.eu=2),1)
  
#Example 1
rr.obs=3.9
rr.ud=4
rr.eu=2
bf.dvw.ex1<-bf.dvw(rr.ud=rr.ud,rr.eu=rr.eu)
round(bf.dvw.ex1,2)


#from baros:
pc=.356
#from Victora
pe=.298
#assume
rr.eu=2
p1.est<-pconf.est(pc=pc,pe=pe,rr.eu=rr.eu)
round(p1.est,2)
bf.s.ex1<-bf.s(rr.ud=rr.ud,rr.eu=rr.eu,p1=p1.est)
round(bf.s.ex1,2)
round(rr.obs/bf.dvw.ex1,1)
round(rr.obs/bf.s.ex1,1)
  
  
#boundign factors vs PBA
#bias parameters
#from Baros (table 2): p1=.356 (.203-437)
p1.baros<-pconf.est(pc=.356,pe=pe,rr.eu=rr.eu)
p1.baros.up<-pconf.est(pc=.203,pe=pe,rr.eu=rr.eu)
p1.baros.low<-pconf.est(pc=.437,pe=pe,rr.eu=rr.eu)
round(c(p1.baros,p1.baros.low,p1.baros.up),2)
#smoking-breastfeeding: 4
#formula smoking: 2
rr.eu=2
rr.ud=4
#observed effect
rr.obs=3.9
#se from victoria paper
se.obs=(log(8.7)-log(1.8))/(2*1.96)
#check se formula
c(exp(log(rr.obs)-1.96*se.obs),exp(log(rr.obs)+1.96*se.obs))
  
  
  
#sample Pr(U=1|E=1)
p1.s<-rtriangle(10^6,.255,.717,.55)
#compute bias factor
bf.s.p1=(rr.eu+(rr.ud-1)*rr.eu*p1.s)/(rr.eu+(rr.ud-1)*p1.s)
#dist'n of bias factor
round(quantile(bf.s.p1,c(.025,.5,.975)),2)
  
#Figure 2
plot(1, xlim=c(1,1.7),ylim=c(0,10^5), type='n',xlab=expression('BF'[S]),ylab="Frequency")
abline(v=1.60,lwd=3)
hist(bf.s.p1,xlab="",main="",add=TRUE)
  
#adjusted RR
round(quantile(rr.obs/bf.s.p1,c(.025,.5,.975)),2)
  
#incorporate random errer
rr.adj.all<-exp(rnorm(10^6,log(rr.obs/bf.s.p1),se.obs))
round(quantile(rr.adj.all,c(.025,.5,.975)),2)
  
#e-values
#Figure 3
rr.obs=2
p=seq(0,1,length=1000)
plot(p,gval(rr.obs,p),xlab="Prevalence in the Exposed",ylab=expression(paste('G-value(p'[1],')')),type="l",lty=1, ylim=c(1,30),yaxt="n")
v1=c(1,5,10,15,20,25,30)
v2=c("1","5","10","15","20","25","30")
axis(side=2,at=v1,labels=v2)
eval(2)
gval(2,1)
  

#example 2
eval(1.33)
eval(1.16)
  
  
  
#Table
cbind(1:10/10,round(gval(1.33,1:10/10),1),round((1:10/10)/gval(1.33,1:10/10),2))
cbind(1:10/10,round(gval(1.16,1:10/10),1),round((1:10/10)/gval(1.16,1:10/10),2))
  
  
  
  
  
  
  
  
