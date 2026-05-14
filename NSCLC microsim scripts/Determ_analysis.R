library(tidyverse)
library(ggplot2)

#s.df=read.csv(paste0(getwd(),"/Outputs/Individual State Trace.csv"))
#s.df.comp=read.csv(paste0(getwd(),"/Outputs/Individual State Trace Comp.csv"))
#s.df.comp2=read.csv(paste0(getwd(),"/Outputs/Individual State Trace Comp2.csv"))
#c.df=read.csv(paste0(getwd(),"/Outputs/Individual Cost Trace.csv"))
#c.df.comp=read.csv(paste0(getwd(),"/Outputs/Individual Cost Trace Comp.csv"))
#c.df.comp2=read.csv(paste0(getwd(),"/Outputs/Individual Cost Trace Comp2.csv"))
#u.df=read.csv(paste0(getwd(),"/Outputs/Individual Utility Trace.csv"))
#u.df.comp=read.csv(paste0(getwd(),"/Outputs/Individual Utility Trace Comp.csv"))
#u.df.comp2=read.csv(paste0(getwd(),"/Outputs/Individual Utility Trace Comp2.csv"))

s.df=read.csv(paste0(getwd(),"/Outputs/Individual State Trace LY.csv"))
s.df.comp=read.csv(paste0(getwd(),"/Outputs/Individual State Trace Comp LY.csv"))
s.df.comp2=read.csv(paste0(getwd(),"/Outputs/Individual State Trace Comp2 LY.csv"))
c.df=read.csv(paste0(getwd(),"/Outputs/Individual Cost Trace LY.csv"))
c.df.comp=read.csv(paste0(getwd(),"/Outputs/Individual Cost Trace Comp LY.csv"))
c.df.comp2=read.csv(paste0(getwd(),"/Outputs/Individual Cost Trace Comp2 LY.csv"))
u.df=read.csv(paste0(getwd(),"/Outputs/Individual Utility Trace LY.csv"))
u.df.comp=read.csv(paste0(getwd(),"/Outputs/Individual Utility Trace Comp LY.csv"))
u.df.comp2=read.csv(paste0(getwd(),"/Outputs/Individual Utility Trace Comp2 LY.csv"))

###################### First order Model Convergence ###################### 

costs=as.numeric(c.df[300,1:10000])
mean_vals=sapply(1:length(costs), function(i) mean(costs[1:i]))
plot=data.frame(n_individuals=1:length(costs),mean_cost=mean_vals)

cost.conv.plot=ggplot(plot, aes(x = n_individuals, y = mean_cost)) +
  geom_line(color="blue") +
  labs(
    x="Number of Individuals",
    y="Mean Cost",
    title="Mean Cost vs Number of Individuals"
  )+theme_minimal()
cost.conv.plot

QALY=as.numeric(u.df[300,1:10000])
mean_vals=sapply(1:length(QALY), function(i) mean(QALY[1:i]))
plot=data.frame(n_individuals = 1:length(QALY), mean_cost = mean_vals)

Util.conv.plot=ggplot(plot, aes(x = n_individuals, y = mean_cost)) +
  geom_line(color = "blue") +
  labs(
    x = "Number of Individuals",
    y = "Mean QALY per person",
    title = "Mean QALY per person vs Number of Individuals"
  )+theme_minimal()
Util.conv.plot

############################# Subgroup analysis base model #############################

table2=data.frame("Subgroup"=c("","Overall","Stage I","Stage II","Stage III","Stage IV","EGFR","ALK","PD-L1: >50%","PD-L1: <50%"),
                  "Mean (95% CI)"=c("Cost Per Case $CAD 2025","","","","","","","","",""),
                  "Mean (95% CI)"=c("QALY","","","","","","","","",""))

me.ov.c=mean(unlist(c.df[300,]))
me.ov.c.ci=t.test(unlist(c.df[300,]))$conf.int
table2[2,2]=paste0(prettyNum(round(me.ov.c,digits=0), big.mark = ",")," (",prettyNum(round(me.ov.c.ci[2],digits=0), big.mark = ","),"-",prettyNum(round(me.ov.c.ci[1],digits=0), big.mark = ","),")")

me.ov.u=mean(unlist(u.df[300,]))
me.ov.u.ci=t.test(unlist(u.df[300,]))$conf.int
table2[2,3]=paste0(round(me.ov.u,digits=2)," (",round(me.ov.u.ci[2],digits=2),"-",round(me.ov.u.ci[1],digits=2),")")

u.det.temp=u.df[,s.df[1,]=="S1"]

#"Stage I Overall"
me.s1.u=mean(unlist(u.det.temp[300,]))
me.s1.u.ci=t.test(unlist(u.det.temp[300,]))$conf.int
table2[3,3]=paste0(round(me.s1.u,digits=2)," (",round(me.s1.u.ci[2],digits=2),"-",round(me.s1.u.ci[1],digits=2),")")

u.det.temp=u.df[,s.df[1,]=="S2"]

#"Stage II Overall"
me.s2.u=mean(unlist(u.det.temp[300,]))
me.s2.u.ci=t.test(unlist(u.det.temp[300,]))$conf.int
table2[4,3]=paste0(round(me.s2.u,digits=2)," (",round(me.s2.u.ci[2],digits=2),"-",round(me.s2.u.ci[1],digits=2),")")

u.det.temp=u.df[,s.df[1,]=="S3"]

#"Stage III Overall"
me.s3.u=mean(unlist(u.det.temp[300,]))
me.s3.u.ci=t.test(unlist(u.det.temp[300,]))$conf.int
table2[5,3]=paste0(round(me.s3.u,digits=2)," (",round(me.s3.u.ci[2],digits=2),"-",round(me.s3.u.ci[1],digits=2),")")

c.det.temp=c.df[,s.df[1,]=="S1"]


#"Stage I Overall"
me.s1.c=mean(unlist(c.det.temp[300,]))
me.s1.c.ci=t.test(unlist(c.det.temp[300,]))$conf.int
table2[3,2]=paste0(prettyNum(round(me.s1.c,digits=0), big.mark = ",")," (",prettyNum(round(me.s1.c.ci[2],digits=0), big.mark = ","),"-",prettyNum(round(me.s1.c.ci[1],digits=0), big.mark = ","),")")

c.det.temp=c.df[,s.df[1,]=="S2"]
s.df.s1=s.df[,s.df[1,]=="S2"]

#"Stage II Overall"
me.s2.c=mean(unlist(c.det.temp[300,]))
me.s2.c.ci=t.test(unlist(c.det.temp[300,]))$conf.int
table2[4,2]=paste0(prettyNum(round(me.s2.c,digits=0), big.mark = ",")," (",prettyNum(round(me.s2.c.ci[2],digits=0), big.mark = ","),"-",prettyNum(round(me.s2.c.ci[1],digits=0), big.mark = ","),")")

c.det.temp=c.df[,s.df[1,]=="S3"]
s.df.s1=s.df[,s.df[1,]=="S3"]

#"Stage III Overall"
me.s3.c=mean(unlist(c.det.temp[300,]))
me.s3.c.ci=t.test(unlist(c.det.temp[300,]))$conf.int
table2[5,2]=paste0(prettyNum(round(me.s3.c,digits=0), big.mark = ",")," (",prettyNum(round(me.s3.c.ci[2],digits=0), big.mark = ","),"-",prettyNum(round(me.s3.c.ci[1],digits=0), big.mark = ","),")")

print("EGFR 1L")
c.egfr=c.df[,s.df[1,]=="EGFR 1L"]
u.egfr=u.df[,s.df[1,]=="EGFR 1L"]

me.egfr.u=mean(unlist(u.egfr[300,]))
me.egfr.u.ci=t.test(unlist(u.egfr[300,]))$conf.int
table2[7,3]=paste0(round(me.egfr.u,digits=2)," (",round(me.egfr.u.ci[2],digits=2),"-",round(me.egfr.u.ci[1],digits=2),")")

me.egfr.c=mean(unlist(c.egfr[300,]))
me.egfr.c.ci=t.test(unlist(c.egfr[300,]))$conf.int
table2[7,2]=paste0(prettyNum(round(me.egfr.c,digits=0), big.mark = ",")," (",prettyNum(round(me.egfr.c.ci[2],digits=0), big.mark = ","),"-",prettyNum(round(me.egfr.c.ci[1],digits=0), big.mark = ","),")")

print("ALK 1L")
c.alk=c.df[,s.df[1,]=="ALK 1L"]
u.alk=u.df[,s.df[1,]=="ALK 1L"]

me.alk.u=mean(unlist(u.alk[300,]))
me.alk.u.ci=t.test(unlist(u.alk[300,]))$conf.int
table2[8,3]=paste0(round(me.alk.u,digits=2)," (",round(me.alk.u.ci[2],digits=2),"-",round(me.alk.u.ci[1],digits=2),")")

me.alk.c=mean(unlist(c.alk[300,]))
me.alk.c.ci=t.test(unlist(c.alk[300,]))$conf.int
table2[8,2]=paste0(prettyNum(round(me.alk.c,digits=0), big.mark = ",")," (",prettyNum(round(me.alk.c.ci[2],digits=0), big.mark = ","),"-",prettyNum(round(me.alk.c.ci[1],digits=0), big.mark = ","),")")

print("H 1L")
c.H=c.df[,s.df[1,]=="H 1L"]
u.H=u.df[,s.df[1,]=="H 1L"]

me.H.u=mean(unlist(u.H[300,]))
me.H.u.ci=t.test(unlist(u.H[300,]))$conf.int
table2[9,3]=paste0(round(me.H.u,digits=2)," (",round(me.H.u.ci[2],digits=2),"-",round(me.H.u.ci[1],digits=2),")")

me.H.c=mean(unlist(c.H[300,]))
me.H.c.ci=t.test(unlist(c.H[300,]))$conf.int
table2[9,2]=paste0(prettyNum(round(me.H.c,digits=0), big.mark = ",")," (",prettyNum(round(me.H.c.ci[2],digits=0), big.mark = ","),"-",prettyNum(round(me.H.c.ci[1],digits=0), big.mark = ","),")")

print("L 1L")
c.L=c.df[,s.df[1,]=="L 1L"]
u.L=u.df[,s.df[1,]=="L 1L"]

me.L.u=mean(unlist(u.L[300,]))
me.L.u.ci=t.test(unlist(u.L[300,]))$conf.int
table2[10,3]=paste0(round(me.L.u,digits=2)," (",round(me.L.u.ci[2],digits=2),"-",round(me.L.u.ci[1],digits=2),")")

me.L.c=mean(unlist(c.L[300,]))
me.L.u.ci=t.test(unlist(c.L[300,]))$conf.int
table2[10,2]=paste0(prettyNum(round(me.L.c,digits=0), big.mark = ",")," (",prettyNum(round(me.L.u.ci[2],digits=0), big.mark = ","),"-",prettyNum(round(me.L.u.ci[1],digits=0), big.mark = ","),")")

print("Stage IV")
c.IV=c.df[,s.df[1,]=="L 1L" | s.df[1,]=="H 1L"| s.df[1,]=="ALK 1L" |  s.df[1,]=="EGFR 1L"]
u.IV=u.df[,s.df[1,]=="L 1L" | s.df[1,]=="H 1L"| s.df[1,]=="ALK 1L" |  s.df[1,]=="EGFR 1L"]

me.IV.u=mean(unlist(u.IV[300,]))
me.IV.u.ci=t.test(unlist(u.IV[300,]))$conf.int
table2[6,3]=paste0(round(me.IV.u,digits=2)," (",round(me.IV.u.ci[2],digits=2),"-",round(me.IV.u.ci[1],digits=2),")")

me.IV.c=mean(unlist(c.IV[300,]))
me.IV.u.ci=t.test(unlist(c.IV[300,]))$conf.int
table2[6,2]=paste0(prettyNum(round(me.IV.c,digits=0), big.mark = ",")," (",prettyNum(round(me.IV.u.ci[2],digits=0), big.mark = ","),"-",prettyNum(round(me.IV.u.ci[1],digits=0), big.mark = ","),")")

#write.csv(table2, paste0(getwd(),"/Outputs/Table 2.csv"), row.names = FALSE)
write.csv(table2, paste0(getwd(),"/Outputs/Table 2 LY.csv"), row.names = FALSE)

############################# Deterministic Cost Effectiveness #############################

#Whole disease
IncCosts.ind=(mean(unlist(c.df.comp[300,]))-mean(unlist(c.df[300,])))
IncQALY.ind=(mean(unlist(u.df.comp[300,]))-mean(unlist(u.df[300,])))

IncCosts.ind
IncQALY.ind
IncCosts.ind/IncQALY.ind
(IncQALY.ind*50000)-IncCosts.ind


IncCosts.ind=(mean(unlist(c.df.comp2[300,]))-mean(unlist(c.df[300,])))
IncQALY.ind=(mean(unlist(u.df.comp2[300,]))-mean(unlist(u.df[300,])))

IncCosts.ind
IncQALY.ind
IncCosts.ind/IncQALY.ind
(IncQALY.ind*50000)-IncCosts.ind

#EGFR subgroup

c.egfr.det=c.df[,s.df[1,]=="EGFR 1L"]
u.egfr.det=u.df[,s.df[1,]=="EGFR 1L"]
#u.egfr
mean(unlist(u.egfr.det[300,]))
mean(unlist(c.egfr.det[300,]))


c.egfr.comp.det=c.df.comp[,s.df.comp[1,]=="EGFR 1L"]
u.egfr.comp.det=u.df.comp[,s.df.comp[1,]=="EGFR 1L"]

mean(unlist(u.egfr.comp.det[300,]))
mean(unlist(c.egfr.comp.det[300,]))


IC=mean(unlist(c.egfr.comp.det[300,]))-mean(unlist(c.egfr.det[300,]))
IQ=mean(unlist(u.egfr.comp.det[300,]))-mean(unlist(u.egfr.det[300,]))
IC
IQ
IC/IQ


c.egfr.comp2.det=c.df.comp2[,s.df.comp2[1,]=="EGFR 1L"]
u.egfr.comp2.det=u.df.comp2[,s.df.comp2[1,]=="EGFR 1L"]

mean(unlist(u.egfr.comp2.det[300,]))
mean(unlist(c.egfr.comp2.det[300,]))

IC=mean(unlist(c.egfr.comp2.det[300,]))-mean(unlist(c.egfr.det[300,]))
IQ=mean(unlist(u.egfr.comp2.det[300,]))-mean(unlist(u.egfr.det[300,]))
IC
IQ
IC/IQ

