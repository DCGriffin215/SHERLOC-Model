

incCosts=(c.df.avg.comp-c.df.avg)
incQALY=(u.df.avg.comp-u.df.avg)
ICER=incCosts/incQALY
NMB=(incQALY*50000)-incCosts

PSAresult=data.frame("incCosts"=c(unlist(incCosts)),"incQALY"=c(unlist(incQALY)),"ICER"=c(unlist(ICER)),"NMB"=c(unlist(NMB)))

incCosts=(c.df.avg.comp2-c.df.avg)
incQALY=(u.df.avg.comp2-u.df.avg)
ICER=incCosts/incQALY
NMB=(incQALY*50000)-incCosts

PSAresult2=data.frame("incCosts"=c(unlist(incCosts)),"incQALY"=c(unlist(incQALY)),"ICER"=c(unlist(ICER)),"NMB"=c(unlist(NMB)))


write.csv(c.df.avg,paste0(getwd(),"/Outputs/PSA flaura2 cost avg.csv"), row.names = FALSE)
write.csv(c.df.avg.comp,paste0(getwd(),"/Outputs/PSA flaura2 cost avg comp.csv"), row.names = FALSE)
write.csv(u.df.avg,paste0(getwd(),"/Outputs/PSA flaura2 util avg.csv"), row.names = FALSE)
write.csv(u.df.avg.comp,paste0(getwd(),"/Outputs/PSA flaura2 util avg comp.csv"), row.names = FALSE)
write.csv(c.df.tot,paste0(getwd(),"/Outputs/PSA flaura2 cost tot.csv"), row.names = FALSE)
write.csv(c.df.tot.comp,paste0(getwd(),"/Outputs/PSA flaura2 cost tot comp.csv"), row.names = FALSE)
write.csv(u.df.tot,paste0(getwd(),"/Outputs/PSA flaura2 util tot.csv"), row.names = FALSE)
write.csv(u.df.tot.comp,paste0(getwd(),"/Outputs/PSA flaura2 util tot comp.csv"), row.names = FALSE)
write.csv(PSAresult,paste0(getwd(),"/Outputs/PSA flaura2 inc and ICERs.csv"), row.names = FALSE)
write.csv(c.df.avg.comp2,paste0(getwd(),"/Outputs/PSA MARIPOSA cost avg comp.csv"), row.names = FALSE)
write.csv(u.df.avg.comp2,paste0(getwd(),"/Outputs/PSA MARIPOSA util avg comp.csv"), row.names = FALSE)
write.csv(c.df.tot.comp2,paste0(getwd(),"/Outputs/PSA MARIPOSA cost tot comp.csv"), row.names = FALSE)
write.csv(u.df.tot.comp2,paste0(getwd(),"/Outputs/PSA MARIPOSA util tot comp.csv"), row.names = FALSE)
write.csv(PSAresult2,paste0(getwd(),"/Outputs/PSA MARIPOSA inc and ICERs.csv"), row.names = FALSE)


WTP.thresh=data.frame()
WTP.alt=data.frame()
for (i in 0:6000){
  count=0
  for(j in 1:1000){
    temp=PSAresult[j,"ICER"]<(i*100)
    temp2=PSAresult[j,"incQALY"]>0
    if (temp==TRUE & temp2==TRUE){count=count+1}
  }
  WTP.thresh[i+1,"prob"]=count/1000
  WTP.thresh[i+1,"WTP"]=(i*100)
  WTP.thresh[i+1,"Treatment"]="FLAURA2"
  
  WTP.alt[i+1,"prob"]=1-(count/1000)
  WTP.alt[i+1,"WTP"]=(i*100)
  WTP.alt[i+1,"Treatment"]="Standard Of Care"
}
WTP.thresh=rbind(WTP.thresh,WTP.alt)
WTP.thresh


WTP.thresh2=data.frame()
WTP.alt2=data.frame()
for (i in 0:6000){
  count=0
  for(j in 1:1000){
    temp=PSAresult2[j,"ICER"]<(i*100)
    temp2=PSAresult2[j,"incQALY"]>0
    if (temp==TRUE & temp2==TRUE){count=count+1}
  }
  WTP.thresh2[i+1,"prob"]=count/1000
  WTP.thresh2[i+1,"WTP"]=(i*100)
  WTP.thresh2[i+1,"Treatment"]="MARIPOSA"
  
  WTP.alt2[i+1,"prob"]=1-(count/1000)
  WTP.alt2[i+1,"WTP"]=(i*100)
  WTP.alt2[i+1,"Treatment"]="Standard Of Care"
}
WTP.thresh2=rbind(WTP.thresh2,WTP.alt2)
WTP.thresh2

write.csv(WTP.thresh,paste0(getwd(),"/Outputs/PSA flaura2 WTP.csv"), row.names = FALSE)
write.csv(WTP.thresh2,paste0(getwd(),"/Outputs/PSA MARIPOSA WTP.csv"), row.names = FALSE)


########################### EGFR Results ###########################


incCosts=(c.egfr.comp-c.egfr)
incQALY=(u.egfr.comp-u.egfr)
ICER=incCosts/incQALY
NMB=(incQALY*50000)-incCosts

PSAresult=data.frame("incCosts"=c(unlist(incCosts)),"incQALY"=c(unlist(incQALY)),"ICER"=c(unlist(ICER)),"NMB"=c(unlist(NMB)))

incCosts=(c.egfr.comp2-c.egfr)
incQALY=(u.egfr.comp2-u.egfr)
ICER=incCosts/incQALY
NMB=(incQALY*50000)-incCosts

PSAresult2=data.frame("incCosts"=c(unlist(incCosts)),"incQALY"=c(unlist(incQALY)),"ICER"=c(unlist(ICER)),"NMB"=c(unlist(NMB)))


write.csv(c.egfr,paste0(getwd(),"/Outputs/PSA  flaura2 c_egfr.csv"), row.names = FALSE)
write.csv(c.egfr.comp,paste0(getwd(),"/Outputs/PSA flaura2 c_egfr_comp.csv"), row.names = FALSE)
write.csv(c.egfr.comp2,paste0(getwd(),"/Outputs/PSA MARIPOSA c_egfr_comp2.csv"), row.names = FALSE)
write.csv(u.egfr,paste0(getwd(),"/Outputs/PSA  flaura2 u_egfr.csv"), row.names = FALSE)
write.csv(u.egfr.comp,paste0(getwd(),"/Outputs/PSA flaura2 u_egfr_comp.csv"), row.names = FALSE)
write.csv(u.egfr.comp2,paste0(getwd(),"/Outputs/PSA MARIPOSA u_egfr_comp2.csv"), row.names = FALSE)
write.csv(PSAresult,paste0(getwd(),"/Outputs/PSA flaura2 EGFR inc and ICERs.csv"), row.names = FALSE)
write.csv(PSAresult2,paste0(getwd(),"/Outputs/PSA MARIPOSA EGFR inc and ICERs.csv"), row.names = FALSE)


WTP.thresh=data.frame()
WTP.alt=data.frame()
for (i in 0:6000){
  count=0
  for(j in 1:1000){
    temp=PSAresult[j,"ICER"]<(i*100)
    temp2=PSAresult[j,"incQALY"]>0
    if (temp==TRUE & temp2==TRUE){count=count+1}
  }
  WTP.thresh[i+1,"prob"]=count/1000
  WTP.thresh[i+1,"WTP"]=(i*100)
  WTP.thresh[i+1,"Treatment"]="FLAURA2"
  
  WTP.alt[i+1,"prob"]=1-(count/1000)
  WTP.alt[i+1,"WTP"]=(i*100)
  WTP.alt[i+1,"Treatment"]="Standard Of Care"
}
WTP.thresh=rbind(WTP.thresh,WTP.alt)
WTP.thresh


WTP.thresh2=data.frame()
WTP.alt2=data.frame()
for (i in 0:6000){
  count=0
  for(j in 1:1000){
    temp=PSAresult2[j,"ICER"]<(i*100)
    temp2=PSAresult2[j,"incQALY"]>0
    if (temp==TRUE & temp2==TRUE){count=count+1}
  }
  WTP.thresh2[i+1,"prob"]=count/1000
  WTP.thresh2[i+1,"WTP"]=(i*100)
  WTP.thresh2[i+1,"Treatment"]="MARIPOSA"
  
  WTP.alt2[i+1,"prob"]=1-(count/1000)
  WTP.alt2[i+1,"WTP"]=(i*100)
  WTP.alt2[i+1,"Treatment"]="Standard Of Care"
}
WTP.thresh2=rbind(WTP.thresh2,WTP.alt2)
WTP.thresh2

write.csv(WTP.thresh,paste0(getwd(),"/Outputs/PSA flaura2 EGFR WTP.csv"), row.names = FALSE)
write.csv(WTP.thresh2,paste0(getwd(),"/Outputs/PSA MARIPOSA EGFR WTP.csv"), row.names = FALSE)


########################### Survival Results ###########################

write.csv(surv1y,paste0(getwd(),"/Outputs/PSA Survival 1 Year.csv"), row.names = FALSE)
write.csv(surv3y,paste0(getwd(),"/Outputs/PSA Survival 3 Year.csv"), row.names = FALSE)
write.csv(surv5y,paste0(getwd(),"/Outputs/PSA Survival 5 Year.csv"), row.names = FALSE)
write.csv(surv10y,paste0(getwd(),"/Outputs/PSA Survival 10 Year.csv"), row.names = FALSE)

write.csv(surv1y.comp,paste0(getwd(),"/Outputs/PSA Survival 1 Year Comparitor1.csv"), row.names = FALSE)
write.csv(surv3y.comp,paste0(getwd(),"/Outputs/PSA Survival 3 Year Comparitor1.csv"), row.names = FALSE)
write.csv(surv5y.comp,paste0(getwd(),"/Outputs/PSA Survival 5 Year Comparitor1.csv"), row.names = FALSE)
write.csv(surv10y.comp,paste0(getwd(),"/Outputs/PSA Survival 10 Year Comparitor1.csv"), row.names = FALSE)

write.csv(surv1y.comp2,paste0(getwd(),"/Outputs/PSA Survival 1 Year Comparitor2.csv"), row.names = FALSE)
write.csv(surv3y.comp2,paste0(getwd(),"/Outputs/PSA Survival 3 Year Comparitor2.csv"), row.names = FALSE)
write.csv(surv5y.comp2,paste0(getwd(),"/Outputs/PSA Survival 5 Year Comparitor2.csv"), row.names = FALSE)
write.csv(surv10y.comp2,paste0(getwd(),"/Outputs/PSA Survival 10 Year Comparitor2.csv"), row.names = FALSE)