#rm(list = ls())
library(parallel);
cl_sock <- makePSOCKcluster(detectCores() *0.75);


NSCLCmicro=function(start_i, end_i,comp,mari){
  
  #model parameters
  n.i=end_i-start_i
  n.t=300
  #Costs
  
  Osimertinib=9020
  Alectinib=9918
  Lorlatinib=8982
  Crizotinib=7280.00
  Pembrolizumab=11733
  Durvalumab=11733
  Atezolizumab=9518.91
  Pemetrexed=600
  Carboplatin=987
  Cisplatin=540
  Docetaxel=2837
  Pembrolizumab.Carboplatin.Pemetrexed=Pembrolizumab+Carboplatin+Pemetrexed
  Pembrolizumab.Pemetrexed=Pembrolizumab+Pemetrexed
  Cisplatin.Pemetrexed=Cisplatin+Pemetrexed
  Carbo.Pemetrexed=Carboplatin+Pemetrexed
  Cisplatin.Docetaxel=Cisplatin+Docetaxel
  Nivolumab=9387
  Lazert.Amivan=12135
  
  
  Consultation=207.17
  Office.Visit=49.72
  CT.CAP.with.contrast=1275.00
  MRI.brain=730.00
  PET.Body=999.00
  Chemistry=71.75
  Electrolytes=6.53
  Liver=53.57
  Renal=11.65
  Hematology=18.96
  NGS=1306.75
  Bronchoscopy=2856.77
  Endoscopic.lung.resection=17079.85
  Radiation.1=811.2163743*10
  Radiation.2=811.2163743*10	
  Radiation.3=811.2163743*30
  BSC.cost=5772.608396
  hospital=9959.969718
  
  
  
  
  #Costing
  
  Surveilance=Office.Visit+(CT.CAP.with.contrast*2)
  Surveilance.4=Office.Visit+(CT.CAP.with.contrast*4)
  stage1.cost.misc=(Consultation*1)+(Office.Visit*2)+(PET.Body*1)+(MRI.brain*1)+(Chemistry*2)+(Hematology*2)+(Bronchoscopy*1)
  stage2.cost.misc=(Consultation*1)+(Office.Visit*5)+(CT.CAP.with.contrast*1)+(MRI.brain*1)+(PET.Body*1)+(Chemistry*6)+(Hematology*6)+(Bronchoscopy*1)
  stage3.cost.misc=(Consultation*1)+(Office.Visit*17)+(CT.CAP.with.contrast*4)+(MRI.brain*1)+(PET.Body*1)+(Chemistry*20)+(Hematology*20)+(NGS*1)+(Bronchoscopy*1)
  stage4.cost.misc=(Consultation*1)+(Office.Visit*17)+(CT.CAP.with.contrast*4)+(MRI.brain*1)+(PET.Body*1)+(Chemistry*20)+(Hematology*20)+(NGS*1)+(Bronchoscopy*1)
  stage42L.cost.misc=(Office.Visit*17)+(CT.CAP.with.contrast*4)+(Hematology*17)
  stage42L.cost.misc.act=(Office.Visit*17)+(CT.CAP.with.contrast*4)+(Hematology*17)+(MRI.brain*2)
  
  
  #Utilities
  
  u.s1.s=0.82/12
  u.s1.r=0.76/12
  u.s1.dfs=0.82/12
  
  u.s2.cs=0.77/12
  u.s2.r=0.76/12
  u.s2.dfs=0.80/12
  
  u.s3.cr=0.76/12
  u.s3.cs=0.78/12
  u.s3.dfs=0.79/12
  
  u.H.1L=0.71/12
  u.L.1L=0.71/12
  u.A.1L=0.80/12
  u.H.2L=0.67/12
  u.L.2L=0.67/12
  u.A.2L=0.64/12
  u.BSC=0.473/12
  
  if(mari==FALSE){u.A.1L.comp=(0.8-0.039)/12}
  if(mari==TRUE){u.A.1L.comp=(0.8-0.019)/12}
  
  
  
  #u.s1.s=1/12
  #u.s1.r=1/12
  #u.s1.dfs=1/12
  
  #u.s2.cs=1/12
  #u.s2.r=1/12
  #u.s2.dfs=1/12
  
  #u.s3.cr=1/12
  #u.s3.cs=1/12
  #u.s3.dfs=1/12
  
  #u.H.1L=1/12
  #u.L.1L=1/12
  #u.A.1L=1/12
  #u.H.2L=1/12
  #u.L.2L=1/12
  #u.A.2L=1/12
  #u.BSC=1/12
  
  #u.A.1L.comp=1/12
  
  #transition probabilities
  s4=(1193/2410)
  s1=(633/2410)
  s2=(170/2410)
  s3=(414/2410)
  
  
  #internal cohort staging & 
  #s4=(534/868)
  #s1=(142/868)
  #s2=(137/868)
  #s3=(55/868)
  
  #act=(207/534)
  #low=((82/534)+(87/534))
  #high=(158/534)
  #alk=act*(72/203)
  #egfr=act*(131/203)
  
  #real world mutation and 
  alk=0.27*0.14
  egfr=0.27*0.77
  low=(1-(egfr+alk))*0.6
  high=(1-(egfr+alk))*0.4
  
  #alk=0.27*0.14
  #egfr=0.5
  #low=(1-(egfr+alk))*0.6
  #high=(1-(egfr+alk))*0.4
  
  parms=read.csv(paste0(getwd(),"/NSCLC microsim data/parms_comb.csv"))
  
  df.list=Probs.determ(low,high,alk,egfr,parms)
  
  sp.1=df.list$sp.1
  sp.2=df.list$sp.2
  sp.3=df.list$sp.3
  sp.alk.1L=df.list$sp.alk.1L
  
  if (comp==FALSE){sp.egfr.1L=df.list$sp.egfr.1L}
  if (comp==TRUE & mari==FALSE){sp.egfr.1L=df.list$sp.FL2_HR.1L}
  if (comp==TRUE & mari==TRUE){sp.egfr.1L=df.list$sp.MARIP_HR.1L}
  
  sp.low.1L=df.list$sp.low.1L
  sp.high.1L=df.list$sp.high.1L
  sp.low.2L=df.list$sp.low.2L
  sp.high.2L=df.list$sp.high.2L
  sp.egfr.2L=df.list$sp.egfr.2L
  sp.alk.2L=df.list$sp.alk.2L
  
  #checks
  s4+s1+s2+s3==1
  alk+egfr+low+high==1
  
  #dataframes
  
  
  
  t.df=data.frame()
  
  n.rows=301
  s.mat=matrix(NA, nrow = n.rows, ncol = n.i+1)
  c.mat=matrix(NA, nrow = n.rows-1, ncol = n.i+1)
  u.mat=matrix(NA, nrow = n.rows-1, ncol = n.i+1)
  
  #discount function
  d.rt=0.015
  discount=function(unit,t,d.rt){unit/((1+d.rt)^((t-1)/12))}
  
  c.tot=0
  u.tot=0
  
  #begin simulation
  
  for (i in start_i:end_i){
    set.seed(1 + i) #set seed for each individual
    i=(i+1)-start_i
    #stage counters
    s1.ct=0
    s2.ct=0
    s3.ct=0
    alk.1L.ct=0
    egfr.1L.ct=0
    high.1L.ct=0
    low.1L.ct=0
    alk.2L.ct=0
    egfr.2L.ct=0
    high.2L.ct=0
    low.2L.ct=0
    
    #costs and utility counters
    
    n.s=sample(c("S1","S2","S3","S4"),prob=c(s1,s2,s3,s4),size=1)
    n.s
    if (n.s=="S4"){
      n.s=sample(c("ALK 1L","H 1L","L 1L","EGFR 1L"),prob=c(alk,high,low,egfr),size=1)
    }
    s.mat[1,i]=n.s
    
    c.ind=0
    u.ind=0
    
    #begin sim
    
    for (t in 1:n.t){
      c.s=n.s
      c.temp=0
      u.temp=0
      if (c.s=="S1"){
        if (s1.ct==0){s1.treat=sample(c("Surg","Rad"),prob=c(0.52349,0.47651),size=1)}
        s1.ct=s1.ct+1
        if (s1.treat=="Surg"){
          if (s1.ct==1){c.temp=c.temp+stage1.cost.misc+Endoscopic.lung.resection}
          if (s1.ct>1 & s1.ct<=61){c.temp=c.temp+(Surveilance/12)}
          
          if (s1.ct<=4){u.temp=u.temp+u.s1.s}
          if (s1.ct>4){u.temp=u.temp+u.s1.dfs}
          
        }
        if (s1.treat=="Rad"){
          if (s1.ct==1){c.temp=c.temp+((stage1.cost.misc+Radiation.1)/2)}
          if (s1.ct==2){c.temp=c.temp+((stage1.cost.misc+Radiation.1)/2)}
          if (s1.ct>2 & s1.ct<=62){c.temp=c.temp+(Surveilance/12)}
          
          if (s1.ct<=3){u.temp=u.temp+u.s1.r}
          if (s1.ct>3){u.temp=u.temp+u.s1.dfs}
          
        }
        
        
        c.ind=c.ind+discount(c.temp,t,d.rt)
        u.ind=u.ind+discount(u.temp,t,d.rt)
        
        n.s=sample(c("S1","Dead","ALK 1L","L 1L","H 1L","EGFR 1L"),prob=unlist(sp.1[s1.ct,]),size=1)
      }
      
      # Stage II
      
      if (c.s=="S2"){
        if (s2.ct==0){s2.treat=sample(c("Neoadj","Rad","Surg+adj"),prob=c((1-0.3252)*0.3,0.3252,(1-0.3252)*0.7),size=1)}
        s2.ct=s2.ct+1
        
        if (s2.treat=="Neoadj"){
          if (s2.ct==1){c.temp=c.temp+(stage2.cost.misc/3)+(Nivolumab+Cisplatin.Pemetrexed)*1.15}
          if (s2.ct==2){c.temp=c.temp+(stage2.cost.misc/3)+(Nivolumab+Cisplatin.Pemetrexed)*1.15}
          if (s2.ct==3){c.temp=c.temp+(stage2.cost.misc/3)+Endoscopic.lung.resection}
          if (s2.ct>3 & s2.ct<=63){c.temp=c.temp+(Surveilance/12)}
          
          if (s2.ct<=6){u.temp=u.temp+u.s2.cs}
          if (s2.ct>6){u.temp=u.temp+u.s2.dfs}
        }
        if (s2.treat=="Rad"){
          if (s2.ct==1){c.temp=c.temp+((stage2.cost.misc+Radiation.2)/2)}
          if (s2.ct==2){c.temp=c.temp+((stage2.cost.misc+Radiation.2)/2)}
          if (s2.ct>2 & s2.ct<=62){c.temp=c.temp+(Surveilance/12)}
          
          if (s2.ct<=3){u.temp=u.temp+u.s2.r}
          if (s2.ct>3){u.temp=u.temp+u.s2.dfs}
        }
        if (s2.treat=="Surg+adj"){
          if (s2.ct==1){c.temp=c.temp+(stage2.cost.misc/12)+Cisplatin.Pemetrexed+Endoscopic.lung.resection}
          if (s2.ct==2){c.temp=c.temp+(stage2.cost.misc/12)+Cisplatin.Pemetrexed}
          if (s2.ct==3){c.temp=c.temp+(stage2.cost.misc/12)+Cisplatin.Pemetrexed}
          if (s2.ct>3 & s2.ct<=15){c.temp=c.temp+(stage2.cost.misc/12)+(Atezolizumab*0.5)+(Surveilance*0.5)+((hospital/12)*0.18*0.5)}
          if (s2.ct>15 & s2.ct<=75){c.temp=c.temp+(Surveilance/12)}
          
          if (s2.ct<=6){u.temp=u.temp+u.s2.cs}
          if (s2.ct>6 & s2.ct<=18){u.temp=u.temp+(u.s2.dfs*0.5)+(u.s2.cs*0.5)}
          if (s2.ct>18){u.temp=u.temp+u.s2.dfs}
        }
        
        c.ind=c.ind+discount(c.temp,t,d.rt)
        u.ind=u.ind+discount(u.temp,t,d.rt)
        
        n.s=sample(c("S2","Dead","ALK 1L","L 1L","H 1L","EGFR 1L"),prob=unlist(sp.2[s2.ct,]),size=1)
      }
      
      # Stage III
      
      if (c.s=="S3"){
        if (s3.ct==0){s3.treat=sample(c("Unresect","Neoadj","Surg+adj"),prob=c(0.7,(0.3*0.3),(0.3*0.7)),size=1)}
        s3.ct=s3.ct+1
        
        if (s3.treat=="Unresect"){
          if (s3.ct==1){c.temp=c.temp+(stage3.cost.misc/12)+Cisplatin.Pemetrexed+(Radiation.3/2)}
          if (s3.ct==2){c.temp=c.temp+(stage3.cost.misc/12)+Cisplatin.Pemetrexed+(Radiation.3/2)}
          if (s3.ct==3){c.temp=c.temp+(stage3.cost.misc/12)+Cisplatin.Pemetrexed}
          if (s3.ct>3 & s3.ct<=15){c.temp=c.temp+(stage3.cost.misc/12)+Durvalumab}
          if (s3.ct>15 & s3.ct<=75){c.temp=c.temp+(Surveilance/12)}
          
          if (s3.ct<=18){u.temp=u.temp+u.s3.cr}
          if (s3.ct>18){u.temp=u.temp+u.s3.dfs}
        }
        if (s3.treat=="Neoadj"){
          if (s3.ct==1){c.temp=c.temp+(stage3.cost.misc/3)+(Nivolumab+Cisplatin.Pemetrexed)*1.15}
          if (s3.ct==2){c.temp=c.temp+(stage3.cost.misc/3)+(Nivolumab+Cisplatin.Pemetrexed)*1.15}
          if (s3.ct==3){c.temp=c.temp+(stage3.cost.misc/3)+Endoscopic.lung.resection}
          if (s3.ct>3 & s3.ct<=63){c.temp=c.temp+(Surveilance/12)}
          
          if (s3.ct<=6){u.temp=u.temp+u.s3.cs}
          if (s3.ct>6){u.temp=u.temp+u.s3.dfs}
        }
        if (s3.treat=="Surg+adj"){
          if (s3.ct==1){c.temp=c.temp+(stage3.cost.misc/12)+Cisplatin.Pemetrexed+Endoscopic.lung.resection}
          if (s3.ct==2){c.temp=c.temp+(stage3.cost.misc/12)+Cisplatin.Pemetrexed}
          if (s3.ct==3){c.temp=c.temp+(stage3.cost.misc/12)+Cisplatin.Pemetrexed}
          if (s3.ct>3 & s3.ct<=15){c.temp=c.temp+(stage3.cost.misc/12)+(Atezolizumab*0.5)+(Surveilance*0.5)+((hospital/12)*0.18*0.5)}
          if (s3.ct>15 & s3.ct<=75){c.temp=c.temp+(Surveilance/12)}
          
          if (s3.ct<=6){u.temp=u.temp+u.s3.cs}
          if (s3.ct>6 & s3.ct<=18){u.temp=u.temp+(u.s3.dfs*0.5)+(u.s3.cs*0.5)}
          if (s3.ct>18){u.temp=u.temp+u.s3.dfs}
        }
        
        c.ind=c.ind+discount(c.temp,t,d.rt)
        u.ind=u.ind+discount(u.temp,t,d.rt)
        
        n.s=sample(c("S3","Dead","ALK 1L","L 1L","H 1L","EGFR 1L"),prob=unlist(sp.3[s3.ct,]),size=1)
      }
      
      if (c.s=="ALK 1L"){
        alk.1L.ct=alk.1L.ct+1
        
        c.temp=c.temp+(stage4.cost.misc/12)+Alectinib
        
        u.temp=u.temp+u.A.1L
        
        
        c.ind=c.ind+discount(c.temp,t,d.rt)
        u.ind=u.ind+discount(u.temp,t,d.rt)
        
        n.s=sample(c("ALK 1L","ALK 2L","BSC"),prob=unlist(sp.alk.1L[alk.1L.ct,]),size=1)
      }
      
      if (c.s=="EGFR 1L" & comp==FALSE){
        egfr.1L.ct=egfr.1L.ct+1
        
        c.temp=c.temp+(stage4.cost.misc/12)+Osimertinib
        
        u.temp=u.temp+u.A.1L
        
        
        c.ind=c.ind+discount(c.temp,t,d.rt)
        u.ind=u.ind+discount(u.temp,t,d.rt)
        
        n.s=sample(c("EGFR 1L","EGFR 2L","BSC"),prob=unlist(sp.egfr.1L[egfr.1L.ct,]),size=1)
      }
      
      if (c.s=="EGFR 1L" & comp==TRUE & mari==FALSE){
        egfr.1L.ct=egfr.1L.ct+1
        
        if (egfr.1L.ct==1){c.temp=c.temp+(stage4.cost.misc/12)+Cisplatin.Pemetrexed+Osimertinib+(hospital*0.14)}
        if (egfr.1L.ct==2){c.temp=c.temp+(stage4.cost.misc/12)+Cisplatin.Pemetrexed+Osimertinib}
        if (egfr.1L.ct==3){c.temp=c.temp+(stage4.cost.misc/12)+Cisplatin.Pemetrexed+Osimertinib}
        if (egfr.1L.ct>3){c.temp=c.temp+(stage4.cost.misc/12)+Pemetrexed+Osimertinib}
        
        
        u.temp=u.temp+u.A.1L.comp
        
        
        c.ind=c.ind+discount(c.temp,t,d.rt)
        u.ind=u.ind+discount(u.temp,t,d.rt)
        
        n.s=sample(c("EGFR 1L","EGFR 2L","BSC"),prob=unlist(sp.egfr.1L[egfr.1L.ct,]),size=1)
      }
      
      if (c.s=="EGFR 1L" & comp==TRUE & mari==TRUE){
        egfr.1L.ct=egfr.1L.ct+1
        
        c.temp=c.temp+(stage4.cost.misc/12)+Lazert.Amivan
        
        u.temp=u.temp+u.A.1L.comp
        
        
        c.ind=c.ind+discount(c.temp,t,d.rt)
        u.ind=u.ind+discount(u.temp,t,d.rt)
        
        n.s=sample(c("EGFR 1L","EGFR 2L","BSC"),prob=unlist(sp.egfr.1L[egfr.1L.ct,]),size=1)
      }
      
      if (c.s=="H 1L"){
        high.1L.ct= high.1L.ct+1
        
        if (high.1L.ct<=24){c.temp=c.temp+(stage4.cost.misc/12)+Pembrolizumab+((hospital/24)*0.15)}
        if (high.1L.ct>24){c.temp=c.temp+(Surveilance.4/12)}
        
        u.temp=u.temp+u.H.1L
        
        c.ind=c.ind+discount(c.temp,t,d.rt)
        u.ind=u.ind+discount(u.temp,t,d.rt)
        
        n.s=sample(c("H 1L","H 2L","BSC"),prob=unlist(sp.high.1L[high.1L.ct,]),size=1)
      }
      
      if (c.s=="L 1L"){
        low.1L.ct=low.1L.ct+1
        
        if (low.1L.ct<=3){c.temp=c.temp+(stage4.cost.misc/12)+Pembrolizumab.Carboplatin.Pemetrexed+((hospital/24)*0.15)}
        if (low.1L.ct>3 &low.1L.ct<=24){c.temp=c.temp+(stage4.cost.misc/12)+Pembrolizumab.Pemetrexed+((hospital/24)*0.15)}
        if (low.1L.ct>24){c.temp=c.temp+(Surveilance.4/12)}
        
        u.temp=u.temp+u.L.1L
        
        c.ind=c.ind+discount(c.temp,t,d.rt)
        u.ind=u.ind+discount(u.temp,t,d.rt)
        
        n.s=sample(c("L 1L","L 2L","BSC"),prob=unlist(sp.low.1L[low.1L.ct,]),size=1)
      }
      
      if (c.s=="ALK 2L"){
        alk.2L.ct=alk.2L.ct+1
        
        c.temp=c.temp+(stage42L.cost.misc.act/12)+Lorlatinib
        
        u.temp=u.temp+u.A.2L
        
        c.ind=c.ind+discount(c.temp,t,d.rt)
        u.ind=u.ind+discount(u.temp,t,d.rt)
        
        n.s=sample(c("ALK 2L","BSC"),prob=unlist(sp.alk.2L[alk.2L.ct,]),size=1)
      }
      
      if (c.s=="EGFR 2L" & comp==FALSE){
        egfr.2L.ct=egfr.2L.ct+1
        
        c.temp=c.temp+(stage42L.cost.misc.act/12)+Carboplatin+Pemetrexed
        
        u.temp=u.temp+u.A.2L
        
        c.ind=c.ind+discount(c.temp,t,d.rt)
        u.ind=u.ind+discount(u.temp,t,d.rt)
        
        n.s=sample(c("EGFR 2L","BSC"),prob=unlist(sp.egfr.2L[egfr.2L.ct,]),size=1)
      }
      
      if (c.s=="EGFR 2L" & mari==TRUE){
        egfr.2L.ct=egfr.2L.ct+1
        
        c.temp=c.temp+(stage42L.cost.misc.act/12)+Carboplatin+Pemetrexed
        
        u.temp=u.temp+u.A.2L
        
        c.ind=c.ind+discount(c.temp,t,d.rt)
        u.ind=u.ind+discount(u.temp,t,d.rt)
        
        n.s=sample(c("EGFR 2L","BSC"),prob=unlist(sp.egfr.2L[egfr.2L.ct,]),size=1)
      }
      
      if (c.s=="EGFR 2L" & comp==TRUE & mari==FALSE){
        egfr.2L.ct=egfr.2L.ct+1
        
        c.temp=c.temp+(stage42L.cost.misc.act/12)+(Carboplatin+Pemetrexed*0.5)+(Docetaxel*0.5)
        
        u.temp=u.temp+u.A.2L
        
        c.ind=c.ind+discount(c.temp,t,d.rt)
        u.ind=u.ind+discount(u.temp,t,d.rt)
        
        n.s=sample(c("EGFR 2L","BSC"),prob=unlist(sp.egfr.2L[egfr.2L.ct,]),size=1)
      }
      
      if (c.s=="H 2L"){
        high.2L.ct=high.2L.ct+1
        
        c.temp=c.temp+(stage42L.cost.misc/12)+Carboplatin+Pemetrexed
        
        u.temp=u.temp+u.H.2L
        
        c.ind=c.ind+discount(c.temp,t,d.rt)
        u.ind=u.ind+discount(u.temp,t,d.rt)
        
        n.s=sample(c("H 2L","BSC"),prob=unlist(sp.high.2L[high.2L.ct,]),size=1)
      }
      
      if (c.s=="L 2L"){
        low.2L.ct=low.2L.ct+1
        
        c.temp=c.temp+(stage42L.cost.misc/12)+Docetaxel
        
        u.temp=u.temp+u.L.2L
        
        c.ind=c.ind+discount(c.temp,t,d.rt)
        u.ind=u.ind+discount(u.temp,t,d.rt)
        
        n.s=sample(c("L 2L","BSC"),prob=unlist(sp.low.2L[low.2L.ct,]),size=1)
      }
      
      if (c.s=="BSC"){
        
        c.temp=c.temp+BSC.cost
        
        u.temp=u.temp+u.BSC
        
        c.ind=c.ind+discount(c.temp,t,d.rt)
        u.ind=u.ind+discount(u.temp,t,d.rt)
        
        n.s="Dead"
      }
      s.mat[t+1,i]=n.s
      c.mat[t,i]=c.ind
      u.mat[t,i]=u.ind
    }
    c.tot=c.tot+c.ind
  }
  
  return(list(s.mat=s.mat,c.mat=c.mat,u.mat=u.mat))
  
}


source(paste0(getwd(),"/Probs_determ.R"))


#export functions
clusterExport(cl_sock, varlist = c("NSCLCmicro","Probs.determ"))

#divide simulation into chunks
n=10000
n.chunks=15
chunk_size=n / n.chunks
chunks=split(1:n, ceiling(seq_along(1:n) / chunk_size))

#run simulations
system.time({
  soc.output=parLapply(cl_sock, chunks, function(x) {
    
    library(readxl)
    
    NSCLCmicro(min(x), max(x),comp=FALSE,mari=FALSE)
  })
  
  comp.output=parLapply(cl_sock, chunks, function(x) {
    
    library(readxl)
    
    NSCLCmicro(min(x), max(x),comp=TRUE,mari=FALSE)
  })

  comp.output2=parLapply(cl_sock, chunks, function(x) {
  
    library(readxl)
  
    NSCLCmicro(min(x), max(x),comp=TRUE,mari=TRUE)
  })
})
#stop cluster
stopCluster(cl = cl_sock);

#extract matrixes
soc.list1=lapply(soc.output, function(res) res$s.mat)
soc.list2=lapply(soc.output, function(res) res$c.mat)
soc.list3=lapply(soc.output, function(res) res$u.mat)

# Convert to data frames
s.df=as.data.frame(do.call(cbind, soc.list1))
c.df=as.data.frame(do.call(cbind, soc.list2))
u.df=as.data.frame(do.call(cbind, soc.list3))
print(1-(table(unlist(s.df[13,]))[4]/10000))
print(1-(table(unlist(s.df[37,]))[4]/10000))
print(1-(table(unlist(s.df[61,]))[4]/10000))
print(1-(table(unlist(s.df[121,]))[4]/10000))
mean(unlist(c.df[300,]))
median(unlist(c.df[300,]))
t.test(unlist(c.df[300,]))$conf.int
sum(unlist(u.df[300,]))
mean(unlist(u.df[300,]))
t.test(unlist(u.df[300,]))$conf.int
median(unlist(u.df[300,]))

#extract matrixes
comp.list1=lapply(comp.output, function(res) res$s.mat)
comp.list2=lapply(comp.output, function(res) res$c.mat)
comp.list3=lapply(comp.output, function(res) res$u.mat)

# Convert to data frames
s.df.comp=as.data.frame(do.call(cbind, comp.list1))
c.df.comp=as.data.frame(do.call(cbind, comp.list2))
u.df.comp=as.data.frame(do.call(cbind, comp.list3))

#extract matrixes
comp2.list1=lapply(comp.output2, function(res) res$s.mat)
comp2.list2=lapply(comp.output2, function(res) res$c.mat)
comp2.list3=lapply(comp.output2, function(res) res$u.mat)

# Convert to data frames
s.df.comp2=as.data.frame(do.call(cbind, comp2.list1))
c.df.comp2=as.data.frame(do.call(cbind, comp2.list2))
u.df.comp2=as.data.frame(do.call(cbind, comp2.list3))

print(table(unlist(s.df.comp[1,])))  #check initial distribution
print(table(unlist(s.df.comp[61,]))) #check 5 year survival
print(table(unlist(s.df.comp[121,]))) #check 10 year survival

#summary results
mean(unlist(c.df.comp[300,]))
t.test(unlist(c.df.comp[300,]))$conf.int
median(unlist(c.df.comp[300,]))
sum(unlist(u.df.comp[300,]))
mean(unlist(u.df.comp[300,]))
t.test(unlist(u.df.comp[300,]))$conf.int
median(unlist(u.df.comp[300,]))

IncCosts.ind=(mean(unlist(c.df.comp[300,]))-mean(unlist(c.df[300,])))
IncQALY.ind=(mean(unlist(u.df.comp[300,]))-mean(unlist(u.df[300,])))
IncCosts.tot=(sum(unlist(c.df.comp[300,]))-sum(unlist(c.df[300,])))
IncQALY.tot=(sum(unlist(u.df.comp[300,]))-sum(unlist(u.df[300,])))

#cost effectiveness results
IncCosts.ind
IncQALY.ind
IncCosts.ind/IncQALY.ind
(IncQALY.ind*50000)-IncCosts.ind

IncCosts.tot
IncQALY.tot
IncCosts.tot/IncQALY.tot