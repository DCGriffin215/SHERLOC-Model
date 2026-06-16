rm(list = ls())
library(readxl)
library(parallel);
cl_sock=makePSOCKcluster(detectCores() *0.75);

#setup




NSCLCmicro=function(start_i, end_i,comp, mari){
  j.ind=1
  set.seed(1)
  
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
  
  
  Osimertinib.L=qgamma(runif(1000),shape=((Osimertinib^2)/((((Osimertinib*1.2)-Osimertinib)/1.96)^2)),scale=(((((Osimertinib*1.2)-Osimertinib)/1.96)^2)/Osimertinib))
  Alectinib.L=qgamma(runif(1000),shape=((Alectinib^2)/((((Alectinib*1.2)-Alectinib)/1.96)^2)),scale=(((((Alectinib*1.2)-Alectinib)/1.96)^2)/Alectinib))
  Lorlatinib.L=qgamma(runif(1000),shape=((Lorlatinib^2)/((((Lorlatinib*1.2)-Lorlatinib)/1.96)^2)),scale=(((((Lorlatinib*1.2)-Lorlatinib)/1.96)^2)/Lorlatinib))
  Crizotinib.L=qgamma(runif(1000),shape=((Crizotinib^2)/((((Crizotinib*1.2)-Crizotinib)/1.96)^2)),scale=(((((Crizotinib*1.2)-Crizotinib)/1.96)^2)/Crizotinib))
  Pembrolizumab.L=qgamma(runif(1000),shape=((Pembrolizumab^2)/((((Pembrolizumab*1.2)-Pembrolizumab)/1.96)^2)),scale=(((((Pembrolizumab*1.2)-Pembrolizumab)/1.96)^2)/Pembrolizumab))
  Durvalumab.L=qgamma(runif(1000),shape=((Durvalumab^2)/((((Durvalumab*1.2)-Durvalumab)/1.96)^2)),scale=(((((Durvalumab*1.2)-Durvalumab)/1.96)^2)/Durvalumab))
  Atezolizumab.L=qgamma(runif(1000),shape=((Atezolizumab^2)/((((Atezolizumab*1.2)-Atezolizumab)/1.96)^2)),scale=(((((Atezolizumab*1.2)-Atezolizumab)/1.96)^2)/Atezolizumab))
  Pemetrexed.L=qgamma(runif(1000),shape=((Pemetrexed^2)/((((Pemetrexed*1.2)-Pemetrexed)/1.96)^2)),scale=(((((Pemetrexed*1.2)-Pemetrexed)/1.96)^2)/Pemetrexed))
  Carboplatin.L=qgamma(runif(1000),shape=((Carboplatin^2)/((((Carboplatin*1.2)-Carboplatin)/1.96)^2)),scale=(((((Carboplatin*1.2)-Carboplatin)/1.96)^2)/Carboplatin))
  Cisplatin.L=qgamma(runif(1000),shape=((Cisplatin^2)/((((Cisplatin*1.2)-Cisplatin)/1.96)^2)),scale=(((((Cisplatin*1.2)-Cisplatin)/1.96)^2)/Cisplatin))
  Docetaxel.L=qgamma(runif(1000),shape=((Docetaxel^2)/((((Docetaxel*1.2)-Docetaxel)/1.96)^2)),scale=(((((Docetaxel*1.2)-Docetaxel)/1.96)^2)/Docetaxel))
  Nivolumab.L=qgamma(runif(1000),shape=((Nivolumab^2)/((((Nivolumab*1.2)-Nivolumab)/1.96)^2)),scale=(((((Nivolumab*1.2)-Nivolumab)/1.96)^2)/Nivolumab))
  Lazert.Amivan.L=qgamma(runif(1000),shape=((Lazert.Amivan^2)/((((Lazert.Amivan*1.2)-Lazert.Amivan)/1.96)^2)),scale=(((((Lazert.Amivan*1.2)-Lazert.Amivan)/1.96)^2)/Lazert.Amivan))
  
  Consultation.L=qgamma(runif(1000),shape=((Consultation^2)/((((Consultation*1.2)-Consultation)/1.96)^2)),scale=(((((Consultation*1.2)-Consultation)/1.96)^2)/Consultation))
  Office.Visit.L=qgamma(runif(1000),shape=((Office.Visit^2)/((((Office.Visit*1.2)-Office.Visit)/1.96)^2)),scale=(((((Office.Visit*1.2)-Office.Visit)/1.96)^2)/Office.Visit))
  CT.CAP.with.contrast.L=qgamma(runif(1000),shape=((CT.CAP.with.contrast^2)/((((CT.CAP.with.contrast*1.2)-CT.CAP.with.contrast)/1.96)^2)),scale=(((((CT.CAP.with.contrast*1.2)-CT.CAP.with.contrast)/1.96)^2)/CT.CAP.with.contrast))
  MRI.brain.L=qgamma(runif(1000),shape=((MRI.brain^2)/((((MRI.brain*1.2)-MRI.brain)/1.96)^2)),scale=(((((MRI.brain*1.2)-MRI.brain)/1.96)^2)/MRI.brain))
  PET.Body.L=qgamma(runif(1000),shape=((PET.Body^2)/((((PET.Body*1.2)-PET.Body)/1.96)^2)),scale=(((((PET.Body*1.2)-PET.Body)/1.96)^2)/PET.Body))
  Electrolytes.L=qgamma(runif(1000),shape=((Electrolytes^2)/((((Electrolytes*1.2)-Electrolytes)/1.96)^2)),scale=(((((Electrolytes*1.2)-Electrolytes)/1.96)^2)/Electrolytes))
  Liver.L=qgamma(runif(1000),shape=((Liver^2)/((((Liver*1.2)-Liver)/1.96)^2)),scale=(((((Liver*1.2)-Liver)/1.96)^2)/Liver))
  Renal.L=qgamma(runif(1000),shape=((Renal^2)/((((Renal*1.2)-Renal)/1.96)^2)),scale=(((((Renal*1.2)-Renal)/1.96)^2)/Renal))
  Hematology.L=qgamma(runif(1000),shape=((Hematology^2)/((((Hematology*1.2)-Hematology)/1.96)^2)),scale=(((((Hematology*1.2)-Hematology)/1.96)^2)/Hematology))
  NGS.L=qgamma(runif(1000),shape=((NGS^2)/((((NGS*1.2)-NGS)/1.96)^2)),scale=(((((NGS*1.2)-NGS)/1.96)^2)/NGS))
  Bronchoscopy.L=qgamma(runif(1000),shape=((Bronchoscopy^2)/((((Bronchoscopy*1.2)-Bronchoscopy)/1.96)^2)),scale=(((((Bronchoscopy*1.2)-Bronchoscopy)/1.96)^2)/Bronchoscopy))
  Endoscopic.lung.resection.L=qgamma(runif(1000),shape=((Endoscopic.lung.resection^2)/((((Endoscopic.lung.resection*1.2)-Endoscopic.lung.resection)/1.96)^2)),scale=(((((Endoscopic.lung.resection*1.2)-Endoscopic.lung.resection)/1.96)^2)/Endoscopic.lung.resection))
  hospital.L=qgamma(runif(1000),shape=((hospital^2)/((((hospital*1.2)-hospital)/1.96)^2)),scale=(((((hospital*1.2)-hospital)/1.96)^2)/hospital))
  Radiation.1.L=qgamma(runif(1000),shape=((Radiation.1^2)/((((Radiation.1*1.2)-Radiation.1)/1.96)^2)),scale=(((((Radiation.1*1.2)-Radiation.1)/1.96)^2)/Radiation.1))		
  Radiation.2.L=qgamma(runif(1000),shape=((Radiation.2^2)/((((Radiation.2*1.2)-Radiation.2)/1.96)^2)),scale=(((((Radiation.2*1.2)-Radiation.2)/1.96)^2)/Radiation.2))		
  Radiation.3.L=qgamma(runif(1000),shape=((Radiation.3^2)/((((Radiation.3*1.2)-Radiation.3)/1.96)^2)),scale=(((((Radiation.3*1.2)-Radiation.3)/1.96)^2)/Radiation.3))
  BSC.cost.L=qgamma(runif(1000),shape=((BSC.cost^2)/((((BSC.cost*1.2)-BSC.cost)/1.96)^2)),scale=(((((BSC.cost*1.2)-BSC.cost)/1.96)^2)/BSC.cost))
  
  
  u.s1.s=0.82
  u.s1.r=0.76
  u.s1.dfs=0.82
  
  u.s2.cs=0.77
  u.s2.r=0.76
  u.s2.dfs=0.80
  
  u.s3.cr=0.76
  u.s3.cs=0.78
  u.s3.dfs=0.79
  
  u.H.1L=0.71
  u.L.1L=0.71
  u.A.1L=0.80
  u.H.2L=0.67
  u.L.2L=0.67
  u.A.2L=0.64
  u.BSC=0.473
  
  
  
  u.s1.s.sd=0.16
  u.s1.r.sd=0.20
  u.s1.dfs.sd=0.17
  
  u.s2.cs.sd=0.18
  u.s2.r.sd=0.19
  u.s2.dfs.sd=0.17
  
  u.s3.cr.sd=0.18
  u.s3.cs.sd=0.18
  u.s3.dfs.sd=0.19
  
  u.H.1L.sd=0.24
  u.L.1L.sd=0.24
  u.A.1L.sd=(0.07/(2*1.96))*sqrt(81)
  u.H.2L.sd=0.2
  u.L.2L.sd=0.2
  u.A.2L.sd=(0.64*1.2-0.64)
  u.BSC.sd=0.02169
  
  a.s1.s=u.s1.s*(((u.s1.s*(1-u.s1.s))/(u.s1.s.sd^2))-1)
  b.s1.s=a.s1.s*((1/u.s1.s)-1)
  
  a.s1.r=u.s1.r*(((u.s1.r*(1-u.s1.r))/(u.s1.r.sd^2))-1)
  b.s1.r=a.s1.r*((1/u.s1.r)-1)
  
  a.s1.d=u.s1.dfs*(((u.s1.dfs*(1-u.s1.dfs))/(u.s1.dfs.sd^2))-1)
  b.s1.d=a.s1.d*((1/u.s1.dfs)-1)
  
  
  
  a.s2.cs=u.s2.cs*(((u.s2.cs*(1-u.s2.cs))/(u.s2.cs.sd^2))-1)
  b.s2.cs=a.s2.cs*((1/u.s2.cs)-1)
  
  a.s2.r=u.s2.r*(((u.s2.r*(1-u.s2.r))/(u.s2.r.sd^2))-1)
  b.s2.r=a.s2.r*((1/u.s2.r)-1)
  
  a.s2.d=u.s2.dfs*(((u.s2.dfs*(1-u.s2.dfs))/(u.s2.dfs.sd^2))-1)
  b.s2.d=a.s2.d*((1/u.s2.dfs)-1)
  
  
  
  a.s3.cr=u.s3.cr*(((u.s3.cr*(1-u.s3.cr))/(u.s3.cr.sd^2))-1)
  b.s3.cr=a.s3.cr*((1/u.s3.cr)-1)
  
  a.s3.cs=u.s3.cs*(((u.s3.cs*(1-u.s3.cs))/(u.s3.cs.sd^2))-1)
  b.s3.cs=a.s3.cs*((1/u.s3.cs)-1)
  
  a.s3.d=u.s3.dfs*(((u.s3.dfs*(1-u.s3.dfs))/(u.s3.dfs.sd^2))-1)
  b.s3.d=a.s3.d*((1/u.s3.dfs)-1)
  
  a.A.1L=u.A.1L*(((u.A.1L*(1-u.A.1L))/(u.A.1L.sd^2))-1)
  b.A.1L=a.A.1L*((1/u.A.1L)-1)
  
  a.H.1L=u.H.1L*(((u.H.1L*(1-u.H.1L))/(u.H.1L.sd^2))-1)
  b.H.1L=a.H.1L*((1/u.H.1L)-1)
  
  a.L.1L=u.L.1L*(((u.L.1L*(1-u.L.1L))/(u.L.1L.sd^2))-1)
  b.L.1L=a.L.1L*((1/u.L.1L)-1)
  
  a.A.2L=u.A.2L*(((u.A.2L*(1-u.A.2L))/(u.A.2L.sd^2))-1)
  b.A.2L=a.A.2L*((1/u.A.2L)-1)
  
  a.H.2L=u.H.2L*(((u.H.2L*(1-u.H.2L))/(u.H.2L.sd^2))-1)
  b.H.2L=a.H.2L*((1/u.H.2L)-1)
  
  a.L.2L=u.L.2L*(((u.L.2L*(1-u.L.2L))/(u.L.2L.sd^2))-1)
  b.L.2L=a.L.2L*((1/u.L.2L)-1)
  
  a.BSC=u.BSC*(((u.BSC*(1-u.BSC))/(u.BSC.sd^2))-1)
  b.BSC=a.BSC*((1/u.BSC)-1)
  
  
  
  u.s1.s.L=qbeta(runif(1000), shape1=a.s1.s, shape2=b.s1.s)/12
  u.s1.r.L=qbeta(runif(1000), shape1=a.s1.r, shape2=b.s1.r)/12
  u.s1.dfs.L=qbeta(runif(1000), shape1=a.s1.d, shape2=b.s1.d)/12
  u.s2.cs.L=qbeta(runif(1000), shape1=a.s2.cs, shape2=b.s2.cs)/12
  u.s2.r.L=qbeta(runif(1000), shape1=a.s2.r, shape2=b.s2.r)/12
  u.s2.dfs.L=qbeta(runif(1000), shape1=a.s2.d, shape2=b.s2.d)/12
  u.s3.cr.L=qbeta(runif(1000), shape1=a.s3.cr, shape2=b.s3.cr)/12
  u.s3.cs.L=qbeta(runif(1000), shape1=a.s3.cs, shape2=b.s3.cs)/12
  u.s3.dfs.L=qbeta(runif(1000), shape1=a.s3.d, shape2=b.s3.d)/12
  u.H.1L.L=qbeta(runif(1000), shape1=a.H.1L, shape2=b.H.1L)/12
  u.L.1L.L=qbeta(runif(1000), shape1=a.L.1L, shape2=b.L.1L)/12
  u.A.1L.L=qbeta(runif(1000), shape1=a.A.1L, shape2=b.A.1L)/12
  u.H.2L.L=qbeta(runif(1000), shape1=a.H.2L, shape2=b.H.2L)/12
  u.L.2L.L=qbeta(runif(1000), shape1=a.L.2L, shape2=b.L.2L)/12
  u.A.2L.L=qbeta(runif(1000), shape1=a.A.2L, shape2=b.A.2L)/12
  u.BSC.L=qbeta(runif(1000), shape1=a.BSC, shape2=b.BSC)/12
  
  if (mari==FALSE){
    u.A.1L.disut=0.039
    u.A.1L.disut.sd=(((0.039*1.2)-0.039)/1.96)}
  
  if (mari==TRUE){
    u.A.1L.disut=0.019
    u.A.1L.disut.sd=(((0.019*1.2)-0.019)/1.96)}
  
  
  a.disut=u.A.1L.disut*(((u.A.1L.disut*(1-u.A.1L.disut))/(u.A.1L.disut.sd^2))-1)
  b.disut=a.disut*((1/u.A.1L.disut)-1)
  
  
  u.A.1L.comp.L=((u.A.1L.L*12)-qbeta(runif(1000), shape1=a.disut, shape2=b.disut))/12
  
  
  parms=read.csv(paste0(getwd(),"/NSCLC microsim data/parms_comb.csv"))
  parms.samp.m=Probs.sample(1000,parms)
  
  
  mean.os=0.75
  lse.os=(log(0.92)-log(0.61))/(2*1.96)
  
  mean.pfs=0.76
  lse.pfs=(log(0.90)-log(0.64))/(2*1.96)
  
  corr = 0.641
  
  mu=c(log(mean.pfs), log(mean.os))
  
  Sigma=matrix(c(lse.pfs^2, corr* lse.pfs * lse.os, corr * lse.pfs * lse.os,  lse.os^2), nrow = 2, byrow = TRUE)
  
  samp.HR.MARIP=exp(mvrnorm(n = 1000, mu = mu, Sigma = Sigma))
  
  hr.prob.MARIP=data.frame(os.hr=samp.HR.MARIP[,2],pfs.hr=samp.HR.MARIP[,1])
  
  
  mean.os=0.77
  lse.os=(log(0.96)-log(0.61))/(2*1.96)
  
  mean.pfs=0.73
  lse.pfs=(log(0.89)-log(0.60))/(2*1.96)
  
  corr = 0.641
  
  mu=c(log(mean.pfs), log(mean.os))
  
  Sigma=matrix(c(lse.pfs^2, corr* lse.pfs * lse.os, corr * lse.pfs * lse.os,  lse.os^2), nrow = 2, byrow = TRUE)
  
  samp.HR.FLAURA2=exp(mvrnorm(n = 1000, mu = mu, Sigma = Sigma))
  
  hr.prob.FLAURA2=data.frame(os.hr=samp.HR.FLAURA2[,2],pfs.hr=samp.HR.FLAURA2[,1])
  
  
  c.mat.avg=data.frame()
  u.mat.avg=data.frame()
  c.mat.tot=data.frame()
  u.mat.tot=data.frame()
  check=data.frame()
  surv1y=data.frame()
  surv3y=data.frame()
  surv5y=data.frame()
  surv10y=data.frame()
  c.egfr=data.frame()
  u.egfr=data.frame()
  for (j in start_i:end_i){
    #model parameters
    
    n.m=end_i-start_i
    n.i=10000
    n.t=300
    
    
    #Base Costs
    Osimertinib=Osimertinib.L[j]
    Alectinib=Alectinib.L[j]
    Lorlatinib=Lorlatinib.L[j]
    Crizotinib=Crizotinib.L[j]
    Pembrolizumab=Pembrolizumab.L[j]
    Durvalumab=Durvalumab.L[j]
    Atezolizumab=Atezolizumab.L[j]
    Pemetrexed=Pemetrexed.L[j]
    Carboplatin=Carboplatin.L[j]
    Cisplatin=Cisplatin.L[j]
    Docetaxel=Docetaxel.L[j]
    Pembrolizumab.Carboplatin.Pemetrexed=Pembrolizumab+Carboplatin+Pemetrexed
    Pembrolizumab.Pemetrexed=Pembrolizumab+Pemetrexed
    Cisplatin.Pemetrexed=Cisplatin+Pemetrexed
    Carbo.Pemetrexed=Carboplatin+Pemetrexed
    Cisplatin.Docetaxel=Cisplatin+Docetaxel
    Nivolumab=Nivolumab.L[j]
    Lazert.Amivan=Lazert.Amivan.L[j]
    
    Consultation=Consultation.L[j]
    Office.Visit=Office.Visit.L[j]
    CT.CAP.with.contrast=CT.CAP.with.contrast.L[j]
    MRI.brain=MRI.brain.L[j]
    PET.Body=PET.Body.L[j]
    Chemistry=Electrolytes.L[j]+Liver.L[j]+Renal.L[j]
    Hematology=Hematology.L[j]
    NGS=NGS.L[j]
    Bronchoscopy=Bronchoscopy.L[j]
    Endoscopic.lung.resection=Endoscopic.lung.resection.L[j]
    Radiation.1=Radiation.1.L[j]	
    Radiation.2=Radiation.2.L[j]		
    Radiation.3=Radiation.3.L[j]
    BSC.cost=BSC.cost.L[j]
    hospital=hospital.L[j]
    
    #Costing
    
    Surveilance=Office.Visit+(CT.CAP.with.contrast*2)
    Surveilance.4=Office.Visit+(CT.CAP.with.contrast*4)
    stage1.cost.misc=(Consultation*1)+(Office.Visit*2)+(CT.CAP.with.contrast*1)+(MRI.brain*1)+(Chemistry*2)+(Hematology*2)+(Bronchoscopy*1)
    stage2.cost.misc=(Consultation*1)+(Office.Visit*5)+(CT.CAP.with.contrast*1)+(MRI.brain*1)+(PET.Body*1)+(Chemistry*6)+(Hematology*6)+(Bronchoscopy*1)
    stage3.cost.misc=(Consultation*1)+(Office.Visit*17)+(CT.CAP.with.contrast*4)+(MRI.brain*1)+(PET.Body*1)+(Chemistry*20)+(Hematology*20)+(NGS*1)+(Bronchoscopy*1)
    stage4.cost.misc=(Consultation*1)+(Office.Visit*17)+(CT.CAP.with.contrast*4)+(MRI.brain*1)+(PET.Body*1)+(Chemistry*20)+(Hematology*20)+(NGS*1)+(Bronchoscopy*1)
    stage42L.cost.misc=(Office.Visit*17)+(CT.CAP.with.contrast*4)+(Hematology*17)
    stage42L.cost.misc.act=(Office.Visit*17)+(CT.CAP.with.contrast*4)+(Hematology*17)+(MRI.brain*2)
    #Utilities
    
    u.s1.s=u.s1.s.L[j]
    u.s1.r=u.s1.r.L[j]
    u.s1.dfs=u.s1.dfs.L[j]
    
    u.s2.cs=u.s2.cs.L[j]
    u.s2.r=u.s2.r.L[j]
    u.s2.dfs=u.s2.dfs.L[j]
    
    u.s3.cr=u.s3.cr.L[j]
    u.s3.cs=u.s3.cs.L[j]
    u.s3.dfs=u.s3.dfs.L[j]
    
    u.H.1L=u.H.1L.L[j]
    u.L.1L=u.L.1L.L[j]
    u.A.1L=u.A.1L.L[j]
    u.H.2L=u.H.2L.L[j]
    u.L.2L=u.L.2L.L[j]
    u.A.2L=u.A.2L.L[j]
    u.BSC=u.BSC.L[j]
    u.A.1L.comp=u.A.1L.comp.L[j]
    
    #u.s1=1/12
    #u.s1.dfs=1/12
    #u.s2=1/12
    #u.s2.dfs=1/12
    #u.s3=1/12
    #u.s3.dfs=1/12
    #u.H.1L=1/12
    #u.L.1L=1/12
    #u.A.1L=1/12
    #u.H.2L=1/12
    #u.L.2L=1/12
    #u.A.2L=1/12
    #u.BSC=1/12
    
    
    
    
    
    #transition probabilities
    # normal staging
    s4=(1193/2410)
    s1=(633/2410)
    s2=(170/2410)
    s3=(414/2410)
    
    # pancam staging
    #s1=100/152
    #s2=14/152
    #s3=21/152
    #s4=17/152
    
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
    #checks
    s4+s1+s2+s3==1
    alk+egfr+low+high==1
    
    
    df.list=Probs.PSA(low, high, alk, egfr, parms.samp.m, hr.prob.FLAURA2, hr.prob.MARIP, j)
    
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
    
    #dataframes
    
    t.df=data.frame()
    
    n.rows=301
    s.mat=matrix(NA, nrow = n.rows, ncol = n.i)
    c.mat=matrix(NA, nrow = n.rows, ncol = n.i)
    u.mat=matrix(NA, nrow = n.rows, ncol = n.i)
    c.mat.alt=matrix(NA, nrow = n.rows, ncol = n.i)
    u.mat.alt=matrix(NA, nrow = n.rows, ncol = n.i)
    
    #discount function
    d.rt=0.015
    discount=function(unit,t,d.rt){unit/((1+d.rt)^((t-1)/12))}
    
    c.tot=0
    u.tot=0
    
    #begin simulation
    
    for (i in 1:n.i){
      set.seed(1 + i) #set seed for each individual
      
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
            if (s1.ct>1 & s1.ct<=62){c.temp=c.temp+(Surveilance/12)}
            
            if (s1.ct<=4){u.temp=u.temp+u.s1.s}
            if (s1.ct>4){u.temp=u.temp+u.s1.dfs}
            
          }
          if (s1.treat=="Rad"){
            if (s1.ct==1){c.temp=c.temp+((stage1.cost.misc+Radiation.1)/2)}
            if (s1.ct==2){c.temp=c.temp+((stage1.cost.misc+Radiation.1)/2)}
            if (s1.ct>2& s1.ct<=63){c.temp=c.temp+(Surveilance/12)}
            
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
        #c.mat[t,i]=c.ind
        #u.mat[t,i]=u.ind
      }
      c.mat.alt[1,i]=c.ind
      u.mat.alt[1,i]=u.ind
      
      
      c.tot=c.tot+c.ind
      u.tot=u.tot+u.ind
    }
    #if (j.ind==1){
    #s.mat.comb=s.mat
    #c.mat.comb=c.mat
    #u.mat.comb=u.mat
    #}
    #if (j.ind>1){
    #s.mat.comb=cbind(s.mat.comb,s.mat)
    #c.mat.comb=cbind(c.mat.comb,c.mat)
    #u.mat.comb=cbind(u.mat.comb,u.mat)
    #}
    c.mat.avg[1,j.ind]=mean(unlist(c.mat.alt[1,]))
    u.mat.avg[1,j.ind]=mean(unlist(u.mat.alt[1,]))
    c.mat.tot[1,j.ind]=c.tot
    u.mat.tot[1,j.ind]=u.tot
    surv1y[1,j.ind]=mean(s.mat[13,] != "Dead")
    surv3y[1,j.ind]=mean(s.mat[37,] != "Dead")
    surv5y[1,j.ind]=mean(s.mat[61,] != "Dead")
    surv10y[1,j.ind]=mean(s.mat[121,] != "Dead")
    c.egfr[1,j.ind]=mean(unlist(c.mat.alt[1,s.mat[1,]=="EGFR 1L"]))
    u.egfr[1,j.ind]=mean(unlist(u.mat.alt[1,s.mat[1,]=="EGFR 1L"]))
    
    check[1,j.ind]=u.A.1L*12
    check[2,j.ind]=u.A.1L.comp*12
    j.ind=j.ind+1
    
  }
  return(list(c.mat.avg=c.mat.avg,u.mat.avg=u.mat.avg,c.mat.tot=c.mat.tot,u.mat.tot=u.mat.tot,
              surv5y=surv5y,surv10y=surv10y,c.egfr=c.egfr,u.egfr=u.egfr,surv1y=surv1y,surv3y=surv3y,s.mat=s.mat))
}


source(paste0(getwd(),"/Parms_sample.R"))
source(paste0(getwd(),"/Probs_PSA.R"))


#export functions
clusterExport(cl_sock, varlist = c("NSCLCmicro","Probs.sample","Probs.PSA"))

#divide second order simulations into chunks
n=1000
n.chunks=15
chunk_size= n / n.chunks
chunks=split(1:n, ceiling(seq_along(1:n) / chunk_size))

#run simulations
system.time({
  soc.output=parLapply(cl_sock, chunks, function(x) {
    
    NSCLCmicro(min(x), max(x),comp=FALSE,mari=FALSE)
  })
}) 

system.time({
  comp.output=parLapply(cl_sock, chunks, function(x) {
    
    NSCLCmicro(min(x), max(x),comp=TRUE,mari=FALSE)
  })
})
system.time({
  comp.output2=parLapply(cl_sock, chunks, function(x) {
    
    NSCLCmicro(min(x), max(x),comp=TRUE,mari=TRUE)
  })
})

#stop cluster
stopCluster(cl = cl_sock);

#extract matrices

#soc.list1=lapply(soc.output, function(res) res$s.mat)
#soc.list2=lapply(soc.output, function(res) res$c.mat.comb)
#soc.list3=lapply(soc.output, function(res) res$u.mat.comb)
soc.list4=lapply(soc.output, function(res) res$c.mat.avg)
soc.list5=lapply(soc.output, function(res) res$u.mat.avg)
soc.list6=lapply(soc.output, function(res) res$c.mat.tot)
soc.list7=lapply(soc.output, function(res) res$u.mat.tot)
soc.list8=lapply(soc.output, function(res) res$check)
soc.list9=lapply(soc.output, function(res) res$surv5y)
soc.list10=lapply(soc.output, function(res) res$surv10y)
soc.list11=lapply(soc.output, function(res) res$c.egfr)
soc.list12=lapply(soc.output, function(res) res$u.egfr)
soc.list13=lapply(soc.output, function(res) res$surv1y)
soc.list14=lapply(soc.output, function(res) res$surv3y)

comp.list1=lapply(comp.output, function(res) res$c.mat.avg)
comp.list2=lapply(comp.output, function(res) res$u.mat.avg)
comp.list3=lapply(comp.output, function(res) res$c.mat.tot)
comp.list4=lapply(comp.output, function(res) res$u.mat.tot)
comp.list5=lapply(comp.output, function(res) res$check)
comp.list9=lapply(comp.output, function(res) res$surv5y)
comp.list10=lapply(comp.output, function(res) res$surv10y)
comp.list11=lapply(comp.output, function(res) res$c.egfr)
comp.list12=lapply(comp.output, function(res) res$u.egfr)
comp.list13=lapply(comp.output, function(res) res$surv1y)
comp.list14=lapply(comp.output, function(res) res$surv3y)


comp2.list1=lapply(comp.output2, function(res) res$c.mat.avg)
comp2.list2=lapply(comp.output2, function(res) res$u.mat.avg)
comp2.list3=lapply(comp.output2, function(res) res$c.mat.tot)
comp2.list4=lapply(comp.output2, function(res) res$u.mat.tot)
comp2.list5=lapply(comp.output2, function(res) res$check)
comp2.list9=lapply(comp.output2, function(res) res$surv5y)
comp2.list10=lapply(comp.output2, function(res) res$surv10y)
comp2.list11=lapply(comp.output2, function(res) res$c.egfr)
comp2.list12=lapply(comp.output2, function(res) res$u.egfr)
comp2.list13=lapply(comp.output2, function(res) res$surv1y)
comp2.list14=lapply(comp.output2, function(res) res$surv3y)

# Convert to data frames
#s.df=as.data.frame(do.call(cbind, soc.list1))
#c.df=as.data.frame(do.call(cbind, soc.list2))
#u.df=as.data.frame(do.call(cbind, soc.list3))
c.df.avg=as.data.frame(do.call(cbind, soc.list4))
u.df.avg=as.data.frame(do.call(cbind, soc.list5))
c.df.tot=as.data.frame(do.call(cbind, soc.list6))
u.df.tot=as.data.frame(do.call(cbind, soc.list7))
check=as.data.frame(do.call(cbind, soc.list8))
surv5y=as.data.frame(do.call(cbind, soc.list9))
surv10y=as.data.frame(do.call(cbind, soc.list10))
c.egfr=as.data.frame(do.call(cbind, soc.list11))
u.egfr=as.data.frame(do.call(cbind, soc.list12))
surv1y=as.data.frame(do.call(cbind, soc.list13))
surv3y=as.data.frame(do.call(cbind, soc.list14))

c.df.avg.comp=as.data.frame(do.call(cbind, comp.list1))
u.df.avg.comp=as.data.frame(do.call(cbind, comp.list2))
c.df.tot.comp=as.data.frame(do.call(cbind, comp.list3))
u.df.tot.comp=as.data.frame(do.call(cbind, comp.list4))
ckeck.comp=as.data.frame(do.call(cbind, comp.list5))
surv5y.comp=as.data.frame(do.call(cbind, comp.list9))
surv10y.comp=as.data.frame(do.call(cbind, comp.list10))
c.egfr.comp=as.data.frame(do.call(cbind, comp.list11))
u.egfr.comp=as.data.frame(do.call(cbind, comp.list12))
surv1y.comp=as.data.frame(do.call(cbind, comp.list13))
surv3y.comp=as.data.frame(do.call(cbind, comp.list14))

c.df.avg.comp2=as.data.frame(do.call(cbind, comp2.list1))
u.df.avg.comp2=as.data.frame(do.call(cbind, comp2.list2))
c.df.tot.comp2=as.data.frame(do.call(cbind, comp2.list3))
u.df.tot.comp2=as.data.frame(do.call(cbind, comp2.list4))
ckeck.comp2=as.data.frame(do.call(cbind, comp2.list5))
surv5y.comp2=as.data.frame(do.call(cbind, comp2.list9))
surv10y.comp2=as.data.frame(do.call(cbind, comp2.list10))
c.egfr.comp2=as.data.frame(do.call(cbind, comp2.list11))
u.egfr.comp2=as.data.frame(do.call(cbind, comp2.list12))
surv1y.comp2=as.data.frame(do.call(cbind, comp2.list13))
surv3y.comp2=as.data.frame(do.call(cbind, comp2.list14))




