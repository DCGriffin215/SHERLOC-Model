Probs.sample.het=function(n_sim,parms){
  library(MASS)
  parms.samp=matrix(nrow=n_sim,ncol=47)
  
  ############################ Stage I ############################
  
  temp.cov.s1D=read.csv(paste0(getwd(),"/SHERLOC microsim Het Data/cov s1_D.csv"))
  temp.cov.s1s4=read.csv(paste0(getwd(),"/SHERLOC microsim Het Data/cov s1_s4.csv"))
  
  log.s1_D.p1=log(parms[1,2])
  
  op.par=c(log.s1_D.p1)
  
  mvnorm.res=mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.s1D)
  
  v.s1_D.r=exp(mvnorm.res[,1])
  
  parms.samp[,1]=v.s1_D.r
  ####
  s1_s4.p1=parms[2,2]
  log.s1_s4.p2=log(parms[3,2])
  B=parms[4,2]
  
  op.par=c(s1_s4.p1,log.s1_s4.p2, B)
  
  mvnorm.res=mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.s1s4)
  
  v.s1_s4.p1=mvnorm.res[,1]
  v.s1_s4.p2=exp(mvnorm.res[,2])
  v.s1_s4.p3=mvnorm.res[,3]
  
  parms.samp[,2]=v.s1_s4.p1
  parms.samp[,3]=v.s1_s4.p2
  parms.samp[,4]=v.s1_s4.p3
  
  ############################ Stage II ############################
  
  temp.cov.s2D=read.csv(paste0(getwd(),"/SHERLOC microsim Het Data/cov s2_D.csv"))
  temp.cov.s2s4=read.csv(paste0(getwd(),"/SHERLOC microsim Het Data/cov s2_s4.csv"))
  
  log.s2_D.p1=log(parms[5,2])
  op.par=c(log.s2_D.p1)
  
  mvnorm.res=mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.s2D)
  
  v.s2_D=exp(mvnorm.res[,1])
  
  parms.samp[,5]=v.s2_D
  ###
  s2_s4.p1=parms[6,2]
  log.s2_s4.p2=log(parms[7,2])
  s2_s4.p3=parms[8,2]
  B=parms[9:10,2]
  
  op.par=c(s2_s4.p1,log.s2_s4.p2,s2_s4.p3, B)
  
  mvnorm.res=mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.s2s4)
  
  v.s1_s4.p1=mvnorm.res[,1]
  v.s1_s4.p2=exp(mvnorm.res[,2])
  v.s1_s4.p3=mvnorm.res[,3]
  v.s1_s4.p4=mvnorm.res[,4]
  v.s1_s4.p5=mvnorm.res[,5]
  
  parms.samp[,6]=v.s1_s4.p1
  parms.samp[,7]=v.s1_s4.p2
  parms.samp[,8]=v.s1_s4.p3
  parms.samp[,9]=v.s1_s4.p4
  parms.samp[,10]=v.s1_s4.p5
  
  ############################ Stage III ############################
  
  temp.cov.s3D=read.csv(paste0(getwd(),"/SHERLOC microsim Het Data/cov s3_D.csv"))
  temp.cov.s3s4=read.csv(paste0(getwd(),"/SHERLOC microsim Het Data/cov s3_s4.csv"))
  
  s3_D.p1=parms[11,2]
  log.s3_D.p2=log(parms[12,2])
  
  op.par=c(s3_D.p1,log.s3_D.p2)
  
  mvnorm.res=mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.s3D)
  
  v.s3_D.p1=mvnorm.res[,1]
  v.s3_D.p2=exp(mvnorm.res[,2])
  
  parms.samp[,11]=v.s3_D.p1
  parms.samp[,12]=v.s3_D.p2
  ###
  s3_s4.p1=parms[13,2]
  log.s3_s4.p2=log(parms[14,2])
  s3_s4.p3=parms[15,2]
  
  op.par=c(s3_s4.p1,log.s3_s4.p2,s3_s4.p3)
  
  mvnorm.res=mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.s3s4)
  
  v.s3_s4.p1=mvnorm.res[,1]
  v.s3_s4.p2=exp(mvnorm.res[,2])
  v.s3_s4.p3=mvnorm.res[,3]
  
  parms.samp[,13]=v.s3_s4.p1
  parms.samp[,14]=v.s3_s4.p2
  parms.samp[,15]=v.s3_s4.p3
  
  
  ############################ EGFR 1L ############################
  
  temp.cov.e1LD=read.csv(paste0(getwd(),"/SHERLOC microsim Het Data/cov EGFR-1L_D.csv"))
  temp.cov.e1L2L=read.csv(paste0(getwd(),"/SHERLOC microsim Het Data/cov EGFR-1L_2L.csv"))
  
  log.egfr_D.p1=log(parms[16,2])
  B=parms[17,2]
  op.par=c(log.egfr_D.p1, B)
  
  mvnorm.res=mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.e1LD)
  
  v.egfr_D=exp(mvnorm.res[,1])
  v.egfr_D.2=mvnorm.res[,2]
  
  parms.samp[,16]=v.egfr_D
  parms.samp[,17]=v.egfr_D.2
  ###
  
  log.egfr_2L.p1=log(parms[18,2])
  B=parms[19,2]
  
  op.par=c(log.egfr_2L.p1, B)
  
  mvnorm.res=mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.e1L2L)
  
  v.egfr_2L=exp(mvnorm.res[,1])
  v.egfr_2L.2=mvnorm.res[,2]
  
  parms.samp[,18]=v.egfr_2L
  parms.samp[,19]=v.egfr_2L.2
  
  ############################ ALK 1L ############################
  
  temp.cov.a1LD=read.csv(paste0(getwd(),"/SHERLOC microsim Het Data/cov ALK-1L_D.csv"))
  temp.cov.a1L2L=read.csv(paste0(getwd(),"/SHERLOC microsim Het Data/cov ALK-1L_2L.csv"))
  
  log.alk_D.p1=log(parms[20,2])
  op.par=c(log.alk_D.p1)
  
  mvnorm.res=mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.a1LD)
  
  v.alk_D=exp(mvnorm.res[,1])
  
  parms.samp[,20]=v.alk_D
  ####
  alk_2L.p1=parms[21,2]
  log.alk_2L.p2=log(parms[22,2])
  
  op.par=c(alk_2L.p1,log.alk_2L.p2)
  
  mvnorm.res=mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.a1L2L)
  
  v.alk_2L.p1=mvnorm.res[,1]
  v.alk_2L.p2=exp(mvnorm.res[,2])
  
  parms.samp[,21]=v.alk_2L.p1
  parms.samp[,22]=v.alk_2L.p2
  
  
  ############################ PDL1 High 1L ############################
  
  temp.cov.high1LD=read.csv(paste0(getwd(),"/SHERLOC microsim Het Data/cov PDL1H-1L_D.csv"))
  temp.cov.high1L2L=read.csv(paste0(getwd(),"/SHERLOC microsim Het Data/cov PDL1H-1L_2L.csv"))
  
  high_D.p1=parms[23,2]
  log.high_D.p2=log(parms[24,2])
  B=parms[25:27,2]
  
  op.par=c(high_D.p1,log.high_D.p2, B)
  
  mvnorm.res=mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.high1LD)
  
  v.high_D.p1=mvnorm.res[,1]
  v.high_D.p2=exp(mvnorm.res[,2])
  v.high_D.p3=mvnorm.res[,3]
  v.high_D.p4=mvnorm.res[,4]
  v.high_D.p5=mvnorm.res[,5]
  
  parms.samp[,23]=v.high_D.p1
  parms.samp[,24]=v.high_D.p2
  parms.samp[,25]=v.high_D.p3
  parms.samp[,26]=v.high_D.p4
  parms.samp[,27]=v.high_D.p5
  ####
  high_2L.p1=parms[28,2]
  log.high_2L.p2=log(parms[29,2])
  B=parms[30:31,2]
  
  op.par=c(high_2L.p1,log.high_2L.p2, B)
  
  mvnorm.res=mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.high1L2L)
  
  v.high_2L.p1=mvnorm.res[,1]
  v.high_2L.p2=exp(mvnorm.res[,2])
  v.high_2L.p3=mvnorm.res[,3]
  v.high_2L.p4=mvnorm.res[,4]
  
  parms.samp[,28]=v.high_2L.p1
  parms.samp[,29]=v.high_2L.p2
  parms.samp[,30]=v.high_2L.p3
  parms.samp[,31]=v.high_2L.p4
  
  
  ############################ PDL1 low 1L ############################
  
  temp.cov.low1LD=read.csv(paste0(getwd(),"/SHERLOC microsim Het Data/cov PDL1L-1L_D.csv"))
  temp.cov.low1L2L=read.csv(paste0(getwd(),"/SHERLOC microsim Het Data/cov PDL1L-1L_2L.csv"))
  
  low_D.p1=parms[32,2]
  log.low_D.p2=log(parms[33,2])
  B=parms[34:36,2]
  
  op.par=c(low_D.p1,log.low_D.p2, B)
  
  mvnorm.res=mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.low1LD)
  
  v.low_D.p1=mvnorm.res[,1]
  v.low_D.p2=exp(mvnorm.res[,2])
  v.low_D.p3=mvnorm.res[,3]
  v.low_D.p4=mvnorm.res[,4]
  v.low_D.p5=mvnorm.res[,5]
  
  parms.samp[,32]=v.low_D.p1
  parms.samp[,33]=v.low_D.p2
  parms.samp[,34]=v.low_D.p3
  parms.samp[,35]=v.low_D.p4
  parms.samp[,36]=v.low_D.p5
  
  #####
  low_2L.p1=parms[37,2]
  log.low_2L.p2=log(parms[38,2])
  B=parms[39,2]
  
  op.par=c(low_2L.p1,log.low_2L.p2, B)
  
  mvnorm.res=mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.low1L2L)
  
  v.low_2L.p1=mvnorm.res[,1]
  v.low_2L.p2=exp(mvnorm.res[,2])
  v.low_2L.p3=mvnorm.res[,3]
  
  parms.samp[,37]=v.low_2L.p1
  parms.samp[,38]=v.low_2L.p2
  parms.samp[,39]=v.low_2L.p3
  
  
  ############################ PDL1 low 2L ############################
  
  temp.cov.lowD=read.csv(paste0(getwd(),"/SHERLOC microsim Het Data/cov low-2L_D.csv"))
  
  log.low2L_D.p1=parms[40,2]
  log.low2L_D.p2=log(parms[41,2])
  B=parms[42:43,2]
  
  op.par=c(log.low2L_D.p1,log.low2L_D.p2, B)
  
  mvnorm.res=mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.lowD)
  
  v.low2L_D.p1=mvnorm.res[,1]
  v.low2L_D.p2=exp(mvnorm.res[,2])
  v.low2L_D.p3=mvnorm.res[,3]
  v.low2L_D.p4=mvnorm.res[,4]
  
  parms.samp[,40]=v.low2L_D.p1
  parms.samp[,41]=v.low2L_D.p2
  parms.samp[,42]=v.low2L_D.p3
  parms.samp[,43]=v.low2L_D.p4
  
  ############################ PDL1 high 2L ############################
  
  temp.cov.highD=read.csv(paste0(getwd(),"/SHERLOC microsim Het Data/cov high-2L_D.csv"))
  
  log.high2L_D.p1=log(parms[44,2])
  log.high2L_D.p2=log(parms[45,2])
  
  op.par=c(log.high2L_D.p1,log.high2L_D.p2)
  
  mvnorm.res=mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.highD)
  
  v.high2L_D.p1=exp(mvnorm.res[,1])
  v.high2L_D.p2=exp(mvnorm.res[,2])
  
  parms.samp[,44]=v.high2L_D.p1
  parms.samp[,45]=v.high2L_D.p2
  
  
  ############################ EGFR 2L ############################
  
  temp.cov.egfrD=read.csv(paste0(getwd(),"/SHERLOC microsim Het Data/cov EGFR-2L_D.csv"))
  
  log.egfr2L_D.p1=log(parms[46,2])
  
  op.par=c(log.egfr2L_D.p1)
  
  mvnorm.res=mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.egfrD)
  
  v.egfr2L_D.p1=exp(mvnorm.res[,1])
  
  parms.samp[,46]=v.egfr2L_D.p1
  
  
  ############################ ALK 2L ############################
  
  temp.cov.alkD=read.csv(paste0(getwd(),"/SHERLOC microsim Het Data/cov ALK-2L_D.csv"))
  
  log.alk2L_D.p1=log(parms[47,2])
  
  op.par=c(log.alk2L_D.p1)
  
  mvnorm.res=mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.alkD)
  
  v.alk2L_D.p1=exp(mvnorm.res[,1])
  
  parms.samp[,47]=v.alk2L_D.p1
  
  
  parms.samp=as.data.frame(t(parms.samp))
  
  return(parms.samp)
}




