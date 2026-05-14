


Probs.sample=function(n_sim,parms){
  library(MASS)
  
  parms.samp=matrix(nrow=n_sim,ncol=34)
  
  ############################ Stage I ############################
  
  temp.cov.s1D=read.csv(paste0(getwd(),"/NSCLC microsim data/cov s1_D.csv"))
  temp.cov.s1s4=read.csv(paste0(getwd(),"/NSCLC microsim data/cov s1_s4.csv"))
  
  log.s1_D.p1=log(parms[1,2])
  op.par=c(log.s1_D.p1)
  
  mvnorm.res <- mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.s1D)
  
  v.s1_D=exp(mvnorm.res[,1])
  
  parms.samp[,1]=v.s1_D
  
  
  s1_s4.p1=parms[2,2]
  log.s1_s4.p2=log(parms[3,2])
  
  
  op.par=c(s1_s4.p1,log.s1_s4.p2)
  
  mvnorm.res <- mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.s1s4)
  
  v.s1_s4.p1=mvnorm.res[,1]
  v.s1_s4.p2=exp(mvnorm.res[,2])
  
  parms.samp[,2]=v.s1_s4.p1
  parms.samp[,3]=v.s1_s4.p2
  
  ############################ Stage II ############################
  
  temp.cov.s2D=read.csv(paste0(getwd(),"/NSCLC microsim data/cov s2_D.csv"))
  temp.cov.s2s4=read.csv(paste0(getwd(),"/NSCLC microsim data/cov s2_s4.csv"))
  
  log.s2_D.p1=log(parms[4,2])
  op.par=c(log.s2_D.p1)
  
  mvnorm.res <- mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.s2D)
  
  v.s2_D=exp(mvnorm.res[,1])
  
  parms.samp[,4]=v.s2_D
  
  
  
  s2_s4.p1=parms[5,2]
  log.s2_s4.p2=log(parms[6,2])
  s2_s4.p3=parms[7,2]
  
  
  op.par=c(s2_s4.p1,log.s2_s4.p2,s2_s4.p3)
  
  mvnorm.res <- mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.s2s4)
  
  v.s1_s4.p1=mvnorm.res[,1]
  v.s1_s4.p2=exp(mvnorm.res[,2])
  v.s1_s4.p3=mvnorm.res[,3]
  
  parms.samp[,5]=v.s1_s4.p1
  parms.samp[,6]=v.s1_s4.p2
  parms.samp[,7]=v.s1_s4.p3
  
  
  
  ############################ Stage III ############################
  
  temp.cov.s3D=read.csv(paste0(getwd(),"/NSCLC microsim data/cov s3_D.csv"))
  temp.cov.s3s4=read.csv(paste0(getwd(),"/NSCLC microsim data/cov s3_s4.csv"))
  
  s3_D.p1=parms[8,2]
  log.s3_D.p2=log(parms[9,2])
  
  
  op.par=c(s3_D.p1,log.s3_D.p2)
  
  mvnorm.res <- mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.s3D)
  
  v.s3_D.p1=mvnorm.res[,1]
  v.s3_D.p2=exp(mvnorm.res[,2])
  
  parms.samp[,8]=v.s3_D.p1
  parms.samp[,9]=v.s3_D.p2
  
  
  
  
  s3_s4.p1=parms[10,2]
  log.s3_s4.p2=log(parms[11,2])
  s3_s4.p3=parms[12,2]
  
  
  op.par=c(s3_s4.p1,log.s3_s4.p2,s3_s4.p3)
  
  mvnorm.res <- mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.s3s4)
  
  v.s3_s4.p1=mvnorm.res[,1]
  v.s3_s4.p2=exp(mvnorm.res[,2])
  v.s3_s4.p3=mvnorm.res[,3]
  
  parms.samp[,10]=v.s3_s4.p1
  parms.samp[,11]=v.s3_s4.p2
  parms.samp[,12]=v.s3_s4.p3
  
  
  
  ############################ EGFR 1L ############################
  
  temp.cov.e1LD=read.csv(paste0(getwd(),"/NSCLC microsim data/cov EGFR-1L_D.csv"))
  temp.cov.e1L2L=read.csv(paste0(getwd(),"/NSCLC microsim data/cov EGFR-1L_2L.csv"))
  
  log.egfr_D.p1=log(parms[13,2])
  op.par=c(log.egfr_D.p1)
  
  mvnorm.res <- mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.e1LD)
  
  v.egfr_D=exp(mvnorm.res[,1])
  
  parms.samp[,13]=v.egfr_D
  
  
  
  log.egfr_2L.p1=log(parms[14,2])
  op.par=c(log.egfr_2L.p1)
  
  mvnorm.res <- mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.e1L2L)
  
  v.egfr_2L=exp(mvnorm.res[,1])
  
  parms.samp[,14]=v.egfr_2L
  
  
  
  ############################ ALK 1L ############################
  
  temp.cov.a1LD=read.csv(paste0(getwd(),"/NSCLC microsim data/cov ALK-1L_D.csv"))
  temp.cov.a1L2L=read.csv(paste0(getwd(),"/NSCLC microsim data/cov ALK-1L_2L.csv"))
  
  log.alk_D.p1=log(parms[15,2])
  op.par=c(log.alk_D.p1)
  
  mvnorm.res <- mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.a1LD)
  
  v.alk_D=exp(mvnorm.res[,1])
  
  parms.samp[,15]=v.alk_D
  
  
  
  alk_2L.p1=parms[16,2]
  log.alk_2L.p2=log(parms[17,2])
  
  op.par=c(alk_2L.p1,log.alk_2L.p2)
  
  mvnorm.res <- mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.a1L2L)
  
  v.alk_2L.p1=mvnorm.res[,1]
  v.alk_2L.p2=exp(mvnorm.res[,2])
  
  parms.samp[,16]=v.alk_2L.p1
  parms.samp[,17]=v.alk_2L.p2
  
  
  ############################ PDL1 High 1L ############################
  
  temp.cov.high1LD=read.csv(paste0(getwd(),"/NSCLC microsim data/cov PDL1H-1L_D.csv"))
  temp.cov.high1L2L=read.csv(paste0(getwd(),"/NSCLC microsim data/cov PDL1H-1L_2L.csv"))
  
  high_D.p1=parms[18,2]
  log.high_D.p2=log(parms[19,2])
  high_D.p3=parms[20,2]
  
  
  op.par=c(high_D.p1,log.high_D.p2,high_D.p3)
  
  mvnorm.res <- mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.high1LD)
  
  v.high_D.p1=mvnorm.res[,1]
  v.high_D.p2=exp(mvnorm.res[,2])
  v.high_D.p3=mvnorm.res[,3]
  
  parms.samp[,18]=v.high_D.p1
  parms.samp[,19]=v.high_D.p2
  parms.samp[,20]=v.high_D.p3
  
  
  high_2L.p1=parms[21,2]
  log.high_2L.p2=log(parms[22,2])
  
  op.par=c(high_2L.p1,log.high_2L.p2)
  
  mvnorm.res <- mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.high1L2L)
  
  v.high_2L.p1=mvnorm.res[,1]
  v.high_2L.p2=exp(mvnorm.res[,2])
  
  parms.samp[,21]=v.high_2L.p1
  parms.samp[,22]=v.high_2L.p2
  
  
  
  ############################ PDL1 low 1L ############################
  
  temp.cov.low1LD=read.csv(paste0(getwd(),"/NSCLC microsim data/cov PDL1L-1L_D.csv"))
  temp.cov.low1L2L=read.csv(paste0(getwd(),"/NSCLC microsim data/cov PDL1L-1L_2L.csv"))
  
  low_D.p1=parms[23,2]
  log.low_D.p2=log(parms[24,2])
  
  op.par=c(low_D.p1,log.low_D.p2)
  
  mvnorm.res <- mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.low1LD)
  
  v.low_D.p1=mvnorm.res[,1]
  v.low_D.p2=exp(mvnorm.res[,2])
  
  parms.samp[,23]=v.low_D.p1
  parms.samp[,24]=v.low_D.p2
  
  
  
  low_2L.p1=parms[25,2]
  log.low_2L.p2=log(parms[26,2])
  
  op.par=c(low_2L.p1,log.low_2L.p2)
  
  mvnorm.res <- mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.low1L2L)
  
  v.low_2L.p1=mvnorm.res[,1]
  v.low_2L.p2=exp(mvnorm.res[,2])
  
  parms.samp[,25]=v.low_2L.p1
  parms.samp[,26]=v.low_2L.p2
  
  
  ############################ PDL1 low 2L ############################
  
  temp.cov.lowD=read.csv(paste0(getwd(),"/NSCLC microsim data/cov low-2L_D.csv"))
  
  log.low2L_D.p1=parms[27,2]
  log.low2L_D.p2=log(parms[28,2])
  
  op.par=c(log.low2L_D.p1,log.low2L_D.p2)
  
  mvnorm.res <- mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.lowD)
  
  v.low2L_D.p1=mvnorm.res[,1]
  v.low2L_D.p2=exp(mvnorm.res[,2])
  
  parms.samp[,27]=v.low2L_D.p1
  parms.samp[,28]=v.low2L_D.p2
  
  
  
  ############################ PDL1 high 2L ############################
  
  temp.cov.highD=read.csv(paste0(getwd(),"/NSCLC microsim data/cov high-2L_D.csv"))
  
  log.high2L_D.p1=log(parms[29,2])
  log.high2L_D.p2=log(parms[30,2])
  
  op.par=c(log.high2L_D.p1,log.high2L_D.p2)
  
  mvnorm.res <- mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.highD)
  
  v.high2L_D.p1=exp(mvnorm.res[,1])
  v.high2L_D.p2=exp(mvnorm.res[,2])
  
  parms.samp[,29]=v.high2L_D.p1
  parms.samp[,30]=v.high2L_D.p2
  
  
  ############################ EGFR 2L ############################
  
  temp.cov.egfrD=read.csv(paste0(getwd(),"/NSCLC microsim data/cov EGFR-2L_D.csv"))
  
  log.egfr2L_D.p1=log(parms[31,2])
  log.egfr2L_D.p2=log(parms[32,2])
  
  op.par=c(log.egfr2L_D.p1,log.egfr2L_D.p2)
  
  mvnorm.res <- mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.egfrD)
  
  v.egfr2L_D.p1=exp(mvnorm.res[,1])
  v.egfr2L_D.p2=exp(mvnorm.res[,2])
  
  parms.samp[,31]=v.egfr2L_D.p1
  parms.samp[,32]=v.egfr2L_D.p2
  
  
  
  ############################ ALK 2L ############################
  
  temp.cov.alkD=read.csv(paste0(getwd(),"/NSCLC microsim data/cov ALK-2L_D.csv"))
  
  log.alk2L_D.p1=log(parms[33,2])
  log.alk2L_D.p2=log(parms[34,2])
  
  op.par=c(log.alk2L_D.p1,log.alk2L_D.p2)
  
  mvnorm.res <- mvrnorm(n=n_sim, mu=op.par, Sigma=temp.cov.alkD)
  
  v.alk2L_D.p1=exp(mvnorm.res[,1])
  v.alk2L_D.p2=exp(mvnorm.res[,2])
  
  parms.samp[,33]=v.alk2L_D.p1
  parms.samp[,34]=v.alk2L_D.p2
  
  
  return(parms.samp)
}



