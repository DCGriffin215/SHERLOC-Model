


Probs.determ=function(low,high,alk,egfr,parms){
  library('flexsurv')
  library('survival')
  
  ############################ Stage I ############################
  
  p.tmp=pexp(seq(30, 9030, by=30), rate=parms[1,2], lower.tail=F)
  p.tmp.1=pexp(seq(0, 9000, by=30), rate=parms[1,2], lower.tail=F)
  
  death.p.s1=1-(p.tmp/p.tmp.1)
  
  p.tmp=plnorm(seq(30, 9030, by=30), meanlog=parms[2,2], sdlog=parms[3,2], lower.tail=F)
  p.tmp.1=plnorm(seq(0, 9000, by=30), meanlog=parms[2,2], sdlog=parms[3,2], lower.tail=F)
  
  rec.p.s1=1-(p.tmp/p.tmp.1)
  
  alk.p.s1=rec.p.s1*alk
  low.p.s1=rec.p.s1*low
  high.p.s1=rec.p.s1*high
  egfr.p.s1=rec.p.s1*egfr
  
  stay.p.s1=1-(death.p.s1+rec.p.s1)
  
  sp.1=data.frame(stay.p.s1,death.p.s1,alk.p.s1,low.p.s1,high.p.s1,egfr.p.s1)
  
  
  ############################ Stage II ############################
  
  p.tmp=pexp(seq(30, 9030, by=30), rate=parms[4,2], lower.tail=F)
  p.tmp.1=pexp(seq(0, 9000, by=30), rate=parms[4,2], lower.tail=F)
  
  death.p.s2=1-(p.tmp/p.tmp.1)
  
  p.tmp=pgengamma(seq(30, 9030, by=30), mu=parms[5,2], sigma=parms[6,2], Q=parms[7,2], lower.tail=F)
  p.tmp.1=pgengamma(seq(0, 9000, by=30), mu=parms[5,2], sigma=parms[6,2], Q=parms[7,2], lower.tail=F)
  
  rec.p.s2=1-(p.tmp/p.tmp.1)
  
  alk.p.s2=rec.p.s2*alk
  low.p.s2=rec.p.s2*low
  high.p.s2=rec.p.s2*high
  egfr.p.s2=rec.p.s2*egfr
  
  stay.p.s2=1-(death.p.s2+rec.p.s2)
  
  sp.2=data.frame(stay.p.s2,death.p.s2,alk.p.s2,low.p.s2,high.p.s2,egfr.p.s2)
  
  
  ############################ Stage III ############################
  
  p.tmp=plnorm(seq(30, 9030, by=30), meanlog=parms[8,2], sdlog=parms[9,2], lower.tail=F)
  p.tmp.1=plnorm(seq(0, 9000, by=30), meanlog=parms[8,2], sdlog=parms[9,2], lower.tail=F)
  
  death.p.s3=1-(p.tmp/p.tmp.1)
  
  p.tmp=pgengamma(seq(30, 9030, by=30), mu=parms[10,2], sigma=parms[11,2], Q=parms[12,2], lower.tail=F)
  p.tmp.1=pgengamma(seq(0, 9000, by=30), mu=parms[10,2], sigma=parms[11,2], Q=parms[12,2], lower.tail=F)
  
  rec.p.s3=1-(p.tmp/p.tmp.1)
  
  alk.p.s3=rec.p.s3*alk
  low.p.s3=rec.p.s3*low
  high.p.s3=rec.p.s3*high
  egfr.p.s3=rec.p.s3*egfr
  
  stay.p.s3=1-(death.p.s3+rec.p.s3)
  
  sp.3=data.frame(stay.p.s3,death.p.s3,alk.p.s3,low.p.s3,high.p.s3,egfr.p.s3)
  
  
  ############################ EGFR 1L ############################
  
  p.tmp=pexp(seq(30, 9030, by=30), rate=parms[13,2], lower.tail=F)
  p.tmp.1=pexp(seq(0, 9000, by=30), rate=parms[13,2], lower.tail=F)
  
  death.p.egfr1L=1-(p.tmp/p.tmp.1)
  
  p.tmp=pexp(seq(30, 9030, by=30), rate=parms[14,2], lower.tail=F)
  p.tmp.1=pexp(seq(0, 9000, by=30), rate=parms[14,2], lower.tail=F)
  
  prog.p.egfr1L=1-(p.tmp/p.tmp.1)
  
  stay.p.egfr1L=1-(death.p.egfr1L+prog.p.egfr1L)
  
  sp.egfr.1L=data.frame(stay.p.egfr1L, prog.p.egfr1L, death.p.egfr1L)
  
  
  ############################ ALK 1L ############################
  
  p.tmp=pexp(seq(30, 9030, by=30), rate=parms[15,2], lower.tail=F)
  p.tmp.1=pexp(seq(0, 9000, by=30), rate=parms[15,2], lower.tail=F)
  
  death.p.alk1L=1-(p.tmp/p.tmp.1)
  
  p.tmp=plnorm(seq(30, 9030, by=30), meanlog=parms[16,2], sdlog=parms[17,2], lower.tail=F)
  p.tmp.1=plnorm(seq(0, 9000, by=30), meanlog=parms[16,2], sdlog=parms[17,2], lower.tail=F)
  
  prog.p.alk1L=1-(p.tmp/p.tmp.1)
  
  stay.p.alk1L=1-(death.p.alk1L+prog.p.alk1L)
  
  sp.alk.1L=data.frame(stay.p.alk1L, prog.p.alk1L, death.p.alk1L)
  
  
  ############################ PDL1 High 1L ############################
  
  p.tmp=pgengamma(seq(30, 9030, by=30), mu=parms[18,2], sigma=parms[19,2], Q=parms[20,2], lower.tail=F)
  p.tmp.1=pgengamma(seq(0, 9000, by=30), mu=parms[18,2], sigma=parms[19,2], Q=parms[20,2], lower.tail=F)
  
  death.p.high1L=1-(p.tmp/p.tmp.1)
  
  p.tmp=plnorm(seq(30, 9030, by=30), meanlog=parms[21,2], sdlog=parms[22,2], lower.tail=F)
  p.tmp.1=plnorm(seq(0, 9000, by=30), meanlog=parms[21,2], sdlog=parms[22,2], lower.tail=F)
  
  prog.p.high1L=1-(p.tmp/p.tmp.1)
  
  stay.p.high1L=1-(death.p.high1L+prog.p.high1L)
  
  sp.high.1L=data.frame(stay.p.high1L, prog.p.high1L, death.p.high1L)
  
  
  ############################ PDL1 low 1L ############################
  
  p.tmp=plnorm(seq(30, 9030, by=30), meanlog=parms[23,2], sdlog=parms[24,2], lower.tail=F)
  p.tmp.1=plnorm(seq(0, 9000, by=30), meanlog=parms[23,2], sdlog=parms[24,2], lower.tail=F)
  
  death.p.low1L=1-(p.tmp/p.tmp.1)
  
  p.tmp=plnorm(seq(30, 9030, by=30), meanlog=parms[25,2], sdlog=parms[26,2], lower.tail=F)
  p.res.1=plnorm(seq(0, 9000, by=30), meanlog=parms[25,2], sdlog=parms[26,2], lower.tail=F)
  
  prog.p.low1L=1-(p.tmp/p.res.1)
  
  stay.p.low1L=1-(death.p.low1L+prog.p.low1L)
  
  sp.low.1L=data.frame(stay.p.low1L, prog.p.low1L, death.p.low1L)
  
  
  ############################ PDL1 Low 2L ############################
  
  p.tmp=plnorm(seq(30, 9030, by=30), meanlog=parms[27,2], sdlog=parms[28,2], lower.tail=F)
  p.tmp.1=plnorm(seq(0, 9000, by=30), meanlog=parms[27,2], sdlog=parms[28,2], lower.tail=F)
  
  death.PDL1low2L=1-(p.tmp/p.tmp.1)
  
  stay.p.PDL1low2=1-(death.PDL1low2L)
  
  sp.low.2L=data.frame(stay.p.PDL1low2,death.PDL1low2L)
  
  
  ############################ PDL1 high 2L ############################
  
  p.tmp=pllogis(seq(30, 9030, by=30), shape=parms[29,2], scale=parms[30,2], lower.tail=F)
  p.tmp.1=pllogis(seq(0, 9000, by=30), shape=parms[29,2], scale=parms[30,2], lower.tail=F)
  
  death.PDL1high2L=1-(p.tmp/p.tmp.1)
  
  stay.p.PDL1high2L=1-(death.PDL1high2L)
  
  sp.high.2L=data.frame(stay.p.PDL1high2L,death.PDL1high2L)
  
  
  ############################ EGFR 2L ############################
  
  p.tmp=pllogis(seq(30, 9030, by=30), shape=parms[31,2], scale=parms[32,2], lower.tail=F)
  p.tmp.1=pllogis(seq(0, 9000, by=30), shape=parms[31,2], scale=parms[32,2], lower.tail=F)
  
  death.egfr2L=1-(p.tmp/p.tmp.1)
  
  stay.p.egfr2L=1-(death.egfr2L)
  
  sp.egfr.2L=data.frame(stay.p.egfr2L,death.egfr2L)
  
  
  ############################ ALK 2L ############################
  
  p.tmp=pgamma(seq(30, 9030, by=30), shape=parms[33,2], rate=parms[34,2], lower.tail=F)
  p.tmp.1=pgamma(seq(0, 9000, by=30), shape=parms[33,2], rate=parms[34,2], lower.tail=F)
  
  death.alk2L=1-(p.tmp/p.tmp.1)
  
  stay.p.alk2L=1-(death.alk2L)
  
  sp.alk.2L=data.frame(stay.p.alk2L,death.alk2L)
  
  
  ############################ FLAURA2 ############################
  
  death_FL2=1-exp(log(1-death.p.egfr1L)* 0.77)
  
  prog_FL2=1-exp(log(1-prog.p.egfr1L)* 0.73)
  
  stay.p.FL2=1-(death_FL2+prog_FL2)
  
  sp.FL2_HR.1L=data.frame(stay.p.FL2, prog_FL2, death_FL2)
  
  
  ############################ MARIPOSA ############################
  
  death_MARIP=1-exp(log(1-death.p.egfr1L)* 0.75)
  
  prog_MARIP=1-exp(log(1-prog.p.egfr1L)* 0.76)
  
  stay.p.MARIP=1-(death_MARIP+prog_MARIP)
  
  sp.MARIP_HR.1L=data.frame(stay.p.MARIP, prog_MARIP, death_MARIP)
  
  return(list(
    sp.1=sp.1,
    sp.2=sp.2,
    sp.3=sp.3,
    sp.egfr.1L=sp.egfr.1L,
    sp.alk.1L=sp.alk.1L,
    sp.high.1L=sp.high.1L,
    sp.low.1L=sp.low.1L,
    sp.low.2L=sp.low.2L,
    sp.high.2L=sp.high.2L,
    sp.egfr.2L=sp.egfr.2L,
    sp.alk.2L=sp.alk.2L,
    sp.MARIP_HR.1L=sp.MARIP_HR.1L,
    sp.FL2_HR.1L=sp.FL2_HR.1L
  ))
}







