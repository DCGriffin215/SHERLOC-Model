Probs.determ.het.late=function(parms,ind.df, hr.prob.FLAURA2, hr.prob.MARIP){
  library('flexsurv')
  library('survival')
  
  
  ############################ EGFR 1L ############################
  
  pred.raw=log(parms[16,2])+(parms[17,2]*ind.df$ecog)  #ecog 
  pred.death=exp(pred.raw)
  
  p.tmp=pllogis(seq(30, 9030, by=30), shape=parms[15,2], scale=pred.death, lower.tail=F)
  p.tmp.1=pllogis(seq(0, 9000, by=30), shape=parms[15,2], scale=pred.death, lower.tail=F)
  
  death.p.egfr1L=1-(p.tmp/p.tmp.1)
  
  
  pred.raw=log(parms[18,2])+(parms[19,2]*ind.df$ecog)  #ecog 
  pred.prog=exp(pred.raw)
  
  p.tmp=pexp(seq(30, 9030, by=30), rate=pred.prog, lower.tail=F)
  p.tmp.1=pexp(seq(0, 9000, by=30), rate=pred.prog, lower.tail=F)
  
  prog.p.egfr1L=1-(p.tmp/p.tmp.1)
  
  stay.p.egfr1L=1-(death.p.egfr1L+prog.p.egfr1L)
  
  sp.egfr.1L=data.frame(stay.p.egfr1L, prog.p.egfr1L, death.p.egfr1L)
  
  
  ############################ ALK 1L ############################
  
  p.tmp=pexp(seq(30, 9030, by=30), rate=parms[20,2], lower.tail=F)
  p.tmp.1=pexp(seq(0, 9000, by=30), rate=parms[20,2], lower.tail=F)
  
  death.p.alk1L=1-(p.tmp/p.tmp.1)
  
  
  p.tmp=pexp(seq(30, 9030, by=30), rate=parms[21,2], lower.tail=F)
  p.tmp.1=pexp(seq(0, 9000, by=30), rate=parms[21,2], lower.tail=F)
  
  prog.p.alk1L=1-(p.tmp/p.tmp.1)
  
  stay.p.alk1L=1-(death.p.alk1L+prog.p.alk1L)
  
  sp.alk.1L=data.frame(stay.p.alk1L, prog.p.alk1L, death.p.alk1L)
  
  
  ############################ PDL1 High 1L ############################
  
  pred.raw=parms[22,2]+ (parms[24,2]*ind.df$sex+parms[25,2]*ind.df$age+parms[26,2]*ind.df$ecog)    #sex+age+ecog
  pred.death=pred.raw  
  
  p.tmp=plnorm(seq(30, 9030, by=30), meanlog=pred.death, sdlog=parms[23,2], lower.tail=F)
  p.tmp.1=plnorm(seq(0, 9000, by=30), meanlog=pred.death, sdlog=parms[23,2], lower.tail=F)
  
  death.p.high1L=1-(p.tmp/p.tmp.1)
  
  pred.raw=parms[27,2]+ (parms[29,2]*ind.df$sex+parms[30,2]*ind.df$ecog)    #sex+ecog
  pred.prog=pred.raw 
  
  p.tmp=plnorm(seq(30, 9030, by=30), meanlog=pred.prog, sdlog=parms[28,2], lower.tail=F)
  p.tmp.1=plnorm(seq(0, 9000, by=30), meanlog=pred.prog, sdlog=parms[28,2], lower.tail=F)
  
  prog.p.high1L=1-(p.tmp/p.tmp.1)
  
  stay.p.high1L=1-(death.p.high1L+prog.p.high1L)
  
  sp.high.1L=data.frame(stay.p.high1L, prog.p.high1L, death.p.high1L)
  
  
  ############################ PDL1 low 1L ############################
  
  pred.raw=parms[31,2]+ (parms[33,2]*ind.df$sex+parms[34,2]*ind.df$age+parms[35,2]*ind.df$ecog)    #sex+age+ecog
  pred.death=pred.raw  
  
  p.tmp=plnorm(seq(30, 9030, by=30), meanlog=pred.death, sdlog=parms[32,2], lower.tail=F)
  p.tmp.1=plnorm(seq(0, 9000, by=30), meanlog=pred.death, sdlog=parms[32,2], lower.tail=F)
  
  death.p.low1L=1-(p.tmp/p.tmp.1)
  
  
  pred.raw=parms[36,2]+ (parms[38,2]*ind.df$age)    #age
  pred.prog=pred.raw 
  
  p.tmp=plnorm(seq(30, 9030, by=30), meanlog=pred.prog, sdlog=parms[37,2], lower.tail=F)
  p.res.1=plnorm(seq(0, 9000, by=30), meanlog=pred.prog, sdlog=parms[37,2], lower.tail=F)
  
  prog.p.low1L=1-(p.tmp/p.res.1)
  
  stay.p.low1L=1-(death.p.low1L+prog.p.low1L)
  
  sp.low.1L=data.frame(stay.p.low1L, prog.p.low1L, death.p.low1L)
  
  
  ############################ PDL1 Low 2L ############################
  
  pred.raw=parms[39,2]+ (parms[41,2]*ind.df$age+parms[42,2]*ind.df$ecog)    #age+ecog
  pred.death=pred.raw 
  
  p.tmp=plnorm(seq(30, 9030, by=30), meanlog=pred.death, sdlog=parms[40,2], lower.tail=F)
  p.tmp.1=plnorm(seq(0, 9000, by=30), meanlog=pred.death, sdlog=parms[40,2], lower.tail=F)
  
  death.PDL1low2L=1-(p.tmp/p.tmp.1)
  
  stay.p.PDL1low2=1-(death.PDL1low2L)
  
  sp.low.2L=data.frame(stay.p.PDL1low2,death.PDL1low2L)
  
  
  ############################ PDL1 high 2L ############################
  
  p.tmp=pllogis(seq(30, 9030, by=30), shape=parms[43,2], scale=parms[44,2], lower.tail=F)
  p.tmp.1=pllogis(seq(0, 9000, by=30), shape=parms[43,2], scale=parms[44,2], lower.tail=F)
  
  death.PDL1high2L=1-(p.tmp/p.tmp.1)
  
  stay.p.PDL1high2L=1-(death.PDL1high2L)
  
  sp.high.2L=data.frame(stay.p.PDL1high2L,death.PDL1high2L)
  
  
  ############################ EGFR 2L ############################
  
  p.tmp=pexp(seq(30, 9030, by=30), rate=parms[45,2], lower.tail=F)
  p.tmp.1=pexp(seq(0, 9000, by=30), rate=parms[45,2], lower.tail=F)
  
  death.egfr2L=1-(p.tmp/p.tmp.1)
  
  stay.p.egfr2L=1-(death.egfr2L)
  
  sp.egfr.2L=data.frame(stay.p.egfr2L,death.egfr2L)
  
  
  ############################ ALK 2L ############################
  
  p.tmp=pexp(seq(30, 9030, by=30),rate=parms[46,2], lower.tail=F)
  p.tmp.1=pexp(seq(0, 9000, by=30),rate=parms[46,2], lower.tail=F)
  
  death.alk2L=1-(p.tmp/p.tmp.1)
  
  stay.p.alk2L=1-(death.alk2L)
  
  sp.alk.2L=data.frame(stay.p.alk2L,death.alk2L)
  
  
  ############################ FLAURA2 ############################
  
  death_FL2=1-exp(log(1-death.p.egfr1L)* hr.prob.FLAURA2[,1])
  
  prog_FL2=1-exp(log(1-prog.p.egfr1L)* hr.prob.FLAURA2[,2])
  
  stay.p.FL2=1-(death_FL2+prog_FL2)
  
  sp.FL2_HR.1L=data.frame(stay.p.FL2, prog_FL2, death_FL2)
  
  
  ############################ MARIPOSA ############################
  
  death_MARIP=1-exp(log(1-death.p.egfr1L)* hr.prob.MARIP[,1])
  
  prog_MARIP=1-exp(log(1-prog.p.egfr1L)* hr.prob.MARIP[,2])
  
  stay.p.MARIP=1-(death_MARIP+prog_MARIP)
  
  sp.MARIP_HR.1L=data.frame(stay.p.MARIP, prog_MARIP, death_MARIP)
  
  return(list(
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