


Probs.PSA=function(low, high, alk, egfr, parms.samp.m, hr.prob.FLAURA2, hr.prob.MARIP, j){
  library('flexsurv')
  library('survival')
  
  ############################ Stage I ############################
  
  p.res   <- pexp(seq(30, 9030, by = 30), rate =parms.samp.m[j,1], lower.tail = F)
  p.res.1 <- pexp(seq(0, 9000, by = 30), rate = parms.samp.m[j,1], lower.tail = F)
  
  death.p.s1=1 - (p.res/p.res.1)
  
  p.res   <- plnorm(seq(30, 9030, by = 30), meanlog = parms.samp.m[j,2], sdlog = parms.samp.m[j,3], lower.tail = F)
  p.res.1 <- plnorm(seq(0, 9000, by = 30), meanlog = parms.samp.m[j,2], sdlog = parms.samp.m[j,3], lower.tail = F)
  
  rec.p.s1=1 - (p.res/p.res.1)
  
  alk.p.s1=rec.p.s1*alk
  low.p.s1=rec.p.s1*low
  high.p.s1=rec.p.s1*high
  egfr.p.s1=rec.p.s1*egfr
  
  stay.p.s1=1-(death.p.s1+rec.p.s1)
  
  sp.1=data.frame(stay.p.s1,death.p.s1,alk.p.s1,low.p.s1,high.p.s1,egfr.p.s1)
  
  
  
  ############################ Stage II ############################
  
  p.res   <- pexp(seq(30, 9030, by = 30), rate = parms.samp.m[j,4], lower.tail = F)
  p.res.1 <- pexp(seq(0, 9000, by = 30), rate = parms.samp.m[j,4], lower.tail = F)
  
  death.p.s2=1 - (p.res/p.res.1)
  
  
  p.res   <- pgengamma(seq(30, 9030, by = 30), mu = parms.samp.m[j,5], sigma = parms.samp.m[j,6], Q=parms.samp.m[j,7], lower.tail = F)
  p.res.1 <- pgengamma(seq(0, 9000, by = 30), mu = parms.samp.m[j,5], sigma = parms.samp.m[j,6], Q=parms.samp.m[j,7], lower.tail = F)
  
  rec.p.s2=1 - (p.res/p.res.1)
  
  alk.p.s2=rec.p.s2*alk
  low.p.s2=rec.p.s2*low
  high.p.s2=rec.p.s2*high
  egfr.p.s2=rec.p.s2*egfr
  
  stay.p.s2=1-(death.p.s2+rec.p.s2)
  
  sp.2=data.frame(stay.p.s2,death.p.s2,alk.p.s2,low.p.s2,high.p.s2,egfr.p.s2)

  
  
  
  ############################ Stage III ############################
  
  p.res   <- plnorm(seq(30, 9030, by = 30), meanlog = parms.samp.m[j,8], sdlog = parms.samp.m[j,9], lower.tail = F)
  p.res.1 <- plnorm(seq(0, 9000, by = 30), meanlog = parms.samp.m[j,8], sdlog = parms.samp.m[j,9], lower.tail = F)
  
  death.p.s3=1 - (p.res/p.res.1)
  
  
  p.res   <- pgengamma(seq(30, 9030, by = 30), mu = parms.samp.m[j,10], sigma = parms.samp.m[j,11], Q=parms.samp.m[j,12], lower.tail = F)
  p.res.1 <- pgengamma(seq(0, 9000, by = 30), mu = parms.samp.m[j,10], sigma = parms.samp.m[j,11], Q=parms.samp.m[j,12], lower.tail = F)
  
  rec.p.s3=1 - (p.res/p.res.1)
  
  alk.p.s3=rec.p.s3*alk
  low.p.s3=rec.p.s3*low
  high.p.s3=rec.p.s3*high
  egfr.p.s3=rec.p.s3*egfr
  
  stay.p.s3=1-(death.p.s3+rec.p.s3)
  
  sp.3=data.frame(stay.p.s3,death.p.s3,alk.p.s3,low.p.s3,high.p.s3,egfr.p.s3)

  
  
  ############################ EGFR 1L ############################
  
  p.res   <- pexp(seq(30, 9030, by = 30), rate = parms.samp.m[j,13], lower.tail = F)
  p.res.1 <- pexp(seq(0, 9000, by = 30), rate = parms.samp.m[j,13], lower.tail = F)
  
  death.p.egfr1L=1 - (p.res/p.res.1)
  
  
  p.res   <- pexp(seq(30, 9030, by = 30), rate = parms.samp.m[j,14], lower.tail = F)
  p.res.1 <- pexp(seq(0, 9000, by = 30), rate = parms.samp.m[j,14], lower.tail = F)
  
  prog.p.egfr1L=1 - (p.res/p.res.1)
  
  
  stay.p.egfr1L=1-(death.p.egfr1L+prog.p.egfr1L)
  
  sp.egfr.1L=data.frame(stay.p.egfr1L, prog.p.egfr1L, death.p.egfr1L)

  
  
  ############################ ALK 1L ############################
  
  p.res   <- pexp(seq(30, 9030, by = 30), rate = parms.samp.m[j,15], lower.tail = F)
  p.res.1 <- pexp(seq(0, 9000, by = 30), rate = parms.samp.m[j,15], lower.tail = F)
  
  death.p.alk1L=1 - (p.res/p.res.1)
  
  p.res   <- plnorm(seq(30, 9030, by = 30), meanlog = parms.samp.m[j,16], sdlog = parms.samp.m[j,17], lower.tail = F)
  p.res.1 <- plnorm(seq(0, 9000, by = 30), meanlog = parms.samp.m[j,16], sdlog = parms.samp.m[j,17], lower.tail = F)
  
  prog.p.alk1L=1 - (p.res/p.res.1)
  
  
  stay.p.alk1L=1-(death.p.alk1L+prog.p.alk1L)
  
  sp.alk.1L=data.frame(stay.p.alk1L, prog.p.alk1L, death.p.alk1L)
  
  
  ############################ PDL1 High 1L ############################
  
  p.res   <- pgengamma(seq(30, 9030, by = 30), mu = parms.samp.m[j,18], sigma = parms.samp.m[j,19], Q=parms.samp.m[j,20], lower.tail = F)
  p.res.1 <- pgengamma(seq(0, 9000, by = 30), mu = parms.samp.m[j,18], sigma = parms.samp.m[j,19], Q=parms.samp.m[j,20], lower.tail = F)
  
  death.p.high1L=1 - (p.res/p.res.1)
  
  p.res   <- plnorm(seq(30, 9030, by = 30), meanlog = parms.samp.m[j,21], sdlog = parms.samp.m[j,22], lower.tail = F)
  p.res.1 <- plnorm(seq(0, 9000, by = 30), meanlog = parms.samp.m[j,21], sdlog = parms.samp.m[j,22], lower.tail = F)
  
  prog.p.high1L=1 - (p.res/p.res.1)
  
  
  stay.p.high1L=1-(death.p.high1L+prog.p.high1L)
  
  sp.high.1L=data.frame(stay.p.high1L, prog.p.high1L, death.p.high1L)
  
  
  
  ############################ PDL1 low 1L ############################
  
  p.res   <- plnorm(seq(30, 9030, by = 30), meanlog = parms.samp.m[j,23], sdlog = parms.samp.m[j,24], lower.tail = F)
  p.res.1 <- plnorm(seq(0, 9000, by = 30), meanlog = parms.samp.m[j,23], sdlog = parms.samp.m[j,24], lower.tail = F)
  
  death.p.low1L=1 - (p.res/p.res.1)
  
  p.res   <- plnorm(seq(30, 9030, by = 30), meanlog = parms.samp.m[j,25], sdlog = parms.samp.m[j,26], lower.tail = F)
  p.res.1 <- plnorm(seq(0, 9000, by = 30), meanlog = parms.samp.m[j,25], sdlog = parms.samp.m[j,26], lower.tail = F)
  
  prog.p.low1L=1 - (p.res/p.res.1)
  
  
  stay.p.low1L=1-(death.p.low1L+prog.p.low1L)
  
  sp.low.1L=data.frame(stay.p.low1L, prog.p.low1L, death.p.low1L)
  
  
  
  ############################ PDL1 Low 2L ############################
  
  p.res   <- plnorm(seq(30, 9030, by = 30), meanlog = parms.samp.m[j,27], sdlog = parms.samp.m[j,28], lower.tail = F)
  p.res.1 <- plnorm(seq(0, 9000, by = 30), meanlog = parms.samp.m[j,27], sdlog = parms.samp.m[j,28], lower.tail = F)
  
  death.PDL1low2L=1-(p.res/p.res.1)
  
  stay.p.PDL1low2=1-(death.PDL1low2L)
  
  sp.low.2L=data.frame(stay.p.PDL1low2,death.PDL1low2L)
  
  
  
  ############################ PDL1 High 2L ############################
  
  p.res   <- pllogis(seq(30, 9030, by = 30), shape = parms.samp.m[j,29], scale = parms.samp.m[j,30], lower.tail = F)
  p.res.1 <- pllogis(seq(0, 9000, by = 30), shape = parms.samp.m[j,29], scale = parms.samp.m[j,30], lower.tail = F)
  
  death.PDL1high2L=1-(p.res/p.res.1)
  
  stay.p.PDL1high2L=1-(death.PDL1high2L)
  
  sp.high.2L=data.frame(stay.p.PDL1high2L,death.PDL1high2L)
  
  
  
  ############################ EGFR 2L ############################
  
  p.res   <- pllogis(seq(30, 9030, by = 30), shape = parms.samp.m[j,31], scale = parms.samp.m[j,32], lower.tail = F)
  p.res.1 <- pllogis(seq(0, 9000, by = 30), shape = parms.samp.m[j,31], scale = parms.samp.m[j,32], lower.tail = F)
  
  death.egfr2L=1-(p.res/p.res.1)
  
  stay.p.egfr2L=1-(death.egfr2L)
  
  sp.egfr.2L=data.frame(stay.p.egfr2L,death.egfr2L)
  
  
  ############################ ALK 2L ############################
  
  p.res   <- pgamma(seq(30, 9030, by = 30), shape = parms.samp.m[j,33], rate = parms.samp.m[j,34], lower.tail = F)
  p.res.1 <- pgamma(seq(0, 9000, by = 30), shape = parms.samp.m[j,33], rate = parms.samp.m[j,34], lower.tail = F)
  
  death.alk2L=1-(p.res/p.res.1)
  
  stay.p.alk2L=1-(death.alk2L)
  
  sp.alk.2L=data.frame(stay.p.alk2L,death.alk2L)
  
  
  
  ############################ FLAURA2 ############################
  
  death_FL2 <- 1 - exp(log(1 - death.p.egfr1L)* hr.prob.FLAURA2[j,1])
  
  prog_FL2 <- 1 - exp(log(1 - prog.p.egfr1L)* hr.prob.FLAURA2[j,2])
  
  stay.p.FL2=1-(death_FL2+prog_FL2)
  
  sp.FL2_HR.1L=data.frame(stay.p.FL2, prog_FL2, death_FL2)
  
  
  
  ############################ MARIPOSA ############################
  
  death_MARIP <- 1 - exp(log(1 - death.p.egfr1L)* hr.prob.MARIP[j,1])
  
  prog_MARIP <- 1 - exp(log(1 - prog.p.egfr1L)* hr.prob.MARIP[j,2])
  
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



