Probs.determ.het.early=function(parms,ind.df){
  library('flexsurv')
  library('survival')
  
  ############################ Stage I ############################
  
  p.tmp=pexp(seq(30, 9030, by=30), rate=parms[1,2], lower.tail=F)
  p.tmp.1=pexp(seq(0, 9000, by=30), rate=parms[1,2], lower.tail=F)
  
  death.p.s1=1-(p.tmp/p.tmp.1)
  
  
  pred.raw=parms[2,2]+ (parms[4,2]*ind.df$age)    #age
  pred.rec=pred.raw              
  
  p.tmp=plnorm(seq(30, 9030, by=30), meanlog=pred.rec, sdlog=parms[3,2], lower.tail=F)
  p.tmp.1=plnorm(seq(0, 9000, by=30), meanlog=pred.rec, sdlog=parms[3,2], lower.tail=F)
  
  rec.p.s1=1-(p.tmp/p.tmp.1)
  
  
  stay.p.s1=1-(death.p.s1+rec.p.s1)
  
  sp.1=data.frame(stay.p.s1,death.p.s1,rec.p.s1)
  
  
  ############################ Stage II ############################
  
  p.tmp=pexp(seq(30, 9030, by=30), rate=parms[5,2], lower.tail=F)
  p.tmp.1=pexp(seq(0, 9000, by=30), rate=parms[5,2], lower.tail=F)
  
  death.p.s2=1-(p.tmp/p.tmp.1)
  
  
  pred.raw=parms[6,2]+(parms[9,2]*ind.df$sex+parms[10,2]*ind.df$age)  #sex+age 
  pred.rec=pred.raw
  
  p.tmp=pgengamma(seq(30, 9030, by=30), mu=pred.rec, sigma=parms[7,2], Q=parms[8,2], lower.tail=F)
  p.tmp.1=pgengamma(seq(0, 9000, by=30), mu=pred.rec, sigma=parms[7,2], Q=parms[8,2], lower.tail=F)
  
  rec.p.s2=1-(p.tmp/p.tmp.1)
  
  
  stay.p.s2=1-(death.p.s2+rec.p.s2)
  
  sp.2=data.frame(stay.p.s2,death.p.s2,rec.p.s2)
  
  
  ############################ Stage III ############################
  
  p.tmp=plnorm(seq(30, 9030, by=30), meanlog=parms[11,2], sdlog=parms[12,2], lower.tail=F)
  p.tmp.1=plnorm(seq(0, 9000, by=30), meanlog=parms[11,2], sdlog=parms[12,2], lower.tail=F)
  
  death.p.s3=1-(p.tmp/p.tmp.1)
  
  p.tmp=pgengamma(seq(30, 9030, by=30), mu=parms[13,2], sigma=parms[14,2], Q=parms[15,2], lower.tail=F)
  p.tmp.1=pgengamma(seq(0, 9000, by=30), mu=parms[13,2], sigma=parms[14,2], Q=parms[15,2], lower.tail=F)
  
  rec.p.s3=1-(p.tmp/p.tmp.1)
  
  
  stay.p.s3=1-(death.p.s3+rec.p.s3)
  
  sp.3=data.frame(stay.p.s3,death.p.s3,rec.p.s3)
  
  return(list(
    sp.1=sp.1,
    sp.2=sp.2,
    sp.3=sp.3
  ))
}


