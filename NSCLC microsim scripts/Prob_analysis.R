library(ggplot2)
library(scales)
library(patchwork)

###################### Cost Effectiveness Plots ######################

WTP.thresh=read.csv(paste0(getwd(),"/Outputs/PSA flaura2 EGFR WTP.csv"))
WTP.thresh2=read.csv(paste0(getwd(),"/Outputs/PSA MARIPOSA EGFR WTP.csv"))

temp=WTP.thresh[WTP.thresh$Treatment=="FLAURA2",]
temp2=WTP.thresh2[WTP.thresh2$Treatment=="MARIPOSA",]

WTP.bind=rbind(temp,temp2)

WTP.bind

ceac.plot=ggplot(WTP.bind, aes(x=WTP, y=prob,group=Treatment))+
  geom_line(aes(color=Treatment))+
  geom_vline(xintercept = 50000,linetype="dashed")+ 
  geom_vline(xintercept = 100000,linetype="dashed")+
  scale_x_continuous(name="Willingness-to-Pay",labels=dollar)+
  scale_y_continuous(name="Probability of Cost-Effectiveness",labels=percent)+
  annotate("text", x = 50000, y = 0.9, label = "$50k/QALY", hjust = -0.05, vjust = -0.5,size=3)+
  annotate("text", x = 100000, y = 1, label = "$100k/QALY", hjust = -0.05, vjust = -0.5,size=3)+
  scale_color_manual(
    values = c(#F8766D
      "FLAURA2" = "#F8766D",
      "MARIPOSA" = "#619CFF"))+theme(aspect.ratio = 1)+theme_bw()
ceac.plot

PSAresult=read.csv(paste0(getwd(),"/Outputs/PSA flaura2 EGFR inc and ICERs.csv"))
PSAresult2=read.csv(paste0(getwd(),"/Outputs/PSA MARIPOSA EGFR inc and ICERs.csv"))

WTP50=data.frame("x"=c(-0.5,2),"y"=c(-25000,100000))
WTP100=data.frame("x"=c(-0.5,2),"y"=c(-50000,200000))

ceap.plot=ggplot()+
  geom_point(data = PSAresult2, aes(x = incQALY, y = incCosts, colour = "MARIPOSA"), size = 0.7, shape = 16, alpha=0.8)+
  stat_ellipse(data = PSAresult2, aes(x = incQALY, y = incCosts, colour = "MARIPOSA"),level = 0.95,type = "t",linewidth = 0.6)+
  geom_point(data = PSAresult,aes(x = incQALY, y = incCosts, colour = "FLAURA2"),size = 0.7,shape = 16, alpha=0.8)+
  stat_ellipse(data = PSAresult,aes(x = incQALY, y = incCosts, colour = "FLAURA2"),level = 0.95,type = "t",linewidth = 0.6)+
  geom_line(data=WTP50, aes(x = x, y = y),color = "blue", linetype = "dashed", linewidth = 1)+
  geom_line(data=WTP100, aes(x = x, y = y),color = "red", linetype = "dashed", linewidth = 1)+
  geom_segment(aes(x = 0, y = -50000, xend = 0, yend = 500000))+
  geom_segment(aes(x = -0.5, y = 0, xend = 2, yend = 0))+
  geom_segment(data = data.frame(y = seq(-50000,500000, by=50000 )),aes(x = -0.05, xend = 0.05, y = y, yend = y))+
  geom_segment(data = data.frame(x = seq(-0.5, 2, by = 0.5)),aes(x = x, xend = x, y = -10000, yend = 10000))+
  annotate("text", x = 1.2, y = 75000,label = "WTP $50,000/QALY", angle = 13, hjust = 0)+
  annotate("text", x = 1.1, y = 127000,label = "WTP $100,000/QALY", angle = 23, hjust = 0)+
  scale_x_continuous(name = "Incremental QALYs")+
  scale_y_continuous(name = "Incremental Costs (2025 $CAD)",labels = scales::dollar_format())+
  scale_colour_manual(name = "PSA Results",values = c("MARIPOSA" = "#619CFF","FLAURA2" = "#F8766D"),guide = "none")+
  theme(aspect.ratio = 1)+theme_classic()
ceap.plot

comb_case=(ceac.plot | ceap.plot) + plot_layout(guides="collect") & theme(legend.position="bottom")
comb_case

ggsave("C:/Users/Owner/Documents/case study combined.pdf", plot = comb_case, width = 10, height = 5, units = "in", dpi = 400)
ggsave("C:/Users/Owner/Documents/case study combined.jpg", plot = comb_case, width = 10, height = 5, units = "in", dpi = 400)


c.egfr.comp=read.csv(paste0(getwd(),"/Outputs/PSA flaura2 c_egfr_comp.csv"))
c.egfr.comp2=read.csv(paste0(getwd(),"/Outputs/PSA MARIPOSA c_egfr_comp2.csv"))

u.egfr.comp=read.csv(paste0(getwd(),"/Outputs/PSA flaura2 u_egfr_comp.csv"))
u.egfr.comp2=read.csv(paste0(getwd(),"/Outputs/PSA MARIPOSA u_egfr_comp2.csv"))


incCosts.tmp=(c.egfr.comp2-c.egfr.comp)
incQALY.tmp=(u.egfr.comp2-u.egfr.comp)
NMB.tmp=(incQALY.tmp*50000)-incCosts.tmp


sum(NMB.tmp <0)

temp[501,]
temp[1001,]
temp2[501,]
temp2[1001,]

mean(PSAresult$incCosts)
t.test(PSAresult$incCosts)$conf.int
mean(PSAresult$incQALY)
t.test(PSAresult$incQALY)$conf.int
mean(PSAresult$NMB)
t.test(PSAresult$NMB)$conf.int

mean(PSAresult$incCosts)/mean(PSAresult$incQALY)


mean(PSAresult2$incCosts)
t.test(PSAresult2$incCosts)$conf.int
mean(PSAresult2$incQALY)
t.test(PSAresult2$incQALY)$conf.int
mean(PSAresult2$NMB)
t.test(PSAresult2$NMB)$conf.int


mean(PSAresult2$incCosts)/mean(PSAresult2$incQALY)






