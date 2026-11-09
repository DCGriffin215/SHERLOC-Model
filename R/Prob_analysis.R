library(ggplot2)
library(scales)
library(patchwork)
library(ggprism)

###################### Cost Effectiveness Plots ######################

WTP.thresh=read.csv(paste0(getwd(),"/Outputs/PSA flaura2 EGFR WTP alt.csv"))
WTP.thresh2=read.csv(paste0(getwd(),"/Outputs/PSA MARIPOSA EGFR WTP alt.csv"))

temp=WTP.thresh[WTP.thresh$Treatment=="FLAURA2",]
temp2=WTP.thresh2[WTP.thresh2$Treatment=="MARIPOSA",]

WTP.bind=rbind(temp,temp2)

ceac.plot=ggplot(WTP.bind, aes(x=WTP, y=prob,group=Treatment))+
  geom_line(aes(color=Treatment),linewidth = 0.3)+
  geom_vline(xintercept = 50000,linetype="dashed",linewidth = 0.2)+ 
  geom_vline(xintercept = 100000,linetype="dashed",linewidth = 0.2)+
  scale_x_continuous(name="Willingness-to-Pay",labels=dollar)+
  scale_y_continuous(name="Probability of Cost-Effectiveness",labels=percent)+
  annotate("text", x = 50000, y = 0.9, label = "$50k/QALY", hjust = -0.05, vjust = -0.5,size=1.5)+
  annotate("text", x = 100000, y = 1, label = "$100k/QALY", hjust = -0.05, vjust = -0.5,size=1.5)+
  annotate("text", x=0.25, y=1,label="A",size=2, fontface = "bold")+
  scale_color_prism(labels=c("FLAURA2","MARIPOSA"),palette = "winter_bright")+theme_prism()+
  theme(aspect.ratio = 1,
        axis.title.x= element_text(size = 6,margin = margin(t = 3, 0, 0, 0),unit='in'),
        axis.title.y = element_text(size = 6,margin = margin(r = 0.5, 0, 0, 0),unit='in'),
        axis.text = element_text(size = 4,margin = margin(0, 0, 0, 0)),
        legend.text = element_text(size = 4),
        legend.key.size = unit(0.05, "in"),
        legend.spacing.x = unit(0.05, "in"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0),plot.margin = margin(0, 0, 0, 0, unit = "in"),
        axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.5),
        axis.ticks.length = unit(0.05, "in"))
ceac.plot

PSAresult=read.csv(paste0(getwd(),"/Outputs/PSA flaura2 EGFR inc and ICERs alt.csv"))
PSAresult2=read.csv(paste0(getwd(),"/Outputs/PSA MARIPOSA EGFR inc and ICERs alt.csv"))

WTP50=data.frame("x"=c(-0.5,2),"y"=c(-25000,100000))
WTP100=data.frame("x"=c(-0.5,2),"y"=c(-50000,200000))

ceap.plot=ggplot()+
  geom_point(data = PSAresult2, aes(x = incQALY, y = incCosts, colour = "MARIPOSA"), size = 0.15, shape = 16, alpha=0.8)+
  stat_ellipse(data = PSAresult2, aes(x = incQALY, y = incCosts, colour = "MARIPOSA"),level = 0.95,type = "t",linewidth = 0.3)+
  geom_point(data = PSAresult,aes(x = incQALY, y = incCosts, colour = "FLAURA2"),size = 0.15,shape = 16, alpha=0.8)+
  stat_ellipse(data = PSAresult,aes(x = incQALY, y = incCosts, colour = "FLAURA2"),level = 0.95,type = "t",linewidth = 0.3)+
  geom_line(data=WTP50, aes(x = x, y = y),color = "#0072B2", linetype = "dashed", linewidth = 0.5)+
  geom_line(data=WTP100, aes(x = x, y = y),color = "#F8766D", linetype = "dashed", linewidth = 0.5)+
  geom_segment(aes(x = 0, y = -50000, xend = 0, yend = 500000),linewidth = 0.2)+
  geom_segment(aes(x = -0.5, y = 0, xend = 2, yend = 0),linewidth = 0.2)+
  geom_segment(data = data.frame(y = seq(-50000,500000, by=50000 )),aes(x = -0.05, xend = 0.05, y = y, yend = y),linewidth = 0.2)+
  geom_segment(data = data.frame(x = seq(-0.5, 2, by = 0.5)),aes(x = x, xend = x, y = -10000, yend = 10000),linewidth = 0.2)+
  annotate("text", x = 1.2, y = 75000,label = "WTP $50,000/QALY", angle = 13, hjust = 0,size=1.5)+
  annotate("text", x = 1.1, y = 127000,label = "WTP $100,000/QALY", angle = 23, hjust = 0,size=1.5)+
  annotate("text", x=-0.45, y=500000,label="B",size=2, fontface = "bold")+
  scale_x_continuous(name = "Incremental QALYs")+
  scale_y_continuous(name = "Incremental Costs (2025 $CAD)",labels = scales::dollar_format())+
  scale_colour_prism(name = "PSA Results",labels=c("FLAURA2","MARIPOSA"),palette = "winter_bright",guide = "none")+theme_prism()+
  theme(aspect.ratio = 1,
        axis.title.x= element_text(size = 6,margin = margin(t = 3, 0, 0, 0),unit='in'),
        axis.title.y = element_text(size = 6,margin = margin(r = 0.5, 0, 0, 0),unit='in'),
        axis.text = element_text(size = 4,margin = margin(0, 0, 0, 0)),
        legend.text = element_text(size = 4),
        legend.key.size = unit(0.05, "in"),
        legend.spacing.x = unit(0.05, "in"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0),plot.margin = margin(0, 0, 0, 0, unit = "in"),
        axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.5),
        axis.ticks.length = unit(0.05, "in"))
ceap.plot

comb_case=(ceac.plot | ceap.plot) + plot_layout(guides="collect") & theme(legend.position="bottom")
comb_case
paste0(getwd(),"/Outputs/case study combined alt.jpg")
#ggsave(paste0(getwd(),"/Outputs/case study combined alt.jpg"), plot = comb_case, width = 5, height = 2.5, units = "in", dpi = 400)


c.egfr.comp=read.csv(paste0(getwd(),"/Outputs/PSA flaura2 c_egfr_comp alt.csv"))
c.egfr.comp2=read.csv(paste0(getwd(),"/Outputs/PSA MARIPOSA c_egfr_comp2 alt.csv"))

u.egfr.comp=read.csv(paste0(getwd(),"/Outputs/PSA flaura2 u_egfr_comp alt.csv"))
u.egfr.comp2=read.csv(paste0(getwd(),"/Outputs/PSA MARIPOSA u_egfr_comp2 alt.csv"))


incCosts.tmp=(c.egfr.comp2-c.egfr.comp)
incQALY.tmp=(u.egfr.comp2-u.egfr.comp)
NMB.tmp=(incQALY.tmp*50000)-incCosts.tmp


sum(PSAresult$NMB>0)/500

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


sum(PSAresult2$NMB <PSAresult$NMB)/500

sum(((PSAresult2$incQALY*100000)-PSAresult2$incCosts)<((PSAresult$incQALY*100000)-PSAresult$incCosts))/500

incCosts.tmp=(c.egfr.comp2-c.egfr.comp)
incQALY.tmp=(u.egfr.comp2-u.egfr.comp)
NMB.tmp=(incQALY.tmp*50000)-incCosts.tmp
NMB.tmp.100=(incQALY.tmp*100000)-incCosts.tmp
ICER=(incCosts.tmp/incQALY.tmp)

count=0
for(j in 1:500){
  temp=ICER[j]<(50000)
  temp2=incQALY.tmp[j]>0
  temp3=ICER[j]>(50000)
  temp4=incQALY.tmp[j]<0
  if (temp==TRUE & temp2==TRUE){count=count+1}
  if (temp3==TRUE & temp4==TRUE){count=count+1}
}
count/500

for(j in 1:500){
  temp=ICER[j]<(100000)
  temp2=incQALY.tmp[j]>0
  temp3=ICER[j]>(100000)
  temp4=incQALY.tmp[j]<0
  if (temp==TRUE & temp2==TRUE){count=count+1}
  if (temp3==TRUE & temp4==TRUE){count=count+1}
}
count/500

sum(NMB.tmp>0)/500
sum(NMB.tmp.100>0)/500

mean(unlist(incCosts.tmp))
t.test(unlist(incCosts.tmp))$conf.int
mean(unlist(incQALY.tmp))
t.test(unlist(incQALY.tmp))$conf.int
mean(unlist(NMB.tmp))
t.test(unlist(NMB.tmp))$conf.int


mean(unlist(incCosts.tmp))/mean(unlist(incQALY.tmp))

