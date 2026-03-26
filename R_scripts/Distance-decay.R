library(ggplot2)

WT = read.delim(file = "wt_intrainteractions.txt",sep="\t",h=F)
WT$Sample = "WT"
WT$distance = WT$V3 - WT$V2
WT$logvalue = log10(WT$V4)
WT$xvalue = WT$distance/5000
WT1 = aggregate(V4 ~ xvalue, data=WT, median)
WT1$logvalue = log10(WT1$V4)
WT1$Sample = "WT"


Top2 = read.delim(file = "top2_intrainteractions.txt",sep="\t",h=F)
Top2$Sample = "Top2-1"
Top2$distance = Top2$V3 - Top2$V2
Top2$logvalue = log10(Top2$V4)
Top2$xvalue = Top2$distance/5000
Top2_1 = aggregate(V4 ~ xvalue, data=Top2, median)
Top2_1$logvalue = log10(Top2_1$V4)
Top2_1$Sample = "Top2-1"


hmo1 = read.delim(file = "hmo1_intrainteractions.txt",sep="\t",h=F)
hmo1$Sample = "Hmo1d"
hmo1$distance = hmo1$V3 - hmo1$V2
hmo1$logvalue = log10(hmo1$V4)
hmo1$xvalue = hmo1$distance/5000
hmo1_1 = aggregate(V4 ~ xvalue, data=hmo1, median)
hmo1_1$logvalue = log10(hmo1_1$V4)
hmo1_1$Sample = "Hmo1d"


####

alldata = rbind(WT1,Top2_1,hmo1_1)

alldata$Sample = factor(alldata$Sample, levels=c("WT","Top2-1","Hmo1d"))

alldata1= alldata[(alldata$xvalue >=0 & alldata$xvalue <= 300),]
dp <- ggplot(alldata1, aes(x=xvalue, y=logvalue, colour=Sample)) + 
  geom_smooth(method = lm, se = FALSE) +
  labs(title="Contact probability with genomic distance",x="Genomic Distance (5Kb)", y = "Log10 Contact probability") + ylim(0,1.6)
ggsave(file = "distance1.jpeg",dp,height=6,width=8,dpi=200)
