library(ggplot2)
library(dplyr)
library(bedr)
options(scipen=999)
library(gtools)
library(plyr)

#### Bins file in ref_files folder

tenkbbins = read.delim("sc3_10kb_bins.txt",sep="\t",h=F)
colnames(tenkbbins) = c("gchr","gstart","gend","bins")
tenkbbins[which(tenkbbins$gstart == 0),"gstart"] = 1

wt_ten = read.delim(file = "wt_inter_telomere.txt",sep="\t",h=F)
telomere = read.delim(file = "telomere.bed",sep="\t",h=F)
wt_ten = wt_ten[,c(4:7,11)]
colnames(wt_ten) = c("chr","start","end","score","telomere")
wt_ten[which(wt_ten$start == 0),"start"] = 1
sc = tenkbbins
alltelo = unique(wt_ten$telomere)
listtelo = list()

bfunc = function(sc,findata,telo){
sc$gid = paste(sc$gchr,":",sc$gstart,"-",sc$gend,sep="")
scv = sc$gid
scv.sort = bedr.sort.region(scv)
findata$pid = paste(findata$chr,":",findata$start,"-",findata$end,sep="")
dtv = findata$pid
dtv.sort = bedr.sort.region(dtv)
dtv.int <- bedr(input = list(a = scv.sort, b = dtv.sort), engine="/home/hpc-344/miniconda3/bin/bedtools", method="intersect",params ="-wao -sorted",tmpDir="/tmp",deleteTmpDir=FALSE)
dtv.int$pid = paste(dtv.int$V4,":",dtv.int$V5,"-",dtv.int$V6,sep="")
colnames(dtv.int)[1] = "gid"
findata1 = merge(findata,dtv.int,by=c("pid"),all=TRUE)
findata1 = merge(findata1,sc,by=c("gid"),all=TRUE)
findata1 = findata1[,c("gchr","gstart","gend","bins","score","telomere")]
findata1[is.na(findata1)] = 0
findata1$telomere = telo
findata1 = findata1[mixedorder(findata1$bins),]
findata1$id = as.numeric(gsub('bid','',findata1$bins))
findata1$score = as.numeric(findata1$score)
return (findata1)
}

for (i in 1:length(alltelo)){
telo = alltelo[i]
findata = wt_ten[which(wt_ten$telomere == telo),]
listtelo[[telo]] = bfunc(sc,findata,telo)
}

df <- ldply(listtelo, data.frame)
head(df)
breaksx = as.numeric(gsub("bid","",tenkbbins[tenkbbins$gstart == 1,"bins"]))

df$telomere = factor(df$telomere, levels = telomere$V4)
head(df)


df_TEL03L <- df[df$telomere == 'TEL03L',]
df_TEL11L <- df[df$telomere == 'TEL11L',]
DF_1 <- rbind(df_TEL03L,df_TEL11L)


p <- ggplot(DF_1,aes(id, score)) + geom_line(aes(colour=telomere), size =0.5) + scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, by = 50)) + scale_x_continuous(breaks=breaksx,labels=paste("chr",c(1:16),sep="")) + scale_colour_manual(values=c(TEL11L="#339999",TEL03L="#CC0033")) + theme_classic() + geom_vline(xintercept = c(breaksx,dim(tenkbbins)[1]), colour = "grey48", linetype = "dotted") + theme(strip.text = element_text(size = 12, colour = "red",face = "bold")) + theme(strip.background = element_blank()) + theme(axis.title.x = element_text( face="bold", colour="black", size=10)) + xlab("TEL03L & TEL11L vs all Telomere Interactions") + theme(axis.title.y = element_text( face="bold", colour="black", size=10))
ggsave(file = "wt_TEL03L-TEL11L.jpeg", p, height = 6, width = 18, dpi = 400)
