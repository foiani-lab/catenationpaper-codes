library("strawr")
library("reshape2")

wt_tel = read.delim("telomere_interactions_50kb.txt",sep="\t",h=F)
tel_armlength = tel_armlength = read.delim(file = "chromosomearm_shortandlong.txt",sep="\t",h=F)
tnames = paste0(tel_armlength$V1,"(",tel_armlength$V2,")")
wt_tel = wt_tel[,c(4,8,9)]
vecmat = aggregate(V9 ~ V4, wt_tel, median)[,2]
wt_tel_mat = as.matrix(acast(wt_tel, V4~V8, value.var="V9"))
wt_tel_mat = round(wt_tel_mat/vecmat,2)*100
wt_tel_mat[upper.tri(wt_tel_mat)] <- -1
wt_tel_mat[is.na(wt_tel_mat)] = -1
wt_tel_mat[wt_tel_mat <  as.numeric(summary(as.numeric(wt_tel_mat))[5])] = -1
rownames(wt_tel_mat) = tnames
colnames(wt_tel_mat) = tnames
d <- dist(wt_tel_mat)
hc <- hclust(d)

jpeg(filename = "telomere_dendogram.jpeg",height=2500,units = "px",width=3400,res=300)
plot(hc,  ylim=c(0,1500),cex=1.2, ylab="Interaction frequency")
dev.off()
