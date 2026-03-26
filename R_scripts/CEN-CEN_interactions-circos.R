library("circlize")
library(data.table)

chrrom <- function(data){
  chrrom = paste0("chr",c("X","XI","XII","XIII","XIV","XV","XVI","I","II","III","IV","V","VI","VII","VIII","IX"))
  chrnum = c("10","11","12","13","14","15","16","1","2","3","4","5","6","7","8","9")
  data$V1 = as.vector(data$V1)
  data$V7 = as.vector(data$V7)
  for (i in 1:length(chrnum)){
    data$V1 = gsub(chrnum[i], chrrom[i], data$V1)
    data$V7 = gsub(chrnum[i], chrrom[i], data$V7)
  }
  return(data)
}

circosallcenplot <- function(allinterdata, cent, outdir, tag){
  #allinterwt = allinterwt[allinterwt$V1 < allinterwt$V7,]
  allinterdata = allinterdata[((allinterdata$V6 != ".") & (allinterdata$V6 != "-1")),]
  allinterdata = chrrom(allinterdata)
  allinterdata = allinterdata[order(-allinterdata$V4),]
  allinterdataf = allinterdata[c(1:round(dim(allinterdata)[1]*0.1)),]
  allinterdataf = allinterdataf[allinterdataf$V4 >= 25,]
  coldat = data.frame(V1 = c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI"), colr = c(1:16))
  allinterdataf = merge(allinterdataf,coldat,by=c("V1"))
  bed1 = allinterdataf[,c(1,2,3,4)]
  colnames(bed1) = c("chr","start","end","value1")
  #bed1$end = bed1$start + 1000
  bed2 = allinterdataf[,c(7,8,9,10)]
  colnames(bed2) = c("chr","start","end","value1")
  filename1 = paste(outdir,tag,".jpeg",sep="")
  jpeg(filename1,height = 4000, width = 6000, res=800)
  circos.initializeWithIdeogram(species = "sacCer3",chromosome.index = paste0("chr",c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI")))
  circos.genomicTrack(cent, ylim = c(0, 1), panel.fun = function(region, value, ...) {
    circos.genomicText(region, value, y = 0.6, cex = 0.25, facing = "reverse.clockwise", labels.column = 1, ...)
  })
  circos.genomicLink(bed1, bed2, col = allinterdataf$colr, border = NA)
  dev.off()
}

circoscenplot <- function(allinterdata, cent, outdir, tag){
  #allinterwt = allinterwt[allinterwt$V1 < allinterwt$V7,]
  allinterdata = allinterdata[((allinterdata$V1 == "4") & (allinterdata$V6 == "CEN4")),]
  allinterdata = chrrom(allinterdata)
  allinterdata = allinterdata[order(-allinterdata$V4),]
  allinterdataf = allinterdata[c(1:round(dim(allinterdata)[1]*0.1)),]
  allinterdataf = allinterdataf[allinterdataf$V4 >= 25,]
  bed1 = allinterdataf[,c(1,2,3,4)]
  colnames(bed1) = c("chr","start","end","value1")
  #bed1$end = bed1$start + 1000
  bed2 = allinterdataf[,c(7,8,9,10)]
  colnames(bed2) = c("chr","start","end","value1")
  filename1 = paste(outdir,tag,".jpeg",sep="")
  jpeg(filename1,height = 4000, width = 6000, res=800)
  circos.initializeWithIdeogram(species = "sacCer3",chromosome.index = paste0("chr",c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI")))
  circos.genomicTrack(cent, ylim = c(0, 1), panel.fun = function(region, value, ...) {
    circos.genomicText(region, value, y = 0.6, cex = 0.25, facing = "reverse.clockwise", labels.column = 1, ...)
  })
  circos.genomicLink(bed1, bed2, col = "red", border = NA)
  dev.off()
}

cent = read.delim(file = "centromere.bed", h = F,skip = 1)
cent = cent[,c(1:4)]
colnames(cent) = c("chr","start","end","label")

outdir = "/home/hpc-344/Data/2_Analysis/Plots/Centromere_clustering/"
allinterwt = fread(file = "/home/hpc-344/Data/2_Analysis/Data/wt/allinter_obsvc_5kb/wt_interinteractions_withid.txt",sep="\t",h=F)

circoscenplot(allinterwt, cent, outdir, "WT_cen4")
circosallcenplot(allinterwt, cent, outdir, "WT_allcen")
