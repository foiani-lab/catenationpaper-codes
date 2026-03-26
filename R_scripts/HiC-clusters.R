library(HiCExperiment)
library(HiCool)
library(HiContacts)
library(purrr)
library(GenomicRanges)
library(dplyr)
library(cowplot)
library(ggplot2)

hics <- import("wt.mcool",format = 'mcool', resolution = 1000)


#### centromere clustering

centros_yeast <- makeGRangesFromDataFrame(data.frame(read.table("centromere.bed",header=T,sep="\t")),seqnames.field="chrom",start.field="start",end.field="end")

centros_pairs <- lapply(1:length(centros_yeast), function(i) {                                                                
    lapply(1:length(centros_yeast), function(j) {                                                                               
        S4Vectors::Pairs(centros_yeast[i], centros_yeast[j])                                                                    
    })                                                                                                                          
}) |>
	do.call(c, args = _) |> 
	do.call(c, args = _) |>
	InteractionSet::makeGInteractionsFromGRangesPairs()

centros_pairs <- centros_pairs[anchors(centros_pairs, 'first') != anchors(centros_pairs, 'second')]                               

aggr_maps <- purrr::imap(hics, ~ {
    aggr <- aggregate(.x, centros_pairs, maxDistance = 1e999)
    plotMatrix(
        aggr, use.scores = 'balanced', limits = c(-6, -1), 
        cmap = HiContacts::rainbowColors(), 
        caption = FALSE
    ) 
})

jpeg("centromere-interactions-wt.jpeg",width=10,height=8,unit="in",res=300)
cowplot::plot_grid(plotlist = aggr_maps, nrow = 1) 
dev.off()



#### rDNA locus interactions with chrIV

rDNA_locus <- makeGRangesFromDataFrame(data.frame(read.table("rDNA-locus.bed",header=T,sep="\t")),seqnames.field="chrom",start.field="start",end.field="end")


                                                                  
v4c_rDNA <- imap_dfr(hics, ~ virtual4C(.x, GenomicRanges::resize(rDNA_locus[1], 10)) |> 
    as_tibble() |> 
    mutate(background = .y) |> 
    filter(seqnames == 'chrIV')
)


jpeg("rDNAclusters-chrIV-WT.jpeg",width=10,height=6,unit="in",res=300)
ggplot(v4c_rDNA, aes(x = start, y = score, fill = background)) + geom_area() + theme_bw() + labs( x = "chrIV position",y = "Contacts from rDNA locus",title = "Interaction profile of rDNA clusters in chrIV") + coord_cartesian(ylim = c(0, 0.001))
dev.off()
