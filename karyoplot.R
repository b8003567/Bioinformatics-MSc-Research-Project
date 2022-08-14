library(karyoploteR)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

MYC.region = GRanges("chr8:127000000-130000000")
kp <- plotKaryotype(zoom = MYC.region)

genes.data <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg38.knownGene,
                                    karyoplot=kp,
                                    plot.transcripts = TRUE, 
                                    plot.transcripts.structure = TRUE)

genes.data <- addGeneNames(genes.data)
genes.data <- mergeTranscripts(genes.data)

pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 15
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10

kp <- plotKaryotype(zoom = MYC.region, plot.params = pp)
kpAddBaseNumbers(kp, tick.dist = 1000000, minor.tick.dist = 500000,
                 add.units = TRUE, digits = 3)

myc_enhancers = as.data.frame(matrix(nrow = 21, ncol = 3))
names(myc_enhancers) = c("CHR", "START", "END")
myc_enhancers$CHR = rep("chr8", 21)
myc_enhancers$START = c("127117755", "127917270", "128251139", "129546728", "129149529")
myc_enhancers$END = c("127287755", "127918105", "128251614", "129712179", "129189766")
myc_enhancers_gr = makeGRangesFromDataFrame(myc_enhancers, seqnames.field = "CHR", start.field = "START", end.field = "END")

colours_main = c("#629363", "#8D5B96", "#8D5B96", "#8D5B96", "#EFDF0F")

kpPlotRegions(kp, data = myc_enhancers_gr, r0=0.15, r1 = 0.2,
              col = c(colours_main[1], colours_main[1], colours_main[1], colours_main[2], colours_main[3], colours_main[3], colours_main[4], colours_main[5], colours_main[6], colours_main[7], colours_main[8], colours_main[9], colours_main[9], colours_main[9], colours_main[10], colours_main[10], colours_main[11], colours_main[12], colours_main[13], colours_main[14], colours_main[15]))

kpPlotGenes(kp, data=genes.data, r0=0, r1=0.15, gene.name.cex = 0.3, gene.name.position = "left")
