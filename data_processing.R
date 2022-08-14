library(dplyr)
library(tidyverse)
library(GenomicRanges)

PEAK_FILES = list.files("/media/aneta/Data/students/jack/haem_bed_files_v2", recursive = T, full.names = T, pattern = ".bed")

MYC = data.frame(CHR="chr8", START="127000000", END="130000000")
MYC_gr = GRanges(MYC, seqnames.field = 'CHR', start.field = 'START', end.field = 'END')
MYC_gr = tile(MYC_gr, width=200)
MYC_df = as.data.frame(MYC_gr)[, 3:7]

MYC_H3K9me3 = MYC_df
MYC_H3K4me3 = MYC_df
MYC_H3K4me1 = MYC_df
MYC_H3K27me3 = MYC_df
MYC_H3K36me3 = MYC_df
MYC_H3K27ac = MYC_df

for (sample in PEAK_FILES){
  NAME           <- basename(sample) %>% gsub(".bed","",.)
  if (file.info(sample)$size == 0){
    combined = as.data.frame(rep(0, nrow(MYC_df)))
    colnames(combined) = "Overlap"
  } else {
    DATASET        <- data.frame(NA,NA,NA,NA,NA)
    names(DATASET) <- c("CHR", "START", "END", "PEAK", "ID")
    
    FILE          <- read.table(sample) %>% dplyr::select(V1, V2, V3, V4) 
    names(FILE)   <- c("CHR", "START", "END", "PEAK")
    FILE          <- mutate(FILE, ID = NAME)
    
    DATASET       <- rbind(DATASET, FILE)
    DATASET       <- filter(DATASET, ID !="NA")
    DATASET       <- arrange(DATASET, CHR, START, END, PEAK, ID)
    
    gr = makeGRangesFromDataFrame(DATASET, seqnames.field = "CHR", start.field = "START", end.field = "END", keep.extra.columns = TRUE)
    
    overlap = subsetByOverlaps(MYC_gr[[1]], gr, minoverlap = 100L)
    no_overlap = MYC_gr[[1]][!MYC_gr[[1]] %in% overlap]
    overlap = as.data.frame(overlap)
    overlap$Overlap = rep(1, nrow(overlap))
    no_overlap = as.data.frame(no_overlap)
    no_overlap$Overlap = rep(0, nrow(no_overlap))
    
    combined = rbind(overlap, no_overlap)
    combined = combined[order(combined$start), ]
  }
  
  if (grepl("H3K9me3", NAME) == TRUE){
    MYC_H3K9me3 = cbind(MYC_H3K9me3, H3K9me3 = combined$Overlap)
    names(MYC_H3K9me3)[names(MYC_H3K9me3) == "H3K9me3"] = NAME
  } else if (grepl("H3K4me3", NAME) == TRUE){
    MYC_H3K4me3 = cbind(MYC_H3K4me3, H3K4me3 = combined$Overlap)
    names(MYC_H3K4me3)[names(MYC_H3K4me3) == "H3K4me3"] = NAME
  } else if (grepl("H3K4me1", NAME) == TRUE){
    MYC_H3K4me1 = cbind(MYC_H3K4me1, H3K4me1 = combined$Overlap)
    names(MYC_H3K4me1)[names(MYC_H3K4me1) == "H3K4me1"] = NAME
  } else if (grepl("H3K27me3", NAME) == TRUE){
    MYC_H3K27me3 = cbind(MYC_H3K27me3, H3K27me3 = combined$Overlap)
    names(MYC_H3K27me3)[names(MYC_H3K27me3) == "H3K27me3"] = NAME
  } else if (grepl("H3K36me3", NAME) == TRUE){
    MYC_H3K36me3 = cbind(MYC_H3K36me3, H3K36me3 = combined$Overlap)
    names(MYC_H3K36me3)[names(MYC_H3K36me3) == "H3K36me3"] = NAME
  } else if (grepl("H3K27ac", NAME) == TRUE){
    MYC_H3K27ac = cbind(MYC_H3K27ac, H3K27ac = combined$Overlap)
    names(MYC_H3K27ac)[names(MYC_H3K27ac) == "H3K27ac"] = NAME
  }
}

MYC_H3K27ac <- MYC_H3K27ac[, !duplicated(colnames(MYC_H3K27ac))]
MYC_H3K27me3 <- MYC_H3K27me3[, !duplicated(colnames(MYC_H3K27me3))]
MYC_H3K36me3 <- MYC_H3K36me3[, !duplicated(colnames(MYC_H3K36me3))]
MYC_H3K4me1 <- MYC_H3K4me1[, !duplicated(colnames(MYC_H3K4me1))]
MYC_H3K4me3 <- MYC_H3K4me3[, !duplicated(colnames(MYC_H3K4me3))]
MYC_H3K9me3 <- MYC_H3K9me3[, !duplicated(colnames(MYC_H3K9me3))]

MYC_H3K4me1 = MYC_H3K4me1[, 6:ncol(MYC_H3K4me1)]
MYC_H3K4me3 = MYC_H3K4me3[, 6:ncol(MYC_H3K4me3)]
MYC_H3K9me3 = MYC_H3K9me3[, 6:ncol(MYC_H3K9me3)]
MYC_H3K27ac = MYC_H3K27ac[, 6:ncol(MYC_H3K27ac)]
MYC_H3K27me3 = MYC_H3K27me3[, 6:ncol(MYC_H3K27me3)]
MYC_H3K36me3 = MYC_H3K36me3[, 6:ncol(MYC_H3K36me3)]

MYC_H3K4me1 = MYC_H3K4me1[, grepl('_C_|_H_', colnames(MYC_H3K4me1))]
MYC_H3K4me3 = MYC_H3K4me3[, grepl('_C_|_H_', colnames(MYC_H3K4me3))]
MYC_H3K9me3 = MYC_H3K9me3[, grepl('_C_|_H_', colnames(MYC_H3K9me3))]
MYC_H3K27ac = MYC_H3K27ac[, grepl('_C_|_H_', colnames(MYC_H3K27ac))]
MYC_H3K27me3 = MYC_H3K27me3[, grepl('_C_|_H_', colnames(MYC_H3K27me3))]
MYC_H3K36me3 = MYC_H3K36me3[, grepl('_C_|_H_', colnames(MYC_H3K36me3))]

write_csv(MYC_H3K9me3, "/media/aneta/Data/students/jack/haem_histone_marks_127_130/MYC_H3K9me3.csv")
write_csv(MYC_H3K4me3, "/media/aneta/Data/students/jack/haem_histone_marks_127_130/MYC_H3K4me3.csv")
write_csv(MYC_H3K4me1, "/media/aneta/Data/students/jack/haem_histone_marks_127_130/MYC_H3K4me1.csv")
write_csv(MYC_H3K36me3, "/media/aneta/Data/students/jack/haem_histone_marks_127_130/MYC_H3K36me3.csv")
write_csv(MYC_H3K27me3, "/media/aneta/Data/students/jack/haem_histone_marks_127_130/MYC_H3K27me3.csv")
write_csv(MYC_H3K27ac, "/media/aneta/Data/students/jack/haem_histone_marks_127_130/MYC_H3K27ac.csv")
