library(dplyr)
library(tidyverse)

PEAK_FILES = list.files("/media/aneta/Data/students/data/IHEC/ChIP-Seq_narrowPeak_edit_v2_220504")
OUTDIR = "/media/aneta/Data/students/jack/bed_files_v2/"

for (sample in PEAK_FILES){
  NAME          <- basename(sample) %>% gsub(".narrowPeak","",.)
  FILE          <- read.table(sample) %>% dplyr::select(V1, V2, V3, V4)
  names(FILE)   <- c("CONTIG", "START", "END", "CHRSTAT")
  FILE          <- mutate(FILE, ID = NAME)
  FILE          <- subset(FILE, (FILE$CONTIG == "chr8" & FILE$START >= 126000000 & FILE$END <= 131000000))
  write_tsv(FILE, paste0(OUTDIR, NAME, ".bed"), col_names = F)
}
