data = read_csv("/media/aneta/Data/students/jack/chrom_states/chrom_states_127_130_v2_healthy.csv")

H3K27ac = read_csv("/media/aneta/Data/students/jack/haem_histone_marks_127_130/MYC_H3K27ac.csv")
H3K27me3 = read_csv("/media/aneta/Data/students/jack/haem_histone_marks_127_130/MYC_H3K27me3.csv")
H3K4me1 = read_csv("/media/aneta/Data/students/jack/haem_histone_marks_127_130/MYC_H3K4me1.csv")
H3K4me3 = read_csv("/media/aneta/Data/students/jack/haem_histone_marks_127_130/MYC_H3K4me3.csv")

H3K27ac = H3K27ac[, grepl("_H_", colnames(H3K27ac))]
H3K27me3 = H3K27me3[, grepl("_H_", colnames(H3K27me3))]
H3K4me1 = H3K4me1[, grepl("_H_", colnames(H3K4me1))]
H3K4me3 = H3K4me3[, grepl("_H_", colnames(H3K4me3))]

H3K27ac_cols = colnames(H3K27ac)
H3K27me3_cols = colnames(H3K27me3)
H3K4me1_cols = colnames(H3K4me1)
H3K4me3_cols = colnames(H3K4me3)

for(i in 1:length(H3K27ac_cols)){
  temp = sub(".*v1.1.4.", "", H3K27ac_cols[i])
  H3K27ac_cols[i] = sub("\\..*", "", temp)
}

for(i in 1:length(H3K27me3_cols)){
  temp = sub(".*v1.1.4.", "", H3K27me3_cols[i])
  H3K27me3_cols[i] = sub("\\..*", "", temp)
}

for(i in 1:length(H3K4me1_cols)){
  temp = sub(".*v1.1.4.", "", H3K4me1_cols[i])
  H3K4me1_cols[i] = sub("\\..*", "", temp)
}

for(i in 1:length(H3K4me3_cols)){
  temp = sub(".*v1.1.4.", "", H3K4me3_cols[i])
  H3K4me3_cols[i] = sub("\\..*", "", temp)
}

shared_cols = Reduce(intersect, list(H3K27ac_cols, H3K27me3_cols, H3K4me1_cols, H3K4me3_cols))

H3K27ac_df = rep(NA, nrow(H3K27ac))
for(i in 1:length(H3K27ac_cols)){
  for(j in 1:length(shared_cols)){
    if(grepl(shared_cols[j], H3K27ac_cols[i]) == TRUE){
      H3K27ac_df = cbind(H3K27ac_df, H3K27ac[i])
    }
  }
}
H3K27ac = H3K27ac_df[2:ncol(H3K27ac_df)]
remove(H3K27ac_df)

H3K27me3_df = rep(NA, nrow(H3K27me3))
for(i in 1:length(H3K27me3_cols)){
  for(j in 1:length(shared_cols)){
    if(grepl(shared_cols[j], H3K27me3_cols[i]) == TRUE){
      H3K27me3_df = cbind(H3K27me3_df, H3K27me3[i])
    }
  }
}
H3K27me3 = H3K27me3_df[2:ncol(H3K27me3_df)]
remove(H3K27me3_df)

H3K4me1_df = rep(NA, nrow(H3K4me1))
for(i in 1:length(H3K4me1_cols)){
  for(j in 1:length(shared_cols)){
    if(grepl(shared_cols[j], H3K4me1_cols[i]) == TRUE){
      H3K4me1_df = cbind(H3K4me1_df, H3K4me1[i])
    }
  }
}
H3K4me1 = H3K4me1_df[2:ncol(H3K4me1_df)]
remove(H3K4me1_df)

H3K4me3_df = rep(NA, nrow(H3K4me3))
for(i in 1:length(H3K4me3_cols)){
  for(j in 1:length(shared_cols)){
    if(grepl(shared_cols[j], H3K4me3_cols[i]) == TRUE){
      H3K4me3_df = cbind(H3K4me3_df, H3K4me3[i])
    }
  }
}
H3K4me3 = H3K4me3_df[2:ncol(H3K4me3_df)]
remove(H3K4me3_df)

metadata = data.frame(matrix(ncol=4, nrow=ncol(H3K27ac)))
colnames(metadata) = c("ID", "Health", "Cell Lineage", "Cell Type")

for(i in 1:length(colnames(H3K27ac))){
  temp = sub(".*v1.1.4.", "", colnames(H3K27ac)[i])
  metadata$ID[i] = sub("\\..*", "", temp)
}

for(i in 1:nrow(metadata)){
  metadata$Health[i] = sub("_.*", "", substring(colnames(H3K27ac)[i], 6))
}

for(i in 1:nrow(metadata)){
  metadata$`Cell Type`[i] = sub("_.*", "", substring(colnames(H3K27ac)[i], 10))
}

for(i in 1:nrow(metadata)){
  metadata$`Cell Lineage`[i] = sub("_.*", "", substring(colnames(H3K27ac)[i], 8))
}

health_status = c("C" = "Cancer", 
                  "H" = "Healthy")

cell_type = c("Bcell"                              = "B cell", 
              "eosinophil"                         = "Eosinophil", 
              "erythroblast"                       = "Erythroblast", 
              "erythroid-progenitor"               = "Erythroid progenitor", 
              "granulocyte-monocyte-progenitor"    = "Granulocyte-monocyte progenitor", 
              "macrophage"                         = "Macrophage", 
              "megakaryocyte"                      = "Megakaryocyte", 
              "megakaryocyte-erythroid-progenitor" = "Megakaryocyte-erythroid progenitor", 
              "monocyte"                           = "Monocyte", 
              "myeloid"                            = "Myeloid", 
              "neutrophil"                         = "Neutrophil", 
              "neutrophilic-metamyelocyte"         = "Neutrophilic metamyelocyte", 
              "neutrophilic-myelocyte"             = "Neutrophilic myelocyte", 
              "NK"                                 = "NK cell", 
              "progenitor"                         = "Progenitor", 
              "Tcell"                              = "T cell",
              "MNC"                                = "Mononuclear cell")

cell_lineage = c("L" = "Lymphoid",
                 "M" = "Myeloid",
                 "MNC" = "Mononuclear cell")

metadata$Health = health_status[metadata$Health]
metadata$`Cell Lineage` = cell_lineage[metadata$`Cell Lineage`]
metadata$`Cell Type` = cell_type[metadata$`Cell Type`]

write_csv(metadata, "/media/aneta/Data/students/jack/chrom_states/chrom_state_meta_127_130_v2_healthy.csv")
