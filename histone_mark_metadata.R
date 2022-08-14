library(readr)

H3K4me1 = read_csv("/media/aneta/Data/students/jack/haem_histone_marks_127_130/MYC_H3K4me1.csv", show_col_types = F)
H3K4me3 = read_csv("/media/aneta/Data/students/jack/haem_histone_marks_127_130/MYC_H3K4me3.csv", show_col_types = F)
H3K9me3 = read_csv("/media/aneta/Data/students/jack/haem_histone_marks_127_130/MYC_H3K9me3.csv", show_col_types = F)
H3K27ac = read_csv("/media/aneta/Data/students/jack/haem_histone_marks_127_130/MYC_H3K27ac.csv", show_col_types = F)
H3K27me3 = read_csv("/media/aneta/Data/students/jack/haem_histone_marks_127_130/MYC_H3K27me3.csv", show_col_types = F)
H3K36me3 = read_csv("/media/aneta/Data/students/jack/haem_histone_marks_127_130/MYC_H3K36me3.csv", show_col_types = F)

meta_H3K4me1 = as.data.frame(matrix(ncol=4, nrow = ncol(H3K4me1)))
meta_H3K4me3 = as.data.frame(matrix(ncol=4, nrow = ncol(H3K4me3)))
meta_H3K9me3 = as.data.frame(matrix(ncol=4, nrow = ncol(H3K9me3)))
meta_H3K27ac = as.data.frame(matrix(ncol=4, nrow = ncol(H3K27ac)))
meta_H3K27me3 = as.data.frame(matrix(ncol=4, nrow = ncol(H3K27me3)))
meta_H3K36me3 = as.data.frame(matrix(ncol=4, nrow = ncol(H3K36me3)))

colnames(meta_H3K4me1) = c("ID", "Health", "Cell Type", "Cell Lineage")
colnames(meta_H3K4me3) = c("ID", "Health", "Cell Type", "Cell Lineage")
colnames(meta_H3K9me3) = c("ID", "Health", "Cell Type", "Cell Lineage")
colnames(meta_H3K27ac) = c("ID", "Health", "Cell Type", "Cell Lineage")
colnames(meta_H3K27me3) = c("ID", "Health", "Cell Type", "Cell Lineage")
colnames(meta_H3K36me3) = c("ID", "Health", "Cell Type", "Cell Lineage")

meta_H3K4me1$ID = colnames(H3K4me1)
meta_H3K4me3$ID = colnames(H3K4me3)
meta_H3K9me3$ID = colnames(H3K9me3)
meta_H3K27ac$ID = colnames(H3K27ac)
meta_H3K27me3$ID = colnames(H3K27me3)
meta_H3K36me3$ID = colnames(H3K36me3)

#===============================================================================

for(i in 1:nrow(meta_H3K4me1)){
  temp = sub(".*v1.1.4.", "", meta_H3K4me1$ID[i])
  meta_H3K4me1$ID[i] = sub("\\..*", "", temp)
}

for(i in 1:nrow(meta_H3K4me3)){
  temp = sub(".*v1.1.4.", "", meta_H3K4me3$ID[i])
  meta_H3K4me3$ID[i] = sub("\\..*", "", temp)
}

for(i in 1:nrow(meta_H3K9me3)){
  temp = sub(".*v1.1.4.", "", meta_H3K9me3$ID[i])
  meta_H3K9me3$ID[i] = sub("\\..*", "", temp)
}

for(i in 1:nrow(meta_H3K27ac)){
  temp = sub(".*v1.1.4.", "", meta_H3K27ac$ID[i])
  meta_H3K27ac$ID[i] = sub("\\..*", "", temp)
}

for(i in 1:nrow(meta_H3K27me3)){
  temp = sub(".*v1.1.4.", "", meta_H3K27me3$ID[i])
  meta_H3K27me3$ID[i] = sub("\\..*", "", temp)
}

for(i in 1:nrow(meta_H3K36me3)){
  temp = sub(".*v1.1.4.", "", meta_H3K36me3$ID[i])
  meta_H3K36me3$ID[i] = sub("\\..*", "", temp)
}

#===============================================================================

for(i in 1:nrow(meta_H3K4me1)){
  meta_H3K4me1$Health[i] = sub("_.*", "", substring(colnames(H3K4me1)[i], 6))
}

for(i in 1:nrow(meta_H3K4me3)){
  meta_H3K4me3$Health[i] = sub("_.*", "", substring(colnames(H3K4me3)[i], 6))
}

for(i in 1:nrow(meta_H3K9me3)){
  meta_H3K9me3$Health[i] = sub("_.*", "", substring(colnames(H3K9me3)[i], 6))
}

for(i in 1:nrow(meta_H3K27ac)){
  meta_H3K27ac$Health[i] = sub("_.*", "", substring(colnames(H3K27ac)[i], 6))
}

for(i in 1:nrow(meta_H3K27me3)){
  meta_H3K27me3$Health[i] = sub("_.*", "", substring(colnames(H3K27me3)[i], 6))
}

for(i in 1:nrow(meta_H3K36me3)){
  meta_H3K36me3$Health[i] = sub("_.*", "", substring(colnames(H3K36me3)[i], 6))
}

#===============================================================================

for(i in 1:nrow(meta_H3K4me1)){
  meta_H3K4me1$`Cell Type`[i] = sub("_.*", "", substring(colnames(H3K4me1)[i], 10))
}

for(i in 1:nrow(meta_H3K4me3)){
  meta_H3K4me3$`Cell Type`[i] = sub("_.*", "", substring(colnames(H3K4me3)[i], 10))
}

for(i in 1:nrow(meta_H3K9me3)){
  meta_H3K9me3$`Cell Type`[i] = sub("_.*", "", substring(colnames(H3K9me3)[i], 10))
}

for(i in 1:nrow(meta_H3K27ac)){
  meta_H3K27ac$`Cell Type`[i] = sub("_.*", "", substring(colnames(H3K27ac)[i], 10))
}

for(i in 1:nrow(meta_H3K27me3)){
  meta_H3K27me3$`Cell Type`[i] = sub("_.*", "", substring(colnames(H3K27me3)[i], 10))
}

for(i in 1:nrow(meta_H3K36me3)){
  meta_H3K36me3$`Cell Type`[i] = sub("_.*", "", substring(colnames(H3K36me3)[i], 10))
}

#===============================================================================

for(i in 1:nrow(meta_H3K4me1)){
  meta_H3K4me1$`Cell Lineage`[i] = sub("_.*", "", substring(colnames(H3K4me1)[i], 8))
}

for(i in 1:nrow(meta_H3K4me3)){
  meta_H3K4me3$`Cell Lineage`[i] = sub("_.*", "", substring(colnames(H3K4me3)[i], 8))
}

for(i in 1:nrow(meta_H3K9me3)){
  meta_H3K9me3$`Cell Lineage`[i] = sub("_.*", "", substring(colnames(H3K9me3)[i], 8))
}

for(i in 1:nrow(meta_H3K27ac)){
  meta_H3K27ac$`Cell Lineage`[i] = sub("_.*", "", substring(colnames(H3K27ac)[i], 8))
}

for(i in 1:nrow(meta_H3K27me3)){
  meta_H3K27me3$`Cell Lineage`[i] = sub("_.*", "", substring(colnames(H3K27me3)[i], 8))
}

for(i in 1:nrow(meta_H3K36me3)){
  meta_H3K36me3$`Cell Lineage`[i] = sub("_.*", "", substring(colnames(H3K36me3)[i], 8))
}

#===============================================================================

for(i in 1:nrow(meta_H3K4me1)){
  if(meta_H3K4me1$`Cell Lineage`[i] == "MNC"){
    meta_H3K4me1$`Cell Type`[i] = "MNC"
  }
}

for(i in 1:nrow(meta_H3K4me3)){
  if(meta_H3K4me3$`Cell Lineage`[i] == "MNC"){
    meta_H3K4me3$`Cell Type`[i] = "MNC"
  }
}

for(i in 1:nrow(meta_H3K9me3)){
  if(meta_H3K9me3$`Cell Lineage`[i] == "MNC"){
    meta_H3K9me3$`Cell Type`[i] = "MNC"
  }
}

for(i in 1:nrow(meta_H3K27ac)){
  if(meta_H3K27ac$`Cell Lineage`[i] == "MNC"){
    meta_H3K27ac$`Cell Type`[i] = "MNC"
  }
}

for(i in 1:nrow(meta_H3K27me3)){
  if(meta_H3K27me3$`Cell Lineage`[i] == "MNC"){
    meta_H3K27me3$`Cell Type`[i] = "MNC"
  }
}

for(i in 1:nrow(meta_H3K36me3)){
  if(meta_H3K36me3$`Cell Lineage`[i] == "MNC"){
    meta_H3K36me3$`Cell Type`[i] = "MNC"
  }
}

#===============================================================================

health_map = c("C" = "Cancer",
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

#===============================================================================

meta_H3K4me1$Health = health_map[meta_H3K4me1$Health]
meta_H3K4me3$Health = health_map[meta_H3K4me3$Health]
meta_H3K9me3$Health = health_map[meta_H3K9me3$Health]
meta_H3K27ac$Health = health_map[meta_H3K27ac$Health]
meta_H3K27me3$Health = health_map[meta_H3K27me3$Health]
meta_H3K36me3$Health = health_map[meta_H3K36me3$Health]

meta_H3K4me1$`Cell Type` = cell_type[meta_H3K4me1$`Cell Type`]
meta_H3K4me3$`Cell Type` = cell_type[meta_H3K4me3$`Cell Type`]
meta_H3K9me3$`Cell Type` = cell_type[meta_H3K9me3$`Cell Type`]
meta_H3K27ac$`Cell Type` = cell_type[meta_H3K27ac$`Cell Type`]
meta_H3K27me3$`Cell Type` = cell_type[meta_H3K27me3$`Cell Type`]
meta_H3K36me3$`Cell Type` = cell_type[meta_H3K36me3$`Cell Type`]

meta_H3K4me1$`Cell Lineage` = cell_lineage[meta_H3K4me1$`Cell Lineage`]
meta_H3K4me3$`Cell Lineage` = cell_lineage[meta_H3K4me3$`Cell Lineage`]
meta_H3K9me3$`Cell Lineage` = cell_lineage[meta_H3K9me3$`Cell Lineage`]
meta_H3K27ac$`Cell Lineage` = cell_lineage[meta_H3K27ac$`Cell Lineage`]
meta_H3K27me3$`Cell Lineage` = cell_lineage[meta_H3K27me3$`Cell Lineage`]
meta_H3K36me3$`Cell Lineage` = cell_lineage[meta_H3K36me3$`Cell Lineage`]

write_csv(meta_H3K4me1, file = "/media/aneta/Data/students/jack/metadata/meta_H3K4me1_v3.csv")
write_csv(meta_H3K4me3, file = "/media/aneta/Data/students/jack/metadata/meta_H3K4me3_v3.csv")
write_csv(meta_H3K9me3, file = "/media/aneta/Data/students/jack/metadata/meta_H3K9me3_v3.csv")
write_csv(meta_H3K27ac, file = "/media/aneta/Data/students/jack/metadata/meta_H3K27ac_v3.csv")
write_csv(meta_H3K27me3, file = "/media/aneta/Data/students/jack/metadata/meta_H3K27me3_v3.csv")
write_csv(meta_H3K36me3, file = "/media/aneta/Data/students/jack/metadata/meta_H3K36me3_v3.csv")
