library(pheatmap)
library(readr)
library(dplyr)

data = read_csv("/media/aneta/Data/students/jack/haem_histone_marks_127_130/MYC_H3K4me1.csv", show_col_types = F)
metadata = read_csv("/media/aneta/Data/students/jack/metadata/meta_H3K4me1_v3.csv", show_col_types = F)
colnames(data) = metadata$ID
data = t(data)

heat = pheatmap(data, cluster_cols = F, show_rownames = F)

order = heat$tree_row$order

metadata$ID = metadata$ID[order]
metadata$Health = metadata$Health[order]
metadata$`Cell Type` = metadata$`Cell Type`[order]
metadata$`Cell Lineage` = metadata$`Cell Lineage`[order]

metadata$Group = metadata$`Cell Type`

for(i in 1:nrow(metadata)){
  if(metadata$`Cell Lineage`[i] == "Myeloid"){
    metadata$Group[i] = "Myeloid"
  }
}

id = metadata$ID

metadata_2 = as.data.frame(metadata$Health) %>% cbind(metadata$Group)
colnames(metadata_2) = c("Health", "Group")
rownames(metadata_2) = id

colours = list("Group" = c("T cell" = "darkorange2",
                           "B cell" = "darkgoldenrod4",
                           "NK cell" = "darkgoldenrod1",
                           "Myeloid" = "deepskyblue"),
               "Health" = c("Cancer" = "red", 
                            "Healthy" = "chartreuse2"))
heat = pheatmap(data, 
                cluster_cols = F,
                show_rownames = F,
                annotation_row = metadata_2,
                annotation_colors = colours,
                clustering_distance_rows = "euclidean",
                legend = F)
