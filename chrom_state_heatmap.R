library(readr)
library(dplyr)

data = read_csv("/media/aneta/Data/students/jack/chrom_states/chrom_states_127_130_v2.csv", show_col_types = F) %>% as.data.frame()

H3K27ac = read_csv("/media/aneta/Data/students/jack/haem_histone_marks_127_130/MYC_H3K27ac.csv", show_col_types = F)
H3K27me3 = read_csv("/media/aneta/Data/students/jack/haem_histone_marks_127_130/MYC_H3K27me3.csv", show_col_types = F)
H3K4me1 = read_csv("/media/aneta/Data/students/jack/haem_histone_marks_127_130/MYC_H3K4me1.csv", show_col_types = F)
H3K4me3 = read_csv("/media/aneta/Data/students/jack/haem_histone_marks_127_130/MYC_H3K4me3.csv", show_col_types = F)

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

colnames(data) = shared_cols

metadata = read_csv("/media/aneta/Data/students/jack/chrom_states/chrom_state_meta_127_130_v2.csv", show_col_types = F) %>% as.data.frame()

for(i in 1:nrow(metadata)){
  if(metadata$`Cell Lineage`[i] == "Myeloid"){
    metadata$Group[i] = "Myeloid"
  }
  else{
    metadata$Group[i] = metadata$`Cell Type`[i]
  }
}

metadata = subset(metadata, select = c("ID", "Health", "Group"))

state_1 = data
state_2 = data
state_3 = data
state_4 = data

for(i in 1:ncol(state_1)){
  state_1[,i][state_1[,i] == 2] = 0
  state_1[,i][state_1[,i] == 3] = 0
  state_1[,i][state_1[,i] == 4] = 0
}

for(i in 1:ncol(state_2)){
  state_2[,i][state_2[,i] == 1] = 0
  state_2[,i][state_2[,i] == 2] = 1
  state_2[,i][state_2[,i] == 3] = 0
  state_2[,i][state_2[,i] == 4] = 0
}

for(i in 1:ncol(state_3)){
  state_3[,i][state_3[,i] == 1] = 0
  state_3[,i][state_3[,i] == 2] = 0
  state_3[,i][state_3[,i] == 3] = 1
  state_3[,i][state_3[,i] == 4] = 0
}

for(i in 1:ncol(state_4)){
  state_4[,i][state_4[,i] == 1] = 0
  state_4[,i][state_4[,i] == 2] = 0
  state_4[,i][state_4[,i] == 3] = 0
  state_4[,i][state_4[,i] == 4] = 1
}

state_1_mat = t(state_1)
state_2_mat = t(state_2)
state_3_mat = t(state_3)
state_4_mat = t(state_4)

state_1_heat = pheatmap(state_1_mat, cluster_cols = F, show_rownames = F)
state_2_heat = pheatmap(state_2_mat, cluster_cols = F, show_rownames = F)
state_3_heat = pheatmap(state_3_mat, cluster_cols = F, show_rownames = F)
state_4_heat = pheatmap(state_4_mat, cluster_cols = F, show_rownames = F)

clustering_order_1 = state_1_heat$tree_row$order
metadata_1 = metadata
metadata_1$ID = metadata_1$ID[clustering_order_1]
metadata_1$Health = metadata_1$Health[clustering_order_1]
metadata_1$Group = metadata_1$Group[clustering_order_1]
id_1 = metadata_1$ID
metadata_1 = subset(metadata_1, select = c("Health", "Group"))
rownames(metadata_1) = id_1

clustering_order_2 = state_2_heat$tree_row$order
metadata_2 = metadata
metadata_2$ID = metadata_2$ID[clustering_order_2]
metadata_2$Health = metadata_2$Health[clustering_order_2]
metadata_2$Group = metadata_2$Group[clustering_order_2]
id_2 = metadata_2$ID
metadata_2 = subset(metadata_2, select = c("Health", "Group"))
rownames(metadata_2) = id_2

clustering_order_3 = state_3_heat$tree_row$order
metadata_3 = metadata
metadata_3$ID = metadata_3$ID[clustering_order_3]
metadata_3$Health = metadata_3$Health[clustering_order_3]
metadata_3$Group = metadata_3$Group[clustering_order_3]
id_3 = metadata_3$ID
metadata_3 = subset(metadata_3, select = c("Health", "Group"))
rownames(metadata_3) = id_3

clustering_order_4 = state_4_heat$tree_row$order
metadata_4 = metadata
metadata_4$ID = metadata_4$ID[clustering_order_4]
metadata_4$Health = metadata_4$Health[clustering_order_4]
metadata_4$Group = metadata_4$Group[clustering_order_4]
id_4 = metadata_4$ID
metadata_4 = subset(metadata_4, select = c("Health", "Group"))
rownames(metadata_4) = id_4

colours = list("Group" = c("T cell" = "darkorange2",
                           "B cell" = "darkgoldenrod4",
                           "NK cell" = "darkgoldenrod1",
                           "Myeloid" = "deepskyblue"),
               "Health" = c("Cancer" = "red", 
                            "Healthy" = "chartreuse2"))

state_1_heat = pheatmap(state_1_mat, 
                        cluster_cols = F,
                        show_rownames = F,
                        annotation_row = metadata_1,
                        annotation_colors = colours,
                        clustering_distance_rows = "euclidean",
                        legend = F)

state_2_heat = pheatmap(state_2_mat, 
                        cluster_cols = F,
                        show_rownames = F,
                        annotation_row = metadata_2,
                        annotation_colors = colours,
                        clustering_distance_rows = "euclidean",
                        legend = F)

state_3_heat = pheatmap(state_3_mat, 
                        cluster_cols = F,
                        show_rownames = F,
                        annotation_row = metadata_3,
                        annotation_colors = colours,
                        clustering_distance_rows = "euclidean",
                        legend = F)

state_4_heat = pheatmap(state_4_mat, 
                        cluster_cols = F,
                        show_rownames = F,
                        annotation_row = metadata_4,
                        annotation_colors = colours,
                        clustering_distance_rows = "euclidean",
                        legend = F)