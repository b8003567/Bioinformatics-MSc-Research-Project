library(FactoMineR)
library(factoextra)

data = read_csv("/media/aneta/Data/students/jack/chrom_states_v2/polycomb_cancer.csv", col_types = cols(col_factor()))

metadata = read_csv("/media/aneta/Data/students/jack/chrom_states/chrom_state_meta_127_130_v2_cancer.csv", show_col_types = F)
data = t(data)

metadata$Group = metadata$`Cell Type`

for(i in 1:nrow(metadata)){
  if(metadata$`Cell Lineage`[i] == "Myeloid"){
    metadata$Group[i] = "Myeloid"
  }
}

metadata = subset(metadata, select = "Group")

data = cbind(data, metadata)

res = MCA(data)

fviz_mca_ind(res,
             label = "none",
             habillage = 15002,
             ggtheme = theme_minimal(),
             palette = c("#FF7256", "#A2CD5A", "#BF3EFF"),
             xlim = c(-0.5,0.5),
             ylim = c(-0.5,0.5))

