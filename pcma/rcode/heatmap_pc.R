if (!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2",repos='https://cloud.r-project.org/',dependencies=TRUE)
    library(ggplot2)
}
if (!require(dplyr, quietly = TRUE)) {
    install.packages("dplyr",repos='https://cloud.r-project.org/', dependencies=TRUE)
    library(dplyr)
}
if (!require(tidyr, quietly = TRUE)) {
    install.packages("tidyr",repos='https://cloud.r-project.org/', dependencies=TRUE)
    library(tidyr)
}
args = commandArgs(trailingOnly = TRUE)
data_dir = args[1]
plot_dir = args[2]

CD = read.csv(data_dir)
CD = CD[,c(2,1)]
colnames(CD) = c("Bacteria","PC_note")
CD = CD[which(CD$Bacteria %in% unique(CD$Bacteria)[1:10]),]
CD[which(CD$PC_note == "no_sig"),2] = ""
pc_CD = expand.grid(Bacteria = unique(CD$Bacteria), PC = paste0("PC", 1:50))

# 标记PC_note中出现的PC和no
pc_CD$Present = with(CD, mapply(function(bact, pc) {
  any(Bacteria == bact & PC_note == pc)
}, pc_CD$Bacteria, pc_CD$PC))
pc_CD$IsNo = with(CD, mapply(function(bact, pc) {
  any(Bacteria == bact & PC_note == "no")
}, pc_CD$Bacteria, pc_CD$PC))

# 颜色转换
pc_CD$Color = ifelse(pc_CD$IsNo, "lightgrey", ifelse(pc_CD$Present, "red", "white"))
plot_dir_all = file.path(plot_dir, "heatmap.pdf")

p = ggplot(pc_CD, aes(x = Bacteria, y = PC, fill = Color)) + 
    geom_tile(color = "grey") +
    scale_fill_manual(values = c("red" = "red", "white" = "white", "lightgrey" = "lightgrey")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_text(size = 9, color = "black"),
          legend.position = "none") +
    coord_fixed(ratio = 1)

if (nrow(pc_CD) > 0) {
  pdf(plot_dir_all, width=8.3, height=11.7) 
  print(p)
  dev.off()
}else{
  print("there is no significant PC")
}


