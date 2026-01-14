# heatmap
if (!require(pheatmap, quietly = TRUE)) {
  install.packages("pheatmap",repos='https://cloud.r-project.org/',dependencies=TRUE)
  library(pheatmap)
}

args = commandArgs(trailingOnly = TRUE)
PCA_data = args[1]  
plot_dir = args[2]

meta = read.csv(PCA_data,header = T,check.names = F)
rownames(meta) <- paste0("PC", 1:min(50,nrow(meta)))
select_df = function(meta){
  num_cols = ncol(meta)
  selected_cols = if (num_cols > 25) {
    sample(seq_len(num_cols), 20)
  } else {
    seq_len(num_cols)
  }
  num_rows = nrow(meta)
  selected_rows = 1:15
  return(meta[selected_rows, selected_cols])
}

meta = select_df(meta)

pdf(plot_dir,width=6, height=5)
pheatmap(meta,
         cluster_rows = F,
         cluster_cols = F,
         border = NA,
         legend = T,
         cellwidth = 10, 
         cellheight = 10,
         fontsize = 7)
dev.off()