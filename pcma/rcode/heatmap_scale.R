# heatmap
if (!require(pheatmap, quietly = TRUE)) {
  install.packages("pheatmap",repos='https://cloud.r-project.org/',dependencies=TRUE)
  library(pheatmap)
}

args = commandArgs(trailingOnly = TRUE)
raw_data = args[1]
scale_data = args[2]
plot_dir_1 = args[3]
plot_dir_2 = args[4]
meta = read.csv(raw_data,row.names = 1,header = T,check.names = F)
meta_scale = read.csv(scale_data,row.names = 1,header = T,check.names = F)
select_df = function(meta){
  num_cols = ncol(meta)
  selected_cols = if (num_cols > 25) {
    sample(seq_len(num_cols), 20)
  } else {
    seq_len(num_cols)
  }
  num_rows = nrow(meta)
  selected_rows = if (num_rows > 15) {
    sample(seq_len(num_rows), 15)
  } else {
    seq_len(num_rows)
  }
  return(meta[selected_rows, selected_cols])
}

meta_scale = select_df(meta_scale)
meta = meta[rownames(meta_scale),colnames(meta_scale)]
pdf(plot_dir_1,width=6, height=5)
pheatmap(meta,
         cluster_rows = F,
         cluster_cols = F,
         border = NA,
         legend = T,
         cellwidth = 10, 
         cellheight = 10,
         fontsize = 7)
dev.off()
pdf(plot_dir_2, width=6, height=5)
pheatmap(meta_scale,
         cluster_rows = F,
         cluster_cols = F,
         border = NA,
         legend = T,
         cellwidth = 10, 
         cellheight = 10,
         fontsize = 7)
dev.off()
