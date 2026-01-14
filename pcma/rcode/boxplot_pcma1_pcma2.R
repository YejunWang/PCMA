# boxplot
if (!require(ggplot2, quietly = TRUE)) {
  install.packages("ggplot2",repos='https://cloud.r-project.org/',dependencies=TRUE)
  library(ggplot2)
}

if (!require(dplyr, quietly = TRUE)) {
  install.packages("dplyr", repos = 'https://cloud.r-project.org/', dependencies = TRUE)
  library(dplyr)
}
if (!require(ggpubr, quietly = TRUE)) {
  install.packages("ggpubr", repos = 'https://cloud.r-project.org/', dependencies = TRUE)
  library(ggpubr)
}
if (!require(gridExtra, quietly = TRUE)) {
  install.packages("gridExtra", repos = 'https://cloud.r-project.org/', dependencies = TRUE)
  library(gridExtra)
}

# read bacteria 
args = commandArgs(trailingOnly = TRUE)


pathway = args[1]
PC_res = args[2]

bacteria = args[3]
metabolite = args[4]
diagnosis = args[5]

PC_components = args[6]

output_dir = args[7]
designate_pathway_num = as.numeric(args[8])


pathway = read.csv(pathway)
PC_res = read.csv(PC_res)
bacteria = read.csv(bacteria,check.names = F)
metabolite = read.csv(metabolite,check.names = F)
diagnosis = read.csv(diagnosis,check.names = F)
PC_components = read.csv(PC_components,check.names = F)

setwd(output_dir)

# fix
rownames(PC_components) = paste0("PC",1:nrow(PC_components))
designate_pathway_num = min(nrow(pathway),designate_pathway_num)

pathway['Diagnosis'] = 'Diagnosis'

pathway_sample = pathway[sample(nrow(pathway), designate_pathway_num), ]


# Bacteria compare
bacteria = merge(bacteria,diagnosis,by.x=colnames(bacteria)[1],by.y=colnames(diagnosis)[1])
plot_bact = list()
for(i in 1:nrow(pathway_sample)){
  bact = pathway_sample[i,'Bacteria_PC']
  meta = pathway_sample[i,'Significant_PC']
  p = ggplot(bacteria, aes_string(x = factor(bacteria[,ncol(bacteria)]), y = bact)) +
    geom_boxplot() +
    labs(x = "Diagnosis", y = "Relative abundance") +
    ggtitle(paste('Pathway:',bact,'-',meta,'-','Diagnosis')) +
    theme_bw() +
    stat_compare_means(method = "wilcox.test", 
                       label = "p.format", 
                       label.y = max(bacteria[[bact]]) * 0.9,
                       comparisons = list(as.character(unique(bacteria[,ncol(bacteria)])))) +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5,face = "bold"),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10)
    )
  plot_bact[[i]] = p
}



pdf(file = "boxplot_bacteria.pdf", width = 16, height = 12)
do.call(grid.arrange, c(plot_bact, ncol = 3))
dev.off()

# PC compare
PC_res = merge(PC_res,diagnosis,by.x = colnames(PC_res)[1], by.y = colnames(diagnosis)[1])
plot_pc = list()
for(i in 1:nrow(pathway_sample)){
  bact = pathway_sample[i,'Bacteria_PC']
  meta = pathway_sample[i,'Significant_PC']
  p = ggplot(PC_res, aes_string(x = factor(PC_res[,ncol(PC_res)]), y = meta)) +
    geom_boxplot() +
    labs(x = "Diagnosis", y = "PC") +
    ggtitle(paste('Pathway:',bact,'-',meta,'-','Diagnosis')) +
    theme_bw() +
    stat_compare_means(method = "wilcox.test", 
                       label = "p.format", 
                       label.y = max(PC_res[[meta]]) * 0.9,
                       comparisons = list(as.character(unique(PC_res[,ncol(PC_res)])))) +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5,face = "bold"),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10)
    )
  plot_pc[[i]] = p
}


pdf(file = "boxplot_metabolite_pc.pdf", width = 16, height = 12)
do.call(grid.arrange, c(plot_pc, ncol = 3))
dev.off()

# metabolite
metabolite = merge(metabolite,diagnosis,by.x = colnames(metabolite)[1], by.y = colnames(diagnosis)[1])
plot_meta = list()
for(i in 1:nrow(pathway_sample)){
  bact = pathway_sample[i,'Bacteria_PC']
  meta = pathway_sample[i,'Significant_PC']
  
  meta_real = colnames(PC_components)[which.max(abs(PC_components[meta,]))] # get the significant metabolite
  
  p = ggplot(metabolite, aes(x = factor(metabolite[,ncol(metabolite)]), y = .data[[meta_real]])) +
    geom_boxplot() +
    labs(x = "Diagnosis", y = `meta_real`) +
    ggtitle(paste('Pathway:',bact,'-',meta,'-','Diagnosis')) +
    theme_bw() +
    stat_compare_means(method = "wilcox.test", 
                       label = "p.format", 
                       label.y = max(metabolite[[`meta_real`]]) * 0.9,
                       comparisons = list(as.character(unique(metabolite[,ncol(metabolite)])))) +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5,face = "bold"),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10)
    )
  plot_meta[[i]] = p
}


pdf(file = "boxplot_metabolite.pdf", width = 16, height = 12)
do.call(grid.arrange, c(plot_meta, ncol = 3))
dev.off()