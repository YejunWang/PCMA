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

meta_PC_res = args[2]
bact_PC_res = args[3]

bacteria = args[4]
metabolite = args[5]
diagnosis = args[6]

meta_PC_components = args[7]
bact_PC_components = args[8]

output_dir = args[9]

designate_pathway_num = as.numeric(args[10])

# read data
pathway = read.csv(pathway)

meta_PC_res = read.csv(meta_PC_res,check.names = F)
bact_PC_res = read.csv(bact_PC_res,check.names = F)

bacteria = read.csv(bacteria,check.names = F)
metabolite = read.csv(metabolite,check.names = F)
diagnosis = read.csv(diagnosis,check.names = F)

meta_PC_components = read.csv(meta_PC_components,check.names = F)
bact_PC_components = read.csv(bact_PC_components,check.names = F)


# fix rownames and data
rownames(meta_PC_components) = paste0("PC",1:nrow(meta_PC_components))
rownames(bact_PC_components) = paste0("PC",1:nrow(bact_PC_components))

designate_pathway_num = min(nrow(pathway),designate_pathway_num)

pathway['Diagnosis'] = 'Diagnosis'

setwd(output_dir)

# sampling
if (designate_pathway_num > nrow(pathway)) {
  stop("Sample size exceeds the total number of rows in the data.")
}
pathway_sample = pathway[sample(nrow(pathway), designate_pathway_num), ]


# Bacteria_PC compare
bact_PC_res = merge(bact_PC_res,diagnosis,by.x=colnames(bact_PC_res)[1],by.y=colnames(diagnosis)[1])
plot_bact_pc = list()
for(i in 1:nrow(pathway_sample)){
  bact = pathway_sample[i,'Bacteria_PC']
  meta = pathway_sample[i,'Significant_PC']
  p = ggplot(bact_PC_res, aes(x = factor(bact_PC_res[,ncol(bact_PC_res)]), y = .data[[bact]])) +
    geom_boxplot() +
    labs(x = "Diagnosis", y = "PC") +
    ggtitle(paste('Pathway:',bact,'-',meta,'-','Diagnosis')) +
    theme_bw() +
    stat_compare_means(method = "wilcox.test", 
                       label = "p.format", 
                       label.y = max(bact_PC_res[[bact]]) * 0.9,
                       comparisons = list(as.character(unique(bact_PC_res[,ncol(bact_PC_res)])))) +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5,face = "bold"),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10)
    )
  plot_bact_pc[[i]] = p
}

pdf(file = "boxplot_bacteria_pc.pdf", width = 12, height = 8)
do.call(grid.arrange, c(plot_bact_pc, ncol = 3))
dev.off()

# metabolite_PC compare
meta_PC_res = merge(meta_PC_res,diagnosis,by.x = colnames(meta_PC_res)[1], by.y = colnames(diagnosis)[1])
plot_meta_pc = list()

for(i in 1:nrow(pathway_sample)){
  bact = pathway_sample[i,'Bacteria_PC']
  meta = pathway_sample[i,'Significant_PC']
  p = ggplot(meta_PC_res, aes(x = factor(.data[[colnames(meta_PC_res)[ncol(meta_PC_res)]]]), y = .data[[meta]])) +
    geom_boxplot() +
    labs(x = "Diagnosis", y = "PC") +
    ggtitle(paste('Pathway:',bact,'-',meta,'-','Diagnosis')) +
    theme_bw() +
    stat_compare_means(method = "wilcox.test", 
                       label = "p.format", 
                       label.y = max(meta_PC_res[[meta]]) * 0.9,
                       comparisons = list(as.character(unique(meta_PC_res[,ncol(meta_PC_res)])))) +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5,face = "bold"),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10)
    )
  plot_meta_pc[[i]] = p
}
pdf(file = "boxplot_metabolite_pc.pdf", width = 12, height = 8)
do.call(grid.arrange, c(plot_meta_pc, ncol = 3))
dev.off()

# metabolite
metabolite = merge(metabolite,diagnosis,by.x = colnames(metabolite)[1], by.y = colnames(diagnosis)[1])
plot_meta = list()
for(i in 1:nrow(pathway_sample)){
  bact = pathway_sample[i,'Bacteria_PC']
  meta = pathway_sample[i,'Significant_PC']
  
  meta_real = colnames(meta_PC_components)[which.max(abs(meta_PC_components[meta,]))] # get the significant metabolite
  
  p = ggplot(metabolite, aes(x = factor(metabolite[,ncol(metabolite)]), y = .data[[rlang::sym(meta_real)]])) +
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
pdf(file = "boxplot_metabolite.pdf", width = 12, height = 8)
do.call(grid.arrange, c(plot_meta, ncol = 3))
dev.off()

# bacteria
bacteria = merge(bacteria,diagnosis,by.x = colnames(bacteria)[1], by.y = colnames(diagnosis)[1])
plot_bact = list()
for(i in 1:nrow(pathway_sample)){
  bact = pathway_sample[i,'Bacteria_PC']
  meta = pathway_sample[i,'Significant_PC']
  
  bact_real = colnames(bact_PC_components)[which.max(abs(bact_PC_components[bact,]))] # get the significant bacteria
  
  p = ggplot(bacteria, aes(x = factor(bacteria[,ncol(bacteria)]), y = .data[[bact_real]])) +
    geom_boxplot() +
    labs(x = "Diagnosis", y = `bact_real`) +
    ggtitle(paste('Pathway:',bact,'-',meta,'-','Diagnosis')) +
    theme_bw() +
    stat_compare_means(method = "wilcox.test", 
                       label = "p.format", 
                       label.y = max(bacteria[[`bact_real`]]) * 0.9,
                       comparisons = list(as.character(unique(bacteria[,ncol(bacteria)])))) +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5,face = "bold"),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10)
    )
  plot_meta[[i]] = p
}
pdf(file = "boxplot_bacteria.pdf", width = 12, height = 8)
do.call(grid.arrange, c(plot_meta, ncol = 3))
dev.off()