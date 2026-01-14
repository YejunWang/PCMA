if (!require(tidyverse, quietly = TRUE)) {
  install.packages("tidyverse",repos='https://cloud.r-project.org/',dependencies=TRUE)
  library(tidyverse)
}

args = commandArgs(trailingOnly = TRUE)
Bacteria_dir = args[1]
Metabolite_pca_dir = args[2]
Diag_dir = args[3]
Sig_dir = args[4]
output_data_dir = args[5]

Bacteria = read.csv(Bacteria_dir, check.names=F)
Metabolite = read.csv(Metabolite_pca_dir, check.names=F)
Diag = read.csv(Diag_dir, check.names=F)
Sig = read.csv(Sig_dir, check.names=F)

Sig$Significant_PC_bact = paste(Sig$Significant_PC, Sig$Bacteria_PC, sep = "_")
links = data.frame(source = character(), target = character(), Spearman_Correlation = numeric(), stringsAsFactors = FALSE)

for (i in 1:nrow(Sig)) {
  Sig_pc = Sig$Significant_PC[i]
  bacteria = Sig$Bacteria_PC[i]
  Sig_pc_bact = Sig$Significant_PC_bact[i]
  
  bacteria_values = Bacteria[bacteria]
  metabolite_values = Metabolite[Sig_pc]
  correlation_bacteria_metabolite = c(cor(bacteria_values, metabolite_values, method = "spearman"))
  
  diagnosis_values = Diag[2]
  correlation_metabolite_diagnosis = c(cor(metabolite_values, diagnosis_values, method = "spearman"))
  
  links = rbind(links, data.frame(source = bacteria, target = Sig_pc_bact, Spearman_Correlation = correlation_bacteria_metabolite))
  links = rbind(links, data.frame(source = Sig_pc_bact, target = "Disease Phenotype", Spearman_Correlation = correlation_metabolite_diagnosis))
}
write.csv(links,output_data_dir)