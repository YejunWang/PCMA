# Volcano plot
if (!require(ggplot2, quietly = TRUE)) {
  install.packages("ggplot2",repos='https://cloud.r-project.org/',dependencies=TRUE)
  library(ggplot2)
}
args = commandArgs(trailingOnly = TRUE)
resource_data = args[1]
plot_dir = args[2]
y_value = args[3]
y_threshold = as.numeric(args[4])
SCC_threshold = as.numeric(args[5])

volcano = read.csv(resource_data)

# compute fdr
volcano['FDR'] = p.adjust(volcano$P_value, method = "fdr")

if(y_value == "p"){
  volcano$threshold="normal"
  volcano$threshold[volcano$P_value<y_threshold & (volcano$Correlation >= SCC_threshold)] = "up"
  volcano$threshold[volcano$P_value<y_threshold & (volcano$Correlation <= -SCC_threshold)] = "down"
  volcano$threshold=factor(volcano$threshold, levels=c("up","down","normal"), order=T)
  p = ggplot(volcano,aes(x=Correlation,y=-log10(P_value),color=threshold)) + 
    theme_bw() + 
    labs(x="Correlation",y="-log10(P_value)",title="Correlation Analysis Volcano Plot")
}else{
  if(y_value == 'fdr'){
    volcano$threshold="normal"
    volcano$threshold[volcano$FDR<y_threshold & (volcano$Correlation >= SCC_threshold)] = "up"
    volcano$threshold[volcano$FDR<y_threshold & (volcano$Correlation <= -SCC_threshold)] = "down"
    volcano$threshold=factor(volcano$threshold, levels=c("up","down","normal"), order=T)
    p = ggplot(volcano,aes(x=Correlation,y=-log10(FDR),color=threshold)) + 
    labs(x="Correlation",y="-log10(FDR)",title="Correlation Analysis Volcano Plot")
  }
}

pdf(plot_dir,width=6, height=4)
p = p + geom_point(alpha = 0.3) + 
  scale_color_manual(values=c("#DC143C","#00008B","#808080")) + 
  geom_hline(yintercept = -log10(y_threshold),lty=3,col="black",lwd=0.7) + 
  geom_vline(xintercept = SCC_threshold,lty=3,col="black",lwd=0.7) + 
  geom_vline(xintercept = -SCC_threshold,lty=3,col="black",lwd=0.7) + 
  theme(plot.title = element_text(size=10, color="black", hjust = 0.5),
        legend.position="none", 
        legend.title = element_blank(),legend.text= element_text(size=12, color="black"),
        axis.text=element_text(size=7, color="black"),
        axis.title=element_text(size=10, color="black", vjust=0.5, hjust=0.5))
p
dev.off()