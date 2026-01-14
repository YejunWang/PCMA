# PCMA

Mediation analysis is an important approach that has promising applications in pathway analysis in microbiota studies. However, thousands of intermediate metabolites and various interactions among them make it a big challenge to identify the real relevant pathways with traditional mediation analysis. In this study, we proposed a new method, PCMA, which couples Principal Component Analysis with mediation analysis and can avoid the interference caused by various noises and cofounding factors. With multiple testing datasets, PCMA could always effectively identify extensively simplified significant pathways connecting microbiota and diseases, which are mediated by Principal Components of metabolites. Extensive annotation or enrichment analysis on the metabolites contributing most to the significant Principal Components could provide functional insights on the connection pathways. Furthermore, application of PCMA disclosed a few common metabolite pathways, each connecting multiple bacterial taxa and inflammatory bowel diseases. PCMA also revealed a specific pathway that connects both Enterococcus faecium and Lactobacillus gasseri, and Crohn's disease. In conclusion, the study proposed a new scheme of mediation analysis, which could be effectively applied in identifying relevant connection pathways and providing insights into the molecular mechanisms in microbiota studies. A webserver to implement PCMA on line or the standalone version of PCMA is accessible through the link: http://www.szu-bioinf.org/PCMA/.

# Install

To install the latest version of PCMA using `uv pip`.
```bash
uv pip install git+https://github.com/YejunWang/PCMA.git
```

# cite us

Weimiao Kong, Yuqi Lin, and Yejun Wang. 2024. PCMA: A New Mediation Analysis Method Effectively Identifies Relevant Pathways Connecting Microbiota and Diseases. In Proceedings of the 15th ACM International Conference on Bioinformatics, Computational Biology and Health Informatics (BCB '24). Association for Computing Machinery, New York, NY, USA, Article 79, 1. https://doi.org/10.1145/3698587.3701428
        
        
        
        

