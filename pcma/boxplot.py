import os
import subprocess
from pkg_resources import resource_filename
'''
Boxplot for pathway in PCMA1/PCMA2/PCPCMA
'''


def boxplot_1(file_dir,
              Bacteria_dir,
              Metabolite_dir,
              Diagnosis_dir,
              output_dir,
              designate_pathway_num=3):

    boxlplot = resource_filename('pcma', 'rcode/boxplot_pcma1_pcma2.R')

    pathway = os.path.join(file_dir, "Siginficant_PC.csv")
    PC_res = os.path.join(file_dir, "Metabolite_pca.csv")

    PC_components = os.path.join(file_dir, "Metabolite_pca_components.csv")

    result = subprocess.run([
        "Rscript", boxlplot, pathway, PC_res, Bacteria_dir, Metabolite_dir,
        Diagnosis_dir, PC_components, output_dir,
        str(designate_pathway_num)
    ],
                            capture_output=True,
                            text=True)
    print(result.stdout)


def boxplot_2(file_dir,
              Bacteria_dir,
              Metabolite_dir,
              Diagnosis_dir,
              output_dir,
              designate_pathway_num=3):

    boxlplot = resource_filename('pcma', 'rcode/boxplot_pcpcma.R')

    pathway = os.path.join(file_dir, "Siginficant_PC.csv")
    meta_PC_res = os.path.join(file_dir, "Metabolite_pca.csv")
    bact_PC_res = os.path.join(file_dir, "Bacteria_pca.csv")

    meta_PC_components = os.path.join(file_dir,
                                      "Metabolite_pca_components.csv")
    bact_PC_components = os.path.join(file_dir, "Bacteria_pca_components.csv")

    result = subprocess.run([
        "Rscript", boxlplot, pathway, meta_PC_res, bact_PC_res, Bacteria_dir,
        Metabolite_dir, Diagnosis_dir, meta_PC_components, bact_PC_components,
        output_dir,
        str(designate_pathway_num)
    ],
                            capture_output=True,
                            text=True)
    print(result.stdout)
