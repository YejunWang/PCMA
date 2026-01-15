#%%
import os
import pkg_resources
import pandas as pd
from jinja2 import Environment, FileSystemLoader
from weasyprint import HTML
from pdf2image import convert_from_path
from PyPDF2 import PdfReader, PdfWriter


#%%
def report(result_dir, file_dir, plot_dir):
    template_path = pkg_resources.resource_filename(
        'pcma', 'report_templete/templete.html')
    parms = pd.read_csv(os.path.join(file_dir, 'paramters_list.csv'))
    output_path = os.path.join(result_dir, "Report.pdf")

    def render_pdf(template_path, output_path, context):
        template_dir = os.path.dirname(template_path)
        template_file = os.path.basename(template_path)
        env = Environment(loader=FileSystemLoader(template_dir))
        template = env.get_template(template_file)
        html_content = template.render(context)
        HTML(string=html_content).write_pdf(output_path)

    def pdf2png(dir, imgName):
        # deal with pdf
        pdf_dir = os.path.join(dir, f'{imgName}.pdf')
        images = convert_from_path(pdf_dir, dpi=600)
        image_dir = os.path.join(dir, f'{imgName}.png')
        images[0].save(image_dir, 'PNG')
        return image_dir

    def get_meta_filter_result(dir):
        meta_filt = pd.read_csv(os.path.join(dir, 'Metabolite_filtered.csv'))
        meta = pd.read_csv(parms.loc[parms['Parameter'] == 'Metabolite_dir',
                                     'Value'].values[0])
        return (len(meta.columns) - 1), (len(meta_filt.columns) - 1)

    main_method = f"{parms.loc[parms['Parameter'] == 'method', 'Value'].values[0]}"

    context = {
        # 1.1
        'main_method':
        main_method,
        # 1.2
        'Input_files':
        [[
            'X(Independent Variable)',
            f"{parms.loc[parms['Parameter'] == 'Bacteria_dir', 'Value'].values[0]}",
            "Discrete variable"
        ],
         [
             'M(mediator)',
             f"{parms.loc[parms['Parameter'] == 'Metabolite_dir', 'Value'].values[0]}",
             "Discrete variable"
         ],
         [
             'Y(Dependent Variable)',
             f"{parms.loc[parms['Parameter'] == 'Diagnosis_dir', 'Value'].values[0]}",
             f"{parms.loc[parms['Parameter'] == 'y_type', 'Value'].values[0]}"
         ]],
        # 1.3
        'feature_selection':
        [[
            "X~Y",
            f"{parms.loc[parms['Parameter'] == 'Selection_X_Y', 'Value'].values[0]}",
            f"{parms.loc[parms['Parameter'] == 'Selection_X_Y_threshold', 'Value'].values[0]}"
        ],
         [
             "M~Y",
             f"SCC; {parms.loc[parms['Parameter'] == 'p_filter_method_Metabolite_Diagnosis', 'Value'].values[0]}",
             f"SCC >= {parms.loc[parms['Parameter'] == 'SCC_threshold_Metabolite_Diagnosis', 'Value'].values[0]}; {parms.loc[parms['Parameter'] == 'p_filter_method_Metabolite_Diagnosis', 'Value'].values[0]} <= {parms.loc[parms['Parameter'] == 'p_threshold_Metabolite_Diagnosis', 'Value'].values[0]}"
         ],
         [
             "X~M",
             f"SCC; {parms.loc[parms['Parameter'] == 'p_filter_method_Bacteria_Metabolite', 'Value'].values[0]}",
             f"SCC >= {parms.loc[parms['Parameter'] == 'SCC_threshold_Bacteria_Metabolite', 'Value'].values[0]}; {parms.loc[parms['Parameter'] == 'p_filter_method_Bacteria_Metabolite', 'Value'].values[0]} <= {parms.loc[parms['Parameter'] == 'p_threshold_Bacteria_Metabolite', 'Value'].values[0]}"
         ]],
        # 1.4
        'normalization':
        [["X", "NA"],
         [
             "M",
             f"{parms.loc[parms['Parameter'] == 'meta_PC_norm_method', 'Value'].values[0]}"
         ]],
        # 2
        # 2.1
        'meta_all_num':
        f"{get_meta_filter_result(file_dir)[0]}",
        'meta_filt_num':
        f"{get_meta_filter_result(file_dir)[1]}",
        'volcano_diagram':
        [['X~Y', 'NA', 'NA'],
         [
             'M~Y', 'Volcano',
             f"{os.path.join(file_dir,'Metabolite_Diagnosis.csv')}"
         ],
         [
             'X~M', 'Volcano',
             f"{os.path.join(file_dir,'Bacteria_Metabolite.csv')}"
         ]],
        'volcano_plot_1':
        "",
        'volcano_plot_1_caption':
        "",
        'volcano_plot_2':
        "",
        'volcano_plot_2_caption':
        "",

        # 2.2
        'heatmap_diagram':
        [['X', 'NA', 'NA'],
         [
             'M', 'Heatmap',
             f"{os.path.join(file_dir, 'Metabolite_scale_data.csv')}"
         ]],
        'heatmap_plot_1':
        f"file:///mnt/md0/Public/PCMA_Web_sever/PCMA/{pdf2png(plot_dir,'heatmap_metabolite_raw')}",
        'heatmap_plot_1_caption':
        "Heatmap (Before Normalization)",
        'heatmap_plot_2':
        f"file:///mnt/md0/Public/PCMA_Web_sever/PCMA/{pdf2png(plot_dir,'heatmap_metabolite_scale')}",
        'heatmap_plot_2_caption':
        "Heatmap (After Normalization)",
        # 3
        # 3.1：X：
        'main_method_all':
        f"{parms.loc[parms['Parameter'] == 'method', 'Value'].values[0]}",
        'PCA_X_diagram':
        [[
            "Result (Values of each sample's PC)",
            "It shows the values of each principal component (PC) in each sample after the independent variable X has been processed by PCA.",
            f"{parms.loc[parms['Parameter'] == 'PCA_X_result', 'Value'].values[0]}"
        ],
         [
             "Components",
             "It shows the principal components (PCs) derived from the independent variable X after PCA processing, with each PC linearly composed of the original variables with certain coefficients.",
             f"{parms.loc[parms['Parameter'] == 'PCA_X_components', 'Value'].values[0]}"
         ],
         [
             "Explained Variance Ratio",
             "-",
             f"{parms.loc[parms['Parameter'] == 'PCA_X_EVR', 'Value'].values[0]}",
         ]],
        'PCA_X_heatmap':
        f"{'file:///mnt/md0/Public/PCMA_Web_sever/PCMA/' + pdf2png(plot_dir, 'PCA_X_heatmap') if parms.loc[parms['Parameter'] == 'PCA_X_heatmap', 'Value'].values[0] == 'True' else 'NA'}",
        'PCA_X_heatmap_caption':
        'PCA Bacteria Heatmap (PC1 to PC15)',

        # 3.2
        'PCA_M_heatmap':
        f"file:///mnt/md0/Public/PCMA_Web_sever/PCMA/{pdf2png(plot_dir,'PCA_M_heatmap')}",
        'PCA_M_diagram':
        [[
            "Result (Values of each sample's PC)",
            "It shows the values of each principal component (PC) in each sample after the mediator variable M has been processed by PCA",
            f"{parms.loc[parms['Parameter'] == 'PCA_M_result', 'Value'].values[0]}"
        ],
         [
             "Components",
             "It shows the principal components (PCs) derived from the mediator variable M after PCA processing, with each PC linearly composed of the original variables with certain coefficients.",
             f"{parms.loc[parms['Parameter'] == 'PCA_M_components', 'Value'].values[0]}"
         ],
         [
             "Explained Variance Ratio", "-",
             f"{parms.loc[parms['Parameter'] =='PCA_M_EVR', 'Value'].values[0]}"
         ]],
        'PCA_M_heatmap_caption':
        'PCA Metabolite Heatmap (PC1 to PC15)',
        # 4
        'heatmap_sig_PC':
        f"file:///mnt/md0/Public/PCMA_Web_sever/PCMA/{pdf2png(plot_dir,'heatmap')}",
        'heatmap_sig_PC_caption':
        'Significant PC Plot. X-axis: bacteria (PC), Y-axis: metabolites (PC)',
        'significant_results_path':
        f"{os.path.join(file_dir,'Siginficant_PC.csv')}",
        'coff_mediation_path':
        f"{os.path.join(file_dir,'coefficient.csv')}",
        'sankey_Plot':
        f"file:///mnt/md0/Public/PCMA_Web_sever/PCMA/{os.path.join(plot_dir,'sankey_plot.png')}",
        'sanekey_Plot_caption':
        'Sankey Plot. From left to right, the three layers are Bacteria (PC), Metabolites (PC), and Disease Phenotype.',
        # 5
        'boxplot_bacteria':
        f"file:///mnt/md0/Public/PCMA_Web_sever/PCMA/{pdf2png(plot_dir,'boxplot_bacteria')}",
        'boxplot_bacteria_caption':
        '',
        'boxplot_metabolite':
        f"file:///mnt/md0/Public/PCMA_Web_sever/PCMA/{pdf2png(plot_dir,'boxplot_metabolite')}",
        'boxplot_metabolite_caption':
        'Boxplot for metabolite in Significant Metabolite PC.',
        'boxplot_metabolite_pc':
        f"file:///mnt/md0/Public/PCMA_Web_sever/PCMA/{pdf2png(plot_dir,'boxplot_metabolite_pc')}",
        'boxplot_metabolite_pc_caption':
        'Boxplot for Metabolite Significant PC'
    }
    if main_method == 'PCMA1':
        context[
            'volcano_plot_1'] = f"file:///mnt/md0/Public/PCMA_Web_sever/PCMA/{pdf2png(plot_dir,'Volcano_plot_Bacteria_Metabolite_correlation')}"
        context[
            'volcano_plot_1_caption'] = "Bacteria-Metabolite Correlation Analysis Volcano Plot"
        context[
            'volcano_plot_2'] = f"file:///mnt/md0/Public/PCMA_Web_sever/PCMA/{pdf2png(plot_dir,'Volcano_plot_Metabolite_Diagnosis_correlation')}"
        context[
            'volcano_plot_2_caption'] = "Metabolite-Diagnosis Correlation Analysis Volcano Plot"
        context[
            'boxplot_bacteria_caption'] = 'Boxplot for Bacteria Relative abundance.'
    elif main_method == 'PCMA2':
        context[
            'boxplot_bacteria_caption'] = 'Boxplot for Bacteria Relative Abundance'
        context['feature_selection'] = [
            [
                "X~Y",
                f"{parms.loc[parms['Parameter'] == 'Selection_X_Y', 'Value'].values[0]}",
                f"{parms.loc[parms['Parameter'] == 'Selection_X_Y_threshold', 'Value'].values[0]}"
            ],
            [
                "M~Y",
                f"SCC; {parms.loc[parms['Parameter'] == 'p_filter_method_Metabolite_Diagnosis', 'Value'].values[0]}",
                f"SCC >= {parms.loc[parms['Parameter'] == 'SCC_threshold_Metabolite_Diagnosis', 'Value'].values[0]}; {parms.loc[parms['Parameter'] == 'p_filter_method_Metabolite_Diagnosis', 'Value'].values[0]} <= {parms.loc[parms['Parameter'] == 'p_threshold_Metabolite_Diagnosis', 'Value'].values[0]}"
            ], ["X~M", "nan", "nan"]
        ]
    elif main_method == 'PCPCMA':
        context[
            'volcano_plot_1'] = f"file:///mnt/md0/Public/PCMA_Web_sever/PCMA/{pdf2png(plot_dir,'Volcano_plot_Bacteria_Diagnosis_correlation')}"
        context[
            'volcano_plot_1_caption'] = "Bacteria-Metabolite Correlation Analysis Volcano Plot"
        context[
            'volcano_plot_2'] = f"file:///mnt/md0/Public/PCMA_Web_sever/PCMA/{pdf2png(plot_dir,'Volcano_plot_Metabolite_Diagnosis_correlation')}"
        context[
            'volcano_plot_2_caption'] = "Metabolite-Diagnosis Correlation Analysis Volcano Plot"
        context[
            'boxplot_bacteria_caption'] = 'Boxplot for Bacteria Relative abundance in Significant Bacteria PC'
        context[
            'boxplot_bacteria_pc'] = f"file:///mnt/md0/Public/PCMA_Web_sever/PCMA/{pdf2png(plot_dir,'boxplot_bacteria_pc')}"
        context[
            'boxplot_bacteria_pc_caption'] = 'Boxplot for Bacteria Significant PC'
        context['feature_selection'] = [
            [
                "X~Y",
                f"{parms.loc[parms['Parameter'] == 'Selection_X_Y', 'Value'].values[0]}",
                f"{parms.loc[parms['Parameter'] == 'Selection_X_Y_threshold', 'Value'].values[0]}"
            ],
            [
                "M~Y",
                f"SCC; {parms.loc[parms['Parameter'] == 'p_filter_method_Metabolite_Diagnosis', 'Value'].values[0]}",
                f"SCC >= {parms.loc[parms['Parameter'] == 'SCC_threshold_Metabolite_Diagnosis', 'Value'].values[0]}; {parms.loc[parms['Parameter'] == 'p_filter_method_Metabolite_Diagnosis', 'Value'].values[0]} <= {parms.loc[parms['Parameter'] == 'p_threshold_Metabolite_Diagnosis', 'Value'].values[0]}"
            ], ["X~M", "nan", "nan"]
        ]
    render_pdf(template_path, output_path, context)
    if parms.loc[parms['Parameter'] == 'is_function_analysis',
                 'Value'].values[0] == 'True':
        merge_pdfs(output_path, os.path.join(plot_dir, 'barplot_fisher.pdf'),
                   os.path.join(result_dir, 'Report_pcma.pdf'))


def merge_pdfs(pdf1_path, pdf2_path, output_path):
    pdf_writer = PdfWriter()

    def add_pdf_to_writer(pdf_path):
        pdf_reader = PdfReader(pdf_path)
        for page_num in range(len(pdf_reader.pages)):
            pdf_writer.add_page(pdf_reader.pages[page_num])

    add_pdf_to_writer(pdf1_path)
    add_pdf_to_writer(pdf2_path)

    with open(output_path, 'wb') as output_pdf:
        pdf_writer.write(output_pdf)

    os.remove(pdf1_path)
    os.remove(pdf2_path)


# %%
