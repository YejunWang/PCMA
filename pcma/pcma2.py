import pandas as pd
import os
from pcma.correlation import compute_correlation
from pcma.pca_analysis import run_pca
from pcma.mediation_pcma2 import mediation_pcma2
from pcma.heatmap_pc import heatmap_pc
from pcma.volcano import volcano
from pcma.heatmap_scale import heatmap_scale
from pcma.heatmap_pca import heatmap_pca
from pcma.sankey import sankey
from pcma.boxplot import boxplot_1
from pcma.metabolite_analysis import metabolite_analysis
from pcma.report import report


def pcma2(output_dir: str,
          Bacteria_dir: str,
          Metabolite_dir: str,
          Diagnosis_dir: str,
          SCC_threshold_Metabolite_Diagnosis=0.2,
          SCC_threshold_Bacteria_Diagnosis=0.2,
          p_threshold_Metabolite_Diagnosis=0.05,
          p_threshold_Bacteria_Diagnosis=0.05,
          p_filter_method_Metabolite_Diagnosis='p',
          p_filter_method_Bacteria_Diagnosis='p',
          meta_PC_num=50,
          meta_PC_norm_method='MinMaxScaler',
          meta_if_whiten=True,
          meta_Norm=True,
          meta_seeds=0,
          n_bootstrap=1000,
          is_func_anal=False,
          func_anal_file=""):
    '''PCMA1 workflow'''
    # set the result dir
    result_dir = os.path.join(output_dir, 'result')
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
    print(f'The result file has been estilablished in {result_dir}')

    file_dir = os.path.join(result_dir, 'data')
    if not os.path.exists(file_dir):
        os.makedirs(file_dir)

    plot_dir = os.path.join(result_dir, 'plot')
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    # read data
    Bacteria = pd.read_csv(Bacteria_dir)
    Metabolite = pd.read_csv(Metabolite_dir)
    Diagnosis = pd.read_csv(Diagnosis_dir)

    # decide y type
    if len(pd.unique(Diagnosis.iloc[:, 1:].values.flatten())) == 2:
        y_type = 'Binary Variable'
    if len(pd.unique(Diagnosis.iloc[:, 1:].values.flatten())) > 2 and len(
            pd.DataFrame(Diagnosis.iloc[:, 1:])) < 5:
        y_type = 'Multi-categorical variable'
    if len(pd.unique(Diagnosis.iloc[:, 1:].values.flatten())) >= 5:
        y_type = 'Discrete variable'
    print(f'All data file has been read')

    # filte data with SCC
    Bacteria_rowname_remove = Bacteria.iloc[:, 1:]
    Metabolite_rowname_remove = Metabolite.iloc[:, 1:]
    Diagnosis_rowname_remove = Diagnosis.iloc[:, 1:]

    # Metaboilite-Diagnosis SCC filiter
    Meta_Diag_corr_origin, Meta_Diag_corr = compute_correlation(
        Metabolite_rowname_remove, Diagnosis_rowname_remove,
        SCC_threshold_Metabolite_Diagnosis, p_threshold_Metabolite_Diagnosis,
        p_filter_method_Metabolite_Diagnosis)
    meta_to_keep = set(Meta_Diag_corr['Column_A']).intersection(
        set(Metabolite_rowname_remove.columns))
    Metabolite_filt_final = Metabolite_rowname_remove[list(meta_to_keep)]

    # Bacteria-Diagnosis SCC filter
    Bact_Diag_corr_origin, Bact_Diag_corr = compute_correlation(
        Bacteria_rowname_remove, Diagnosis_rowname_remove,
        SCC_threshold_Bacteria_Diagnosis, p_threshold_Bacteria_Diagnosis,
        p_filter_method_Bacteria_Diagnosis)
    columns_to_keep = set(Bact_Diag_corr['Column_A']).intersection(
        set(Bacteria.columns))
    Bacteria_filt = Bacteria[list(columns_to_keep)]

    # PCA
    meta_scale_data, meta_pca_df, meta_pca_contribution, meta_pca_components = run_pca(
        Metabolite_filt_final,
        PC_num=meta_PC_num,
        norm_method=meta_PC_norm_method,
        if_whiten=meta_if_whiten,
        Norm=meta_Norm,
        seeds=meta_seeds)

    # %%
    # mediation
    result_mediation, coef, p_value = mediation_pcma2(Bacteria_filt, meta_pca_df,
                                             Metabolite.iloc[:, 0],
                                             Diagnosis.iloc[:, 1])

    # serialize data
    meta_scale_data.insert(0, 'Sample_Name', Metabolite.iloc[:, 0])
    meta_scale_data.to_csv(os.path.join(file_dir, 'Metabolite_scale_data.csv'),
                           index=False)
    meta_pca_df.insert(0, 'Sample_Name', Metabolite.iloc[:, 0])
    meta_pca_df.to_csv(os.path.join(file_dir, 'Metabolite_pca.csv'),
                       index=False)
    meta_pca_contribution.to_csv(os.path.join(
        file_dir, 'Metabolite_pca_contribution.csv'),
                                 index=False)
    meta_pca_components.to_csv(os.path.join(file_dir,
                                            'Metabolite_pca_components.csv'),
                               index=False)
    Meta_Diag_corr_origin.to_csv(os.path.join(
        file_dir, 'Metabolite_Diagnosis_correlation_all.csv'),
                                 index=False)
    Meta_Diag_corr.to_csv(os.path.join(
        file_dir, 'Metabolite_Diagnosis_correlation_filtered.csv'),
                          index=False)
    Bact_Diag_corr_origin.to_csv(os.path.join(
        file_dir, 'Bacteria_Diagnosis_correaltion_all.csv'),
                                 index=False)
    Bact_Diag_corr.to_csv(os.path.join(
        file_dir, 'Bacteria_Diagnosis_correaltion_filtered.csv'),
                          index=False)
    Metabolite_filt_final.insert(0, 'Sample_Name', Metabolite.iloc[:, 0])
    Metabolite_filt_final.to_csv(os.path.join(file_dir,
                                              'Metabolite_filtered.csv'),
                                 index=False)
    Bacteria_filt.to_csv(os.path.join(file_dir, 'Bacteria_filter.csv'),
                         index=False)
    result_mediation.to_csv(os.path.join(file_dir, 'Siginficant_PC.csv'),
                            index=False)
    coef.to_csv(os.path.join(file_dir, 'coefficient.csv'), index=False)
    p_value.to_csv(os.path.join(file_dir, 'p_value.csv'), index=False)

    # output all paraments users choose
    params = {
        'method': 'PCMA2',
        'Selection_X_Y': f'SCC; {p_filter_method_Metabolite_Diagnosis}',
        'Selection_X_Y_threshold':
        f"SCC >= {SCC_threshold_Metabolite_Diagnosis}; {p_filter_method_Metabolite_Diagnosis} <= {p_threshold_Metabolite_Diagnosis}",
        'PCA_X': 'NA',
        'PCA_X_heatmap': 'NA',
        'PCA_X_result': 'NA',
        'PCA_X_components': 'NA',
        'PCA_X_EVR': 'NA',
        'PCA_M_result': f"{os.path.join(file_dir,'Metabolite_pca.csv')}",
        'PCA_M_components':
        f"{os.path.join(file_dir,'Metabolite_pca_components.csv')}",
        'PCA_M_EVR':
        f"{os.path.join(file_dir,'Metabolite_pca_contribution.csv')}",
        'output_dir': output_dir,
        'Bacteria_dir': Bacteria_dir,
        'Metabolite_dir': Metabolite_dir,
        'Diagnosis_dir': Diagnosis_dir,
        'SCC_threshold_Metabolite_Diagnosis':
        SCC_threshold_Metabolite_Diagnosis,
        'SCC_threshold_Bacteria_Diagnosis': SCC_threshold_Bacteria_Diagnosis,
        'SCC_threshold_Bacteria_Metabolite': 'NA',
        'p_threshold_Metabolite_Diagnosis': p_threshold_Metabolite_Diagnosis,
        'p_threshold_Bacteria_Diagnosis': p_threshold_Bacteria_Diagnosis,
        'p_threshold_Bacteria_Metabolite': 'NA',
        'p_filter_method_Metabolite_Diagnosis':
        p_filter_method_Metabolite_Diagnosis,
        'p_filter_method_Bacteria_Diagnosis':
        p_filter_method_Bacteria_Diagnosis,
        'p_filter_method_Bacteria_Metabolite': 'NA',
        'meta_PC_num': meta_PC_num,
        'meta_PC_norm_method': meta_PC_norm_method,
        'meta_if_whiten': meta_if_whiten,
        'meta_Norm': meta_Norm,
        'meta_seeds': meta_seeds,
        'n_bootstrap': n_bootstrap,
        'is_function_analysis': is_func_anal,
        "function_analysis_file": func_anal_file
    }
    params_df = pd.DataFrame(list(params.items()),
                             columns=['Parameter', 'Value'])
    params_df.to_csv(os.path.join(file_dir, 'paramters_list.csv'), index=False)
    print('paprmeters has been stored.')

    # plot
    # heatmap
    result = heatmap_pc(result_mediation, file_dir, plot_dir, method='pcma2')
    # volcano
    # Bact_Diag
    volcano(resource_data=os.path.join(
        file_dir, 'Bacteria_Diagnosis_correaltion_all.csv'),
            plot_dir=os.path.join(
                plot_dir, 'Volcano_plot_Bacteria_Diagnosis_correlation.pdf'),
            y_value=p_filter_method_Bacteria_Diagnosis,
            y_threshold=str(p_threshold_Bacteria_Diagnosis),
            SCC_threshold=str(SCC_threshold_Bacteria_Diagnosis))
    # Meta_Diag
    volcano(resource_data=os.path.join(
        file_dir, 'Metabolite_Diagnosis_correlation_all.csv'),
            plot_dir=os.path.join(
                plot_dir, 'Volcano_plot_Metabolite_Diagnosis_correlation.pdf'),
            y_value=p_filter_method_Metabolite_Diagnosis,
            y_threshold=str(p_threshold_Metabolite_Diagnosis),
            SCC_threshold=str(SCC_threshold_Metabolite_Diagnosis))
    # heatmap_scale
    heatmap_scale(Metabolite_dir,
                  os.path.join(file_dir, 'Metabolite_scale_data.csv'),
                  os.path.join(plot_dir, 'heatmap_metabolite_raw.pdf'),
                  os.path.join(plot_dir, 'heatmap_metabolite_scale.pdf'))
    # heatmap:pca_M
    heatmap_pca(os.path.join(file_dir, 'Metabolite_pca_components.csv'),
                os.path.join(plot_dir, 'PCA_M_heatmap.pdf'))
    # sankey plot
    sankey(Bacteria_dir, os.path.join(file_dir, 'Metabolite_pca.csv'),
           Diagnosis_dir, os.path.join(file_dir, 'Siginficant_PC.csv'),
           os.path.join(file_dir, 'sankey_data_dir.csv'),
           os.path.join(plot_dir, 'sankey_plot.png'))
    # boxplot
    boxplot_1(file_dir, Bacteria_dir, Metabolite_dir, Diagnosis_dir, plot_dir,
              3)
    # functional analysis
    if is_func_anal == False:
        pass
    if is_func_anal == True:
        metabolite_analysis(
            os.path.join(file_dir, 'Metabolite_pca_components.csv'),
            func_anal_file, os.path.join(file_dir, 'Siginficant_PC.csv'),
            os.path.join(plot_dir, 'barplot_fisher.pdf'), file_dir)
    # report
    report(result_dir, file_dir, plot_dir)
