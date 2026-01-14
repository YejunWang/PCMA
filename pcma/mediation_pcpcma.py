import pandas as pd
import numpy as np
import statsmodels.api as sm
from .bootstrap_analysis import bootstrap_p_values


def mediation_PCPCMA(meta_pca_df: pd.DataFrame,
                     bact_pca_df: pd.DataFrame,
                     Sample_Name: pd.Series,
                     Diagnoisis: pd.Series,
                     n_bootstrap=1000):
    # running PCPCMA
    media_df_origin = meta_pca_df
    media_df_origin.insert(0, 'Sample_Name', Sample_Name)
    media_df_origin.insert(1, 'Diagnoisis', Diagnoisis)
    result_all = {}
    result_all_bact_coef = {}

    def mulit_mediation(bact_pc):
        print(f' {bact_pc} is processing...')
        media_df = media_df_origin.copy()
        media_df.insert(2, f'{bact_pc}_bact', bact_pca_df[bact_pc])
        cols = media_df.columns.tolist()
        cols[1], cols[2] = cols[2], cols[1]
        media_df = media_df[cols]

        x_col = media_df.columns[1]
        y_col = media_df.columns[2]
        encoded_y = media_df[y_col]

        # Y = cX + ε1
        model1 = sm.Logit(encoded_y,
                          sm.add_constant(media_df[x_col])).fit(disp=False)
        coef_c = model1.params.iloc[1]

        result_dict_single_mediaton = {}
        result_coef = {}

        def single_mediation(index_media_df):

            mediator_col = media_df.columns[index_media_df]

            try:

                # M=aX+ε2
                model2 = sm.OLS(media_df[mediator_col],
                                sm.add_constant(
                                    media_df[x_col])).fit(disp=False)
                coef_a = model2.params.iloc[1]  # get coef
                # Y=c'X+bM+ε3
                model3 = sm.Logit(
                    encoded_y,
                    sm.add_constant(media_df[[x_col,
                                              mediator_col]])).fit(disp=False)
                coef_c_prime = model3.params.iloc[1]
                coef_b = model3.params.iloc[2]
                #  X→M regression coef of a standard error
                a_se = model2.bse.iloc[1]
                #  M→Y regression coef of a standard error
                b_se = model3.bse.iloc[2]
                # compute t
                t_stat_a = coef_a / a_se  #  X→M regression coef t
                t_stat_b = coef_b / b_se  #  M→Y regression coef t
                # z_ab
                z_ab = t_stat_a * t_stat_b
                # Bootstrap compute p
                bootstrap_z_values = bootstrap_p_values(
                    media_df, encoded_y, x_col, mediator_col, n_bootstrap)
                # p
                p_value = (np.sum(
                    np.abs(np.array(bootstrap_z_values)) >= np.abs(z_ab)) +
                           1) / (n_bootstrap + 1)

                # significant
                significance = "Significant" if p_value < 0.05 else "Insignificant"

                result_dict_single_mediaton[mediator_col] = (p_value,
                                                             significance)
                result_coef[mediator_col] = (coef_a, coef_b, coef_c,
                                             coef_c_prime, coef_a * coef_b)

            except Exception as e:
                print(f"Error in dealing wtih {mediator_col} ：{e}")
                pass

        # loop mediation
        for index_media_df in range(3, len(media_df.columns)):
            single_mediation(index_media_df)

        # get sig PC
        significant_pc = [
            pc for pc, (p_value,
                        significance) in result_dict_single_mediaton.items()
            if significance == 'Significant'
        ]

        # create DataFrame included sig PC
        significant_pc_df = pd.DataFrame(significant_pc,
                                         columns=['Significant_PC'])
        result_all[bact_pc] = significant_pc_df

        # get coef
        result_coef = pd.DataFrame.from_dict(result_coef,
                                             orient='index',
                                             columns=[
                                                 'coef_a', 'coef_b', 'coef_c',
                                                 'coef_c_prime',
                                                 'coef_a * coef_b'
                                             ])
        result_coef.reset_index(inplace=True)
        result_coef.rename(columns={'index': 'Metabolite_PC'}, inplace=True)
        result_all_bact_coef[bact_pc] = result_coef

    for bact_pc in bact_pca_df.columns:
        mulit_mediation(bact_pc)

    result_finall_df = pd.DataFrame()

    for bact_pc, df in result_all.items():
        df['Bacteria_PC'] = bact_pc
        result_finall_df = pd.concat([result_finall_df, df], ignore_index=True)

    result_final_coef = pd.concat(result_all_bact_coef).reset_index(
        level=0).rename(columns={'level_0': 'Bacteria_PC'})

    return result_finall_df, result_final_coef
