#%%
import pandas as pd
from scipy.stats import spearmanr
from itertools import product
from statsmodels.stats.multitest import multipletests


def compute_correlation(df1: pd.DataFrame,
                        df2: pd.DataFrame,
                        threshold=0.2,
                        p_threshold=0.05,
                        filter_method='p'):
    '''compute the correlation between 2 dataframe and filte with threshold'''
    pairs = list(product(df1.columns, df2.columns))

    def calc_correlation(pair):
        column_A, column_B = pair
        correlation, p_value = spearmanr(df1[column_A], df2[column_B])
        return {
            "Column_A": column_A,
            "Column_B": column_B,
            "Correlation": correlation,
            "P_value": p_value
        }

    result_data = pd.Series(pairs).apply(calc_correlation).tolist()
    result_data = pd.DataFrame(result_data)
    # FDR
    p_values = result_data['P_value']
    rejected, fdr_corrected_p_values, _, _ = multipletests(p_values,
                                                           method='fdr_bh')
    result_data['FDR'] = fdr_corrected_p_values
    result_data_filter = result_data[result_data['Correlation'].abs() >=
                                     threshold]
    if filter_method == 'p':
        result_data_filter = result_data_filter[result_data_filter['P_value']
                                                <= p_threshold]
    if filter_method == 'fdr':
        result_data_filter = result_data_filter[result_data_filter['FDR'] <=
                                                p_threshold]
    return result_data, result_data_filter


# %%
