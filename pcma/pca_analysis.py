import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, MinMaxScaler


def run_pca(get_bact_meta: pd.DataFrame,
            PC_num=50,
            norm_method='MinMaxScaler',
            if_whiten=True,
            Norm=True,
            seeds=0):

    # number of PC
    if PC_num != 'mle':
        if len(get_bact_meta.columns) >= 50:
            n_components = int(PC_num)
        else:
            n_components = get_bact_meta.shape[1]
    else:
        n_components = 'mle'

    # if normalization
    normalize_data = Norm

    # if whiten
    whiten = if_whiten

    # svd_solver
    svd_solver = 'auto'  # automatic

    # random seed
    random_state = seeds

    # if Normalization
    if Norm == True:
        if norm_method == 'MinMaxScaler':
            scaled_data = MinMaxScaler().fit_transform(get_bact_meta)
        if norm_method == 'StandardScaler':
            scaled_data = StandardScaler().fit_transform(get_bact_meta)
        if norm_method == 'MeanNormalization':
            mean_vals = get_bact_meta.mean()
            range_vals = get_bact_meta.max() - get_bact_meta.min()
            scaled_data = get_bact_meta.copy()
            # escape from all 0
            for column in get_bact_meta.columns:
                if range_vals[column] != 0:
                    scaled_data[column] = (
                        get_bact_meta[column] -
                        mean_vals[column]) / range_vals[column]
                else:
                    scaled_data[column] = 0
        if norm_method == 'Centralization':
            mean_vals = get_bact_meta.mean()
            scaled_data = get_bact_meta - mean_vals

    else:
        # do not normalization
        scaled_data = get_bact_meta

    # PCA
    pca = PCA(n_components=n_components,
              whiten=whiten,
              svd_solver=svd_solver,
              random_state=random_state)

    pca_result = pca.fit_transform(scaled_data)

    explained_variance_ratio = pca.explained_variance_ratio_
    explained_variance_ratio = pd.DataFrame(explained_variance_ratio)
    components = pca.components_
    components_df = pd.DataFrame(data=components,
                                 columns=get_bact_meta.columns)
    pca_df = pd.DataFrame(
        data=pca_result,
        columns=[f'PC{i}' for i in range(1, pca_result.shape[1] + 1)])
    # dealwhith scale data
    scaled_data = pd.DataFrame(scaled_data)
    scaled_data.columns = get_bact_meta.columns
    return scaled_data, pca_df, explained_variance_ratio, components_df
