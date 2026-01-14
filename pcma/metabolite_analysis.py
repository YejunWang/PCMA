import os
import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def metabolite_analysis(components_file_path,
                        label_file_path,
                        significant_pc_path,
                        fig_output,
                        file_output,
                        pathway_num=3,
                        alpha=0.05,
                        label_column_index=1,
                        coff_left=0.5):
    '''Labels processing'''

    def process_pca_components(data_df,
                               label_df,
                               label_column_index=1,
                               coff_left=0.5):
        feature_names = data_df.columns.to_numpy()
        data_matrix = data_df.to_numpy()
        labels = label_df.set_index(label_df.columns[0])[
            label_df.columns[label_column_index]].to_dict()

        def retain_positive_negative_50_percent(row):
            positive_indices = np.where(row > 0)[0]
            negative_indices = np.where(row < 0)[0]
            sorted_positive_indices = positive_indices[np.argsort(
                row[positive_indices])[::-1]]
            sorted_negative_indices = negative_indices[np.argsort(
                np.abs(row[negative_indices]))[::-1]]
            total_positive_sum = np.sum(row[sorted_positive_indices])
            total_negative_sum = np.sum(np.abs(row[sorted_negative_indices]))
            cumulative_positive_sum = 0
            cumulative_negative_sum = 0
            retained_positive_indices = []
            retained_negative_indices = []
            for idx in sorted_positive_indices:
                cumulative_positive_sum += row[idx]
                retained_positive_indices.append(idx)
                if cumulative_positive_sum / total_positive_sum >= coff_left:
                    break
            for idx in sorted_negative_indices:
                cumulative_negative_sum += np.abs(row[idx])
                retained_negative_indices.append(idx)
                if cumulative_negative_sum / total_negative_sum >= coff_left:
                    break
            retained_indices = retained_positive_indices + retained_negative_indices
            return retained_indices

        feature_coefficient_matrix = []
        for i in range(data_matrix.shape[0]):
            retained_indices = retain_positive_negative_50_percent(
                data_matrix[i])
            retained_features = feature_names[retained_indices]
            retained_labels = [
                labels[feature] for feature in retained_features
            ]
            retained_coefficients = data_matrix[i, retained_indices]
            feature_coefficient_matrix.append(['PC' + str(i + 1)] +
                                              retained_features.tolist())
            feature_coefficient_matrix.append(['PC' + str(i + 1)] +
                                              retained_labels)
            feature_coefficient_matrix.append(['PC' + str(i + 1)] +
                                              retained_coefficients.tolist())

        flattened_data = []
        for trio in feature_coefficient_matrix:
            flattened_data.append(trio)

        result_df = pd.DataFrame(flattened_data)
        return result_df

    '''Labels counting'''

    def count_labels_per_pc(result_df):
        pc_label_counts = {}
        all_labels = set()
        for i in range(1, result_df.shape[0], 3):
            pc_name = result_df.iloc[i, 0]
            labels = result_df.iloc[i, 1:].dropna()
            coefficients = result_df.iloc[i + 1, 1:].dropna().astype(float)
            pos_labels = labels[coefficients > 0]
            neg_labels = labels[coefficients < 0]
            if pc_name not in pc_label_counts:
                pc_label_counts[pc_name] = {'positive': {}, 'negative': {}}
            for label in pos_labels:
                all_labels.add(label)
                if label in pc_label_counts[pc_name]['positive']:
                    pc_label_counts[pc_name]['positive'][label] += 1
                else:
                    pc_label_counts[pc_name]['positive'][label] = 1
            for label in neg_labels:
                all_labels.add(label)
                if label in pc_label_counts[pc_name]['negative']:
                    pc_label_counts[pc_name]['negative'][label] += 1
                else:
                    pc_label_counts[pc_name]['negative'][label] = 1

        pc_names = sorted(pc_label_counts.keys())
        all_labels = sorted(all_labels)
        result_dict = {'PC': pc_names}

        for label in all_labels:
            result_dict[f'positive_{label}'] = []
            result_dict[f'negative_{label}'] = []
            for pc in pc_names:
                pos_count = pc_label_counts[pc]['positive'].get(label, 0)
                neg_count = pc_label_counts[pc]['negative'].get(label, 0)
                result_dict[f'positive_{label}'].append(pos_count)
                result_dict[f'negative_{label}'].append(neg_count)

        result_df = pd.DataFrame(result_dict)
        return result_df

    '''Randomly select rows'''

    def select_random_pcs(file_path, pathway_num=3):
        df = pd.read_csv(file_path)
        pathway_num = min(pathway_num, len(df))
        selected_rows = df.sample(n=pathway_num)
        return selected_rows

    '''Plotting '''

    def plot_pc_data(pc, data1, data2):
        data1_pc = data1[data1['PC'] == pc].drop(columns=['PC'])
        data2_pc = data2[data2['PC'] == pc].drop(columns=['PC'])
        positive_labels = [
            col for col in data1_pc.columns if 'positive' in col
        ]
        negative_labels = [
            col for col in data1_pc.columns if 'negative' in col
        ]
        data1_pos = data1_pc[positive_labels]
        data1_neg = data1_pc[negative_labels]
        data2_pos = data2_pc[positive_labels]
        data2_neg = data2_pc[negative_labels]
        fig1, ax1 = plt.subplots(figsize=(30, 60))
        bar_positions = np.arange(len(positive_labels))
        ax1.barh(bar_positions,
                 data1_pos.values.flatten(),
                 color='#1f77b4',
                 label='N50')
        ax1.barh(bar_positions,
                 data2_pos.values.flatten(),
                 left=data1_pos.values.flatten(),
                 color='#ff7f0e',
                 label='All')
        ax1.set_yticks(bar_positions)
        ax1.set_yticklabels(positive_labels, fontsize=5)
        ax1.set_title(f'barplot_fisher_{pc}_positive')
        ax1.legend()
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.xticks(fontsize=10)
        plt.yticks(rotation=45)
        plt.tight_layout()
        fig2, ax2 = plt.subplots(figsize=(30, 60))
        bar_positions = np.arange(len(negative_labels))
        ax2.barh(bar_positions,
                 data1_neg.values.flatten(),
                 color='#1f77b4',
                 label='enriched')
        ax2.barh(bar_positions,
                 data2_neg.values.flatten(),
                 left=data1_neg.values.flatten(),
                 color='#ff7f0e',
                 label='original')
        ax2.set_yticks(bar_positions)
        ax2.set_yticklabels(negative_labels, fontsize=5)
        ax2.set_title(f'barplot_fisher_{pc}_negative')
        ax2.legend()
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.xticks(fontsize=10)
        plt.yticks(rotation=45)
        plt.tight_layout()
        return fig1, fig2

    '''Fisher'''

    def fisher_test_enrichment(df1, df2, alpha=0.05):
        all_labels = df2.columns[1:]
        pcs = df2['PC']

        # fill labels that are in df2 but not in df1
        for label in all_labels:
            if label not in df1.columns:
                df1[label] = 0

        df1_matrix = df1.set_index('PC').reindex(pcs).fillna(
            0)[all_labels].values
        df2_matrix = df2.set_index('PC')[all_labels].values
        df1_sums = np.tile(df1[all_labels].sum().values,
                           (df1_matrix.shape[0], 1))
        df2_sums = np.tile(df2[all_labels].sum().values,
                           (df2_matrix.shape[0], 1))
        table_matrix_11 = df1_matrix
        table_matrix_12 = df1_sums - df1_matrix
        table_matrix_21 = df2_matrix
        table_matrix_22 = df2_sums - df2_matrix
        results_list = []
        for i in range(df1_matrix.shape[0]):
            for j in range(df1_matrix.shape[1]):
                table = np.array(
                    [[table_matrix_11[i, j], table_matrix_12[i, j]],
                     [table_matrix_21[i, j], table_matrix_22[i, j]]])
                _, p_value = stats.fisher_exact(table)
                results_list.append({
                    'PC': pcs[i],
                    'Label': all_labels[j],
                    'p_value': p_value
                })
        results = pd.DataFrame(results_list)
        significant_results = results[results['p_value'] <= alpha]
        return significant_results, results

    # Main part
    df1 = pd.read_csv(components_file_path)
    df2 = pd.read_csv(label_file_path)
    df_result = process_pca_components(df1, df2, label_column_index, coff_left)
    df_comparison_result = process_pca_components(df1, df2, label_column_index,
                                                  1)
    label_count = count_labels_per_pc(df_result)
    label_count_comparison = count_labels_per_pc(df_comparison_result)
    selected_rows = select_random_pcs(significant_pc_path, pathway_num)

    with PdfPages(fig_output) as pdf:
        for idx, row in selected_rows.iterrows():
            pc = row['Significant_PC']
            bacteria_pc = row['Bacteria_PC']
            fig1, fig2 = plot_pc_data(pc, label_count, label_count_comparison)
            pdf.savefig(fig1)
            pdf.savefig(fig2)
            plt.close(fig1)
            plt.close(fig2)
            print(
                f'Significant metabolic pathways of bacteria {bacteria_pc} for PC {pc}'
            )

    significant_results, fisher_results = fisher_test_enrichment(
        label_count, label_count_comparison)
    final_significant_results = significant_results[
        significant_results['PC'].isin(selected_rows['Significant_PC'])]
    final_significant_results.to_csv(
        os.path.join(file_output, 'significant_pathway.csv'))
    fisher_results.to_csv(os.path.join(file_output, 'Fisher_results.csv'))
