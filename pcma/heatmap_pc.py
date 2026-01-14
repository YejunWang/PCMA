import os
import pandas as pd
import subprocess
from pkg_resources import resource_filename


# %%
def heatmap_pc(PC_df, result_dir, plot_dir, method):
    heatmap_script_path = resource_filename('pcma', 'rcode/heatmap_pc.R')
    # data cleaning
    if method == 'pcpcma':
        all_pcs = PC_df['Bacteria_PC']
        pc_numbers = [int(pc[2:]) for pc in all_pcs if pc.startswith('PC')]
        max_pc = max(pc_numbers)
        all_pcs = ['PC' + str(i) for i in range(1, max_pc + 1)]
        existing_pcs = PC_df['Bacteria_PC'].unique()
        missing_pcs = [pc for pc in all_pcs if pc not in existing_pcs]

        missing_rows = pd.DataFrame({
            'Significant_PC': [None] * len(missing_pcs),
            'Bacteria_PC': missing_pcs
        })

        PC_df = pd.concat([PC_df, missing_rows], ignore_index=True)
        # rank
        PC_df['Sort_Key'] = PC_df['Bacteria_PC'].str.extract('(\d+)').astype(
            int)
        PC_df.sort_values('Sort_Key', inplace=True)
        PC_df.drop('Sort_Key', axis=1, inplace=True)
    if method == 'pcma1':
        pass
    if method == 'pcma2':
        pass

    # plot

    data_dir = os.path.join(result_dir, "heatmap_data.csv")
    PC_df.to_csv(data_dir, index=False)
    r_result = subprocess.run(
        ["Rscript", heatmap_script_path, data_dir, plot_dir],
        capture_output=True,
        text=True)
    print(r_result.stdout)
    return (PC_df)
