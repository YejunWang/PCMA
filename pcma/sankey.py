#%%
import pandas as pd
import plotly.graph_objects as go
import subprocess
from pkg_resources import resource_filename


#%%
# Load the data
def sankey(Bacteria_dir, Metabolite_pca_dir, Diag_dir, Sig_dir,
           output_data_dir, output_png_dir):
    sankey_script_path = resource_filename('pcma', 'rcode/sankey.R')
    result = subprocess.run([
        "Rscript", sankey_script_path, Bacteria_dir, Metabolite_pca_dir,
        Diag_dir, Sig_dir, output_data_dir
    ],
                            capture_output=True,
                            text=True)
    print(result.stdout)
    data = pd.read_csv(output_data_dir)

    # Create node labels
    all_labels = list(pd.concat([data['source'], data['target']]).unique())
    label_dict = {label: idx for idx, label in enumerate(all_labels)}

    # Create node and link data
    source_indices = data['source'].map(label_dict)
    target_indices = data['target'].map(label_dict)
    values = data['Spearman_Correlation'].abs()

    # Create node labels
    all_labels = list(pd.concat([data['source'], data['target']]).unique())
    label_dict = {label: idx for idx, label in enumerate(all_labels)}

    # Create node and link data
    source_indices = data['source'].map(label_dict)
    target_indices = data['target'].map(label_dict)
    values = data['Spearman_Correlation'].abs()

    # Set colors based on Spearman Correlation
    colors = [
        '#F48892' if val > 0 else '#91CAE8'
        for val in data['Spearman_Correlation']
    ]

    # Create the Sankey diagram
    fig = go.Figure(data=[
        go.Sankey(node=dict(pad=15,
                            thickness=20,
                            line=dict(color="black", width=0.5),
                            label=all_labels),
                  link=dict(source=source_indices,
                            target=target_indices,
                            value=values,
                            color=colors))
    ])
    fig.update_layout(title_text="Sankey Plot", font_size=10)
    fig.write_image(output_png_dir, scale=8)
    #fig.show()
