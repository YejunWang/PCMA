import subprocess
from pkg_resources import resource_filename


def heatmap_pca(PCA_data, plot_dir):
    heatmap_scale_script_path = resource_filename('pcma',
                                                  'rcode/heatmap_pca.R')
    result = subprocess.run(
        ["Rscript", heatmap_scale_script_path, PCA_data, plot_dir],
        capture_output=True,
        text=True)
    print(result.stdout)
