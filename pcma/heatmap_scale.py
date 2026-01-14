import subprocess
from pkg_resources import resource_filename


def heatmap_scale(Metabolite_dir, scale_data, plot_dir_1, plot_dir_2):
    heatmap_scale_script_path = resource_filename('pcma',
                                                  'rcode/heatmap_scale.R')
    result = subprocess.run([
        "Rscript", heatmap_scale_script_path, Metabolite_dir, scale_data,
        plot_dir_1, plot_dir_2
    ],
                            capture_output=True,
                            text=True)
    print(result.stdout)
