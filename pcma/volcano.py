import subprocess
from unittest import result
from pkg_resources import resource_filename


def volcano(resource_data,
            plot_dir,
            y_value="p",
            y_threshold=0.05,
            SCC_threshold=0.2):
    volcano_script_path = resource_filename('pcma', 'rcode/volcano.R')
    y_threshold = str(y_threshold)
    SCC_threshold = str(SCC_threshold)
    result = subprocess.run([
        "Rscript", volcano_script_path, resource_data, plot_dir, y_value,
        y_threshold, SCC_threshold
    ],
                            capture_output=True,
                            text=True)
    print(result.stdout)
