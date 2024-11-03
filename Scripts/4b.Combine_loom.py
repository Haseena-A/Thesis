
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
import loompy
import os

# Directory containing the loom files
loom_directory = "/home/vibhor/haseena_ip/mps/loom"

# List all loom files in the directory
loom_files = [os.path.join(loom_directory, f) for f in os.listdir(loom_directory) if f.endswith(".loom")]

# Combine the loom files
loompy.combine(loom_files, output_file="/home/vibhor/haseena_ip/mps/combined_loom/mps.loom", key="Accession")

print("Loom files combined successfully into 'combined.loom'")


