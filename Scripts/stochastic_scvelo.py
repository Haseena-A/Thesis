import sys
import os
import scvelo as scv
import pandas as pd

# Set up scVelo figure parameters
scv.settings.verbosity = 3  # Show detailed output information
scv.settings.presenter_view = True  # Enable presenter view for wider plots
scv.set_figure_params('scvelo')  # For beautified visualization

# Load the Loom file into an AnnData object
loom_file = '/home/haseena/endoderm/velocity/all_endoderm.loom'  # Replace with the actual path
adata = scv.read(loom_file, cache=True)

# Preprocess the data
scv.pp.filter_and_normalize(adata, min_shared_counts=2, n_top_genes=adata.shape[1])
scv.pp.moments(adata)
scv.tl.velocity(adata, mode='stochastic')

# Convert RNA velocity data to a DataFrame with proper row and column names
data1 = pd.DataFrame(adata.layers['velocity'], index=adata.obs_names, columns=adata.var_names)


# Add a new row with a name and values
# Example: Adding a row named 'new_row' with values [1, 2, ..., n] where n is the number of columns
new_row_name = 'new_row'
new_row_values = [1] * data1.shape[1]
data1.loc[new_row_name] = new_row_values

# Save the modified DataFrame to CSV
data1.to_csv('/home/haseena/endoderm/velocity/all_endoderm_stochastic.csv', index=True)

