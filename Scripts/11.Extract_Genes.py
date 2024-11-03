import pandas as pd

# Load the CSV file
file_path = 'Ectoderm/removed_bnplot_new_bn200_ecto_stochastic.csv'
data = pd.read_csv(file_path)

# Filter rows where "from" or "to" columns contain "MESODERMAL_GENES"
filtered_data = data[(data['from'].str.contains('ECTODERM')) | (data['to'].str.contains('ECTODERM'))]

# Save the filtered data to a new CSV file
filtered_file_path = 'Ectoderm/unique_bnplot_new_bn200_ecto_stochastic.csv'
filtered_data.to_csv(filtered_file_path, index=False)


