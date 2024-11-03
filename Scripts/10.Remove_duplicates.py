import pandas as pd
# Add this line at the beginning
input_path = 'Mesoderm/bnplot_new_bn200_meso_stochastic.csv' 

# Load the data
data = pd.read_csv(input_path)

data['pair'] = data.apply(lambda row: tuple(sorted([row['from'], row['to']])), axis=1)

# Remove duplicates based on the pair
unique_data = data.drop_duplicates(subset='pair')

# Drop the helper 'pair' column
unique_data = unique_data.drop(columns=['pair'])

# Save the unique data to a new CSV file
output_path = 'Mesoderm/removed_bnplot_new_bn200_meso_stochastic.csv'
unique_data.to_csv(output_path, index=False)

output_path
