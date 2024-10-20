import pandas as pd

# Load the uploaded files
filtered_mesodermal_genes_df = pd.read_csv('Mesoderm/unique_bnplot_new_bn200_meso_stochastic.csv')
absolute_pathway_df = pd.read_csv('Mesoderm/absolute_meso_stochastic.csv')

# Merge the two dataframes on the Pathways column
merged_df = filtered_mesodermal_genes_df.merge(
    absolute_pathway_df, 
    left_on='from', 
    right_on='Pathways', 
    how='left'
)

# Add the 'Values' column to the 'filtered_mesodermal_genes_df'
merged_df['strength_with_values'] = merged_df['strength'] + merged_df['Values'].fillna(0)

# Save the result to a new CSV file
output_path = 'Mesoderm/abscorp_bnplot_new_bn200_meso_stochastic.csv'
merged_df.to_csv(output_path, index=False)


