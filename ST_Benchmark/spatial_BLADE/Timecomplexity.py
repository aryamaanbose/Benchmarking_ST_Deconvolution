
import pandas as pd
import numpy as np

from BLADE import Framework_Iterative

import os
import time

import matplotlib.pyplot as plt
script_path = os.path.abspath(__file__)

# Manually backtrack to the project root directory
project_root = os.path.join(script_path, *[os.path.pardir]*4)  # Adjust the number of os.path.pardir based on the depth
project_root = os.path.abspath(project_root)

print(f"Project root directory: {project_root}")


start_start_time = time.time()

##Load Data
# 1. Load the Spatial Data
# Load the spatial data
# Construct the path to your data file relative to the project root
data_file_path = os.path.join(project_root,"Processing", "BLADE", "Data","Processed_BLADE", "BLADE_spatial.csv")

# Read the CSV file into a pandas DataFrame
spatial_data = pd.read_csv(data_file_path, sep=",", header=None)

###Mouse olfactory bulb data taken from CARD
# Set the first row as column headers
spatial_data.columns = spatial_data.iloc[0]

# Rename the first column (currently NaN)
spatial_data.columns.values[0] = 'Gene'

# Set the first column as the row index
spatial_data = spatial_data.set_index('Gene')

# Drop the first row as it's now set as the header
spatial_data = spatial_data.drop(spatial_data.index[0])

spatial_data.replace(0, 0.000000001, inplace=True)

data_file_path = os.path.join(project_root, "Processing", "BLADE", "Data","Processed_BLADE", "BLADE_scrna_Mean.csv")

# Read the CSV file into a pandas DataFrame
scrna_mean = pd.read_csv(data_file_path, sep=",", header=None)

# Set the first row as the column headers
scrna_mean.columns = scrna_mean.iloc[0]

# Drop the first row as it is now set as the column header
scrna_mean = scrna_mean.drop(scrna_mean.index[0])
# Set the first column as the row indices
scrna_mean = scrna_mean.set_index(scrna_mean.columns[0])



data_file_path = os.path.join(project_root,"Processing", "BLADE", "Data","Processed_BLADE", "BLADE_scrna_SD.csv")

# Read the CSV file into a pandas DataFrame
scrna_sd = pd.read_csv(data_file_path, sep=",", header=None)

# Set the first row as the column headers
scrna_sd.columns = scrna_sd.iloc[0]

# Drop the first row as it is now set as the column header
scrna_sd = scrna_sd.drop(scrna_sd.index[0])
# Set the first column as the row indices
scrna_sd = scrna_sd.set_index(scrna_sd.columns[0])

scrna_mean.replace(0, 0.000000001, inplace=True)

# Replace all zeros with 0.00001 in the scrna_sd DataFrame
scrna_sd.replace(0, 0.000000001, inplace=True)

###Subset the data

import random

# Assuming spatial_data, scrna_mean, and scrna_sd are pre-defined DataFrames

# Find common genes
common_genes = set(spatial_data.index) & set(scrna_mean.index) & set(scrna_sd.index)
common_genes_list = list(common_genes)

# Subset Data to include only common genes
spatial_data_common = spatial_data.loc[common_genes_list]
scrna_mean_common = scrna_mean.loc[common_genes_list]
scrna_sd_common = scrna_sd.loc[common_genes_list]

# Ensure that the data is in a numec format
spatial_data_common = spatial_data_common.apply(pd.to_numeric, errors='coerce')
scrna_mean_common = scrna_mean_common.apply(pd.to_numeric, errors='coerce')
scrna_sd_common = scrna_sd_common.apply(pd.to_numeric, errors='coerce')

###Run Blade


 # Set hyperparameters for BLADE
hyperpars = {
    'Alpha': [1],
    'Kappa0': [1],
    'SY': [1]}

Nrep = 10
Nrepfinal = 100
Njob = 10

# Prepare the input data for BLADE
mean = scrna_mean_common.values
sd = scrna_sd_common.values
Y = spatial_data_common.values



# Apply the BLADE framework


Y_subset = Y[:, :278]



# Define the sequence of column counts to test
num_columns_list = np.arange(278,8,-27)

# Placeholder for timing results
timing_results = []

for num_columns in num_columns_list:
    # Subset `Y_subset` to the first `num_columns` columns
    Y_current = Y_subset[:, :num_columns]
    print(Y_subset.shape,"Y_subset")
    print(Y_current.shape, "Y_current")
    # Time the operation you're interested in
    start_time = time.time()
    final_obj, best_obj, best_set, outs = Framework_Iterative(
        mean,
        sd,
        Y_current,
        Alpha=hyperpars['Alpha'],
        Nrep=Nrep,
        Njob=Njob,
        Nrepfinal=Nrepfinal
    )
    end_time = time.time()

    # Calculate elapsed time and store it
    elapsed_time = end_time - start_time
    timing_results.append((num_columns, elapsed_time))

# Convert timing results to a DataFrame
df_results = pd.DataFrame(timing_results, columns=['Num_Columns', 'Execution_Time'])

results_dir = os.path.join(project_root,"Processing", "BLADE", "Spatial_Blade","Results","execution_time_vs_samples.png")

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(df_results['Num_Columns'], df_results['Execution_Time'], marker='o')
plt.title('Execution Time vs Number of Columns in Spatial Data')
plt.xlabel('Number of Columns')
plt.ylabel('Execution Time (s)')
plt.grid(True)
plt.savefig(results_dir)

print("Saved file")



end_end_time = time.time()


elapsed_time2 = end_end_time - start_start_time

print(elapsed_time2)
