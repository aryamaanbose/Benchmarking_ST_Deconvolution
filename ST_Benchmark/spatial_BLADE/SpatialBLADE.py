##Import

import pandas as pd
import numpy as np



from BLADE import Framework_Iterative

import os

import matplotlib.pyplot as plt


script_path = os.path.abspath(__file__)

# Manually backtrack to the project root directory
project_root = os.path.join(script_path, *[os.path.pardir]*4)  # Adjust the number of os.path.pardir based on the depth
project_root = os.path.abspath(project_root)

print(f"Project root directory: {project_root}")




##Load Data
# 1. Load the Spatial Data
# Load the spatial data
# Construct the path to your data file relative to the project root
data_file_path = os.path.join(project_root,"Processing", "BLADE", "Data","Processed_BLADE", "SpatialDataMOB.csv")

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

#spatial_data.replace(0, 0.00001, inplace=True)

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

# Remove genes starting with "Rp" from scrna_mean
#scrna_mean = scrna_mean[~scrna_mean.index.str.startswith("Rp")]

# Remove genes starting with "Rp" from scrna_sd
#scrna_sd = scrna_sd[~scrna_sd.index.str.startswith("Rp")]

print(scrna_mean.shape)
print(scrna_sd.shape)

#scrna_mean.replace(0, 0.00001, inplace=True)

# Replace all zeros with 0.00001 in the scrna_sd DataFrame
#scrna_sd.replace(0, 0.00001, inplace=True)

###Subset the data

import random

# Assuming spatial_data, scrna_mean, and scrna_sd are pre-defined DataFrames

# Find common genes
common_genes = set(spatial_data.index) & set(scrna_mean.index) & set(scrna_sd.index)
common_genes_list = list(common_genes)

##Autogene_list
#gene_list2 = ['Gm26901', 'Lman2l', 'Ormdl1', 'Fam117b', 'Fzd5', 'Igfbp5', 'Scg2', 'Dock10', 'Krtap28-13', 'Itm2c', 'Tmem183a', 'B3galt2', 'Vamp4', 'Ndufs2', 'Pex19', 'Opn3', 'Olfm1', 'Zeb2', 'Ttc21b', 'Dlx1', 'Olfr1056', 'C1qtnf4', 'Ano3', 'Rmdn3', 'Exd1', 'Nusap1', 'Casc4', 'Slc20a1', 'Chgb', 'Napb', 'Rbm39', 'Snhg11', 'Rims4', 'Zfos1', 'Syp', 'Pgrmc1', 'Xpnpep2', 'Phf6', 'Dkc1', 'Xist', 'Drp2', 'Bex2', 'Ngfrap1', 'Phf8', 'Ppef1', 'Piga', 'Sec62', 'Ttc14', '4932438A13Rik', 'Ccdc169', 'Fstl5', 'Gria2', 'Tsacc', 'S100a5', 'Gm43062', 'Dennd2c', 'Wdr77', 'Slc16a4', 'Ntng1', 'Plppr4', 'Sep15', 'Gng5', 'Pnisr', 'Gm20878', 'Col15a1', 'Whrn', 'Zmym4', 'Zbtb40', 'Cda', 'Ttc34', 'Mterf1a', 'Reln', 'Lhfpl3', 'Fosl2', 'Epha5', 'Adamts3', 'Pf4', 'Tmem150cos', 'Rabgef1', 'Gnb2', 'Prkar1b', 'Gm43597', 'Met', 'Nap1l5', 'Stambp', 'Pcyox1', 'March8', 'Etnk1', 'Mboat7', 'Chmp2a', 'Zfp658', 'Slc17a7', 'Ube3a', 'Nmb', 'Omp', 'Plekhb1', 'Spon1', 'RP23-44H21.1', 'Grm1', 'Themis', 'Fam229b', 'Serinc1', 'Atp5d', 'Gm15608', 'Syt1', 'Rdh16', 'Wrn', 'Cpe', 'Hmgxb4', 'Nrn1l', 'Calb2', 'Rfwd3', 'Zcchc14', 'Gm26759', 'Gm2974', 'Fhit', 'Synpr', 'Il3ra', 'Cdhr1', 'Ipo4', 'Ppp2r2a', 'Kbtbd3', 'Trpc6', 'Zfp426', 'Slc37a2', '2010007H06Rik', 'Nptn', 'B930082K07Rik', 'Map2k1', 'Rwdd2a', 'Mthfsl', 'Tcta', 'Cck', 'Nktr', 'Rtn4', 'Gabra1', 'Gabrb2', 'Adra1b', 'Rpl26', 'Aurkb', 'Gm11205', 'Crlf3', 'Bzrap1', 'Pdk2', 'Ikzf3', 'Krt222', 'Map3k3', 'Hn1', 'H3f3b', 'Srsf2', 'Tnrc6c', 'Tk1', 'Tha1', 'Elmo1', 'Hist1h3e', 'Nrsn1', 'Nedd9', 'Hk3', 'Zfp808', 'Ankrd55', 'Mrps30', 'Paip1', '4921508M14Rik', 'Tspan13', 'Etv1', 'Rtn1', 'Six1', 'Meg3', 'Amn', 'Rictor', 'Rpl30', 'Oxr1', 'Adcy8', '1700109K24Rik', 'Syngr1', 'Fbln1', 'Selo', '5330439K02Rik', 'Tns2', 'Alg1', 'Gm15738', 'Chrd', 'Gsk3b', 'Tagln3', 'St3gal6', 'Son', 'Sft2d1', 'Slc22a3', 'Ndufb10', 'Gng13', 'Gm28043', 'H2-D1', 'Safb2', 'Wdr43', 'Taf4b', 'Pcdhac2', 'Ablim1', '6430562O15Rik', 'Gm9821', 'Lrrc71']




##CARD_list

# Subset Data to include only common genes
spatial_data_common = spatial_data.loc[common_genes]
scrna_mean_common = scrna_mean.loc[common_genes]
scrna_sd_common = scrna_sd.loc[common_genes]

# Ensure that the data is in a numec format
spatial_data_common = spatial_data_common.apply(pd.to_numeric, errors='coerce')
scrna_mean_common = scrna_mean_common.apply(pd.to_numeric, errors='coerce')
scrna_sd_common = scrna_sd_common.apply(pd.to_numeric, errors='coerce')

# Debug: Print shapes after subsetting
print(f"Spatial data shape after subsetting: {spatial_data_common.shape}")
print(f"scrna_mean shape after subsetting: {scrna_mean_common.shape}")
print(f"scrna_sd shape after subsetting: {scrna_sd_common.shape}")


###Run Blade


 # Set hyperparameters for BLADE
hyperpars = {
    'Alpha': [1],
    'Kappa0': [1],
    'SY': [1]}

Nrep = 10
Nrepfinal = 1000
Njob = 10

# Prepare the input data for BLADE
mean = scrna_mean_common.values
sd = scrna_sd_common.values
Y = spatial_data_common.values

Y = Y[:, :10]

print(f"Dimensions of mean: {mean.shape}")
print(f"Dimensions of sd: {sd.shape}")
print(f"Dimensions of Y: {Y.shape}")

# Apply the BLADE framework



final_obj, best_obj, best_set, outs = Framework_Iterative(
    mean,
    sd,
    Y,
    Alpha=hyperpars['Alpha'],
    Nrep=Nrep,
    Njob=Njob,
    Nrepfinal=Nrepfinal
)





# Assuming you've loaded the scRNA reference dataset as scrna_mean_common
cell_type_labels = scrna_mean_common.columns.tolist()

# Extracting cell type proportions
# Ensure that the BLADE framework has been applied to your spatial data and that final_obj_spatial is defined
cell_fractions_spatial = final_obj.ExpF(final_obj.Beta)
##ExpF is used to convert them into fractions

#df_alpha = pd.DataFrame(Alpha, columns=cell_type_labels)






# Convert to DataFrame
df_cell_fractions_spatial = pd.DataFrame(cell_fractions_spatial, columns=cell_type_labels)

results_dir = os.path.join(project_root,"Processing", "BLADE", "Spatial_Blade", "Results") # Adjusted for a more typical relative path
base_file_name = "BLADE_results"
extension = ".csv"

# Ensure the results directory exists
if not os.path.exists(results_dir):
    os.makedirs(results_dir)

# Initialize the run number and output file name
run_number = 1
output_file = os.path.join(results_dir, f"{base_file_name}{extension}")

# Check if the file already exists
while os.path.exists(output_file):
    # If it exists, update the file name with a new run number
    output_file = os.path.join(results_dir, f"{base_file_name}_{run_number}{extension}")
    run_number += 1

# Save the DataFrame to CSV in the results directory with the new file name
df_cell_fractions_spatial.to_csv(output_file, index=False)




df_alpha = pd.DataFrame(final_obj.Alpha, columns=cell_type_labels)
df_alpha.to_csv(os.path.join(results_dir,"alpha_selectgenes.csv"))


df_beta = pd.DataFrame(final_obj.Beta, columns=cell_type_labels)
df_beta.to_csv(os.path.join(results_dir,"beta_selectgenes.csv"))




print(final_obj.Alpha)
print(final_obj.Beta)
print(final_obj.ExpF(final_obj.Beta))



print(f"DataFrame saved to {output_file}")

plt.scatter(x=range(len(best_obj)), y=best_obj)  # Adapt based on actual data structure
plt.xlabel('Iteration')
plt.ylabel('Objective Function Value')
plt.title('BLADE Convergence Check')
plt.show()



