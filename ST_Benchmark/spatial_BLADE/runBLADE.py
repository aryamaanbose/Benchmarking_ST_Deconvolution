import pandas as pd
import numpy as np



from BLADE import Framework_Iterative

import os





def run_blade(scrna_mean_file, scrna_sd_file,spatial, result_name, small_number=0.00001, expected_file=""):
    # Get the project root directory
    script_path = os.path.abspath(__file__)
    project_root = os.path.join(script_path, *[os.path.pardir]*4)
    project_root = os.path.abspath(project_root)

    print(f"Project root directory: {project_root}")

    def load_data(file_path):
        data = pd.read_csv(file_path, sep=",", header=None)
        data.columns = data.iloc[0]
        data = data.drop(data.index[0])
        data = data.set_index(data.columns[0])
        return data

    # Construct paths to the data files
    spatial_data_path = os.path.join(project_root, "Benchmarking_ST_Deconvolution", "ST_Benchmark", "Data", "Processed_BLADE", spatial)
    scrna_mean_path = os.path.join(project_root,  "Benchmarking_ST_Deconvolution", "ST_Benchmark", "Data", "Processed_BLADE", scrna_mean_file)
    scrna_sd_path = os.path.join(project_root, "Benchmarking_ST_Deconvolution", "ST_Benchmark", "Data", "Processed_BLADE", scrna_sd_file)

    # Load data
    spatial_data = load_data(spatial_data_path)
    scrna_mean = load_data(scrna_mean_path)
    scrna_sd = load_data(scrna_sd_path)

    # Load expected data if provided
    if expected_file:
        expected_path = os.path.join(project_root, "Processing", "BLADE", "Data", "Processed_BLADE", expected_file)
        expected_data = load_data(expected_path)
        expected = np.array(expected_data).astype(float)
    else:
        expected = None

    # Find common genes
    common_genes = list(set(spatial_data.index) & set(scrna_mean.index) & set(scrna_sd.index))

    # Subset data to include only common genes
    spatial_data_common = spatial_data.loc[common_genes]
    scrna_mean_common = scrna_mean.loc[common_genes]
    scrna_sd_common = scrna_sd.loc[common_genes]

    # Ensure that the data is in numeric format
    spatial_data_common = spatial_data_common.apply(pd.to_numeric, errors='coerce')
    scrna_mean_common = scrna_mean_common.apply(pd.to_numeric, errors='coerce')
    scrna_sd_common = scrna_sd_common.apply(pd.to_numeric, errors='coerce')

    # Replace zeros with a small number
    spatial_data_common.replace(0, small_number, inplace=True)
    scrna_mean_common.replace(0, small_number, inplace=True)
    scrna_sd_common.replace(0, small_number, inplace=True)

    # If expected data is provided, ensure spatial_data_common columns match the number of rows in expected
    if expected is not None:
        if spatial_data_common.shape[1] < expected.shape[0]:
            raise ValueError("Spatial data has fewer columns than expected data has rows. Please check your data.")
        spatial_data_common = spatial_data_common.iloc[:, :expected.shape[0]]

    # Debug: Print shapes after subsetting
    print(f"Spatial data shape after subsetting: {spatial_data_common.shape}")
    print(f"scrna_mean shape after subsetting: {scrna_mean_common.shape}")
    print(f"scrna_sd shape after subsetting: {scrna_sd_common.shape}")

    # Set hyperparameters for BLADE
    hyperpars = {
        'Alpha': [1],
        'Kappa0': [1],
        'SY': [1]
    }

    Nrep = 10
    Nrepfinal = 100
    Njob = 10

    # Prepare the input data for BLADE
    mean = scrna_mean_common.values
    sd = scrna_sd_common.values
    Y = spatial_data_common.values

    print(f"Dimensions of mean: {mean.shape}")
    print(f"Dimensions of sd: {sd.shape}")
    print(f"Dimensions of Y: {Y.shape}")

    # Apply the BLADE framework
    if expected is not None:
        final_obj, best_obj, best_set, outs = Framework_Iterative(
            mean,
            sd,
            Y,
            Alpha=hyperpars['Alpha'],
            Nrep=Nrep,
            Njob=Njob,
            Nrepfinal=Nrepfinal,
            Expectation=expected
        )
    else:
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
    cell_fractions_spatial = final_obj.ExpF(final_obj.Beta)

    # Convert to DataFrame
    df_cell_fractions_spatial = pd.DataFrame(cell_fractions_spatial, columns=cell_type_labels)

    results_dir = os.path.join(project_root,  "Benchmarking_ST_Deconvolution", "ST_Benchmark" ,"spatial_Blade", "Results_new")
    base_file_name = result_name
    extension = ".csv"

    # Ensure the results directory exists
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    # Initialize the output file name
    output_file = os.path.join(results_dir, f"{base_file_name}{extension}")

    # Save the DataFrame to CSV in the results directory
    df_cell_fractions_spatial.to_csv(output_file, index=False)

    print(f"Results saved to {output_file}")


script_path = os.path.abspath(__file__)
project_root = os.path.join(script_path, *[os.path.pardir]*4)
project_root = os.path.abspath(project_root)

print(f"Project root directory: {project_root}")

def load_data(file_path):
    data = pd.read_csv(file_path, sep=",", header=None)
    data.columns = data.iloc[0]
    data = data.drop(data.index[0])
    data = data.set_index(data.columns[0])
    return data









#run_blade("deg_mean200.csv", "deg_sd200.csv", "deg200", small_number=0.00001, expected_file="Expected_manual.csv")

#run_blade("deg_mean400.csv", "deg_sd400.csv", "deg400", small_number=0.00001)

#run_blade("deg_mean600.csv", "deg_sd600.csv", "deg600", small_number=0.00001)

#run_blade("deg_mean_norp.csv", "deg_mean_norp.csv", "deg_norp", small_number=0.00001)

#run_blade("deg_mean_ag200.csv", "deg_sd_ag200.csv", "deg_ag200", small_number=0.00001)

#run_blade("deg_mean_ag200_2.csv", "deg_sd_ag200_2.csv", "deg_ag200_2", small_number=0.00001)

#run_blade("deg_mean_ag400_2.csv", "deg_sd_ag400_2.csv", "deg_ag400_2", small_number=0.00001)

#run_blade("deg_mean_ag400.csv", "deg_sd_ag400.csv", "deg_ag400", small_number=0.00001)


#run_blade("deg_mean_norp.csv", "deg_sd_norp.csv", "baseline1.csv","BLADE_baseline1sp", small_number=0.00001, expected_file= "identity_matrix_baseline1.csv")


run_blade("deg_mean_norp.csv", "deg_sd_norp.csv", "baseline2.csv","BLADE_baseline2sp", small_number=0.00001, expected_file= "identity_matrix_baseline2.csv")


run_blade("deg_mean_norp.csv", "deg_sd_norp.csv", "baseline3.csv","BLADE_baseline3sp", small_number=0.00001, expected_file= "identity_matrix_baseline3.csv")


run_blade("deg_mean_norp.csv", "deg_sd_norp.csv", "target_medium.csv","BLADE_mediumsp", small_number=0.00001, expected_file= "identity_matrix_medium.csv")

run_blade("deg_mean_norp.csv", "deg_sd_norp.csv", "target_hard.csv","BLADE_hardsp", small_number=0.00001, expected_file= "identity_matrix_hard.csv")




#run_blade("deg_mean_norp.csv", "deg_sd_norp.csv", "baseline2.csv","BLADE_baseline2_group", small_number=0.00001)


#run_blade("deg_mean_norp.csv", "deg_sd_norp.csv", "baseline3.csv","BLADE_baseline3_group", small_number=0.00001)


#run_blade("deg_mean_ag600.csv", "deg_sd_ag600.csv", "SpatialDataMOB.csv","deg_ag600", small_number=0.00001)


#run_blade("deg_mean_sigdiff_norp_3SD_200.csv", "deg_sd_sigdiff_norp_3SD_200.csv", "deg_sigdiff_norp_3SD_200", small_number=0.00001)




