import os
import pandas as pd

# List of directories
directories = ['maf_low_coverage', 'maf_regular_coverage', 'no_maf_low_coverage', 'no_maf_regular_coverage']
base_path = '/UniqueHaplotypesFractions/'

# Loop over each directory
for dir_name in directories:
    directory_path = os.path.join(base_path, dir_name)
    all_files = os.listdir(directory_path)
    csv_files = [f for f in all_files if f.endswith('.csv')]

    dfs = []
    for file in csv_files:
        file_path = os.path.join(directory_path, file)
        df = pd.read_csv(file_path, index_col=0)  # Assumes "Window Size" is the index
        dfs.append(df)

    panel_500 = []
    panel_2500 = []

    for df in dfs:
        panel_500.append(df['500'])
        panel_2500.append(df['2500'])

    panel_500_concat = pd.concat(panel_500, axis=1)
    panel_2500_concat = pd.concat(panel_2500, axis=1)

    # Median, Min, Max for Panel 500
    panel_500_median = panel_500_concat.median(axis=1)
    panel_500_min = panel_500_concat.min(axis=1)
    panel_500_max = panel_500_concat.max(axis=1)

    # Median, Min, Max for Panel 2500
    panel_2500_median = panel_2500_concat.median(axis=1)
    panel_2500_min = panel_2500_concat.min(axis=1)
    panel_2500_max = panel_2500_concat.max(axis=1)

    stats_df = pd.DataFrame({
        "Window Size": panel_500_median.index,  # The index is the window size
        "500 Median": panel_500_median.values,
        "500 Min": panel_500_min.values,
        "500 Max": panel_500_max.values,
        "2500 Median": panel_2500_median.values,
        "2500 Min": panel_2500_min.values,
        "2500 Max": panel_2500_max.values
    }).set_index("Window Size")

    # Save the stats DataFrame with the directory name in the filename
    output_filename = f'stats_{dir_name.strip("/")}.csv'
    stats_df.to_csv(os.path.join(base_path, output_filename))

    print(f"Stats for {dir_name} saved to {output_filename}")
