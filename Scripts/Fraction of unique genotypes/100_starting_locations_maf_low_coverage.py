# plotting window size vs unique genotypes, but now we are using 100 different starting locations
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


#Generate the 100 different starting locations
# Total number of SNPs
total_snps = 933347

# Number of starting locations to select
num_locations = 100

# Calculate the spacing between starting locations
spacing = total_snps // num_locations

# Generate starting locations
starting_locations = [i * spacing for i in range(num_locations)]

# Define panel sizes and row counts
panel_sizes = [500, None]  # None for full dataset
row_counts = [1] + list(range(5, 205, 5))

for loc in range(0,100):
    results = {}
    chunk_size = 100000
    target_rows = 200
    start_idx = starting_locations[loc]
    selected_rows = pd.DataFrame()  # To store selected rows

    # Initialize the chunk reader
    chunks = pd.read_csv('filtered_genotype_matrix.txt', delimiter='\t', header=None, chunksize=chunk_size)

    # Calculate the starting chunk to process
    start_chunk = start_idx // chunk_size

    # Process chunks
    for i, chunk in enumerate(chunks):
        if i < start_chunk:
            continue  # Skip chunks before the start_chunk
        # Determine the actual row indices in the current chunk if it's the start chunk
        if i == start_chunk:
            chunk = chunk.iloc[start_idx - i * chunk_size:]
        # Randomly select rows with a 10% chance
        selected = chunk[np.random.rand(chunk.shape[0]) < 0.1]
        # Append the selected rows to the accumulator DataFrame
        selected_rows = pd.concat([selected_rows, selected], ignore_index=True)
        # Check if we have accumulated the target number of rows
        if selected_rows.shape[0] >= target_rows:
            selected_rows = selected_rows.iloc[:target_rows]  # Trim excess rows
            break
    # Results now contain exactly 200 rows with a 10% selection probability
    df_subset = selected_rows

    # Process each panel size
    for panel_size in panel_sizes:
        fraction_of_unique_haplotypes_results = []
        for nrows in row_counts:
            df_subset_subset = df_subset.head(nrows).reset_index(drop=True)
            # Extract columns from index 4 to the end
            sample_columns = df_subset_subset.iloc[:, 4:]
            # Create a new DataFrame with the first 4 columns
            new_df = df_subset_subset.iloc[:, :4].copy()
            # Create a dictionary to collect new columns
            columns_dict = {}
            # Iterate through each row to populate columns_dict with sample names and values

            for i, row in df_subset_subset.iterrows():
                for col in sample_columns.columns:
                    # Extract sample name and value
                    sample_name, value = row[col].split('=')
                    sample_name = sample_name.strip()
                    value = value.strip()
                    # Check if the sample name already exists in columns_dict
                    if sample_name not in columns_dict:
                        columns_dict[sample_name] = [None] * len(df_subset_subset)  # Initialize with None values
                    # Set the value for the current row
                    columns_dict[sample_name][i] = value

            # Convert columns_dict to a DataFrame
            new_samples_df = pd.DataFrame(columns_dict)
            # Concatenate the original DataFrame with the new columns
            final_df = pd.concat([new_df, new_samples_df], axis=1)
            
            # Define a function to map genotype values to the desired output
            def map_genotype(genotype):
                if genotype == '0|0':
                    return 0
                elif genotype in ['0|1', '1|0']:
                    return 1
                elif genotype == '1|1':
                    return 2
                else:
                    return None  # Or handle unexpected cases
                
            for col in final_df.columns[4:]:  # Adjust if your columns are 0-indexed
                final_df[col] = final_df[col].apply(map_genotype)

            # Get the joint genotype vectors for the samples
            result = final_df.iloc[:, 4:].apply(lambda x: ' '.join(x.astype(str)), axis=0)
            
            if panel_size is not None:
                # Sample the result if a panel size is specified
                sampled_result = result.sample(n=panel_size, random_state=1)  # Set random_state for reproducibility if needed
            else:
                # Use all individuals (no sampling) for the full dataset case
                sampled_result = result

            # Count the number of times the joint genotype vectors appear
            unique_counts = sampled_result.value_counts()
            # Take the number of unique haplotypes
            count_of_lines_with_single_occurrence = unique_counts[unique_counts == 1].count()
            fraction_of_unique_haplotypes = count_of_lines_with_single_occurrence / len(unique_counts)
            fraction_of_unique_haplotypes_results.append(fraction_of_unique_haplotypes)

        # Store the results for the current panel size
        results[panel_size] = fraction_of_unique_haplotypes_results

    column_names = [panel_size if panel_size is not None else 2500 for panel_size in panel_sizes]
    results_df = pd.DataFrame(results, index=row_counts)
    results_df.columns = column_names
    output_file = f"/maf_low_coverage/haplotype_fraction_starting_location_{starting_locations[loc]}.csv"
    results_df.to_csv(output_file, index_label='Window Size')

