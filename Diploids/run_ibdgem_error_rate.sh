#!/bin/bash

# Define the directories containing the subdirectories
main_directories=(
    "/error_rate/n100_k50_r2_order1"
    "/error_rate/n100_k50_r2_order2"
    "/error_rate/n100_k50_r2_order3"
    "/error_rate/n100_k50_r2_order4"
    "/error_rate/n100_k50_r2_order5"
    "/error_rate/n100_k50_r2_order6"
    "/error_rate/n100_k50_r2_order7"
)

# Define the error rates corresponding to each order number
error_rates=(0.0001 0.001 0.01 0.02 0.05 0.1 0.2)

# Loop through each main directory
for main_directory in "${main_directories[@]}"; do
    # Get the value of n from the directory name
    n_value=$(basename "$main_directory" | cut -d 'n' -f 2 | cut -d '_' -f 1)
    order_number=$(basename "$main_directory" | awk -F'order' '{print $2}')
    order_number=$((order_number + 0))

    # Loop through each subdirectory in the main directory
    for subdir in "$main_directory"/*; do
        # Check if the subdirectory contains the necessary files
        if [[ -e "$subdir/test.hap" && -e "$subdir/test.legend" && -e "$subdir/test.indv" && -e "$subdir/unknowndna.pileup" ]]; then
            # Determine the value of -N based on the value of n
            sample_value="sample$n_value"
            # Get the corresponding error rate for the order number
            error_rate=${error_rates[$((order_number - 1))]}
            # Run the ibdgem command
            ./IBDGem-2.0.2/ibdgem -H "$subdir/test.hap" -L "$subdir/test.legend" -I "$subdir/test.indv" -P "$subdir/unknowndna.pileup" -N "$sample_value" -O "$subdir/output_non_LD" --allele-freqs "$subdir/output.frq" --error-rate "$error_rate"
            ./IBDGem-2.0.2/ibdgem --LD -H "$subdir/test.hap" -L "$subdir/test.legend" -I "$subdir/test.indv" -P "$subdir/unknowndna.pileup" -N "$sample_value" -O "$subdir/output_LD" --allele-freqs "$subdir/output.frq" --error-rate "$error_rate"
        else
            echo "Files not found in $subdir"
        fi
    done
done
