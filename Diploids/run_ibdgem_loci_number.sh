#!/bin/bash

# Define the directories containing the subdirectories
main_directories=(
    "/loci_number/n100_k1_r2_order1"
    "/loci_number/n100_k25_r2_order2"
    "/loci_number/n100_k50_r2_order3"
    "/loci_number/n100_k75_r2_order4"
    "/loci_number/n100_k100_r2_order5"
)

# Loop through each main directory
for main_directory in "${main_directories[@]}"; do
    # Get the value of n from the directory name
    n_value=$(basename "$main_directory" | cut -d 'n' -f 2 | cut -d '_' -f 1)
    # Get the value of k from the directory name
    k_value=$(basename "$main_directory" | cut -d 'k' -f 2 | cut -d '_' -f 1)

    # Loop through each subdirectory in the main directory
    for subdir in "$main_directory"/*; do
        # Check if the subdirectory contains the necessary files
        if [[ -e "$subdir/test.hap" && -e "$subdir/test.legend" && -e "$subdir/test.indv" && -e "$subdir/unknowndna.pileup" ]]; then
            # Determine the value of -N based on the value of n
            sample_value="sample$n_value"

            # Run the ibdgem command
            if [ "$k_value" -eq 1 ]; then
                ./IBDGem-2.0.2/ibdgem -H "$subdir/test.hap" -L "$subdir/test.legend" -I "$subdir/test.indv" -P "$subdir/unknowndna.pileup" -N "$sample_value" -O "$subdir/output_non_LD" --allele-freqs "$subdir/output.frq"
                ./IBDGem-2.0.2/ibdgem --LD -H "$subdir/test.hap" -L "$subdir/test.legend" -I "$subdir/test.indv" -P "$subdir/unknowndna.pileup" -N "$sample_value" -O "$subdir/output_LD" --allele-freqs "$subdir/output.frq"
            else
                ./IBDGem-2.0.2/ibdgem -H "$subdir/test.hap" -L "$subdir/test.legend" -I "$subdir/test.indv" -P "$subdir/unknowndna.pileup" -N "$sample_value" -O "$subdir/output_non_LD" --allele-freqs "$subdir/output.frq" --window-size "$k_value"
                ./IBDGem-2.0.2/ibdgem --LD -H "$subdir/test.hap" -L "$subdir/test.legend" -I "$subdir/test.indv" -P "$subdir/unknowndna.pileup" -N "$sample_value" -O "$subdir/output_LD" --allele-freqs "$subdir/output.frq" --window-size "$k_value"
            fi
        else
            echo "Files not found in $subdir"
        fi
    done
done
